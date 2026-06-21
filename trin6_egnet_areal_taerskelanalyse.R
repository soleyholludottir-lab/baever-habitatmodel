# ============================================================
# TRIN 6 (SAMLET): TÆRSKLER OG EGNET HABITAT
# Castor fiber, Sjælland og Lolland-Falster, Danmark (2026)
# Baseret på 5-variabel modellen med alle arealanvendelsesklasser
# (v5 _alle)
#
# Dette script samler tre tidligere separate scripts:
#   - trin6a (areal pr. vandopland/kommune/Natura 2000, fast tærskel)
#   - trin6b (tærskelberegning + samlet/per-vandopland areal)
#   - Trin6f (tærskelfølsomhed: areal ved flere tærskler)
#
# Formål:
#   1. Find alle relevante tærskler ud fra habitategnethed ved de
#      55 NST-validerede forekomstpunkter (direkte punktekstraktion,
#      ingen buffer, WGS84-raster — matcher QGIS "Sample Raster Values")
#   2. Beregn egnet areal (km²) konsekvent med terra::cellSize() på
#      WGS84-rasteret 
#   3. Beregn areal for både gennemsnitstærsklen (Swinnen et al. 2017) og
#      10. percentilen, samt giv mulighed for en brugerdefineret tærskel
#   4. Bryd egnet areal ned pr. vandopland, kommune og Natura 2000-område
#
# OBS — ændringer i forhold til de oprindelige tre scripts:
#   - EPSG:25832-reprojektion + prod(res()) (brugt i trin6b/Trin6f) er
#     droppet til fordel for cellSize() på WGS84 (brugt i trin6a) —
#     ét konsekvent areal-regnestykke i hele scriptet.
#   - NoData-rensningen fra trin6b (fjerner -3.4e38-flag og falske
#     0-celler) anvendes nu konsekvent, også i polygon-opdelingen,
#     hvor den oprindeligt manglede i trin6a.
#   - occ_coords_wgs84.csv (brugt i Trin6f) er droppet til fordel for
#     den samme gpkg-kilde som trin6b ("Baeverbo til Soley.gpkg").
# ============================================================

# === load packages ===
library(terra)
library(sf)

# === stier ===
lag_sti    <- "C:/qGIS/bachelor/layers_endelig/"
github_sti <- "C:/Users/Soley/Documents/GitHub/baever-habitatmodel/"
  
  # ============================================================
# TRIN 6.1: INDLÆS HABITATKORT (WGS84) OG FORBERED AREALBEREGNING
# ============================================================

hab <- rast(paste0(lag_sti, "habitatkort_v5_5var_alle.tif"))
hab[hab < 0]  <- NA   # fjern NoData-flag (-3.4e+38)
hab[hab == 0] <- NA   # fjern falske 0-celler

# cellSize() håndterer WGS84 geodetisk korrekt og bruges til ALLE
# arealberegninger i dette script
csz <- cellSize(hab, unit = "km")          # km² pr. celle
total_km2_cell <- (!is.na(hab)) * csz      # km² hvor der er data

tot_areal <- global(total_km2_cell, "sum", na.rm = TRUE)[1, 1]
cat("Samlet studieområde (km²):", round(tot_areal, 1), " (forventet: ~9272)\n")

# RESULTAT: Samlet studieområde (km²): 9272.4  (forventet: ~9272)

# === indlæs forekomstpunkter og vandoplande ===
occ <- st_read(paste0(lag_sti, "Baeverbo til Soley.gpkg"),   # Naturstyrelsens data
               layer = "bverbo2026sg")
vandoplande <- st_read(paste0(lag_sti, "hovedvandoplande.gpkg"),
                       layer = "vp3endelig2022vp3e2022_hovedoplande")

# studieområdets afgrænsning - bruges til at filtrere kommuner og
# N2000-områder til kun dem der overlapper studieområdet
studieomr_union <- st_union(vandoplande)

occ_wgs84 <- st_transform(occ, crs = "EPSG:4326")  # til direkte punktekstraktion

cat("Habitatkort CRS:    ", crs(hab, describe = TRUE)$code, "\n")
# RESULTAT: Habitatkort CRS:     4326 
cat("Habitatkort opløsning:", res(hab), "\n")
# RESULTAT: Habitatkort opløsning: 0.000112945 0.000112945 
cat("Antal ikke-NA celler:", global(hab, "notNA")$notNA, "\n")
# RESULTAT: Antal ikke-NA celler: 103069194 
cat("Antal occurrence:    ", nrow(occ_wgs84), "\n")
# RESULTAT: Antal occurrence:     55 

# ============================================================
# TRIN 6.2: UDTRÆK EGNETHED VED FOREKOMSTPUNKTER
# ============================================================
# Direkte punktekstraktion uden buffer fra original WGS84-raster —
# matcher QGIS "Sample Raster Values" og undgår bilineær interpolation

suit_vaerdier <- terra::extract(hab, vect(occ_wgs84), ID = FALSE)
suit_ok <- suit_vaerdier[!is.na(suit_vaerdier[, 1]), 1]

cat("\n--- Habitategnethed ved forekomstpunkter (5-variabel model, alle klasser) ---\n")
cat("  Antal punkter med værdi:", length(suit_ok), "/", nrow(occ_wgs84), "\n")
cat("  NA:                     ", nrow(occ_wgs84) - length(suit_ok), "\n")
print(round(summary(suit_ok), 3))

# ============================================================
# TRIN 6.3 (DIAGNOSTIK): LAV-SUITABILITY PUNKTER
# ============================================================

occ_wgs84$suitability <- suit_vaerdier[, 1]

lav_suit <- occ_wgs84[!is.na(occ_wgs84$suitability) &
                        occ_wgs84$suitability < 0.1, ]

cat("\nAntal punkter med suitability < 0.1:", nrow(lav_suit), "\n")
# RESULTAT: Antal punkter med suitability < 0.1: 1 
# Bæverbo nr. 45 i datasættet med suitability 0.06089008

if (nrow(lav_suit) > 0) {
  cat("Attributter for disse punkter:\n")
  print(st_drop_geometry(lav_suit))
}

# --- hurtigt kort i R ---
plot(hab,
     main = "Egnethed ved forekomstpunkter (5-variabel model, alle klasser)",
     col  = hcl.colors(100, "viridis"))

plot(st_geometry(occ_wgs84), add = TRUE,
     pch = 16, col = "white", cex = 0.8)

plot(st_geometry(lav_suit), add = TRUE,
     pch = 16, col = "red", cex = 1.2)

legend("bottomright",
       legend = c("NST-punkt (ok)", "Suitability < 0.1"),
       col    = c("white", "red"),
       pch    = 16)

# --- eksporter til QGIS ---
st_write(occ_wgs84,
         dsn   = paste0(lag_sti, "occ_med_suitability_5var_alle.gpkg"),
         layer = "occ_suitability",
         delete_layer = TRUE)
cat("Geopackage gemt: occ_med_suitability_5var_alle.gpkg\n")

# ============================================================
# TRIN 6.4: BEREGN ALLE TÆRSKLER
# ============================================================
# Metodisk begrundelse: hovedtærsklen beregnes som gennemsnits-
# egnethed ved forekomstpunkterne (Swinnen et al. 2017, mean = 0.54
# i Dijle-dalen). Her suppleres med 95 % bootstrap-konfidensinterval
# samt median og 10. percentil som referencetærskler.

set.seed(42)   # projektets seed fra trin 3
boot_mean <- replicate(2000, mean(sample(suit_ok, replace = TRUE)))
ci <- quantile(boot_mean, c(0.025, 0.975))

thresholds <- c(
  "CI nedre (2.5%)"      = unname(ci[1]),
  "Gennemsnit (anvendt)" = mean(suit_ok),
  "CI øvre (97.5%)"      = unname(ci[2]),
  "Median"               = median(suit_ok),
  "10. percentil"        = unname(quantile(suit_ok, 0.10))
)

cat("\n--- Tærskelværdier (5-variabel model, alle klasser) ---\n")
print(round(thresholds, 3))

# --- Tærskelværdier (5-variabel model, alle klasser) ---
# CI nedre (2.5%)   Gennemsnit (anvendt)    CI øvre (97.5%)     Median        10. percentil 
#           0.669                  0.740              0.807      0.860                0.305 

# === de to hovedtærskler, som arealberegningerne nedenfor kører på ===
taerskel_gennemsnit <- unname(thresholds["Gennemsnit (anvendt)"])
taerskel_p10    <- unname(thresholds["10. percentil"])

# === valgfri brugerdefineret tærskel ===
# Sæt en ønsket værdi her for at få den regnet med i areal- og
# polygonopgørelserne nedenfor (TRIN 6.5/6.6). Sæt til NA for at
# springe over.
taerskel_egen <- NA   # fx 0.5

# ============================================================
# TRIN 6.5: SAMLET EGNET AREAL (STUDIEOMRÅDE)
# ============================================================

egnet_areal_km2 <- function(taerskel) {
  global((hab >= taerskel) * csz, "sum", na.rm = TRUE)[1, 1]
}

areal_gennemsnit <- egnet_areal_km2(taerskel_gennemsnit)
areal_p10    <- egnet_areal_km2(taerskel_p10)

cat("\n--- Samlet egnet areal (cellSize, WGS84) ---\n")
cat("Studieområde i alt:                  ", round(tot_areal, 1), "km²\n")
cat("Egnet ved gennemsnit (", round(taerskel_gennemsnit, 3), "): ",
    round(areal_gennemsnit, 1), "km² (", round(areal_gennemsnit / tot_areal * 100, 1), "%)\n")
cat("Egnet ved 10. percentil (", round(taerskel_p10, 3), "):    ",
    round(areal_p10, 1), "km² (", round(areal_p10 / tot_areal * 100, 1), "%)\n")

if (!is.na(taerskel_egen)) {
  areal_egen <- egnet_areal_km2(taerskel_egen)
  cat("Egnet ved brugerdefineret tærskel (", round(taerskel_egen, 3), "):",
      round(areal_egen, 1), "km²\n")
}

# === binære kort gemmes for begge hovedtærskler ===
binaert_gennemsnit <- hab >= taerskel_gennemsnit
binaert_p10    <- hab >= taerskel_p10

writeRaster(binaert_gennemsnit,
            filename  = paste0(lag_sti, "habitatkort_binaert_5var_alle_taerskel_gennemsnit.tif"),
            overwrite = TRUE, datatype = "INT1U")
writeRaster(binaert_p10,
            filename  = paste0(lag_sti, "habitatkort_binaert_5var_alle_taerskel_p10.tif"),
            overwrite = TRUE, datatype = "INT1U")
cat("\nBinære kort gemt (gennemsnits og 10.percentil-tærskel).\n")

if (!is.na(taerskel_egen)) {
  binaert_egen <- hab >= taerskel_egen
  writeRaster(binaert_egen,
              filename  = paste0(lag_sti, "habitatkort_binaert_5var_alle_taerskel_egen.tif"),
              overwrite = TRUE, datatype = "INT1U")
  cat("Binært kort gemt (brugerdefineret tærskel).\n")
}

# ============================================================
# TRIN 6.6: EGNET AREAL PR. POLYGON (VANDOPLAND / KOMMUNE / N2000)
# ============================================================
# Generisk funktion, kører for én tærskel ad gangen

areal_pr_polygon <- function(polygon_sf, navn_felt, taerskel) {
  
  # --- diagnostik: findes feltnavnet? ---
  if (!(navn_felt %in% names(polygon_sf))) {
    stop("navn_felt '", navn_felt, "' findes ikke. Tilgaengelige felter: ",
         paste(names(polygon_sf), collapse = ", "))
  }
  
  # --- reparer evt. ugyldige geometrier ---
  ugyldige <- !st_is_valid(polygon_sf)
  if (any(ugyldige, na.rm = TRUE)) {
    cat("OBS:", sum(ugyldige, na.rm = TRUE), "ugyldige geometrier - reparerer med st_make_valid()\n")
    polygon_sf <- st_make_valid(polygon_sf)
  }
  
  # --- filtrer til kun polygoner der overlapper studieomraadet ---
  # (løser "extents do not overlap" - kommuner/N2000 daekker hele DK,
  #  men habitatkortet daekker kun Sjaelland/Lolland-Falster)
  polygon_sf <- st_transform(polygon_sf, crs = st_crs(studieomr_union))
  n_foer <- nrow(polygon_sf)
  polygon_sf <- st_filter(polygon_sf, st_as_sf(studieomr_union), .predicate = st_intersects)
  cat("Filtreret fra", n_foer, "til", nrow(polygon_sf), "polygoner der overlapper studieomraadet\n")
  
  pol <- st_transform(polygon_sf, crs = crs(hab))   # polygoner -> rasterets CRS (WGS84)
  
  egnet_km2_cell <- (hab >= taerskel) * csz
  
  out <- data.frame(navn = character(), total_km2 = numeric(),
                    egnet_km2 = numeric(), andel_pct = numeric(),
                    stringsAsFactors = FALSE)
  for (i in seq_len(nrow(pol))) {
    pv  <- vect(pol[i, ])
    
    # tryCatch som sikkerhedsnet for kantsager (sliver-overlap der
    # bestaar bbox-filteret i st_filter, men reelt ikke overlapper rasteret)
    e <- tryCatch({
      eg  <- mask(crop(egnet_km2_cell, pv), pv)
      val <- global(eg, "sum", na.rm = TRUE)[1, 1]
      if (is.na(val)) 0 else val
    }, error = function(err) {
      cat("  -> ingen overlap for", as.character(pol[[navn_felt]][i]), "- saetter 0\n")
      0
    })
    
    t <- tryCatch({
      tot <- mask(crop(total_km2_cell, pv), pv)
      val <- global(tot, "sum", na.rm = TRUE)[1, 1]
      if (is.na(val)) 0 else val
    }, error = function(err) 0)
    
    out[i, ] <- list(
      navn      = as.character(pol[[navn_felt]][i]),
      total_km2 = round(t, 1),
      egnet_km2 = round(e, 1),
      andel_pct = if (t > 0) round(e / t * 100, 2) else NA_real_
    )
    cat("Faerdig:", out$navn[i], "\n")
  }
  out
}

# === kør for begge hovedtærskler og saml i én tabel pr. polygontype ===
areal_pr_polygon_begge <- function(polygon_sf, navn_felt) {
  res_gennemsnit <- areal_pr_polygon(polygon_sf, navn_felt, taerskel_gennemsnit)
  res_gennemsnit$taerskel_type   <- "gennemsnit"
  res_gennemsnit$taerskel_vaerdi <- round(taerskel_gennemsnit, 3)
  
  res_p10 <- areal_pr_polygon(polygon_sf, navn_felt, taerskel_p10)
  res_p10$taerskel_type   <- "p10"
  res_p10$taerskel_vaerdi <- round(taerskel_p10, 3)
  
  res <- rbind(res_gennemsnit, res_p10)
  
  if (!is.na(taerskel_egen)) {
    res_egen <- areal_pr_polygon(polygon_sf, navn_felt, taerskel_egen)
    res_egen$taerskel_type   <- "egen"
    res_egen$taerskel_vaerdi <- round(taerskel_egen, 3)
    res <- rbind(res, res_egen)
  }
  
  res
}

# --- 1) vandoplande ---
res_vand <- areal_pr_polygon_begge(vandoplande, navn_felt = "hov_na")
cat("\n--- Egnet habitat pr. vandopland ---\n")
print(res_vand, row.names = FALSE)
write.csv(res_vand, paste0(github_sti, "v5/egnet_pr_vandopland.csv"), row.names = FALSE)

# RESULTAT: --- Egnet habitat pr. vandopland ---
# navn                      total_km2 egnet_km2 andel_pct taerskel_type taerskel_vaerdi
# Kalundborg                  985.0       3.4      0.35    gennemsnit           0.740
# Isefjord og Roskilde Fjord 1948.4      14.0      0.72    gennemsnit           0.740
# Øresund                     939.8       6.7      0.72    gennemsnit           0.740
# Køge Bugt                   874.0       2.1      0.25    gennemsnit           0.740
# Smålandsfarvandet          3279.5      18.8      0.57    gennemsnit           0.740
# Østersøen                  1254.5       3.3      0.26    gennemsnit           0.740
#-------------------------------------------------------------------------------------
# Kalundborg                  985.0      18.2      1.85           p10           0.305
# Isefjord og Roskilde Fjord 1948.4      55.2      2.83           p10           0.305
# Øresund                     939.8      32.0      3.40           p10           0.305
# Køge Bugt                   874.0       8.6      0.99           p10           0.305
# Smålandsfarvandet          3279.5      62.3      1.90           p10           0.305
# Østersøen                  1254.5      11.8      0.94           p10           0.305

# --- 2) kommuner ---
kommuner <- st_read(paste0(lag_sti, "kommuner.gpkg"), layer = "matmatrikelkommune_gaeldende_co")
res_komm <- areal_pr_polygon_begge(kommuner, navn_felt = "kommunenavn")
cat("\n--- Egnet habitat pr. kommune ---\n")
print(res_komm, row.names = FALSE)
write.csv(res_komm, paste0(github_sti, "v5/egnet_pr_kommune.csv"), row.names = FALSE)

# RESULTAT: --- Egnet habitat pr. kommune ---

# navn                total_km2 egnet_km2 andel_pct taerskel_type taerskel_vaerdi
# Vallensbæk Kommune       9.3       0.0      0.41    gennemsnit           0.740
# Lejre Kommune          239.6       0.9      0.40    gennemsnit           0.740
# Albertslund Kommune     23.6       0.2      0.66    gennemsnit           0.740
# Lyngby-Taarbæk Kommune  39.0       0.4      1.01    gennemsnit           0.740
# Gladsaxe Kommune        25.1       0.1      0.33    gennemsnit           0.740
# Rødovre Kommune         12.3       0.1      0.90    gennemsnit           0.740
# Frederiksberg Kommune    8.8       0.0      0.01    gennemsnit           0.740
# Tårnby Kommune          65.2       0.0      0.02    gennemsnit           0.740
# Ringsted Kommune       295.3       0.4      0.15    gennemsnit           0.740
# Faxe Kommune           405.6       0.6      0.14    gennemsnit           0.740
# Lolland Kommune        886.2       5.0      0.56    gennemsnit           0.740
# Dragør Kommune          18.2       0.0      0.06    gennemsnit           0.740
# Allerød Kommune         67.7       0.0      0.01    gennemsnit           0.740
# Næstved Kommune        678.0       4.4      0.65    gennemsnit           0.740
# Sorø Kommune           309.3       0.4      0.14    gennemsnit           0.740
# Egedal Kommune         126.2       1.8      1.46    gennemsnit           0.740
# Ballerup Kommune        34.1       0.3      0.75    gennemsnit           0.740
# Hvidovre Kommune        22.9       0.0      0.18    gennemsnit           0.740
# Kalundborg Kommune     576.7       2.7      0.47    gennemsnit           0.740
# Guldborgsund Kommune   901.1       5.7      0.63    gennemsnit           0.740
# Odsherred Kommune      354.9       1.1      0.31    gennemsnit           0.740
# Helsingør Kommune      119.0       0.6      0.54    gennemsnit           0.740
# Rudersdal Kommune       73.6       0.6      0.84    gennemsnit           0.740
# Frederikssund Kommune  249.1       1.3      0.51    gennemsnit           0.740
# Hørsholm Kommune        31.5       0.2      0.60    gennemsnit           0.740
# Slagelse Kommune       569.4       1.5      0.27    gennemsnit           0.740
# Greve Kommune           60.3       0.3      0.56    gennemsnit           0.740
# Ishøj Kommune           26.3       0.2      0.61    gennemsnit           0.740
# Holbæk Kommune         578.1       2.2      0.38    gennemsnit           0.740
# Solrød Kommune          40.3       0.1      0.22    gennemsnit           0.740
# Herlev Kommune          12.2       0.0      0.32    gennemsnit           0.740
# Furesø Kommune          57.1       0.8      1.41    gennemsnit           0.740
# Hillerød Kommune       213.9       2.7      1.26    gennemsnit           0.740
# Brøndby Kommune         21.0       0.1      0.36    gennemsnit           0.740
# Glostrup Kommune        13.4       0.1      0.74    gennemsnit           0.740
# Fredensborg Kommune    112.4       1.0      0.93    gennemsnit           0.740
# Køge Kommune           257.1       0.7      0.25    gennemsnit           0.740
# Roskilde Kommune       212.5       1.1      0.49    gennemsnit           0.740
# Københavns Kommune      89.5       0.2      0.17    gennemsnit           0.740
# Gentofte Kommune        25.7       0.1      0.50    gennemsnit           0.740
# Stevns Kommune         250.2       0.7      0.26    gennemsnit           0.740
# Halsnæs Kommune        121.3       1.2      0.97    gennemsnit           0.740
# Høje-Taastrup Kommune   78.5       0.6      0.83    gennemsnit           0.740
# Gribskov Kommune       279.6       4.4      1.59    gennemsnit           0.740
# Vordingborg Kommune    620.2       1.7      0.28    gennemsnit           0.740
#------------------------------------------------------------------------------
# Vallensbæk Kommune       9.3       0.3      2.89           p10           0.305
# Lejre Kommune          239.6       3.6      1.51           p10           0.305
# Albertslund Kommune     23.6       0.5      1.95           p10           0.305
# Lyngby-Taarbæk Kommune  39.0       1.7      4.39           p10           0.305
# Gladsaxe Kommune        25.1       0.6      2.36           p10           0.305
# Rødovre Kommune         12.3       0.1      1.08           p10           0.305
# Frederiksberg Kommune    8.8       0.0      0.22           p10           0.305
# Tårnby Kommune          65.2       0.2      0.26           p10           0.305
# Ringsted Kommune       295.3       4.1      1.39           p10           0.305
# Faxe Kommune           405.6       2.5      0.63           p10           0.305
# Lolland Kommune        886.2      13.7      1.55           p10           0.305
# Dragør Kommune          18.2       0.1      0.57           p10           0.305
# Allerød Kommune         67.7       0.4      0.59           p10           0.305
# Næstved Kommune        678.0      12.0      1.77           p10           0.305
# Sorø Kommune           309.3       3.5      1.15           p10           0.305
# Egedal Kommune         126.2       6.7      5.27           p10           0.305
# Ballerup Kommune        34.1       1.0      2.82           p10           0.305
# Hvidovre Kommune        22.9       0.3      1.26           p10           0.305
# Kalundborg Kommune     576.7      11.8      2.05           p10           0.305
# Guldborgsund Kommune   901.1      17.2      1.91           p10           0.305
# Odsherred Kommune      354.9       4.8      1.35           p10           0.305
# Helsingør Kommune      119.0       3.5      2.91           p10           0.305
# Rudersdal Kommune       73.6       3.4      4.59           p10           0.305
# Frederikssund Kommune  249.1       4.6      1.83           p10           0.305
# Hørsholm Kommune        31.5       1.2      3.89           p10           0.305
# Slagelse Kommune       569.4       6.0      1.05           p10           0.305
# Greve Kommune           60.3       1.0      1.69           p10           0.305
# Ishøj Kommune           26.3       0.7      2.48           p10           0.305
# Holbæk Kommune         578.1       7.6      1.31           p10           0.305
# Solrød Kommune          40.3       0.5      1.16           p10           0.305
# Herlev Kommune          12.2       0.3      2.15           p10           0.305
# Furesø Kommune          57.1       4.1      7.27           p10           0.305
# Hillerød Kommune       213.9      11.6      5.43           p10           0.305
# Brøndby Kommune         21.0       0.3      1.36           p10           0.305
# Glostrup Kommune        13.4       0.2      1.41           p10           0.305
# Fredensborg Kommune    112.4       4.8      4.30           p10           0.305
# Køge Kommune           257.1       2.7      1.03           p10           0.305
# Roskilde Kommune       212.5       4.7      2.23           p10           0.305
# Københavns Kommune      89.5       1.3      1.46           p10           0.305
# Gentofte Kommune        25.7       0.5      1.84           p10           0.305
# Stevns Kommune         250.2       2.1      0.83           p10           0.305
# Halsnæs Kommune        121.3       6.7      5.50           p10           0.305
# Høje-Taastrup Kommune   78.5       2.4      3.06           p10           0.305
# Gribskov Kommune       279.6      15.9      5.68           p10           0.305
# Vordingborg Kommune    620.2       6.3      1.01           p10           0.305

# --- 3) Natura 2000 ---
n2000 <- st_read(paste0(lag_sti, "natura2000.gpkg"), layer = "natura_2000_omrder")
res_n2000 <- areal_pr_polygon_begge(n2000, navn_felt = "n2000_nr")
cat("\n--- Egnet habitat pr. Natura 2000-område ---\n")
print(res_n2000, row.names = FALSE)
write.csv(res_n2000, paste0(github_sti, "v5/egnet_pr_n2000.csv"), row.names = FALSE)



# ============================================================
# TRIN 6.7: TÆRSKELFØLSOMHED (SAMLET STUDIEOMRÅDE)
# ============================================================
# Kvantificerer hvor følsomt det egnede areal (km²) er over for
# valget af tærskel, ved at anvende alle tærskler fra TRIN 6.4

res_df <- data.frame(
  metode    = names(thresholds),
  taerskel  = round(unname(thresholds), 3),
  egnet_km2 = round(sapply(unname(thresholds), egnet_areal_km2), 1),
  stringsAsFactors = FALSE
)
res_df$andel_pct <- round(res_df$egnet_km2 / tot_areal * 100, 2)

cat("\n--- Tærskelfølsomhed: egnet areal ved forskellige tærskler (cellSize, WGS84) ---\n")
print(res_df, row.names = FALSE)

# RESULTAT: --- Tærskelfølsomhed: egnet areal ved forskellige tærskler (cellSize, WGS84) ---

# metode                taerskel egnet_km2 andel_pct
# CI nedre (2.5%)         0.669      62.4      0.67
# Gennemsnit (anvendt)    0.740      48.4      0.52
# CI øvre (97.5%)         0.807      36.8      0.40
# Median                  0.860      28.0      0.30
# 10. percentil           0.305     188.1      2.03

write.csv(res_df,
          file      = paste0(github_sti, "v5/taerskelfoelsomhed.csv"),
          row.names = FALSE)
cat("\nGemt: v5/taerskelfoelsomhed.csv\n")