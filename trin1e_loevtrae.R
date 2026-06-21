# ============================================================
# TRIN 1f: UDLED DIST_LOEVTRAE.TIF + KORRELATIONSDIAGNOSTIK
# Formål: (a) Beregn afstand til løvtræ (lu_leaf_type_2024, Basemap05) som
#         et kontinuert miljøvariabel, og (b) test korrelation
#         med eksisterende lag inden inkludering i modellen (Dormann et al. 2013).
# Metode: Binær maske (løvtræs-klasse) → terra::distance() →
#         korrelationsmatrix på 10.000 baggrundspunkter.
# Litteratur: Løvtræ (Salix, Populus, Betula) er primær fødekilde
#             for Castor fiber (Jones et al. 2009); vegetation forklarede
#             ~28,8 % af modelydelse i Swinnen et al. (2017).
#             Tilsvarende vegetationsvariabler indgår i bæver-SDM
#             hos Treves et al. (2022).
# Hvad scriptet gør:
## 1. Diagnosticerer lu_leaf_type_2024.tif (klassefordeling, CRS)
## 2. Laver binær maske: klasse 1 (løvtræ) = 1, resten = NA
## 3. Beregner afstand med terra::distance() → dist_loevtrae.tif
## 4. Clipper og alignerer midlertidigt til reference-grid for korrelation
## 5. Sampler 10.000 baggrundspunkter og ekstraherer værdier
## 6. Beregner Pearson-korrelation mellem dist_loevtrae og alle
##    øvrige kontinuerte lag
## 7. Viser middel dist_loevtrae pr. arealanvendelsesklasse
##    (Pearson er ikke meningsfuld for kategorisk variabel)
###
### Basemap05 kilde:
### Levin, G. (2026). Basemap05. Documentation of the data and method for elaboration 
### of a land use and land cover map for Denmark (Scientific Advisory Report Nr. 351-400; 
### Technical Report from DCE – Danish Centre for Environment and Energy No. 363, s. 69). 
### Aarhus University, Department of Environmental Science. https://dce.au.dk/udgivelser/tr/nr-351-400
# ============================================================

# === load packages ===
library(terra)  # raster-håndtering og distance()
library(sf)     # indlæsning af projektområdepolygon

# === stier ===
leaf_sti      <- #INDSÆT STI MED " " PÅ HVER SIDE (der hvor oprindeligt lag er gemt)
lag_sti       <- #INDSÆT STI MED " " PÅ HVER SIDE (mappen hvor alle miljøvariabler samles)
studieomr_sti <- #INDSÆT STI MED " " PÅ HVER SIDE (studieområdet)

# ============================================================
# DEL A: UDLED DIST_LOEVTRAE.TIF
# ============================================================

# === indlæs Copernicus leaf type-lag ===
r_leaf <- rast(leaf_sti)

# === diagnostik: bekræft struktur inden distance-beregning ===
# Forventet: klasse 1 = løvtræ, klasse 2 = nåletræ, klasse 3 = blandet
# Verificer at klasse 1 faktisk findes og udgør en rimelig andel
print(r_leaf)        # opløsning, extent, CRS = ETRS89 / UTM zone 32N (EPSG:25832) 
cats(r_leaf)         # RAT hvis tilstede
freq(r_leaf)         # frekvens pr. pixelværdi — vis klassefordeling

# === fjern RAT så vi arbejder med rå pixelværdier ===
# freq() ovenfor viste labels ("Copernicus leaf type") i stedet for tal,
# hvilket indikerer at en RAT-kategori er aktiv. Strip RAT for at
# garantere at ifel(r_leaf == 92000100) sammenligner mod pixelværdier.
levels(r_leaf) <- NULL
freq(r_leaf)  # bør nu vise 0, 92000100, 92000200 som tal

# === byg binær maske: klasse 92000100 (løvtræ) = 1, alt andet = NA ===
# Copernicus Dominant Leaf Type bruger 8-cifrede koder, ikke 1/2/3 —
# se RAT-kolonne "Value" fra cats(r_leaf)
loev_target <- ifel(r_leaf == 92000100, 1, NA)

# === beregn afstandslag med terra::distance() ===
dist_loev <- distance(loev_target)

# === verificer raw output ===
res(dist_loev)       # bør være 10 10
crs(dist_loev, describe = TRUE)$code  # bør være 25832
summary(dist_loev)   # min bør være 0 (pixels med løvtræ)

# === gem rå version ===
# Clipping og resampling sker centralt i trin1e hvis laget besluttes inkluderet
writeRaster(dist_loev,
            filename  = paste0(lag_sti, "dist_loevtrae.tif"),
            overwrite = TRUE,
            datatype  = "FLT4S")  # kontinuert afstandslag

# ============================================================
# DEL B: KORRELATIONSDIAGNOSTIK
# ============================================================
## Denne del af processen stammer fra at afstand til løv blev inkluderet sent.
## Det vil sige at der allerede var kørt mange modeller igennem med de andre
## miljøvariabler. De øvrige variabler havde derfor allerede kørt igennem trin 
## 1f på dette tidspunkt.

# === indlæs projektområde (samme som trin1f) ===
studieomr       <- st_read(studieomr_sti, layer = "vp3endelig2022vp3e2022_hovedoplande")
studieomr_union <- st_union(studieomr)
studieomr_vect  <- vect(studieomr_union)

# === indlæs alle eksisterende lag (allerede clippet i trin1f) ===
r_anvend   <- rast(paste0(lag_sti, "arealanvendelse.tif"))
r_infra    <- rast(paste0(lag_sti, "dist_infrastruktur.tif"))
r_vandloeb <- rast(paste0(lag_sti, "dist_vandloeb.tif"))
r_soer     <- rast(paste0(lag_sti, "dist_soer.tif"))
r_slope    <- rast(paste0(lag_sti, "slope_10m.tif"))
r_dhm      <- rast(paste0(lag_sti, "dhm_10m.tif"))

# === clip og align dist_loevtrae til reference-grid (midlertidigt) ===
# Brug r_anvend som reference for korrelationstesten — løvtræ er endnu
# ikke kørt gennem trin1f så vi tilpasser det her uden at gemme den
# tilpassede version (rå version er allerede gemt ovenfor).
r_loev_clip <- crop(dist_loev, studieomr_vect)
r_loev_clip <- mask(r_loev_clip, studieomr_vect)
r_loev_clip <- resample(r_loev_clip, r_anvend, method = "bilinear")

# === byg stak af kontinuerte lag til Pearson-korrelation ===
# Arealanvendelse (kategorisk) håndteres separat længere nede 
stak_kont <- c(r_loev_clip, r_infra, r_vandloeb, r_soer, r_slope, r_dhm)
names(stak_kont) <- c("dist_loevtrae", "dist_infrastruktur",
                      "dist_vandloeb", "dist_soer", "slope", "dhm")

# === sample 10.000 baggrundspunkter ===
# Samme stikprøvestørrelse som modelkørslen i trin3_ENMeval —
# sikrer at korrelationen reflekterer den fordeling modellen ser
set.seed(42)  # reproducerbarhed
bg_vals <- spatSample(stak_kont,
                      size    = 10000,
                      method  = "random",
                      na.rm   = TRUE,
                      as.df   = TRUE)

cat("\nAntal komplette baggrundspunkter:", nrow(bg_vals), "\n")

# === beregn Pearson-korrelationsmatrix ===
cor_mat <- cor(bg_vals, use = "complete.obs", method = "pearson")
cat("\n=== Pearson-korrelationsmatrix (kontinuerte lag) ===\n")
print(round(cor_mat, 3))

# === fremhæv kun korrelationer involverende dist_loevtrae ===
# Dette er nøgleinformationen for beslutningen om inkludering
cat("\n=== Korrelationer mellem dist_loevtrae og øvrige lag ===\n")
loev_korr <- cor_mat["dist_loevtrae", -1]
print(round(sort(abs(loev_korr), decreasing = TRUE), 3))

# Tommelfingerregel for SDM (Dormann et al. 2013):
# |r| < 0.7 = uafhængigt nok til at begge lag kan beholdes
# |r| >= 0.7 = redundans, ét lag bør udelades
cat("\nTommelfingerregel (Dormann et al. 2013): |r| < 0.7 = acceptabel\n")

# === dist_loevtrae fordelt på arealanvendelsesklasser ===
# Sampler samtidigt fra arealanvendelse og dist_loevtrae for at
# se om afstandslaget tilføjer information ud over kategorisk skov/ikke-skov
stak_blandet <- c(r_loev_clip, r_anvend)
names(stak_blandet) <- c("dist_loevtrae", "arealanvendelse")

set.seed(42)
bg_blandet <- spatSample(stak_blandet,
                         size    = 10000,
                         method  = "random",
                         na.rm   = TRUE,
                         as.df   = TRUE)

cat("\n=== Middel dist_loevtrae (meter) pr. arealanvendelsesklasse ===\n")
agg_tab <- aggregate(dist_loevtrae ~ arealanvendelse, data = bg_blandet,
                     FUN = function(x) c(middel = mean(x),
                                         median = median(x),
                                         n      = length(x)))
print(agg_tab)

# === fortolkningsguide ===
# Lav middelafstand i skovklasser (4 = tør skov, 5 = våd skov) =
# forventeligt og bekræfter at laget virker som det skal.
#
# Spørgsmålet er hvor stor variation der er INDEN FOR ikke-skovklasser:
#   - Hvis alle ikke-skovklasser har lignende middelafstand til løv,
#     tilføjer laget primært "indenfor skov vs. udenfor skov"-info
#     som arealanvendelse allerede har.
#   - Hvis der er stor variation, tilføjer laget en gradient
#     der ikke fanges af arealanvendelse — det er argumentet
#     for inkludering.