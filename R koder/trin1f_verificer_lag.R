# ============================================================
# TRIN 1e: CLIP OG VERIFICER ALLE MILJØLAG
# Formål: Clip alle miljølag til projektområdet og verificer
#         at alle lag har identisk extent, opløsning og CRS
# OBS: dist_vandloeb og dist_soer er allerede clippet i trin 1d
#      men verificeres her sammen med de øvrige lag.
#      dist_loevtrae er udledt i trin 1e på fuld DK-extent og
#      clippes til projektområdet her.
# ============================================================

# === load packages ===
library(terra)
library(sf)

# === stier ===
lag_sti       <- #INDSÆT STI MED " " PÅ HVER SIDE (mappen hvor alle miljøvariabler samles)
studieomr_sti <- #INDSÆT STI MED " " PÅ HVER SIDE (studieområdet)

# === indlæs projektområde ===
studieomr       <- st_read(studieomr_sti, layer = "vp3endelig2022vp3e2022_hovedoplande")
studieomr_union <- st_union(studieomr)
studieomr_vect  <- vect(studieomr_union)

# === indlæs alle lag ===
r_anvend   <- rast(paste0(lag_sti, "arealanvendelse_alle.tif"))
r_infra    <- rast(paste0(lag_sti, "dist_infrastruktur.tif"))
r_vandloeb <- rast(paste0(lag_sti, "dist_vandloeb.tif"))
r_soer     <- rast(paste0(lag_sti, "dist_soer.tif"))
r_loev     <- rast(paste0(lag_sti, "dist_loevtrae.tif"))
r_slope    <- rast(paste0(lag_sti, "slope_10m.tif"))
r_dhm      <- rast(paste0(lag_sti, "dhm_10m.tif"))

# === clip alle lag til projektområdet ===
# dist_vandloeb og dist_soer er allerede clippet i trin 1d —
# croppes her igen for at sikre identisk extent på tværs af alle lag.
# dist_loevtrae deler 10x10m grid med arealanvendelse (begge fra
# Basemap05_geotiff) og kræver ingen resampling efter crop.
r_anvend   <- crop(r_anvend,   studieomr_vect); r_anvend   <- mask(r_anvend,   studieomr_vect)
r_infra    <- crop(r_infra,    studieomr_vect); r_infra    <- mask(r_infra,    studieomr_vect)
r_vandloeb <- crop(r_vandloeb, studieomr_vect); r_vandloeb <- mask(r_vandloeb, studieomr_vect)
r_soer     <- crop(r_soer,     studieomr_vect); r_soer     <- mask(r_soer,     studieomr_vect)
r_loev     <- crop(r_loev,     studieomr_vect); r_loev     <- mask(r_loev,     studieomr_vect)
r_slope    <- crop(r_slope,    studieomr_vect); r_slope    <- mask(r_slope,    studieomr_vect)
r_dhm      <- crop(r_dhm,      studieomr_vect); r_dhm      <- mask(r_dhm,      studieomr_vect)

# === resample slope og dhm til referencegitter ===
# Opløsningen er marginalt forskellig (10.00017 x 9.999947) —
# resample() tilpasser til præcis 10x10m grid fra de øvrige lag
# method = "bilinear" bruges til kontinuerte lag (glat interpolation)
r_slope <- resample(r_slope, r_anvend, method = "bilinear")
r_dhm   <- resample(r_dhm,   r_anvend, method = "bilinear")

# === gem clippede og tilpassede versioner ===
writeRaster(r_anvend,
            filename  = paste0(lag_sti, "arealanvendelse_alle.tif"),
            overwrite = TRUE,
            datatype  = "INT2U")

writeRaster(r_infra,
            filename  = paste0(lag_sti, "dist_infrastruktur.tif"),
            overwrite = TRUE,
            datatype  = "FLT4S")

writeRaster(r_vandloeb,
            filename  = paste0(lag_sti, "dist_vandloeb.tif"),
            overwrite = TRUE,
            datatype  = "FLT4S")

writeRaster(r_soer,
            filename  = paste0(lag_sti, "dist_soer.tif"),
            overwrite = TRUE,
            datatype  = "FLT4S")

writeRaster(r_loev,
            filename  = paste0(lag_sti, "dist_loevtrae.tif"),
            overwrite = TRUE,
            datatype  = "FLT4S")

writeRaster(r_slope,
            filename  = paste0(lag_sti, "slope_10m.tif"),
            overwrite = TRUE,
            datatype  = "FLT4S")

writeRaster(r_dhm,
            filename  = paste0(lag_sti, "dhm_10m.tif"),
            overwrite = TRUE,
            datatype  = "FLT4S")

# === verificer alle lag ===
# CRS, opløsning og extent skal være identiske på tværs af alle lag
lag_liste <- list(r_anvend, r_infra, r_vandloeb, r_soer, r_loev, r_slope, r_dhm)
lag_navne <- c("arealanvendelse", "dist_infrastruktur", "dist_vandloeb", "dist_soer",
               "dist_loevtrae", "slope", "dhm")

for (i in seq_along(lag_liste)) {
  cat("\n---", lag_navne[i], "---\n")
  cat("CRS:       ", crs(lag_liste[[i]], describe = TRUE)$code, "\n")
  cat("Opløsning: ", res(lag_liste[[i]]), "\n")
  cat("Extent:    ", as.vector(ext(lag_liste[[i]])), "\n")
  cat("is.factor: ", is.factor(lag_liste[[i]]), "\n")
}
