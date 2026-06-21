# ============================================================
# TRIN 1d: BYG VANDLAGSRASTER (dist_vandloeb.tif og dist_soer.tif)
# Formål: Beregn kontinuerte afstandslag til vandløb og søer
#         fra Datafordeler geometrier, clippet til projektområdet
# Metode: st_intersection() → terra::rasterize() → terra::distance()
# Reference-grid: 10x10m, EPSG:25832 (fra lu_aggregated_2024.tif (Levin, G. 2026))
#
##
## Datakilde: GeoDanmark. (2023). Danmarks geografi – GeoDanmark grunddata 
##[Geodatasæt; lag: Vandløbsmidte, Sø]. Klimadatastyrelsen. Hentet 15. marts 2026, 
##fra Datafordeler https://dataforsyningen.dk/data/3563 
# ============================================================

# === load packages ===
library(terra)  # raster-håndtering
library(sf)     # vektor-håndtering

# === stier ===
basemap_sti   <- #INDSÆT STI MED " " PÅ HVER SIDE
lag_sti       <- #INDSÆT STI MED " " PÅ HVER SIDE (mappen hvor alle miljøvariabler samles)
vandloeb_sti  <- #INDSÆT STI MED " " PÅ HVER SIDE (der hvor original-laget er gemt)
soer_sti      <- #INDSÆT STI MED " " PÅ HVER SIDE (der hvor original-laget er gemt)
studieomr_sti <- #INDSÆT STI MED " " PÅ HVER SIDE (der hvor studieområdets afgrænsning er gemt)

# === indlæs reference-raster og projektområde ===
r_ref     <- rast(basemap_sti)
studieomr <- st_read(studieomr_sti, layer = "vp3endelig2022vp3e2022_hovedoplande")

# === indlæs vandgeometrier ===
vandloeb <- st_read(vandloeb_sti, layer = "watercourselink")
soer     <- st_read(soer_sti,     layer = "soe")

# === verificer geometritype ===
st_geometry_type(vandloeb) |> unique()
st_geometry_type(soer)     |> unique()

# === slå geometrier sammen inden clip ===
# st_union() reducerer kompleksiteten fra ~936.000 segmenter til én geometri
# før intersection — markant hurtigere end segment-for-segment
studieomr_union <- st_union(studieomr)
vandloeb_union  <- st_union(vandloeb)
soer_union      <- st_union(soer)

# === clip til projektområdet ===
# st_intersection() clipper vandgeometrierne til studieområdets ydre grænse
# st_union() på studieomr sikrer at intersection sker mod én samlet geometri
vandloeb_clip <- st_intersection(vandloeb_union, studieomr_union)
soer_clip     <- st_intersection(soer_union,     studieomr_union)

# === verificer clip ===
st_is_empty(vandloeb_clip)
st_is_empty(soer_clip)
st_geometry_type(vandloeb_clip)
st_geometry_type(soer_clip)

# === konverter til SpatVector til brug med terra ===
vandloeb_vect <- vect(vandloeb_clip)
soer_vect     <- vect(soer_clip)

# === crop reference-raster til studieområdet ===
r_ref_crop <- crop(r_ref, studieomr_union)
r_ref_crop <- mask(r_ref_crop, vect(studieomr_union))

# === rasteriser vandgeometrier med croppet Basemap05 ===
# Sikrer at output-rasteret har præcis samme extent, opløsning og CRS
# som referencegitteret (10x10m, EPSG:25832)
vandloeb_r <- rasterize(vandloeb_vect, r_ref_crop, field = 1, background = NA)
soer_r     <- rasterize(soer_vect,     r_ref_crop, field = 1, background = NA)

# === beregn kontinuerte afstandslag ===
# terra::distance() beregner afstand i meter fra nærmeste vandpixel
dist_vandloeb <- distance(vandloeb_r)
dist_soer     <- distance(soer_r)

# === verificer ===
summary(dist_vandloeb)  # min bør være 0
summary(dist_soer)      # min bør være 0

# === gem som individuelle .tif-filer ===
writeRaster(dist_vandloeb,
            filename  = paste0(lag_sti, "dist_vandloeb.tif"),
            overwrite = TRUE,
            datatype  = "FLT4S")  # kontinuert afstandslag

writeRaster(dist_soer,
            filename  = paste0(lag_sti, "dist_soer.tif"),
            overwrite = TRUE,
            datatype  = "FLT4S")  # kontinuert afstandslag
