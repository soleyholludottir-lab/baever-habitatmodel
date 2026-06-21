# ============================================================
# TRIN 1c: UDLED DIST_INFRASTRUKTUR.TIF (kontinuert afstandslag)
# Formål: Beregn afstand til klasse 100+101 (Terrestrial not 
#         classified + Settlements/infrastructure) fra Basemap05
# Metode: Binær maske → terra::distance()
# Litteratur: Afstandslag til infrastruktur som forstyrrelsesvariabel
#             brugt i (Swinnen et al. 2017,  Treves et al. 2022)
#
# OBS: dette trin blev udeladt, da det ikke bidrog til modellens performance
#      og bæverhabitater kan findes i både klasse 100+101 
#
# Hvad scriptet gør:
## 1. Reklassificerer Basemap05 til C_05-værdier (samme som Trin 1b)
## 2. Laver en binær maske hvor klasse 100+101 = 1, resten = 0
## 3. Gemmer masken som infra_maske.tif
## 4. Beregner afstand til nærmeste infra-pixel med terra::distance()
## 5. Verificerer og gemmer det færdige lag
###
### Basemap05 kilde:
### Levin, G. (2026). Basemap05. Documentation of the data and method for elaboration 
### of a land use and land cover map for Denmark (Scientific Advisory Report Nr. 351-400; 
### Technical Report from DCE – Danish Centre for Environment and Energy No. 363, s. 69). 
### Aarhus University, Department of Environmental Science. https://dce.au.dk/udgivelser/tr/nr-351-400
# ============================================================

# === load packages ===
library(terra)  # raster-håndtering

# === stier ===
basemap_sti <- #INDSÆT STI MED " " PÅ HVER SIDE
lag_sti     <- #INDSÆT STI MED " " PÅ HVER SIDE (mappen hvor alle miljøvariabler samles)

# === indlæs Basemap05 og reklassificer til C_05 ===
# (samme reklassificering som trin 1b)
r     <- rast(basemap_sti)
rat   <- cats(r)[[1]]
rcl   <- matrix(c(rat$Value, as.integer(rat$C_05)), ncol = 2)
r_c05 <- classify(r, rcl, others = NA)

# Klasse 100 (Terrestrial not classified) udelades —
# er støjklasse uden biologisk betydning og forringede modellen
# ved at tilføje multikollinearitet (jackknife AUC-fald < 0)
## OBS: dette viste sig at være en fejltagelse
infra_maske <- ifel(r_c05 == 101, 1, 0)

# === verificer masken ===
freq(infra_maske)  # bør vise to værdier: 0 og 1

# === gem maske som individuel .tif ===
writeRaster(infra_maske,
            filename  = paste0(lag_sti, "infra_maske.tif"),
            overwrite = TRUE,
            datatype  = "INT1U")  # binær: kun 0 og 1

# === beregn afstandslag med terra::distance() ===
# Sætter infra-pixels (= 1) som mål, beregner afstand i meter
infra_target <- ifel(infra_maske == 1, 1, NA)
dist_infra   <- distance(infra_target)

# === verificer ===
res(dist_infra)      # skal være 10 10
crs(dist_infra)      # skal være EPSG:25832
summary(dist_infra)  # min bør være 0

# === gem ===
writeRaster(dist_infra,
            filename  = paste0(lag_sti, "dist_infrastruktur.tif"),
            overwrite = TRUE,
            datatype  = "FLT4S")  # kontinuert afstandslag