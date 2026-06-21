# ============================================================
# TRIN 2: GENNEMGÅ NATURSTYRELSENS FOREKOMSTDATA
# Formål: Indlæs og verificer NST's validerede bæverbo-lokationer
#         og forbered koordinater til ENMeval
#
# OBS: NST forekomstdata er validerede feltdata med høj
#      positionel præcision (Lise, NST). Ingen filtrering eller tynding nødvendig.
# ============================================================

# === load packages ===
library(terra)
library(sf)

# === stier ===
occ_sti       <- #INDSÆT STI MED " " PÅ HVER SIDE (forekomstdata)
studieomr_sti <- #INDSÆT STI MED " " PÅ HVER SIDE (studieområdet)
lag_sti       <- #INDSÆT STI MED " " PÅ HVER SIDE (mappen hvor alle miljøvariabler samles)

# === tjek lagstruktur ===
st_layers(occ_sti)

# === indlæs occurrence-data ===
occ <- st_read(occ_sti, layer = "bverbo2026sg")

# === grundlæggende verificering ===
nrow(occ)           # skal være 55
st_crs(occ)$epsg    # skal være 25832
names(occ)          # vis kolonnenavne

# === undersøg nøglekolonner ===
table(occ$Activity)             # aktive vs. inaktive lokaliteter
table(occ$Activity.latest.obs)  # seneste aktivitetsstatus
table(occ$Naturtype)            # naturtype ved hvert punkt
range(occ$Årstal.observeret., na.rm = TRUE)  # observationsår

# === alle 55 punkter inkluderes ===
# Både aktive og inaktive lokaliteter repræsenterer bekræftet
# habitatbrug og inkluderes i modellen.
# Inaktive lokaliteter kan skyldes forstyrrelser eller naturlig
# spredning.

# === transformer til WGS84 til ENMeval ===
occ_wgs84 <- st_transform(occ, crs = 4326)

# === udtræk koordinater til ENMeval ===
occ_coords <- data.frame(
  longitude = st_coordinates(occ_wgs84)[, 1],
  latitude  = st_coordinates(occ_wgs84)[, 2]
)

# === verificer ===
head(occ_coords)
nrow(occ_coords)  # skal være 55

# === afstandsdiagnostik ===
# Udtræk afstand fra hvert punkt til nærmeste vandløb og sø
# Bruges til at verificere buffer-beslutningen (for baggrundspunkter, trin 3)
r_vandloeb <- rast(paste0(lag_sti, "dist_vandloeb.tif"))
r_soer     <- rast(paste0(lag_sti, "dist_soer.tif"))

occ_vect <- vect(occ)

afstand_vandloeb <- extract(r_vandloeb, occ_vect)[, 2]
afstand_soer     <- extract(r_soer,     occ_vect)[, 2]

summary(afstand_vandloeb)  # median bør være tæt på 0–50m
summary(afstand_soer)      # forventet højere end vandløb

# === gem occ_coords som CSV til ENMeval ===
write.csv(occ_coords,
          file      = paste0(lag_sti, "occ_coords_wgs84.csv"),
          row.names = FALSE)