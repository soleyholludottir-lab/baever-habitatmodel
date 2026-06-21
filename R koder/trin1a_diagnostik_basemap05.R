# ============================================================
# TRIN 1a: DIAGNOSTICER BASEMAP05-STRUKTUR
# Formål: Forstå pixelværdier og RAT-struktur i 
#         lu_aggregated_2024.tif inden lagudledning
# ============================================================

# === load packages ===
library(terra)  # raster-håndtering

# === stier ===
basemap_sti <- #INDSÆT STI MED " " PÅ HVER SIDE

# === indlæs Basemap05 ===
r <- rast(basemap_sti)

# --- grundlæggende info ---
print(r)           # opløsning, extent, CRS, antal lag
crs(r)             # verificer CRS
res(r)             # pixelstørrelse

# --- RAT-struktur (Raster Attribute Table) ---
cats(r)            # returnerer RAT som data.frame


