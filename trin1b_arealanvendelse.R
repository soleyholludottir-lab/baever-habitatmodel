# ============================================================
# TRIN 1b: UDLED AREALANVENDELSE_ALLE.TIF (kategorisk, klasse 100–112)
# Formål: Reklassificer Basemap05 til EU level 1 klasser
#         Alle klasser inkl. 100 og 101 bevares (blev udeladt i forrige kørsler, men beholdes nu)
# Hvad scriptet gør trin for trin:
##  Udtrækker RAT'en fra cats() og bygger en mapping fra lange pixelkoder → C_05 værdier
##  classify() laver et nyt raster hvor alle pixels har fået deres C_05-værdi (100–112)
##  Rasteret konverteres til kategorisk faktor
##  Der laves danske oversættelser til labels
##  Gemt som arealanvendelse_alle.tif
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
lag_sti     <- #INDSÆT STI MED " " PÅ HVER SIDE


# === indlæs Basemap05 (allerede diagnosticeret i trin 1a) ===
r <- rast(basemap_sti)

# === byg reklassificeringsmatrix fra RAT ===
rat <- cats(r)[[1]]
rcl <- matrix(
  c(rat$Value, as.integer(rat$C_05)),
  ncol = 2
)

# === reklassificer rasteret til EU level 1 koder ===
r_c05 <- classify(r, rcl, others = NA)

# === verificer at reklassificeringen lykkedes ===
freq(r_c05)  # bør vise værdierne 100, 101, 102 ... 112

# === sæt som kategorisk lag ===
r_anvend_alle <- as.factor(r_c05)

# === verificer ===
is.factor(r_anvend_alle)   # skal returnere TRUE
cats(r_anvend_alle)        # bør vise klasse 100–112 med labels
crs(r_anvend_alle)         # skal være EPSG:25832
res(r_anvend_alle)         # skal være 10 10

# === tilføj danske labels ===
# denne oversættelse har forfatter selv lavet ud fra de engelske navne i den 
# tekniske anvisning for Basemap05 (Levin, G. 2026)
levels(r_anvend_alle) <- data.frame(
  ID    = c(100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112),
  label = c("Ikke klassificeret",
            "Bebyggelse og infrastruktur",
            "Landbrugsarealer",
            "Græsland",
            "Skov",
            "Heder og krat",
            "Sparsomt bevoksede arealer",
            "Ferske vådområder",
            "Vandløb",
            "Søer",
            "Strandeng og marine indløb",
            "Kyststrande og klitter",
            "Hav")
)

# === verificer labels ===
cats(r_anvend_alle)

# === gem med labels ===
writeRaster(r_anvend_alle,
            filename  = paste0(lag_sti, "arealanvendelse_alle.tif"),
            overwrite = TRUE,
            datatype  = "INT2U")