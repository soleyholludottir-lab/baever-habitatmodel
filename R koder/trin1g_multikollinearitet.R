# ============================================================
# TRIN 1g: MULTIKOLLINEARITETSTJEK PÅ KONTINUERTE MILJØVARIABLER
# Formål: Verificere at de fire kontinuerte prædiktorer ikke
#         er stærkt indbyrdes korrelerede.
# Metode: Pearson-korrelation beregnet på 10.000 baggrundspunkter
#         (samme set.seed(42) og spatSample-procedure som trin3
#         — sikrer konsistens mellem korrelationstjek og modelfit)
# Tærskel: |r| < 0,7 (Dormann et al. 2013)
# ============================================================

# === load packages ===
library(terra)

# === stier ===
lag_sti <- #INDSÆT STI MED " " PÅ HVER SIDE (mappen hvor alle miljøvariabler samles)

# === indlæs de fire kontinuerte prædiktorlag ===
env_kont <- c(
  rast(paste0(lag_sti, "dist_vandloeb.tif")),
  rast(paste0(lag_sti, "dist_soer.tif")),
  rast(paste0(lag_sti, "dist_loevtrae.tif")),
  rast(paste0(lag_sti, "dhm_10m.tif"))
)
names(env_kont) <- c("dist_vandloeb", "dist_soer",
                     "dist_loevtrae", "dhm")

# === sample 10.000 baggrundspunkter ===
set.seed(42)
bg_pts <- spatSample(env_kont, size = 10000,
                     method = "random", na.rm = TRUE)

# === Pearson-korrelationsmatrix ===
cor_mat <- cor(bg_pts, method = "pearson")
cat("\n=== Pearson-korrelationsmatrix ===\n")
print(round(cor_mat, 3))

# === maksimal absolut korrelation (uden for diagonalen) ===
diag(cor_mat) <- NA
max_abs_r <- max(abs(cor_mat), na.rm = TRUE)

cat("\nMaksimal absolut Pearson r:", round(max_abs_r, 3), "\n")
cat("Tærskel (Dormann et al. 2013): 0,7\n")
cat("Resultat:",
    ifelse(max_abs_r < 0.7,
           "OK — under tærskel",
           "ADVARSEL — over tærskel"), "\n")

# === gem til disk ===
write.csv(round(cor_mat, 3),
          file = paste0(lag_sti, "korrelationsmatrix.csv"))
cat("\nKorrelationsmatrix gemt til:",
    paste0(lag_sti, "korrelationsmatrix.csv"), "\n")

### RESULTATER
### === Pearson-korrelationsmatrix ===
###  > print(round(cor_mat, 3))
###
###                dist_vandloeb dist_soer  dist_loevtrae  dhm
### dist_vandloeb         1.000     0.207         0.135  0.045
### dist_soer             0.207     1.000         0.090 -0.151
### dist_loevtrae         0.135     0.090         1.000 -0.199
### dhm                   0.045    -0.151        -0.199  1.000
