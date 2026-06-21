# ============================================================
# TRIN 3: DEN ENDELIGE MODEL MED 5 MILJØVARIABLER
#
# KØR MAXENT MODEL MED REDUCERET VARIABELSÆT TIL SAMMENLIGNING
# Formål: Empirisk test af om en model med kun de bidragende
#         variabler performer bedre end den fulde 7-variabel model.
#         Resultaterne sammenlignes med trin3_ENMeval.R (7 var).
#         Denne version inkluderer alle klasser også klasse 100 og 101 i
#         arealanvendelse (bebyggelse og ikke klassificeret).
#
# Variabler inkluderet (5):
#   - arealanvendelse_alle  (alle klasser 100–112, inkl. 100 og 101)
#   - dhm                   (AUC-fald 0.0137 — sekundær)
#   - dist_vandloeb         (AUC-fald 0.0025)
#   - dist_soer             (AUC-fald 0.0021)
#   - dist_loevtrae         (AUC-fald 0.00003 — bevares for empirisk
#                           test af om bidraget stiger uden variablerne hældning 
#                           og afstand til infrastruktur)
#
# Variable udeladt fra denne kørsel (begrundelse: støj):
#   - hældning             (AUC-fald -0.0005, negativt bidrag)
#   - dist_infrastruktur   (AUC-fald -0.0012, negativt bidrag)
#   Begge variabler forværrede modelydelse i den fulde jackknife-
#   analyse. — fjernelse er i overenstemmelse med L1-regulariseringens variabel-
#   udvælgelsesprincip (Elith et al. 2011).
#
# Litteratur:
#   - ENMeval 2.0: Kass et al. (2021)
#   - SWD-format: Elith et al. (2010)
#   - AICc modeludvælgelse: Muscarella et al. (2014)
#   - Block partitioning: Elith et al. (2010), Treves et al. (2022)
#   - Feature classes og rm 0.5–5: Swinnen et al. (2017)
#   - Jackknife variabelbidrag: Phillips (2006), Elith et al. (2011)
#   - L1-regularisering som implicit variabel-udvælgelse:
#     Elith et al. (2011)
# ============================================================

# === load packages ===
library(terra)      # raster-håndtering
library(sf)         # vektor-håndtering
library(ENMeval)    # MaxEnt modellering og tuning
library(maxnet)     # maxnet backend til ENMeval
library(ggplot2)    # visualiserer resultater

# === stier ===
lag_sti    <- #INDSÆT STI MED " " PÅ HVER SIDE (mappen hvor alle miljøvariabler samles)
github_sti <- #INDSÆT STI MED " " PÅ HVER SIDE (baever_habitatmodel)

# === indlæs occurrence-koordinater (WGS84) ===
occ_coords <- read.csv(paste0(lag_sti, "occ_coords_wgs84.csv"))

# === byg env_stack (5 variable) ===
# slope og dist_infrastruktur udeladt (se kodens indledning)
# arealanvendelse_alle inkluderer klasse 100 og 101 (bebyggelse
# og ikke klassificeret) — lagnavnet i stack sættes til
# "arealanvendelse" så MaxEnt behandler det som samme variabeltype
env_stack <- c(
  rast(paste0(lag_sti, "arealanvendelse_alle.tif")),  # alle klasser inkl. 100–112
  rast(paste0(lag_sti, "dist_vandloeb.tif")),
  rast(paste0(lag_sti, "dist_soer.tif")),
  rast(paste0(lag_sti, "dist_loevtrae.tif")),
  rast(paste0(lag_sti, "dhm_10m.tif"))
)

# === navngiv lag i stack ===
# arealanvendelse_alle omdøbes til "arealanvendelse" i stack —
# MaxEnt håndterer det som kategorisk variabel via categoricals-
# argumentet i ENMevaluate()
names(env_stack) <- c("arealanvendelse", "dist_vandloeb",
                      "dist_soer", "dist_loevtrae", "dhm")

# === verificer at arealanvendelse er faktor med alle 13 klasser ===
# bør vise klasse 100–112 inkl. 100 og 101
cats(env_stack[[1]])

# === transformer env_stack til WGS84 til ENMeval ===
# ENMeval kræver WGS84 — projected CRS giver intern vect()-fejl
env_stack_wgs84 <- project(env_stack, "EPSG:4326")

# === byg 100m pre-model buffer maske til baggrundspunkter ===
# Biologisk begrundelse: bævere fouragerer primært inden for
# 10–50m fra vandveje (Macdonald et al. 1995; Schley 2004).
# En 100m buffer er konservativ og konsistent med Swinnen et al.
# (2017), som afgrænsede modellen til habitat 50m fra vandveje.
# Masken bygges fra projected CRS (EPSG:25832) for korrekte
# meterdistancer, derefter transformeret til WGS84 inden sampling.
dist_vl_proj <- rast(paste0(lag_sti, "dist_vandloeb.tif"))
dist_so_proj <- rast(paste0(lag_sti, "dist_soer.tif"))

# minimums-afstand til vand (vandløb ELLER sø) — celle beholdes
# hvis den er inden for 100m fra enten vandløb eller sø
min_dist     <- min(dist_vl_proj, dist_so_proj)
buffer_maske <- ifel(min_dist <= 100, 1, NA)

# diagnostik: rapportér andel af studieområde inden for buffer
n_buf <- global(buffer_maske, "sum", na.rm = TRUE)[[1]]
n_tot <- global(!is.na(rast(paste0(lag_sti, "dhm_10m.tif"))),
                "sum", na.rm = TRUE)[[1]]
cat(sprintf("100m buffer maske: %.1f%% af studieomraadet (%d celler)\n",
            n_buf / n_tot * 100, as.integer(n_buf)))

# transformer buffer-maske til WGS84 til brug med env_stack_wgs84
buffer_maske_wgs84 <- project(buffer_maske, "EPSG:4326", method = "near")

# anvend maske på env_stack_wgs84 — baggrundspunkter samples
# kun fra celler inden for 100m fra vandløb eller sø
env_stack_buf <- mask(env_stack_wgs84, buffer_maske_wgs84)

# === sample baggrundspunkter fra 100m buffer-zonen ===
# Samme seed (42) som 7-variabel-modellen — sikrer at sammenligningen
# kører på samme baggrundspunkter, så forskelle skyldes variabelvalg
# og ikke stikprøvevariation
set.seed(42)  # reproducerbarhed
bg_pts    <- spatSample(env_stack_buf, size = 10000,
                        method = "random", na.rm = TRUE, xy = TRUE)
bg_coords <- bg_pts[, c("x", "y")]
bg_env    <- bg_pts[, names(env_stack_buf)]

# === pre-ekstraher miljøværdier ved occurrence-punkter (SWD) ===
# Elith et al. (2010): SWD-format leverer pre-ekstraherede
# miljøværdier som data.frame i stedet for raster —
# modelberegningerne er identiske, men forudsigelsesraster
# laves separat for den bedste model med terra::predict()
occ_env <- extract(env_stack_wgs84,
                   vect(occ_coords, geom = c("longitude", "latitude"),
                        crs = "EPSG:4326"))[, -1]

# === kombiner koordinater og miljøværdier (SWD-format) ===
# Dokumentation: occs og bg skal indeholde koordinater + miljøværdier
# når raster ikke angives (Kass et al. 2021, ENMeval dokumentation)
occs_swd <- cbind(occ_coords, occ_env)
bg_swd   <- cbind(bg_coords,  bg_env)

# === ensret kolonnenavne i bg_swd ===
# spatSample() returnerer "x" og "y" — omdøbes til "longitude" og "latitude"
# så rbind(occs, bg) internt i ENMeval ikke fejler
names(bg_swd)[1:2] <- c("longitude", "latitude")

# === ensret factor levels for arealanvendelse ===
# ENMeval kører rbind(occs, bg) internt — kræver identiske factor levels.
# Med arealanvendelse_alle (13 klasser) dækker baggrundspunkterne
# sandsynligvis alle klasser inkl. 100 og 101, men occurrence-punkter
# dækker ikke nødvendigvis alle — union() sikrer konsistente levels.
# MaxEnt bruger fraværet af occurrence-punkter i visse klasser som
# information om habitatundgåelse — biologisk korrekt.
alle_levels <- union(levels(occs_swd$arealanvendelse),
                     levels(bg_swd$arealanvendelse))

occs_swd$arealanvendelse <- factor(occs_swd$arealanvendelse,
                                   levels = alle_levels)
bg_swd$arealanvendelse   <- factor(bg_swd$arealanvendelse,
                                   levels = alle_levels)

# === verificer SWD ===
nrow(occs_swd)   # skal være 55 (minus evt. NA-punkter)
nrow(bg_swd)     # skal være 10000
head(occs_swd)   # tjek at miljøværdier er korrekt ekstraheret

# === kør ENMeval med SWD-format ===
# 30 modeller: 3 fc-kombinationer x 10 rm-værdier (0.5–5)
# fc: L, LQ, LQH — forsvarligt med 55 punkter
#   L:               altid tilladt
#   LQ:              kræver ≥10 forekomstpunkter
#   LQH:             kræver ≥15 forekomstpunkter
#   T og P udelades: kræver ≥80 forekomstpunkter (Elith et al. 2011)
#
# rm: 0.5–5 i 0.5 intervaller (Swinnen et al. 2017)
# block partitioning: spatial cross-validation
#   (Elith et al. 2010, Muscarella et al. 2014, Treves et al. 2022)
# raster.preds = FALSE: forudsigelsesraster laves separat for
#   den bedste model med terra::predict() — undgår lang køretid
#   for alle 30 modeller (Kass et al. 2021, ENMeval dokumentation)
# envs udelades: SWD-format — miljøværdier er pre-ekstraheret
enmeval_res <- ENMevaluate(
  occs         = occs_swd,
  bg           = bg_swd,
  algorithm    = "maxnet",
  partitions   = "block",
  tune.args    = list(fc = c("L", "LQ", "LQH"),
                      rm = seq(0.5, 5, 0.5)),
  categoricals = "arealanvendelse",
  raster.preds = FALSE,
  parallel     = FALSE
)

# === gem enmeval_res til disk (undgår tab ved R-crash) ===
saveRDS(enmeval_res, file = paste0(github_sti, "v5/enmeval_res_5var_alle.rds"))
cat("enmeval_res_5var_alle gemt til disk.\n")

# === vis modelresultater ===
res <- eval.results(enmeval_res)
res[, c("fc", "rm", "AICc", "delta.AICc", "auc.val.avg", "or.10p.avg")]

# === plot tuning-resultater på tværs af fc/rm-kombinationer ===
# evalplot.stats() er ENMevals indbyggede tuning-figur (Kass et al. 2021).
# Bruger or.10p (omission rate) og cbi.val (Boyce Index) som i vignettens
# sekventielle modelvalgskriterium — viser visuelt hvordan kandidatmodellerne
# rangerer på tværs af fc og rm, og dermed hvorfor den valgte model er optimal.
library(ggplot2)  # ggsave kræver ggplot2 eksplicit indlæst

tuning_plot <- evalplot.stats(
  e          = enmeval_res,
  stats      = c("or.10p", "cbi.val"),
  color      = "fc",
  x.var      = "rm",
  error.bars = TRUE
)

print(tuning_plot)

ggsave(
  filename = paste0(github_sti, "v5/tuning_plot_v5_5var_alle.png"),
  plot     = tuning_plot,
  width    = 9,
  height   = 5,
  dpi      = 300
)
cat("Tuning-plot gemt: tuning_plot_v5_5var_alle.png\n")

# === udtræk bedste model (laveste AICc) ===
# Muscarella et al. (2014): AICc-baseret modeludvælgelse
bedste_model_navn <- res$tune.args[which.min(res$AICc)]
bedste_model      <- eval.models(enmeval_res)[[bedste_model_navn]]
cat("Bedste model:", as.character(bedste_model_navn), "\n")

# ============================================================
# JACKKNIFE VARIABELBIDRAG FOR DEN BEDSTE MODEL
# Metode: ENMevals eval.variable.importance() virker kun med
#         maxent.jar, ikke maxnet (Kass et al. 2021).
#         For maxnet bruges jackknife-test: modellen genoptrænes
#         med én variabel udeladt ad gangen, og faldet i AUC
#         måler variablens bidrag. Dette svarer til jackknife-
#         testen i MaxEnt-softwaren (Phillips 2006; Elith et al.
#         2011) og er den anbefalede metode.
# ============================================================

# hjælpefunktion til AUC-beregning
auc_beregn <- function(pred_occ, pred_bg) {
  u <- sum(outer(pred_occ, pred_bg, ">")) +
    0.5 * sum(outer(pred_occ, pred_bg, "=="))
  u / (length(pred_occ) * length(pred_bg))
}

# fjern koordinatkolonner og NA-rækker til jackknife
occ_env_jk <- occs_swd[complete.cases(occs_swd), -(1:2)]
bg_env_jk  <- bg_swd[, -(1:2)]
var_navne  <- names(occ_env_jk)

# beregn baseline AUC for fuld model
pred_occ_fuld <- as.vector(predict(bedste_model, occ_env_jk, type = "cloglog"))
pred_bg_fuld  <- as.vector(predict(bedste_model, bg_env_jk,  type = "cloglog"))
auc_fuld      <- auc_beregn(pred_occ_fuld, pred_bg_fuld)
cat("Fuld model AUC:", round(auc_fuld, 4), "\n")

# jackknife: genoptræn model uden én variabel ad gangen
cat("\nKører jackknife variabelbidrag...\n")
jackknife_res <- data.frame(variabel = var_navne, auc_uden = NA_real_)

for (i in seq_along(var_navne)) {
  occ_reduced <- occ_env_jk[, -i, drop = FALSE]
  bg_reduced  <- bg_env_jk[,  -i, drop = FALSE]
  
  # ensret factor levels for arealanvendelse hvis den er med
  if ("arealanvendelse" %in% names(occ_reduced)) {
    alle_lev <- union(levels(occ_reduced$arealanvendelse),
                      levels(bg_reduced$arealanvendelse))
    occ_reduced$arealanvendelse <- factor(occ_reduced$arealanvendelse,
                                          levels = alle_lev)
    bg_reduced$arealanvendelse  <- factor(bg_reduced$arealanvendelse,
                                          levels = alle_lev)
  }
  
  # genoptræn model uden variabel i
  p_vec <- c(rep(1, nrow(occ_reduced)), rep(0, nrow(bg_reduced)))
  mod_reduced <- maxnet::maxnet(
    p    = p_vec,
    data = rbind(occ_reduced, bg_reduced)
  )
  
  pred_occ_red <- as.vector(predict(mod_reduced, occ_reduced, type = "cloglog"))
  pred_bg_red  <- as.vector(predict(mod_reduced, bg_reduced,  type = "cloglog"))
  
  jackknife_res$auc_uden[i] <- auc_beregn(pred_occ_red, pred_bg_red)
  cat(" ", var_navne[i], "- AUC uden:", round(jackknife_res$auc_uden[i], 4), "\n")
}

# AUC-fald = fuld model AUC minus AUC uden variabel i
# højere fald = vigtigere variabel
jackknife_res$auc_fald <- auc_fuld - jackknife_res$auc_uden
jackknife_res <- jackknife_res[order(-jackknife_res$auc_fald), ]

cat("\n=== Jackknife variabelbidrag (5 variable — arealanvendelse alle klasser) ===\n")
print(jackknife_res)

# gem jackknife resultater
write.csv(jackknife_res,
          file      = paste0(github_sti, "v5/jackknife_res_5var_alle.csv"),
          row.names = FALSE)

# plot jackknife
barplot(jackknife_res$auc_fald,
        names.arg = jackknife_res$variabel,
        las       = 2,
        main      = "Jackknife variabelbidrag — 5 variable, alle klasser (AUC-fald)",
        ylab      = "AUC-fald ved udeladelse",
        col       = "steelblue",
        cex.names = 0.8)

# === lav forudsigelsesraster for bedste model ===
# terra::predict() bruges til SpatRaster output —
# maxnet:::predict.maxnet() (triple colon) kræves da funktionen
# ikke er eksporteret fra maxnet namespace
habitatkort <- terra::predict(
  env_stack_wgs84,
  bedste_model,
  fun = function(model, ...) maxnet:::predict.maxnet(model, ..., type = "cloglog"),
  na.rm = TRUE
)

# === gem habitatkort til disk ===
# '_alle' suffix skelner fra tidligere kørsel med 7 variabler (originale arealanvendelse)
writeRaster(habitatkort,
            filename  = paste0(lag_sti, "habitatkort_v5_5var_alle.tif"),
            overwrite = TRUE,
            datatype  = "FLT4S")
cat("Habitatkort_5var_alle gemt til disk.\n")

# === verificer habitatkort ===
summary(habitatkort)

