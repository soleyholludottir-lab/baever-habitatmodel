# ============================================================
# TRIN 4: NULL MODEL VALIDERING — 5 VARIABLER, ALLE KLASSER
# Castor fiber, Vanddistrikt Sjælland, Danmark (2026)
#
# Formål: Kvantificere statistisk signifikans
# for den valgte models performancemetrikker (AUC og CBI)
# ved at sammenligne med null-modeller bygget på tilfældige
# forekomstdata med samme rumlige struktur.
#
# Denne version svarer til trin3_ENMeval_alle.R —
# Resultaterne overskriver IKKE tidligere modelkørsel.
#
# Metodisk begrundelse: Kass et al. (2021) argumenterer for at
# AUC alene ikke er tilstrækkeligt som performancemetrik,
# Null-modelvalidering udføres som beskrevet af Bohl et al., (2019), implementeret 
# i ENMnulls() (Kass et al., 2021). 100 null-modeller genereres på fiktive forekomster, 
# men med samme baggrundspunkter (her kommer stikprøve-værdien i spil igen) og med 
# samme tuning-indstillinger som den rigtige model. Modellens AUC og CBI sammenlignes 
# derefter med null-fordelingens. Scorer den højere på præstationsparametrene, 
# tyder det på at dens forudsigelser har fanget en reel økologisk sammenhæng 
# (Bohl et al., 2019). 
#
# ============================================================

library(ENMeval)
library(ggplot2)

# === stier ===
github_sti <- #INDSÆT STI MED " " PÅ HVER SIDE (baever_habitatmodel)
lag_sti    <- #INDSÆT STI MED " " PÅ HVER SIDE (mappen hvor alle miljøvariabler samles)

# ============================================================
# TRIN 4A: INDLÆS ENMeval-OBJEKT FRA 5VAR_ALLE KØRSEL
# ============================================================

enmeval_res <- readRDS(paste0(github_sti, "v5/enmeval_res_5var_alle.rds"))

cat("ENMeval-objekt indlæst.\n")
cat("Klasse:", class(enmeval_res), "\n")
cat("Antal modeller evalueret:", nrow(eval.results(enmeval_res)), "\n")

# identificer bedste model
res         <- eval.results(enmeval_res)
bedste_navn <- as.character(res$tune.args[which.min(res$AICc)])
cat("Bedste model til null-test:", bedste_navn, "\n")
cat("Validation AUC (empirisk):",
    round(res$auc.val.avg[which.min(res$AICc)], 3), "\n")

# ============================================================
# TRIN 4B: KØR ENMnulls()
# OBS: mod.settings skal matche den bedste model fra
#      trin3_ENMeval_5var_alle.R — tjek bedste_navn ovenfor
#      og opdater fc og rm inden kørsel.
#
# 100 iterationer er standard (Bohl et al. 2019; Kass et al. 2021)
# Forventet køretid: 1-3 timer alt efter computerens kraft
# ============================================================

cat("\nStarter ENMnulls() — 100 iterationer...\n")
cat("ADVARSEL: Dette kan tage 1-3 timer.\n\n")

no_iter <- 100

# OBS: opdater fc og rm til at matche bedste_navn fra trin 4A
null_model <- ENMnulls(
  e            = enmeval_res,
  mod.settings = list(fc = "LQH", rm = 2),  # <-- ret til bedste model
  no.iter      = no_iter,
  eval.stats   = c("auc.val", "cbi.val"),
  parallel     = FALSE
)

cat("ENMnulls() fuldført.\n")

# ============================================================
# TRIN 4C: GEM OG INSPICÉR RESULTATER
# ============================================================

saveRDS(null_model,
        file = paste0(github_sti, "v5/null_model_v5_5var_alle.rds"))
cat("Null model-objekt gemt: null_model_v5_5var_alle.rds\n")

null_resultater <- null.results(null_model)
cat("\nNull model resultater:\n")
print(null_resultater)

write.csv(null_resultater,
          paste0(github_sti, "v5/null_model_resultater_v5_5var_alle.csv"),
          row.names = FALSE)

# nøgletal
cat("\n--- FORTOLKNING ---\n")
cat("Empirisk AUC (validation):",
    round(res$auc.val.avg[which.min(res$AICc)], 3), "\n")

if ("auc.val.avg" %in% names(null_resultater)) {
  null_auc_mean <- mean(null_resultater$auc.val.avg, na.rm = TRUE)
  null_auc_sd   <- sd(null_resultater$auc.val.avg,   na.rm = TRUE)
  cat("Null AUC (gennemsnit):", round(null_auc_mean, 3), "\n")
  cat("Null AUC (SD):",         round(null_auc_sd,   3), "\n")
  cat("Effect size (z-score):",
      round((res$auc.val.avg[which.min(res$AICc)] - null_auc_mean) /
              null_auc_sd, 2), "\n")
}

# ============================================================
# TRIN 4D: VISUALISERING
# ============================================================

# --- Plot 1: AUC null-distribution ---
if ("auc.val.avg" %in% names(null_resultater)) {
  
  auc_emp <- res$auc.val.avg[which.min(res$AICc)]
  
  null_auc_plot <- ggplot(null_resultater, aes(x = auc.val.avg)) +
    geom_histogram(binwidth = 0.02, fill = "#abd9e9",
                   color = "white", alpha = 0.9) +
    geom_vline(xintercept = auc_emp,
               color = "#d7191c", linewidth = 1.2) +
    geom_vline(xintercept = mean(null_resultater$auc.val.avg, na.rm = TRUE),
               color = "#4d4d4d", linewidth = 0.8, linetype = "dashed") +
    annotate("text", x = auc_emp + 0.01, y = Inf, vjust = 1.5, hjust = 0,
             label = paste0("Empirisk AUC = ", round(auc_emp, 3)),
             color = "#d7191c", size = 3.5) +
    scale_x_continuous(limits = c(0.4, 1.0), breaks = seq(0.4, 1.0, 0.1)) +
    labs(title    = paste0("Null model — AUC (validation), ", bedste_navn,
                           " (alle klasser)"),
         subtitle = paste0("Rød linje = empirisk model | Stiplet = null-gennemsnit",
                           " | n = ", no_iter, " iterationer\n",
                           "Kass et al. 2021; Bohl et al. 2019"),
         x = "AUC (validation)", y = "Antal null-modeller") +
    theme_bw() +
    theme(plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8,  color = "grey40"))
  
  print(null_auc_plot)
  
  ggsave(filename = paste0(github_sti, "v5/null_model_auc_",
                           bedste_navn, "_v5_5var_alle.png"),
         plot = null_auc_plot, width = 7, height = 5, dpi = 300)
  cat("AUC null model-plot gemt.\n")
  
} else {
  cat("OBS: auc.val.avg ikke fundet — tjek names(null_resultater)\n")
}

# --- Plot 2: CBI null-distribution ---
if ("cbi.val.avg" %in% names(null_resultater)) {
  
  cbi_emp  <- res$cbi.val.avg[which.min(res$AICc)]
  cbi_data <- null_resultater$cbi.val.avg
  cbi_ok   <- !is.nan(cbi_data) & !is.na(cbi_data)
  
  if (sum(cbi_ok) > 10) {
    
    null_cbi_plot <- ggplot(data.frame(cbi = cbi_data[cbi_ok]), aes(x = cbi)) +
      geom_histogram(binwidth = 0.05, fill = "#abd9e9",
                     color = "white", alpha = 0.9) +
      geom_vline(xintercept = cbi_emp,
                 color = "#d7191c", linewidth = 1.2) +
      geom_vline(xintercept = 0,
                 linetype = "dashed", color = "grey50", linewidth = 0.5) +
      annotate("text", x = cbi_emp + 0.03, y = Inf, vjust = 1.5, hjust = 0,
               label = paste0("Empirisk CBI = ", round(cbi_emp, 3)),
               color = "#d7191c", size = 3.5) +
      scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25)) +
      labs(title    = paste0("Null model — Boyce Index (validation), ",
                             bedste_navn, " (alle klasser)"),
           subtitle = paste0("Rød linje = empirisk model | Stiplet = CBI = 0",
                             " | n = ", sum(cbi_ok), " gyldige iterationer"),
           x = "Continuous Boyce Index", y = "Antal null-modeller") +
      theme_bw() +
      theme(plot.title    = element_text(size = 11, face = "bold"),
            plot.subtitle = element_text(size = 8,  color = "grey40"))
    
    print(null_cbi_plot)
    
    ggsave(filename = paste0(github_sti, "v5/null_model_cbi_",
                             bedste_navn, "_v5_5var_alle.png"),
           plot = null_cbi_plot, width = 7, height = 5, dpi = 300)
    cat("CBI null model-plot gemt.\n")
    
  } else {
    cat("OBS: CBI returnerer NaN for de fleste null-iterationer — plot springes over.\n")
  }
} else {
  cat("OBS: cbi.val.avg ikke fundet — tjek names(null_resultater)\n")
}

# ============================================================
# TRIN 4E: NØGLETAL TIL METODESEKTION
# ============================================================

cat("\n============================================================\n")
cat("NØGLETAL TIL METODESEKTION\n")
cat("============================================================\n")
cat("Null model: ENMnulls(), ENMeval 2.0.5 (Kass et al. 2021)\n")
cat("Antal iterationer:", no_iter, "\n")
cat("Testet model:", bedste_navn, "\n")
cat("Partitionering (empirisk model): block (4-fold)\n\n")

# p-værdier
if ("p.auc.val.avg" %in% names(null_resultater)) {
  cat("AUC p-værdi (empirisk > null):",
      round(null_resultater$p.auc.val.avg[1], 4), "\n")
} else {
  cat("OBS: p-værdi kolonne ikke fundet — tjek names(null.results(null_model))\n")
}

# ============================================================
# TRIN 4F: EKSPORTÉR rangeModelMetadata-OBJEKT
# Formål: Eksportér automatisk genereret metadata om modellen
#         som CSV-bilag. Komplementerer ODMAP-protokollen.
# Reference: Merow et al. (2019); Kass et al. (2021)
# ============================================================

library(rangeModelMetadata)

# === hent det auto-byggede rmm-objekt ===
rmm <- enmeval_res@rmm

# === authorship ===
rmm$authorship$rmmName   <- "Habitategnethedsmodel for Castor fiber i Vandistrikt Sjælland"
rmm$authorship$names     <- "Sóley Hölludóttir"
rmm$authorship$contact   <- "prg740@alumni.ku.dk"
rmm$authorship$miscNotes <- paste0(
  "Bachelorprojekt 2026, Skovskolen, Institut For Geovidenskab og Naturforvaltning, ",
  "Københavns Universitet."
)

# === studyObjective ===
rmm$studyObjective$purpose <- paste(
  "At producere et habitategnethedskort for bæver (Castor fiber) i" ,
  "Vandområdedistrikt Sjælland ved hjælp af MaxEnt-modellering. MaxEnt algoritmen udregner ",
  "den relative habitategnethed på en kontinuert skala fra 1 til 0 som kortbaseret rasterlag. ",
  "Rasterlaget behandles og visualiseres derefter som et forståeligt habitatkort over hele ",
  "projektområdet og har til formål at vise hvor habitatforholdene fremstår egnede, også udenfor ",
  "bæverens nuværende udbredelsesområde"
)
rmm$studyObjective$rangeType  <- "potential"
rmm$studyObjective$hypotheses <- paste(
  "Habitategnethed bestemmes af afstand til vand, arealanvendelse",
  ", terraen og afstand til løvtræer."
)

# === model ===
rmm$model$selectionRules     <- "AICc-minimering (Muscarella et al. 2014)"
rmm$model$finalModelSettings <- as.character(bedste_navn)
rmm$model$notes              <- paste0(
  "Bedste model valgt blandt 30 kandidater (3 fc x 10 rm) ",
  "ved laveste AICc. Arealanvendelse inkluderer alle klasser 100–112."
)
rmm$model$justification[[1]] <- paste(
  "Block-partitionering valgt ud fra forekomst-datasættets størrelse og ujævne rumlige fordeling",
  "(Elith et al. 2011)"
)

# === assessment ===
rmm$assessment$notes <- paste0(
  "Null-model: ENMnulls med ", no_iter, " iterationer ",
  "(Kass et al. 2021). Empirisk model testet ",
  "for AUC og CBI mod 100 randomiserede null-modeller."
)

# === prediction ===
rmm$prediction$notes <- paste0(
  "Habitategnethed (cloglog, 0-1). Tærskel for binaer klassifikation: 0.740",
  "Gennemsnits-habitategnethed ved forekomst punkterne (Swinnen et al. 2017).",
  "OPDATERING: Ved nærmere vurdering bør tærsklen fastsættes til 10 percentilen i stedet."
)

# === eksportér til CSV ===
rmmToCSV(rmm, paste0(github_sti, "v5/rmm_v5_5var_alle.csv"))
cat("rangeModelMetadata gemt: rmm_v5_5var_alle.csv\n")