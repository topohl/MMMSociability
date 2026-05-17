# ================================================================
# GAMM trajectory-feature extraction
# MMMSociability
# ================================================================
# Goal:
#   Fit group-blind smooth behavioral trajectories and extract animal-level features:
#   AUC, mean predicted value, peak, trough, dynamic range, time-to-peak,
#   RMSSD of prediction, and ACF1 of prediction.
#
# Input expectation:
#   Run Analysis/03_build_multiscale_behavior_metrics.R first.
#
# Recommended scale:
#   10–30 min bins. Default is 30min_based for smooth trajectory summaries.
# Design:
#   Group/Sex are intentionally not used in the GAMM formula. They are used
#   only after prediction for summaries and downstream descriptive contrasts.
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(mgcv)
  library(pracma)
})

source("C:/Users/topohl/Documents/GitHub/MMMSociability/Functions/behavioral_dynamics_helpers.R")
source("C:/Users/topohl/Documents/GitHub/MMMSociability/Functions/duration_normalization_helpers.R")

# ------------------------------------------------
# USER INPUT
# ------------------------------------------------

bin_level <- "30min_based"
input_file <- file.path("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/analysis_ready/03_derived_metrics", bin_level, "all_behavior_metrics.csv")
output_dir <- file.path("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/analysis_ready/06_behavioral_dynamics/gamm_features", bin_level)

# Use normalized proximity for GAMM trajectories. Raw contact seconds scale with bin size.
proximity_col <- "ProximityFraction"

ensure_dir(output_dir)
ensure_dir(file.path(output_dir, "tables"))
ensure_dir(file.path(output_dir, "figures"))
output_dirs <- analysis_output_dirs(output_dir)
write_output_manifest(
  output_dir,
  script_name = "11_gamm_trajectory_features.R",
  analysis_name = "GAMM trajectory features",
  primary_tables = c(
    "tables/gamm_trajectory_features.csv",
    "tables/gamm_model_qc.csv",
    "tables/gamm_trajectory_robustness_summary.csv",
    "tables/gamm_trajectory_interpretation_guide.csv",
    "tables/epoch_duration_qc.csv"
  ),
  notes = c("Use auc_per_hour rather than raw auc when comparing unequal-duration cage changes.")
)

raw_dat <- read_behavior_table(input_file)
if (!proximity_col %in% names(raw_dat)) proximity_col <- "Proximity"
behav <- standardize_behavior_columns(raw_dat, proximity_col = proximity_col)

epoch_duration_qc <- write_epoch_duration_qc(behav, output_dir, metric_source = "11_gamm_trajectory_features", bin_size_sec = infer_bin_size_sec(behav))

metric_names <- c("Movement", "Entropy", "Proximity")
all_features <- list()
model_qc <- list()

for (metric in metric_names) {
  dat <- behav %>%
    transmute(
      AnimalNum = factor(AnimalNum),
      Group = factor(Group),
      Sex = factor(Sex),
      Phase = factor(Phase),
      CageChange = factor(CageChange),
      TimeIndex,
      Value = .data[[metric]]
    ) %>%
    filter(
      is.finite(Value),
      is.finite(TimeIndex),
      !is.na(AnimalNum),
      !is.na(Group)
    ) %>%
    droplevels()

  if (nrow(dat) < 20 || n_distinct(dat$TimeIndex) < 5 || n_distinct(dat$AnimalNum) < 3) {
    model_qc[[metric]] <- tibble(
      Metric = metric,
      BinLevel = bin_level,
      ProximityInput = proximity_col,
      status = "skipped",
      reason = "Too few rows, time points, or animals",
      n_rows = nrow(dat),
      n_animals = n_distinct(dat$AnimalNum),
      n_timepoints = n_distinct(dat$TimeIndex)
    )
    next
  }

  dat <- dat %>% mutate(TimeScaled = as.numeric(scale(TimeIndex)))
  smooth_k <- min(6, max(3, n_distinct(dat$TimeIndex) - 1))
  animal_smooth_k <- min(4, max(3, n_distinct(dat$TimeIndex) - 1))

  fit <- mgcv::bam(
    Value ~ Phase + CageChange +
      s(TimeScaled, AnimalNum, bs = "fs", k = animal_smooth_k),
    data = dat,
    discrete = TRUE,
    method = "fREML"
  )

  dat$Pred <- predict(fit, newdata = dat)

  feature_tbl <- dat %>%
    group_by(Group, Sex, Phase, CageChange, AnimalNum) %>%
    arrange(TimeIndex, .by_group = TRUE) %>%
    summarise(
      Metric = metric,
      BinLevel = bin_level,
      ProximityInput = proximity_col,
      auc = if_else(n() >= 4, pracma::trapz(TimeIndex, Pred), NA_real_),
      mean_pred = mean(Pred, na.rm = TRUE),
      peak = max(Pred, na.rm = TRUE),
      trough = min(Pred, na.rm = TRUE),
      dynamic_range = peak - trough,
      time_to_peak = TimeIndex[which.max(Pred)][1],
      rmssd_pred = calc_rmssd(Pred),
      acf1_pred = calc_acf1(Pred),
      n_bins = n(),
      .groups = "drop"
    ) %>%
    join_duration_qc(epoch_duration_qc) %>%
    normalize_auc_per_hour("auc")

  all_features[[metric]] <- feature_tbl
  write_table(feature_tbl, file.path(output_dir, "tables", paste0(safe_name(metric), "_gamm_features.csv")))

  model_qc[[metric]] <- tibble(
    Metric = metric,
    BinLevel = bin_level,
    ProximityInput = proximity_col,
    status = "fit",
    reason = NA_character_,
    group_blind = TRUE,
    model_formula = paste(
      "Value ~ Phase + CageChange +",
      "s(TimeScaled, AnimalNum, bs='fs')"
    ),
    n_rows = nrow(dat),
    n_animals = n_distinct(dat$AnimalNum),
    n_timepoints = n_distinct(dat$TimeIndex),
    smooth_k = smooth_k,
    animal_smooth_k = animal_smooth_k,
    edf_total = sum(summary(fit)$s.table[, "edf"], na.rm = TRUE),
    deviance_explained = summary(fit)$dev.expl
  )
}

combined_features <- bind_rows(all_features)
write_table(combined_features, file.path(output_dir, "tables", "combined_gamm_features.csv"))
write_table(bind_rows(model_qc), file.path(output_dir, "tables", "gamm_model_qc.csv"))

trajectory_robustness_summary <- combined_features %>%
  group_by(Metric, BinLevel, ProximityInput, Sex, Phase, CageChange) %>%
  summarise(
    full_data_effect = mean(auc_per_hour, na.rm = TRUE),
    excluding_short_duration_effect = mean(auc_per_hour[!short_epoch %in% TRUE & (is.na(cage_change_duration_class) | cage_change_duration_class != "short")], na.rm = TRUE),
    delta_estimate = excluding_short_duration_effect - full_data_effect,
    delta_cohen_d = NA_real_,
    direction_stable = sign(full_data_effect) == sign(excluding_short_duration_effect) | abs(full_data_effect) < .Machine$double.eps,
    StableForMainText = direction_stable %in% TRUE & (is.na(delta_cohen_d) | abs(delta_cohen_d) < 0.30),
    .groups = "drop"
  )

trajectory_interpretation_guide <- tibble(
  Feature = c("auc_per_hour", "dynamic_range", "time_to_peak", "rmssd_pred", "acf1_pred"),
  BiologicalInterpretation = c(
    "Duration-normalized integrated trajectory magnitude.",
    "Amplitude of longitudinal adaptation trajectory.",
    "Timing of peak response after regrouping.",
    "Smoothed trajectory instability.",
    "Smoothed trajectory inertia/persistence."
  ),
  ReviewerRisk = c("medium", "low", "low", "medium", "medium"),
  PreferredLanguage = c("trajectory burden per hour", "trajectory dynamic range", "peak timing", "trajectory instability", "trajectory persistence")
)

write_table(trajectory_robustness_summary, file.path(output_dir, "tables", "gamm_trajectory_robustness_summary.csv"))
write_table(trajectory_interpretation_guide, file.path(output_dir, "tables", "gamm_trajectory_interpretation_guide.csv"))

message("GAMM trajectory-feature extraction complete.")
