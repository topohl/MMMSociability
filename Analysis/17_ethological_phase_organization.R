# ================================================================
# Ethological Active/Inactive Phase Organization
# MMMSociability
# ================================================================
# Goal:
#   Treat active/inactive structure as a biologically meaningful axis of
#   home-cage adaptation after social instability.
#
# Language guardrail:
#   Use active/inactive, day/night behavioral structure, dark-phase
#   adaptation, and light-phase rest-like organization. Do not claim
#   circadian or sleep disruption without independent validation.
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(purrr)
  library(tibble)
  library(stringr)
})

source_candidates <- c(
  file.path("Functions", "behavioral_dynamics_helpers.R"),
  file.path("..", "Functions", "behavioral_dynamics_helpers.R"),
  "C:/Users/topohl/Documents/GitHub/MMMSociability/Functions/behavioral_dynamics_helpers.R"
)
duration_helper_candidates <- c(
  file.path("Functions", "duration_normalization_helpers.R"),
  file.path("..", "Functions", "duration_normalization_helpers.R"),
  "C:/Users/topohl/Documents/GitHub/MMMSociability/Functions/duration_normalization_helpers.R"
)
helper_path <- source_candidates[file.exists(source_candidates)][1]
duration_helper_path <- duration_helper_candidates[file.exists(duration_helper_candidates)][1]
if (is.na(helper_path)) stop("Could not find behavioral_dynamics_helpers.R", call. = FALSE)
if (is.na(duration_helper_path)) stop("Could not find duration_normalization_helpers.R", call. = FALSE)
source(helper_path)
source(duration_helper_path)

bin_level <- "5min_based"
base_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID"
input_file <- file.path(base_dir, "analysis_ready/03_derived_metrics", bin_level, "all_behavior_metrics.csv")
output_dir <- file.path(base_dir, "analysis_ready/17_ethological_phase_organization", bin_level)
proximity_col <- "ProximityFraction"
phase_metrics <- c("Movement", "Entropy", "Proximity")

output_dirs <- analysis_output_dirs(output_dir)
write_output_manifest(
  output_dir,
  script_name = "17_ethological_phase_organization.R",
  analysis_name = "ethological active/inactive phase organization",
  primary_tables = c(
    "tables/phase_contrast_features.csv",
    "tables/phase_timing_features.csv",
    "tables/phase_fragmentation_features.csv",
    "tables/phase_recovery_kinetics.csv",
    "tables/phase_predictability_features.csv"
  ),
  primary_figures = c("figures/publication_panels/Fig17_phase_contrast_overview.svg"),
  notes = c("Interpret as phase-specific behavioral organization, not validated circadian/sleep disruption.")
)

raw_dat <- read_behavior_table(input_file)
if (!proximity_col %in% names(raw_dat)) proximity_col <- "Proximity"
behav <- standardize_behavior_columns(raw_dat, proximity_col = proximity_col) %>%
  mutate(
    Group = factor(as.character(Group), levels = unique(c(mmm_group_levels, sort(unique(as.character(Group)))))),
    PhaseClass = case_when(
      str_detect(str_to_lower(as.character(Phase)), "inactive|light|day") ~ "Inactive",
      str_detect(str_to_lower(as.character(Phase)), "active|dark|night") ~ "Active",
      TRUE ~ as.character(Phase)
    ),
    CageChangeIndex = suppressWarnings(as.integer(str_extract(as.character(CageChange), "\\d+"))),
    CageChangeIndex = if_else(is.finite(CageChangeIndex), CageChangeIndex, dense_rank(as.character(CageChange))),
    BinLevel = bin_level
  ) %>%
  arrange(AnimalNum, CageChangeIndex, PhaseClass, TimeIndex)

bin_size_sec <- infer_bin_size_sec(behav)
epoch_duration_qc <- write_epoch_duration_qc(behav, output_dir, metric_source = "17_ethological_phase_organization", bin_size_sec = bin_size_sec)

long_dat <- behav %>%
  pivot_longer(all_of(phase_metrics), names_to = "Metric", values_to = "Value") %>%
  group_by(AnimalNum, CageChange, PhaseClass, Metric) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  mutate(
    WithinPhaseRank = row_number(),
    PhaseN = n(),
    PhaseFraction = WithinPhaseRank / pmax(PhaseN, 1),
    PhaseHour = (WithinPhaseRank - 1) * bin_size_sec / 3600
  ) %>%
  ungroup()

phase_summary <- long_dat %>%
  group_by(BinLevel, AnimalNum, Group, Sex, CageChange, CageChangeIndex, PhaseClass, Metric) %>%
  summarise(
    n_bins = n(),
    phase_mean = mean(Value, na.rm = TRUE),
    phase_median = median(Value, na.rm = TRUE),
    phase_rmssd = calc_rmssd(Value),
    phase_acf1 = calc_acf1(Value),
    onset_slope = if (sum(is.finite(Value) & PhaseFraction <= 0.25) >= 4) unname(coef(lm(Value[PhaseFraction <= 0.25] ~ PhaseHour[PhaseFraction <= 0.25]))[2]) else NA_real_,
    offset_slope = if (sum(is.finite(Value) & PhaseFraction >= 0.75) >= 4) unname(coef(lm(Value[PhaseFraction >= 0.75] ~ PhaseHour[PhaseFraction >= 0.75]))[2]) else NA_real_,
    activity_center_of_mass = if (sum(is.finite(Value)) > 0 && sum(abs(Value), na.rm = TRUE) > 0) weighted.mean(PhaseFraction, w = abs(Value), na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  join_duration_qc(epoch_duration_qc, by_cols = c("AnimalNum", "Group", "Sex", "CageChange"))

missing_phase_classes <- setdiff(c("Active", "Inactive"), unique(as.character(phase_summary$PhaseClass)))
if (length(missing_phase_classes) > 0) {
  stop(
    "Phase contrast analysis requires both Active and Inactive phases. Missing: ",
    paste(missing_phase_classes, collapse = ", "),
    call. = FALSE
  )
}

phase_contrast_features <- phase_summary %>%
  select(BinLevel, AnimalNum, Group, Sex, CageChange, CageChangeIndex, PhaseClass, Metric, phase_mean, phase_rmssd, phase_acf1) %>%
  pivot_wider(names_from = PhaseClass, values_from = c(phase_mean, phase_rmssd, phase_acf1)) %>%
  mutate(
    active_minus_inactive_mean = phase_mean_Active - phase_mean_Inactive,
    active_inactive_ratio_mean = phase_mean_Active / pmax(abs(phase_mean_Inactive), 1e-9),
    active_minus_inactive_rmssd = phase_rmssd_Active - phase_rmssd_Inactive,
    active_minus_inactive_acf1 = phase_acf1_Active - phase_acf1_Inactive,
    phase_contrast_strength = abs(active_minus_inactive_mean)
  )

phase_timing_features <- phase_summary %>%
  select(BinLevel, AnimalNum, Group, Sex, CageChange, CageChangeIndex, PhaseClass, Metric, onset_slope, offset_slope, activity_center_of_mass, n_bins)

state_bins <- behav %>%
  group_by(AnimalNum) %>%
  mutate(active_threshold = quantile(Movement, 0.60, na.rm = TRUE, names = FALSE)) %>%
  ungroup() %>%
  mutate(above_activity_threshold = Movement >= active_threshold)

phase_fragmentation_features <- state_bins %>%
  group_by(BinLevel, AnimalNum, Group, Sex, CageChange, CageChangeIndex, PhaseClass) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  summarise(
    n_bins = n(),
    active_fraction = mean(above_activity_threshold, na.rm = TRUE),
    phase_switch_rate = if_else(n() >= 3, mean(above_activity_threshold != lag(above_activity_threshold), na.rm = TRUE), NA_real_),
    inactivity_fragmentation = if_else(n() >= 3, mean(!above_activity_threshold & lag(above_activity_threshold), na.rm = TRUE), NA_real_),
    mean_movement_rmssd = calc_rmssd(Movement),
    .groups = "drop"
  ) %>%
  normalize_counts_to_rates(c(), duration_col = "total_observation_duration_hours")

control_phase_profile <- phase_contrast_features %>%
  filter(as.character(Group) == "CON") %>%
  group_by(CageChange, CageChangeIndex, Metric) %>%
  summarise(control_phase_contrast = mean(active_minus_inactive_mean, na.rm = TRUE), .groups = "drop")

phase_recovery_kinetics <- phase_contrast_features %>%
  left_join(control_phase_profile, by = c("CageChange", "CageChangeIndex", "Metric")) %>%
  group_by(BinLevel, AnimalNum, Group, Sex, Metric) %>%
  arrange(CageChangeIndex, .by_group = TRUE) %>%
  mutate(
    deviation_from_control_phase_profile = abs(active_minus_inactive_mean - control_phase_contrast),
    phase_profile_recovery = lag(deviation_from_control_phase_profile) - deviation_from_control_phase_profile
  ) %>%
  summarise(
    mean_phase_contrast = mean(active_minus_inactive_mean, na.rm = TRUE),
    mean_deviation_from_control = mean(deviation_from_control_phase_profile, na.rm = TRUE),
    phase_stabilization_slope = if (sum(is.finite(deviation_from_control_phase_profile)) >= 3) unname(coef(lm(deviation_from_control_phase_profile ~ CageChangeIndex))[2]) else NA_real_,
    mean_phase_profile_recovery = mean(phase_profile_recovery, na.rm = TRUE),
    .groups = "drop"
  )

phase_predictability_features <- long_dat %>%
  group_by(BinLevel, AnimalNum, Group, Sex, CageChange, CageChangeIndex, Metric) %>%
  summarise(
    active_mean = mean(Value[PhaseClass == "Active"], na.rm = TRUE),
    inactive_mean = mean(Value[PhaseClass == "Inactive"], na.rm = TRUE),
    pooled_sd = sd(Value, na.rm = TRUE),
    phase_predictability_index = abs(active_mean - inactive_mean) / pmax(pooled_sd, 1e-9),
    phase_regularization_index = 1 / (1 + calc_rmssd(Value)),
    .groups = "drop"
  )

interpretation_guide <- tibble(
  Output = c("phase_contrast_features", "phase_timing_features", "phase_fragmentation_features", "phase_recovery_kinetics", "phase_predictability_features"),
  BiologicalQuestion = c(
    "Is active-phase behavior distinct from inactive-phase behavior?",
    "Where within the phase is behavior concentrated or changing?",
    "Is active/inactive structure fragmented?",
    "Does phase organization move toward the control profile across cage changes?",
    "How separable and regular are active versus inactive behavioral profiles?"
  ),
  PreferredLanguage = c(
    "phase contrast",
    "phase-specific timing structure",
    "phase-specific fragmentation",
    "phase organization recovery",
    "phase predictability/regularity"
  ),
  AvoidLanguage = c("circadian amplitude", "circadian phase shift", "sleep fragmentation", "circadian recovery", "sleep/circadian disruption")
)

write_table(phase_contrast_features, file.path(output_dir, "tables", "phase_contrast_features.csv"))
write_table(phase_timing_features, file.path(output_dir, "tables", "phase_timing_features.csv"))
write_table(phase_fragmentation_features, file.path(output_dir, "tables", "phase_fragmentation_features.csv"))
write_table(phase_recovery_kinetics, file.path(output_dir, "tables", "phase_recovery_kinetics.csv"))
write_table(phase_predictability_features, file.path(output_dir, "tables", "phase_predictability_features.csv"))
write_table(interpretation_guide, file.path(output_dir, "tables", "phase_organization_interpretation_guide.csv"))
generate_duration_sensitivity_outputs(phase_contrast_features, epoch_duration_qc, output_dir, "phase_contrast_features", by_cols = c("AnimalNum", "Group", "Sex", "CageChange"))

plot_tbl <- phase_contrast_features %>% filter(Metric == "Movement", is.finite(active_minus_inactive_mean))
if (nrow(plot_tbl) > 0) {
  p <- ggplot(plot_tbl, aes(Group, active_minus_inactive_mean, fill = Group, colour = Group)) +
    geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey55") +
    geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.75, linewidth = 0.25) +
    geom_jitter(width = 0.06, size = 0.9, alpha = 0.75) +
    facet_grid(Sex ~ CageChange) +
    scale_fill_manual(values = mmm_group_colors, drop = FALSE) +
    scale_colour_manual(values = mmm_group_colors, drop = FALSE) +
    labs(title = "Ethological active/inactive organization", subtitle = "Movement active-minus-inactive phase contrast", x = NULL, y = "Active - inactive movement") +
    make_publication_theme(base_size = 7) +
    theme(legend.position = "none")
  save_plot_svg_pdf(p, file.path(output_dir, "figures/publication_panels/Fig17_phase_contrast_overview"), width = 150, height = 90)
}

message("Ethological phase organization complete: ", output_dir)
