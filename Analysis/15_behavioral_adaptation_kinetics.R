# ================================================================
# Behavioral Adaptation / Recovery Kinetics
# MMMSociability
# ================================================================
# Goal:
#   Quantify interpretable longitudinal stabilization after social
#   perturbation without relying on abstract manifold geometry.
#
# Biological framing:
#   Resilience may reflect faster stabilization or adaptive reorganization
#   after regrouping, not simply more or less movement.
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
output_dir <- file.path(base_dir, "analysis_ready/15_behavioral_adaptation_kinetics", bin_level)
proximity_col <- "ProximityFraction"
primary_metrics <- c("Movement", "Entropy", "Proximity")
early_fraction <- 0.25
late_fraction <- 0.25

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  min(x)
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

output_dirs <- analysis_output_dirs(output_dir)
write_output_manifest(
  output_dir,
  script_name = "15_behavioral_adaptation_kinetics.R",
  analysis_name = "behavioral adaptation and recovery kinetics",
  primary_tables = c(
    "tables/adaptation_kinetics_features.csv",
    "tables/distance_to_control_trajectories.csv",
    "tables/adaptation_kinetics_interpretation_guide.csv"
  ),
  primary_figures = c("figures/publication_panels/Fig15_adaptation_kinetics_overview.svg"),
  notes = c("Focus on interpretable stabilization/recovery metrics; avoid causal recovery language without endpoint validation.")
)

raw_dat <- read_behavior_table(input_file)
if (!proximity_col %in% names(raw_dat)) proximity_col <- "Proximity"
behav <- standardize_behavior_columns(raw_dat, proximity_col = proximity_col) %>%
  mutate(
    Group = factor(as.character(Group), levels = unique(c(mmm_group_levels, sort(unique(as.character(Group)))))),
    PhaseClass = case_when(
      str_detect(str_to_lower(as.character(Phase)), "active|dark|night") ~ "Active",
      str_detect(str_to_lower(as.character(Phase)), "inactive|light|day") ~ "Inactive",
      TRUE ~ as.character(Phase)
    ),
    CageChangeIndex = suppressWarnings(as.integer(str_extract(as.character(CageChange), "\\d+"))),
    CageChangeIndex = if_else(is.finite(CageChangeIndex), CageChangeIndex, dense_rank(as.character(CageChange))),
    BinLevel = bin_level
  ) %>%
  arrange(AnimalNum, CageChangeIndex, PhaseClass, TimeIndex)

epoch_duration_qc <- write_epoch_duration_qc(behav, output_dir, metric_source = "15_behavioral_adaptation_kinetics", bin_size_sec = infer_bin_size_sec(behav))

long_dat <- behav %>%
  pivot_longer(all_of(primary_metrics), names_to = "Metric", values_to = "Value") %>%
  group_by(AnimalNum, CageChange, PhaseClass, Metric) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  mutate(
    WithinEpochRank = row_number(),
    EpochN = n(),
    EpochFraction = WithinEpochRank / pmax(EpochN, 1),
    EpochHour = (WithinEpochRank - 1) * infer_bin_size_sec(behav) / 3600
  ) %>%
  ungroup()

control_profile <- long_dat %>%
  filter(as.character(Group) == "CON") %>%
  group_by(CageChange, PhaseClass, Metric, WithinEpochRank) %>%
  summarise(control_mean = safe_mean(Value), .groups = "drop")

distance_to_control <- long_dat %>%
  left_join(control_profile, by = c("CageChange", "PhaseClass", "Metric", "WithinEpochRank")) %>%
  mutate(abs_distance_to_control = abs(Value - control_mean)) %>%
  group_by(BinLevel, AnimalNum, Group, Sex, CageChange, CageChangeIndex, PhaseClass, Metric) %>%
  summarise(
    mean_distance_to_control = safe_mean(abs_distance_to_control),
    early_distance_to_control = safe_mean(abs_distance_to_control[EpochFraction <= early_fraction]),
    late_distance_to_control = safe_mean(abs_distance_to_control[EpochFraction >= 1 - late_fraction]),
    control_convergence = early_distance_to_control - late_distance_to_control,
    .groups = "drop"
  )

adaptation_features <- long_dat %>%
  group_by(BinLevel, AnimalNum, Group, Sex, CageChange, CageChangeIndex, PhaseClass, Metric) %>%
  summarise(
    n_bins = n(),
    mean_value = safe_mean(Value),
    early_mean = safe_mean(Value[EpochFraction <= early_fraction]),
    late_mean = safe_mean(Value[EpochFraction >= 1 - late_fraction]),
    early_late_shift = late_mean - early_mean,
    rebound_magnitude = safe_max(Value) - early_mean,
    volatility_early = calc_rmssd(Value[EpochFraction <= 0.50]),
    volatility_late = calc_rmssd(Value[EpochFraction > 0.50]),
    volatility_decay = volatility_early - volatility_late,
    recovery_slope = if (sum(is.finite(Value)) >= 6) unname(coef(lm(Value ~ EpochHour))[2]) else NA_real_,
    stabilization_time_bins = {
      baseline <- safe_median(Value[EpochFraction >= 0.75])
      tol <- sd(Value[is.finite(Value)]) * 0.50
      stable_idx <- which(abs(Value - baseline) <= tol)
      if (length(stable_idx) == 0 || !is.finite(tol)) NA_real_ else min(stable_idx)
    },
    adaptation_half_life_bins = {
      start <- early_mean
      target <- late_mean
      half_target <- start + 0.5 * (target - start)
      distance_to_half <- abs(Value - half_target)
      min_distance <- safe_min(distance_to_half)
      idx <- which(distance_to_half == min_distance)[1]
      if (length(idx) == 0 || !is.finite(min_distance)) NA_real_ else idx
    },
    .groups = "drop"
  ) %>%
  left_join(distance_to_control, by = c("BinLevel", "AnimalNum", "Group", "Sex", "CageChange", "CageChangeIndex", "PhaseClass", "Metric")) %>%
  join_duration_qc(epoch_duration_qc, by_cols = c("AnimalNum", "Group", "Sex", "CageChange")) %>%
  mutate(
    recovery_slope_per_hour = recovery_slope,
    stabilization_time_hours = stabilization_time_bins * infer_bin_size_sec(behav) / 3600,
    adaptation_half_life_hours = adaptation_half_life_bins * infer_bin_size_sec(behav) / 3600,
    StableForMainText = !short_epoch %in% TRUE & (is.na(cage_change_duration_class) | cage_change_duration_class != "short")
  )

interpretation_guide <- tibble(
  Feature = c("early_late_shift", "volatility_decay", "recovery_slope_per_hour", "stabilization_time_hours", "adaptation_half_life_hours", "control_convergence"),
  BiologicalInterpretation = c(
    "Within-epoch change from early response to late adaptation.",
    "Reduction in local temporal instability across the epoch.",
    "Directional behavioral drift after regrouping.",
    "Time until behavior remains near its late-epoch reference level.",
    "Approximate time to halfway between early response and late response.",
    "Movement toward the contemporaneous control trajectory."
  ),
  ClaimType = c("descriptive", "descriptive", "descriptive", "descriptive", "descriptive", "associative"),
  ReviewerRisk = c("medium", "medium", "medium", "medium", "medium", "medium"),
  PreferredLanguage = c("adaptation shift", "stabilization of volatility", "post-regrouping slope", "stabilization time", "adaptation half-life", "distance-to-control convergence")
)

write_table(adaptation_features, file.path(output_dir, "tables", "adaptation_kinetics_features.csv"))
write_table(distance_to_control, file.path(output_dir, "tables", "distance_to_control_trajectories.csv"))
write_table(interpretation_guide, file.path(output_dir, "tables", "adaptation_kinetics_interpretation_guide.csv"))
generate_duration_sensitivity_outputs(adaptation_features, epoch_duration_qc, output_dir, "adaptation_kinetics_features", by_cols = c("AnimalNum", "Group", "Sex", "CageChange"))

plot_tbl <- adaptation_features %>%
  filter(Metric == "Movement", PhaseClass == "Active", is.finite(volatility_decay))

if (nrow(plot_tbl) > 0) {
  p <- ggplot(plot_tbl, aes(Group, volatility_decay, fill = Group, colour = Group)) +
    geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey55") +
    geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.75, linewidth = 0.25) +
    geom_jitter(width = 0.06, size = 0.9, alpha = 0.75) +
    facet_grid(Sex ~ CageChange) +
    scale_fill_manual(values = mmm_group_colors, drop = FALSE) +
    scale_colour_manual(values = mmm_group_colors, drop = FALSE) +
    labs(title = "Adaptation kinetics after regrouping", subtitle = "Movement volatility decay in active phases", x = NULL, y = "Early RMSSD - late RMSSD") +
    make_publication_theme(base_size = 7) +
    theme(legend.position = "none")
  save_plot_svg_pdf(p, file.path(output_dir, "figures/publication_panels/Fig15_adaptation_kinetics_overview"), width = 140, height = 90)
}

message("Behavioral adaptation kinetics complete: ", output_dir)
