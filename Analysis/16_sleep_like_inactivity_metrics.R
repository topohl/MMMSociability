# ================================================================
# Sleep-like Inactivity / Quiescence Metrics
# MMMSociability
# ================================================================
# Goal:
#   Quantify rest-like inactivity structure from RFID home-cage behavior.
#
# Caveat:
#   These are not EEG-validated sleep metrics. Use "sleep-like",
#   "rest-like" or "quiescence" language only.
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
output_dir <- file.path(base_dir, "analysis_ready/16_sleep_like_inactivity_metrics", bin_level)
proximity_col <- "ProximityFraction"
low_activity_quantile <- 0.20
prolonged_inactivity_min <- 30

output_dirs <- analysis_output_dirs(output_dir)
write_output_manifest(
  output_dir,
  script_name = "16_sleep_like_inactivity_metrics.R",
  analysis_name = "sleep-like inactivity and quiescence fragmentation",
  primary_tables = c(
    "tables/sleep_like_inactivity_features.csv",
    "tables/inactivity_bout_architecture.csv",
    "tables/sleep_like_inactivity_interpretation_guide.csv"
  ),
  primary_figures = c("figures/publication_panels/Fig16_quiescence_fragmentation.svg"),
  notes = c("Do not claim EEG-validated sleep; interpret as rest-like inactivity/quiescence structure.")
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

bin_size_sec <- infer_bin_size_sec(behav)
bin_size_min <- bin_size_sec / 60
prolonged_bins <- ceiling(prolonged_inactivity_min / bin_size_min)
epoch_duration_qc <- write_epoch_duration_qc(behav, output_dir, metric_source = "16_sleep_like_inactivity_metrics", bin_size_sec = bin_size_sec)

animal_thresholds <- behav %>%
  group_by(AnimalNum) %>%
  summarise(
    low_activity_threshold = quantile(Movement, low_activity_quantile, na.rm = TRUE, names = FALSE),
    zero_like_threshold = if_else(any(Movement == 0, na.rm = TRUE), 0, low_activity_threshold),
    .groups = "drop"
  )

inactivity_bins <- behav %>%
  left_join(animal_thresholds, by = "AnimalNum") %>%
  mutate(
    inactive_like = is.finite(Movement) & Movement <= low_activity_threshold,
    zero_like_inactive = is.finite(Movement) & Movement <= zero_like_threshold
  )

make_bout_table <- function(dat) {
  dat <- dat %>% arrange(TimeIndex)
  r <- rle(dat$inactive_like %in% TRUE)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  tibble(
    bout_id = seq_along(r$lengths),
    inactive_like = r$values,
    start_bin = starts,
    end_bin = ends,
    duration_bins = r$lengths,
    duration_min = duration_bins * bin_size_min
  ) %>%
    filter(inactive_like)
}

bout_architecture <- inactivity_bins %>%
  group_by(BinLevel, AnimalNum, Group, Sex, CageChange, CageChangeIndex, PhaseClass) %>%
  group_modify(~make_bout_table(.x)) %>%
  ungroup()

inactivity_features <- inactivity_bins %>%
  group_by(BinLevel, AnimalNum, Group, Sex, CageChange, CageChangeIndex, PhaseClass) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  summarise(
    n_bins = n(),
    inactivity_fraction = mean(inactive_like, na.rm = TRUE),
    zero_like_inactivity_fraction = mean(zero_like_inactive, na.rm = TRUE),
    active_inactive_transition_rate = if_else(n() >= 3, mean(inactive_like != lag(inactive_like), na.rm = TRUE), NA_real_),
    inactive_to_active_rate = if_else(n() >= 3, mean(inactive_like[-1] == FALSE & inactive_like[-n()] == TRUE, na.rm = TRUE), NA_real_),
    active_to_inactive_rate = if_else(n() >= 3, mean(inactive_like[-1] == TRUE & inactive_like[-n()] == FALSE, na.rm = TRUE), NA_real_),
    .groups = "drop"
  ) %>%
  left_join(
    bout_architecture %>%
      group_by(BinLevel, AnimalNum, Group, Sex, CageChange, CageChangeIndex, PhaseClass) %>%
      summarise(
        inactivity_bout_count = n(),
        mean_inactivity_bout_min = mean(duration_min, na.rm = TRUE),
        median_inactivity_bout_min = median(duration_min, na.rm = TRUE),
        max_inactivity_bout_min = max(duration_min, na.rm = TRUE),
        prolonged_inactivity_episodes = sum(duration_bins >= prolonged_bins, na.rm = TRUE),
        prolonged_inactivity_fraction = mean(duration_bins >= prolonged_bins, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("BinLevel", "AnimalNum", "Group", "Sex", "CageChange", "CageChangeIndex", "PhaseClass")
  ) %>%
  mutate(
    inactivity_bout_count = replace_na(inactivity_bout_count, 0L),
    inactivity_fragmentation = active_inactive_transition_rate,
    prolonged_inactivity_episodes = replace_na(prolonged_inactivity_episodes, 0L)
  ) %>%
  join_duration_qc(epoch_duration_qc, by_cols = c("AnimalNum", "Group", "Sex", "CageChange")) %>%
  normalize_counts_to_rates(c("inactivity_bout_count", "prolonged_inactivity_episodes"), duration_col = "total_observation_duration_hours")

interpretation_guide <- tibble(
  Feature = c("inactivity_fraction", "mean_inactivity_bout_min", "inactivity_fragmentation", "prolonged_inactivity_episodes_per_hour", "active_inactive_transition_rate"),
  AllowedInterpretation = c(
    "Fraction of bins in a low-activity/quiescent state.",
    "Continuity of rest-like inactivity bouts.",
    "Fragmentation of inactivity into shorter or more interrupted episodes.",
    "Frequency of prolonged quiescence episodes.",
    "Switching between active and inactive-like states."
  ),
  DisallowedInterpretation = c(
    "Total sleep time",
    "Sleep bout length",
    "Sleep fragmentation",
    "Validated sleep episodes",
    "Sleep-stage transition rate"
  ),
  ReviewerRisk = c("medium", "medium", "medium", "medium", "medium")
)

write_table(inactivity_features, file.path(output_dir, "tables", "sleep_like_inactivity_features.csv"))
write_table(bout_architecture, file.path(output_dir, "tables", "inactivity_bout_architecture.csv"))
write_table(interpretation_guide, file.path(output_dir, "tables", "sleep_like_inactivity_interpretation_guide.csv"))
generate_duration_sensitivity_outputs(inactivity_features, epoch_duration_qc, output_dir, "sleep_like_inactivity_features", by_cols = c("AnimalNum", "Group", "Sex", "CageChange"))

plot_tbl <- inactivity_features %>% filter(is.finite(inactivity_fragmentation))
if (nrow(plot_tbl) > 0) {
  p <- ggplot(plot_tbl, aes(Group, inactivity_fragmentation, fill = Group, colour = Group)) +
    geom_boxplot(width = 0.22, outlier.shape = NA, alpha = 0.75, linewidth = 0.25) +
    geom_jitter(width = 0.06, size = 0.9, alpha = 0.75) +
    facet_grid(PhaseClass ~ Sex) +
    scale_fill_manual(values = mmm_group_colors, drop = FALSE) +
    scale_colour_manual(values = mmm_group_colors, drop = FALSE) +
    labs(title = "Sleep-like inactivity fragmentation", subtitle = "Low-activity bout switching; not EEG-validated sleep", x = NULL, y = "Active-inactive transition rate") +
    make_publication_theme(base_size = 7) +
    theme(legend.position = "none")
  save_plot_svg_pdf(p, file.path(output_dir, "figures/publication_panels/Fig16_quiescence_fragmentation"), width = 120, height = 90)
}

message("Sleep-like inactivity metrics complete: ", output_dir)
