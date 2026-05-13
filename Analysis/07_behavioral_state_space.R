# ================================================================
# Behavioral State-Space Analysis
# MMMSociability
# ================================================================
# Goal:
#   Infer latent behavioral states from movement / entropy /
#   proximity dynamics without requiring pose estimation, then test
#   whether group, sex, phase, or cage-change context shifts state
#   occupancy and switching.
#
# Input expectation:
#   Run Analysis/03_build_multiscale_behavior_metrics.R first.
#
# Recommended scale:
#   5-10 min bins. Phase-level data are too coarse for state switching.
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cluster)
  library(purrr)
  library(readr)
  library(tibble)
})

source("C:/Users/topohl/Documents/GitHub/MMMSociability/Functions/behavioral_dynamics_helpers.R")
source("C:/Users/topohl/Documents/GitHub/MMMSociability/Functions/behavioral_dynamics_stats_helpers.R")

# ------------------------------------------------
# USER INPUT
# ------------------------------------------------

bin_level <- "1min_based"
input_file <- file.path("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/analysis_ready/03_derived_metrics", bin_level, "all_behavior_metrics.csv")
output_dir <- file.path("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/analysis_ready/06_behavioral_dynamics/state_space", bin_level)
n_states <- 4

# Use normalized proximity for state-space analyses. Raw proximity seconds scale
# with bin size and can dominate k-means/PCA for trivial duration reasons.
proximity_col <- "ProximityFraction"

group_colors <- c(
  "CON" = "#3d3b6e",
  "SUS" = "#e63947",
  "RES" = "#C6C3BB"
)
group_levels <- c("CON", "RES", "SUS")
pairwise_contrasts <- c("RES-CON", "SUS-CON", "SUS-RES")
state_colors <- c("#2F4858", "#7E9F35", "#F2A65A", "#B23A48", "#6D597A", "#4D908E")

# ------------------------------------------------
# HELPERS
# ------------------------------------------------

parse_cage_change_index <- function(x) {
  idx <- suppressWarnings(as.integer(gsub("\\D+", "", as.character(x))))
  fallback <- match(as.character(x), sort(unique(as.character(x))))
  ifelse(is.na(idx), fallback, idx)
}

sig_from_p <- function(p) {
  case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    p < 0.10 ~ "trend",
    TRUE ~ "ns"
  )
}

format_p_label <- function(p) {
  case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "p<.001",
    TRUE ~ paste0("p=", sub("^0", "", formatC(p, format = "f", digits = 3)))
  )
}

add_stat_columns <- function(tbl) {
  if (nrow(tbl) == 0) return(tbl)
  tbl %>%
    mutate(
      p.sign = sig_from_p(p.adjust_bh_family),
      Significant = p.adjust_bh_family < 0.05,
      ReportingP = p.adjust_bh_family,
      ReportingSignificance = p.sign,
      ReportingCorrection = "BH across all pairwise group contrasts"
    ) %>%
    select(-any_of(c("is_primary_contrast", "p.adjust_bh_primary", "p.sign_primary")))
}

write_stats_package <- function(dat, analysis_name, value_cols, by_cols, summary_group_cols) {
  summary_tbl <- make_dynamics_group_summary(dat, value_cols = value_cols, group_cols = summary_group_cols)
  contrast_tbl <- make_dynamics_group_contrasts(dat, value_cols = value_cols, by_cols = by_cols) %>%
    add_stat_columns()

  write_table(summary_tbl, file.path(output_dir, "stats_tables", paste0(analysis_name, "_group_summary.csv")))
  write_table(contrast_tbl, file.path(output_dir, "stats_tables", paste0(analysis_name, "_group_contrasts.csv")))

  invisible(list(summary = summary_tbl, contrasts = contrast_tbl))
}

make_publication_theme <- function(base_size = 7) {
  make_nature_theme(base_size = base_size) +
    theme(
      legend.position = "top",
      legend.key.width = unit(8, "mm"),
      panel.grid.major.y = element_line(linewidth = 0.15, colour = "grey92"),
      panel.grid.major.x = element_blank(),
      panel.spacing = unit(1.0, "lines"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1),
      plot.subtitle = element_text(hjust = 0, size = base_size)
    )
}

summarise_mean_ci <- function(dat, value_col, group_cols) {
  dat %>%
    filter(is.finite(.data[[value_col]])) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n_animals = n_distinct(AnimalNum),
      n_obs = n(),
      mean = mean(.data[[value_col]], na.rm = TRUE),
      sd = sd(.data[[value_col]], na.rm = TRUE),
      sem = sd / sqrt(n_obs),
      ci95 = 1.96 * sem,
      ymin = mean - ci95,
      ymax = mean + ci95,
      .groups = "drop"
    )
}

make_pairwise_brackets <- function(data,
                                   contrast_tbl,
                                   outcome,
                                   facet_cols,
                                   value_col = NULL,
                                   p_threshold = 0.05,
                                   include_trends = FALSE) {
  if (is.null(value_col)) value_col <- outcome
  empty_brackets <- tibble(
    contrast = character(), x_start = numeric(), x_end = numeric(),
    y = numeric(), y_tip = numeric(), label = character()
  )
  if (length(facet_cols) > 0) {
    for (nm in facet_cols) empty_brackets[[nm]] <- data[[nm]][0]
  }
  if (nrow(contrast_tbl) == 0 || !value_col %in% names(data)) return(empty_brackets)

  threshold <- if (include_trends) 0.10 else p_threshold

  brackets <- contrast_tbl %>%
    filter(
      Outcome == outcome,
      contrast %in% pairwise_contrasts,
      !is.na(ReportingP),
      ReportingP < threshold
    ) %>%
    mutate(
      x_start = match(group_ref, group_levels),
      x_end = match(group_comp, group_levels),
      label = format_p_label(ReportingP)
    ) %>%
    filter(is.finite(x_start), is.finite(x_end))

  if (nrow(brackets) == 0) return(empty_brackets)

  y_tbl <- data %>%
    group_by(across(all_of(facet_cols))) %>%
    summarise(
      y_max = max(.data[[value_col]], na.rm = TRUE),
      y_min = min(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y_range = if_else(is.finite(y_max - y_min) & y_max > y_min, y_max - y_min, abs(y_max)),
      y_step = if_else(is.finite(y_range) & y_range > 0, 0.18 * y_range, 0.1),
      y_tip_size = 0.035 * y_step
    )

  brackets %>%
    group_by(across(all_of(facet_cols))) %>%
    arrange(ReportingP, .by_group = TRUE) %>%
    mutate(bracket_rank = row_number()) %>%
    ungroup() %>%
    left_join(y_tbl, by = facet_cols) %>%
    mutate(
      y = y_max + bracket_rank * y_step,
      y_tip = y - y_tip_size
    ) %>%
    select(any_of(facet_cols), contrast, x_start, x_end, y, y_tip, label, ReportingP, ReportingSignificance)
}

add_pairwise_brackets <- function(plot, bracket_tbl, text_size = 1.8) {
  if (nrow(bracket_tbl) == 0) return(plot)
  plot +
    geom_segment(
      data = bracket_tbl,
      aes(x = x_start, xend = x_end, y = y, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.25
    ) +
    geom_segment(
      data = bracket_tbl,
      aes(x = x_start, xend = x_start, y = y, yend = y_tip),
      inherit.aes = FALSE,
      linewidth = 0.25
    ) +
    geom_segment(
      data = bracket_tbl,
      aes(x = x_end, xend = x_end, y = y, yend = y_tip),
      inherit.aes = FALSE,
      linewidth = 0.25
    ) +
    geom_text(
      data = bracket_tbl,
      aes(x = (x_start + x_end) / 2, y = y, label = label),
      inherit.aes = FALSE,
      vjust = -0.25,
      size = text_size
    )
}

make_state_profile_labels <- function(centroids) {
  centroids %>%
    rowwise() %>%
    mutate(
      dominant_feature = c("Movement", "Entropy", "Proximity")[which.max(c(Movement_z, Entropy_z, Proximity_z))],
      dominant_direction = if_else(max(c(Movement_z, Entropy_z, Proximity_z)) >= abs(min(c(Movement_z, Entropy_z, Proximity_z))), "high", "low"),
      dominant_feature = if_else(dominant_direction == "low", c("Movement", "Entropy", "Proximity")[which.min(c(Movement_z, Entropy_z, Proximity_z))], dominant_feature),
      StateLabel = paste0("S", State, " (", dominant_direction, " ", tolower(dominant_feature), ")")
    ) %>%
    ungroup() %>%
    select(State, StateLabel)
}

# ------------------------------------------------
# LOAD + STANDARDIZE
# ------------------------------------------------

raw_dat <- read_behavior_table(input_file)
if (!proximity_col %in% names(raw_dat)) proximity_col <- "Proximity"

behav <- standardize_behavior_columns(raw_dat, proximity_col = proximity_col) %>%
  mutate(
    Group = factor(as.character(Group), levels = group_levels),
    CageChangeIndex = parse_cage_change_index(CageChange)
  )

ensure_dir(output_dir)
ensure_dir(file.path(output_dir, "tables"))
ensure_dir(file.path(output_dir, "stats_tables"))
ensure_dir(file.path(output_dir, "figures"))
ensure_dir(file.path(output_dir, "figures", "publication"))
ensure_dir(file.path(output_dir, "figures", "publication", "panels"))
ensure_dir(file.path(output_dir, "figures", "publication", "overview"))

# ------------------------------------------------
# STATE FEATURES
# ------------------------------------------------

state_dat <- behav %>%
  mutate(
    Movement_z = z_within_metric(Movement),
    Entropy_z = z_within_metric(Entropy),
    Proximity_z = z_within_metric(Proximity),
    BinLevel = bin_level,
    ProximityInput = proximity_col
  ) %>%
  drop_na(Movement_z, Entropy_z, Proximity_z)

feature_matrix <- state_dat %>%
  select(Movement_z, Entropy_z, Proximity_z) %>%
  as.matrix()

set.seed(123)
km <- kmeans(feature_matrix, centers = n_states, nstart = 100, iter.max = 100)
state_dat$State <- factor(km$cluster, levels = seq_len(n_states))

set.seed(124)
silhouette_n <- min(10000L, nrow(state_dat))
silhouette_idx <- sample(seq_len(nrow(state_dat)), size = silhouette_n)
sil <- cluster::silhouette(as.integer(state_dat$State[silhouette_idx]), dist(feature_matrix[silhouette_idx, , drop = FALSE]))
state_dat$SilhouetteWidth <- NA_real_
state_dat$SilhouetteWidth[silhouette_idx] <- sil[, "sil_width"]

model_quality_tbl <- tibble(
  BinLevel = bin_level,
  ProximityInput = proximity_col,
  n_states = n_states,
  n_bins = nrow(state_dat),
  silhouette_sample_n = silhouette_n,
  total_withinss = km$tot.withinss,
  between_total_ss_ratio = km$betweenss / km$totss,
  mean_silhouette_width = mean(state_dat$SilhouetteWidth, na.rm = TRUE),
  median_silhouette_width = median(state_dat$SilhouetteWidth, na.rm = TRUE)
)

write_table(model_quality_tbl, file.path(output_dir, "tables", "state_model_quality.csv"))

# ------------------------------------------------
# STATE CENTROIDS + PROFILE LABELS
# ------------------------------------------------

centroids <- state_dat %>%
  group_by(BinLevel, ProximityInput, State) %>%
  summarise(
    Movement_z = mean(Movement_z, na.rm = TRUE),
    Entropy_z = mean(Entropy_z, na.rm = TRUE),
    Proximity_z = mean(Proximity_z, na.rm = TRUE),
    mean_silhouette_width = mean(SilhouetteWidth, na.rm = TRUE),
    n_bins = n(),
    .groups = "drop"
  )

state_labels <- make_state_profile_labels(centroids)
state_levels <- state_labels$StateLabel[match(levels(state_dat$State), state_labels$State)]

state_dat <- state_dat %>%
  left_join(state_labels, by = "State") %>%
  mutate(
    StateLabel = factor(StateLabel, levels = state_levels),
    State = factor(State, levels = seq_len(n_states))
  )

centroids <- centroids %>%
  left_join(state_labels, by = "State") %>%
  mutate(StateLabel = factor(StateLabel, levels = state_levels))

centroid_long <- centroids %>%
  pivot_longer(
    cols = c(Movement_z, Entropy_z, Proximity_z),
    names_to = "Feature",
    values_to = "CentroidZ"
  ) %>%
  mutate(
    Feature = recode(Feature, Movement_z = "Movement", Entropy_z = "Entropy", Proximity_z = "Proximity"),
    Feature = factor(Feature, levels = c("Movement", "Entropy", "Proximity"))
  )

write_table(centroids, file.path(output_dir, "tables", "state_centroids.csv"))
write_table(centroid_long, file.path(output_dir, "tables", "state_centroid_feature_profiles.csv"))

# ------------------------------------------------
# OCCUPANCY + STATE DIVERSITY
# ------------------------------------------------

occupancy_tbl <- state_dat %>%
  count(BinLevel, ProximityInput, Group, Sex, Phase, CageChange, CageChangeIndex, AnimalNum, State, StateLabel, name = "n_bins_state") %>%
  group_by(BinLevel, ProximityInput, Group, Sex, Phase, CageChange, CageChangeIndex, AnimalNum) %>%
  mutate(frac_time = n_bins_state / sum(n_bins_state)) %>%
  ungroup()

state_diversity_tbl <- occupancy_tbl %>%
  group_by(BinLevel, ProximityInput, Group, Sex, Phase, CageChange, CageChangeIndex, AnimalNum) %>%
  summarise(
    state_entropy = -sum(frac_time * log(frac_time + 1e-9)),
    normalized_state_entropy = state_entropy / log(n_states),
    n_observed_states = sum(frac_time > 0),
    max_state_fraction = max(frac_time, na.rm = TRUE),
    dominant_state = as.character(StateLabel[which.max(frac_time)]),
    .groups = "drop"
  )

write_table(occupancy_tbl, file.path(output_dir, "tables", "state_occupancy.csv"))
write_table(state_diversity_tbl, file.path(output_dir, "tables", "state_diversity_metrics.csv"))

occupancy_stats <- write_stats_package(
  occupancy_tbl,
  analysis_name = "state_occupancy",
  value_cols = "frac_time",
  by_cols = c("BinLevel", "ProximityInput", "StateLabel", "Sex", "Phase", "CageChange"),
  summary_group_cols = c("BinLevel", "ProximityInput", "StateLabel", "Sex", "Phase", "CageChange", "Group")
)

diversity_stats <- write_stats_package(
  state_diversity_tbl,
  analysis_name = "state_diversity",
  value_cols = c("state_entropy", "normalized_state_entropy", "n_observed_states", "max_state_fraction"),
  by_cols = c("BinLevel", "ProximityInput", "Sex", "Phase", "CageChange"),
  summary_group_cols = c("BinLevel", "ProximityInput", "Sex", "Phase", "CageChange", "Group")
)

# ------------------------------------------------
# TRANSITIONS + SWITCHING RATE
# ------------------------------------------------

transition_events <- state_dat %>%
  group_by(BinLevel, ProximityInput, Group, Sex, Phase, CageChange, CageChangeIndex, AnimalNum) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  mutate(
    NextState = lead(State),
    NextStateLabel = lead(StateLabel),
    StateSwitch = State != lead(State)
  ) %>%
  filter(!is.na(NextState)) %>%
  ungroup()

transition_tbl <- transition_events %>%
  count(BinLevel, ProximityInput, Group, Sex, Phase, CageChange, CageChangeIndex, AnimalNum, State, StateLabel, NextState, NextStateLabel, name = "Transitions") %>%
  group_by(BinLevel, ProximityInput, Group, Sex, Phase, CageChange, CageChangeIndex, AnimalNum, State, StateLabel) %>%
  mutate(transition_probability = Transitions / sum(Transitions)) %>%
  ungroup()

switch_tbl <- transition_events %>%
  group_by(BinLevel, ProximityInput, Group, Sex, Phase, CageChange, CageChangeIndex, AnimalNum) %>%
  summarise(
    switch_rate = mean(StateSwitch, na.rm = TRUE),
    stay_rate = 1 - switch_rate,
    n_switches = sum(StateSwitch, na.rm = TRUE),
    n_transitions = n(),
    .groups = "drop"
  )

transition_summary_tbl <- transition_tbl %>%
  group_by(BinLevel, ProximityInput, Group, Sex, Phase, CageChange, CageChangeIndex, StateLabel, NextStateLabel) %>%
  summarise(
    n_animals = n_distinct(AnimalNum),
    total_transitions = sum(Transitions, na.rm = TRUE),
    mean_transition_probability = mean(transition_probability, na.rm = TRUE),
    .groups = "drop"
  )

write_table(transition_tbl, file.path(output_dir, "tables", "state_transition_counts.csv"))
write_table(transition_summary_tbl, file.path(output_dir, "tables", "state_transition_probability_summary.csv"))
write_table(switch_tbl, file.path(output_dir, "tables", "state_switching_metrics.csv"))

switch_stats <- write_stats_package(
  switch_tbl,
  analysis_name = "state_switching",
  value_cols = c("switch_rate", "stay_rate", "n_switches"),
  by_cols = c("BinLevel", "ProximityInput", "Sex", "Phase", "CageChange"),
  summary_group_cols = c("BinLevel", "ProximityInput", "Sex", "Phase", "CageChange", "Group")
)

# ------------------------------------------------
# PCA LANDSCAPE
# ------------------------------------------------

pca_fit <- prcomp(feature_matrix, center = TRUE, scale. = TRUE)
state_dat$PC1 <- pca_fit$x[, 1]
state_dat$PC2 <- pca_fit$x[, 2]

pca_var_tbl <- tibble(
  PC = paste0("PC", seq_along(pca_fit$sdev)),
  variance_explained = pca_fit$sdev^2 / sum(pca_fit$sdev^2),
  cumulative_variance_explained = cumsum(variance_explained)
)

pca_loadings_tbl <- as_tibble(pca_fit$rotation, rownames = "Feature") %>%
  mutate(Feature = recode(Feature, Movement_z = "Movement", Entropy_z = "Entropy", Proximity_z = "Proximity"))

write_table(pca_var_tbl, file.path(output_dir, "tables", "state_space_pca_variance.csv"))
write_table(pca_loadings_tbl, file.path(output_dir, "tables", "state_space_pca_loadings.csv"))

pc1_lab <- paste0("PC1 (", round(100 * pca_var_tbl$variance_explained[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(100 * pca_var_tbl$variance_explained[2], 1), "%)")

p_pca <- state_dat %>%
  ggplot(aes(PC1, PC2, colour = StateLabel)) +
  geom_point(alpha = 0.35, size = 0.55) +
  stat_ellipse(type = "norm", linewidth = 0.25, alpha = 0.75, show.legend = FALSE) +
  facet_grid(Group ~ Phase) +
  labs(
    title = "Behavioral state-space landscape",
    subtitle = paste0("K-means states from z-scored movement, entropy, and proximity; mean silhouette = ", round(model_quality_tbl$mean_silhouette_width, 2)),
    x = pc1_lab,
    y = pc2_lab,
    colour = NULL
  ) +
  scale_colour_manual(values = state_colors[seq_len(n_states)], drop = FALSE) +
  make_publication_theme(base_size = 6) +
  theme(panel.grid.major = element_blank())

save_plot_svg_pdf(p_pca, file.path(output_dir, "figures", "behavioral_state_space_pca"), width = 190, height = 135)

# ------------------------------------------------
# STATE PROFILE PLOTS
# ------------------------------------------------

p_centroids <- centroid_long %>%
  ggplot(aes(Feature, CentroidZ, fill = StateLabel)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey45") +
  geom_col(width = 0.72, colour = "grey20", linewidth = 0.18, show.legend = FALSE) +
  facet_wrap(~StateLabel, nrow = 1) +
  labs(
    title = "Behavioral state feature profiles",
    subtitle = "Centroids are z-scored within the full analysis set",
    x = NULL,
    y = "Centroid z-score"
  ) +
  scale_fill_manual(values = state_colors[seq_len(n_states)], drop = FALSE) +
  make_publication_theme(base_size = 6) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid.major.x = element_blank())

save_plot_svg_pdf(p_centroids, file.path(output_dir, "figures", "state_centroid_profiles"), width = 180, height = 70)

# ------------------------------------------------
# STATE OCCUPANCY PLOTS
# ------------------------------------------------

occupancy_brackets <- make_pairwise_brackets(
  occupancy_tbl,
  occupancy_stats$contrasts,
  outcome = "frac_time",
  facet_cols = c("Sex", "Phase", "StateLabel"),
  value_col = "frac_time"
)

p_occ <- occupancy_tbl %>%
  ggplot(aes(Group, frac_time, fill = Group)) +
  geom_violin(alpha = 0.55, linewidth = 0.25, trim = FALSE, colour = "grey25") +
  geom_jitter(width = 0.08, size = 0.75, alpha = 0.70, shape = 21, colour = "grey20", stroke = 0.12, show.legend = FALSE) +
  facet_grid(Sex + Phase ~ StateLabel, scales = "free_y") +
  labs(
    title = "Behavioral state occupancy",
    subtitle = "Each point is one animal within a cage-change and phase; brackets show BH-adjusted pairwise tests",
    y = "Fraction of time",
    x = NULL
  ) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  make_publication_theme(base_size = 6)

p_occ <- add_pairwise_brackets(p_occ, occupancy_brackets, text_size = 1.55)

save_plot_svg_pdf(p_occ, file.path(output_dir, "figures", "state_occupancy"), width = 210, height = 150)

occupancy_summary <- summarise_mean_ci(
  occupancy_tbl,
  "frac_time",
  c("Group", "Sex", "Phase", "CageChange", "CageChangeIndex", "StateLabel")
) %>%
  mutate(CageChange = factor(as.character(CageChange), levels = unique(as.character(CageChange[order(CageChangeIndex)]))))

write_table(occupancy_summary, file.path(output_dir, "tables", "state_occupancy_mean_ci_summary.csv"))

p_occ_traj <- ggplot(occupancy_summary, aes(CageChange, mean, colour = Group, group = Group)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.08, linewidth = 0.28) +
  geom_line(linewidth = 0.55) +
  geom_point(aes(fill = Group), size = 1.45, shape = 21, stroke = 0.22) +
  facet_grid(Sex + Phase ~ StateLabel, scales = "free_y") +
  labs(
    title = "State occupancy trajectories across cage changes",
    subtitle = "Group mean +/- 95% CI",
    x = "Cage change",
    y = "Fraction of time"
  ) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  make_publication_theme(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

save_plot_svg_pdf(p_occ_traj, file.path(output_dir, "figures", "state_occupancy_trajectories"), width = 210, height = 150)

# ------------------------------------------------
# SWITCHING + DIVERSITY PLOTS
# ------------------------------------------------

switch_brackets <- make_pairwise_brackets(
  switch_tbl,
  switch_stats$contrasts,
  outcome = "switch_rate",
  facet_cols = c("Sex", "Phase"),
  value_col = "switch_rate"
)

p_switch <- switch_tbl %>%
  ggplot(aes(Group, switch_rate, fill = Group)) +
  geom_violin(alpha = 0.55, linewidth = 0.25, trim = FALSE, colour = "grey25") +
  geom_jitter(width = 0.08, size = 0.75, alpha = 0.70, shape = 21, colour = "grey20", stroke = 0.12, show.legend = FALSE) +
  facet_grid(Sex ~ Phase, scales = "free_y") +
  labs(
    title = "Behavioral state switching",
    subtitle = "Probability that consecutive bins occupy different latent states",
    y = "Switching rate",
    x = NULL
  ) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  make_publication_theme(base_size = 6)

p_switch <- add_pairwise_brackets(p_switch, switch_brackets, text_size = 1.8)

save_plot_svg_pdf(p_switch, file.path(output_dir, "figures", "state_switching_rate"), width = 160, height = 110)

diversity_brackets <- make_pairwise_brackets(
  state_diversity_tbl,
  diversity_stats$contrasts,
  outcome = "normalized_state_entropy",
  facet_cols = c("Sex", "Phase"),
  value_col = "normalized_state_entropy"
)

p_diversity <- state_diversity_tbl %>%
  ggplot(aes(Group, normalized_state_entropy, fill = Group)) +
  geom_violin(alpha = 0.55, linewidth = 0.25, trim = FALSE, colour = "grey25") +
  geom_jitter(width = 0.08, size = 0.75, alpha = 0.70, shape = 21, colour = "grey20", stroke = 0.12, show.legend = FALSE) +
  facet_grid(Sex ~ Phase, scales = "free_y") +
  labs(
    title = "Behavioral state diversity",
    subtitle = "Normalized entropy of state occupancy within animal, phase, and cage change",
    y = "Normalized state entropy",
    x = NULL
  ) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  make_publication_theme(base_size = 6)

p_diversity <- add_pairwise_brackets(p_diversity, diversity_brackets, text_size = 1.8)

save_plot_svg_pdf(p_diversity, file.path(output_dir, "figures", "state_diversity_entropy"), width = 160, height = 110)

# ------------------------------------------------
# EFFECT-SIZE MAPS + TRANSITION MAP
# ------------------------------------------------

occupancy_effect_tbl <- occupancy_stats$contrasts %>%
  filter(contrast %in% pairwise_contrasts, status == "tested") %>%
  mutate(
    StateLabel = factor(StateLabel, levels = state_levels),
    CagePhase = paste(as.character(CageChange), Phase, sep = "\n"),
    CagePhase = factor(CagePhase, levels = unique(CagePhase[order(CageChange, Phase)])),
    contrast = factor(contrast, levels = pairwise_contrasts),
    SigLabel = if_else(ReportingP < 0.05, "*", "")
  )

write_table(occupancy_effect_tbl, file.path(output_dir, "tables", "state_occupancy_effect_size_heatmap_data.csv"))

p_occ_heat <- ggplot(occupancy_effect_tbl, aes(CagePhase, StateLabel, fill = cohen_d)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = SigLabel), size = 2.3, colour = "black") +
  facet_grid(Sex ~ contrast, scales = "free_x", space = "free_x") +
  labs(
    title = "State occupancy effect-size map",
    subtitle = "Cohen's d for pairwise contrasts; asterisk marks BH-adjusted FDR < 0.05",
    x = NULL,
    y = NULL,
    fill = "Cohen's d"
  ) +
  scale_fill_gradient2(low = "#3d3b6e", mid = "white", high = "#e63947", midpoint = 0, na.value = "grey90") +
  make_publication_theme(base_size = 6) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  )

save_plot_svg_pdf(p_occ_heat, file.path(output_dir, "figures", "publication", "overview", "state_occupancy_effect_size_heatmap"), width = 190, height = 135)

switch_effect_tbl <- switch_stats$contrasts %>%
  filter(contrast %in% pairwise_contrasts, status == "tested") %>%
  mutate(
    Outcome = factor(Outcome, levels = c("switch_rate", "stay_rate", "n_switches"), labels = c("Switch rate", "Stay rate", "N switches")),
    CagePhase = paste(as.character(CageChange), Phase, sep = "\n"),
    CagePhase = factor(CagePhase, levels = unique(CagePhase[order(CageChange, Phase)])),
    contrast = factor(contrast, levels = pairwise_contrasts),
    SigLabel = if_else(ReportingP < 0.05, "*", "")
  )

write_table(switch_effect_tbl, file.path(output_dir, "tables", "state_switching_effect_size_heatmap_data.csv"))

p_switch_heat <- ggplot(switch_effect_tbl, aes(CagePhase, Outcome, fill = cohen_d)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = SigLabel), size = 2.3, colour = "black") +
  facet_grid(Sex ~ contrast, scales = "free_x", space = "free_x") +
  labs(
    title = "State-transition effect-size map",
    subtitle = "Cohen's d for pairwise contrasts; asterisk marks BH-adjusted FDR < 0.05",
    x = NULL,
    y = NULL,
    fill = "Cohen's d"
  ) +
  scale_fill_gradient2(low = "#3d3b6e", mid = "white", high = "#e63947", midpoint = 0, na.value = "grey90") +
  make_publication_theme(base_size = 6) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  )

save_plot_svg_pdf(p_switch_heat, file.path(output_dir, "figures", "publication", "overview", "state_switching_effect_size_heatmap"), width = 190, height = 120)

transition_plot_tbl <- transition_summary_tbl %>%
  group_by(Group, Sex, Phase, StateLabel, NextStateLabel) %>%
  summarise(
    mean_transition_probability = weighted.mean(mean_transition_probability, w = pmax(total_transitions, 1), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    StateLabel = factor(StateLabel, levels = state_levels),
    NextStateLabel = factor(NextStateLabel, levels = state_levels)
  )

p_transition <- ggplot(transition_plot_tbl, aes(NextStateLabel, StateLabel, fill = mean_transition_probability)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  facet_grid(Sex + Phase ~ Group) +
  labs(
    title = "Behavioral state transition probabilities",
    subtitle = "Rows are current states; columns are next-bin states",
    x = "Next state",
    y = "Current state",
    fill = "Probability"
  ) +
  scale_fill_gradient(low = "grey95", high = "#2F4858", limits = c(0, NA)) +
  make_publication_theme(base_size = 6) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  )

save_plot_svg_pdf(p_transition, file.path(output_dir, "figures", "state_transition_probability_heatmap"), width = 190, height = 140)

message("Behavioral state-space analysis complete. State occupancy, diversity, switching, PCA, and effect-size outputs written.")
