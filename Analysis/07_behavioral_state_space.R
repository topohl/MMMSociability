# ================================================================
# Behavioral State-Space Analysis
# MMMSociability
# ================================================================
# Goal:
#   Infer latent behavioral states from movement / entropy /
#   proximity dynamics without requiring pose estimation.
#
# Input expectation:
#   Run Analysis/03_build_multiscale_behavior_metrics.R first.
#
# Recommended scale:
#   5–10 min bins. Phase-level data are too coarse for state switching.
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cluster)
  library(purrr)
  library(readr)
})

source("Functions/behavioral_dynamics_helpers.R")

# ------------------------------------------------
# USER INPUT
# ------------------------------------------------

bin_level <- "5min_based"
input_file <- file.path("analysis_ready/03_derived_metrics", bin_level, "all_behavior_metrics.csv")
output_dir <- file.path("analysis_ready/06_behavioral_dynamics/state_space", bin_level)
n_states <- 4

# Use normalized proximity for state-space analyses. Raw proximity seconds scale
# with bin size and can dominate k-means/PCA for trivial duration reasons.
proximity_col <- "ProximityFraction"

# ------------------------------------------------
# LOAD
# ------------------------------------------------

raw_dat <- read_behavior_table(input_file)
if (!proximity_col %in% names(raw_dat)) proximity_col <- "Proximity"
behav <- standardize_behavior_columns(raw_dat, proximity_col = proximity_col)

ensure_dir(output_dir)
ensure_dir(file.path(output_dir, "tables"))
ensure_dir(file.path(output_dir, "figures"))

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
km <- kmeans(feature_matrix, centers = n_states, nstart = 50)
state_dat$State <- factor(km$cluster)

# ------------------------------------------------
# STATE CENTROIDS
# ------------------------------------------------

centroids <- state_dat %>%
  group_by(BinLevel, ProximityInput, State) %>%
  summarise(
    Movement_z = mean(Movement_z, na.rm = TRUE),
    Entropy_z = mean(Entropy_z, na.rm = TRUE),
    Proximity_z = mean(Proximity_z, na.rm = TRUE),
    n_bins = n(),
    .groups = "drop"
  )

write_table(centroids, file.path(output_dir, "tables", "state_centroids.csv"))

# ------------------------------------------------
# OCCUPANCY
# ------------------------------------------------

occupancy_tbl <- state_dat %>%
  count(BinLevel, ProximityInput, Group, Sex, Phase, CageChange, AnimalNum, State) %>%
  group_by(BinLevel, ProximityInput, Group, Sex, Phase, CageChange, AnimalNum) %>%
  mutate(
    frac_time = n / sum(n),
    state_entropy = -sum(frac_time * log(frac_time + 1e-9))
  ) %>%
  ungroup()

write_table(occupancy_tbl, file.path(output_dir, "tables", "state_occupancy.csv"))

# ------------------------------------------------
# TRANSITIONS
# ------------------------------------------------

transition_tbl <- state_dat %>%
  group_by(BinLevel, ProximityInput, Group, Sex, Phase, CageChange, AnimalNum) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  mutate(
    NextState = lead(State),
    StateSwitch = State != lead(State)
  ) %>%
  filter(!is.na(NextState)) %>%
  count(BinLevel, ProximityInput, Group, Sex, Phase, CageChange, AnimalNum, State, NextState, name = "Transitions")

write_table(transition_tbl, file.path(output_dir, "tables", "state_transition_counts.csv"))

# ------------------------------------------------
# SWITCHING RATE
# ------------------------------------------------

switch_tbl <- state_dat %>%
  group_by(BinLevel, ProximityInput, Group, Sex, Phase, CageChange, AnimalNum) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  summarise(
    switch_rate = mean(State != lag(State), na.rm = TRUE),
    n_switches = sum(State != lag(State), na.rm = TRUE),
    n_bins = n(),
    .groups = "drop"
  )

write_table(switch_tbl, file.path(output_dir, "tables", "state_switching_metrics.csv"))

# ------------------------------------------------
# PCA LANDSCAPE
# ------------------------------------------------

pca_fit <- prcomp(feature_matrix, center = TRUE, scale. = TRUE)
state_dat$PC1 <- pca_fit$x[, 1]
state_dat$PC2 <- pca_fit$x[, 2]

p_pca <- state_dat %>%
  ggplot(aes(PC1, PC2, colour = State)) +
  geom_point(alpha = 0.4, size = 0.8) +
  facet_grid(Group ~ Phase) +
  labs(
    title = "Behavioral state-space landscape",
    subtitle = paste0("Bin level: ", bin_level, "; proximity input: ", proximity_col),
    x = "PC1",
    y = "PC2"
  ) +
  make_nature_theme()

save_plot_svg_pdf(p_pca, file.path(output_dir, "figures", "behavioral_state_space_pca"), width = 180, height = 140)

# ------------------------------------------------
# STATE OCCUPANCY PLOT
# ------------------------------------------------

p_occ <- occupancy_tbl %>%
  ggplot(aes(State, frac_time, fill = Group)) +
  geom_violin(alpha = 0.5, linewidth = 0.2, trim = FALSE) +
  geom_jitter(width = 0.1, size = 0.8, alpha = 0.7) +
  facet_grid(Phase ~ State) +
  labs(
    title = "Behavioral state occupancy",
    subtitle = paste0("Bin level: ", bin_level),
    y = "Fraction of time",
    x = NULL
  ) +
  make_nature_theme()

save_plot_svg_pdf(p_occ, file.path(output_dir, "figures", "state_occupancy"), width = 180, height = 120)

message("Behavioral state-space analysis complete.")
