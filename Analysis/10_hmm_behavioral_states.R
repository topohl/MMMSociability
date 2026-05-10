# ================================================================
# Hidden Markov Behavioral States
# MMMSociability
# ================================================================
# Goal:
#   Infer latent behavioral states and transitions using HMMs.
#
# Compared with k-means:
#   - explicitly models temporal persistence
#   - estimates transition probabilities
#   - estimates dwell times
#   - more biologically interpretable dynamics
#
# Requires:
#   depmixS4
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
})

source("Functions/behavioral_dynamics_helpers.R")

input_file <- "analysis_ready/03_derived_metrics/phase_based/all_behavior_metrics.csv"
output_dir <- "analysis_ready/06_behavioral_dynamics/hmm_states"
n_states <- 4

ensure_dir(output_dir)
ensure_dir(file.path(output_dir, "tables"))
ensure_dir(file.path(output_dir, "figures"))

if (!requireNamespace("depmixS4", quietly = TRUE)) {
  stop("Please install depmixS4 to run HMM behavioral states.")
}

raw_dat <- read_behavior_table(input_file)
behav <- standardize_behavior_columns(raw_dat)

hmm_dat <- behav %>%
  mutate(
    Movement_z = z_within_metric(Movement),
    Entropy_z = z_within_metric(Entropy),
    Proximity_z = z_within_metric(Proximity)
  ) %>%
  arrange(AnimalNum, TimeIndex)

mod <- depmixS4::depmix(
  list(
    Movement_z ~ 1,
    Entropy_z ~ 1,
    Proximity_z ~ 1
  ),
  data = hmm_dat,
  nstates = n_states,
  family = list(gaussian(), gaussian(), gaussian())
)

fit_mod <- depmixS4::fit(mod, verbose = FALSE)

post <- depmixS4::posterior(fit_mod)

hmm_dat$State <- factor(post$state)

state_summary <- hmm_dat %>%
  group_by(State) %>%
  summarise(
    Movement_z = mean(Movement_z, na.rm = TRUE),
    Entropy_z = mean(Entropy_z, na.rm = TRUE),
    Proximity_z = mean(Proximity_z, na.rm = TRUE),
    n_bins = n(),
    .groups = "drop"
  )

write_table(state_summary, file.path(output_dir, "tables", "hmm_state_summary.csv"))

transition_tbl <- hmm_dat %>%
  group_by(Group, Sex, Phase, CageChange, AnimalNum) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  mutate(NextState = lead(State)) %>%
  filter(!is.na(NextState)) %>%
  count(State, NextState, Group, Phase, name = "Transitions")

write_table(transition_tbl, file.path(output_dir, "tables", "hmm_transition_counts.csv"))

occupancy_tbl <- hmm_dat %>%
  count(Group, Sex, Phase, CageChange, AnimalNum, State) %>%
  group_by(Group, Sex, Phase, CageChange, AnimalNum) %>%
  mutate(frac_time = n / sum(n)) %>%
  ungroup()

write_table(occupancy_tbl, file.path(output_dir, "tables", "hmm_state_occupancy.csv"))

p_occ <- occupancy_tbl %>%
  ggplot(aes(State, frac_time, fill = Group)) +
  geom_violin(alpha = 0.5, linewidth = 0.2, trim = FALSE) +
  geom_jitter(width = 0.1, size = 0.8, alpha = 0.7) +
  facet_grid(Phase ~ State, scales = "free_y") +
  labs(
    title = "Hidden Markov behavioral state occupancy",
    y = "Fraction of time",
    x = NULL
  ) +
  make_nature_theme()

save_plot_svg_pdf(p_occ, file.path(output_dir, "figures", "hmm_state_occupancy"), width = 180, height = 120)

message("HMM behavioral-state analysis complete.")
