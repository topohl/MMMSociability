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
#   - group-blind: Group/Sex are used only after state inference for summaries
#
# Input expectation:
#   Run Analysis/03_build_multiscale_behavior_metrics.R first.
#
# Recommended scale:
#   5–10 min bins. Phase-level data are too coarse for HMMs.
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

source("C:/Users/topohl/Documents/GitHub/MMMSociability/Functions/behavioral_dynamics_helpers.R")

# ------------------------------------------------
# USER INPUT
# ------------------------------------------------

bin_level <- "10min_based"
input_file <- file.path("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/analysis_ready/03_derived_metrics", bin_level, "all_behavior_metrics.csv")
output_dir <- file.path("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/analysis_ready/06_behavioral_dynamics/hmm_states", bin_level)
n_states <- 4

# Use normalized proximity because raw contact seconds scale with bin size.
proximity_col <- "ProximityFraction"

ensure_dir(output_dir)
ensure_dir(file.path(output_dir, "tables"))
ensure_dir(file.path(output_dir, "figures"))

if (!requireNamespace("depmixS4", quietly = TRUE)) {
  stop("Please install depmixS4 to run HMM behavioral states.")
}

raw_dat <- read_behavior_table(input_file)
if (!proximity_col %in% names(raw_dat)) proximity_col <- "Proximity"
behav <- standardize_behavior_columns(raw_dat, proximity_col = proximity_col)

hmm_dat <- behav %>%
  mutate(
    Movement_z = z_within_metric(Movement),
    Entropy_z = z_within_metric(Entropy),
    Proximity_z = z_within_metric(Proximity),
    BinLevel = bin_level,
    ProximityInput = proximity_col
  ) %>%
  filter(is.finite(Movement_z), is.finite(Entropy_z), is.finite(Proximity_z)) %>%
  arrange(AnimalNum, CageChange, Phase, TimeIndex) %>%
  mutate(SequenceID = interaction(AnimalNum, CageChange, Phase, drop = TRUE, sep = "__"))

sequence_tbl <- hmm_dat %>%
  count(SequenceID, AnimalNum, Group, Sex, CageChange, Phase, name = "n_bins") %>%
  arrange(AnimalNum, CageChange, Phase) %>%
  filter(n_bins >= 2)

hmm_dat <- hmm_dat %>%
  semi_join(sequence_tbl %>% select(SequenceID), by = "SequenceID") %>%
  mutate(SequenceID = factor(SequenceID, levels = sequence_tbl$SequenceID)) %>%
  arrange(SequenceID, TimeIndex)

ntimes <- sequence_tbl$n_bins

if (nrow(hmm_dat) < n_states * 10) {
  stop("Too few usable rows for HMM. Use a finer bin level or fewer states.", call. = FALSE)
}

set.seed(1)
mod <- depmixS4::depmix(
  list(
    Movement_z ~ 1,
    Entropy_z ~ 1,
    Proximity_z ~ 1
  ),
  data = hmm_dat,
  ntimes = ntimes,
  nstates = n_states,
  family = list(gaussian(), gaussian(), gaussian())
)

fit_warnings <- character()
fit_mod <- withCallingHandlers(
  depmixS4::fit(mod, verbose = FALSE),
  warning = function(w) {
    fit_warnings <<- c(fit_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
post <- depmixS4::posterior(fit_mod, type = "viterbi")
hmm_dat$State <- factor(post$state)
fit_loglik <- tryCatch(as.numeric(fit_mod@logLike), error = function(e) NA_real_)

hmm_qc <- tibble(
  BinLevel = bin_level,
  ProximityInput = proximity_col,
  GroupBlind = TRUE,
  n_states = n_states,
  n_rows = nrow(hmm_dat),
  n_animals = n_distinct(hmm_dat$AnimalNum),
  n_sequences = length(ntimes),
  min_sequence_bins = min(ntimes),
  median_sequence_bins = median(ntimes),
  max_sequence_bins = max(ntimes),
  logLik = fit_loglik,
  fit_warning = ifelse(length(fit_warnings) == 0, NA_character_, paste(unique(fit_warnings), collapse = " | "))
)

write_table(sequence_tbl, file.path(output_dir, "tables", "hmm_sequence_design.csv"))
write_table(hmm_qc, file.path(output_dir, "tables", "hmm_model_qc.csv"))
write_table(
  hmm_dat %>%
    select(BinLevel, ProximityInput, SequenceID, AnimalNum, Group, Sex, Phase, CageChange, TimeIndex, Movement_z, Entropy_z, Proximity_z, State),
  file.path(output_dir, "tables", "hmm_state_assignments.csv")
)

state_summary <- hmm_dat %>%
  group_by(BinLevel, ProximityInput, State) %>%
  summarise(
    Movement_z = mean(Movement_z, na.rm = TRUE),
    Entropy_z = mean(Entropy_z, na.rm = TRUE),
    Proximity_z = mean(Proximity_z, na.rm = TRUE),
    n_bins = n(),
    .groups = "drop"
  )

write_table(state_summary, file.path(output_dir, "tables", "hmm_state_summary.csv"))

transition_tbl <- hmm_dat %>%
  group_by(BinLevel, ProximityInput, Group, Sex, Phase, CageChange, AnimalNum) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  mutate(NextState = lead(State)) %>%
  filter(!is.na(NextState)) %>%
  count(BinLevel, ProximityInput, Group, Sex, Phase, CageChange, AnimalNum, State, NextState, name = "Transitions")

write_table(transition_tbl, file.path(output_dir, "tables", "hmm_transition_counts.csv"))

transition_prob_tbl <- transition_tbl %>%
  group_by(BinLevel, ProximityInput, Group, Sex, Phase, CageChange, AnimalNum, State) %>%
  mutate(TransitionProbability = Transitions / sum(Transitions)) %>%
  ungroup()

write_table(transition_prob_tbl, file.path(output_dir, "tables", "hmm_transition_probabilities.csv"))

dwell_tbl <- hmm_dat %>%
  group_by(BinLevel, ProximityInput, Group, Sex, Phase, CageChange, AnimalNum) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  mutate(
    StateRun = cumsum(State != lag(State, default = first(State))) + 1L
  ) %>%
  group_by(BinLevel, ProximityInput, Group, Sex, Phase, CageChange, AnimalNum, State, StateRun) %>%
  summarise(dwell_bins = n(), .groups = "drop") %>%
  group_by(BinLevel, ProximityInput, Group, Sex, Phase, CageChange, AnimalNum, State) %>%
  summarise(
    mean_dwell_bins = mean(dwell_bins, na.rm = TRUE),
    median_dwell_bins = median(dwell_bins, na.rm = TRUE),
    max_dwell_bins = max(dwell_bins, na.rm = TRUE),
    n_bouts = n(),
    .groups = "drop"
  )

write_table(dwell_tbl, file.path(output_dir, "tables", "hmm_state_dwell_times.csv"))

occupancy_tbl <- hmm_dat %>%
  count(BinLevel, ProximityInput, Group, Sex, Phase, CageChange, AnimalNum, State) %>%
  group_by(BinLevel, ProximityInput, Group, Sex, Phase, CageChange, AnimalNum) %>%
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
    subtitle = paste0("Bin level: ", bin_level, "; proximity input: ", proximity_col),
    y = "Fraction of time",
    x = NULL
  ) +
  make_nature_theme()

save_plot_svg_pdf(p_occ, file.path(output_dir, "figures", "hmm_state_occupancy"), width = 180, height = 120)

message("HMM behavioral-state analysis complete.")
