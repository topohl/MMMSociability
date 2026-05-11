# ================================================================
# Burstiness / Temporal Instability Analysis
# MMMSociability
# ================================================================
# Goal:
#   Quantify whether stress changes the temporal structure of
#   movement / entropy / proximity beyond mean-level effects.
#
# Input expectation:
#   Run Analysis/03_build_multiscale_behavior_metrics.R first.
#
# Recommended scale:
#   1–5 min bins. Phase-level data are too coarse for RMSSD/ACF/burstiness.
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(zoo)
})

source("Functions/behavioral_dynamics_helpers.R")

# ------------------------------------------------
# USER INPUT
# ------------------------------------------------

bin_level <- "5min_based"
input_file <- file.path("analysis_ready/03_derived_metrics", bin_level, "all_behavior_metrics.csv")
output_dir <- file.path("analysis_ready/06_behavioral_dynamics/burstiness", bin_level)

# Prefer normalized proximity for temporal instability because raw contact
# seconds scale with bin size. Falls back to Proximity if ProximityFraction is absent.
proximity_col <- "ProximityFraction"

# ------------------------------------------------
# LOAD + STANDARDIZE
# ------------------------------------------------

raw_dat <- read_behavior_table(input_file)
if (!proximity_col %in% names(raw_dat)) proximity_col <- "Proximity"

behav <- standardize_behavior_columns(raw_dat, proximity_col = proximity_col)

ensure_dir(output_dir)
ensure_dir(file.path(output_dir, "tables"))
ensure_dir(file.path(output_dir, "figures"))

# ------------------------------------------------
# LONG FORMAT
# ------------------------------------------------

long_dat <- behav %>%
  pivot_longer(
    cols = c(Movement, Entropy, Proximity),
    names_to = "Metric",
    values_to = "Value"
  )

# ------------------------------------------------
# INSTABILITY METRICS
# ------------------------------------------------

instability_tbl <- long_dat %>%
  group_by(Group, Sex, Phase, CageChange, AnimalNum, Metric) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  summarise(calc_instability_metrics(Value), .groups = "drop") %>%
  mutate(BinLevel = bin_level, ProximityInput = proximity_col)

write_table(instability_tbl, file.path(output_dir, "tables", "burstiness_metrics_per_animal.csv"))

# ------------------------------------------------
# GROUP SUMMARY
# ------------------------------------------------

group_summary <- instability_tbl %>%
  group_by(BinLevel, ProximityInput, Group, Sex, Phase, CageChange, Metric) %>%
  summarise(
    across(
      c(mean, cv, fano, rmssd, acf1),
      list(mean = ~mean(.x, na.rm = TRUE),
           sem = ~sd(.x, na.rm = TRUE) / sqrt(sum(is.finite(.x)))),
      .names = "{.col}_{.fn}"
    ),
    n_animals = n_distinct(AnimalNum),
    .groups = "drop"
  )

write_table(group_summary, file.path(output_dir, "tables", "burstiness_group_summary.csv"))

# ------------------------------------------------
# PLOTS
# ------------------------------------------------

plot_metrics <- c("rmssd", "cv", "fano", "acf1")

for (metric_name in plot_metrics) {
  p <- instability_tbl %>%
    ggplot(aes(x = Group, y = .data[[metric_name]], fill = Group)) +
    geom_violin(alpha = 0.5, linewidth = 0.2, trim = FALSE) +
    geom_jitter(width = 0.1, size = 0.8, alpha = 0.7) +
    facet_grid(Metric ~ Phase, scales = "free_y") +
    labs(
      title = paste0(toupper(metric_name), " behavioral instability"),
      subtitle = paste0("Bin level: ", bin_level, "; proximity input: ", proximity_col),
      y = metric_name,
      x = NULL
    ) +
    make_nature_theme()

  save_plot_svg_pdf(
    p,
    file.path(output_dir, "figures", paste0(metric_name, "_group_comparison")),
    width = 180,
    height = 120
  )
}

# ------------------------------------------------
# TIME-RESOLVED VOLATILITY
# ------------------------------------------------

rolling_tbl <- long_dat %>%
  group_by(Group, Sex, Phase, CageChange, AnimalNum, Metric) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  mutate(
    rolling_sd = zoo::rollapply(Value, width = 3, FUN = sd, align = "right", fill = NA, na.rm = TRUE),
    rolling_rmssd = zoo::rollapply(Value, width = 4, FUN = calc_rmssd, align = "right", fill = NA)
  ) %>%
  ungroup() %>%
  mutate(BinLevel = bin_level, ProximityInput = proximity_col)

write_table(rolling_tbl, file.path(output_dir, "tables", "rolling_instability_metrics.csv"))

p_roll <- rolling_tbl %>%
  group_by(Group, Metric, Phase, TimeIndex) %>%
  summarise(
    mean_rmssd = mean(rolling_rmssd, na.rm = TRUE),
    sem_rmssd = sd(rolling_rmssd, na.rm = TRUE) / sqrt(sum(is.finite(rolling_rmssd))),
    .groups = "drop"
  ) %>%
  ggplot(aes(TimeIndex, mean_rmssd, colour = Group, fill = Group)) +
  geom_line(linewidth = 0.5) +
  geom_ribbon(aes(ymin = mean_rmssd - sem_rmssd, ymax = mean_rmssd + sem_rmssd), alpha = 0.2, linewidth = 0) +
  facet_grid(Metric ~ Phase, scales = "free_y") +
  labs(
    title = "Temporal evolution of behavioral instability",
    subtitle = paste0("Bin level: ", bin_level),
    y = "Rolling RMSSD",
    x = "Time bin"
  ) +
  make_nature_theme()

save_plot_svg_pdf(p_roll, file.path(output_dir, "figures", "rolling_rmssd_timecourse"), width = 180, height = 140)

message("Burstiness analysis complete.")
