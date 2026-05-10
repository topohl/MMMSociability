# ================================================================
# Burstiness / Temporal Instability Analysis
# MMMSociability
# ================================================================
# Goal:
#   Quantify whether stress changes the temporal structure of
#   movement / entropy / proximity beyond mean-level effects.
#
# Key outputs:
#   - RMSSD
#   - coefficient of variation
#   - Fano factor
#   - lag-1 autocorrelation
#   - burst-amplitude metrics
#
# Recommended interpretation:
#   Higher RMSSD/Fano/CV with lower ACF1 may indicate more
#   fragmented or unstable behavior.
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(readr)
})

source("Functions/behavioral_dynamics_helpers.R")

# ------------------------------------------------
# USER INPUT
# ------------------------------------------------

input_file <- "analysis_ready/03_derived_metrics/phase_based/all_behavior_metrics.csv"
output_dir <- "analysis_ready/06_behavioral_dynamics/burstiness"

# ------------------------------------------------
# LOAD + STANDARDIZE
# ------------------------------------------------

raw_dat <- read_behavior_table(input_file)
behav <- standardize_behavior_columns(raw_dat)

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
  summarise(calc_instability_metrics(Value), .groups = "drop")

write_table(
  instability_tbl,
  file.path(output_dir, "tables", "burstiness_metrics_per_animal.csv")
)

# ------------------------------------------------
# GROUP SUMMARY
# ------------------------------------------------

group_summary <- instability_tbl %>%
  group_by(Group, Sex, Phase, CageChange, Metric) %>%
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

write_table(
  group_summary,
  file.path(output_dir, "tables", "burstiness_group_summary.csv")
)

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
    rolling_sd = zoo::rollapply(Value, width = 3, FUN = sd,
                                align = "right", fill = NA, na.rm = TRUE),
    rolling_rmssd = zoo::rollapply(
      Value,
      width = 4,
      FUN = function(x) calc_rmssd(x),
      align = "right",
      fill = NA
    )
  ) %>%
  ungroup()

write_table(
  rolling_tbl,
  file.path(output_dir, "tables", "rolling_instability_metrics.csv")
)

p_roll <- rolling_tbl %>%
  group_by(Group, Metric, Phase, TimeIndex) %>%
  summarise(
    mean_rmssd = mean(rolling_rmssd, na.rm = TRUE),
    sem_rmssd = sd(rolling_rmssd, na.rm = TRUE) / sqrt(sum(is.finite(rolling_rmssd))),
    .groups = "drop"
  ) %>%
  ggplot(aes(TimeIndex, mean_rmssd, colour = Group, fill = Group)) +
  geom_line(linewidth = 0.5) +
  geom_ribbon(aes(ymin = mean_rmssd - sem_rmssd,
                  ymax = mean_rmssd + sem_rmssd),
              alpha = 0.2,
              linewidth = 0) +
  facet_grid(Metric ~ Phase, scales = "free_y") +
  labs(
    title = "Temporal evolution of behavioral instability",
    y = "Rolling RMSSD",
    x = "Time"
  ) +
  make_nature_theme()

save_plot_svg_pdf(
  p_roll,
  file.path(output_dir, "figures", "rolling_rmssd_timecourse"),
  width = 180,
  height = 140
)

message("Burstiness analysis complete.")
