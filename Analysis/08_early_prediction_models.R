# ================================================================
# Early Behavioral Prediction Models
# MMMSociability
# ================================================================
# Goal:
#   Test whether early movement / entropy / proximity dynamics predict
#   later stress burden or phenotype labels.
#
# Input expectation:
#   Run Analysis/03_build_multiscale_behavior_metrics.R first.
#
# Recommended scale:
#   10 min default. This balances temporal resolution with noise.
#   Use 5min_based as sensitivity analysis.
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

source("Functions/behavioral_dynamics_helpers.R")

# ------------------------------------------------
# USER INPUT
# ------------------------------------------------

bin_level <- "10min_based"
input_file <- file.path("analysis_ready/03_derived_metrics", bin_level, "all_behavior_metrics.csv")
output_dir <- file.path("analysis_ready/06_behavioral_dynamics/early_prediction", bin_level)

# Optional endpoint file. If NULL or missing, the script will try to use endpoints
# already present in input_file.
endpoint_file <- NULL

# Endpoint column to predict. Change this to your real stress score column.
outcome_col <- "stress_z_score"

# Early window definition. Default: first active-phase bins per animal/cage-change.
early_phase_pattern <- "active|dark|night"
max_early_bins_per_animal <- 4

# Use normalized proximity when available.
animal_col <- NULL
time_col <- NULL
group_col <- NULL
sex_col <- NULL
phase_col <- NULL
cage_col <- NULL
movement_col <- NULL
entropy_col <- NULL
proximity_col <- "ProximityFraction"

# ------------------------------------------------
# LOAD + STANDARDIZE
# ------------------------------------------------

raw_dat <- read_behavior_table(input_file)
if (!proximity_col %in% names(raw_dat)) proximity_col <- "Proximity"

behav <- standardize_behavior_columns(
  raw_dat,
  animal_col = animal_col,
  time_col = time_col,
  group_col = group_col,
  sex_col = sex_col,
  phase_col = phase_col,
  cage_col = cage_col,
  movement_col = movement_col,
  entropy_col = entropy_col,
  proximity_col = proximity_col
)

ensure_dir(output_dir)
ensure_dir(file.path(output_dir, "tables"))
ensure_dir(file.path(output_dir, "figures"))

# ------------------------------------------------
# EARLY WINDOW SELECTION
# ------------------------------------------------

has_active_phase <- any(str_detect(str_to_lower(as.character(behav$Phase)), early_phase_pattern))

early_dat <- if (has_active_phase) {
  behav %>% filter(str_detect(str_to_lower(as.character(Phase)), early_phase_pattern))
} else {
  behav
}

early_dat <- early_dat %>%
  group_by(AnimalNum, CageChange, Phase) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  mutate(early_rank = row_number()) %>%
  filter(early_rank <= max_early_bins_per_animal) %>%
  ungroup() %>%
  mutate(BinLevel = bin_level, ProximityInput = proximity_col)

write_table(early_dat, file.path(output_dir, "tables", "early_window_rows_used.csv"))

# ------------------------------------------------
# FEATURE EXTRACTION
# ------------------------------------------------

metric_long <- early_dat %>%
  pivot_longer(cols = c(Movement, Entropy, Proximity), names_to = "Metric", values_to = "Value")

feature_long <- metric_long %>%
  group_by(AnimalNum, Group, Sex, Metric) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  summarise(calc_instability_metrics(Value), .groups = "drop")

feature_wide <- feature_long %>%
  pivot_wider(
    id_cols = c(AnimalNum, Group, Sex),
    names_from = Metric,
    values_from = c(mean, sd, cv, fano, rmssd, acf1, p95, max),
    names_glue = "{Metric}_{.value}"
  ) %>%
  mutate(
    BinLevel = bin_level,
    ProximityInput = proximity_col,
    SocialWithdrawal_mean = -scale(Proximity_mean)[, 1] + scale(Movement_mean)[, 1],
    PassiveIsolation_mean = -scale(Proximity_mean)[, 1] - scale(Movement_mean)[, 1],
    SocialEngagement_mean = scale(Proximity_mean)[, 1] + scale(Movement_mean)[, 1]
  )

write_table(feature_wide, file.path(output_dir, "tables", "early_behavior_features.csv"))

# ------------------------------------------------
# ENDPOINT HANDLING
# ------------------------------------------------

endpoint_dat <- NULL

if (!is.null(endpoint_file) && file.exists(endpoint_file)) {
  endpoint_raw <- read_behavior_table(endpoint_file)
  endpoint_animal_col <- first_existing_col(endpoint_raw, c("AnimalNum", "Animal", "MouseID", "Mouse", "ID", "RFID", "animal_id"), TRUE, "endpoint animal column")
  endpoint_dat <- endpoint_raw %>%
    transmute(AnimalNum = .data[[endpoint_animal_col]], outcome = suppressWarnings(as.numeric(.data[[outcome_col]])))
} else if (outcome_col %in% names(raw_dat)) {
  endpoint_animal_col <- first_existing_col(raw_dat, c("AnimalNum", "Animal", "MouseID", "Mouse", "ID", "RFID", "animal_id"), TRUE, "endpoint animal column")
  endpoint_dat <- raw_dat %>%
    group_by(AnimalNum = .data[[endpoint_animal_col]]) %>%
    summarise(outcome = first(na.omit(suppressWarnings(as.numeric(.data[[outcome_col]])))), .groups = "drop")
}

if (is.null(endpoint_dat)) {
  message("No endpoint column found. Feature extraction finished, but predictive modeling skipped.")
  message("Set outcome_col and/or endpoint_file in Analysis/08_early_prediction_models.R.")
  quit(save = "no", status = 0)
}

model_dat <- feature_wide %>%
  left_join(endpoint_dat, by = "AnimalNum") %>%
  filter(is.finite(outcome))

write_table(model_dat, file.path(output_dir, "tables", "early_prediction_model_input.csv"))

if (nrow(model_dat) < 8) {
  warning("Fewer than 8 animals with outcome data. Modeling skipped; feature and input tables were written.")
  quit(save = "no", status = 0)
}

# ------------------------------------------------
# UNIVARIATE ASSOCIATIONS
# ------------------------------------------------

feature_cols <- model_dat %>%
  select(where(is.numeric)) %>%
  select(-outcome) %>%
  names()

association_tbl <- map_dfr(feature_cols, function(fc) {
  tibble(
    feature = fc,
    spearman_rho = safe_cor(model_dat[[fc]], model_dat$outcome, method = "spearman"),
    pearson_r = safe_cor(model_dat[[fc]], model_dat$outcome, method = "pearson"),
    n = sum(is.finite(model_dat[[fc]]) & is.finite(model_dat$outcome))
  )
}) %>%
  arrange(desc(abs(spearman_rho))) %>%
  mutate(BinLevel = bin_level, ProximityInput = proximity_col)

write_table(association_tbl, file.path(output_dir, "tables", "early_feature_outcome_associations.csv"))

# ------------------------------------------------
# LEAVE-ONE-ANIMAL-OUT RIDGE/ELASTIC-NET IF glmnet AVAILABLE
# ------------------------------------------------

if (requireNamespace("glmnet", quietly = TRUE)) {
  x <- model_dat %>%
    select(all_of(feature_cols)) %>%
    mutate(across(everything(), ~replace_na(.x, median(.x, na.rm = TRUE)))) %>%
    as.matrix()

  y <- model_dat$outcome
  loo_pred <- rep(NA_real_, length(y))

  set.seed(123)
  for (i in seq_along(y)) {
    train_idx <- setdiff(seq_along(y), i)
    fit <- glmnet::cv.glmnet(
      x[train_idx, , drop = FALSE],
      y[train_idx],
      alpha = 0.5,
      standardize = TRUE,
      nfolds = min(5, length(train_idx))
    )
    loo_pred[i] <- as.numeric(predict(fit, newx = x[i, , drop = FALSE], s = "lambda.min"))
  }

  pred_tbl <- model_dat %>%
    transmute(AnimalNum, Group, Sex, BinLevel, ProximityInput, observed = outcome, predicted = loo_pred, residual = observed - predicted)

  perf_tbl <- tibble(
    BinLevel = bin_level,
    ProximityInput = proximity_col,
    n = nrow(pred_tbl),
    loo_pearson_r = safe_cor(pred_tbl$observed, pred_tbl$predicted, "pearson"),
    loo_spearman_rho = safe_cor(pred_tbl$observed, pred_tbl$predicted, "spearman"),
    loo_rmse = sqrt(mean((pred_tbl$observed - pred_tbl$predicted)^2, na.rm = TRUE)),
    baseline_rmse_mean_only = sqrt(mean((pred_tbl$observed - mean(pred_tbl$observed, na.rm = TRUE))^2, na.rm = TRUE))
  ) %>%
    mutate(relative_rmse_vs_baseline = loo_rmse / baseline_rmse_mean_only)

  write_table(pred_tbl, file.path(output_dir, "tables", "leave_one_animal_out_predictions.csv"))
  write_table(perf_tbl, file.path(output_dir, "tables", "leave_one_animal_out_performance.csv"))

  p_pred <- pred_tbl %>%
    ggplot(aes(observed, predicted, colour = Group)) +
    geom_point(size = 1.6, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.5) +
    labs(
      title = "Early behavior predicts later endpoint",
      subtitle = paste0("LOO elastic-net; bin level: ", bin_level),
      x = paste0("Observed ", outcome_col),
      y = paste0("Predicted ", outcome_col)
    ) +
    make_nature_theme()

  save_plot_svg_pdf(p_pred, file.path(output_dir, "figures", "early_prediction_observed_vs_predicted"), width = 85, height = 75)
} else {
  message("Package glmnet not installed. Wrote feature associations, but skipped cross-validated elastic-net modeling.")
}

# ------------------------------------------------
# TOP ASSOCIATION PLOT
# ------------------------------------------------

p_assoc <- association_tbl %>%
  slice_head(n = 12) %>%
  mutate(feature = reorder(feature, abs(spearman_rho))) %>%
  ggplot(aes(feature, spearman_rho)) +
  geom_col(width = 0.7) +
  coord_flip() +
  labs(
    title = "Top early behavioral correlates",
    subtitle = paste0("Bin level: ", bin_level),
    y = "Spearman rho with endpoint",
    x = NULL
  ) +
  make_nature_theme()

save_plot_svg_pdf(p_assoc, file.path(output_dir, "figures", "top_early_feature_associations"), width = 90, height = 80)

message("Early prediction analysis complete.")
