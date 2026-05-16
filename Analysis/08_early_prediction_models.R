# ================================================================
# Early Behavioral Prediction Models
# MMMSociability
# ================================================================
# Goal:
#   Test whether early movement / entropy / proximity dynamics predict
#   later stress burden or phenotype labels, with explicit feature QC,
#   corrected univariate screening, and cross-validated prediction.
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

source("C:/Users/topohl/Documents/GitHub/MMMSociability/Functions/behavioral_dynamics_helpers.R")
source("C:/Users/topohl/Documents/GitHub/MMMSociability/Functions/behavioral_dynamics_stats_helpers.R")
source("C:/Users/topohl/Documents/GitHub/MMMSociability/Functions/duration_normalization_helpers.R")

# ------------------------------------------------
# USER INPUT
# ------------------------------------------------

bin_level <- "10min_based"
input_file <- file.path("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/analysis_ready/03_derived_metrics", bin_level, "all_behavior_metrics.csv")
output_dir <- file.path("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/analysis_ready/06_behavioral_dynamics/early_prediction", bin_level)

# Optional endpoint file. If NULL or missing, the script will try to use endpoints
# already present in input_file.
endpoint_file <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/SIS_Analysis/E9_Behavior_Data.xlsx"
endpoint_sheet <- "zScore"
endpoint_cols <- c(
  "CombZ", "stress_z_score", "NOR", "sucrose_pref", "weight_dev", "delta_cort", "adrenal_weight", "spleen_weight",
  "SucrosePreference", "Corticosterone", "DeltaCorticosterone"
)
primary_outcome <- "CombZ"
primary_outcome_label <- "CombZ (lower = worse depressive-like endpoint; higher = more resilient-like endpoint)"
primary_outcome_direction <- "lower_worse"

# Endpoint column to predict. Uses primary_outcome by default.
outcome_col <- primary_outcome

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

group_colors <- c(
  "CON" = "#3d3b6e",
  "SUS" = "#e63947",
  "RES" = "#C6C3BB",
  "All" = "grey55"
)
group_levels <- c("CON", "RES", "SUS", "All")
metric_levels <- c("Movement", "Entropy", "Proximity", "SocialWithdrawal", "PassiveIsolation", "SocialEngagement")
measure_levels <- c("mean", "sd", "cv", "fano", "rmssd", "acf1", "p95", "max")

# Cheap permutation test for the observed-vs-predicted association. This keeps
# training fixed and asks whether prediction alignment exceeds random pairing.
n_prediction_permutations <- 5000

# ------------------------------------------------
# HELPERS
# ------------------------------------------------

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

parse_feature_name <- function(feature) {
  tibble(feature = feature) %>%
    tidyr::extract(
      feature,
      into = c("Metric", "Measure"),
      regex = "^(Movement|Entropy|Proximity|SocialWithdrawal|PassiveIsolation|SocialEngagement)_(.*)$",
      remove = FALSE
    ) %>%
    mutate(
      Metric = factor(Metric, levels = metric_levels),
      Measure = factor(Measure, levels = measure_levels),
      FeatureLabel = case_when(
        is.na(Metric) ~ feature,
        TRUE ~ paste(as.character(Metric), as.character(Measure), sep = " ")
      )
    )
}

safe_scale <- function(x) {
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

safe_cor_test <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4 || sd(x[ok]) == 0 || sd(y[ok]) == 0) {
    return(tibble(
      n = sum(ok), estimate = NA_real_, statistic = NA_real_, p.value = NA_real_,
      conf.low = NA_real_, conf.high = NA_real_
    ))
  }
  test <- try(stats::cor.test(x[ok], y[ok], method = method, exact = FALSE), silent = TRUE)
  if (inherits(test, "try-error")) {
    return(tibble(
      n = sum(ok), estimate = NA_real_, statistic = NA_real_, p.value = NA_real_,
      conf.low = NA_real_, conf.high = NA_real_
    ))
  }
  ci <- if (!is.null(test$conf.int)) unname(test$conf.int) else c(NA_real_, NA_real_)
  tibble(
    n = sum(ok),
    estimate = unname(test$estimate),
    statistic = unname(test$statistic),
    p.value = test$p.value,
    conf.low = ci[1],
    conf.high = ci[2]
  )
}

make_feature_qc <- function(dat, feature_cols) {
  map_dfr(feature_cols, function(fc) {
    x <- dat[[fc]]
    tibble(
      feature = fc,
      n_animals = nrow(dat),
      n_finite = sum(is.finite(x)),
      missing_fraction = mean(!is.finite(x)),
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      q25 = quantile(x, 0.25, na.rm = TRUE, names = FALSE),
      q75 = quantile(x, 0.75, na.rm = TRUE, names = FALSE),
      min = min(x, na.rm = TRUE),
      max = max(x, na.rm = TRUE),
      zero_variance = sd(x, na.rm = TRUE) == 0
    )
  }) %>%
    left_join(parse_feature_name(feature_cols), by = "feature") %>%
    arrange(Metric, Measure)
}

impute_feature_matrix <- function(dat, feature_cols) {
  dat %>%
    select(all_of(feature_cols)) %>%
    mutate(across(everything(), ~{
      med <- median(.x, na.rm = TRUE)
      if (!is.finite(med)) med <- 0
      replace_na(.x, med)
    })) %>%
    as.matrix()
}

standardize_matrix <- function(x) {
  x_center <- colMeans(x, na.rm = TRUE)
  x_scale <- apply(x, 2, sd, na.rm = TRUE)
  x_scale[!is.finite(x_scale) | x_scale == 0] <- 1
  sweep(sweep(x, 2, x_center, "-"), 2, x_scale, "/")
}

prediction_permutation_p <- function(observed, predicted, n_perm = 5000, seed = 123) {
  ok <- is.finite(observed) & is.finite(predicted)
  if (sum(ok) < 4 || sd(observed[ok]) == 0 || sd(predicted[ok]) == 0) {
    return(tibble(observed_r = NA_real_, permutation_p = NA_real_, n_permutations = n_perm))
  }
  observed_r <- suppressWarnings(cor(observed[ok], predicted[ok], method = "pearson"))
  set.seed(seed)
  null_r <- replicate(n_perm, suppressWarnings(cor(observed[ok], sample(predicted[ok]), method = "pearson")))
  tibble(
    observed_r = observed_r,
    permutation_p = (sum(abs(null_r) >= abs(observed_r), na.rm = TRUE) + 1) / (sum(is.finite(null_r)) + 1),
    n_permutations = n_perm
  )
}

make_prediction_summary <- function(pred_tbl, outcome_col, n_perm = 5000) {
  baseline_pred <- mean(pred_tbl$observed, na.rm = TRUE)
  baseline_rmse <- sqrt(mean((pred_tbl$observed - baseline_pred)^2, na.rm = TRUE))
  loo_rmse <- sqrt(mean((pred_tbl$observed - pred_tbl$predicted)^2, na.rm = TRUE))
  perm_tbl <- prediction_permutation_p(pred_tbl$observed, pred_tbl$predicted, n_perm = n_perm)

  tibble(
    Outcome = outcome_col,
    n = nrow(pred_tbl),
    loo_pearson_r = safe_cor(pred_tbl$observed, pred_tbl$predicted, "pearson"),
    loo_spearman_rho = safe_cor(pred_tbl$observed, pred_tbl$predicted, "spearman"),
    loo_rmse = loo_rmse,
    loo_mae = mean(abs(pred_tbl$observed - pred_tbl$predicted), na.rm = TRUE),
    baseline_rmse_mean_only = baseline_rmse,
    relative_rmse_vs_baseline = loo_rmse / baseline_rmse,
    cross_validated_r2_vs_mean = 1 - sum((pred_tbl$observed - pred_tbl$predicted)^2, na.rm = TRUE) /
      sum((pred_tbl$observed - baseline_pred)^2, na.rm = TRUE)
  ) %>%
    bind_cols(perm_tbl %>% select(prediction_permutation_p = permutation_p, prediction_permutations = n_permutations))
}

summarise_model_coefficients <- function(coef_tbl) {
  coef_tbl %>%
    filter(feature != "(Intercept)") %>%
    group_by(feature) %>%
    summarise(
      nonzero_frequency = mean(coefficient != 0, na.rm = TRUE),
      median_coefficient = median(coefficient, na.rm = TRUE),
      mean_coefficient = mean(coefficient, na.rm = TRUE),
      q25_coefficient = quantile(coefficient, 0.25, na.rm = TRUE, names = FALSE),
      q75_coefficient = quantile(coefficient, 0.75, na.rm = TRUE, names = FALSE),
      max_abs_coefficient = max(abs(coefficient), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(parse_feature_name(unique(coef_tbl$feature)), by = "feature") %>%
    arrange(desc(nonzero_frequency), desc(max_abs_coefficient))
}

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
) %>%
  mutate(Group = factor(as.character(Group), levels = unique(c(group_levels, sort(unique(as.character(Group)))))))

ensure_dir(output_dir)
ensure_dir(file.path(output_dir, "tables"))
ensure_dir(file.path(output_dir, "stats_tables"))
ensure_dir(file.path(output_dir, "figures"))
output_dirs <- analysis_output_dirs(output_dir)
ensure_dir(file.path(output_dir, "figures", "publication"))
ensure_dir(file.path(output_dir, "figures", "publication", "overview"))
ensure_dir(file.path(output_dir, "figures", "publication", "panels"))
write_output_manifest(
  output_dir,
  script_name = "08_early_prediction_models.R",
  analysis_name = "early prediction models",
  primary_tables = c(
    "tables/early_behavior_features.csv",
    "tables/early_behavior_features_excluding_short_duration.csv",
    "tables/leave_one_animal_out_performance.csv",
    "tables/leave_one_animal_out_performance_duration_sensitivity.csv"
  ),
  primary_figures = c(
    "figures/publication/panels",
    "figures/publication/overview"
  ),
  notes = c("Prediction figures should report whether duration-sensitive rows were excluded or retained.")
)

epoch_duration_qc <- write_epoch_duration_qc(behav, output_dir, metric_source = "08_early_prediction_models", bin_size_sec = infer_bin_size_sec(behav))

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

early_window_summary <- early_dat %>%
  group_by(BinLevel, ProximityInput, AnimalNum, Group, Sex, CageChange, Phase) %>%
  summarise(
    n_early_bins = n(),
    first_time_index = min(TimeIndex, na.rm = TRUE),
    last_time_index = max(TimeIndex, na.rm = TRUE),
    .groups = "drop"
  )

early_design_tbl <- tibble(
  BinLevel = bin_level,
  ProximityInput = proximity_col,
  EarlyPhasePattern = early_phase_pattern,
  ActivePhaseDetected = has_active_phase,
  MaxEarlyBinsPerAnimalCagePhase = max_early_bins_per_animal,
  n_animals = n_distinct(early_dat$AnimalNum),
  n_rows = nrow(early_dat),
  n_cage_changes = n_distinct(early_dat$CageChange),
  n_phases = n_distinct(early_dat$Phase)
)

write_table(early_dat, file.path(output_dir, "tables", "early_window_rows_used.csv"))
write_table(early_window_summary, file.path(output_dir, "tables", "early_window_summary_by_animal.csv"))
write_table(early_design_tbl, file.path(output_dir, "tables", "early_window_design_summary.csv"))
write_table(filter_short_duration_epochs(early_dat, epoch_duration_qc), file.path(output_dir, "tables", "early_window_rows_used_excluding_short_duration.csv"))

# ------------------------------------------------
# FEATURE EXTRACTION
# ------------------------------------------------

metric_long <- early_dat %>%
  pivot_longer(cols = c(Movement, Entropy, Proximity), names_to = "Metric", values_to = "Value")

feature_long <- metric_long %>%
  group_by(AnimalNum, Group, Sex, Metric) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  summarise(calc_instability_metrics(Value), .groups = "drop")

early_duration_by_animal <- early_dat %>%
  join_duration_qc(epoch_duration_qc) %>%
  distinct(AnimalNum, Group, Sex, CageChange, Phase, .keep_all = TRUE) %>%
  group_by(AnimalNum, Group, Sex) %>%
  summarise(
    early_observed_bins = sum(observed_bins, na.rm = TRUE),
    early_observation_hours = sum(total_observation_duration_hours, na.rm = TRUE),
    min_duration_completeness_fraction = min(duration_completeness_fraction, na.rm = TRUE),
    contains_short_duration_epoch = any(short_epoch %in% TRUE | cage_change_duration_class == "short", na.rm = TRUE),
    .groups = "drop"
  )

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
    SocialWithdrawal_mean = -safe_scale(Proximity_mean) + safe_scale(Movement_mean),
    PassiveIsolation_mean = -safe_scale(Proximity_mean) - safe_scale(Movement_mean),
    SocialEngagement_mean = safe_scale(Proximity_mean) + safe_scale(Movement_mean)
  ) %>%
  left_join(early_duration_by_animal, by = c("AnimalNum", "Group", "Sex"))

write_table(feature_long, file.path(output_dir, "tables", "early_behavior_features_long.csv"))
write_table(feature_wide, file.path(output_dir, "tables", "early_behavior_features.csv"))
write_table(feature_wide %>% filter(!contains_short_duration_epoch %in% TRUE), file.path(output_dir, "tables", "early_behavior_features_excluding_short_duration.csv"))

feature_value_cols <- feature_wide %>%
  select(where(is.numeric)) %>%
  names() %>%
  setdiff(c("early_observed_bins", "early_observation_hours", "min_duration_completeness_fraction"))

feature_group_summary <- make_dynamics_group_summary(
  feature_wide,
  value_cols = feature_value_cols,
  group_cols = c("BinLevel", "ProximityInput", "Sex", "Group")
) %>%
  left_join(parse_feature_name(unique(.$Outcome)), by = c("Outcome" = "feature"))

feature_group_contrasts <- make_dynamics_group_contrasts(
  feature_wide,
  value_cols = feature_value_cols,
  by_cols = c("BinLevel", "ProximityInput", "Sex")
) %>%
  mutate(
    ReportingP = p.adjust_bh_family,
    ReportingSignificance = sig_from_p(ReportingP),
    ReportingCorrection = "BH across pairwise group contrasts within sex"
  ) %>%
  left_join(parse_feature_name(unique(.$Outcome)), by = c("Outcome" = "feature"))

write_table(feature_group_summary, file.path(output_dir, "stats_tables", "early_feature_group_summary.csv"))
write_table(feature_group_contrasts, file.path(output_dir, "stats_tables", "early_feature_group_contrasts.csv"))

# ------------------------------------------------
# EARLY FEATURE QC PLOTS
# ------------------------------------------------

early_feature_overview <- feature_wide %>%
  select(AnimalNum, Group, Sex, any_of(c(
    "Movement_mean", "Movement_rmssd", "Movement_acf1",
    "Entropy_mean", "Entropy_rmssd",
    "Proximity_mean", "Proximity_rmssd",
    "SocialWithdrawal_mean", "PassiveIsolation_mean", "SocialEngagement_mean"
  ))) %>%
  pivot_longer(
    cols = -c(AnimalNum, Group, Sex),
    names_to = "feature",
    values_to = "Value"
  ) %>%
  left_join(parse_feature_name(unique(.$feature)), by = "feature") %>%
  filter(is.finite(Value))

if (nrow(early_feature_overview) > 0) {
  p_early_features <- early_feature_overview %>%
    mutate(FeatureLabel = factor(FeatureLabel, levels = unique(FeatureLabel))) %>%
    ggplot(aes(Group, Value, fill = Group)) +
    geom_violin(alpha = 0.55, linewidth = 0.25, trim = FALSE, colour = "grey25") +
    geom_jitter(width = 0.08, size = 0.75, alpha = 0.70, shape = 21, colour = "grey20", stroke = 0.12, show.legend = FALSE) +
    facet_grid(Sex ~ FeatureLabel, scales = "free_y") +
    labs(
      title = "Early-window behavioral feature profiles",
      subtitle = "Endpoint-independent QC view of primary predictive features",
      x = NULL,
      y = "Feature value"
    ) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    make_publication_theme(base_size = 6)

  save_plot_svg_pdf(p_early_features, file.path(output_dir, "figures", "early_window_feature_profiles"), width = 200, height = 120)
}

feature_group_effect_tbl <- feature_group_contrasts %>%
  filter(status == "tested", contrast %in% c("RES-CON", "SUS-CON", "SUS-RES"), !is.na(Metric), !is.na(Measure)) %>%
  mutate(
    Metric = factor(Metric, levels = metric_levels),
    Measure = factor(Measure, levels = measure_levels),
    SigLabel = if_else(ReportingP < 0.05, "*", ""),
    contrast = factor(contrast, levels = c("RES-CON", "SUS-CON", "SUS-RES"))
  )

if (nrow(feature_group_effect_tbl) > 0) {
  p_feature_group_heat <- ggplot(feature_group_effect_tbl, aes(Measure, Metric, fill = cohen_d)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = SigLabel), size = 2.3, colour = "black") +
    facet_grid(Sex ~ contrast) +
    labs(
      title = "Early feature group-difference map",
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

  save_plot_svg_pdf(
    p_feature_group_heat,
    file.path(output_dir, "figures", "publication", "overview", "early_feature_group_effect_size_heatmap"),
    width = 180,
    height = 115
  )
}

# ------------------------------------------------
# ENDPOINT HANDLING
# ------------------------------------------------

endpoint_dat <- NULL
endpoint_source_col <- NA_character_

if (!is.null(endpoint_file) && file.exists(endpoint_file)) {
  endpoint_ext <- tolower(tools::file_ext(endpoint_file))
  endpoint_raw <- if (endpoint_ext %in% c("xlsx", "xls") && !is.null(endpoint_sheet)) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Install readxl to read endpoint Excel files.", call. = FALSE)
    }
    readxl::read_excel(endpoint_file, sheet = endpoint_sheet)
  } else {
    read_behavior_table(endpoint_file)
  }

  endpoint_animal_col <- first_existing_col(endpoint_raw, c("AnimalNum", "Animal", "MouseID", "Mouse", "ID", "RFID", "animal_id"), TRUE, "endpoint animal column")
  endpoint_outcome_col <- first_existing_col(endpoint_raw, unique(c(outcome_col, endpoint_cols)), TRUE, "endpoint outcome column")
  endpoint_source_col <- endpoint_outcome_col

  endpoint_dat <- endpoint_raw %>%
    transmute(
      AnimalNum = as.character(.data[[endpoint_animal_col]]),
      outcome = suppressWarnings(as.numeric(.data[[endpoint_outcome_col]]))
    ) %>%
    group_by(AnimalNum) %>%
    summarise(outcome = first(na.omit(outcome)), .groups = "drop")
} else {
  endpoint_animal_col <- first_existing_col(raw_dat, c("AnimalNum", "Animal", "MouseID", "Mouse", "ID", "RFID", "animal_id"), TRUE, "endpoint animal column")
  endpoint_outcome_col <- first_existing_col(raw_dat, unique(c(outcome_col, endpoint_cols)), FALSE, "endpoint outcome column")
  if (!is.na(endpoint_outcome_col)) {
    endpoint_source_col <- endpoint_outcome_col
    endpoint_dat <- raw_dat %>%
      group_by(AnimalNum = as.character(.data[[endpoint_animal_col]])) %>%
      summarise(outcome = first(na.omit(suppressWarnings(as.numeric(.data[[endpoint_outcome_col]])))), .groups = "drop")
  }
}

if (is.null(endpoint_dat)) {
  message("No endpoint column found. Feature extraction and group summaries finished, but predictive modeling skipped.")
  message("Set endpoint_file/endpoint_sheet/endpoint_cols and/or outcome_col in Analysis/08_early_prediction_models.R.")
  quit(save = "no", status = 0)
}

message("Endpoint source column used for modeling: ", endpoint_source_col)

model_dat <- feature_wide %>%
  mutate(AnimalNum = as.character(AnimalNum)) %>%
  left_join(endpoint_dat, by = "AnimalNum") %>%
  filter(is.finite(outcome))

write_table(model_dat, file.path(output_dir, "tables", "early_prediction_model_input.csv"))

feature_cols <- model_dat %>%
  select(where(is.numeric)) %>%
  select(-outcome) %>%
  names() %>%
  setdiff(c("early_observed_bins", "early_observation_hours", "min_duration_completeness_fraction"))

feature_qc_tbl <- make_feature_qc(model_dat, feature_cols) %>%
  mutate(BinLevel = bin_level, ProximityInput = proximity_col)

endpoint_summary_tbl <- model_dat %>%
  summarise(
    Outcome = outcome_col,
    OutcomeLabel = primary_outcome_label,
    OutcomeDirection = primary_outcome_direction,
    EndpointSourceColumn = endpoint_source_col,
    n_animals = n(),
    n_groups = n_distinct(Group),
    mean = mean(outcome, na.rm = TRUE),
    sd = sd(outcome, na.rm = TRUE),
    median = median(outcome, na.rm = TRUE),
    q25 = quantile(outcome, 0.25, na.rm = TRUE, names = FALSE),
    q75 = quantile(outcome, 0.75, na.rm = TRUE, names = FALSE),
    min = min(outcome, na.rm = TRUE),
    max = max(outcome, na.rm = TRUE)
  ) %>%
  mutate(BinLevel = bin_level, ProximityInput = proximity_col)

write_table(feature_qc_tbl, file.path(output_dir, "tables", "early_feature_qc.csv"))
write_table(endpoint_summary_tbl, file.path(output_dir, "tables", "endpoint_summary.csv"))

if (nrow(model_dat) < 8) {
  warning("Fewer than 8 animals with outcome data. Modeling skipped; feature and input tables were written.")
  quit(save = "no", status = 0)
}

usable_feature_cols <- feature_qc_tbl %>%
  filter(n_finite >= 4, !zero_variance, missing_fraction < 1) %>%
  pull(feature)

if (length(usable_feature_cols) < 2) {
  warning("Fewer than 2 usable features. Modeling skipped; feature and input tables were written.")
  quit(save = "no", status = 0)
}

# ------------------------------------------------
# UNIVARIATE ASSOCIATIONS
# ------------------------------------------------

association_tbl <- map_dfr(usable_feature_cols, function(fc) {
  spearman <- safe_cor_test(model_dat[[fc]], model_dat$outcome, method = "spearman")
  pearson <- safe_cor_test(model_dat[[fc]], model_dat$outcome, method = "pearson")
  tibble(
    feature = fc,
    n = spearman$n,
    spearman_rho = spearman$estimate,
    spearman_statistic = spearman$statistic,
    spearman_p = spearman$p.value,
    pearson_r = pearson$estimate,
    pearson_p = pearson$p.value,
    pearson_conf.low = pearson$conf.low,
    pearson_conf.high = pearson$conf.high
  )
}) %>%
  mutate(
    spearman_p_bh = p.adjust(spearman_p, method = "BH"),
    pearson_p_bh = p.adjust(pearson_p, method = "BH"),
    ReportingP = spearman_p_bh,
    ReportingSignificance = sig_from_p(ReportingP),
    ReportingCorrection = "BH across all early features",
    BinLevel = bin_level,
    ProximityInput = proximity_col
  ) %>%
  left_join(parse_feature_name(usable_feature_cols), by = "feature") %>%
  arrange(ReportingP, desc(abs(spearman_rho)))

write_table(association_tbl, file.path(output_dir, "tables", "early_feature_outcome_associations.csv"))

association_shortlist <- association_tbl %>%
  filter(!is.na(ReportingP)) %>%
  mutate(
    EvidenceClass = case_when(
      ReportingP < 0.05 & abs(spearman_rho) >= 0.6 ~ "Strong association",
      ReportingP < 0.05 ~ "FDR-supported",
      ReportingP < 0.10 ~ "Trend",
      abs(spearman_rho) >= 0.6 ~ "Large rho, uncertain",
      TRUE ~ "No clear evidence"
    ),
    Direction = if_else(spearman_rho >= 0, "Higher early feature predicts higher endpoint", "Higher early feature predicts lower endpoint"),
    PriorityScore = if_else(ReportingP > 0, -log10(ReportingP) * abs(spearman_rho), 0)
  ) %>%
  arrange(desc(EvidenceClass %in% c("Strong association", "FDR-supported")), desc(PriorityScore)) %>%
  slice_head(n = 40)

write_table(association_shortlist, file.path(output_dir, "tables", "early_feature_association_shortlist.csv"))

# ------------------------------------------------
# FEATURE SPACE SUMMARIES
# ------------------------------------------------

x_imputed <- impute_feature_matrix(model_dat, usable_feature_cols)
x_scaled <- standardize_matrix(x_imputed)

feature_pca <- stats::prcomp(x_scaled, center = FALSE, scale. = FALSE)
pca_score_tbl <- as_tibble(feature_pca$x[, seq_len(min(3, ncol(feature_pca$x))), drop = FALSE]) %>%
  bind_cols(model_dat %>% select(AnimalNum, Group, Sex, outcome)) %>%
  mutate(BinLevel = bin_level, ProximityInput = proximity_col)

pca_variance_tbl <- tibble(
  PC = paste0("PC", seq_along(feature_pca$sdev)),
  variance_explained = feature_pca$sdev^2 / sum(feature_pca$sdev^2),
  cumulative_variance_explained = cumsum(variance_explained)
)

pca_loading_tbl <- as_tibble(feature_pca$rotation, rownames = "feature") %>%
  left_join(parse_feature_name(usable_feature_cols), by = "feature") %>%
  arrange(desc(abs(PC1)))

write_table(pca_score_tbl, file.path(output_dir, "tables", "early_feature_pca_scores.csv"))
write_table(pca_variance_tbl, file.path(output_dir, "tables", "early_feature_pca_variance.csv"))
write_table(pca_loading_tbl, file.path(output_dir, "tables", "early_feature_pca_loadings.csv"))

# ------------------------------------------------
# LEAVE-ONE-ANIMAL-OUT ELASTIC-NET IF glmnet AVAILABLE
# ------------------------------------------------

if (requireNamespace("glmnet", quietly = TRUE)) {
  run_duration_subset_elastic_net <- function(dat, feature_cols, label) {
    if (nrow(dat) < 8 || length(feature_cols) < 2) {
      return(list(predictions = tibble(), performance = tibble(DurationAnalysisSet = label, status = "skipped_low_n")))
    }
    x_subset <- impute_feature_matrix(dat, feature_cols)
    y_subset <- dat$outcome
    loo_pred_subset <- rep(NA_real_, length(y_subset))
    set.seed(123)
    for (i in seq_along(y_subset)) {
      train_idx <- setdiff(seq_along(y_subset), i)
      fit <- glmnet::cv.glmnet(
        x_subset[train_idx, , drop = FALSE],
        y_subset[train_idx],
        alpha = 0.5,
        standardize = TRUE,
        nfolds = min(5, length(train_idx))
      )
      loo_pred_subset[i] <- as.numeric(predict(fit, newx = x_subset[i, , drop = FALSE], s = "lambda.min"))
    }
    pred_subset <- dat %>%
      transmute(
        AnimalNum, Group, Sex, BinLevel, ProximityInput,
        observed = outcome,
        predicted = loo_pred_subset,
        residual = observed - predicted,
        abs_residual = abs(residual),
        DurationAnalysisSet = label
      )
    perf_subset <- make_prediction_summary(pred_subset, outcome_col = outcome_col, n_perm = n_prediction_permutations) %>%
      mutate(BinLevel = bin_level, ProximityInput = proximity_col, Model = "LOO elastic-net alpha=0.5", DurationAnalysisSet = label, status = "fit")
    list(predictions = pred_subset, performance = perf_subset)
  }

  x <- x_imputed
  y <- model_dat$outcome
  loo_pred <- rep(NA_real_, length(y))
  coef_tbl <- tibble()
  lambda_tbl <- tibble()

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

    this_coef <- as.matrix(stats::coef(fit, s = "lambda.min"))[, 1]
    coef_tbl <- bind_rows(
      coef_tbl,
      tibble(held_out_animal = model_dat$AnimalNum[i], feature = names(this_coef), coefficient = as.numeric(this_coef))
    )
    lambda_tbl <- bind_rows(
      lambda_tbl,
      tibble(held_out_animal = model_dat$AnimalNum[i], lambda_min = fit$lambda.min, lambda_1se = fit$lambda.1se, cvm_min = min(fit$cvm, na.rm = TRUE))
    )
  }

  pred_tbl <- model_dat %>%
    transmute(
      AnimalNum, Group, Sex, BinLevel, ProximityInput,
      observed = outcome,
      predicted = loo_pred,
      residual = observed - predicted,
      abs_residual = abs(residual)
    )

  perf_tbl <- make_prediction_summary(pred_tbl, outcome_col = outcome_col, n_perm = n_prediction_permutations) %>%
    mutate(BinLevel = bin_level, ProximityInput = proximity_col, Model = "LOO elastic-net alpha=0.5")

  coefficient_summary_tbl <- summarise_model_coefficients(coef_tbl) %>%
    mutate(BinLevel = bin_level, ProximityInput = proximity_col)

  write_table(pred_tbl, file.path(output_dir, "tables", "leave_one_animal_out_predictions.csv"))
  write_table(perf_tbl, file.path(output_dir, "tables", "leave_one_animal_out_performance.csv"))
  write_table(coef_tbl, file.path(output_dir, "tables", "leave_one_animal_out_model_coefficients.csv"))
  write_table(lambda_tbl, file.path(output_dir, "tables", "leave_one_animal_out_lambda_diagnostics.csv"))
  write_table(coefficient_summary_tbl, file.path(output_dir, "tables", "elastic_net_coefficient_stability.csv"))

  no_short_model_dat <- model_dat %>% filter(!contains_short_duration_epoch %in% TRUE)
  duration_prediction_full <- list(
    predictions = pred_tbl %>% mutate(DurationAnalysisSet = "full"),
    performance = perf_tbl %>% mutate(DurationAnalysisSet = "full", status = "fit")
  )
  duration_prediction_no_short <- run_duration_subset_elastic_net(no_short_model_dat, usable_feature_cols, "excluding_short_duration")
  duration_prediction_perf <- bind_rows(duration_prediction_full$performance, duration_prediction_no_short$performance) %>%
    mutate(
      delta_cv_r2_vs_full = cross_validated_r2_vs_mean - cross_validated_r2_vs_mean[DurationAnalysisSet == "full"][1],
      delta_pearson_r_vs_full = loo_pearson_r - loo_pearson_r[DurationAnalysisSet == "full"][1],
      delta_rmse_vs_full = loo_rmse - loo_rmse[DurationAnalysisSet == "full"][1]
    )
  write_table(bind_rows(duration_prediction_full$predictions, duration_prediction_no_short$predictions), file.path(output_dir, "tables", "leave_one_animal_out_predictions_duration_sensitivity.csv"))
  write_table(duration_prediction_perf, file.path(output_dir, "tables", "leave_one_animal_out_performance_duration_sensitivity.csv"))

  p_pred <- pred_tbl %>%
    ggplot(aes(observed, predicted, colour = Group, fill = Group)) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey45") +
    geom_point(size = 1.8, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.55, alpha = 0.14, colour = "grey20") +
    facet_grid(. ~ Sex) +
    labs(
      title = "Early behavior predicts later endpoint",
      subtitle = paste0("LOO elastic-net; r=", round(perf_tbl$loo_pearson_r, 2), ", CV R2=", round(perf_tbl$cross_validated_r2_vs_mean, 2), ", permutation ", format_p_label(perf_tbl$prediction_permutation_p)),
      x = paste0("Observed ", outcome_col),
      y = paste0("Predicted ", outcome_col)
    ) +
    scale_colour_manual(values = group_colors, drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    make_publication_theme()

  save_plot_svg_pdf(p_pred, file.path(output_dir, "figures", "early_prediction_observed_vs_predicted"), width = 120, height = 80)

  p_resid <- pred_tbl %>%
    ggplot(aes(predicted, residual, colour = Group, fill = Group)) +
    geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey45") +
    geom_point(size = 1.8, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.5, alpha = 0.14, colour = "grey20") +
    facet_grid(. ~ Sex) +
    labs(
      title = "Prediction residuals",
      subtitle = "Residual = observed - predicted",
      x = paste0("Predicted ", outcome_col),
      y = "Residual"
    ) +
    scale_colour_manual(values = group_colors, drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    make_publication_theme()

  save_plot_svg_pdf(p_resid, file.path(output_dir, "figures", "early_prediction_residuals"), width = 120, height = 80)

  p_coef <- coefficient_summary_tbl %>%
    slice_max(order_by = nonzero_frequency * abs(median_coefficient), n = 18, with_ties = FALSE) %>%
    mutate(FeatureLabel = factor(FeatureLabel, levels = rev(FeatureLabel[order(nonzero_frequency * abs(median_coefficient))]))) %>%
    ggplot(aes(median_coefficient, FeatureLabel, colour = Metric, size = nonzero_frequency)) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey45") +
    geom_errorbar(aes(xmin = q25_coefficient, xmax = q75_coefficient), orientation = "y", width = 0, linewidth = 0.28, alpha = 0.85) +
    geom_point(alpha = 0.9) +
    labs(
      title = "Elastic-net coefficient stability",
      subtitle = "Top features by selection frequency and coefficient magnitude across LOO fits",
      x = "Median elastic-net coefficient",
      y = NULL,
      colour = NULL,
      size = "Selection frequency"
    ) +
    make_publication_theme(base_size = 6) +
    theme(panel.grid.major.y = element_blank())

  save_plot_svg_pdf(p_coef, file.path(output_dir, "figures", "elastic_net_coefficient_stability"), width = 145, height = 95)
} else {
  message("Package glmnet not installed. Wrote feature associations, but skipped cross-validated elastic-net modeling.")
}

# ------------------------------------------------
# PUBLICATION-STYLE SUMMARY FIGURES
# ------------------------------------------------

p_assoc <- association_tbl %>%
  slice_min(order_by = ReportingP, n = 18, with_ties = FALSE) %>%
  mutate(
    FeatureLabel = factor(FeatureLabel, levels = rev(FeatureLabel[order(abs(spearman_rho))])),
    SignificantLabel = if_else(ReportingP < 0.05, "FDR < 0.05", "n.s.")
  ) %>%
  ggplot(aes(spearman_rho, FeatureLabel, fill = Metric, alpha = SignificantLabel)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey45") +
  geom_col(width = 0.68, colour = "grey20", linewidth = 0.15) +
  labs(
    title = "Top early behavioral correlates",
    subtitle = "Ranked by BH-adjusted Spearman association with endpoint",
    x = paste0("Spearman rho with ", outcome_col),
    y = NULL,
    fill = NULL,
    alpha = NULL
  ) +
  scale_alpha_manual(values = c("FDR < 0.05" = 1, "n.s." = 0.55), drop = FALSE) +
  make_publication_theme(base_size = 6) +
  theme(panel.grid.major.y = element_blank())

save_plot_svg_pdf(p_assoc, file.path(output_dir, "figures", "top_early_feature_associations"), width = 130, height = 95)

assoc_heat_tbl <- association_tbl %>%
  filter(!is.na(Metric), !is.na(Measure)) %>%
  mutate(
    Metric = factor(Metric, levels = metric_levels),
    Measure = factor(Measure, levels = measure_levels),
    SigLabel = if_else(ReportingP < 0.05, "*", "")
  )

p_assoc_heat <- ggplot(assoc_heat_tbl, aes(Measure, Metric, fill = spearman_rho)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = SigLabel), size = 2.5, colour = "black") +
  labs(
    title = "Early feature association map",
    subtitle = "Spearman rho; asterisk marks BH-adjusted FDR < 0.05",
    x = NULL,
    y = NULL,
    fill = "rho"
  ) +
  scale_fill_gradient2(low = "#3d3b6e", mid = "white", high = "#e63947", midpoint = 0, na.value = "grey90") +
  make_publication_theme(base_size = 6) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  )

save_plot_svg_pdf(p_assoc_heat, file.path(output_dir, "figures", "publication", "overview", "early_feature_association_heatmap"), width = 150, height = 85)

pc1_lab <- paste0("PC1 (", round(100 * pca_variance_tbl$variance_explained[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(100 * pca_variance_tbl$variance_explained[2], 1), "%)")

p_pca <- pca_score_tbl %>%
  ggplot(aes(PC1, PC2, colour = Group, fill = Group)) +
  geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey88") +
  geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey88") +
  geom_point(aes(size = outcome), alpha = 0.82) +
  facet_grid(. ~ Sex) +
  labs(
    title = "Early behavior feature space",
    subtitle = "Point size encodes endpoint value",
    x = pc1_lab,
    y = pc2_lab,
    size = outcome_col
  ) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  make_publication_theme()

save_plot_svg_pdf(p_pca, file.path(output_dir, "figures", "early_feature_space_pca"), width = 120, height = 80)

top_feature_values <- association_tbl %>%
  slice_min(order_by = ReportingP, n = 6, with_ties = FALSE) %>%
  select(feature, FeatureLabel)

feature_scatter_tbl <- model_dat %>%
  select(AnimalNum, Group, Sex, outcome, all_of(top_feature_values$feature)) %>%
  pivot_longer(cols = all_of(top_feature_values$feature), names_to = "feature", values_to = "FeatureValue") %>%
  left_join(top_feature_values, by = "feature") %>%
  mutate(FeatureLabel = factor(FeatureLabel, levels = top_feature_values$FeatureLabel))

p_feature_scatter <- feature_scatter_tbl %>%
  ggplot(aes(FeatureValue, outcome, colour = Group, fill = Group)) +
  geom_point(size = 1.5, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.45, alpha = 0.12, colour = "grey20") +
  facet_wrap(~FeatureLabel, scales = "free_x", ncol = 3) +
  labs(
    title = "Top early features versus endpoint",
    subtitle = "Panels show the strongest BH-ranked univariate associations",
    x = "Early feature value",
    y = outcome_col
  ) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  make_publication_theme(base_size = 6)

save_plot_svg_pdf(p_feature_scatter, file.path(output_dir, "figures", "top_feature_endpoint_scatterplots"), width = 160, height = 110)

message("Early prediction analysis complete. Feature QC, associations, prediction diagnostics, and publication figures written.")
