# ================================================================
# Early Behavioral Prediction Model Ladder
# MMMSociability
# ================================================================
# Goal:
#   Test whether early first-active-phase behavior after the first cage
#   change predicts later composite stress burden (CombZ), with explicit
#   model comparison between conventional magnitude features and temporal
#   organization features.
#
# Biological use case:
#   P25 first cage change, first 12 h active phase, 5-min bins.
#   Primary features: Movement_mean and Entropy_acf1.
#
# Main question:
#   Does early entropy persistence explain later CombZ beyond movement,
#   group labels, sex, and batch/cage-change covariates?
#
# Input expectation:
#   Run Analysis/03_build_multiscale_behavior_metrics.R first.
#   The endpoint file must contain AnimalNum and CombZ or your selected
#   outcome column.
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

bin_level <- "5min_based"
base_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID"
input_file <- file.path(base_dir, "analysis_ready/03_derived_metrics", bin_level, "all_behavior_metrics.csv")
output_dir <- file.path(base_dir, "analysis_ready/06_behavioral_dynamics/early_prediction_model_ladder", bin_level)

# Endpoint file should contain one row per animal or repeated rows with a stable endpoint.
# If NULL, the script tries to read the outcome from input_file.
endpoint_file <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/SIS_Analysis/E9_Behavior_Data.xlsx"
endpoint_excel_sheet <- "zScore"  # sheet name for Excel endpoint files; NULL = first sheet
outcome_col <- "CombZ"

# Primary prospective window: first 12 h active phase after the first cage change.
early_phase_pattern <- "active|dark|night"
first_cage_change_only <- TRUE
early_window_hours <- 12
bin_size_min <- 5
max_early_bins_per_animal <- early_window_hours * 60 / bin_size_min

# Optional column overrides. Leave NULL for automatic detection by helper functions.
animal_col <- NULL
time_col <- NULL
group_col <- NULL
sex_col <- NULL
phase_col <- NULL
cage_col <- NULL
movement_col <- NULL
entropy_col <- NULL
proximity_col <- "ProximityFraction"

# Model options.
n_prediction_permutations <- 5000
n_bootstrap <- 5000
set.seed(123)

group_colors <- c(
  "CON" = "#3d3b6e",
  "RES" = "#C6C3BB",
  "SUS" = "#e63947",
  "All" = "grey55"
)
group_levels <- c("CON", "RES", "SUS")

# ------------------------------------------------
# SMALL HELPERS
# ------------------------------------------------

ensure_dir_safe <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

safe_scale <- function(x) {
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

format_p <- function(p) {
  case_when(
    is.na(p) ~ "p=NA",
    p < 0.001 ~ "p<.001",
    TRUE ~ paste0("p=", sub("^0", "", formatC(p, format = "f", digits = 3)))
  )
}

make_publication_theme <- function(base_size = 7) {
  make_nature_theme(base_size = base_size) +
    theme(
      legend.position = "top",
      panel.grid.major.y = element_line(linewidth = 0.15, colour = "grey92"),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1),
      plot.subtitle = element_text(hjust = 0, size = base_size)
    )
}

get_first_cage_change <- function(x) {
  ux <- unique(as.character(x))
  cc_num <- suppressWarnings(as.numeric(str_extract(ux, "\\d+")))
  if (any(is.finite(cc_num))) ux[which.min(ifelse(is.finite(cc_num), cc_num, Inf))] else sort(ux)[1]
}

impute_numeric <- function(dat, cols) {
  dat %>%
    mutate(across(all_of(cols), ~{
      med <- median(.x, na.rm = TRUE)
      if (!is.finite(med)) med <- 0
      replace_na(.x, med)
    }))
}

prediction_metrics <- function(observed, predicted) {
  ok <- is.finite(observed) & is.finite(predicted)
  y <- observed[ok]
  p <- predicted[ok]
  baseline <- mean(y, na.rm = TRUE)
  tibble(
    n = length(y),
    pearson_r = safe_cor(y, p, "pearson"),
    spearman_rho = safe_cor(y, p, "spearman"),
    rmse = sqrt(mean((y - p)^2, na.rm = TRUE)),
    mae = mean(abs(y - p), na.rm = TRUE),
    baseline_rmse = sqrt(mean((y - baseline)^2, na.rm = TRUE)),
    cv_r2_vs_mean = 1 - sum((y - p)^2, na.rm = TRUE) / sum((y - baseline)^2, na.rm = TRUE)
  )
}

permutation_prediction_p <- function(observed, predicted, n_perm = 5000, seed = 123) {
  ok <- is.finite(observed) & is.finite(predicted)
  y <- observed[ok]
  p <- predicted[ok]
  if (length(y) < 4 || sd(y) == 0 || sd(p) == 0) return(NA_real_)
  r_obs <- suppressWarnings(cor(y, p, method = "pearson"))
  set.seed(seed)
  r_null <- replicate(n_perm, suppressWarnings(cor(y, sample(p), method = "pearson")))
  (sum(abs(r_null) >= abs(r_obs), na.rm = TRUE) + 1) / (sum(is.finite(r_null)) + 1)
}

make_formula <- function(outcome = "outcome", predictors) {
  if (length(predictors) == 0) as.formula(paste0(outcome, " ~ 1")) else as.formula(paste(outcome, "~", paste(predictors, collapse = " + ")))
}

loo_lm_predict <- function(dat, predictors, model_name) {
  pred <- rep(NA_real_, nrow(dat))
  coef_rows <- list()

  for (i in seq_len(nrow(dat))) {
    train <- dat[-i, , drop = FALSE]
    test <- dat[i, , drop = FALSE]
    form <- make_formula("outcome", predictors)
    fit <- try(lm(form, data = train), silent = TRUE)
    if (inherits(fit, "try-error")) next
    pred[i] <- tryCatch(as.numeric(predict(fit, newdata = test)), error = function(e) NA_real_)
    coef_rows[[i]] <- broom_like_coef(fit, held_out_animal = dat$AnimalNum[i], model_name = model_name)
  }

  pred_tbl <- dat %>%
    transmute(
      AnimalNum, Group, Sex, BinLevel, ProximityInput,
      observed = outcome,
      predicted = pred,
      residual = observed - predicted,
      abs_residual = abs(residual),
      Model = model_name
    )

  list(predictions = pred_tbl, coefficients = bind_rows(coef_rows))
}

broom_like_coef <- function(fit, held_out_animal, model_name) {
  sm <- try(summary(fit)$coefficients, silent = TRUE)
  if (inherits(sm, "try-error")) return(tibble())
  tibble(
    held_out_animal = held_out_animal,
    Model = model_name,
    term = rownames(sm),
    estimate = sm[, 1],
    std_error = sm[, 2],
    statistic = sm[, 3],
    p_value = sm[, 4]
  )
}

bootstrap_correlation_ci <- function(dat, x_col, y_col = "outcome", n_boot = 5000, method = "spearman", seed = 123) {
  if (!x_col %in% names(dat)) {
    return(tibble(feature = x_col, n = nrow(dat), estimate = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }
  ok_dat <- dat %>% filter(is.finite(.data[[x_col]]), is.finite(.data[[y_col]]))
  if (nrow(ok_dat) < 4 || sd(ok_dat[[x_col]]) == 0 || sd(ok_dat[[y_col]]) == 0) {
    return(tibble(feature = x_col, n = nrow(ok_dat), estimate = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }
  set.seed(seed)
  boots <- replicate(n_boot, {
    idx <- sample(seq_len(nrow(ok_dat)), replace = TRUE)
    safe_cor(ok_dat[[x_col]][idx], ok_dat[[y_col]][idx], method = method)
  })
  tibble(
    feature = x_col,
    n = nrow(ok_dat),
    estimate = safe_cor(ok_dat[[x_col]], ok_dat[[y_col]], method = method),
    ci_low = quantile(boots, 0.025, na.rm = TRUE, names = FALSE),
    ci_high = quantile(boots, 0.975, na.rm = TRUE, names = FALSE),
    n_bootstrap = n_boot,
    method = method
  )
}

partial_r_from_lm <- function(dat, feature, covariates, outcome = "outcome") {
  needed <- c(outcome, feature, covariates)
  needed <- needed[needed %in% names(dat)]
  d <- dat %>% select(all_of(needed)) %>% drop_na()
  if (nrow(d) < 6 || !feature %in% names(d)) return(NA_real_)
  covariates <- covariates[covariates %in% names(d)]
  if (length(covariates) == 0) return(safe_cor(d[[feature]], d[[outcome]], "pearson"))
  ry <- residuals(lm(make_formula(outcome, covariates), data = d))
  rx <- residuals(lm(make_formula(feature, covariates), data = d))
  safe_cor(rx, ry, "pearson")
}

# ------------------------------------------------
# LOAD DATA
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
  mutate(
    Group = factor(as.character(Group), levels = unique(c(group_levels, sort(unique(as.character(Group))))))
  )

ensure_dir_safe(output_dir)
ensure_dir_safe(file.path(output_dir, "tables"))
ensure_dir_safe(file.path(output_dir, "figures"))
ensure_dir_safe(file.path(output_dir, "figures", "publication"))
output_dirs <- analysis_output_dirs(output_dir)
write_output_manifest(
  output_dir,
  script_name = "08b_early_prediction_model_ladder.R",
  analysis_name = "early prediction model ladder",
  primary_tables = c(
    "tables/model_ladder_performance.csv",
    "tables/model_ladder_performance_duration_sensitivity.csv",
    "tables/model_ladder_repeated_grouped_kfold_performance.csv",
    "tables/prediction_interpretation_constraints.csv",
    "tables/model_ladder_incremental_summary.csv",
    "tables/primary_movement_entropyacf1_associations.csv"
  ),
  primary_figures = c(
    "figures/publication/model_ladder_cv_r2.svg",
    "figures/publication/primary_movement_entropyacf1_vs_combz.svg"
  ),
  notes = c("Main prediction claim should use the ladder performance plus duration-sensitivity companion table.")
)

epoch_duration_qc <- write_epoch_duration_qc(behav, output_dir, metric_source = "08b_early_prediction_model_ladder", bin_size_sec = infer_bin_size_sec(behav))

# ------------------------------------------------
# EARLY WINDOW: FIRST CAGE CHANGE, FIRST 12 h ACTIVE PHASE
# ------------------------------------------------

if (first_cage_change_only && "CageChange" %in% names(behav)) {
  first_cc <- get_first_cage_change(behav$CageChange)
  behav <- behav %>% filter(as.character(CageChange) == first_cc)
} else {
  first_cc <- "all"
}

has_active_phase <- any(str_detect(str_to_lower(as.character(behav$Phase)), early_phase_pattern))
early_dat <- if (has_active_phase) {
  behav %>% filter(str_detect(str_to_lower(as.character(Phase)), early_phase_pattern))
} else {
  behav
}

early_dat <- early_dat %>%
  group_by(AnimalNum, Phase) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  mutate(early_rank = row_number()) %>%
  filter(early_rank <= max_early_bins_per_animal) %>%
  ungroup() %>%
  mutate(BinLevel = bin_level, ProximityInput = proximity_col)

write_table(early_dat, file.path(output_dir, "tables", "early_window_rows_used.csv"))
write_table(filter_short_duration_epochs(early_dat, epoch_duration_qc), file.path(output_dir, "tables", "early_window_rows_used_excluding_short_duration.csv"))

window_design_tbl <- early_dat %>%
  group_by(AnimalNum, Group, Sex, Phase) %>%
  summarise(
    n_bins = n(),
    approx_hours = n_bins * bin_size_min / 60,
    first_time_index = min(TimeIndex, na.rm = TRUE),
    last_time_index = max(TimeIndex, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    BinLevel = bin_level,
    FirstCageChangeOnly = first_cage_change_only,
    FirstCageChange = first_cc,
    TargetWindowHours = early_window_hours
  )

write_table(window_design_tbl, file.path(output_dir, "tables", "early_window_design_by_animal.csv"))

# ------------------------------------------------
# FEATURE EXTRACTION
# ------------------------------------------------

feature_long <- early_dat %>%
  pivot_longer(cols = c(Movement, Entropy, Proximity), names_to = "Metric", values_to = "Value") %>%
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
    EarlyWindow = paste0("first_", early_window_hours, "h_active_first_cage_change"),
    Movement_z = safe_scale(Movement_mean),
    EntropyACF1_z = safe_scale(Entropy_acf1),
    Movement_x_EntropyACF1 = Movement_z * EntropyACF1_z
  ) %>%
  left_join(early_duration_by_animal, by = c("AnimalNum", "Group", "Sex"))

write_table(feature_long, file.path(output_dir, "tables", "early_behavior_features_long.csv"))
write_table(feature_wide, file.path(output_dir, "tables", "early_behavior_features_wide.csv"))
write_table(feature_wide %>% filter(!contains_short_duration_epoch %in% TRUE), file.path(output_dir, "tables", "early_behavior_features_wide_excluding_short_duration.csv"))

# ------------------------------------------------
# ENDPOINT HANDLING
# ------------------------------------------------

endpoint_dat <- NULL
if (!is.null(endpoint_file) && file.exists(endpoint_file)) {
  ext <- tools::file_ext(endpoint_file) %>% tolower()
  endpoint_raw <- if (ext %in% c("xlsx", "xls") && !is.null(endpoint_excel_sheet)) {
    readxl::read_excel(endpoint_file, sheet = endpoint_excel_sheet)
  } else {
    read_behavior_table(endpoint_file)
  }
  endpoint_animal_col <- first_existing_col(endpoint_raw, c("AnimalNum", "Animal", "MouseID", "Mouse", "ID", "RFID", "animal_id"), TRUE, "endpoint animal column")
  endpoint_dat <- endpoint_raw %>%
    transmute(AnimalNum = .data[[endpoint_animal_col]], outcome = suppressWarnings(as.numeric(.data[[outcome_col]]))) %>%
    group_by(AnimalNum) %>%
    summarise(outcome = first(na.omit(outcome)), .groups = "drop")
} else if (outcome_col %in% names(raw_dat)) {
  endpoint_animal_col <- first_existing_col(raw_dat, c("AnimalNum", "Animal", "MouseID", "Mouse", "ID", "RFID", "animal_id"), TRUE, "endpoint animal column")
  endpoint_dat <- raw_dat %>%
    group_by(AnimalNum = .data[[endpoint_animal_col]]) %>%
    summarise(outcome = first(na.omit(suppressWarnings(as.numeric(.data[[outcome_col]])))), .groups = "drop")
}

if (is.null(endpoint_dat)) {
  stop("No endpoint data found. Set endpoint_file or ensure outcome_col is present in the input file.")
}

model_dat <- feature_wide %>%
  left_join(endpoint_dat, by = "AnimalNum") %>%
  filter(is.finite(outcome)) %>%
  mutate(
    Sex = factor(as.character(Sex)),
    Group = factor(as.character(Group), levels = unique(c(group_levels, sort(unique(as.character(Group))))))
  )

if (nrow(model_dat) < 8) stop("Fewer than 8 animals with endpoint data. Model ladder not reliable.")

write_table(model_dat, file.path(output_dir, "tables", "model_ladder_input.csv"))

# ------------------------------------------------
# PRIMARY FEATURE ASSOCIATIONS
# ------------------------------------------------

primary_features <- c("Movement_mean", "Entropy_acf1", "Movement_x_EntropyACF1")
primary_features <- primary_features[primary_features %in% names(model_dat)]

primary_assoc <- map_dfr(primary_features, function(fc) {
  cor_s <- suppressWarnings(cor.test(model_dat[[fc]], model_dat$outcome, method = "spearman", exact = FALSE))
  cor_p <- suppressWarnings(cor.test(model_dat[[fc]], model_dat$outcome, method = "pearson"))
  boot <- bootstrap_correlation_ci(model_dat, fc, "outcome", n_bootstrap, method = "spearman")
  tibble(
    feature = fc,
    n = sum(is.finite(model_dat[[fc]]) & is.finite(model_dat$outcome)),
    spearman_rho = unname(cor_s$estimate),
    spearman_p = cor_s$p.value,
    pearson_r = unname(cor_p$estimate),
    pearson_p = cor_p$p.value,
    spearman_boot_ci_low = boot$ci_low,
    spearman_boot_ci_high = boot$ci_high,
    partial_r_controlling_movement = if_else(fc == "Entropy_acf1", partial_r_from_lm(model_dat, fc, c("Movement_mean", "Sex", "Group")), NA_real_),
    BinLevel = bin_level,
    Outcome = outcome_col
  )
}) %>%
  mutate(
    spearman_p_bh = p.adjust(spearman_p, method = "BH"),
    Evidence = case_when(
      spearman_p_bh < 0.05 & sign(spearman_boot_ci_low) == sign(spearman_boot_ci_high) ~ "FDR-supported; bootstrap CI excludes zero",
      spearman_p < 0.05 ~ "nominal",
      TRUE ~ "uncertain"
    )
  ) %>%
  arrange(spearman_p_bh)

write_table(primary_assoc, file.path(output_dir, "tables", "primary_movement_entropyacf1_associations.csv"))

# Sex-specific associations, because the biological effect appears female-specific.
sex_specific_assoc <- model_dat %>%
  group_by(Sex) %>%
  group_modify(~{
    map_dfr(primary_features, function(fc) {
      if (nrow(.x) < 4 || sd(.x[[fc]], na.rm = TRUE) == 0 || sd(.x$outcome, na.rm = TRUE) == 0) {
        return(tibble(feature = fc, n = nrow(.x), spearman_rho = NA_real_, spearman_p = NA_real_))
      }
      ct <- suppressWarnings(cor.test(.x[[fc]], .x$outcome, method = "spearman", exact = FALSE))
      tibble(feature = fc, n = nrow(.x), spearman_rho = unname(ct$estimate), spearman_p = ct$p.value)
    })
  }) %>%
  ungroup() %>%
  group_by(Sex) %>%
  mutate(spearman_p_bh_within_sex = p.adjust(spearman_p, method = "BH")) %>%
  ungroup() %>%
  mutate(BinLevel = bin_level, Outcome = outcome_col)

write_table(sex_specific_assoc, file.path(output_dir, "tables", "sex_specific_primary_associations.csv"))

# ------------------------------------------------
# MODEL LADDER: EXPLICIT INCREMENTAL PREDICTION
# ------------------------------------------------

candidate_covars <- c("Sex", "Group")
candidate_covars <- candidate_covars[candidate_covars %in% names(model_dat)]

model_specs <- list(
  "Mean only" = character(0),
  "Sex + Group" = candidate_covars,
  "Movement only" = c(candidate_covars, "Movement_mean"),
  "Entropy ACF1 only" = c(candidate_covars, "Entropy_acf1"),
  "Movement + Entropy ACF1" = c(candidate_covars, "Movement_mean", "Entropy_acf1"),
  "Movement x Entropy ACF1" = c(candidate_covars, "Movement_mean", "Entropy_acf1", "Movement_x_EntropyACF1"),
  "Full behavior compact" = c(candidate_covars, "Movement_mean", "Movement_rmssd", "Movement_acf1", "Entropy_mean", "Entropy_rmssd", "Entropy_acf1", "Proximity_mean", "Proximity_rmssd", "Proximity_acf1")
)

model_specs <- map(model_specs, ~.x[.x %in% names(model_dat)])

numeric_predictors <- unique(unlist(model_specs))
numeric_predictors <- numeric_predictors[numeric_predictors %in% names(model_dat) & sapply(model_dat[numeric_predictors], is.numeric)]
model_dat_imputed <- impute_numeric(model_dat, numeric_predictors)

ladder_results <- imap(model_specs, ~loo_lm_predict(model_dat_imputed, .x, .y))
ladder_predictions <- map_dfr(ladder_results, "predictions")
ladder_coefficients <- map_dfr(ladder_results, "coefficients")

ladder_performance <- ladder_predictions %>%
  group_by(Model) %>%
  group_modify(~prediction_metrics(.x$observed, .x$predicted)) %>%
  ungroup() %>%
  mutate(
    prediction_permutation_p = map_dbl(Model, ~{
      pdat <- ladder_predictions %>% filter(Model == .x)
      permutation_prediction_p(pdat$observed, pdat$predicted, n_prediction_permutations, seed = 123)
    }),
    BinLevel = bin_level,
    Outcome = outcome_col
  ) %>%
  arrange(desc(cv_r2_vs_mean), rmse)

run_ladder_duration_set <- function(dat, analysis_set) {
  if (nrow(dat) < 8 || sum(is.finite(dat$outcome)) < 8) {
    return(list(
      predictions = tibble(),
      coefficients = tibble(),
      performance = tibble(
        Model = names(model_specs),
        n = sum(is.finite(dat$outcome)),
        pearson_r = NA_real_,
        spearman_rho = NA_real_,
        rmse = NA_real_,
        mae = NA_real_,
        cv_r2_vs_mean = NA_real_,
        prediction_permutation_p = NA_real_,
        BinLevel = bin_level,
        Outcome = outcome_col,
        DurationAnalysisSet = analysis_set,
        DurationSensitivityStatus = "skipped_too_few_animals"
      )
    ))
  }

  dat_imp <- impute_numeric(dat, numeric_predictors)
  fits <- imap(model_specs, ~loo_lm_predict(dat_imp, .x, .y))
  preds <- map_dfr(fits, "predictions") %>% mutate(DurationAnalysisSet = analysis_set)
  coefs <- map_dfr(fits, "coefficients") %>% mutate(DurationAnalysisSet = analysis_set)
  perf <- preds %>%
    group_by(Model) %>%
    group_modify(~prediction_metrics(.x$observed, .x$predicted)) %>%
    ungroup() %>%
    mutate(
      prediction_permutation_p = map_dbl(Model, ~{
        pdat <- preds %>% filter(Model == .x)
        permutation_prediction_p(pdat$observed, pdat$predicted, n_prediction_permutations, seed = 123)
      }),
      BinLevel = bin_level,
      Outcome = outcome_col,
      DurationAnalysisSet = analysis_set,
      DurationSensitivityStatus = "fit"
    ) %>%
    arrange(desc(cv_r2_vs_mean), rmse)

  list(predictions = preds, coefficients = coefs, performance = perf)
}

full_duration_ladder <- list(
  predictions = ladder_predictions %>% mutate(DurationAnalysisSet = "full"),
  coefficients = ladder_coefficients %>% mutate(DurationAnalysisSet = "full"),
  performance = ladder_performance %>%
    mutate(DurationAnalysisSet = "full", DurationSensitivityStatus = "fit")
)

no_short_model_dat <- model_dat %>%
  filter(!contains_short_duration_epoch %in% TRUE)
no_short_duration_ladder <- run_ladder_duration_set(no_short_model_dat, "excluding_short_duration")

ladder_predictions_duration_sensitivity <- bind_rows(
  full_duration_ladder$predictions,
  no_short_duration_ladder$predictions
)
ladder_coefficients_duration_sensitivity <- bind_rows(
  full_duration_ladder$coefficients,
  no_short_duration_ladder$coefficients
)
ladder_performance_duration_sensitivity <- bind_rows(
  full_duration_ladder$performance,
  no_short_duration_ladder$performance
) %>%
  group_by(Model) %>%
  mutate(
    full_cv_r2 = cv_r2_vs_mean[DurationAnalysisSet == "full"][1],
    full_pearson_r = pearson_r[DurationAnalysisSet == "full"][1],
    full_rmse = rmse[DurationAnalysisSet == "full"][1],
    delta_cv_r2_vs_full = cv_r2_vs_mean - full_cv_r2,
    delta_pearson_r_vs_full = pearson_r - full_pearson_r,
    delta_rmse_vs_full = rmse - full_rmse
  ) %>%
  ungroup() %>%
  select(-full_cv_r2, -full_pearson_r, -full_rmse)

baseline_rmse <- ladder_performance %>% filter(Model == "Mean only") %>% pull(rmse)
movement_rmse <- ladder_performance %>% filter(Model == "Movement only") %>% pull(rmse)
combined_rmse <- ladder_performance %>% filter(Model == "Movement + Entropy ACF1") %>% pull(rmse)

incremental_summary <- ladder_performance %>%
  mutate(
    delta_rmse_vs_mean_only = rmse - baseline_rmse,
    delta_rmse_vs_movement_only = rmse - movement_rmse,
    relative_rmse_vs_movement_only = rmse / movement_rmse,
    Interpretation = case_when(
      Model == "Movement + Entropy ACF1" & rmse < movement_rmse ~ "Entropy ACF1 improves prediction beyond movement",
      Model == "Movement + Entropy ACF1" & rmse >= movement_rmse ~ "No incremental prediction beyond movement",
      TRUE ~ NA_character_
    )
  )

write_table(ladder_predictions, file.path(output_dir, "tables", "model_ladder_loo_predictions.csv"))
write_table(ladder_coefficients, file.path(output_dir, "tables", "model_ladder_loo_coefficients.csv"))
write_table(ladder_performance, file.path(output_dir, "tables", "model_ladder_performance.csv"))
write_table(incremental_summary, file.path(output_dir, "tables", "model_ladder_incremental_summary.csv"))
write_table(ladder_predictions_duration_sensitivity, file.path(output_dir, "tables", "model_ladder_loo_predictions_duration_sensitivity.csv"))
write_table(ladder_coefficients_duration_sensitivity, file.path(output_dir, "tables", "model_ladder_loo_coefficients_duration_sensitivity.csv"))
write_table(ladder_performance_duration_sensitivity, file.path(output_dir, "tables", "model_ladder_performance_duration_sensitivity.csv"))

# ------------------------------------------------
# GROUPED K-FOLD COMPANION: BEHAVIOR-ONLY PRIMARY CLAIM
# ------------------------------------------------

make_grouped_folds <- function(dat, k = 5, repeats = 100, group_col = "AnimalNum", seed = 123) {
  set.seed(seed)
  ids <- unique(dat[[group_col]])
  k <- min(k, length(ids))
  map_dfr(seq_len(repeats), function(rep_i) {
    shuffled <- sample(ids)
    fold_id <- rep(seq_len(k), length.out = length(shuffled))
    tibble(Repeat = rep_i, Fold = fold_id, !!group_col := shuffled)
  })
}

kfold_lm_predict <- function(dat, predictors, model_name, fold_map) {
  pred_rows <- vector("list", nrow(fold_map))
  for (i in seq_len(nrow(fold_map))) {
    held_out_id <- fold_map$AnimalNum[i]
    repeat_i <- fold_map$Repeat[i]
    fold_i <- fold_map$Fold[i]
    train <- dat %>% filter(AnimalNum != held_out_id)
    test <- dat %>% filter(AnimalNum == held_out_id)
    fit <- try(lm(make_formula("outcome", predictors), data = train), silent = TRUE)
    pred <- if (inherits(fit, "try-error")) NA_real_ else tryCatch(as.numeric(predict(fit, newdata = test)), error = function(e) NA_real_)
    pred_rows[[i]] <- test %>%
      transmute(
        Repeat = repeat_i,
        Fold = fold_i,
        AnimalNum, Group, Sex,
        observed = outcome,
        predicted = pred,
        Model = model_name
      )
  }
  bind_rows(pred_rows)
}

summarise_repeated_cv <- function(pred_tbl, analysis_set) {
  if (nrow(pred_tbl) == 0) return(tibble())
  pred_tbl %>%
    group_by(Model, Repeat) %>%
    group_modify(~prediction_metrics(.x$observed, .x$predicted)) %>%
    ungroup() %>%
    group_by(Model) %>%
    summarise(
      n_repeats = n_distinct(Repeat),
      mean_pearson_r = mean(pearson_r, na.rm = TRUE),
      median_pearson_r = median(pearson_r, na.rm = TRUE),
      mean_spearman_rho = mean(spearman_rho, na.rm = TRUE),
      mean_rmse = mean(rmse, na.rm = TRUE),
      mean_mae = mean(mae, na.rm = TRUE),
      mean_cv_r2 = mean(cv_r2_vs_mean, na.rm = TRUE),
      cv_r2_ci_low = quantile(cv_r2_vs_mean, 0.025, na.rm = TRUE, names = FALSE),
      cv_r2_ci_high = quantile(cv_r2_vs_mean, 0.975, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    ) %>%
    mutate(
      BinLevel = bin_level,
      Outcome = outcome_col,
      CVScheme = "repeated_grouped_kfold_leave_animals_intact",
      DurationAnalysisSet = analysis_set
    )
}

behavior_predictors <- c("Movement_mean", "Movement_rmssd", "Entropy_acf1")
compact_behavior_predictors <- c(
  "Movement_mean", "Movement_rmssd", "Movement_acf1",
  "Entropy_mean", "Entropy_rmssd", "Entropy_acf1",
  "Proximity_mean"
)
behavior_only_specs <- list(
  "Behavior-only: mean" = character(0),
  "Behavior-only: movement mean" = "Movement_mean",
  "Behavior-only: movement + entropy persistence" = c("Movement_mean", "Entropy_acf1"),
  "Behavior-only: primary family" = behavior_predictors,
  "Behavior-only: compact dynamics" = compact_behavior_predictors
) %>%
  map(~.x[.x %in% names(model_dat)])

behavior_group_specs <- map(behavior_only_specs, ~unique(c(.x, candidate_covars))) %>%
  set_names(str_replace(names(behavior_only_specs), "Behavior-only", "Behavior + group/sex"))

cv_specs <- c(behavior_only_specs, behavior_group_specs)
cv_predictors <- unique(unlist(cv_specs))
cv_predictors <- cv_predictors[cv_predictors %in% names(model_dat) & sapply(model_dat[cv_predictors], is.numeric)]
cv_model_dat <- impute_numeric(model_dat, cv_predictors)
fold_map <- make_grouped_folds(cv_model_dat, k = 5, repeats = 100, seed = 321)
repeated_cv_predictions <- imap_dfr(cv_specs, ~kfold_lm_predict(cv_model_dat, .x, .y, fold_map))
repeated_cv_performance <- summarise_repeated_cv(repeated_cv_predictions, "full")

repeated_cv_no_short_predictions <- tibble()
repeated_cv_no_short_performance <- tibble()
if (nrow(no_short_model_dat) >= 8) {
  no_short_cv_dat <- impute_numeric(no_short_model_dat, cv_predictors)
  no_short_fold_map <- make_grouped_folds(no_short_cv_dat, k = 5, repeats = 100, seed = 322)
  repeated_cv_no_short_predictions <- imap_dfr(cv_specs, ~kfold_lm_predict(no_short_cv_dat, .x, .y, no_short_fold_map))
  repeated_cv_no_short_performance <- summarise_repeated_cv(repeated_cv_no_short_predictions, "excluding_short_duration")
}

repeated_cv_performance_all <- bind_rows(repeated_cv_performance, repeated_cv_no_short_performance) %>%
  group_by(Model) %>%
  mutate(
    full_mean_cv_r2 = mean_cv_r2[DurationAnalysisSet == "full"][1],
    delta_mean_cv_r2_vs_full = mean_cv_r2 - full_mean_cv_r2
  ) %>%
  ungroup() %>%
  select(-full_mean_cv_r2)

prediction_interpretation_constraints <- tibble(
  Constraint = c(
    "Primary evidence",
    "Group-label circularity",
    "Behavior + group models",
    "Cross-validation unit",
    "Permutation testing",
    "Duration robustness",
    "Allowed language"
  ),
  Interpretation = c(
    "Use behavior-only models to support early behavior predicting later stress burden.",
    "RES/SUS group labels may be derived from CombZ, so group terms should not be treated as independent prospective predictors.",
    "Use behavior + group/sex models as descriptive adjustment/sensitivity analyses, not as the central claim.",
    "Grouped folds keep all observations from an animal together; the animal is the biological unit.",
    "Permutation p-values test full-pipeline prediction strength for the final observed predictions.",
    "Main-text claims require consistent full-data and excluding-short-duration performance.",
    "Use predictive/associative wording; avoid causal or biomarker language without external validation."
  ),
  ManuscriptUse = c(
    "Main Results",
    "Methods/Limitations",
    "Supplementary",
    "Methods",
    "Methods/Statistics",
    "Reviewer robustness",
    "Discussion"
  )
)

write_table(repeated_cv_predictions, file.path(output_dir, "tables", "model_ladder_repeated_grouped_kfold_predictions.csv"))
write_table(repeated_cv_performance_all, file.path(output_dir, "tables", "model_ladder_repeated_grouped_kfold_performance.csv"))
write_table(prediction_interpretation_constraints, file.path(output_dir, "tables", "prediction_interpretation_constraints.csv"))

# ------------------------------------------------
# FIGURES
# ------------------------------------------------

p_primary <- model_dat %>%
  select(AnimalNum, Group, Sex, outcome, all_of(primary_features)) %>%
  pivot_longer(cols = all_of(primary_features), names_to = "Feature", values_to = "Value") %>%
  ggplot(aes(Value, outcome, colour = Group, fill = Group)) +
  geom_point(size = 1.7, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.45, alpha = 0.12, colour = "grey20") +
  facet_grid(Sex ~ Feature, scales = "free_x") +
  labs(
    title = "Early movement and entropy persistence predict later stress burden",
    subtitle = paste0("First ", early_window_hours, " h active phase after first cage change; ", bin_level),
    x = "Early feature value",
    y = outcome_col
  ) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  make_publication_theme(base_size = 6)

save_plot_svg_pdf(p_primary, file.path(output_dir, "figures", "publication", "primary_movement_entropyacf1_vs_combz"), width = 170, height = 105)

p_ladder <- ladder_performance %>%
  mutate(Model = factor(Model, levels = rev(Model[order(cv_r2_vs_mean)]))) %>%
  ggplot(aes(cv_r2_vs_mean, Model)) +
  geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey55") +
  geom_col(width = 0.65, fill = "grey70", colour = "grey20", linewidth = 0.2) +
  labs(
    title = "Incremental predictive value of early behavioral dynamics",
    subtitle = "Leave-one-animal-out linear model ladder",
    x = "Cross-validated R2 versus mean-only baseline",
    y = NULL
  ) +
  make_publication_theme(base_size = 7) +
  theme(panel.grid.major.y = element_blank())

save_plot_svg_pdf(p_ladder, file.path(output_dir, "figures", "publication", "model_ladder_cv_r2"), width = 130, height = 80)

best_model <- ladder_performance %>% slice_max(cv_r2_vs_mean, n = 1, with_ties = FALSE) %>% pull(Model)
best_pred <- ladder_predictions %>% filter(Model == best_model)
best_perf <- ladder_performance %>% filter(Model == best_model)

p_pred <- best_pred %>%
  ggplot(aes(observed, predicted, colour = Group, fill = Group)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey45") +
  geom_point(size = 1.8, alpha = 0.85) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.5, alpha = 0.12, colour = "grey20") +
  facet_grid(. ~ Sex) +
  labs(
    title = paste0("Best early prediction model: ", best_model),
    subtitle = paste0("LOO r=", round(best_perf$pearson_r, 2), ", CV R2=", round(best_perf$cv_r2_vs_mean, 2), ", ", format_p(best_perf$prediction_permutation_p)),
    x = paste0("Observed ", outcome_col),
    y = paste0("Predicted ", outcome_col)
  ) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  make_publication_theme(base_size = 7)

save_plot_svg_pdf(p_pred, file.path(output_dir, "figures", "publication", "best_model_observed_vs_predicted"), width = 120, height = 80)

# ------------------------------------------------
# TEXT SUMMARY FOR RESULTS WRITING
# ------------------------------------------------

combined_row <- incremental_summary %>% filter(Model == "Movement + Entropy ACF1")
movement_row <- incremental_summary %>% filter(Model == "Movement only")
entropy_row <- primary_assoc %>% filter(feature == "Entropy_acf1")
movement_assoc_row <- primary_assoc %>% filter(feature == "Movement_mean")

results_summary <- tibble(
  Result = c(
    "Primary prospective window",
    "Movement association",
    "Entropy ACF1 association",
    "Incremental model comparison"
  ),
  Text = c(
    paste0("Features were extracted from the first ", early_window_hours, " h of the active phase after the first cage change using ", bin_level, " bins."),
    if (nrow(movement_assoc_row) > 0) paste0("Early movement correlated with later ", outcome_col, " (Spearman rho=", round(movement_assoc_row$spearman_rho, 3), ", BH ", format_p(movement_assoc_row$spearman_p_bh), ").") else "Movement association unavailable.",
    if (nrow(entropy_row) > 0) paste0("Early entropy ACF1 correlated with later ", outcome_col, " (Spearman rho=", round(entropy_row$spearman_rho, 3), ", BH ", format_p(entropy_row$spearman_p_bh), "; bootstrap CI ", round(entropy_row$spearman_boot_ci_low, 3), " to ", round(entropy_row$spearman_boot_ci_high, 3), ").") else "Entropy ACF1 association unavailable.",
    if (nrow(combined_row) > 0 && nrow(movement_row) > 0) paste0("Adding entropy ACF1 to movement changed LOO RMSE from ", round(movement_row$rmse, 3), " to ", round(combined_row$rmse, 3), " and CV R2 from ", round(movement_row$cv_r2_vs_mean, 3), " to ", round(combined_row$cv_r2_vs_mean, 3), ".") else "Incremental comparison unavailable."
  )
)

write_table(results_summary, file.path(output_dir, "tables", "results_summary_text.csv"))

message("Early prediction model ladder complete. Tables and publication figures written to: ", output_dir)
