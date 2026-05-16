# ================================================================
# Nonlinear Systems Dynamics Analysis
# MMMSociability
# ================================================================
# Adds a modern nonlinear systems-neuroscience layer for RFID/home-cage
# stress-resilience data: nonlinear manifolds, latent trajectories,
# unsupervised state transitions, recurrence maps, and nonlinear endpoint
# response curves.
#
# Run after:
#   Analysis/03_build_multiscale_behavior_metrics.R
# Optional but useful:
#   Analysis/08_early_prediction_models.R
#   Analysis/12_systems_neuroscience_summary.R
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

# -----------------------------
# User settings
# -----------------------------
project_root <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID"
repo_root <- "C:/Users/topohl/Documents/GitHub/MMMSociability"
bin_level <- "5min_based"

input_file <- file.path(project_root, "analysis_ready/03_derived_metrics", bin_level, "all_behavior_metrics.csv")
output_dir <- file.path(project_root, "analysis_ready/13_nonlinear_systems_dynamics", bin_level)
data_dir <- file.path(output_dir, "derived_data")
stats_dir <- file.path(output_dir, "statistical_results")
figure_dir <- file.path(output_dir, "figures/manuscript_panels")
overview_figure_dir <- file.path(output_dir, "figures/overview")
first_active_figure_dir <- file.path(output_dir, "figures/first_cage_change_active_12h")
longitudinal_figure_dir <- file.path(output_dir, "figures/phase_sex_cage_change_trajectories")

endpoint_file <- NULL
endpoint_cols <- c("CombZ", "stress_z_score", "SucrosePreference", "Corticosterone", "DeltaCorticosterone")
primary_outcome <- "CombZ"

early_phase_pattern <- "active|dark|night"
early_window_hours <- 12
n_behavioral_states <- 4
state_labels <- c("Inactive", "Explore", "Social", "Burst")
random_seed <- 123

group_levels <- c("CON", "RES", "SUS")
group_colors <- c("CON" = "#3d3b6e", "RES" = "#C6C3BB", "SUS" = "#e63947")
sex_levels <- c("Female", "Male")

# -----------------------------
# Helpers
# -----------------------------
source_if_exists <- function(path) if (file.exists(path)) source(path)
source_if_exists(file.path(repo_root, "Functions/behavioral_dynamics_helpers.R"))
source_if_exists(file.path(repo_root, "Functions/behavioral_dynamics_stats_helpers.R"))

if (!exists("ensure_dir")) ensure_dir <- function(path) { if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE); invisible(path) }
if (!exists("write_table")) write_table <- function(x, path) { ensure_dir(dirname(path)); readr::write_csv(x, path); invisible(path) }
if (!exists("save_plot_svg_pdf")) save_plot_svg_pdf <- function(plot, filename_base, width = 85, height = 65, units = "mm") { ensure_dir(dirname(filename_base)); ggsave(paste0(filename_base, ".svg"), plot, width = width, height = height, units = units); ggsave(paste0(filename_base, ".pdf"), plot, width = width, height = height, units = units); invisible(filename_base) }
if (!exists("first_existing_col")) first_existing_col <- function(dat, candidates, required = TRUE, label = "column") { hit <- candidates[candidates %in% names(dat)][1]; if (is.na(hit) && required) stop("Missing ", label, call. = FALSE); hit }
if (!exists("make_nature_theme")) make_nature_theme <- function(base_size = 7) theme_classic(base_size = base_size) + theme(axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.2), strip.background = element_blank(), strip.text = element_text(face = "bold"), legend.title = element_blank(), legend.position = "top", plot.title = element_text(face = "bold", hjust = 0), plot.subtitle = element_text(hjust = 0, colour = "grey35"))

safe_num <- function(x) suppressWarnings(as.numeric(x))
safe_scale <- function(x) { s <- sd(x, na.rm = TRUE); m <- mean(x, na.rm = TRUE); if (!is.finite(s) || s == 0) return(rep(0, length(x))); (x - m) / s }
safe_cor <- function(x, y, method = "spearman") { ok <- is.finite(x) & is.finite(y); if (sum(ok) < 4 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_); suppressWarnings(cor(x[ok], y[ok], method = method)) }
calc_rmssd <- function(x) { x <- x[is.finite(x)]; if (length(x) < 3) return(NA_real_); sqrt(mean(diff(x)^2, na.rm = TRUE)) }
calc_acf1 <- function(x) { x <- x[is.finite(x)]; if (length(x) < 4 || sd(x) == 0) return(NA_real_); suppressWarnings(acf(x, plot = FALSE, lag.max = 1)$acf[2]) }

panel_theme <- function(base_size = 7) make_nature_theme(base_size) + theme(panel.grid = element_blank(), legend.position = "right")

p_label <- function(p) {
  case_when(
    is.na(p) ~ "p=NA",
    p < 0.001 ~ "p<0.001",
    TRUE ~ paste0("p=", formatC(p, format = "f", digits = 3))
  )
}

p_stars <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    p < 0.10 ~ ".",
    TRUE ~ "ns"
  )
}

classify_phase <- function(x) {
  lx <- str_to_lower(as.character(x))
  case_when(
    str_detect(lx, "active|dark|night") & !str_detect(lx, "inactive") ~ "Active",
    str_detect(lx, "inactive|light|day") ~ "Inactive",
    TRUE ~ as.character(x)
  )
}

standard_error <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x) / sqrt(length(x))
}

write_professional_table <- function(x, filename, subdir = data_dir) {
  write_table(x, file.path(subdir, filename))
}

save_manuscript_plot <- function(plot, filename, width = 120, height = 80) {
  save_plot_svg_pdf(plot, file.path(figure_dir, filename), width = width, height = height)
}
save_first_active_plot <- function(plot, filename, width = 120, height = 80) {
  save_plot_svg_pdf(plot, file.path(first_active_figure_dir, filename), width = width, height = height)
}
save_longitudinal_plot <- function(plot, filename, width = 120, height = 80) {
  save_plot_svg_pdf(plot, file.path(longitudinal_figure_dir, filename), width = width, height = height)
}
save_overview_plot <- function(plot, filename, width = 120, height = 80) {
  save_plot_svg_pdf(plot, file.path(overview_figure_dir, filename), width = width, height = height)
}

make_group_stats <- function(dat, value_cols, by_cols) {
  by_cols <- intersect(by_cols, names(dat))
  value_cols <- intersect(value_cols, names(dat))
  if (nrow(dat) == 0 || length(value_cols) == 0) return(tibble())
  make_dynamics_group_contrasts(dat, value_cols = value_cols, by_cols = by_cols) %>%
    mutate(
      ReportingP = p.adjust_bh_family,
      ReportingSignificance = p_stars(ReportingP),
      ReportingCorrection = "BH within outcome/facet family",
      Significant = ReportingP < 0.05,
      p_label = p_label(p.adjust_bh_family),
      p_stars = p_stars(p.adjust_bh_family),
      effect_label = paste0(contrast, ": d=", round(cohen_d, 2), ", ", p_label)
    )
}

make_group_summary <- function(dat, value_cols, group_cols) {
  make_dynamics_group_summary(dat, value_cols = value_cols, group_cols = group_cols) %>%
    mutate(mean_sem = paste0(round(mean, 3), " +/- ", round(sem, 3)))
}

make_pairwise_brackets <- function(data, contrast_tbl, outcome, facet_cols, value_col = NULL, include_trends = FALSE) {
  if (is.null(value_col)) value_col <- outcome
  empty <- tibble(contrast = character(), x_start = numeric(), x_end = numeric(), y = numeric(), y_tip = numeric(), label = character())
  for (nm in facet_cols) empty[[nm]] <- data[[nm]][0]
  if (nrow(contrast_tbl) == 0 || !value_col %in% names(data)) return(empty)
  threshold <- if (include_trends) 0.10 else 0.05
  brackets <- contrast_tbl %>%
    filter(
      Outcome == outcome,
      contrast %in% c("RES-CON", "SUS-CON", "SUS-RES"),
      !is.na(ReportingP),
      ReportingP < threshold
    ) %>%
    mutate(
      x_start = match(group_ref, group_levels),
      x_end = match(group_comp, group_levels),
      label = p_label(ReportingP)
    ) %>%
    filter(is.finite(x_start), is.finite(x_end))
  if (nrow(brackets) == 0) return(empty)
  y_tbl <- data %>%
    group_by(across(all_of(facet_cols))) %>%
    summarise(y_max = max(.data[[value_col]], na.rm = TRUE), y_min = min(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
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
    mutate(y = y_max + bracket_rank * y_step, y_tip = y - y_tip_size) %>%
    select(any_of(facet_cols), contrast, x_start, x_end, y, y_tip, label, ReportingP, ReportingSignificance)
}

add_pairwise_brackets <- function(plot, bracket_tbl, text_size = 1.8) {
  if (nrow(bracket_tbl) == 0) return(plot)
  plot +
    geom_segment(data = bracket_tbl, aes(x = x_start, xend = x_end, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.25) +
    geom_segment(data = bracket_tbl, aes(x = x_start, xend = x_start, y = y, yend = y_tip), inherit.aes = FALSE, linewidth = 0.25) +
    geom_segment(data = bracket_tbl, aes(x = x_end, xend = x_end, y = y, yend = y_tip), inherit.aes = FALSE, linewidth = 0.25) +
    geom_text(data = bracket_tbl, aes(x = (x_start + x_end) / 2, y = y, label = label), inherit.aes = FALSE, vjust = -0.25, size = text_size)
}

make_longitudinal_sig_labels <- function(summary_tbl, contrast_tbl, outcomes) {
  if (nrow(contrast_tbl) == 0) return(tibble())
  contrast_tbl %>%
    filter(Outcome %in% outcomes, !is.na(ReportingP), ReportingP < 0.10, contrast %in% c("RES-CON", "SUS-CON", "SUS-RES")) %>%
    group_by(Sex, PhaseClass, CageChangeLabel, CageChangeIndex, Outcome) %>%
    summarise(label = paste0(contrast[which.min(ReportingP)], " ", p_stars(min(ReportingP, na.rm = TRUE))), .groups = "drop") %>%
    left_join(
      summary_tbl %>%
        group_by(Sex, PhaseClass, CageChangeLabel, CageChangeIndex, Outcome) %>%
        summarise(y = max(ymax, na.rm = TRUE), y_min = min(ymin, na.rm = TRUE), .groups = "drop") %>%
        mutate(y = y + if_else(is.finite(y - y_min) & y > y_min, 0.12 * (y - y_min), 0.1)),
      by = c("Sex", "PhaseClass", "CageChangeLabel", "CageChangeIndex", "Outcome")
    )
}

fit_time_model <- function(dat, value_col) {
  dat <- dat %>% filter(is.finite(.data[[value_col]]), !is.na(Group), !is.na(Sex), !is.na(PhaseClass), !is.na(CageChangeIndex))
  if (nrow(dat) < 8 || n_distinct(dat$Group) < 2 || n_distinct(dat$CageChangeIndex) < 2) {
    return(tibble(outcome = value_col, term = NA_character_, estimate = NA_real_, statistic = NA_real_, p.value = NA_real_, model = "skipped_low_information"))
  }
  dat <- dat %>% mutate(CageChangeIndex = as.numeric(CageChangeIndex), Group = droplevels(factor(Group)), Sex = droplevels(factor(Sex)), PhaseClass = droplevels(factor(PhaseClass)))
  if (requireNamespace("lme4", quietly = TRUE) && requireNamespace("lmerTest", quietly = TRUE)) {
    form <- as.formula(paste0(value_col, " ~ Group * CageChangeIndex * PhaseClass * Sex + (1|AnimalNum)"))
    fit <- try(lmerTest::lmer(form, data = dat), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      return(as_tibble(coef(summary(fit)), rownames = "term") %>%
               transmute(outcome = value_col, term, estimate = Estimate, statistic = `t value`, p.value = `Pr(>|t|)`, model = "lmer_random_intercept"))
    }
  }
  form <- as.formula(paste0(value_col, " ~ Group * CageChangeIndex * PhaseClass * Sex + AnimalNum"))
  fit <- try(lm(form, data = dat), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(tibble(outcome = value_col, term = NA_character_, estimate = NA_real_, statistic = NA_real_, p.value = NA_real_, model = "failed"))
  }
  as_tibble(coef(summary(fit)), rownames = "term") %>%
    transmute(outcome = value_col, term, estimate = Estimate, statistic = `t value`, p.value = `Pr(>|t|)`, model = "lm_animal_fixed_effect")
}

read_any <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") return(read_csv(path, show_col_types = FALSE))
  if (ext %in% c("tsv", "txt")) return(read_tsv(path, show_col_types = FALSE))
  if (ext == "rds") return(readRDS(path))
  if (ext %in% c("xlsx", "xls")) { if (!requireNamespace("readxl", quietly = TRUE)) stop("Install readxl"); return(readxl::read_excel(path)) }
  NULL
}

standardize_input <- function(dat) {
  animal_col <- first_existing_col(dat, c("AnimalNum", "Animal", "AnimalID", "MouseID", "Mouse", "ID", "RFID"), TRUE, "animal ID")
  time_col <- first_existing_col(dat, c("TimeIndex", "BinStart", "HalfHourElapsed", "HalfHourWithinCC0", "Time", "DateTime"), TRUE, "time")
  group_col <- first_existing_col(dat, c("Group", "Phenotype", "Condition", "Treatment", "StressGroup"), FALSE, "group")
  sex_col <- first_existing_col(dat, c("Sex", "sex"), FALSE, "sex")
  phase_col <- first_existing_col(dat, c("Phase", "phase", "LightDark", "DayNight"), FALSE, "phase")
  cage_col <- first_existing_col(dat, c("CageChange", "CC", "CageChangeNum", "Regrouping", "Batch", "Cage"), FALSE, "cage")
  movement_col <- first_existing_col(dat, c("Movement", "movement", "Distance", "Activity"), TRUE, "movement")
  entropy_col <- first_existing_col(dat, c("Entropy", "entropy", "ShannonEntropy", "PositionEntropy"), TRUE, "entropy")
  proximity_col <- first_existing_col(dat, c("ProximityFraction", "Proximity", "proximity", "MeanProximity", "SocialProximity"), TRUE, "proximity")

  dat %>%
    mutate(
      AnimalNum = .data[[animal_col]],
      TimeIndex = .data[[time_col]],
      Group = if (!is.na(group_col)) as.character(.data[[group_col]]) else "All",
      Sex = if (!is.na(sex_col)) as.character(.data[[sex_col]]) else "All",
      Phase = if (!is.na(phase_col)) as.character(.data[[phase_col]]) else "All",
      CageChange = if (!is.na(cage_col)) as.character(.data[[cage_col]]) else "All",
      Movement = safe_num(.data[[movement_col]]),
      Entropy = safe_num(.data[[entropy_col]]),
      Proximity = safe_num(.data[[proximity_col]])
    ) %>%
    mutate(
      Group = factor(Group, levels = unique(c(group_levels, sort(unique(Group))))),
      Sex = factor(Sex, levels = unique(c(sex_levels, sort(unique(Sex)))))
    ) %>%
    filter(!is.na(AnimalNum), !is.na(TimeIndex)) %>%
    arrange(AnimalNum, CageChange, Phase, TimeIndex)
}

impute_scale <- function(dat, cols) {
  if (length(cols) == 0) stop("No usable nonlinear features after QC; inspect statistical_results/animal_level_nonlinear_feature_quality_control.csv.", call. = FALSE)
  x <- dat %>%
    select(all_of(cols)) %>%
    mutate(across(everything(), ~ {
      med <- median(.x, na.rm = TRUE)
      if (!is.finite(med)) med <- 0
      replace_na(.x, med)
    })) %>%
    as.matrix()

  z <- as.matrix(scale(x))
  bad_cols <- !is.finite(colSums(z))
  if (any(bad_cols)) z[, bad_cols] <- 0
  z[!is.finite(z)] <- 0
  z
}

# -----------------------------
# Load data
# -----------------------------
ensure_dir(output_dir)
output_dirs <- analysis_output_dirs(output_dir)
ensure_dir(file.path(output_dir, "tables"))
ensure_dir(file.path(output_dir, "stats_tables"))
ensure_dir(file.path(output_dir, "figures/publication_panels"))
ensure_dir(data_dir)
ensure_dir(stats_dir)
ensure_dir(figure_dir)
ensure_dir(overview_figure_dir)
ensure_dir(first_active_figure_dir)
ensure_dir(longitudinal_figure_dir)
write_output_manifest(
  output_dir,
  script_name = "13_nonlinear_systems_dynamics.R",
  analysis_name = "nonlinear systems dynamics",
  primary_tables = c(
    "derived_data/early_active_group_recurrence_similarity_maps.csv",
    "derived_data/distance_to_sex_matched_control_manifold.csv",
    "statistical_results",
    "tables"
  ),
  primary_figures = c(
    "figures/manuscript_panels",
    "figures/overview",
    "figures/first_cage_change_active_12h",
    "figures/phase_sex_cage_change_trajectories"
  ),
  notes = c("Nonlinear visualizations are exploratory; script 12 imports only interpretable robust features.")
)

raw <- read_any(input_file)
if (is.null(raw)) stop("Missing input file. Run 03_build_multiscale_behavior_metrics.R first.", call. = FALSE)
behav <- standardize_input(raw) %>%
  mutate(
    CageChangeLabel = as.character(CageChange),
    PhaseLabel = as.character(Phase),
    PhaseClass = factor(classify_phase(Phase), levels = unique(c("Active", "Inactive", classify_phase(Phase)))),
    CageChangeIndex = suppressWarnings(as.integer(str_extract(CageChangeLabel, "\\d+"))),
    CageChangeIndex = if_else(is.na(CageChangeIndex), as.integer(factor(CageChangeLabel, levels = unique(CageChangeLabel))), CageChangeIndex),
    CageChangeLabel = factor(CageChangeLabel, levels = unique(CageChangeLabel[order(CageChangeIndex)])),
    Movement_z_global = safe_scale(Movement),
    Entropy_z_global = safe_scale(Entropy),
    Proximity_z_global = safe_scale(Proximity)
  )

endpoint_dat <- read_any(endpoint_file)
if (!is.null(endpoint_dat)) {
  animal_col <- first_existing_col(endpoint_dat, c("AnimalNum", "Animal", "AnimalID", "MouseID", "Mouse", "ID", "RFID"), TRUE, "endpoint animal ID")
  keep <- intersect(endpoint_cols, names(endpoint_dat))
  endpoint_dat <- endpoint_dat %>% transmute(AnimalNum = .data[[animal_col]], across(all_of(keep), safe_num))
  behav <- left_join(behav, endpoint_dat, by = "AnimalNum")
}

available_outcomes <- intersect(endpoint_cols, names(behav))
available_outcomes <- available_outcomes[map_lgl(available_outcomes, ~ any(is.finite(safe_num(behav[[.x]]))))]
outcome_to_use <- if (primary_outcome %in% available_outcomes) primary_outcome else ifelse(length(available_outcomes) > 0, available_outcomes[1], NA_character_)
bins_per_hour <- case_when(str_detect(bin_level, "10sec") ~ 360, str_detect(bin_level, "5min") ~ 12, str_detect(bin_level, "10min") ~ 6, str_detect(bin_level, "30min") ~ 2, TRUE ~ 12)
max_early_bins <- early_window_hours * bins_per_hour

# -----------------------------
# First-window and longitudinal phase dynamics
# -----------------------------
first_cage_change <- min(behav$CageChangeIndex, na.rm = TRUE)
first_active_phase <- behav %>%
  filter(CageChangeIndex == first_cage_change, PhaseClass == "Active") %>%
  distinct(PhaseLabel) %>%
  slice_head(n = 1) %>%
  pull(PhaseLabel)

if (length(first_active_phase) == 0) {
  first_active_phase <- behav %>%
    filter(CageChangeIndex == first_cage_change) %>%
    distinct(PhaseLabel) %>%
    slice_head(n = 1) %>%
    pull(PhaseLabel)
}
first_active_phase <- first_active_phase[1]

summarise_window_features <- function(dat) {
  dat %>%
    summarise(
      Movement_mean = mean(Movement, na.rm = TRUE),
      Movement_rmssd = calc_rmssd(Movement),
      Entropy_mean = mean(Entropy, na.rm = TRUE),
      Entropy_rmssd = calc_rmssd(Entropy),
      Proximity_mean = mean(Proximity, na.rm = TRUE),
      Proximity_rmssd = calc_rmssd(Proximity),
      SocialWithdrawal = mean(Movement_z_global, na.rm = TRUE) - mean(Proximity_z_global, na.rm = TRUE),
      BehavioralInstability = mean(c(calc_rmssd(Movement_z_global), calc_rmssd(Entropy_z_global), calc_rmssd(Proximity_z_global)), na.rm = TRUE),
      n_bins = n(),
      .groups = "drop"
    )
}

first_active_12h <- behav %>%
  filter(CageChangeIndex == first_cage_change, PhaseLabel == first_active_phase) %>%
  group_by(AnimalNum, Group, Sex, CageChangeLabel, CageChangeIndex, PhaseLabel, PhaseClass) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  mutate(WindowBin = row_number(), WindowHour = (WindowBin - 1) / bins_per_hour) %>%
  filter(WindowBin <= max_early_bins) %>%
  summarise_window_features() %>%
  mutate(Window = "First 12h of first active phase, first cage change")

phase_12h_features <- behav %>%
  filter(PhaseClass %in% c("Active", "Inactive")) %>%
  group_by(AnimalNum, Group, Sex, CageChangeLabel, CageChangeIndex, PhaseLabel, PhaseClass) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  mutate(WindowBin = row_number(), WindowHour = (WindowBin - 1) / bins_per_hour) %>%
  filter(WindowBin <= max_early_bins) %>%
  summarise_window_features() %>%
  mutate(Window = "First 12h of each phase")

window_metrics <- c("Movement_mean", "Entropy_mean", "Proximity_mean", "SocialWithdrawal", "BehavioralInstability")

first_active_summary <- make_group_summary(first_active_12h, window_metrics, c("Window", "Sex", "Group"))
first_active_contrasts <- make_group_stats(first_active_12h, window_metrics, c("Window", "Sex"))
first_active_significant <- first_active_contrasts %>%
  filter(status == "tested", is.finite(ReportingP), ReportingP < 0.10) %>%
  arrange(ReportingP, Outcome, Sex, contrast)

phase_summary <- make_group_summary(phase_12h_features, window_metrics, c("Window", "CageChangeIndex", "CageChangeLabel", "PhaseClass", "Sex", "Group"))
phase_contrasts <- make_group_stats(phase_12h_features, window_metrics, c("Window", "CageChangeIndex", "CageChangeLabel", "PhaseClass", "Sex"))
phase_significant <- phase_contrasts %>%
  filter(status == "tested", is.finite(ReportingP), ReportingP < 0.10) %>%
  arrange(ReportingP, Outcome, PhaseClass, Sex, CageChangeIndex, contrast)
phase_mixed_models <- map_dfr(window_metrics, ~ fit_time_model(phase_12h_features, .x)) %>%
  group_by(outcome) %>%
  mutate(p.adjust_bh_outcome = p.adjust(p.value, method = "BH"), p_label = p_label(p.adjust_bh_outcome), p_stars = p_stars(p.adjust_bh_outcome)) %>%
  ungroup()

write_professional_table(first_active_12h, "first_active_phase_cage_change_1_first_12h_animal_features.csv")
write_professional_table(phase_12h_features, "phase_specific_first_12h_by_cage_change_animal_features.csv")
write_professional_table(first_active_summary, "first_active_phase_cage_change_1_first_12h_group_summary.csv", stats_dir)
write_professional_table(first_active_contrasts, "first_active_phase_cage_change_1_first_12h_group_contrasts.csv", stats_dir)
write_professional_table(first_active_significant, "first_active_phase_cage_change_1_first_12h_significant_and_trend_contrasts.csv", stats_dir)
write_professional_table(phase_summary, "phase_specific_first_12h_by_cage_change_group_summary.csv", stats_dir)
write_professional_table(phase_contrasts, "phase_specific_first_12h_by_cage_change_group_contrasts.csv", stats_dir)
write_professional_table(phase_significant, "phase_specific_first_12h_by_cage_change_significant_and_trend_contrasts.csv", stats_dir)
write_professional_table(phase_mixed_models, "phase_specific_first_12h_longitudinal_mixed_models.csv", stats_dir)

first_active_long <- first_active_12h %>%
  pivot_longer(all_of(window_metrics), names_to = "Metric", values_to = "Value") %>%
  filter(is.finite(Value)) %>%
  mutate(Metric = factor(Metric, levels = window_metrics))

first_active_brackets <- map_dfr(window_metrics, ~ make_pairwise_brackets(
  first_active_long,
  first_active_contrasts %>% mutate(Metric = factor(Outcome, levels = window_metrics)),
  outcome = .x,
  facet_cols = c("Metric", "Sex"),
  value_col = "Value",
  include_trends = TRUE
))

p_first_active <- first_active_long %>%
  ggplot(aes(Group, Value, fill = Group, colour = Group)) +
  geom_violin(alpha = 0.24, colour = NA, trim = FALSE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.25, alpha = 0.70) +
  geom_jitter(width = 0.07, size = 0.9, alpha = 0.78, show.legend = FALSE) +
  facet_grid(Metric ~ Sex, scales = "free_y") +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  labs(
    title = "A  First active-phase response during cage change 1",
    subtitle = paste0("Animal-level first ", early_window_hours, "h summaries; annotations show BH-adjusted Welch contrasts"),
    x = NULL,
    y = "Window summary"
  ) +
  panel_theme(6) +
  theme(legend.position = "none")
p_first_active <- add_pairwise_brackets(p_first_active, first_active_brackets, text_size = 1.6)
save_first_active_plot(p_first_active, "first_active_phase_cage_change_1_first_12h_group_differences", width = 150, height = 145)

phase_long <- phase_12h_features %>%
  pivot_longer(all_of(window_metrics), names_to = "Metric", values_to = "Value") %>%
  filter(is.finite(Value)) %>%
  mutate(Metric = factor(Metric, levels = window_metrics))

phase_traj <- phase_long %>%
  group_by(Metric, PhaseClass, Sex, Group, CageChangeIndex) %>%
  summarise(mean = mean(Value, na.rm = TRUE), sem = standard_error(Value), n_animals = n_distinct(AnimalNum), .groups = "drop")

phase_plot_stats <- phase_contrasts %>%
  filter(status == "tested", contrast == "SUS-CON", ReportingP < 0.10) %>%
  transmute(Metric = Outcome, PhaseClass, Sex, CageChangeIndex, label = p_stars) %>%
  left_join(phase_long %>% group_by(Metric, PhaseClass, Sex, CageChangeIndex) %>% summarise(y = max(Value, na.rm = TRUE), .groups = "drop"), by = c("Metric", "PhaseClass", "Sex", "CageChangeIndex")) %>%
  mutate(y = ifelse(is.finite(y), y, 0), y = y + 0.08 * abs(y + 1))

p_phase_time <- phase_traj %>%
  ggplot(aes(CageChangeIndex, mean, colour = Group, fill = Group, group = Group)) +
  geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.16, colour = NA) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 1.5) +
  geom_text(data = phase_plot_stats, aes(CageChangeIndex, y, label = label), inherit.aes = FALSE, size = 2.2, colour = "grey15") +
  facet_grid(Metric + PhaseClass ~ Sex, scales = "free_y") +
  scale_x_continuous(breaks = sort(unique(phase_traj$CageChangeIndex))) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  labs(
    title = "B  Phase-specific behavioral adaptation across cage changes",
    subtitle = "Mean +/- SEM per animal-level first-12h phase window; stars mark SUS-CON BH-adjusted contrasts",
    x = "Cage change",
    y = "Window summary"
  ) +
  panel_theme(6)
save_longitudinal_plot(p_phase_time, "phase_specific_first_12h_group_trajectories_by_sex_and_cage_change", width = 180, height = 190)

phase_effect_heatmap <- phase_contrasts %>%
  filter(status == "tested", contrast %in% c("SUS-CON", "SUS-RES")) %>%
  mutate(
    Metric = factor(Outcome, levels = window_metrics),
    CageChangeIndex = factor(CageChangeIndex),
    label = paste0(round(cohen_d, 2), "\n", p_stars)
  )

p_phase_effect <- phase_effect_heatmap %>%
  ggplot(aes(CageChangeIndex, PhaseClass, fill = cohen_d)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = label), size = 2.0, lineheight = 0.86) +
  facet_grid(Metric + contrast ~ Sex) +
  scale_fill_gradient2(low = "#3d3b6e", mid = "white", high = "#e63947", midpoint = 0, na.value = "grey92") +
  labs(
    title = "C  Where group differences emerge over time",
    subtitle = "Tiles show Cohen's d for primary contrasts; text shows d and BH-adjusted significance",
    x = "Cage change",
    y = "Phase",
    fill = "Cohen's d"
  ) +
  panel_theme(6)
save_longitudinal_plot(p_phase_effect, "phase_specific_first_12h_longitudinal_effect_size_heatmap", width = 170, height = 210)

# -----------------------------
# Animal-level nonlinear feature matrix
# -----------------------------
animal_features <- behav %>%
  group_by(AnimalNum, Group, Sex) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  summarise(
    Movement_mean = mean(Movement, na.rm = TRUE),
    Movement_sd = sd(Movement, na.rm = TRUE),
    Movement_rmssd = calc_rmssd(Movement),
    Movement_acf1 = calc_acf1(Movement),
    Movement_p95 = quantile(Movement, 0.95, na.rm = TRUE, names = FALSE),
    Entropy_mean = mean(Entropy, na.rm = TRUE),
    Entropy_sd = sd(Entropy, na.rm = TRUE),
    Entropy_rmssd = calc_rmssd(Entropy),
    Entropy_acf1 = calc_acf1(Entropy),
    Entropy_p95 = quantile(Entropy, 0.95, na.rm = TRUE, names = FALSE),
    Proximity_mean = mean(Proximity, na.rm = TRUE),
    Proximity_sd = sd(Proximity, na.rm = TRUE),
    Proximity_rmssd = calc_rmssd(Proximity),
    Proximity_acf1 = calc_acf1(Proximity),
    Proximity_p95 = quantile(Proximity, 0.95, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  ) %>%
  mutate(
    SocialWithdrawal = safe_scale(Movement_mean) - safe_scale(Proximity_mean),
    BehavioralInstability = rowMeans(cbind(safe_scale(Movement_rmssd), safe_scale(Entropy_rmssd), safe_scale(Proximity_rmssd)), na.rm = TRUE),
    BehavioralInertia = rowMeans(cbind(safe_scale(Movement_acf1), safe_scale(Entropy_acf1), safe_scale(Proximity_acf1)), na.rm = TRUE)
  )

if (!is.na(outcome_to_use)) {
  outcome_tbl <- behav %>% group_by(AnimalNum) %>% summarise(outcome = first(na.omit(safe_num(.data[[outcome_to_use]]))), .groups = "drop")
  animal_features <- left_join(animal_features, outcome_tbl, by = "AnimalNum")
}

feature_cols <- animal_features %>% select(where(is.numeric)) %>% names() %>% setdiff(c("AnimalNum", "outcome"))
feature_qc <- map_dfr(feature_cols, function(fc) tibble(feature = fc, n_finite = sum(is.finite(animal_features[[fc]])), missing_fraction = mean(!is.finite(animal_features[[fc]])), sd = sd(animal_features[[fc]], na.rm = TRUE)))
usable_features <- feature_qc %>% filter(n_finite >= 4, missing_fraction <= 0.5, is.finite(sd), sd > 0) %>% pull(feature)

write_professional_table(animal_features, "animal_level_nonlinear_feature_matrix.csv")
write_professional_table(feature_qc, "animal_level_nonlinear_feature_quality_control.csv", stats_dir)

x_scaled <- impute_scale(animal_features, usable_features)

# -----------------------------
# 1. Nonlinear manifolds
# -----------------------------
set.seed(random_seed)
embed <- animal_features %>% select(AnimalNum, Group, Sex, any_of("outcome"))

if (requireNamespace("uwot", quietly = TRUE) && nrow(x_scaled) >= 8) {
  nn <- max(3, min(10, floor(nrow(x_scaled) / 2)))
  um <- uwot::umap(x_scaled, n_neighbors = nn, min_dist = 0.25, metric = "euclidean", init = "spectral")
  embed <- embed %>% mutate(UMAP1 = um[, 1], UMAP2 = um[, 2])
}

# Diffusion-like geometry fallback: classical MDS on distances. If destiny is installed, use true diffusion maps.
if (requireNamespace("destiny", quietly = TRUE) && nrow(x_scaled) >= 8) {
  dm <- destiny::DiffusionMap(x_scaled, k = max(3, min(10, nrow(x_scaled) - 1)))
  ev <- destiny::eigenvectors(dm)[, 1:2, drop = FALSE]
  embed <- embed %>% mutate(Diffusion1 = ev[, 1], Diffusion2 = ev[, 2])
} else {
  md <- cmdscale(dist(x_scaled), k = 2)
  embed <- embed %>% mutate(Diffusion1 = md[, 1], Diffusion2 = md[, 2])
}

if (requireNamespace("vegan", quietly = TRUE) && nrow(x_scaled) >= 8) {
  iso <- vegan::isomap(dist(x_scaled), ndim = 2, k = max(3, min(10, nrow(x_scaled) - 1)))
  embed <- bind_cols(embed, as_tibble(iso$points[, 1:2]) %>% setNames(c("Isomap1", "Isomap2")))
}

write_professional_table(embed, "animal_level_nonlinear_embedding_scores.csv")

plot_embedding <- function(dat, x, y, title, subtitle, filename) {
  if (!all(c(x, y) %in% names(dat))) return(NULL)
  p <- ggplot(dat, aes(.data[[x]], .data[[y]], colour = Group, fill = Group, shape = Sex)) +
    stat_ellipse(aes(group = Group), geom = "polygon", alpha = 0.10, linewidth = 0.2, show.legend = FALSE) +
    geom_point(size = 2.2, alpha = 0.92, stroke = 0.25) +
    scale_colour_manual(values = group_colors, drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    labs(title = title, subtitle = subtitle, x = x, y = y) +
    panel_theme(7)
  save_manuscript_plot(p, filename, width = 95, height = 78)
  p
}

p_umap <- plot_embedding(embed, "UMAP1", "UMAP2", "D  Nonlinear behavioral manifold", "UMAP on animal-level movement, entropy, proximity and temporal dynamics", "Fig13D_umap_behavioral_manifold")
p_diff <- plot_embedding(embed, "Diffusion1", "Diffusion2", "E  Diffusion-like state geometry", "Geometry-preserving behavioral organization", "Fig13E_diffusion_behavioral_geometry")
p_iso <- plot_embedding(embed, "Isomap1", "Isomap2", "F  Curved behavioral manifold", "Geodesic embedding for nonlinear group structure", "Fig13F_isomap_behavioral_manifold")

# -----------------------------
# 2. Early active latent trajectories
# -----------------------------
early <- behav %>%
  filter(str_detect(str_to_lower(as.character(Phase)), early_phase_pattern)) %>%
  group_by(AnimalNum, CageChange, Phase) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  mutate(EarlyBin = row_number(), EarlyHour = (EarlyBin - 1) / bins_per_hour) %>%
  filter(EarlyBin <= max_early_bins) %>%
  ungroup()
if (nrow(early) == 0) early <- behav %>% group_by(AnimalNum, CageChange, Phase) %>% arrange(TimeIndex, .by_group = TRUE) %>% mutate(EarlyBin = row_number(), EarlyHour = (EarlyBin - 1) / bins_per_hour) %>% filter(EarlyBin <= max_early_bins) %>% ungroup()

early_mat <- early %>% transmute(Movement = safe_scale(Movement), Entropy = safe_scale(Entropy), Proximity = safe_scale(Proximity)) %>% as.matrix()
early_mat[!is.finite(early_mat)] <- 0
early_pca <- prcomp(early_mat, center = FALSE, scale. = FALSE)
early_scores <- early %>% select(AnimalNum, Group, Sex, CageChange, Phase, EarlyBin, EarlyHour) %>% bind_cols(as_tibble(early_pca$x[, 1:2]) %>% rename(Latent1 = PC1, Latent2 = PC2))

traj <- early_scores %>%
  mutate(HourBin = floor(EarlyHour * 2) / 2) %>%
  group_by(Sex, Group, HourBin) %>%
  summarise(Latent1 = mean(Latent1, na.rm = TRUE), Latent2 = mean(Latent2, na.rm = TRUE), n_animals = n_distinct(AnimalNum), .groups = "drop") %>%
  arrange(Sex, Group, HourBin)

write_professional_table(early_scores, "early_active_latent_state_scores.csv")
write_professional_table(traj, "early_active_latent_state_group_trajectories.csv")

p_traj <- ggplot(traj, aes(Latent1, Latent2, colour = Group, group = Group)) +
  geom_path(linewidth = 0.85, alpha = 0.86, arrow = grid::arrow(length = grid::unit(2.2, "mm"), type = "closed")) +
  geom_point(data = traj %>% group_by(Sex, Group) %>% slice_min(HourBin, n = 1) %>% ungroup(), size = 1.9) +
  geom_point(data = traj %>% group_by(Sex, Group) %>% slice_max(HourBin, n = 1) %>% ungroup(), size = 2.8, shape = 8) +
  facet_grid(. ~ Sex) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  labs(title = "G  Early-active-phase trajectories through behavioral state space", subtitle = "Circle = start; star = end; arrows show group-average adaptation", x = "Latent state 1", y = "Latent state 2") +
  panel_theme(7)
save_manuscript_plot(p_traj, "Fig13G_early_active_latent_state_trajectories", width = 160, height = 70)

# -----------------------------
# 3. Unsupervised states and Markov transitions
# -----------------------------
state_input <- early %>% transmute(AnimalNum, Group, Sex, CageChange, Phase, EarlyBin, EarlyHour, Movement_z = safe_scale(Movement), Entropy_z = safe_scale(Entropy), Proximity_z = safe_scale(Proximity))
state_input <- state_input %>% mutate(across(c(Movement_z, Entropy_z, Proximity_z), ~ replace_na(.x, 0)))
set.seed(random_seed)
km <- kmeans(state_input %>% select(Movement_z, Entropy_z, Proximity_z), centers = n_behavioral_states, nstart = 50)
centres <- as_tibble(km$centers, rownames = "Cluster") %>% mutate(Cluster = as.integer(Cluster))

# Heuristic state labels: low movement = inactive; high proximity = social; high movement = burst; remaining = explore.
label_tbl <- centres %>% mutate(State = NA_character_)
label_tbl$State[which.min(label_tbl$Movement_z)] <- "Inactive"
label_tbl$State[which.max(label_tbl$Proximity_z)] <- "Social"
label_tbl$State[which.max(label_tbl$Movement_z)] <- "Burst"
label_tbl$State[is.na(label_tbl$State)] <- setdiff(state_labels, label_tbl$State[!is.na(label_tbl$State)])[seq_len(sum(is.na(label_tbl$State)))]

state_tbl <- state_input %>% mutate(Cluster = km$cluster) %>% left_join(label_tbl %>% select(Cluster, State), by = "Cluster") %>% mutate(State = factor(State, levels = state_labels))
state_occupancy <- state_tbl %>% count(AnimalNum, Group, Sex, State, name = "n") %>% group_by(AnimalNum, Group, Sex) %>% mutate(Occupancy = n / sum(n)) %>% ungroup()
transition_tbl <- state_tbl %>% arrange(AnimalNum, CageChange, Phase, EarlyBin) %>% group_by(AnimalNum, Group, Sex, CageChange, Phase) %>% mutate(From = State, To = lead(State)) %>% ungroup() %>% filter(!is.na(To)) %>% count(Group, Sex, From, To, name = "n") %>% group_by(Group, Sex, From) %>% mutate(Probability = n / sum(n)) %>% ungroup()

write_professional_table(label_tbl, "unsupervised_behavioral_state_cluster_labels.csv")
write_professional_table(state_tbl, "early_active_unsupervised_behavioral_states.csv")
write_professional_table(state_occupancy, "early_active_state_occupancy_by_animal.csv")
write_professional_table(transition_tbl, "early_active_state_transition_probabilities_by_sex.csv")

p_trans <- transition_tbl %>%
  mutate(From = factor(From, levels = state_labels), To = factor(To, levels = state_labels)) %>%
  ggplot(aes(To, From, fill = Probability)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  facet_grid(Sex ~ Group) +
  scale_fill_gradient(low = "white", high = "#333333") +
  labs(title = "H  Markov transition architecture of early behavior", subtitle = "Unsupervised states from movement, entropy and proximity", x = "To state", y = "From state", fill = "P") +
  panel_theme(6) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_manuscript_plot(p_trans, "Fig13H_early_active_markov_transition_heatmap", width = 165, height = 95)

p_occ <- state_occupancy %>%
  ggplot(aes(Group, Occupancy, fill = Group, colour = Group)) +
  geom_violin(alpha = 0.32, colour = NA, trim = FALSE) +
  geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.25, alpha = 0.72) +
  geom_jitter(width = 0.07, size = 0.8, alpha = 0.70, show.legend = FALSE) +
  facet_grid(Sex ~ State) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  labs(title = "I  Occupancy of unsupervised behavioral states", subtitle = "State dwell-time distribution in the early active window", x = NULL, y = "State occupancy") +
  panel_theme(6) + theme(legend.position = "none")
save_manuscript_plot(p_occ, "Fig13I_early_active_state_occupancy", width = 170, height = 95)

# -----------------------------
# 4. Recurrence maps
# -----------------------------
make_recurrence <- function(mat, q = 0.35) {
  d <- as.matrix(dist(mat))
  eps <- quantile(d[upper.tri(d)], q, na.rm = TRUE, names = FALSE)
  if (!is.finite(eps) || eps == 0) eps <- median(d[upper.tri(d)], na.rm = TRUE)
  exp(-d / eps)
}

recurrence_long <- traj %>%
  group_by(Group, HourBin) %>%
  summarise(Latent1 = mean(Latent1, na.rm = TRUE), Latent2 = mean(Latent2, na.rm = TRUE), .groups = "drop") %>%
  group_split(Group) %>%
  map_dfr(function(g) {
    g <- arrange(g, HourBin)
    rec <- make_recurrence(g %>% select(Latent1, Latent2) %>% as.matrix())
    hrs <- g$HourBin
    idx <- expand.grid(i = seq_along(hrs), j = seq_along(hrs))
    tibble(Group = unique(g$Group), Hour1 = hrs[idx$i], Hour2 = hrs[idx$j], Similarity = as.vector(rec))
  })
write_professional_table(recurrence_long, "early_active_group_recurrence_similarity_maps.csv")

p_rec <- recurrence_long %>%
  ggplot(aes(Hour1, Hour2, fill = Similarity)) +
  geom_raster() +
  facet_grid(. ~ Group) +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1)) +
  coord_equal(expand = FALSE) +
  labs(title = "J  Recurrence structure of early behavioral dynamics", subtitle = "Repeated motifs and state rigidity appear as time-time similarity structure", x = "Hour", y = "Hour", fill = "State\nsimilarity") +
  panel_theme(6)
save_manuscript_plot(p_rec, "Fig13J_early_active_recurrence_similarity_maps", width = 160, height = 62)

# -----------------------------
# 5. Nonlinear endpoint prediction and response curves
# -----------------------------
if (!is.na(outcome_to_use) && "outcome" %in% names(animal_features) && requireNamespace("ranger", quietly = TRUE)) {
  model_dat <- animal_features %>% filter(is.finite(outcome))
  if (nrow(model_dat) >= 8 && length(usable_features) >= 2) {
    model_x <- model_dat %>% select(all_of(usable_features)) %>% mutate(across(everything(), ~ replace_na(.x, median(.x, na.rm = TRUE))))
    rf_dat <- bind_cols(tibble(outcome = model_dat$outcome), model_x)
    set.seed(random_seed)
    rf <- ranger::ranger(outcome ~ ., data = rf_dat, num.trees = 1000, min.node.size = max(2, floor(nrow(rf_dat) / 10)), importance = "permutation", seed = random_seed)
    pred_tbl <- model_dat %>% transmute(AnimalNum, Group, Sex, outcome, predicted = rf$predictions, residual = outcome - predicted)
    imp <- tibble(feature = names(rf$variable.importance), importance = as.numeric(rf$variable.importance)) %>% arrange(desc(importance))
    perf <- tibble(model = "ranger_random_forest_in_sample", outcome = outcome_to_use, n = nrow(pred_tbl), pearson_r = safe_cor(pred_tbl$outcome, pred_tbl$predicted, "pearson"), spearman_rho = safe_cor(pred_tbl$outcome, pred_tbl$predicted, "spearman"), rmse = sqrt(mean((pred_tbl$outcome - pred_tbl$predicted)^2, na.rm = TRUE)), r2_in_sample = 1 - sum((pred_tbl$outcome - pred_tbl$predicted)^2, na.rm = TRUE) / sum((pred_tbl$outcome - mean(pred_tbl$outcome, na.rm = TRUE))^2, na.rm = TRUE))
    write_professional_table(pred_tbl, "nonlinear_endpoint_predictions.csv")
    write_professional_table(imp, "nonlinear_endpoint_feature_importance.csv", stats_dir)
    write_professional_table(perf, "nonlinear_endpoint_prediction_performance.csv", stats_dir)

    p_imp <- imp %>% slice_head(n = min(12, n())) %>% mutate(feature = factor(feature, levels = rev(feature))) %>%
      ggplot(aes(importance, feature)) + geom_col(fill = "#555555", width = 0.72) +
      labs(title = "K  Nonlinear feature importance", subtitle = paste0("Random forest endpoint model for ", outcome_to_use, "; screening, not causal inference"), x = "Permutation importance", y = NULL) + panel_theme(7) + theme(legend.position = "none")
    save_manuscript_plot(p_imp, "Fig13K_nonlinear_endpoint_feature_importance", width = 90, height = 78)

    p_pred <- pred_tbl %>% ggplot(aes(outcome, predicted, colour = Group, fill = Group, shape = Sex)) +
      geom_abline(slope = 1, intercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey45") +
      geom_point(size = 2.1, alpha = 0.92, stroke = 0.25) +
      geom_smooth(method = "lm", se = FALSE, linewidth = 0.55, colour = "grey20") +
      scale_colour_manual(values = group_colors, drop = FALSE) + scale_fill_manual(values = group_colors, drop = FALSE) +
      labs(title = "L  Nonlinear endpoint prediction", subtitle = paste0("In-sample nonlinear model: r=", round(perf$pearson_r, 2), ", R2=", round(perf$r2_in_sample, 2), "; validate by CV before reporting"), x = paste0("Observed ", outcome_to_use), y = paste0("Predicted ", outcome_to_use)) + panel_theme(7)
    save_manuscript_plot(p_pred, "Fig13L_nonlinear_endpoint_prediction", width = 90, height = 78)

    pdp_features <- imp$feature[seq_len(min(6, nrow(imp)))]
    pdp <- map_dfr(pdp_features, function(fc) {
      vals <- model_x[[fc]]
      grid <- seq(quantile(vals, 0.05, na.rm = TRUE), quantile(vals, 0.95, na.rm = TRUE), length.out = 60)
      yhat <- map_dbl(grid, function(v) { nd <- model_x; nd[[fc]] <- v; mean(predict(rf, data = nd)$predictions, na.rm = TRUE) })
      tibble(feature = fc, value = grid, predicted = yhat)
    })
    write_professional_table(pdp, "nonlinear_endpoint_partial_dependence_curves.csv")
    p_pdp <- pdp %>% ggplot(aes(value, predicted)) + geom_line(linewidth = 0.75, colour = "#333333") + facet_wrap(~ feature, scales = "free_x", ncol = 3) +
      labs(title = "M  Nonlinear feature-response curves", subtitle = "Partial dependence reveals threshold, saturation or U-shaped effects", x = "Feature value", y = paste0("Predicted ", outcome_to_use)) + panel_theme(6) + theme(legend.position = "none")
    save_manuscript_plot(p_pdp, "Fig13M_nonlinear_endpoint_partial_dependence_curves", width = 170, height = 90)
  }
}

# -----------------------------
# 6. Distance from nonlinear control geometry
# -----------------------------
if (all(c("Diffusion1", "Diffusion2") %in% names(embed))) {
  con_centres <- embed %>% filter(Group == "CON") %>% group_by(Sex) %>% summarise(CON1 = mean(Diffusion1, na.rm = TRUE), CON2 = mean(Diffusion2, na.rm = TRUE), .groups = "drop")
  dist_tbl <- embed %>% left_join(con_centres, by = "Sex") %>% mutate(DistanceToSexMatchedCON = sqrt((Diffusion1 - CON1)^2 + (Diffusion2 - CON2)^2))
  write_professional_table(dist_tbl, "distance_to_sex_matched_control_manifold.csv")
  p_dist <- dist_tbl %>% ggplot(aes(Group, DistanceToSexMatchedCON, fill = Group, colour = Group)) +
    geom_violin(alpha = 0.32, colour = NA, trim = FALSE) + geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.25, alpha = 0.72) + geom_jitter(width = 0.07, size = 1.0, alpha = 0.78, show.legend = FALSE) + facet_grid(. ~ Sex) +
    scale_fill_manual(values = group_colors, drop = FALSE) + scale_colour_manual(values = group_colors, drop = FALSE) +
    labs(title = "N  Nonlinear distance from sex-matched control geometry", subtitle = "Quantifies departure from CON behavioral manifold", x = NULL, y = "Distance to CON centre") + panel_theme(7) + theme(legend.position = "none")
  save_manuscript_plot(p_dist, "Fig13N_distance_to_sex_matched_control_manifold", width = 120, height = 70)
}

summary_tbl <- tibble(Item = c("Input", "Bin level", "Animals", "Usable features", "Early window", "States", "Outcome"), Value = c(input_file, bin_level, n_distinct(behav$AnimalNum), length(usable_features), paste0(early_window_hours, " h"), n_behavioral_states, ifelse(is.na(outcome_to_use), "none", outcome_to_use)))
write_professional_table(summary_tbl, "nonlinear_systems_analysis_summary.csv")

message("Nonlinear systems dynamics analysis complete: ", output_dir)
