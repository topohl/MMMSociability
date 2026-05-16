# ================================================================
# Duration-normalization helper functions
# MMMSociability
# ================================================================
# Purpose:
#   Shared safeguards for SIS cage-change analyses where some epochs have
#   shorter observation duration. Helpers export epoch-level duration QC,
#   detect short epochs/cage changes without hard-coding CC labels, and add
#   rate-normalized versions of count/cumulative metrics.
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(stringr)
})

if (!exists("%||%")) `%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

duration_safe_num <- function(x) suppressWarnings(as.numeric(x))

duration_first_existing_col <- function(dat, candidates, required = TRUE, label = "column") {
  hit <- candidates[candidates %in% names(dat)][1]
  if ((length(hit) == 0 || is.na(hit)) && required) {
    stop("Could not find ", label, ". Tried: ", paste(candidates, collapse = ", "), call. = FALSE)
  }
  if (length(hit) == 0) NA_character_ else hit
}

infer_bin_size_sec <- function(dat, default = NA_real_) {
  if ("BinSizeSec" %in% names(dat)) {
    x <- duration_safe_num(dat$BinSizeSec)
    x <- x[is.finite(x) & x > 0]
    if (length(x) > 0) return(stats::median(x, na.rm = TRUE))
  }
  if ("BinLevel" %in% names(dat)) {
    lvl <- unique(as.character(dat$BinLevel))[1]
    if (str_detect(lvl, "10sec")) return(10)
    if (str_detect(lvl, "1min")) return(60)
    if (str_detect(lvl, "5min")) return(300)
    if (str_detect(lvl, "10min")) return(600)
    if (str_detect(lvl, "30min")) return(1800)
  }
  if ("TimeIndex" %in% names(dat)) {
    ti <- sort(unique(duration_safe_num(dat$TimeIndex)))
    step <- stats::median(diff(ti), na.rm = TRUE)
    if (is.finite(step) && step > 0 && step < 100) return(step * 60)
  }
  default
}

calculate_observation_duration <- function(dat,
                                           metric_source = "behavior",
                                           bin_size_sec = NULL,
                                           short_epoch_fraction = 0.75,
                                           grouping_cols = c("AnimalNum", "Group", "Sex", "Batch", "System", "BinLevel", "CageChange", "Phase")) {
  if (is.null(dat) || nrow(dat) == 0) return(tibble())
  grouping_cols <- intersect(grouping_cols, names(dat))
  if (!"CageChange" %in% grouping_cols) return(tibble())
  if (!"Phase" %in% grouping_cols && "Phase" %in% names(dat)) grouping_cols <- c(grouping_cols, "Phase")

  bin_size_sec <- bin_size_sec %||% infer_bin_size_sec(dat)
  if (!is.finite(bin_size_sec) || bin_size_sec <= 0) bin_size_sec <- NA_real_
  obs_col <- duration_first_existing_col(dat, c("observation_seconds", "ObservationSeconds", "duration_seconds", "DurationSec"), FALSE, "observation seconds")
  dyad_col <- duration_first_existing_col(dat, c("dyadic_observation_seconds", "DyadicObservationSeconds"), FALSE, "dyadic observation seconds")
  phase_block_col <- duration_first_existing_col(dat, c("PhaseBlock", "PhaseNumber", "phase_block", "phase_number"), FALSE, "phase block")
  if (is.na(phase_block_col)) dat$.duration_phase_block <- "1" else dat$.duration_phase_block <- as.character(dat[[phase_block_col]])

  epoch_tbl <- dat %>%
    mutate(
      .obs_seconds = if (!is.na(obs_col)) duration_safe_num(.data[[obs_col]]) else NA_real_,
      .dyad_seconds = if (!is.na(dyad_col)) duration_safe_num(.data[[dyad_col]]) else NA_real_,
      .bin_seconds = if_else(is.finite(.obs_seconds) & .obs_seconds > 0, .obs_seconds, bin_size_sec),
      .bin_seconds = if_else(is.finite(.bin_seconds) & .bin_seconds > 0, .bin_seconds, NA_real_),
      Phase = if ("Phase" %in% names(.)) as.character(Phase) else "All",
      .phase_block = .duration_phase_block
    ) %>%
    group_by(across(all_of(grouping_cols))) %>%
    summarise(
      metric_source = metric_source,
      n_bins = n(),
      observed_bins = n(),
      n_phase_blocks = n_distinct(.phase_block),
      bin_size_sec = bin_size_sec,
      total_observation_duration_sec = sum(.bin_seconds, na.rm = TRUE),
      total_observation_duration_hours = total_observation_duration_sec / 3600,
      regrouping_exposure_dose_hours = total_observation_duration_hours,
      active_duration_hours = if_else(str_detect(str_to_lower(first(Phase)), "active|dark|night"),
                                      total_observation_duration_hours, 0),
      inactive_duration_hours = if_else(str_detect(str_to_lower(first(Phase)), "inactive|light|day"),
                                        total_observation_duration_hours, 0),
      dyadic_observation_duration_hours = sum(.dyad_seconds, na.rm = TRUE) / 3600,
      .groups = "drop"
    )

  expected_tbl <- epoch_tbl %>%
    group_by(across(any_of(c("metric_source", "BinLevel", "Batch", "System", "Sex", "Phase")))) %>%
    summarise(
      expected_bins = stats::median(observed_bins, na.rm = TRUE),
      expected_phase_blocks = stats::median(n_phase_blocks, na.rm = TRUE),
      expected_duration_hours = stats::median(total_observation_duration_hours, na.rm = TRUE),
      .groups = "drop"
    )

  cage_tbl <- epoch_tbl %>%
    group_by(metric_source, CageChange) %>%
    summarise(
      cage_change_total_bins = sum(observed_bins, na.rm = TRUE),
      cage_change_total_phase_blocks = sum(n_phase_blocks, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(metric_source) %>%
    mutate(
      median_cage_change_bins = stats::median(cage_change_total_bins, na.rm = TRUE),
      median_cage_change_phase_blocks = stats::median(cage_change_total_phase_blocks, na.rm = TRUE),
      cage_change_duration_fraction = cage_change_total_bins / median_cage_change_bins,
      cage_change_phase_block_fraction = cage_change_total_phase_blocks / median_cage_change_phase_blocks,
      cage_change_duration_class = if_else(cage_change_duration_fraction < short_epoch_fraction |
                                             cage_change_phase_block_fraction < short_epoch_fraction,
                                           "short", "standard")
    ) %>%
    ungroup()

  epoch_tbl %>%
    left_join(expected_tbl, by = intersect(names(epoch_tbl), names(expected_tbl))) %>%
    left_join(cage_tbl, by = c("metric_source", "CageChange")) %>%
    mutate(
      duration_completeness_fraction = observed_bins / expected_bins,
      phase_block_completeness_fraction = n_phase_blocks / expected_phase_blocks,
      short_epoch = duration_completeness_fraction < short_epoch_fraction |
        phase_block_completeness_fraction < short_epoch_fraction,
      short_regrouping_epoch = cage_change_duration_class == "short",
      reviewer_risk = case_when(
        short_epoch | cage_change_duration_class == "short" ~ "high",
        duration_completeness_fraction < 0.90 ~ "medium",
        TRUE ~ "low"
      )
    ) %>%
    arrange(CageChange, Phase, across(any_of("AnimalNum")))
}

write_epoch_duration_qc <- function(dat, output_dir, metric_source = "behavior", bin_size_sec = NULL) {
  qc <- calculate_observation_duration(dat, metric_source = metric_source, bin_size_sec = bin_size_sec)
  if (nrow(qc) > 0) {
    if (!dir.exists(file.path(output_dir, "tables"))) dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(qc, file.path(output_dir, "tables", "epoch_duration_qc.csv"))
  }
  qc
}

join_duration_qc <- function(dat, qc, by_cols = c("AnimalNum", "Group", "Sex", "CageChange", "Phase")) {
  if (is.null(qc) || nrow(qc) == 0) return(dat)
  by_cols <- intersect(by_cols, intersect(names(dat), names(qc)))
  duration_cols <- c(
    "observed_bins", "expected_bins", "total_observation_duration_hours",
    "regrouping_exposure_dose_hours", "active_duration_hours", "inactive_duration_hours",
    "n_phase_blocks", "expected_phase_blocks", "phase_block_completeness_fraction",
    "duration_completeness_fraction", "cage_change_duration_fraction",
    "cage_change_phase_block_fraction", "cage_change_duration_class",
    "short_epoch", "short_regrouping_epoch"
  )
  dat <- dat %>% select(-any_of(duration_cols))
  qc_join <- qc %>%
    select(any_of(by_cols), any_of(duration_cols)) %>%
    group_by(across(all_of(by_cols))) %>%
    summarise(
      observed_bins = sum(observed_bins, na.rm = TRUE),
      expected_bins = sum(expected_bins, na.rm = TRUE),
      total_observation_duration_hours = sum(total_observation_duration_hours, na.rm = TRUE),
      regrouping_exposure_dose_hours = sum(regrouping_exposure_dose_hours, na.rm = TRUE),
      active_duration_hours = sum(active_duration_hours, na.rm = TRUE),
      inactive_duration_hours = sum(inactive_duration_hours, na.rm = TRUE),
      n_phase_blocks = sum(n_phase_blocks, na.rm = TRUE),
      expected_phase_blocks = sum(expected_phase_blocks, na.rm = TRUE),
      phase_block_completeness_fraction = n_phase_blocks / expected_phase_blocks,
      duration_completeness_fraction = observed_bins / expected_bins,
      cage_change_duration_fraction = min(cage_change_duration_fraction, na.rm = TRUE),
      cage_change_phase_block_fraction = min(cage_change_phase_block_fraction, na.rm = TRUE),
      cage_change_duration_class = if_else(any(cage_change_duration_class == "short", na.rm = TRUE), "short", "standard"),
      short_epoch = any(short_epoch %in% TRUE, na.rm = TRUE),
      short_regrouping_epoch = any(short_regrouping_epoch %in% TRUE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      across(c(cage_change_duration_fraction, cage_change_phase_block_fraction), ~ if_else(is.infinite(.x), NA_real_, .x))
    )
  dat %>%
    left_join(qc_join, by = by_cols)
}

filter_short_duration_epochs <- function(dat, qc, by_cols = c("AnimalNum", "Group", "Sex", "CageChange", "Phase")) {
  joined <- join_duration_qc(dat, qc, by_cols = by_cols)
  joined %>% filter(!short_epoch %in% TRUE & (is.na(cage_change_duration_class) | cage_change_duration_class != "short"))
}

normalize_counts_to_rates <- function(dat, count_cols, duration_col = "total_observation_duration_hours", suffix = "_per_hour") {
  count_cols <- intersect(count_cols, names(dat))
  if (length(count_cols) == 0 || !duration_col %in% names(dat)) return(dat)
  for (cc in count_cols) {
    dat[[paste0(cc, suffix)]] <- ifelse(is.finite(dat[[duration_col]]) & dat[[duration_col]] > 0,
                                        duration_safe_num(dat[[cc]]) / dat[[duration_col]], NA_real_)
  }
  dat
}

normalize_auc_per_hour <- function(dat, auc_cols = "auc", duration_col = "total_observation_duration_hours") {
  auc_cols <- intersect(auc_cols, names(dat))
  if (length(auc_cols) == 0 || !duration_col %in% names(dat)) return(dat)
  for (cc in auc_cols) {
    dat[[paste0(cc, "_per_hour")]] <- ifelse(is.finite(dat[[duration_col]]) & dat[[duration_col]] > 0,
                                             duration_safe_num(dat[[cc]]) / dat[[duration_col]], NA_real_)
  }
  dat
}

filter_minimum_bins <- function(x, min_bins = 4) {
  if (length(x[is.finite(x)]) < min_bins) return(rep(NA_real_, length(x)))
  x
}

generate_duration_sensitivity_outputs <- function(full_tbl,
                                                  qc,
                                                  output_dir,
                                                  file_stem,
                                                  by_cols = c("AnimalNum", "Group", "Sex", "CageChange", "Phase")) {
  if (!dir.exists(file.path(output_dir, "tables"))) dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  full_tagged <- join_duration_qc(full_tbl, qc, by_cols = by_cols) %>% mutate(DurationAnalysisSet = "full")
  no_short <- full_tagged %>%
    filter(!short_epoch %in% TRUE, is.na(cage_change_duration_class) | cage_change_duration_class != "short") %>%
    mutate(DurationAnalysisSet = "excluding_short_duration")
  readr::write_csv(full_tagged, file.path(output_dir, "tables", paste0(file_stem, "_duration_tagged_full.csv")))
  readr::write_csv(no_short, file.path(output_dir, "tables", paste0(file_stem, "_excluding_short_duration.csv")))
  bind_rows(full_tagged, no_short)
}

make_duration_analysis_sets <- function(dat,
                                        qc,
                                        by_cols = c("AnimalNum", "Group", "Sex", "CageChange", "Phase")) {
  full_tagged <- join_duration_qc(dat, qc, by_cols = by_cols) %>%
    mutate(DurationAnalysisSet = "full")
  no_short <- full_tagged %>%
    filter(!short_epoch %in% TRUE, is.na(cage_change_duration_class) | cage_change_duration_class != "short") %>%
    mutate(DurationAnalysisSet = "excluding_short_duration")
  bind_rows(full_tagged, no_short)
}

write_duration_sensitivity_statistics <- function(dat,
                                                  qc,
                                                  output_dir,
                                                  analysis_name,
                                                  value_cols,
                                                  by_cols,
                                                  summary_group_cols,
                                                  join_by_cols = c("AnimalNum", "Group", "Sex", "CageChange", "Phase")) {
  if (!exists("make_dynamics_group_summary") || !exists("make_dynamics_group_contrasts")) {
    return(invisible(NULL))
  }
  ensure_dir <- if (exists("ensure_dir")) get("ensure_dir") else function(path) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  ensure_dir(file.path(output_dir, "stats_tables"))
  duration_sets <- make_duration_analysis_sets(dat, qc, by_cols = join_by_cols)
  value_cols <- intersect(value_cols, names(duration_sets))
  if (length(value_cols) == 0) return(invisible(NULL))

  summary_tbl <- make_dynamics_group_summary(
    duration_sets,
    value_cols = value_cols,
    group_cols = unique(c("DurationAnalysisSet", summary_group_cols))
  )
  contrast_tbl <- make_dynamics_group_contrasts(
    duration_sets,
    value_cols = value_cols,
    by_cols = unique(c("DurationAnalysisSet", by_cols))
  )
  if ("p.adjust_bh_family" %in% names(contrast_tbl)) {
    contrast_tbl <- contrast_tbl %>%
      mutate(
        ReportingP = p.adjust_bh_family,
        ReportingCorrection = "BH FDR within duration-analysis set and requested comparison family"
      )
  }

  readr::write_csv(summary_tbl, file.path(output_dir, "stats_tables", paste0(analysis_name, "_duration_sensitivity_group_summary.csv")))
  readr::write_csv(contrast_tbl, file.path(output_dir, "stats_tables", paste0(analysis_name, "_duration_sensitivity_group_contrasts.csv")))

  if (all(c("DurationAnalysisSet", "Outcome", "contrast", "cohen_d") %in% names(contrast_tbl))) {
    consistency_tbl <- contrast_tbl %>%
      select(any_of(c(by_cols, "Outcome", "contrast", "DurationAnalysisSet", "estimate", "cohen_d", "p.value", "ReportingP"))) %>%
      pivot_wider(
        names_from = DurationAnalysisSet,
        values_from = c(estimate, cohen_d, p.value, ReportingP),
        names_sep = "__"
      )
    for (needed_col in c("cohen_d__full", "cohen_d__excluding_short_duration", "estimate__full", "estimate__excluding_short_duration")) {
      if (!needed_col %in% names(consistency_tbl)) consistency_tbl[[needed_col]] <- NA_real_
    }
    consistency_tbl <- consistency_tbl %>%
      mutate(
        effect_direction_consistent = sign(cohen_d__full) == sign(cohen_d__excluding_short_duration),
        delta_cohen_d_excluding_short = cohen_d__excluding_short_duration - cohen_d__full,
        sensitivity_flag = case_when(
          is.na(effect_direction_consistent) ~ "not_comparable",
          !effect_direction_consistent ~ "direction_changed",
          abs(delta_cohen_d_excluding_short) >= 0.30 ~ "effect_size_shift",
          TRUE ~ "stable"
        )
      )
    readr::write_csv(consistency_tbl, file.path(output_dir, "stats_tables", paste0(analysis_name, "_duration_sensitivity_consistency.csv")))
  }

  invisible(list(summary = summary_tbl, contrasts = contrast_tbl))
}

make_duration_sensitivity_audit <- function(script, rows) {
  tibble(script = script) %>%
    bind_cols(as_tibble(rows)) %>%
    mutate(
      duration_sensitive = as.logical(duration_sensitive),
      cc4_exclusion_sensitivity_available = as.logical(cc4_exclusion_sensitivity_available),
      reviewer_risk = factor(reviewer_risk, levels = c("low", "medium", "high"))
    )
}
