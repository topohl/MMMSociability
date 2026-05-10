# ================================================================
# Structure check for behavioral dynamics analyses
# MMMSociability
# ================================================================
# Purpose:
#   Preflight-check an input table before running:
#     Analysis/06_burstiness_temporal_instability.R
#     Analysis/07_behavioral_state_space.R
#
# Checks:
#   - file can be read
#   - required columns can be detected
#   - movement / entropy / proximity are numeric or coercible
#   - one animal ID and one time column exist
#   - duplicate animal-time rows are flagged
#   - enough bins per animal exist for RMSSD/autocorrelation/state analyses
#   - phase / cage-change / group / sex columns are reported if available
#
# Output:
#   analysis_ready/00_structure_checks/behavioral_dynamics_structure_report.csv
#   analysis_ready/00_structure_checks/behavioral_dynamics_column_mapping.csv
#   analysis_ready/00_structure_checks/behavioral_dynamics_bins_per_animal.csv
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
})

source("Functions/behavioral_dynamics_helpers.R")

# ------------------------------------------------
# USER INPUT
# ------------------------------------------------

input_file <- "analysis_ready/03_derived_metrics/phase_based/all_behavior_metrics.csv"
output_dir <- "analysis_ready/00_structure_checks"

# Optional manual overrides. Leave NULL for auto-detection.
animal_col <- NULL
time_col <- NULL
group_col <- NULL
sex_col <- NULL
phase_col <- NULL
cage_col <- NULL
movement_col <- NULL
entropy_col <- NULL
proximity_col <- NULL

# ------------------------------------------------
# HELPERS
# ------------------------------------------------

status_row <- function(check, status, detail = NA_character_, severity = "info") {
  tibble(
    check = check,
    status = status,
    severity = severity,
    detail = as.character(detail)
  )
}

soft_find_col <- function(dat, candidates) {
  hit <- candidates[candidates %in% names(dat)][1]
  ifelse(is.na(hit), NA_character_, hit)
}

numeric_quality <- function(x) {
  raw_non_na <- sum(!is.na(x))
  coerced <- suppressWarnings(as.numeric(x))
  numeric_non_na <- sum(is.finite(coerced))
  tibble(
    raw_non_na = raw_non_na,
    numeric_non_na = numeric_non_na,
    coercion_success_fraction = ifelse(raw_non_na > 0, numeric_non_na / raw_non_na, NA_real_),
    n_negative = sum(coerced < 0, na.rm = TRUE),
    n_zero = sum(coerced == 0, na.rm = TRUE),
    min = suppressWarnings(min(coerced, na.rm = TRUE)),
    median = suppressWarnings(median(coerced, na.rm = TRUE)),
    max = suppressWarnings(max(coerced, na.rm = TRUE))
  )
}

# ------------------------------------------------
# LOAD
# ------------------------------------------------

ensure_dir(output_dir)

report <- list()

if (!file.exists(input_file)) {
  report <- list(status_row(
    "input_file_exists",
    "FAIL",
    paste0("File not found: ", input_file),
    "error"
  ))
  final_report <- bind_rows(report)
  write_table(final_report, file.path(output_dir, "behavioral_dynamics_structure_report.csv"))
  stop("Input file not found. Report written to: ", output_dir, call. = FALSE)
}

raw_dat <- read_behavior_table(input_file)

report <- append(report, list(status_row(
  "input_file_readable",
  "PASS",
  paste0("Rows: ", nrow(raw_dat), "; columns: ", ncol(raw_dat)),
  "info"
)))

# ------------------------------------------------
# COLUMN MAPPING
# ------------------------------------------------

mapping <- tibble(
  role = c("animal", "time", "group", "sex", "phase", "cage_change", "movement", "entropy", "proximity"),
  detected_column = c(
    animal_col %||% soft_find_col(raw_dat, c("AnimalNum", "Animal", "MouseID", "Mouse", "ID", "RFID", "animal_id")),
    time_col %||% soft_find_col(raw_dat, c("HalfHourElapsed", "HalfHourWithinCC0", "HalfHour", "Time", "TimeBin", "ZeitgeberTime", "ZT", "datetime", "DateTime")),
    group_col %||% soft_find_col(raw_dat, c("Group", "Phenotype", "Condition", "Treatment", "StressGroup")),
    sex_col %||% soft_find_col(raw_dat, c("Sex", "sex")),
    phase_col %||% soft_find_col(raw_dat, c("Phase", "phase", "LightDark", "DayNight", "CircadianPhase")),
    cage_col %||% soft_find_col(raw_dat, c("CageChange", "CC", "CageChangeNum", "Regrouping", "Batch", "Cage")),
    movement_col %||% soft_find_col(raw_dat, c("Movement", "movement", "Distance", "distance", "Activity", "activity")),
    entropy_col %||% soft_find_col(raw_dat, c("Entropy", "entropy", "ShannonEntropy", "shannon_entropy", "PositionEntropy")),
    proximity_col %||% soft_find_col(raw_dat, c("Proximity", "proximity", "MeanProximity", "SocialProximity", "CloseProximity"))
  ),
  required_for_dynamics = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE)
)

write_table(mapping, file.path(output_dir, "behavioral_dynamics_column_mapping.csv"))

missing_required <- mapping %>%
  filter(required_for_dynamics, is.na(detected_column))

if (nrow(missing_required) > 0) {
  report <- append(report, list(status_row(
    "required_columns_detected",
    "FAIL",
    paste0("Missing required roles: ", paste(missing_required$role, collapse = ", ")),
    "error"
  )))
} else {
  report <- append(report, list(status_row(
    "required_columns_detected",
    "PASS",
    paste0("Detected: ", paste(mapping$role[!is.na(mapping$detected_column)], mapping$detected_column[!is.na(mapping$detected_column)], sep = "=", collapse = "; ")),
    "info"
  )))
}

optional_missing <- mapping %>%
  filter(!required_for_dynamics, is.na(detected_column))

if (nrow(optional_missing) > 0) {
  report <- append(report, list(status_row(
    "optional_columns_detected",
    "WARN",
    paste0("Optional roles not found: ", paste(optional_missing$role, collapse = ", "), ". The analyses will use 'All' for missing stratifiers."),
    "warning"
  )))
} else {
  report <- append(report, list(status_row(
    "optional_columns_detected",
    "PASS",
    "Group/sex/phase/cage-change stratifiers detected.",
    "info"
  )))
}

# Stop only after writing useful mapping.
if (nrow(missing_required) > 0) {
  final_report <- bind_rows(report)
  write_table(final_report, file.path(output_dir, "behavioral_dynamics_structure_report.csv"))
  stop("Required columns missing. See structure report and column mapping.", call. = FALSE)
}

# ------------------------------------------------
# STANDARDIZE + CHECK NUMERIC QUALITY
# ------------------------------------------------

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

metric_quality <- bind_rows(
  numeric_quality(behav$Movement) %>% mutate(metric = "Movement"),
  numeric_quality(behav$Entropy) %>% mutate(metric = "Entropy"),
  numeric_quality(behav$Proximity) %>% mutate(metric = "Proximity")
) %>%
  relocate(metric)

write_table(metric_quality, file.path(output_dir, "behavioral_dynamics_numeric_quality.csv"))

bad_numeric <- metric_quality %>%
  filter(is.na(coercion_success_fraction) | coercion_success_fraction < 0.95)

if (nrow(bad_numeric) > 0) {
  report <- append(report, list(status_row(
    "metric_numeric_quality",
    "WARN",
    paste0("Potential numeric coercion problem for: ", paste(bad_numeric$metric, collapse = ", ")),
    "warning"
  )))
} else {
  report <- append(report, list(status_row(
    "metric_numeric_quality",
    "PASS",
    "Movement, Entropy, and Proximity are numeric/coercible for >=95% of non-missing values.",
    "info"
  )))
}

# ------------------------------------------------
# DUPLICATE ANIMAL-TIME ROWS
# ------------------------------------------------

dup_tbl <- behav %>%
  count(AnimalNum, CageChange, Phase, TimeIndex, name = "n") %>%
  filter(n > 1)

write_table(dup_tbl, file.path(output_dir, "behavioral_dynamics_duplicate_animal_time_rows.csv"))

if (nrow(dup_tbl) > 0) {
  report <- append(report, list(status_row(
    "duplicate_animal_time_rows",
    "WARN",
    paste0(nrow(dup_tbl), " duplicated AnimalNum x CageChange x Phase x TimeIndex combinations detected."),
    "warning"
  )))
} else {
  report <- append(report, list(status_row(
    "duplicate_animal_time_rows",
    "PASS",
    "No duplicate AnimalNum x CageChange x Phase x TimeIndex combinations detected.",
    "info"
  )))
}

# ------------------------------------------------
# BINS PER ANIMAL / PERIOD
# ------------------------------------------------

bins_per_animal <- behav %>%
  group_by(Group, Sex, Phase, CageChange, AnimalNum) %>%
  summarise(
    n_bins = n(),
    n_movement = sum(is.finite(Movement)),
    n_entropy = sum(is.finite(Entropy)),
    n_proximity = sum(is.finite(Proximity)),
    time_min = suppressWarnings(min(TimeIndex, na.rm = TRUE)),
    time_max = suppressWarnings(max(TimeIndex, na.rm = TRUE)),
    .groups = "drop"
  )

write_table(bins_per_animal, file.path(output_dir, "behavioral_dynamics_bins_per_animal.csv"))

low_bins <- bins_per_animal %>%
  filter(n_bins < 4 | n_movement < 4 | n_entropy < 4 | n_proximity < 4)

if (nrow(low_bins) > 0) {
  report <- append(report, list(status_row(
    "minimum_bins_for_temporal_metrics",
    "WARN",
    paste0(nrow(low_bins), " animal-periods have <4 usable bins. RMSSD/ACF/state switching may be unstable there."),
    "warning"
  )))
} else {
  report <- append(report, list(status_row(
    "minimum_bins_for_temporal_metrics",
    "PASS",
    "All animal-periods have >=4 bins for temporal metrics.",
    "info"
  )))
}

# ------------------------------------------------
# STRATIFIER LEVELS
# ------------------------------------------------

level_summary <- behav %>%
  summarise(
    n_animals = n_distinct(AnimalNum),
    n_groups = n_distinct(Group),
    groups = paste(sort(unique(as.character(Group))), collapse = "; "),
    n_sexes = n_distinct(Sex),
    sexes = paste(sort(unique(as.character(Sex))), collapse = "; "),
    n_phases = n_distinct(Phase),
    phases = paste(sort(unique(as.character(Phase))), collapse = "; "),
    n_cage_changes = n_distinct(CageChange),
    cage_changes = paste(sort(unique(as.character(CageChange))), collapse = "; ")
  )

write_table(level_summary, file.path(output_dir, "behavioral_dynamics_level_summary.csv"))

report <- append(report, list(status_row(
  "level_summary",
  "INFO",
  paste0(
    "Animals=", level_summary$n_animals,
    "; Groups=", level_summary$groups,
    "; Sexes=", level_summary$sexes,
    "; Phases=", level_summary$phases,
    "; Cage changes=", level_summary$cage_changes
  ),
  "info"
)))

# ------------------------------------------------
# FINAL REPORT
# ------------------------------------------------

final_report <- bind_rows(report)
write_table(final_report, file.path(output_dir, "behavioral_dynamics_structure_report.csv"))

print(final_report)

if (any(final_report$status == "FAIL")) {
  stop("Structure check failed. See report in: ", output_dir, call. = FALSE)
}

message("Structure check complete. Reports written to: ", output_dir)
