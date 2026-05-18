# ================================================================
# RFID CHIP LOSS / TRACKING INTEGRITY QC
# MMMSociability
# ================================================================
# Goal:
#   Detect animals whose RFID-based trajectories are likely invalid because
#   the chip was lost, stationary in the bedding, or otherwise no longer
#   reflects the animal's true position.
#
# Important principle:
#   This script is diagnostic and non-destructive. It does NOT remove animals
#   from the data. It produces QC tables, plots, and suggested flags that can
#   be reviewed before updating excluded_animals.csv or a future
#   chip_loss_qc.csv / valid_until table.
#
# Recommended use:
#   1. Run after preprocessing / multiscale metric generation.
#   2. Inspect QC tables and plots.
#   3. Decide whether animals should be globally excluded or censored after a
#      conservative valid_until time.
#   4. Keep final exclusion decisions in a small manually reviewed CSV.
#
# Why this matters:
#   Unknown RFID chip loss can create pseudo-data: the chip may still be read,
#   but movement, entropy, and proximity no longer represent the mouse.
#   This is different from ordinary missingness.
# ================================================================

# -----------------------------
# 0. USER SETTINGS
# -----------------------------

# Leave NULL to auto-detect likely metric files. Otherwise provide paths.
INPUT_FILES <- NULL

# Optional: restrict auto-detection to these directories if they exist.
SEARCH_DIRS <- c(
  "analysis_ready",
  "analysis_ready/03_derived_metrics",
  "analysis_ready/03_derived_metrics/phase_based",
  "analysis_ready/03_derived_metrics/halfhour_based",
  "Results",
  "Formatting",
  "."
)

OUT_DIR <- file.path("analysis_ready", "00_tracking_qc_rfid_loss")

# Metrics used for RFID-loss suspicion.
# Adjust after inspecting the output distributions.
QC_THRESHOLDS <- list(
  zero_movement_fraction_high = 0.95,
  dominant_position_fraction_high = 0.95,
  mean_entropy_low = 0.05,
  mean_positions_visited_low = 1.05,
  transition_rate_low = 0.02,
  longest_zero_run_bins_high = 200,
  rolling_zero_run_bins_high = 100
)

# Rolling window length in rows within each animal x batch x cage change x phase.
# For 10-sec data, 360 rows = 60 min. For 5-min data, 12 rows = 60 min.
# The script also reports the inferred BinSizeSec where available.
ROLLING_WINDOW_ROWS <- 120

# If TRUE, writes a suggested list of animals for manual review.
WRITE_SUGGESTED_EXCLUSIONS <- TRUE

# -----------------------------
# 1. PACKAGES
# -----------------------------

required_packages <- c(
  "dplyr", "tidyr", "readr", "stringr", "purrr", "ggplot2",
  "openxlsx", "lubridate", "scales", "forcats"
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

invisible(lapply(required_packages, install_if_missing))
invisible(lapply(required_packages, library, character.only = TRUE))

# -----------------------------
# 2. HELPERS
# -----------------------------

message_header <- function(x) {
  message("\n", paste(rep("=", nchar(x) + 8), collapse = ""))
  message("=== ", x, " ===")
  message(paste(rep("=", nchar(x) + 8), collapse = ""))
}

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

standardize_names <- function(df) {
  nm <- names(df)

  rename_if_present <- function(from, to) {
    idx <- which(names(df) %in% from)
    if (length(idx) > 0 && !(to %in% names(df))) names(df)[idx[1]] <<- to
  }

  rename_if_present(c("Animal", "animal", "AnimalID", "Animal_Id", "animal_id"), "AnimalID")
  rename_if_present(c("AnimalNum", "animal_num", "Animal_Number"), "AnimalNum")
  rename_if_present(c("DateTime", "Datetime", "datetime", "Time", "Timestamp", "BinStart"), "DateTime")
  rename_if_present(c("CageChange", "CC", "Cage_Change"), "CageChange")
  rename_if_present(c("Phase", "phase"), "Phase")
  rename_if_present(c("Batch", "batch"), "Batch")
  rename_if_present(c("System", "system"), "System")
  rename_if_present(c("Movement", "movement"), "Movement")
  rename_if_present(c("MovementDistance", "movement_distance"), "MovementDistance")
  rename_if_present(c("Entropy", "entropy"), "Entropy")
  rename_if_present(c("DominantPosition", "dominant_position", "Position", "position"), "DominantPosition")
  rename_if_present(c("n_positions_visited", "NPositionsVisited", "positions_visited"), "n_positions_visited")
  rename_if_present(c("BinSizeSec", "bin_size_sec"), "BinSizeSec")
  rename_if_present(c("TimeIndex", "time_index"), "TimeIndex")

  df
}

read_one_file <- function(path) {
  message("Reading: ", path)
  ext <- tolower(tools::file_ext(path))

  df <- switch(
    ext,
    "csv" = readr::read_csv(path, show_col_types = FALSE),
    "tsv" = readr::read_tsv(path, show_col_types = FALSE),
    "txt" = readr::read_delim(path, delim = "\t", show_col_types = FALSE),
    "rds" = readRDS(path),
    stop("Unsupported file type: ", path)
  )

  df <- standardize_names(as.data.frame(df))
  df$SourceFile <- basename(path)
  df$SourcePath <- path
  df
}

find_candidate_files <- function(search_dirs = SEARCH_DIRS) {
  existing_dirs <- search_dirs[dir.exists(search_dirs)]
  if (length(existing_dirs) == 0) stop("None of SEARCH_DIRS exist. Set INPUT_FILES manually.")

  files <- unlist(lapply(existing_dirs, function(d) {
    list.files(d, pattern = "\\.(csv|tsv|txt|rds)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  }))

  files <- unique(files)

  # Prefer derived metric files and avoid previous QC outputs.
  keep <- stringr::str_detect(
    basename(files),
    regex("metric|multiscale|halfhour|phase|AnimalPos|preprocessed|movement|entropy|proximity", ignore_case = TRUE)
  ) & !stringr::str_detect(files, regex("tracking_qc|rfid|suggested_excluded|excluded_animals", ignore_case = TRUE))

  files[keep]
}

parse_datetime_safe <- function(x) {
  if (inherits(x, "POSIXct")) return(x)
  if (inherits(x, "Date")) return(as.POSIXct(x))
  suppressWarnings(lubridate::ymd_hms(x, tz = "UTC", quiet = TRUE)) %>%
    { ifelse(is.na(.), suppressWarnings(lubridate::ymd_hm(x, tz = "UTC", quiet = TRUE)), .) } %>%
    as.POSIXct(origin = "1970-01-01", tz = "UTC")
}

longest_true_run <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0 || !any(x)) return(0L)
  r <- rle(x)
  max(r$lengths[r$values], na.rm = TRUE)
}

rolling_true_count <- function(x, window) {
  x <- as.integer(replace_na(as.logical(x), FALSE))
  n <- length(x)
  if (n == 0) return(integer(0))
  out <- rep(NA_integer_, n)
  cs <- c(0L, cumsum(x))
  for (i in seq_len(n)) {
    lo <- max(1L, i - window + 1L)
    out[i] <- cs[i + 1L] - cs[lo]
  }
  out
}

first_sustained_flag_time <- function(df, flag_col = "rolling_zero_flag") {
  idx <- which(df[[flag_col]] %in% TRUE)
  if (length(idx) == 0) return(as.POSIXct(NA))
  df$DateTime[idx[1]]
}

save_svg <- function(filename, plot, width = 12, height = 7) {
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = "svg"
  )
}

# -----------------------------
# 3. LOAD DATA
# -----------------------------

message_header("RFID tracking QC")
safe_dir_create(OUT_DIR)
safe_dir_create(file.path(OUT_DIR, "figures"))
safe_dir_create(file.path(OUT_DIR, "tables"))

if (is.null(INPUT_FILES)) {
  INPUT_FILES <- find_candidate_files()
}

if (length(INPUT_FILES) == 0) {
  stop("No input files found. Set INPUT_FILES manually to derived metric/preprocessed files.")
}

message("Candidate input files: ", length(INPUT_FILES))

raw_list <- purrr::map(INPUT_FILES, read_one_file)
data <- dplyr::bind_rows(raw_list)
data <- standardize_names(data)

# Required identifiers.
if (!("AnimalID" %in% names(data)) && "AnimalNum" %in% names(data)) {
  data$AnimalID <- as.character(data$AnimalNum)
}

required_cols <- c("AnimalID", "Movement")
missing_required <- setdiff(required_cols, names(data))
if (length(missing_required) > 0) {
  stop("Missing required columns: ", paste(missing_required, collapse = ", "))
}

# Create robust grouping columns where absent.
if (!("Batch" %in% names(data))) data$Batch <- "UnknownBatch"
if (!("CageChange" %in% names(data))) data$CageChange <- "UnknownCC"
if (!("Phase" %in% names(data))) data$Phase <- "UnknownPhase"
if (!("System" %in% names(data))) data$System <- "UnknownSystem"
if (!("DateTime" %in% names(data))) {
  if ("TimeIndex" %in% names(data)) {
    data$DateTime <- as.POSIXct(as.numeric(data$TimeIndex), origin = "1970-01-01", tz = "UTC")
  } else {
    data$DateTime <- as.POSIXct(NA)
  }
} else {
  data$DateTime <- parse_datetime_safe(data$DateTime)
}

# Ensure numeric where expected.
num_cols <- intersect(c("Movement", "MovementDistance", "Entropy", "n_positions_visited", "BinSizeSec", "TimeIndex"), names(data))
data <- data %>% mutate(across(all_of(num_cols), ~ suppressWarnings(as.numeric(.x))))

data <- data %>%
  mutate(
    AnimalID = as.character(AnimalID),
    Batch = as.character(Batch),
    CageChange = as.character(CageChange),
    Phase = as.character(Phase),
    System = as.character(System),
    zero_movement = !is.na(Movement) & Movement == 0,
    valid_movement = !is.na(Movement)
  )

# -----------------------------
# 4. ANIMAL x WINDOW QC METRICS
# -----------------------------

message_header("Calculating QC metrics")

qc_group_cols <- c("AnimalID", "Batch", "CageChange", "Phase", "System")
optional_group_cols <- intersect(c("Sex", "Group"), names(data))
qc_group_cols <- c(qc_group_cols, optional_group_cols)

# Dominant position fraction only if a position-like variable exists.
has_position <- "DominantPosition" %in% names(data)
has_entropy <- "Entropy" %in% names(data)
has_positions_visited <- "n_positions_visited" %in% names(data)
has_binsize <- "BinSizeSec" %in% names(data)

qc_by_animal_phase <- data %>%
  arrange(across(all_of(qc_group_cols)), DateTime, .by_group = FALSE) %>%
  group_by(across(all_of(qc_group_cols))) %>%
  summarise(
    n_rows = n(),
    n_valid_movement = sum(valid_movement, na.rm = TRUE),
    first_time = suppressWarnings(min(DateTime, na.rm = TRUE)),
    last_time = suppressWarnings(max(DateTime, na.rm = TRUE)),
    inferred_bin_size_sec = if (has_binsize) median(BinSizeSec, na.rm = TRUE) else NA_real_,
    mean_movement = mean(Movement, na.rm = TRUE),
    median_movement = median(Movement, na.rm = TRUE),
    sd_movement = sd(Movement, na.rm = TRUE),
    zero_movement_fraction = mean(zero_movement, na.rm = TRUE),
    longest_zero_run_bins = longest_true_run(zero_movement),
    mean_entropy = if (has_entropy) mean(Entropy, na.rm = TRUE) else NA_real_,
    median_entropy = if (has_entropy) median(Entropy, na.rm = TRUE) else NA_real_,
    mean_positions_visited = if (has_positions_visited) mean(n_positions_visited, na.rm = TRUE) else NA_real_,
    dominant_position_fraction = if (has_position) {
      pos <- DominantPosition[!is.na(DominantPosition)]
      if (length(pos) == 0) NA_real_ else max(table(pos)) / length(pos)
    } else NA_real_,
    n_unique_positions = if (has_position) n_distinct(DominantPosition, na.rm = TRUE) else NA_integer_,
    transition_rate = if (has_position) {
      pos <- DominantPosition[!is.na(DominantPosition)]
      if (length(pos) < 2) NA_real_ else mean(pos[-1] != pos[-length(pos)], na.rm = TRUE)
    } else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    first_time = ifelse(is.infinite(first_time), NA, first_time) %>% as.POSIXct(origin = "1970-01-01", tz = "UTC"),
    last_time = ifelse(is.infinite(last_time), NA, last_time) %>% as.POSIXct(origin = "1970-01-01", tz = "UTC"),
    flag_high_zero_movement = zero_movement_fraction >= QC_THRESHOLDS$zero_movement_fraction_high,
    flag_stationary_position = !is.na(dominant_position_fraction) & dominant_position_fraction >= QC_THRESHOLDS$dominant_position_fraction_high,
    flag_low_entropy = !is.na(mean_entropy) & mean_entropy <= QC_THRESHOLDS$mean_entropy_low,
    flag_low_positions = !is.na(mean_positions_visited) & mean_positions_visited <= QC_THRESHOLDS$mean_positions_visited_low,
    flag_low_transition_rate = !is.na(transition_rate) & transition_rate <= QC_THRESHOLDS$transition_rate_low,
    flag_long_zero_run = longest_zero_run_bins >= QC_THRESHOLDS$longest_zero_run_bins_high,
    n_qc_flags = rowSums(across(starts_with("flag_")), na.rm = TRUE),
    tracking_qc_class = case_when(
      n_qc_flags >= 3 ~ "high_suspicion",
      n_qc_flags == 2 ~ "moderate_suspicion",
      n_qc_flags == 1 ~ "review",
      TRUE ~ "pass"
    )
  )

# -----------------------------
# 5. ROLLING COLLAPSE DETECTION
# -----------------------------

message_header("Calculating rolling collapse flags")

rolling_qc <- data %>%
  arrange(across(all_of(qc_group_cols)), DateTime, .by_group = FALSE) %>%
  group_by(across(all_of(qc_group_cols))) %>%
  mutate(
    rolling_zero_count = rolling_true_count(zero_movement, ROLLING_WINDOW_ROWS),
    rolling_zero_fraction = rolling_zero_count / pmin(row_number(), ROLLING_WINDOW_ROWS),
    rolling_zero_flag = rolling_zero_count >= QC_THRESHOLDS$rolling_zero_run_bins_high
  ) %>%
  ungroup()

rolling_onsets <- rolling_qc %>%
  group_by(across(all_of(qc_group_cols))) %>%
  summarise(
    first_rolling_zero_flag_time = first_sustained_flag_time(cur_data_all(), "rolling_zero_flag"),
    max_rolling_zero_fraction = max(rolling_zero_fraction, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    first_rolling_zero_flag_time = as.POSIXct(first_rolling_zero_flag_time, origin = "1970-01-01", tz = "UTC"),
    max_rolling_zero_fraction = ifelse(is.infinite(max_rolling_zero_fraction), NA_real_, max_rolling_zero_fraction)
  )

qc_by_animal_phase <- qc_by_animal_phase %>%
  left_join(rolling_onsets, by = qc_group_cols) %>%
  mutate(
    flag_rolling_zero_collapse = !is.na(first_rolling_zero_flag_time),
    n_qc_flags = rowSums(across(starts_with("flag_")), na.rm = TRUE),
    tracking_qc_class = case_when(
      n_qc_flags >= 3 ~ "high_suspicion",
      n_qc_flags == 2 ~ "moderate_suspicion",
      n_qc_flags == 1 ~ "review",
      TRUE ~ "pass"
    )
  )

qc_by_animal <- qc_by_animal_phase %>%
  group_by(AnimalID, across(any_of(optional_group_cols))) %>%
  summarise(
    n_windows = n(),
    n_pass = sum(tracking_qc_class == "pass", na.rm = TRUE),
    n_review = sum(tracking_qc_class == "review", na.rm = TRUE),
    n_moderate = sum(tracking_qc_class == "moderate_suspicion", na.rm = TRUE),
    n_high = sum(tracking_qc_class == "high_suspicion", na.rm = TRUE),
    max_flags = max(n_qc_flags, na.rm = TRUE),
    worst_zero_fraction = max(zero_movement_fraction, na.rm = TRUE),
    worst_dominant_position_fraction = max(dominant_position_fraction, na.rm = TRUE),
    lowest_mean_entropy = suppressWarnings(min(mean_entropy, na.rm = TRUE)),
    earliest_possible_collapse = suppressWarnings(min(first_rolling_zero_flag_time, na.rm = TRUE)),
    suggested_decision = case_when(
      n_high > 0 ~ "manual_review_high_suspicion_possible_chip_loss",
      n_moderate > 0 ~ "manual_review_moderate_suspicion",
      n_review > 0 ~ "review_single_qc_flag",
      TRUE ~ "pass"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    lowest_mean_entropy = ifelse(is.infinite(lowest_mean_entropy), NA_real_, lowest_mean_entropy),
    earliest_possible_collapse = as.POSIXct(
      ifelse(is.infinite(earliest_possible_collapse), NA, earliest_possible_collapse),
      origin = "1970-01-01",
      tz = "UTC"
    )
  )

# -----------------------------
# 6. SAVE TABLES
# -----------------------------

message_header("Writing QC tables")

tables_dir <- file.path(OUT_DIR, "tables")
figures_dir <- file.path(OUT_DIR, "figures")

readr::write_csv(qc_by_animal_phase, file.path(tables_dir, "tracking_qc_by_animal_phase.csv"))
readr::write_csv(qc_by_animal, file.path(tables_dir, "tracking_qc_by_animal.csv"))

if (WRITE_SUGGESTED_EXCLUSIONS) {
  suggested <- qc_by_animal %>%
    filter(suggested_decision != "pass") %>%
    transmute(
      AnimalID,
      across(any_of(optional_group_cols)),
      suggested_decision,
      earliest_possible_collapse,
      comment = "Manual review required. Do not automatically exclude solely from this table."
    )
  readr::write_csv(suggested, file.path(tables_dir, "suggested_animals_for_manual_tracking_review.csv"))
}

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "animal_phase_qc")
openxlsx::writeData(wb, "animal_phase_qc", qc_by_animal_phase)
openxlsx::addWorksheet(wb, "animal_summary")
openxlsx::writeData(wb, "animal_summary", qc_by_animal)
if (WRITE_SUGGESTED_EXCLUSIONS) {
  openxlsx::addWorksheet(wb, "suggested_review")
  openxlsx::writeData(wb, "suggested_review", suggested)
}
openxlsx::addWorksheet(wb, "thresholds")
openxlsx::writeData(wb, "thresholds", data.frame(
  threshold = names(QC_THRESHOLDS),
  value = unlist(QC_THRESHOLDS),
  row.names = NULL
))
openxlsx::saveWorkbook(wb, file.path(tables_dir, "tracking_qc_rfid_loss_report.xlsx"), overwrite = TRUE)

# -----------------------------
# 7. PLOTS
# -----------------------------

message_header("Generating QC plots")

# Animal-level zero movement heatmap by window.
plot_data_phase <- qc_by_animal_phase %>%
  mutate(
    Window = paste(Batch, CageChange, Phase, sep = "_"),
    AnimalID = forcats::fct_reorder(AnimalID, worst <- zero_movement_fraction, .fun = max, .desc = TRUE)
  )

p_zero_heatmap <- ggplot(plot_data_phase, aes(x = Window, y = AnimalID, fill = zero_movement_fraction)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1), labels = scales::percent_format()) +
  labs(
    title = "RFID tracking QC: zero-movement fraction",
    subtitle = "High sustained zero movement can indicate detached or stationary RFID chips",
    x = "Batch / cage change / phase",
    y = "Animal",
    fill = "Zero movement"
  ) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    legend.position = "right"
  )

save_svg(file.path(figures_dir, "qc_zero_movement_fraction_heatmap.svg"), p_zero_heatmap, width = 14, height = 8)

# QC flag count heatmap.
p_flags_heatmap <- ggplot(plot_data_phase, aes(x = Window, y = AnimalID, fill = n_qc_flags)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(option = "plasma", breaks = 0:6) +
  labs(
    title = "RFID tracking QC: number of suspicious tracking flags",
    x = "Batch / cage change / phase",
    y = "Animal",
    fill = "QC flags"
  ) +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))

save_svg(file.path(figures_dir, "qc_flag_count_heatmap.svg"), p_flags_heatmap, width = 14, height = 8)

# Scatter of zero movement vs entropy, if entropy is available.
if (has_entropy) {
  p_entropy <- qc_by_animal_phase %>%
    ggplot(aes(x = zero_movement_fraction, y = mean_entropy, label = AnimalID)) +
    geom_point(aes(shape = tracking_qc_class), size = 2, alpha = 0.8) +
    geom_vline(xintercept = QC_THRESHOLDS$zero_movement_fraction_high, linetype = "dashed", linewidth = 0.3) +
    geom_hline(yintercept = QC_THRESHOLDS$mean_entropy_low, linetype = "dashed", linewidth = 0.3) +
    labs(
      title = "RFID tracking QC: movement collapse versus entropy collapse",
      x = "Zero-movement fraction",
      y = "Mean entropy",
      shape = "QC class"
    ) +
    theme_classic(base_size = 9)

  save_svg(file.path(figures_dir, "qc_zero_movement_vs_entropy.svg"), p_entropy, width = 8, height = 6)
}

# Dominant position fraction, if available.
if (has_position) {
  p_dompos <- qc_by_animal_phase %>%
    ggplot(aes(x = zero_movement_fraction, y = dominant_position_fraction, shape = tracking_qc_class)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_vline(xintercept = QC_THRESHOLDS$zero_movement_fraction_high, linetype = "dashed", linewidth = 0.3) +
    geom_hline(yintercept = QC_THRESHOLDS$dominant_position_fraction_high, linetype = "dashed", linewidth = 0.3) +
    labs(
      title = "RFID tracking QC: stationary movement versus dominant position",
      x = "Zero-movement fraction",
      y = "Dominant position fraction",
      shape = "QC class"
    ) +
    theme_classic(base_size = 9)

  save_svg(file.path(figures_dir, "qc_zero_movement_vs_dominant_position.svg"), p_dompos, width = 8, height = 6)
}

# Rolling zero-movement time series for flagged animals only.
flagged_animals <- qc_by_animal %>%
  filter(suggested_decision != "pass") %>%
  pull(AnimalID) %>%
  unique()

if (length(flagged_animals) > 0 && any(!is.na(rolling_qc$DateTime))) {
  p_roll <- rolling_qc %>%
    filter(AnimalID %in% flagged_animals) %>%
    mutate(Window = paste(Batch, CageChange, Phase, sep = "_")) %>%
    ggplot(aes(x = DateTime, y = rolling_zero_fraction)) +
    geom_line(linewidth = 0.3, alpha = 0.8) +
    geom_hline(yintercept = QC_THRESHOLDS$rolling_zero_run_bins_high / ROLLING_WINDOW_ROWS,
               linetype = "dashed", linewidth = 0.3) +
    facet_grid(AnimalID ~ Window, scales = "free_x") +
    labs(
      title = "RFID tracking QC: rolling zero-movement fraction in flagged animals",
      x = "Time",
      y = "Rolling zero-movement fraction"
    ) +
    theme_classic(base_size = 7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_svg(file.path(figures_dir, "qc_flagged_animals_rolling_zero_movement.svg"), p_roll, width = 16, height = max(5, 1.2 * length(flagged_animals)))
}

# Summary barplot.
p_summary <- qc_by_animal %>%
  count(suggested_decision) %>%
  ggplot(aes(x = suggested_decision, y = n)) +
  geom_col(width = 0.7) +
  labs(
    title = "RFID tracking QC summary",
    x = "Suggested decision",
    y = "Number of animals"
  ) +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_svg(file.path(figures_dir, "qc_suggested_decision_summary.svg"), p_summary, width = 8, height = 5)

# -----------------------------
# 8. MANUSCRIPT / METHODS NOTE
# -----------------------------

methods_note <- c(
  "RFID tracking integrity QC",
  "",
  "Animals were screened for potential RFID-chip loss or detached-chip artifacts using a non-destructive QC layer applied before downstream behavioral analyses.",
  "QC metrics included the fraction of zero-movement bins, longest zero-movement run, dominant-position occupancy, transition rate, mean entropy, and rolling zero-movement collapse within animal-by-cage-change-by-phase windows.",
  "Animals or windows exceeding multiple QC criteria were flagged for manual review rather than automatically excluded.",
  "Where the timing of RFID loss could not be reliably established, whole-animal exclusion is conservative because detached chips can remain detectable and produce non-biological pseudo-movement.",
  "Where a defensible last valid interval can be established, animals can instead be retained up to that point and censored thereafter."
)
writeLines(methods_note, file.path(OUT_DIR, "rfid_tracking_qc_methods_note.txt"))

message_header("Done")
message("Outputs written to: ", normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE))
message("Review the Excel report and SVG heatmaps before changing excluded_animals.csv.")
