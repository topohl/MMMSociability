# ================================================================
# RAW RFID TRACKING QC FOR SUSPECTED CHIP LOSS
# MMMSociability / E9 SIS AnimalPos files
# ================================================================
# Goal:
#   Detect likely RFID chip-loss / detached-chip artifacts directly from the
#   raw AnimalPos files used by E9_SIS_AnimalPos-preprocessing.R.
#
# Why raw-level QC:
#   If a chip is lost but remains in the bedding, it can still generate reads.
#   These reads are not missing data; they are pseudo-data. Therefore QC should
#   be performed before movement, entropy, proximity, network, HMM, GAMM, and
#   prediction features are derived.
#
# What this script does:
#   - reads raw_data/B*/E9_SIS_B*_CC*_AnimalPos.csv
#   - parses Animal into AnimalID and System
#   - maps xPos/yPos into PositionID using the same 8-position cage map as the
#     preprocessing script
#   - quantifies stationary / low-transition / position-fixation patterns
#   - estimates conservative candidate valid-until timestamps
#   - writes review tables and SVG plots
#
# What this script does NOT do:
#   - it does not edit excluded_animals.csv
#   - it does not remove data
#   - it does not automatically decide biological exclusion
#
# Recommended workflow:
#   1. Run this script before E9_SIS_AnimalPos-preprocessing.R.
#   2. Review tables/plots in raw_tracking_qc_rfid_loss/.
#   3. Decide manually:
#        a) global exclusion -> raw_data/excluded_animals.csv
#        b) censoring after valid_until -> future chip_loss_qc.csv handling
#        c) retain animal -> no action
# ================================================================

# -----------------------------
# 0. SETTINGS
# -----------------------------

WORKING_DIR <- getwd()
RAW_DIR <- file.path(WORKING_DIR, "raw_data")
OUT_DIR <- file.path(WORKING_DIR, "raw_tracking_qc_rfid_loss")

BATCHES <- c("B1", "B2", "B3", "B4", "B5", "B6")
CAGE_CHANGES <- c("CC1", "CC2", "CC3", "CC4")

# Same position map as E9_SIS_AnimalPos-preprocessing.R
POSITION_MAP <- tibble::tibble(
  PositionID = 1:8,
  xPos = c(0, 100, 200, 300, 0, 100, 200, 300),
  yPos = c(0,   0,   0,   0, 116, 116, 116, 116)
)

# Main QC thresholds. Treat as review thresholds, not automatic exclusion rules.
QC_THRESHOLDS <- list(
  dominant_position_fraction_high = 0.95,
  transition_rate_low = 0.01,
  longest_stationary_run_high = 500,
  rolling_stationary_fraction_high = 0.98,
  rolling_transition_rate_low = 0.005,
  min_reads_for_window_qc = 100,
  max_inter_read_gap_sec_high = 3600
)

# Rolling window in raw rows per animal. This is row-based because raw read rates
# may not be perfectly regular. Increase for stricter detection.
ROLLING_WINDOW_ROWS <- 250

# Optional: exclude known globally excluded animals from plots? For diagnosing
# why they were excluded, FALSE is usually better.
DROP_EXISTING_EXCLUDED_ANIMALS <- FALSE

# -----------------------------
# 1. PACKAGES
# -----------------------------

required_packages <- c(
  "readr", "dplyr", "tidyr", "stringr", "lubridate", "purrr",
  "ggplot2", "openxlsx", "scales", "forcats", "tibble"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)
}

# -----------------------------
# 2. HELPERS
# -----------------------------

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

longest_true_run <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0 || !any(x)) return(0L)
  r <- rle(x)
  max(r$lengths[r$values], na.rm = TRUE)
}

rolling_sum <- function(x, window) {
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

rolling_mean_numeric <- function(x, window) {
  x <- as.numeric(x)
  n <- length(x)
  if (n == 0) return(numeric(0))
  out <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    lo <- max(1L, i - window + 1L)
    out[i] <- mean(x[lo:i], na.rm = TRUE)
  }
  out
}

parse_datetime_raw <- function(x) {
  if (inherits(x, "POSIXct")) return(x)
  out <- suppressWarnings(as.POSIXct(x, format = "%d.%m.%Y %H:%M:%S", tz = "UTC"))
  if (all(is.na(out))) {
    out <- suppressWarnings(lubridate::ymd_hms(x, tz = "UTC", quiet = TRUE))
  }
  out
}

map_position_id <- function(df) {
  df %>%
    mutate(
      xPos = suppressWarnings(as.numeric(xPos)),
      yPos = suppressWarnings(as.numeric(yPos))
    ) %>%
    left_join(POSITION_MAP, by = c("xPos", "yPos"))
}

read_raw_file <- function(batch, cage_change) {
  filename <- paste0("E9_SIS_", batch, "_", cage_change, "_AnimalPos.csv")
  path <- file.path(RAW_DIR, batch, filename)

  if (!file.exists(path)) {
    warning("Missing file: ", path)
    return(NULL)
  }

  message("Reading raw AnimalPos file: ", path)

  dat <- readr::read_delim(path, delim = ";", show_col_types = FALSE, progress = FALSE)

  required <- c("DateTime", "Animal", "xPos", "yPos")
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0) {
    warning("Skipping ", path, "; missing columns: ", paste(missing, collapse = ", "))
    return(NULL)
  }

  dat %>%
    mutate(
      Batch = batch,
      CageChange = cage_change,
      SourceFile = filename,
      DateTime = parse_datetime_raw(DateTime),
      Animal = as.character(Animal)
    ) %>%
    tidyr::separate(Animal, into = c("AnimalID", "System"), sep = "[_-]", remove = FALSE, fill = "right", extra = "merge") %>%
    mutate(
      AnimalID = as.character(AnimalID),
      System = as.character(System),
      Phase = if_else(
        format(DateTime, "%H:%M", tz = "UTC") >= "18:30" |
          format(DateTime, "%H:%M", tz = "UTC") < "06:30",
        "Active", "Inactive"
      )
    ) %>%
    map_position_id() %>%
    arrange(AnimalID, System, DateTime)
}

first_candidate_valid_until <- function(df) {
  idx <- which(df$candidate_chip_loss_onset %in% TRUE)
  if (length(idx) == 0) return(as.POSIXct(NA))

  # conservative choice: last time before the sustained suspicious window starts
  i <- max(1L, idx[1] - 1L)
  df$DateTime[i]
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
# 3. LOAD RAW DATA
# -----------------------------

safe_dir_create(OUT_DIR)
safe_dir_create(file.path(OUT_DIR, "tables"))
safe_dir_create(file.path(OUT_DIR, "figures"))

raw_data <- purrr::map_dfr(BATCHES, function(b) {
  purrr::map_dfr(CAGE_CHANGES, function(cc) read_raw_file(b, cc))
})

if (nrow(raw_data) == 0) {
  stop("No raw AnimalPos rows loaded. Check RAW_DIR, BATCHES, and CAGE_CHANGES.")
}

excluded_path <- file.path(RAW_DIR, "excluded_animals.csv")
existing_excluded <- character(0)
if (file.exists(excluded_path)) {
  existing_excluded <- readLines(excluded_path, warn = FALSE)
  existing_excluded <- existing_excluded[nzchar(existing_excluded)]
}

raw_data <- raw_data %>%
  mutate(
    is_existing_excluded = AnimalID %in% existing_excluded,
    valid_position = !is.na(PositionID)
  )

if (DROP_EXISTING_EXCLUDED_ANIMALS && length(existing_excluded) > 0) {
  raw_data <- raw_data %>% filter(!is_existing_excluded)
}

# -----------------------------
# 4. RAW SEQUENCE FEATURES
# -----------------------------

raw_seq <- raw_data %>%
  arrange(Batch, CageChange, Phase, System, AnimalID, DateTime) %>%
  group_by(Batch, CageChange, Phase, System, AnimalID) %>%
  mutate(
    prev_position = lag(PositionID),
    position_changed = !is.na(PositionID) & !is.na(prev_position) & PositionID != prev_position,
    stationary_read = !is.na(PositionID) & !is.na(prev_position) & PositionID == prev_position,
    inter_read_gap_sec = as.numeric(difftime(DateTime, lag(DateTime), units = "secs")),
    rolling_stationary_count = rolling_sum(stationary_read, ROLLING_WINDOW_ROWS),
    rolling_stationary_fraction = rolling_stationary_count / pmax(1, pmin(row_number(), ROLLING_WINDOW_ROWS)),
    rolling_transition_rate = rolling_mean_numeric(position_changed, ROLLING_WINDOW_ROWS),
    candidate_chip_loss_onset =
      row_number() >= QC_THRESHOLDS$min_reads_for_window_qc &
      rolling_stationary_fraction >= QC_THRESHOLDS$rolling_stationary_fraction_high &
      rolling_transition_rate <= QC_THRESHOLDS$rolling_transition_rate_low
  ) %>%
  ungroup()

# -----------------------------
# 5. SUMMARY TABLES
# -----------------------------

qc_by_animal_phase <- raw_seq %>%
  group_by(Batch, CageChange, Phase, System, AnimalID, is_existing_excluded) %>%
  summarise(
    n_raw_reads = n(),
    n_valid_position_reads = sum(valid_position, na.rm = TRUE),
    first_time = min(DateTime, na.rm = TRUE),
    last_time = max(DateTime, na.rm = TRUE),
    median_inter_read_gap_sec = median(inter_read_gap_sec, na.rm = TRUE),
    max_inter_read_gap_sec = max(inter_read_gap_sec, na.rm = TRUE),
    n_unique_positions = n_distinct(PositionID, na.rm = TRUE),
    dominant_position_fraction = {
      pos <- PositionID[!is.na(PositionID)]
      if (length(pos) == 0) NA_real_ else max(table(pos)) / length(pos)
    },
    transition_rate = mean(position_changed, na.rm = TRUE),
    stationary_fraction = mean(stationary_read, na.rm = TRUE),
    longest_stationary_run = longest_true_run(stationary_read),
    max_rolling_stationary_fraction = max(rolling_stationary_fraction, na.rm = TRUE),
    min_rolling_transition_rate = min(rolling_transition_rate, na.rm = TRUE),
    candidate_valid_until = first_candidate_valid_until(cur_data_all()),
    .groups = "drop"
  ) %>%
  mutate(
    max_inter_read_gap_sec = if_else(is.infinite(max_inter_read_gap_sec), NA_real_, max_inter_read_gap_sec),
    min_rolling_transition_rate = if_else(is.infinite(min_rolling_transition_rate), NA_real_, min_rolling_transition_rate),
    candidate_valid_until = as.POSIXct(candidate_valid_until, origin = "1970-01-01", tz = "UTC"),
    flag_position_fixation = !is.na(dominant_position_fraction) &
      dominant_position_fraction >= QC_THRESHOLDS$dominant_position_fraction_high,
    flag_low_transition_rate = !is.na(transition_rate) &
      transition_rate <= QC_THRESHOLDS$transition_rate_low,
    flag_long_stationary_run = !is.na(longest_stationary_run) &
      longest_stationary_run >= QC_THRESHOLDS$longest_stationary_run_high,
    flag_rolling_stationary_collapse = !is.na(candidate_valid_until),
    flag_large_inter_read_gap = !is.na(max_inter_read_gap_sec) &
      max_inter_read_gap_sec >= QC_THRESHOLDS$max_inter_read_gap_sec_high,
    n_raw_qc_flags = rowSums(across(starts_with("flag_")), na.rm = TRUE),
    raw_tracking_qc_class = case_when(
      n_raw_qc_flags >= 3 ~ "high_suspicion",
      n_raw_qc_flags == 2 ~ "moderate_suspicion",
      n_raw_qc_flags == 1 ~ "review",
      TRUE ~ "pass"
    )
  )

qc_by_animal <- qc_by_animal_phase %>%
  group_by(AnimalID, is_existing_excluded) %>%
  summarise(
    n_windows = n(),
    n_pass = sum(raw_tracking_qc_class == "pass", na.rm = TRUE),
    n_review = sum(raw_tracking_qc_class == "review", na.rm = TRUE),
    n_moderate = sum(raw_tracking_qc_class == "moderate_suspicion", na.rm = TRUE),
    n_high = sum(raw_tracking_qc_class == "high_suspicion", na.rm = TRUE),
    max_raw_qc_flags = max(n_raw_qc_flags, na.rm = TRUE),
    worst_dominant_position_fraction = max(dominant_position_fraction, na.rm = TRUE),
    lowest_transition_rate = min(transition_rate, na.rm = TRUE),
    longest_stationary_run_any_window = max(longest_stationary_run, na.rm = TRUE),
    earliest_candidate_valid_until = min(candidate_valid_until, na.rm = TRUE),
    suggested_decision = case_when(
      n_high > 0 ~ "manual_review_high_suspicion_possible_chip_loss",
      n_moderate > 0 ~ "manual_review_moderate_suspicion",
      n_review > 0 ~ "review_single_raw_qc_flag",
      TRUE ~ "pass"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    lowest_transition_rate = if_else(is.infinite(lowest_transition_rate), NA_real_, lowest_transition_rate),
    earliest_candidate_valid_until = as.POSIXct(
      ifelse(is.infinite(earliest_candidate_valid_until), NA, earliest_candidate_valid_until),
      origin = "1970-01-01",
      tz = "UTC"
    )
  )

suggested_valid_until <- qc_by_animal_phase %>%
  filter(!is.na(candidate_valid_until)) %>%
  transmute(
    AnimalID,
    Batch,
    CageChange,
    Phase,
    System,
    valid_until = candidate_valid_until,
    raw_tracking_qc_class,
    n_raw_qc_flags,
    reason = "candidate sustained stationary/low-transition raw RFID pattern; manual review required"
  ) %>%
  arrange(AnimalID, valid_until)

suggested_review <- qc_by_animal %>%
  filter(suggested_decision != "pass") %>%
  arrange(desc(max_raw_qc_flags), AnimalID)

# -----------------------------
# 6. WRITE TABLES
# -----------------------------

tables_dir <- file.path(OUT_DIR, "tables")
figures_dir <- file.path(OUT_DIR, "figures")

readr::write_csv(qc_by_animal_phase, file.path(tables_dir, "raw_tracking_qc_by_animal_phase.csv"))
readr::write_csv(qc_by_animal, file.path(tables_dir, "raw_tracking_qc_by_animal.csv"))
readr::write_csv(suggested_review, file.path(tables_dir, "suggested_animals_for_manual_raw_tracking_review.csv"))
readr::write_csv(suggested_valid_until, file.path(tables_dir, "suggested_valid_until_candidates.csv"))

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "animal_phase_raw_qc")
openxlsx::writeData(wb, "animal_phase_raw_qc", qc_by_animal_phase)
openxlsx::addWorksheet(wb, "animal_raw_summary")
openxlsx::writeData(wb, "animal_raw_summary", qc_by_animal)
openxlsx::addWorksheet(wb, "suggested_review")
openxlsx::writeData(wb, "suggested_review", suggested_review)
openxlsx::addWorksheet(wb, "suggested_valid_until")
openxlsx::writeData(wb, "suggested_valid_until", suggested_valid_until)
openxlsx::addWorksheet(wb, "thresholds")
openxlsx::writeData(wb, "thresholds", data.frame(
  threshold = names(QC_THRESHOLDS),
  value = unlist(QC_THRESHOLDS),
  row.names = NULL
))
openxlsx::saveWorkbook(wb, file.path(tables_dir, "raw_tracking_qc_rfid_loss_report.xlsx"), overwrite = TRUE)

# -----------------------------
# 7. PLOTS
# -----------------------------

plot_phase <- qc_by_animal_phase %>%
  mutate(
    Window = paste(Batch, CageChange, Phase, System, sep = "_"),
    AnimalID = forcats::fct_reorder(AnimalID, n_raw_qc_flags, .fun = max, .desc = TRUE)
  )

p_dom <- ggplot(plot_phase, aes(x = Window, y = AnimalID, fill = dominant_position_fraction)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1), labels = scales::percent_format()) +
  labs(
    title = "Raw RFID QC: dominant-position fraction",
    subtitle = "Detached chips often show persistent occupancy at one antenna/position",
    x = "Batch / cage change / phase / system",
    y = "Animal",
    fill = "Dominant\nposition"
  ) +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))

save_svg(file.path(figures_dir, "raw_qc_dominant_position_fraction_heatmap.svg"), p_dom, width = 16, height = 8)

p_trans <- ggplot(plot_phase, aes(x = Window, y = AnimalID, fill = transition_rate)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(option = "viridis") +
  labs(
    title = "Raw RFID QC: antenna/position transition rate",
    subtitle = "Very low transition rates suggest stationary-chip artifacts or tracking collapse",
    x = "Batch / cage change / phase / system",
    y = "Animal",
    fill = "Transition\nrate"
  ) +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))

save_svg(file.path(figures_dir, "raw_qc_transition_rate_heatmap.svg"), p_trans, width = 16, height = 8)

p_flags <- ggplot(plot_phase, aes(x = Window, y = AnimalID, fill = n_raw_qc_flags)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(option = "plasma", breaks = 0:5) +
  labs(
    title = "Raw RFID QC: number of tracking-failure flags",
    x = "Batch / cage change / phase / system",
    y = "Animal",
    fill = "Raw QC\nflags"
  ) +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))

save_svg(file.path(figures_dir, "raw_qc_flag_count_heatmap.svg"), p_flags, width = 16, height = 8)

p_scatter <- qc_by_animal_phase %>%
  ggplot(aes(x = dominant_position_fraction, y = transition_rate, shape = raw_tracking_qc_class)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_vline(xintercept = QC_THRESHOLDS$dominant_position_fraction_high, linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = QC_THRESHOLDS$transition_rate_low, linetype = "dashed", linewidth = 0.3) +
  labs(
    title = "Raw RFID QC: position fixation versus transition collapse",
    x = "Dominant-position fraction",
    y = "Transition rate",
    shape = "QC class"
  ) +
  theme_classic(base_size = 9)

save_svg(file.path(figures_dir, "raw_qc_position_fixation_vs_transition_rate.svg"), p_scatter, width = 8, height = 6)

flagged_animals <- suggested_review$AnimalID %>% unique()

if (length(flagged_animals) > 0) {
  p_roll <- raw_seq %>%
    filter(AnimalID %in% flagged_animals) %>%
    mutate(Window = paste(Batch, CageChange, Phase, System, sep = "_")) %>%
    ggplot(aes(x = DateTime, y = rolling_stationary_fraction)) +
    geom_line(linewidth = 0.25, alpha = 0.8) +
    geom_hline(yintercept = QC_THRESHOLDS$rolling_stationary_fraction_high, linetype = "dashed", linewidth = 0.3) +
    facet_grid(AnimalID ~ Window, scales = "free_x") +
    labs(
      title = "Raw RFID QC: rolling stationary-read fraction in flagged animals",
      x = "Time",
      y = "Rolling stationary fraction"
    ) +
    theme_classic(base_size = 7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_svg(
    file.path(figures_dir, "raw_qc_flagged_animals_rolling_stationary_fraction.svg"),
    p_roll,
    width = 16,
    height = max(5, 1.2 * length(flagged_animals))
  )
}

p_summary <- qc_by_animal %>%
  count(suggested_decision, is_existing_excluded) %>%
  ggplot(aes(x = suggested_decision, y = n, fill = is_existing_excluded)) +
  geom_col(position = "stack", width = 0.7) +
  labs(
    title = "Raw RFID QC summary",
    subtitle = "Existing excluded_animals.csv status is shown for comparison",
    x = "Suggested decision",
    y = "Number of animals",
    fill = "Already in\nexcluded_animals.csv"
  ) +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_svg(file.path(figures_dir, "raw_qc_suggested_decision_summary.svg"), p_summary, width = 9, height = 5)

# -----------------------------
# 8. METHODS NOTE
# -----------------------------

methods_note <- c(
  "Raw RFID tracking integrity QC",
  "",
  "Raw AnimalPos files were screened before behavioral feature extraction for patterns consistent with RFID-chip loss or detached-chip artifacts.",
  "For each animal within batch, cage-change, system, and light/dark phase, we quantified dominant-position occupancy, antenna/position transition rate, stationary-read fraction, longest stationary run, inter-read gaps, and rolling stationary collapse.",
  "Candidate chip-loss windows were flagged when sustained raw-position stationarity coincided with near-zero transition rate across a rolling read window.",
  "Flagged animals or intervals were not automatically excluded; they were exported for manual review together with conservative candidate valid-until timestamps.",
  "This approach distinguishes true missing data from pseudo-data generated by detached RFID chips that remain detectable in the cage bedding."
)
writeLines(methods_note, file.path(OUT_DIR, "raw_rfid_tracking_qc_methods_note.txt"))

message("\nRaw RFID QC complete.")
message("Output directory: ", normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE))
message("Review suggested_animals_for_manual_raw_tracking_review.csv and suggested_valid_until_candidates.csv before changing preprocessing exclusions.")
