# ================================================================
# Build dyadic RFID contact table from preprocessed position data
# MMMSociability
# ================================================================
# Goal:
#   Convert long-term RFID position events into true dyadic social-contact
#   summaries that can feed dynamic social-network analyses.
#
# Input:
#   preprocessed_data/*_preprocessed.csv
#   Expected core columns:
#     DateTime, AnimalID, System, PositionID
#   Optional columns used when present:
#     Phase, CageChange, Batch, HalfHoursElapsed, Group, Sex
#
# Output:
#   analysis_ready/06_behavioral_dynamics/dyadic_contacts/tables/
#     dyadic_contacts_by_bin.csv
#     dyadic_contacts_interval_level.csv
#     dyadic_contact_qc_by_file.csv
#
# Interpretation:
#   same_position_fraction = fraction of observed seconds where a dyad is
#                            detected in the same RFID position.
#   adjacent_fraction      = fraction of observed seconds where a dyad is in
#                            neighbouring RFID positions on the 2 x 4 grid.
#   Weight                 = default edge weight for graph analysis; currently
#                            same_position_fraction.
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(lubridate)
  library(parallel)
})

source("C:/Users/topohl/Documents/GitHub/MMMSociability/Functions/behavioral_dynamics_helpers.R")
source("C:/Users/topohl/Documents/GitHub/MMMSociability/Functions/duration_normalization_helpers.R")

# ------------------------------------------------
# USER INPUT
# ------------------------------------------------

input_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability/preprocessed_data"
output_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/analysis_ready/06_behavioral_dynamics/dyadic_contacts"

# The aggregation bin determines the temporal resolution of the dyadic table.
# 1800 sec = 30 min, matching the current main analysis scale.
# Use 600, 300, 60, or 10 for 10 min, 5 min, 1 min, or 10 sec analyses.
bin_size_sec <- 1800

# A dyad is considered direct contact when both animals occupy the same RFID
# PositionID. Adjacent contact is computed separately using the 2 x 4 cage grid.
contact_definition <- "same_position"

# Long gaps are retained in the interval-level export, but excluded from the
# network-ready bin aggregation by default.
long_gap_threshold_sec <- 3600
exclude_long_gaps_from_aggregation <- TRUE

# Parallelize file reading and per-system interval reconstruction. Set
# use_multicore <- FALSE for easier debugging or very small datasets.
use_multicore <- TRUE
n_cores <- max(1L, parallel::detectCores(logical = TRUE) - 1L)

if (use_multicore && requireNamespace("furrr", quietly = TRUE) && requireNamespace("future", quietly = TRUE)) {
  future::plan(future::multisession, workers = n_cores)
  map_dfr_parallel <- furrr::future_map_dfr
  message("Parallel processing enabled with ", n_cores, " workers.")
} else {
  map_dfr_parallel <- purrr::map_dfr
  message("Parallel processing disabled; using sequential purrr.")
}

# Optional metadata file with AnimalID/AnimalNum plus Group/Sex columns.
# Leave NULL if Group and Sex are unavailable at this stage.
metadata_file <- NULL

output_dirs <- analysis_output_dirs(output_dir)
write_output_manifest(
  output_dir,
  script_name = "05_build_dyadic_rfid_contacts.R",
  analysis_name = "dyadic RFID contact table",
  primary_tables = c(
    "tables/dyadic_contacts_by_bin.csv",
    "tables/dyadic_contacts_interval_level.csv",
    "tables/dyadic_network_ready.csv",
    "tables/dyadic_contact_qc_by_file.csv",
    "tables/epoch_duration_qc.csv"
  ),
  notes = c("Canonical output folders are tables/, stats_tables/, qc/, and figures/{publication_panels,supplementary,qc,exploratory}.")
)

# ------------------------------------------------
# POSITION GEOMETRY
# ------------------------------------------------

position_grid <- tibble(
  PositionID = 1:8,
  x_grid = rep(1:4, times = 2),
  y_grid = rep(1:2, each = 4)
)

pair_distance_lookup <- tidyr::crossing(
  PositionID_A = position_grid$PositionID,
  PositionID_B = position_grid$PositionID
) %>%
  left_join(position_grid %>% rename(PositionID_A = PositionID, x_A = x_grid, y_A = y_grid), by = "PositionID_A") %>%
  left_join(position_grid %>% rename(PositionID_B = PositionID, x_B = x_grid, y_B = y_grid), by = "PositionID_B") %>%
  mutate(
    grid_distance = abs(x_A - x_B) + abs(y_A - y_B),
    same_position = PositionID_A == PositionID_B,
    adjacent_position = grid_distance == 1
  ) %>%
  select(PositionID_A, PositionID_B, grid_distance, same_position, adjacent_position)

# ------------------------------------------------
# HELPERS
# ------------------------------------------------

read_preprocessed_position_file <- function(path) {
  dat <- readr::read_csv(path, show_col_types = FALSE)

  required <- c("DateTime", "AnimalID", "System", "PositionID")
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0) {
    stop("Missing required columns in ", path, ": ", paste(missing, collapse = ", "), call. = FALSE)
  }

  dat %>%
    mutate(
      DateTime = suppressWarnings(as.POSIXct(DateTime, tz = "UTC")),
      AnimalID = as.character(AnimalID),
      System = as.character(System),
      PositionID = suppressWarnings(as.integer(PositionID)),
      SourceFile = basename(path),
      Batch = if ("Batch" %in% names(.)) as.character(Batch) else str_extract(basename(path), "B[0-9]+"),
      CageChange = if ("CageChange" %in% names(.)) as.character(CageChange) else str_extract(basename(path), "CC[0-9]+"),
      Phase = if ("Phase" %in% names(.)) as.character(Phase) else if_else(
        format(DateTime, "%H:%M", tz = "UTC") >= "18:30" | format(DateTime, "%H:%M", tz = "UTC") < "06:30",
        "Active", "Inactive"
      ),
      Group = if ("Group" %in% names(.)) as.character(Group) else NA_character_,
      Sex = if ("Sex" %in% names(.)) as.character(Sex) else NA_character_,
      HalfHoursElapsed = if ("HalfHoursElapsed" %in% names(.)) suppressWarnings(as.numeric(HalfHoursElapsed)) else NA_real_
    ) %>%
    filter(!is.na(DateTime), !is.na(AnimalID), !is.na(System), is.finite(PositionID), PositionID > 0) %>%
    arrange(SourceFile, System, DateTime, AnimalID)
}

make_interval_dyads_one_system <- function(dat_sys) {
  dat_sys <- dat_sys %>% arrange(DateTime, AnimalID)
  animals <- sort(unique(dat_sys$AnimalID))

  if (length(animals) < 2) {
    return(tibble())
  }

  event_times <- sort(unique(dat_sys$DateTime))
  if (length(event_times) < 2) {
    return(tibble())
  }

  animal_pairs <- combn(animals, 2, simplify = FALSE)
  out <- vector("list", length(event_times) - 1)

  current_pos <- rep(NA_integer_, length(animals))
  names(current_pos) <- animals
  current_group <- rep(NA_character_, length(animals))
  names(current_group) <- animals
  current_sex <- rep(NA_character_, length(animals))
  names(current_sex) <- animals

  for (i in seq_len(length(event_times) - 1)) {
    t0 <- event_times[i]
    t1 <- event_times[i + 1]
    duration_sec <- as.numeric(difftime(t1, t0, units = "secs"))

    updates <- dat_sys %>% filter(DateTime == t0)
    if (nrow(updates) > 0) {
      for (j in seq_len(nrow(updates))) {
        a <- updates$AnimalID[j]
        current_pos[a] <- updates$PositionID[j]
        if (!is.na(updates$Group[j])) current_group[a] <- updates$Group[j]
        if (!is.na(updates$Sex[j])) current_sex[a] <- updates$Sex[j]
      }
    }

    if (!is.finite(duration_sec) || duration_sec <= 0) next

    # Guard against long gaps caused by acquisition breaks. These intervals are
    # kept but flagged so they can be excluded later if needed.
    long_gap <- duration_sec > long_gap_threshold_sec

    interval_meta <- updates %>%
      slice_tail(n = 1) %>%
      transmute(
        SourceFile = first(SourceFile),
        Batch = first(Batch),
        CageChange = first(CageChange),
        System = first(System),
        Phase = first(Phase),
        HalfHoursElapsed = first(HalfHoursElapsed)
      )

    if (nrow(interval_meta) == 0) {
      interval_meta <- dat_sys %>%
        filter(DateTime <= t0) %>%
        slice_tail(n = 1) %>%
        transmute(
          SourceFile = first(SourceFile),
          Batch = first(Batch),
          CageChange = first(CageChange),
          System = first(System),
          Phase = first(Phase),
          HalfHoursElapsed = first(HalfHoursElapsed)
        )
    }

    pair_tbl <- map_dfr(animal_pairs, function(p) {
      pos_a <- current_pos[p[1]]
      pos_b <- current_pos[p[2]]
      if (!is.finite(pos_a) || !is.finite(pos_b)) return(tibble())
      tibble(
        Focal = p[1],
        Partner = p[2],
        PositionID_A = as.integer(pos_a),
        PositionID_B = as.integer(pos_b),
        Group_Focal = current_group[p[1]],
        Group_Partner = current_group[p[2]],
        Sex_Focal = current_sex[p[1]],
        Sex_Partner = current_sex[p[2]]
      )
    })

    if (nrow(pair_tbl) == 0) next

    out[[i]] <- pair_tbl %>%
      mutate(
        IntervalStart = t0,
        IntervalEnd = t1,
        DurationSec = duration_sec,
        LongGap = long_gap
      ) %>%
      bind_cols(interval_meta[rep(1, nrow(pair_tbl)), ]) %>%
      left_join(pair_distance_lookup, by = c("PositionID_A", "PositionID_B"))
  }

  bind_rows(out)
}

filter_aggregation_intervals <- function(dat) {
  if (exclude_long_gaps_from_aggregation && "LongGap" %in% names(dat)) {
    dat <- dat %>% filter(!LongGap)
  }
  dat
}

split_intervals_to_bins <- function(dat) {
  dat <- dat %>%
    filter_aggregation_intervals() %>%
    filter(
      !is.na(IntervalStart),
      !is.na(IntervalEnd),
      is.finite(DurationSec),
      DurationSec > 0
    )

  if (nrow(dat) == 0) {
    return(dat %>%
             mutate(
               BinStart = as.POSIXct(NA_real_, origin = "1970-01-01", tz = "UTC")
             ) %>%
             slice(0))
  }

  interval_start_num <- as.numeric(dat$IntervalStart)
  interval_end_num <- as.numeric(dat$IntervalEnd)
  bin_start_num <- floor(interval_start_num / bin_size_sec) * bin_size_sec
  bin_end_num <- floor((interval_end_num - 1e-7) / bin_size_sec) * bin_size_sec
  n_bins <- pmax(0L, as.integer((bin_end_num - bin_start_num) / bin_size_sec) + 1L)

  keep <- n_bins > 0
  if (!all(keep)) {
    dat <- dat[keep, , drop = FALSE]
    interval_start_num <- interval_start_num[keep]
    interval_end_num <- interval_end_num[keep]
    bin_start_num <- bin_start_num[keep]
    n_bins <- n_bins[keep]
  }

  row_idx <- rep(seq_len(nrow(dat)), n_bins)
  bin_offsets <- sequence(n_bins) - 1L
  split_bin_start_num <- rep(bin_start_num, n_bins) + bin_offsets * bin_size_sec
  split_start_num <- pmax(rep(interval_start_num, n_bins), split_bin_start_num)
  split_end_num <- pmin(rep(interval_end_num, n_bins), split_bin_start_num + bin_size_sec)
  split_duration_sec <- split_end_num - split_start_num

  out <- dat[row_idx, , drop = FALSE]
  out$IntervalStart <- as.POSIXct(split_start_num, origin = "1970-01-01", tz = "UTC")
  out$IntervalEnd <- as.POSIXct(split_end_num, origin = "1970-01-01", tz = "UTC")
  out$DurationSec <- split_duration_sec
  out$BinStart <- as.POSIXct(split_bin_start_num, origin = "1970-01-01", tz = "UTC")

  out %>% filter(is.finite(DurationSec), DurationSec > 0)
}

aggregate_dyads_by_bin <- function(interval_tbl) {
  if (nrow(interval_tbl) == 0) return(tibble())

  interval_tbl %>%
    split_intervals_to_bins() %>%
    mutate(
      same_position_seconds = if_else(same_position, DurationSec, 0),
      adjacent_seconds = if_else(adjacent_position, DurationSec, 0),
      weighted_grid_distance = grid_distance * DurationSec
    ) %>%
    group_by(SourceFile, Batch, CageChange, System) %>%
    mutate(TimeIndex = as.numeric(difftime(BinStart, min(BinStart, na.rm = TRUE), units = "secs")) / bin_size_sec) %>%
    ungroup() %>%
    group_by(SourceFile, Batch, CageChange, System, Phase, BinStart, TimeIndex, Focal, Partner) %>%
    summarise(
      observation_seconds = sum(DurationSec, na.rm = TRUE),
      same_position_seconds = sum(same_position_seconds, na.rm = TRUE),
      adjacent_seconds = sum(adjacent_seconds, na.rm = TRUE),
      mean_grid_distance = sum(weighted_grid_distance, na.rm = TRUE) / observation_seconds,
      n_intervals = n(),
      n_long_gaps = sum(LongGap, na.rm = TRUE),
      Group = first(na.omit(c(Group_Focal, Group_Partner, NA_character_))),
      Sex = first(na.omit(c(Sex_Focal, Sex_Partner, NA_character_))),
      .groups = "drop"
    ) %>%
    mutate(
      BinSizeSec = bin_size_sec,
      same_position_fraction = same_position_seconds / observation_seconds,
      adjacent_fraction = adjacent_seconds / observation_seconds,
      same_position_seconds_per_hour = same_position_seconds / pmax(observation_seconds / 3600, 1e-9),
      adjacent_seconds_per_hour = adjacent_seconds / pmax(observation_seconds / 3600, 1e-9),
      contact_fraction = case_when(
        contact_definition == "same_position" ~ same_position_fraction,
        contact_definition == "same_or_adjacent" ~ (same_position_seconds + adjacent_seconds) / observation_seconds,
        TRUE ~ same_position_fraction
      ),
      Weight = contact_fraction,
      Contact = Weight > 0,
      Group = if_else(is.na(Group), "All", Group),
      Sex = if_else(is.na(Sex), "All", Sex)
    )
}

# ------------------------------------------------
# LOAD FILES
# ------------------------------------------------

files <- list.files(input_dir, pattern = "_preprocessed\\.csv$", full.names = TRUE)
if (length(files) == 0) {
  stop("No preprocessed files found in ", input_dir, call. = FALSE)
}

message("Found ", length(files), " preprocessed files.")

all_pos <- map_dfr_parallel(files, read_preprocessed_position_file)

# Optional metadata merge
if (!is.null(metadata_file) && file.exists(metadata_file)) {
  meta <- read_behavior_table(metadata_file)
  meta_animal_col <- first_existing_col(meta, c("AnimalID", "AnimalNum", "Animal", "MouseID", "Mouse", "ID", "RFID"), TRUE, "metadata animal column")
  meta_group_col <- first_existing_col(meta, c("Group", "Phenotype", "Condition", "Treatment", "StressGroup"), FALSE, "metadata group column")
  meta_sex_col <- first_existing_col(meta, c("Sex", "sex"), FALSE, "metadata sex column")

  meta_small <- meta %>%
    transmute(
      AnimalID = as.character(.data[[meta_animal_col]]),
      MetaGroup = if (!is.na(meta_group_col)) as.character(.data[[meta_group_col]]) else NA_character_,
      MetaSex = if (!is.na(meta_sex_col)) as.character(.data[[meta_sex_col]]) else NA_character_
    ) %>%
    distinct(AnimalID, .keep_all = TRUE)

  all_pos <- all_pos %>%
    left_join(meta_small, by = "AnimalID") %>%
    mutate(
      Group = coalesce(Group, MetaGroup),
      Sex = coalesce(Sex, MetaSex)
    ) %>%
    select(-MetaGroup, -MetaSex)
}

qc_by_file <- all_pos %>%
  group_by(SourceFile, Batch, CageChange, System) %>%
  summarise(
    n_rows = n(),
    n_animals = n_distinct(AnimalID),
    start_time = min(DateTime, na.rm = TRUE),
    end_time = max(DateTime, na.rm = TRUE),
    duration_hours = as.numeric(difftime(end_time, start_time, units = "hours")),
    missing_position_rows = sum(!is.finite(PositionID)),
    .groups = "drop"
  )

write_table(qc_by_file, file.path(output_dir, "tables", "dyadic_contact_qc_by_file.csv"))

# ------------------------------------------------
# BUILD INTERVAL-LEVEL DYADS
# ------------------------------------------------

message("Building interval-level dyads...")

interval_tbl <- all_pos %>%
  group_by(SourceFile, Batch, CageChange, System) %>%
  group_split() %>%
  map_dfr_parallel(make_interval_dyads_one_system)

if (use_multicore && requireNamespace("future", quietly = TRUE)) {
  future::plan(future::sequential)
  message("Parallel interval reconstruction complete; aggregating bins sequentially.")
}

if (nrow(interval_tbl) == 0) {
  stop("No dyadic intervals could be derived. Check AnimalID/System/PositionID structure.", call. = FALSE)
}

write_table(
  interval_tbl,
  file.path(output_dir, "tables", "dyadic_contacts_interval_level.csv")
)

# ------------------------------------------------
# AGGREGATE TO NETWORK-READY TABLE
# ------------------------------------------------

message("Aggregating dyads to ", bin_size_sec, "-second bins...")

dyad_bin_tbl <- aggregate_dyads_by_bin(interval_tbl)

write_epoch_duration_qc(
  dyad_bin_tbl %>% rename(AnimalNum = Focal),
  output_dir,
  metric_source = "05_dyadic_rfid_contacts",
  bin_size_sec = bin_size_sec
)

write_table(
  dyad_bin_tbl,
  file.path(output_dir, "tables", "dyadic_contacts_by_bin.csv")
)

# A compact file with the exact columns expected by 09_dynamic_social_networks.R
network_ready_tbl <- dyad_bin_tbl %>%
  transmute(
    Focal,
    Partner,
    TimeIndex,
    BinStart,
    Weight,
    Proximity = Weight,
    Contact,
    Group,
    Sex,
    Phase,
    CageChange,
    Batch,
    System,
    BinSizeSec,
    observation_seconds,
    same_position_fraction,
    adjacent_fraction,
    same_position_seconds_per_hour,
    adjacent_seconds_per_hour,
    mean_grid_distance,
    n_intervals,
    n_long_gaps
  )

write_table(
  network_ready_tbl,
  file.path(output_dir, "tables", "dyadic_network_ready.csv")
)

message("Dyadic RFID contact extraction complete.")
message("Network-ready dyad file: ", file.path(output_dir, "tables", "dyadic_network_ready.csv"))
