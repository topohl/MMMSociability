# ================================================================
# Build multiscale animal-level behavior metrics from RFID position data
# MMMSociability
# ================================================================
# Goal:
#   Produce canonical all_behavior_metrics.csv files used by downstream
#   behavioral dynamics scripts.
#
# Input:
#   preprocessed_data/*_preprocessed.csv
#
# Output:
#   analysis_ready/03_derived_metrics/{10sec,1min,5min,10min,30min}_based/all_behavior_metrics.csv
#   analysis_ready/03_derived_metrics/phase_based/all_behavior_metrics.csv
#   analysis_ready/03_derived_metrics/qc/multiscale_behavior_metrics_qc.csv
#   analysis_ready/03_derived_metrics/qc/animal_group_sex_assignment_qc.csv
#   analysis_ready/03_derived_metrics/qc/reference_ids_not_found_in_preprocessed_data.csv
#   analysis_ready/03_derived_metrics/qc/group_sex_assignment_summary.csv
#
# Key definitions:
#   Movement                  = number of RFID position transitions per animal/bin
#   MovementDistance          = summed Manhattan grid distance moved per animal/bin
#   Entropy                   = Shannon entropy of position occupancy seconds
#   ProximitySeconds          = summed same-position dyadic contact seconds
#   ProximityFraction         = ProximitySeconds / dyadic_observation_seconds
#   Proximity                 = backward-compatible alias for ProximityFraction
#   AdjacentProximitySeconds  = summed adjacent-position dyadic seconds
#   AdjacentProximityFraction = AdjacentProximitySeconds / dyadic_observation_seconds
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(lubridate)
})

source("Functions/behavioral_dynamics_helpers.R")

# ------------------------------------------------
# USER INPUT
# ------------------------------------------------

input_dir <- "preprocessed_data"
output_root <- "analysis_ready/03_derived_metrics"

# Optional animal reference lists. These are one-ID-per-line CSV/text files.
# Matching is done after robust character normalization, so mixed numeric/string
# IDs such as 4, 303, OQ754, OR111 are handled consistently.
sus_animals_file <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/sus_animals.csv"
con_animals_file <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/con_animals.csv"

# If TRUE, animals not listed as CON or SUS are assigned RES. This is appropriate
# only when all remaining animals in the preprocessed files are stress-exposed
# resilient animals.
assign_unlisted_animals_as_res <- TRUE

bin_specs <- tibble::tribble(
  ~bin_label, ~bin_size_sec,
  "10sec",    10,
  "1min",     60,
  "5min",     300,
  "10min",    600,
  "30min",    1800
)

metadata_file <- NULL
long_gap_threshold_sec <- 3600

ensure_dir(output_root)
ensure_dir(file.path(output_root, "qc"))

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

calc_entropy <- function(seconds_by_position) {
  x <- seconds_by_position[is.finite(seconds_by_position) & seconds_by_position > 0]
  if (length(x) == 0 || sum(x) <= 0) return(NA_real_)
  p <- x / sum(x)
  -sum(p * log2(p))
}

safe_divide <- function(num, den) {
  ifelse(is.finite(den) & den > 0, num / den, NA_real_)
}

normalize_animal_id <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_trim() %>%
    stringr::str_replace_all("\\s+", "") %>%
    toupper()
}

read_animal_id_list <- function(path, label) {
  if (is.null(path) || !file.exists(path)) {
    warning("Animal reference file not found for ", label, ": ", path, call. = FALSE)
    return(character())
  }

  readr::read_lines(path, progress = FALSE) %>%
    normalize_animal_id() %>%
    purrr::discard(~ is.na(.x) || .x == "") %>%
    unique()
}

assign_batch_sex <- function(batch) {
  batch_norm <- toupper(stringr::str_trim(as.character(batch)))
  dplyr::case_when(
    batch_norm %in% c("B1", "B2", "B5") ~ "Male",
    batch_norm %in% c("B3", "B4", "B6") ~ "Female",
    TRUE ~ NA_character_
  )
}

# ------------------------------------------------
# LOAD PREPROCESSED POSITION DATA
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
      AnimalNum = AnimalID,
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
      ConsecActive = if ("ConsecActive" %in% names(.)) suppressWarnings(as.integer(ConsecActive)) else NA_integer_,
      ConsecInactive = if ("ConsecInactive" %in% names(.)) suppressWarnings(as.integer(ConsecInactive)) else NA_integer_,
      HalfHoursElapsed = if ("HalfHoursElapsed" %in% names(.)) suppressWarnings(as.numeric(HalfHoursElapsed)) else NA_real_
    ) %>%
    filter(!is.na(DateTime), !is.na(AnimalID), !is.na(System), is.finite(PositionID), PositionID > 0) %>%
    arrange(SourceFile, Batch, CageChange, System, DateTime, AnimalID)
}

files <- list.files(input_dir, pattern = "_preprocessed\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("No preprocessed files found in ", input_dir, call. = FALSE)

message("Found ", length(files), " preprocessed files.")
all_pos <- map_dfr(files, read_preprocessed_position_file)

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
    mutate(Group = coalesce(Group, MetaGroup), Sex = coalesce(Sex, MetaSex)) %>%
    select(-MetaGroup, -MetaSex)
}

# ------------------------------------------------
# ASSIGN GROUP AND SEX METADATA
# ------------------------------------------------

sus_animals <- read_animal_id_list(sus_animals_file, "SUS")
con_animals <- read_animal_id_list(con_animals_file, "CON")

reference_overlap <- intersect(sus_animals, con_animals)
if (length(reference_overlap) > 0) {
  stop(
    "The following normalized AnimalIDs occur in both SUS and CON reference files: ",
    paste(reference_overlap, collapse = ", "),
    call. = FALSE
  )
}

all_pos <- all_pos %>%
  mutate(
    AnimalID_raw = AnimalID,
    AnimalID_norm = normalize_animal_id(AnimalID),
    Batch_norm = toupper(str_trim(Batch)),
    ReferenceGroup = case_when(
      AnimalID_norm %in% sus_animals ~ "SUS",
      AnimalID_norm %in% con_animals ~ "CON",
      assign_unlisted_animals_as_res ~ "RES",
      TRUE ~ NA_character_
    ),
    ReferenceSex = assign_batch_sex(Batch),
    Group = coalesce(ReferenceGroup, Group),
    Sex = coalesce(ReferenceSex, Sex)
  )

animal_group_sex_qc <- all_pos %>%
  distinct(AnimalID_raw, AnimalID_norm, Batch, Batch_norm, Sex, Group, ReferenceGroup, ReferenceSex) %>%
  arrange(Batch_norm, Group, AnimalID_norm)

write_table(
  animal_group_sex_qc,
  file.path(output_root, "qc", "animal_group_sex_assignment_qc.csv")
)

reference_ids_not_found <- tibble(
  AnimalID_norm = c(sus_animals, con_animals),
  ReferenceGroup = c(rep("SUS", length(sus_animals)), rep("CON", length(con_animals)))
) %>%
  distinct() %>%
  anti_join(all_pos %>% distinct(AnimalID_norm), by = "AnimalID_norm") %>%
  arrange(ReferenceGroup, AnimalID_norm)

write_table(
  reference_ids_not_found,
  file.path(output_root, "qc", "reference_ids_not_found_in_preprocessed_data.csv")
)

group_sex_assignment_summary <- animal_group_sex_qc %>%
  count(Batch_norm, Sex, Group, name = "n_animals") %>%
  arrange(Batch_norm, Sex, Group)

write_table(
  group_sex_assignment_summary,
  file.path(output_root, "qc", "group_sex_assignment_summary.csv")
)

if (nrow(reference_ids_not_found) > 0) {
  warning(
    nrow(reference_ids_not_found),
    " reference AnimalIDs were not found in the preprocessed data. Check qc/reference_ids_not_found_in_preprocessed_data.csv",
    call. = FALSE
  )
}

if (any(is.na(animal_group_sex_qc$Sex))) {
  warning("Some animals could not be assigned Sex from Batch. Check qc/animal_group_sex_assignment_qc.csv", call. = FALSE)
}

# ------------------------------------------------
# INTERVAL RECONSTRUCTION
# ------------------------------------------------

last_meta_before <- function(dat_sys, t0) {
  dat_sys %>%
    filter(DateTime <= t0) %>%
    slice_tail(n = 1) %>%
    transmute(
      SourceFile = first(SourceFile),
      Batch = first(Batch),
      CageChange = first(CageChange),
      System = first(System),
      Phase = first(Phase),
      ConsecActive = first(ConsecActive),
      ConsecInactive = first(ConsecInactive)
    )
}

make_occupancy_intervals_one_system <- function(dat_sys) {
  dat_sys <- dat_sys %>% arrange(DateTime, AnimalID)
  animals <- sort(unique(dat_sys$AnimalID))
  event_times <- sort(unique(dat_sys$DateTime))
  if (length(animals) == 0 || length(event_times) < 2) return(tibble())

  out <- vector("list", length(event_times) - 1)
  current_pos <- rep(NA_integer_, length(animals)); names(current_pos) <- animals
  current_group <- rep(NA_character_, length(animals)); names(current_group) <- animals
  current_sex <- rep(NA_character_, length(animals)); names(current_sex) <- animals

  for (i in seq_len(length(event_times) - 1)) {
    t0 <- event_times[i]
    t1 <- event_times[i + 1]
    duration_sec <- as.numeric(difftime(t1, t0, units = "secs"))
    if (!is.finite(duration_sec) || duration_sec <= 0) next

    updates <- dat_sys %>% filter(DateTime == t0)
    if (nrow(updates) > 0) {
      for (j in seq_len(nrow(updates))) {
        a <- updates$AnimalID[j]
        current_pos[a] <- updates$PositionID[j]
        if (!is.na(updates$Group[j])) current_group[a] <- updates$Group[j]
        if (!is.na(updates$Sex[j])) current_sex[a] <- updates$Sex[j]
      }
    }

    valid <- is.finite(current_pos) & current_pos > 0
    if (!any(valid)) next

    interval_meta <- last_meta_before(dat_sys, t0)
    out[[i]] <- tibble(
      AnimalNum = names(current_pos)[valid],
      AnimalID = names(current_pos)[valid],
      PositionID = as.integer(current_pos[valid]),
      Group = current_group[names(current_pos)[valid]],
      Sex = current_sex[names(current_pos)[valid]],
      IntervalStart = t0,
      IntervalEnd = t1,
      DurationSec = duration_sec,
      LongGap = duration_sec > long_gap_threshold_sec
    ) %>%
      bind_cols(interval_meta[rep(1, sum(valid)), ])
  }

  bind_rows(out)
}

make_movement_events <- function(pos_tbl) {
  pos_tbl %>%
    arrange(SourceFile, Batch, CageChange, System, AnimalID, DateTime) %>%
    group_by(SourceFile, Batch, CageChange, System, AnimalID) %>%
    mutate(
      PrevPositionID = lag(PositionID),
      PositionChanged = is.finite(PrevPositionID) & PositionID != PrevPositionID
    ) %>%
    ungroup() %>%
    filter(PositionChanged) %>%
    left_join(
      pair_distance_lookup %>% rename(PrevPositionID = PositionID_A, PositionID = PositionID_B),
      by = c("PrevPositionID", "PositionID")
    ) %>%
    transmute(
      SourceFile, Batch, CageChange, System, Phase, ConsecActive, ConsecInactive,
      AnimalNum = AnimalID,
      AnimalID,
      Group,
      Sex,
      EventTime = DateTime,
      MovementEvent = 1,
      MovementDistanceEvent = grid_distance
    )
}

make_dyadic_intervals_one_system <- function(dat_sys) {
  dat_sys <- dat_sys %>% arrange(DateTime, AnimalID)
  animals <- sort(unique(dat_sys$AnimalID))
  event_times <- sort(unique(dat_sys$DateTime))
  if (length(animals) < 2 || length(event_times) < 2) return(tibble())

  animal_pairs <- combn(animals, 2, simplify = FALSE)
  out <- vector("list", length(event_times) - 1)
  current_pos <- rep(NA_integer_, length(animals)); names(current_pos) <- animals

  for (i in seq_len(length(event_times) - 1)) {
    t0 <- event_times[i]
    t1 <- event_times[i + 1]
    duration_sec <- as.numeric(difftime(t1, t0, units = "secs"))
    if (!is.finite(duration_sec) || duration_sec <= 0) next

    updates <- dat_sys %>% filter(DateTime == t0)
    if (nrow(updates) > 0) {
      for (j in seq_len(nrow(updates))) {
        current_pos[updates$AnimalID[j]] <- updates$PositionID[j]
      }
    }

    pair_tbl <- map_dfr(animal_pairs, function(p) {
      pos_a <- current_pos[p[1]]
      pos_b <- current_pos[p[2]]
      if (!is.finite(pos_a) || !is.finite(pos_b)) return(tibble())
      tibble(Focal = p[1], Partner = p[2], PositionID_A = as.integer(pos_a), PositionID_B = as.integer(pos_b))
    })
    if (nrow(pair_tbl) == 0) next

    interval_meta <- last_meta_before(dat_sys, t0)
    out[[i]] <- pair_tbl %>%
      mutate(
        IntervalStart = t0,
        IntervalEnd = t1,
        DurationSec = duration_sec,
        LongGap = duration_sec > long_gap_threshold_sec
      ) %>%
      bind_cols(interval_meta[rep(1, nrow(pair_tbl)), ]) %>%
      left_join(pair_distance_lookup, by = c("PositionID_A", "PositionID_B"))
  }

  bind_rows(out)
}

message("Building occupancy intervals...")
occupancy_intervals <- all_pos %>%
  group_by(SourceFile, Batch, CageChange, System) %>%
  group_split() %>%
  map_dfr(make_occupancy_intervals_one_system)

message("Building movement events...")
movement_events <- make_movement_events(all_pos)

message("Building dyadic proximity intervals...")
dyadic_intervals <- all_pos %>%
  group_by(SourceFile, Batch, CageChange, System) %>%
  group_split() %>%
  map_dfr(make_dyadic_intervals_one_system)

# ------------------------------------------------
# AGGREGATION HELPERS
# ------------------------------------------------

add_time_bin <- function(dat, time_col, bin_size_sec) {
  dat %>%
    mutate(
      BinStart = as.POSIXct(
        floor(as.numeric(.data[[time_col]]) / bin_size_sec) * bin_size_sec,
        origin = "1970-01-01",
        tz = "UTC"
      )
    )
}

make_time_index <- function(dat) {
  dat %>%
    group_by(SourceFile, Batch, CageChange, System) %>%
    mutate(TimeIndex = as.numeric(difftime(BinStart, min(BinStart, na.rm = TRUE), units = "secs")) / first(BinSizeSec)) %>%
    ungroup()
}

summarise_occupancy_by_bin <- function(bin_size_sec) {
  occupancy_intervals %>%
    add_time_bin("IntervalStart", bin_size_sec) %>%
    mutate(BinSizeSec = bin_size_sec) %>%
    group_by(SourceFile, Batch, CageChange, System, Phase, BinSizeSec, BinStart, AnimalNum, AnimalID, Group, Sex, PositionID) %>%
    summarise(PositionSeconds = sum(DurationSec, na.rm = TRUE), .groups = "drop") %>%
    group_by(SourceFile, Batch, CageChange, System, Phase, BinSizeSec, BinStart, AnimalNum, AnimalID, Group, Sex) %>%
    summarise(
      observation_seconds = sum(PositionSeconds, na.rm = TRUE),
      Entropy = calc_entropy(PositionSeconds),
      DominantPosition = PositionID[which.max(PositionSeconds)],
      n_positions_visited = sum(PositionSeconds > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    make_time_index()
}

summarise_movement_by_bin <- function(bin_size_sec) {
  if (nrow(movement_events) == 0) return(tibble())
  movement_events %>%
    add_time_bin("EventTime", bin_size_sec) %>%
    mutate(BinSizeSec = bin_size_sec) %>%
    group_by(SourceFile, Batch, CageChange, System, Phase, BinSizeSec, BinStart, AnimalNum, AnimalID) %>%
    summarise(
      Movement = sum(MovementEvent, na.rm = TRUE),
      MovementDistance = sum(MovementDistanceEvent, na.rm = TRUE),
      .groups = "drop"
    )
}

summarise_proximity_by_bin <- function(bin_size_sec) {
  if (nrow(dyadic_intervals) == 0) return(tibble())

  dyad_long <- dyadic_intervals %>%
    mutate(
      same_position_seconds = if_else(same_position, DurationSec, 0),
      adjacent_seconds = if_else(adjacent_position, DurationSec, 0),
      weighted_grid_distance = grid_distance * DurationSec
    ) %>%
    select(SourceFile, Batch, CageChange, System, Phase, IntervalStart, DurationSec,
           same_position_seconds, adjacent_seconds, weighted_grid_distance, Focal, Partner) %>%
    bind_rows(
      dyadic_intervals %>%
        mutate(
          same_position_seconds = if_else(same_position, DurationSec, 0),
          adjacent_seconds = if_else(adjacent_position, DurationSec, 0),
          weighted_grid_distance = grid_distance * DurationSec
        ) %>%
        transmute(SourceFile, Batch, CageChange, System, Phase, IntervalStart, DurationSec,
                  same_position_seconds, adjacent_seconds, weighted_grid_distance,
                  Focal = Partner, Partner = Focal)
    )

  dyad_long %>%
    add_time_bin("IntervalStart", bin_size_sec) %>%
    mutate(BinSizeSec = bin_size_sec) %>%
    group_by(SourceFile, Batch, CageChange, System, Phase, BinSizeSec, BinStart, AnimalNum = Focal) %>%
    summarise(
      dyadic_observation_seconds = sum(DurationSec, na.rm = TRUE),
      ProximitySeconds = sum(same_position_seconds, na.rm = TRUE),
      AdjacentProximitySeconds = sum(adjacent_seconds, na.rm = TRUE),
      ProximityFraction = safe_divide(ProximitySeconds, dyadic_observation_seconds),
      AdjacentProximityFraction = safe_divide(AdjacentProximitySeconds, dyadic_observation_seconds),
      Proximity = ProximityFraction,
      AdjacentProximity = AdjacentProximityFraction,
      MeanGridDistanceToOthers = safe_divide(sum(weighted_grid_distance, na.rm = TRUE), dyadic_observation_seconds),
      n_dyadic_intervals = n(),
      .groups = "drop"
    )
}

standardize_metric_output <- function(out, bin_label) {
  out %>%
    mutate(
      Movement = replace_na(Movement, 0),
      MovementDistance = replace_na(MovementDistance, 0),
      ProximitySeconds = replace_na(ProximitySeconds, NA_real_),
      ProximityFraction = replace_na(ProximityFraction, NA_real_),
      AdjacentProximitySeconds = replace_na(AdjacentProximitySeconds, NA_real_),
      AdjacentProximityFraction = replace_na(AdjacentProximityFraction, NA_real_),
      Proximity = ProximityFraction,
      AdjacentProximity = AdjacentProximityFraction,
      MeanGridDistanceToOthers = replace_na(MeanGridDistanceToOthers, NA_real_),
      Group = if_else(is.na(Group) | Group == "", "All", Group),
      Sex = if_else(is.na(Sex) | Sex == "", "All", Sex),
      BinLabel = bin_label
    )
}

build_bin_metrics <- function(bin_label, bin_size_sec) {
  message("Aggregating animal-level metrics for ", bin_label, " bins...")

  occ <- summarise_occupancy_by_bin(bin_size_sec)
  mov <- summarise_movement_by_bin(bin_size_sec)
  prox <- summarise_proximity_by_bin(bin_size_sec)

  out <- occ %>%
    left_join(mov, by = c("SourceFile", "Batch", "CageChange", "System", "Phase", "BinSizeSec", "BinStart", "AnimalNum", "AnimalID")) %>%
    left_join(prox, by = c("SourceFile", "Batch", "CageChange", "System", "Phase", "BinSizeSec", "BinStart", "AnimalNum")) %>%
    standardize_metric_output(bin_label) %>%
    select(
      AnimalNum, AnimalID, Batch, CageChange, System, Group, Sex, Phase,
      BinLabel, BinSizeSec, BinStart, TimeIndex,
      Movement, MovementDistance,
      Proximity, ProximitySeconds, ProximityFraction,
      AdjacentProximity, AdjacentProximitySeconds, AdjacentProximityFraction,
      MeanGridDistanceToOthers,
      Entropy, DominantPosition, n_positions_visited,
      observation_seconds, dyadic_observation_seconds, n_dyadic_intervals,
      SourceFile
    ) %>%
    arrange(Batch, CageChange, System, AnimalNum, BinStart)

  out_dir <- file.path(output_root, paste0(bin_label, "_based"))
  ensure_dir(out_dir)
  write_table(out, file.path(out_dir, "all_behavior_metrics.csv"))
  out
}

build_phase_metrics <- function() {
  message("Aggregating animal-level metrics by phase...")

  occ <- occupancy_intervals %>%
    mutate(PhaseNumber = case_when(Phase == "Active" ~ ConsecActive, Phase == "Inactive" ~ ConsecInactive, TRUE ~ NA_integer_)) %>%
    group_by(SourceFile, Batch, CageChange, System, Phase, PhaseNumber, AnimalNum, AnimalID, Group, Sex, PositionID) %>%
    summarise(PositionSeconds = sum(DurationSec, na.rm = TRUE), .groups = "drop") %>%
    group_by(SourceFile, Batch, CageChange, System, Phase, PhaseNumber, AnimalNum, AnimalID, Group, Sex) %>%
    summarise(
      observation_seconds = sum(PositionSeconds, na.rm = TRUE),
      Entropy = calc_entropy(PositionSeconds),
      DominantPosition = PositionID[which.max(PositionSeconds)],
      n_positions_visited = sum(PositionSeconds > 0, na.rm = TRUE),
      .groups = "drop"
    )

  mov <- movement_events %>%
    mutate(PhaseNumber = case_when(Phase == "Active" ~ ConsecActive, Phase == "Inactive" ~ ConsecInactive, TRUE ~ NA_integer_)) %>%
    group_by(SourceFile, Batch, CageChange, System, Phase, PhaseNumber, AnimalNum, AnimalID) %>%
    summarise(Movement = sum(MovementEvent, na.rm = TRUE), MovementDistance = sum(MovementDistanceEvent, na.rm = TRUE), .groups = "drop")

  prox <- dyadic_intervals %>%
    mutate(
      PhaseNumber = case_when(Phase == "Active" ~ ConsecActive, Phase == "Inactive" ~ ConsecInactive, TRUE ~ NA_integer_),
      same_position_seconds = if_else(same_position, DurationSec, 0),
      adjacent_seconds = if_else(adjacent_position, DurationSec, 0),
      weighted_grid_distance = grid_distance * DurationSec
    ) %>%
    select(SourceFile, Batch, CageChange, System, Phase, PhaseNumber, DurationSec,
           same_position_seconds, adjacent_seconds, weighted_grid_distance, Focal, Partner) %>%
    bind_rows(
      dyadic_intervals %>%
        mutate(
          PhaseNumber = case_when(Phase == "Active" ~ ConsecActive, Phase == "Inactive" ~ ConsecInactive, TRUE ~ NA_integer_),
          same_position_seconds = if_else(same_position, DurationSec, 0),
          adjacent_seconds = if_else(adjacent_position, DurationSec, 0),
          weighted_grid_distance = grid_distance * DurationSec
        ) %>%
        transmute(SourceFile, Batch, CageChange, System, Phase, PhaseNumber, DurationSec,
                  same_position_seconds, adjacent_seconds, weighted_grid_distance,
                  Focal = Partner, Partner = Focal)
    ) %>%
    group_by(SourceFile, Batch, CageChange, System, Phase, PhaseNumber, AnimalNum = Focal) %>%
    summarise(
      dyadic_observation_seconds = sum(DurationSec, na.rm = TRUE),
      ProximitySeconds = sum(same_position_seconds, na.rm = TRUE),
      AdjacentProximitySeconds = sum(adjacent_seconds, na.rm = TRUE),
      ProximityFraction = safe_divide(ProximitySeconds, dyadic_observation_seconds),
      AdjacentProximityFraction = safe_divide(AdjacentProximitySeconds, dyadic_observation_seconds),
      Proximity = ProximityFraction,
      AdjacentProximity = AdjacentProximityFraction,
      MeanGridDistanceToOthers = safe_divide(sum(weighted_grid_distance, na.rm = TRUE), dyadic_observation_seconds),
      n_dyadic_intervals = n(),
      .groups = "drop"
    )

  out <- occ %>%
    left_join(mov, by = c("SourceFile", "Batch", "CageChange", "System", "Phase", "PhaseNumber", "AnimalNum", "AnimalID")) %>%
    left_join(prox, by = c("SourceFile", "Batch", "CageChange", "System", "Phase", "PhaseNumber", "AnimalNum")) %>%
    group_by(SourceFile, Batch, CageChange, System, AnimalNum) %>%
    arrange(Phase, PhaseNumber, .by_group = TRUE) %>%
    mutate(TimeIndex = row_number() - 1L) %>%
    ungroup() %>%
    mutate(BinSizeSec = NA_real_, BinStart = as.POSIXct(NA_real_, origin = "1970-01-01", tz = "UTC")) %>%
    standardize_metric_output("phase") %>%
    select(
      AnimalNum, AnimalID, Batch, CageChange, System, Group, Sex, Phase,
      PhaseNumber, BinLabel, BinSizeSec, BinStart, TimeIndex,
      Movement, MovementDistance,
      Proximity, ProximitySeconds, ProximityFraction,
      AdjacentProximity, AdjacentProximitySeconds, AdjacentProximityFraction,
      MeanGridDistanceToOthers,
      Entropy, DominantPosition, n_positions_visited,
      observation_seconds, dyadic_observation_seconds, n_dyadic_intervals,
      SourceFile
    ) %>%
    arrange(Batch, CageChange, System, AnimalNum, TimeIndex)

  out_dir <- file.path(output_root, "phase_based")
  ensure_dir(out_dir)
  write_table(out, file.path(out_dir, "all_behavior_metrics.csv"))
  out
}

# ------------------------------------------------
# WRITE OUTPUTS
# ------------------------------------------------

all_bin_outputs <- pmap(list(bin_specs$bin_label, bin_specs$bin_size_sec), build_bin_metrics)
phase_output <- build_phase_metrics()

qc_tbl <- bind_rows(
  map2_dfr(bin_specs$bin_label, all_bin_outputs, function(label, dat) {
    tibble(
      output = paste0(label, "_based"),
      n_rows = nrow(dat),
      n_animals = n_distinct(dat$AnimalNum),
      n_batches = n_distinct(dat$Batch),
      n_cage_changes = n_distinct(dat$CageChange),
      n_systems = n_distinct(dat$System),
      total_observation_hours = sum(dat$observation_seconds, na.rm = TRUE) / 3600,
      missing_proximity_fraction = mean(is.na(dat$ProximityFraction)),
      missing_proximity_seconds = mean(is.na(dat$ProximitySeconds))
    )
  }),
  tibble(
    output = "phase_based",
    n_rows = nrow(phase_output),
    n_animals = n_distinct(phase_output$AnimalNum),
    n_batches = n_distinct(phase_output$Batch),
    n_cage_changes = n_distinct(phase_output$CageChange),
    n_systems = n_distinct(phase_output$System),
    total_observation_hours = sum(phase_output$observation_seconds, na.rm = TRUE) / 3600,
    missing_proximity_fraction = mean(is.na(phase_output$ProximityFraction)),
    missing_proximity_seconds = mean(is.na(phase_output$ProximitySeconds))
  )
)

write_table(qc_tbl, file.path(output_root, "qc", "multiscale_behavior_metrics_qc.csv"))

message("Multiscale behavior metric export complete.")
message("Primary downstream file pattern: ", file.path(output_root, "<scale>_based", "all_behavior_metrics.csv"))
