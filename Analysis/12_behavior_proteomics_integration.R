# ================================================================
# Behavior ↔ Proteomics Integration
# MMMSociability
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(purrr)
})

source("C:/Users/topohl/Documents/GitHub/MMMSociability/Functions/behavioral_dynamics_helpers.R")

# ------------------------------------------------
# 1) Paths
# ------------------------------------------------

behavior_file <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/analysis_ready/06_behavioral_dynamics/early_prediction/5min_based/tables/early_behavior_features.csv"

project_root <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID"

analysis_ready_dir <- file.path(project_root, "analysis_ready")

proteomics_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/analysis_ready/proteomics"

base_output_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/analysis_ready/06_behavioral_dynamics"

behavior_bin_level <- "5min_based"

canonical_behavior_bin_level <- "phase_based"

trajectory_bin_level <- "30min_based"

optional_behavior_bin_levels <- unique(c(behavior_bin_level, trajectory_bin_level, "10min_based", "1min_based"))

proteomics_files <- list.files(
  proteomics_dir,
  pattern = "^module_scores_.*\\.csv$",
  full.names = TRUE
)

proteomics_files <- proteomics_files[
  !stringr::str_detect(basename(proteomics_files), "_metadata\\.csv$")
]

if (length(proteomics_files) == 0) {
  stop("No module-score proteomics files found in: ", proteomics_dir)
}

# ------------------------------------------------
# 2) Helpers
# ------------------------------------------------

normalize_animal_id <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x <- stringr::str_replace(x, "\\.0$", "")

  digits <- stringr::str_extract(x, "\\d+$")
  digits_num <- suppressWarnings(as.integer(digits))

  dplyr::if_else(
    !is.na(digits_num),
    sprintf("%04d", digits_num),
    x
  )
}

standardize_sex <- function(x) {
  x <- as.character(x)
  x <- stringr::str_to_lower(stringr::str_trim(x))

  dplyr::case_when(
    x %in% c("m", "male", "männlich") ~ "male",
    x %in% c("f", "female", "weiblich") ~ "female",
    TRUE ~ x
  )
}

read_optional_table <- function(path) {
  if (is.null(path) || is.na(path) || !file.exists(path)) return(NULL)
  read_behavior_table(path)
}

first_existing_path <- function(candidates) {
  candidates <- candidates[!is.na(candidates)]
  hit <- candidates[file.exists(candidates)][1]
  if (length(hit) == 0 || is.na(hit)) NA_character_ else hit
}

make_feature_name <- function(source, domain, scale, metric, statistic, context = "global") {
  paste(
    safe_name(source),
    safe_name(domain),
    safe_name(scale),
    safe_name(metric),
    safe_name(statistic),
    safe_name(context),
    sep = "__"
  )
}

feature_entropy <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0 || sum(x) <= 0) return(NA_real_)
  p <- x / sum(x)
  -sum(p * log2(p))
}

standardize_feature_table_ids <- function(dat) {
  animal_col <- first_existing_col(
    dat,
    c("AnimalNum", "AnimalID", "Animal", "MouseID", "Mouse", "ID", "RFID", "animal_id"),
    TRUE,
    "behavior animal column"
  )

  group_col <- first_existing_col(
    dat,
    c("Group", "Phenotype", "Condition", "Treatment", "StressGroup"),
    FALSE,
    "behavior group column"
  )

  sex_col <- first_existing_col(
    dat,
    c("Sex", "sex"),
    FALSE,
    "behavior sex column"
  )

  out <- dat %>%
    rename(AnimalNum = all_of(animal_col)) %>%
    mutate(AnimalNum = normalize_animal_id(AnimalNum))

  out$Group <- if (!is.na(group_col)) as.character(dat[[group_col]]) else NA_character_
  out$Sex <- if (!is.na(sex_col)) standardize_sex(dat[[sex_col]]) else NA_character_

  out
}

collapse_behavior_feature_rows <- function(feature_long) {
  if (nrow(feature_long) == 0) {
    return(tibble(AnimalNum = character(), Group = character(), Sex = character()))
  }

  first_informative_label <- function(x) {
    vals <- as.character(x)
    vals <- vals[!is.na(vals) & vals != "" & !str_to_lower(vals) %in% c("all", "na", "nan")]
    if (length(vals) == 0) return(NA_character_)
    vals[1]
  }

  feature_long %>%
    group_by(AnimalNum, Group, Sex, feature) %>%
    summarise(FeatureValue = mean(FeatureValue, na.rm = TRUE), .groups = "drop") %>%
    mutate(FeatureValue = if_else(is.nan(FeatureValue), NA_real_, FeatureValue)) %>%
    filter(is.finite(FeatureValue)) %>%
    pivot_wider(names_from = feature, values_from = FeatureValue) %>%
    group_by(AnimalNum) %>%
    summarise(
      Group = first_informative_label(Group),
      Sex = first_informative_label(Sex),
      across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(across(where(is.numeric), ~ if_else(is.nan(.x), NA_real_, .x)))
}

load_behavior_feature_table <- function(path, source_label, domain_label, scale_label = behavior_bin_level) {
  dat <- read_optional_table(path)

  if (is.null(dat) || nrow(dat) == 0) {
    return(list(
      features = tibble(),
      inventory = tibble(
        Source = source_label,
        Domain = domain_label,
        Scale = scale_label,
        Path = path,
        Status = "missing_or_empty",
        n_rows = 0L,
        n_animals = 0L,
        n_numeric_features = 0L
      )
    ))
  }

  dat <- standardize_feature_table_ids(dat)
  scale_vec <- if ("BinLevel" %in% names(dat)) as.character(dat$BinLevel) else rep(scale_label, nrow(dat))

  descriptor_cols <- intersect(
    c(
      "Metric", "StateLabel", "State", "NextState", "Phase", "PhaseClass",
      "CageChange", "CageChangeIndex", "Window", "System", "BinLevel",
      "ProximityInput", "Context", "Source"
    ),
    names(dat)
  )

  numeric_cols <- dat %>%
    select(where(is.numeric)) %>%
    names()

  numeric_cols <- setdiff(
    numeric_cols,
    c(
      "TimeIndex", "WithinEpochRank", "WithinPhaseRank", "EpochN", "PhaseN",
      "BinSizeSec", "bin_size_sec", "DurationSec", "observation_seconds",
      "total_observation_duration_hours", "dyadic_observation_seconds",
      "CageChangeIndex", "n", "n_bins", "n_obs", "n_animals", "n_rows",
      "n_intervals", "n_dyadic_intervals", "first_time_index", "last_time_index"
    )
  )

  if (length(numeric_cols) == 0) {
    return(list(
      features = tibble(),
      inventory = tibble(
        Source = source_label,
        Domain = domain_label,
        Scale = scale_label,
        Path = path,
        Status = "no_numeric_features",
        n_rows = nrow(dat),
        n_animals = n_distinct(dat$AnimalNum),
        n_numeric_features = 0L
      )
    ))
  }

  out <- dat %>%
    select(AnimalNum, Group, Sex, any_of(descriptor_cols), all_of(numeric_cols)) %>%
    mutate(ScaleName = scale_vec) %>%
    mutate(across(any_of(descriptor_cols), ~ safe_name(as.character(.x)))) %>%
    pivot_longer(all_of(numeric_cols), names_to = "MetricName", values_to = "FeatureValue")

  if (length(descriptor_cols) > 0) {
    out <- out %>%
      tidyr::unite("ContextTag", all_of(descriptor_cols), sep = "_", remove = FALSE, na.rm = TRUE)
  } else {
    out <- out %>%
      mutate(ContextTag = "animal_level")
  }

  out <- out %>%
    mutate(
      ContextTag = if_else(is.na(ContextTag) | ContextTag == "", "animal_level", ContextTag),
      feature = make_feature_name(source_label, domain_label, ScaleName, MetricName, "mean", ContextTag),
      FeatureValue = suppressWarnings(as.numeric(FeatureValue))
    ) %>%
    filter(is.finite(FeatureValue)) %>%
    group_by(AnimalNum, Group, Sex, feature) %>%
    summarise(FeatureValue = mean(FeatureValue, na.rm = TRUE), .groups = "drop") %>%
    select(AnimalNum, Group, Sex, feature, FeatureValue)

  list(
    features = out,
    inventory = tibble(
      Source = source_label,
      Domain = domain_label,
      Scale = scale_label,
      Path = path,
      Status = "loaded",
      n_rows = nrow(dat),
      n_animals = n_distinct(dat$AnimalNum),
      n_numeric_features = length(numeric_cols),
      n_exported_features = n_distinct(out$feature)
    )
  )
}

load_curated_behavior_table <- function(path,
                                        source_label,
                                        domain_label,
                                        scale_label = behavior_bin_level,
                                        numeric_keep,
                                        descriptor_keep = character()) {
  dat <- read_optional_table(path)

  if (is.null(dat) || nrow(dat) == 0) {
    return(list(
      features = tibble(),
      inventory = tibble(
        Source = source_label,
        Domain = domain_label,
        Scale = scale_label,
        Path = path,
        Status = "missing_or_empty",
        n_rows = 0L,
        n_animals = 0L,
        n_numeric_features = 0L,
        n_exported_features = 0L
      )
    ))
  }

  dat <- standardize_feature_table_ids(dat)

  descriptor_cols <- intersect(descriptor_keep, names(dat))
  numeric_cols <- intersect(numeric_keep, names(dat))

  if (length(numeric_cols) == 0) {
    return(list(
      features = tibble(),
      inventory = tibble(
        Source = source_label,
        Domain = domain_label,
        Scale = scale_label,
        Path = path,
        Status = "curated_features_absent",
        n_rows = nrow(dat),
        n_animals = n_distinct(dat$AnimalNum),
        n_numeric_features = 0L,
        n_exported_features = 0L
      )
    ))
  }

  out <- dat %>%
    select(AnimalNum, Group, Sex, any_of(descriptor_cols), all_of(numeric_cols)) %>%
    mutate(across(any_of(descriptor_cols), ~ safe_name(as.character(.x)))) %>%
    pivot_longer(all_of(numeric_cols), names_to = "MetricName", values_to = "FeatureValue") %>%
    mutate(FeatureValue = suppressWarnings(as.numeric(FeatureValue)))

  if (length(descriptor_cols) > 0) {
    out <- out %>%
      tidyr::unite("ContextTag", all_of(descriptor_cols), sep = "_", remove = FALSE, na.rm = TRUE)
  } else {
    out <- out %>%
      mutate(ContextTag = "animal_level")
  }

  out <- out %>%
    mutate(
      ContextTag = if_else(is.na(ContextTag) | ContextTag == "", "animal_level", ContextTag),
      feature = make_feature_name(source_label, domain_label, scale_label, MetricName, "animal_mean", ContextTag)
    ) %>%
    filter(is.finite(FeatureValue)) %>%
    group_by(AnimalNum, Group, Sex, feature) %>%
    summarise(FeatureValue = mean(FeatureValue, na.rm = TRUE), .groups = "drop") %>%
    mutate(FeatureValue = if_else(is.nan(FeatureValue), NA_real_, FeatureValue)) %>%
    filter(is.finite(FeatureValue)) %>%
    select(AnimalNum, Group, Sex, feature, FeatureValue)

  list(
    features = out,
    inventory = tibble(
      Source = source_label,
      Domain = domain_label,
      Scale = scale_label,
      Path = path,
      Status = "loaded_curated",
      n_rows = nrow(dat),
      n_animals = n_distinct(dat$AnimalNum),
      n_numeric_features = length(numeric_cols),
      n_exported_features = n_distinct(out$feature),
      BiologicalCollapse = if_else(
        length(descriptor_cols) == 0,
        "animal mean across phases/cage changes/systems",
        paste0("animal mean within: ", paste(descriptor_cols, collapse = ", "))
      )
    )
  )
}

load_hmm_summary_features <- function(scale_label = behavior_bin_level) {
  hmm_dir <- file.path(analysis_ready_dir, "06_behavioral_dynamics", "hmm_states", scale_label, "tables")
  occ <- read_optional_table(file.path(hmm_dir, "hmm_state_occupancy.csv"))
  dwell <- read_optional_table(file.path(hmm_dir, "hmm_state_dwell_times.csv"))
  trans <- read_optional_table(file.path(hmm_dir, "hmm_transition_probabilities.csv"))

  parts <- list()

  if (!is.null(occ) && nrow(occ) > 0) {
    parts$occupancy <- occ %>%
      standardize_feature_table_ids() %>%
      mutate(frac_time = suppressWarnings(as.numeric(frac_time))) %>%
      group_by(AnimalNum, Group, Sex) %>%
      summarise(
        state_occupancy_entropy = feature_entropy(frac_time),
        max_state_fraction = max(frac_time, na.rm = TRUE),
        inactive_state_fraction = mean(frac_time[str_detect(str_to_lower(as.character(State)), "inactive|rest|low")], na.rm = TRUE),
        .groups = "drop"
      )
  }

  if (!is.null(dwell) && nrow(dwell) > 0) {
    parts$dwell <- dwell %>%
      standardize_feature_table_ids() %>%
      group_by(AnimalNum, Group, Sex) %>%
      summarise(
        mean_hmm_dwell_hours = mean(suppressWarnings(as.numeric(mean_dwell_hours)), na.rm = TRUE),
        max_hmm_dwell_hours = max(suppressWarnings(as.numeric(max_dwell_hours)), na.rm = TRUE),
        mean_hmm_dwell_bins = mean(suppressWarnings(as.numeric(mean_dwell_bins)), na.rm = TRUE),
        .groups = "drop"
      )
  }

  if (!is.null(trans) && nrow(trans) > 0) {
    parts$transition <- trans %>%
      standardize_feature_table_ids() %>%
      mutate(
        TransitionProbability = suppressWarnings(as.numeric(TransitionProbability)),
        Transitions = suppressWarnings(as.numeric(Transitions)),
        State = as.character(State),
        NextState = as.character(NextState)
      ) %>%
      group_by(AnimalNum, Group, Sex) %>%
      summarise(
        hmm_self_transition_probability = mean(TransitionProbability[State == NextState], na.rm = TRUE),
        hmm_transition_entropy = feature_entropy(TransitionProbability),
        hmm_state_switch_rate = sum(Transitions[State != NextState], na.rm = TRUE) / sum(Transitions, na.rm = TRUE),
        .groups = "drop"
      )
  }

  parts <- parts[vapply(parts, nrow, integer(1)) > 0]
  if (length(parts) == 0) return(tibble())

  reduce(parts, full_join, by = c("AnimalNum", "Group", "Sex")) %>%
    pivot_longer(-c(AnimalNum, Group, Sex), names_to = "MetricName", values_to = "FeatureValue") %>%
    mutate(
      feature = make_feature_name("hmm", "latent_state", scale_label, MetricName, "animal_level", "summary"),
      FeatureValue = suppressWarnings(as.numeric(FeatureValue))
    ) %>%
    filter(is.finite(FeatureValue)) %>%
    select(AnimalNum, Group, Sex, feature, FeatureValue)
}

select_curated_behavior_feature_long <- function(feature_long, max_features_per_axis = 25) {
  if (nrow(feature_long) == 0) return(feature_long)

  axis_matches <- map_dfr(seq_len(nrow(behavior_axis_specs)), function(i) {
    axis_name <- behavior_axis_specs$Axis[[i]]
    axis_regex <- behavior_axis_specs$Regex[[i]]

    feature_long %>%
      filter(str_detect(feature, regex(axis_regex, ignore_case = TRUE))) %>%
      mutate(BehaviorAxis = axis_name)
  })

  if (nrow(axis_matches) == 0) return(feature_long[0, ])

  feature_rank <- axis_matches %>%
    group_by(BehaviorAxis, feature) %>%
    summarise(
      n_animals = n_distinct(AnimalNum[is.finite(FeatureValue)]),
      feature_sd = sd(FeatureValue, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_animals >= 3, is.finite(feature_sd), feature_sd > 0) %>%
    group_by(BehaviorAxis) %>%
    slice_max(feature_sd, n = max_features_per_axis, with_ties = FALSE) %>%
    ungroup()

  axis_matches %>%
    semi_join(feature_rank, by = c("BehaviorAxis", "feature")) %>%
    distinct(AnimalNum, Group, Sex, feature, FeatureValue, .keep_all = TRUE)
}

build_behavior_feature_matrix <- function() {
  source_specs <- tibble(
    Source = c(
      "raw_multiscale",
      "temporal_instability",
      "state_space",
      "state_space",
      "early_prediction",
      "social_networks",
      "social_networks",
      "hmm_states",
      "hmm_states",
      "hmm_states",
      "gamm_trajectory",
      "nonlinear_systems",
      "adaptation_kinetics",
      "adaptation_kinetics",
      "sleep_like_inactivity",
      "phase_organization",
      "phase_organization",
      "phase_organization",
      "phase_organization",
      "phase_organization"
    ),
    Domain = c(
      "canonical_behavior",
      "temporal_dynamics",
      "latent_state",
      "latent_state",
      "early_adaptation",
      "social_topology",
      "social_topology",
      "latent_state",
      "latent_state",
      "latent_state",
      "trajectory_geometry",
      "nonlinear_dynamics",
      "adaptation_recovery",
      "adaptation_recovery",
      "quiescence_inactivity",
      "phase_organization",
      "phase_organization",
      "phase_organization",
      "phase_organization",
      "phase_organization"
    ),
    Scale = c(
      canonical_behavior_bin_level,
      behavior_bin_level,
      behavior_bin_level,
      behavior_bin_level,
      behavior_bin_level,
      behavior_bin_level,
      behavior_bin_level,
      behavior_bin_level,
      behavior_bin_level,
      behavior_bin_level,
      trajectory_bin_level,
      behavior_bin_level,
      behavior_bin_level,
      behavior_bin_level,
      behavior_bin_level,
      behavior_bin_level,
      behavior_bin_level,
      behavior_bin_level,
      behavior_bin_level,
      behavior_bin_level
    ),
    Path = c(
      file.path(analysis_ready_dir, "03_derived_metrics", canonical_behavior_bin_level, "all_behavior_metrics.csv"),
      first_existing_path(file.path(analysis_ready_dir, "06_behavioral_dynamics", "temporal_instability", optional_behavior_bin_levels, "tables", "temporal_instability_metrics_per_animal_all_metrics.csv")),
      file.path(analysis_ready_dir, "06_behavioral_dynamics", "state_space", behavior_bin_level, "tables", "state_diversity_metrics.csv"),
      file.path(analysis_ready_dir, "06_behavioral_dynamics", "state_space", behavior_bin_level, "tables", "state_switching_metrics.csv"),
      first_existing_path(c(
        behavior_file,
        file.path(analysis_ready_dir, "06_behavioral_dynamics", "early_prediction", behavior_bin_level, "tables", "early_behavior_features_wide.csv"),
        file.path(analysis_ready_dir, "06_behavioral_dynamics", "early_prediction_model_ladder", behavior_bin_level, "tables", "early_behavior_features_wide.csv")
      )),
      file.path(analysis_ready_dir, "06_behavioral_dynamics", "social_networks", behavior_bin_level, "tables", "animal_level_social_dynamics.csv"),
      file.path(analysis_ready_dir, "06_behavioral_dynamics", "social_networks", behavior_bin_level, "tables", "dyadic_node_summary.csv"),
      file.path(analysis_ready_dir, "06_behavioral_dynamics", "hmm_states", behavior_bin_level, "tables", "hmm_state_occupancy.csv"),
      file.path(analysis_ready_dir, "06_behavioral_dynamics", "hmm_states", behavior_bin_level, "tables", "hmm_state_dwell_times.csv"),
      file.path(analysis_ready_dir, "06_behavioral_dynamics", "hmm_states", behavior_bin_level, "tables", "hmm_transition_probabilities.csv"),
      first_existing_path(c(
        file.path(analysis_ready_dir, "06_behavioral_dynamics", "gamm_features", trajectory_bin_level, "tables", "combined_gamm_features.csv"),
        file.path(analysis_ready_dir, "06_behavioral_dynamics", "gamm_trajectory_features", trajectory_bin_level, "tables", "combined_gamm_features.csv")
      )),
      file.path(analysis_ready_dir, "13_nonlinear_systems_dynamics", behavior_bin_level, "derived_data", "animal_level_nonlinear_feature_matrix.csv"),
      file.path(analysis_ready_dir, "15_behavioral_adaptation_kinetics", behavior_bin_level, "tables", "adaptation_kinetics_features.csv"),
      file.path(analysis_ready_dir, "15_behavioral_adaptation_kinetics", behavior_bin_level, "tables", "distance_to_control_trajectories.csv"),
      file.path(analysis_ready_dir, "16_sleep_like_inactivity_metrics", behavior_bin_level, "tables", "sleep_like_inactivity_features.csv"),
      file.path(analysis_ready_dir, "17_ethological_phase_organization", behavior_bin_level, "tables", "phase_contrast_features.csv"),
      file.path(analysis_ready_dir, "17_ethological_phase_organization", behavior_bin_level, "tables", "phase_timing_features.csv"),
      file.path(analysis_ready_dir, "17_ethological_phase_organization", behavior_bin_level, "tables", "phase_fragmentation_features.csv"),
      file.path(analysis_ready_dir, "17_ethological_phase_organization", behavior_bin_level, "tables", "phase_recovery_kinetics.csv"),
      file.path(analysis_ready_dir, "17_ethological_phase_organization", behavior_bin_level, "tables", "phase_predictability_features.csv")
    )
  )

  curated_numeric <- list(
    raw_multiscale = c(
      "MovementPerHour", "MovementDistancePerHour", "ProximityFraction",
      "AdjacentProximityFraction", "MeanGridDistanceToOthers"
    ),
    temporal_instability = c("mean", "cv", "fano", "rmssd", "acf1"),
    state_space_diversity = c("state_entropy", "normalized_state_entropy", "entropy_rate", "max_state_fraction"),
    state_space_switching = c("switch_rate", "stay_rate", "n_switches_per_hour"),
    early_prediction = c(
      "Movement_mean", "Entropy_mean", "Proximity_mean",
      "Movement_rmssd", "Entropy_rmssd", "Proximity_rmssd",
      "Movement_acf1", "Entropy_acf1", "Proximity_acf1",
      "SocialWithdrawal_mean", "PassiveIsolation_mean", "SocialEngagement_mean"
    ),
    social_animal = c(
      "mean_proximity", "proximity_rmssd", "proximity_acf1",
      "contact_fraction", "contact_switch_rate", "contact_entropy",
      "active_isolation", "passive_isolation", "social_engagement",
      "contact_switch_count_per_hour"
    ),
    social_node = c(
      "mean_degree", "mean_strength", "mean_betweenness", "mean_closeness",
      "degree_rmssd", "strength_rmssd", "centrality_instability"
    ),
    hmm_occupancy = c("frac_time"),
    hmm_dwell = c("mean_dwell_hours", "median_dwell_hours", "max_dwell_hours", "n_bouts"),
    hmm_transition = c("TransitionProbability", "Transitions_per_hour"),
    gamm_trajectory = c("auc_per_hour", "dynamic_range", "time_to_peak", "rmssd_pred", "acf1_pred"),
    nonlinear_systems = c(
      "Movement_mean", "Entropy_mean", "Proximity_mean",
      "SocialWithdrawal", "BehavioralInstability"
    ),
    adaptation_kinetics = c(
      "early_late_shift", "volatility_decay", "recovery_slope_per_hour",
      "stabilization_time_hours", "adaptation_half_life_hours", "control_convergence"
    ),
    distance_to_control = c(
      "mean_distance_to_control", "early_distance_to_control",
      "late_distance_to_control", "control_convergence"
    ),
    sleep_like_inactivity = c(
      "inactivity_fraction", "zero_like_inactivity_fraction",
      "active_inactive_transition_rate", "inactive_to_active_rate",
      "active_to_inactive_rate", "mean_inactivity_bout_min",
      "max_inactivity_bout_min", "prolonged_inactivity_fraction",
      "inactivity_fragmentation", "inactivity_bout_count_per_hour",
      "prolonged_inactivity_episodes_per_hour"
    ),
    phase_contrast = c(
      "active_minus_inactive_mean", "active_inactive_ratio_mean",
      "active_minus_inactive_rmssd", "active_minus_inactive_acf1",
      "phase_contrast_strength"
    ),
    phase_timing = c("onset_slope", "offset_slope", "activity_center_of_mass"),
    phase_fragmentation = c(
      "active_fraction", "phase_switch_rate", "inactivity_fragmentation",
      "mean_movement_rmssd"
    ),
    phase_recovery = c(
      "mean_phase_contrast", "mean_deviation_from_control",
      "phase_stabilization_slope", "mean_phase_profile_recovery"
    ),
    phase_predictability = c("phase_predictability_index", "phase_regularization_index")
  )

  source_specs <- source_specs %>%
    mutate(
      CuratedKey = c(
        "raw_multiscale",
        "temporal_instability",
        "state_space_diversity",
        "state_space_switching",
        "early_prediction",
        "social_animal",
        "social_node",
        "hmm_occupancy",
        "hmm_dwell",
        "hmm_transition",
        "gamm_trajectory",
        "nonlinear_systems",
        "adaptation_kinetics",
        "distance_to_control",
        "sleep_like_inactivity",
        "phase_contrast",
        "phase_timing",
        "phase_fragmentation",
        "phase_recovery",
        "phase_predictability"
      ),
      DescriptorKeep = list(
        character(),
        "Metric",
        character(),
        character(),
        character(),
        character(),
        character(),
        "State",
        "State",
        c("State", "NextState"),
        "Metric",
        character(),
        c("Metric", "PhaseClass"),
        c("Metric", "PhaseClass"),
        "PhaseClass",
        "Metric",
        c("Metric", "PhaseClass"),
        "PhaseClass",
        "Metric",
        "Metric"
      )
    )

  loaded <- pmap(
    source_specs,
    function(Source, Domain, Scale, Path, CuratedKey, DescriptorKeep) {
      load_curated_behavior_table(
        path = Path,
        source_label = Source,
        domain_label = Domain,
        scale_label = Scale,
        numeric_keep = curated_numeric[[CuratedKey]],
        descriptor_keep = DescriptorKeep
      )
    }
  )

  feature_long <- map_dfr(loaded, "features") %>%
    bind_rows(map_dfr(optional_behavior_bin_levels, load_hmm_summary_features))

  inventory <- map_dfr(loaded, "inventory") %>%
    bind_rows(tibble(
      Source = "hmm_states",
      Domain = "latent_state",
      Scale = paste(optional_behavior_bin_levels, collapse = ";"),
      Path = "hmm_state_occupancy/dwell/transition summary",
      Status = if_else(any(str_detect(feature_long$feature %||% character(), "^hmm__")), "loaded", "missing_or_empty"),
      n_rows = NA_integer_,
      n_animals = n_distinct(feature_long$AnimalNum[str_detect(feature_long$feature %||% character(), "^hmm__")]),
      n_numeric_features = NA_integer_,
      n_exported_features = sum(str_detect(unique(feature_long$feature %||% character()), "^hmm__"))
    ))

  behavior_matrix <- collapse_behavior_feature_rows(feature_long)

  if (nrow(behavior_matrix) == 0) {
    stop(
      "No behavior features could be loaded for proteomics integration. Checked source specs under: ",
      analysis_ready_dir,
      call. = FALSE
    )
  }

  list(
    matrix = behavior_matrix,
    feature_long = feature_long,
    feature_long_all = feature_long,
    inventory = inventory
  )
}

make_axis <- function(dat, regex_pattern, axis_name) {
  cols <- names(dat)[
    stringr::str_detect(
      names(dat),
      stringr::regex(regex_pattern, ignore_case = TRUE)
    )
  ]

  cols <- cols[vapply(dat[cols], is.numeric, logical(1))]

  if (length(cols) == 0) {
    return(tibble(
      AnimalNum = dat$AnimalNum,
      !!axis_name := NA_real_,
      n_features = 0L
    ))
  }

  tibble(
    AnimalNum = dat$AnimalNum,
    !!axis_name := rowMeans(
      as.data.frame(lapply(dat[cols], z_within_metric)),
      na.rm = TRUE
    ),
    n_features = length(cols)
  )
}

safe_slice_max_abs <- function(df, value_col, n = 80) {
  if (nrow(df) == 0) return(df)
  df %>%
    mutate(.abs_value = abs(.data[[value_col]])) %>%
    slice_max(.abs_value, n = min(n, nrow(.)), with_ties = FALSE) %>%
    select(-.abs_value)
}

short_label <- function(x) {
  x %>%
    stringr::str_replace_all("m_neuron_neuropil__", "") %>%
    stringr::str_replace_all("Neuropil_", "") %>%
    stringr::str_replace_all("RNP_RNA_processing_full", "RNP/RNA full") %>%
    stringr::str_replace_all("RNP_RNA_processing_main", "RNP/RNA main") %>%
    stringr::str_replace_all("mito_bioenergetics", "mito") %>%
    stringr::str_replace_all("ribosome_translation", "ribosome") %>%
    stringr::str_replace_all("chromatin_RNP_related_exploratory", "chrom/RNP") %>%
    stringr::str_replace_all("synaptic_cytoskeleton", "syn/cyto") %>%
    stringr::str_replace_all("__", " | ")
}

pretty_axis_label <- function(x) {
  x %>%
    str_replace_all("_axis$", "") %>%
    str_replace_all("_", " ") %>%
    str_to_sentence()
}

pretty_behavior_label <- function(x) {
  x %>%
    str_replace_all("^[^|]+ \\| [^|]+ \\| [^|]+ \\| ", "") %>%
    str_replace_all(" \\| animal_mean \\| animal_level$", "") %>%
    str_replace_all(" \\| animal_level \\| summary$", "") %>%
    str_replace_all(" \\| animal_mean \\| ", " | ") %>%
    str_replace_all("_", " ") %>%
    str_squish() %>%
    str_to_sentence()
}

safe_cor_p_local <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE)$p.value)
}

cor_ci_approx <- function(r, n, conf = 0.95) {
  if (!is.finite(r) || !is.finite(n) || n < 4 || abs(r) >= 1) {
    return(c(low = NA_real_, high = NA_real_))
  }
  z <- atanh(r)
  se <- 1 / sqrt(n - 3)
  crit <- qnorm(1 - (1 - conf) / 2)
  tanh(c(low = z - crit * se, high = z + crit * se))
}

# ------------------------------------------------
# 3) Axis definitions
# ------------------------------------------------

behavior_axis_specs <- tibble(
  Axis = c(
    "locomotor_activity_axis",
    "temporal_instability_axis",
    "latent_state_flexibility_axis",
    "social_engagement_axis",
    "social_withdrawal_isolation_axis",
    "adaptation_recovery_axis",
    "phase_organization_axis",
    "quiescence_inactivity_axis"
  ),
  Regex = c(
    "movement.*mean|movement.*per_hour|movementdistance.*per_hour|locomotor|activity",
    "rmssd|acf1|cv|fano|instability|volatility|entropy|switch_rate|transition_entropy|persistence",
    "state_occupancy|state_switch|dwell|latent|hmm|flexibility|rigidity|transition_probability",
    "proximity|social_engagement|contact_fraction|mean_contact|degree|strength|partner_entropy|partner_evenness|social_state_fraction",
    "socialwithdrawal|social_withdrawal|passive_isolation|active_isolation|fragmentation|inactive_state_fraction|isolation",
    "adaptation|recovery|stabilization|half_life|early_late_shift|control_convergence|volatility_decay|slope|trajectory",
    "phase|active_minus_inactive|active_inactive|phase_contrast|phase_predictability|phase_regularization|timing",
    "inactivity|quiescence|rest_like|prolonged|bout|zero_like|sleep_like"
  )
)

proteomic_axis_specs <- tibble(
  Axis = c(
    "RNA_RNP_splicing_axis",
    "translation_ribosome_axis",
    "mitochondrial_OXPHOS_axis",
    "proteostasis_endolysosomal_axis",
    "synaptic_plasticity_axis"
  ),
  Regex = c(
    "rna|rnp|splice|splicing",
    "translation|ribosome|rpl|rps|eif",
    "mitochond|oxphos|respirat|electron",
    "proteostasis|lysosom|endosom|ubiquitin|proteasom",
    "synap|plastic|vesicle|glutamate|gaba"
  )
)

# ------------------------------------------------
# 4) Core analysis function
# ------------------------------------------------

run_crossmodal_analysis <- function(
    beh,
    prot,
    subset_label,
    output_dir,
    proteomics_label,
    behavior_axis_specs,
    proteomic_axis_specs
) {

  subset_dir <- file.path(output_dir, subset_label)
  ensure_dir(subset_dir)
  ensure_dir(file.path(subset_dir, "tables"))
  ensure_dir(file.path(subset_dir, "figures"))

  message("\n--- Running subset: ", subset_label, " ---")

  merged <- beh %>%
    inner_join(prot, by = "AnimalNum")

  matched_ids <- sort(unique(merged$AnimalNum))

  cat("Matched animals in ", subset_label, ": ", length(matched_ids), "\n", sep = "")
  print(matched_ids)

  if (length(matched_ids) < 3) {
    warning(
      "Skipping subset ",
      subset_label,
      ": fewer than 3 matched animals."
    )
    return(NULL)
  }

  write_table(
    merged,
    file.path(subset_dir, "tables", paste0("behavior_proteomics_merged_", subset_label, ".csv"))
  )

  behavioral_latent_axes <- reduce(
    purrr::map2(
      behavior_axis_specs$Regex,
      behavior_axis_specs$Axis,
      ~make_axis(beh, .x, .y) %>% select(-n_features)
    ),
    full_join,
    by = "AnimalNum"
  )

  proteomic_latent_axes <- reduce(
    purrr::map2(
      proteomic_axis_specs$Regex,
      proteomic_axis_specs$Axis,
      ~make_axis(prot, .x, .y) %>% select(-n_features)
    ),
    full_join,
    by = "AnimalNum"
  )

  axis_feature_inventory <- bind_rows(
    map2_dfr(
      behavior_axis_specs$Regex,
      behavior_axis_specs$Axis,
      ~tibble(
        Modality = "behavior",
        Axis = .y,
        Feature = names(beh)[
          str_detect(names(beh), regex(.x, ignore_case = TRUE))
        ],
        Subset = subset_label
      )
    ),
    map2_dfr(
      proteomic_axis_specs$Regex,
      proteomic_axis_specs$Axis,
      ~tibble(
        Modality = "proteomics",
        Axis = .y,
        Feature = names(prot)[
          str_detect(names(prot), regex(.x, ignore_case = TRUE))
        ],
        Subset = subset_label
      )
    )
  )

  axis_merged <- behavioral_latent_axes %>%
    inner_join(proteomic_latent_axes, by = "AnimalNum")

  behavior_axis_cols <- setdiff(names(behavioral_latent_axes), "AnimalNum")
  proteomic_axis_cols <- setdiff(names(proteomic_latent_axes), "AnimalNum")

  latent_axis_associations <- expand_grid(
    BehaviorAxis = behavior_axis_cols,
    ProteomicAxis = proteomic_axis_cols
  ) %>%
    mutate(
      n = map2_int(
        BehaviorAxis,
        ProteomicAxis,
        ~sum(is.finite(axis_merged[[.x]]) & is.finite(axis_merged[[.y]]))
      ),
      spearman_rho = map2_dbl(
        BehaviorAxis,
        ProteomicAxis,
        ~safe_cor(axis_merged[[.x]], axis_merged[[.y]], "spearman")
      ),
      pearson_r = map2_dbl(
        BehaviorAxis,
        ProteomicAxis,
        ~safe_cor(axis_merged[[.x]], axis_merged[[.y]], "pearson")
      ),
      spearman_p = map2_dbl(
        BehaviorAxis,
        ProteomicAxis,
        ~safe_cor_p_local(axis_merged[[.x]], axis_merged[[.y]], "spearman")
      ),
      spearman_ci_low = map2_dbl(spearman_rho, n, ~cor_ci_approx(.x, .y)[["low"]]),
      spearman_ci_high = map2_dbl(spearman_rho, n, ~cor_ci_approx(.x, .y)[["high"]]),
      Subset = subset_label,
      ProteomicsLabel = proteomics_label,
      EvidenceUse = if_else(
        n < 12,
        "exploratory_effect_size_small_n",
        "exploratory_association"
      )
    ) %>%
    mutate(
      spearman_fdr = p.adjust(spearman_p, method = "BH"),
      BehaviorAxisLabel = pretty_axis_label(BehaviorAxis),
      ProteomicAxisLabel = pretty_axis_label(ProteomicAxis)
    ) %>%
    arrange(desc(abs(spearman_rho)))

  curated_behavior_proteomics_models <- latent_axis_associations %>%
    mutate(
      BiologicalModel = paste(BehaviorAxis, "aligned_with", ProteomicAxis),
      ClaimType = "associative",
      AllowedInterpretation = paste(
        "Behavioral axis covaries with a curated hippocampal proteomic module axis.",
        "Proteomics scope:",
        proteomics_label,
        "Subset:",
        subset_label
      ),
      ReviewerRisk = if_else(n < 12, "high_small_n", "medium_exploratory"),
      StableForMainText = FALSE
    )

  write_table(
    behavioral_latent_axes,
    file.path(subset_dir, "tables", paste0("behavioral_latent_axes_", subset_label, ".csv"))
  )

  write_table(
    proteomic_latent_axes,
    file.path(subset_dir, "tables", paste0("proteomic_latent_axes_", subset_label, ".csv"))
  )

  write_table(
    latent_axis_associations,
    file.path(subset_dir, "tables", paste0("latent_axis_associations_", subset_label, ".csv"))
  )

  write_table(
    curated_behavior_proteomics_models,
    file.path(subset_dir, "tables", paste0("curated_behavior_proteomics_models_", subset_label, ".csv"))
  )

  write_table(
    axis_feature_inventory,
    file.path(subset_dir, "tables", paste0("latent_axis_feature_inventory_", subset_label, ".csv"))
  )

  axis_plot_tbl <- latent_axis_associations %>%
    filter(is.finite(spearman_rho)) %>%
    mutate(
      BehaviorAxisLabel = factor(BehaviorAxisLabel, levels = rev(unique(BehaviorAxisLabel))),
      ProteomicAxisLabel = factor(ProteomicAxisLabel, levels = unique(ProteomicAxisLabel)),
      EvidenceLabel = case_when(
        !is.na(spearman_fdr) & spearman_fdr < 0.10 ~ "*",
        TRUE ~ ""
      )
    )

  if (nrow(axis_plot_tbl) > 0) {
    p_axis_heat <- axis_plot_tbl %>%
      ggplot(aes(ProteomicAxisLabel, BehaviorAxisLabel, fill = spearman_rho)) +
      geom_tile(color = "white", linewidth = 0.35) +
      geom_text(aes(label = EvidenceLabel), size = 2.4, color = "black") +
      scale_fill_gradient2(
        low = "#457B9D",
        mid = "white",
        high = "#E63946",
        midpoint = 0,
        limits = c(-1, 1),
        name = "Spearman\nrho"
      ) +
      labs(
        title = paste0("Curated behavior-proteomics axes: ", subset_label),
        subtitle = paste0(
          "n = ",
          n_distinct(axis_merged$AnimalNum),
          "; asterisk marks BH FDR < 0.10 across axis pairs"
        ),
        x = "Proteomic axis",
        y = "Behavioral axis"
      ) +
      theme_classic(base_size = 7) +
      theme(
        plot.title = element_text(size = 8.5, face = "bold"),
        plot.subtitle = element_text(size = 6.5),
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_text(size = 7),
        legend.position = "right"
      )

    save_plot_svg_pdf(
      p_axis_heat,
      file.path(subset_dir, "figures", paste0("behavior_proteomics_axis_heatmap_", subset_label)),
      width = 115,
      height = 75
    )
  }

  top_axis <- latent_axis_associations %>%
    filter(is.finite(spearman_rho), n >= 4) %>%
    slice_max(abs(spearman_rho), n = 1, with_ties = FALSE)

  if (nrow(top_axis) == 1) {
    top_axis_plot_tbl <- axis_merged %>%
      select(
        AnimalNum,
        BehaviorScore = all_of(top_axis$BehaviorAxis[[1]]),
        ProteomicScore = all_of(top_axis$ProteomicAxis[[1]])
      ) %>%
      left_join(
        merged %>%
          select(AnimalNum, any_of(c("Group", "StressGroup"))) %>%
          distinct(AnimalNum, .keep_all = TRUE),
        by = "AnimalNum"
      ) %>%
      {
        if (!"Group" %in% names(.)) .$Group <- NA_character_
        if (!"StressGroup" %in% names(.)) .$StressGroup <- NA_character_
        .
      } %>%
      mutate(
        PlotGroup = coalesce(as.character(.data[["Group"]]), as.character(.data[["StressGroup"]]), "All")
      )

    p_top_axis <- top_axis_plot_tbl %>%
      ggplot(aes(BehaviorScore, ProteomicScore, fill = PlotGroup)) +
      geom_hline(yintercept = 0, linewidth = 0.2, color = "grey80") +
      geom_vline(xintercept = 0, linewidth = 0.2, color = "grey80") +
      geom_point(shape = 21, size = 2.4, alpha = 0.9, color = "black", stroke = 0.25) +
      geom_smooth(aes(group = 1), method = "lm", se = FALSE, linewidth = 0.4, color = "black") +
      scale_fill_manual(values = mmm_group_colors, drop = FALSE) +
      labs(
        title = "Strongest curated axis association",
        subtitle = paste0(
          pretty_axis_label(top_axis$BehaviorAxis[[1]]),
          " vs ",
          pretty_axis_label(top_axis$ProteomicAxis[[1]]),
          "; n = ",
          top_axis$n[[1]],
          "; rho = ",
          round(top_axis$spearman_rho[[1]], 2),
          "; FDR = ",
          if_else(is.na(top_axis$spearman_fdr[[1]]), "NA", formatC(top_axis$spearman_fdr[[1]], format = "f", digits = 3))
        ),
        x = pretty_axis_label(top_axis$BehaviorAxis[[1]]),
        y = pretty_axis_label(top_axis$ProteomicAxis[[1]])
      ) +
      theme_classic(base_size = 7) +
      theme(
        plot.title = element_text(size = 8.5, face = "bold"),
        plot.subtitle = element_text(size = 6.5),
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 6),
        legend.position = "top",
        legend.title = element_blank()
      )

    save_plot_svg_pdf(
      p_top_axis,
      file.path(subset_dir, "figures", paste0("strongest_axis_behavior_proteomics_relationship_", subset_label)),
      width = 90,
      height = 72
    )
  }

  beh_num <- merged %>%
    select(where(is.numeric))

  if (ncol(beh_num) < 2) {
    warning("Skipping correlation matrix for ", subset_label, ": fewer than two numeric columns.")
    return(list(latent = latent_axis_associations, cross = tibble()))
  }

  cor_mat <- suppressWarnings(
    cor(beh_num, use = "pairwise.complete.obs", method = "spearman")
  )

  cor_tbl <- as.data.frame(as.table(cor_mat)) %>%
    rename(
      Feature1 = Var1,
      Feature2 = Var2,
      SpearmanRho = Freq
    ) %>%
    filter(Feature1 != Feature2) %>%
    mutate(
      Subset = subset_label,
      ProteomicsLabel = proteomics_label
    ) %>%
    arrange(desc(abs(SpearmanRho)))

  write_table(
    cor_tbl,
    file.path(subset_dir, "tables", paste0("behavior_proteomics_correlations_", subset_label, ".csv"))
  )

  beh_vars <- names(beh %>% select(where(is.numeric)))
  prot_vars <- setdiff(names(prot %>% select(where(is.numeric))), beh_vars)

  cross_tbl <- cor_tbl %>%
    filter(
      (Feature1 %in% beh_vars & Feature2 %in% prot_vars) |
        (Feature2 %in% beh_vars & Feature1 %in% prot_vars)
    ) %>%
    arrange(desc(abs(SpearmanRho)))

  write_table(
    cross_tbl,
    file.path(subset_dir, "tables", paste0("crossmodal_behavior_proteomics_correlations_", subset_label, ".csv"))
  )

  if (nrow(cross_tbl) > 0) {

    plot_tbl <- cross_tbl %>%
      filter(abs(SpearmanRho) >= 0.7) %>%
      mutate(
        BehaviorFeature = if_else(Feature1 %in% beh_vars, Feature1, Feature2),
        ProteomicsFeature = if_else(Feature1 %in% prot_vars, Feature1, Feature2),
        BehaviorFeature = pretty_behavior_label(short_label(BehaviorFeature)),
        ProteomicsFeature = short_label(ProteomicsFeature)
      ) %>%
      distinct(BehaviorFeature, ProteomicsFeature, .keep_all = TRUE) %>%
      filter(is.finite(SpearmanRho)) %>%
      safe_slice_max_abs("SpearmanRho", n = 45)

    p_heat <- plot_tbl %>%
      ggplot(aes(ProteomicsFeature, BehaviorFeature, fill = SpearmanRho)) +
      geom_tile(color = "white", linewidth = 0.25) +
      scale_fill_gradient2(
        low = "#457B9D",
        mid = "white",
        high = "#E63946",
        midpoint = 0,
        limits = c(-1, 1),
        name = "Spearman\nrho"
      ) +
      labs(
        title = paste0(
          "Behavior–proteomics correlations: ",
          subset_label,
          " neuropil (n = ",
          n_distinct(merged$AnimalNum),
          ")"
        ),
        subtitle = paste0("Proteomics: ", proteomics_label, "; exploratory cross-modal correlations"),
        x = "Proteomic module",
        y = "Behavioral feature"
      ) +
      theme_classic(base_size = 7) +
      theme(
        plot.title = element_text(size = 8.5, face = "bold"),
        plot.subtitle = element_text(size = 6.5),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5.5),
        axis.text.y = element_text(size = 5.5),
        axis.title = element_text(size = 7),
        legend.position = "right"
      )

    save_plot_svg_pdf(
      p_heat,
      file.path(subset_dir, "figures", paste0("behavior_proteomics_heatmap_", subset_label)),
      width = 180,
      height = 140
    )

  } else {
    warning("No crossmodal correlations found for subset: ", subset_label)
  }

  if (requireNamespace("mixOmics", quietly = TRUE) && n_distinct(merged$AnimalNum) >= 12) {

    x <- merged %>%
      select(any_of(beh_vars)) %>%
      select(where(is.numeric)) %>%
      as.matrix()

    y <- merged %>%
      select(any_of(prot_vars)) %>%
      select(where(is.numeric)) %>%
      as.matrix()

    x <- x[, colSums(is.finite(x)) > 2, drop = FALSE]
    y <- y[, colSums(is.finite(y)) > 2, drop = FALSE]

    if (ncol(x) > 1 && ncol(y) > 1 && nrow(x) >= 3) {

      pls_fit <- tryCatch(
        mixOmics::pls(x, y, ncomp = 2),
        error = function(e) {
          warning("PLS failed for ", subset_label, ": ", conditionMessage(e))
          NULL
        }
      )

      if (!is.null(pls_fit)) {
        latent_tbl <- tibble(
          AnimalNum = merged$AnimalNum,
          comp1_behavior = pls_fit$variates$X[, 1],
          comp1_proteomics = pls_fit$variates$Y[, 1],
          Subset = subset_label,
          ProteomicsLabel = proteomics_label
        )

        write_table(
          latent_tbl,
          file.path(subset_dir, "tables", paste0("pls_latent_scores_", subset_label, ".csv"))
        )

        rho_pls <- suppressWarnings(cor(
          latent_tbl$comp1_behavior,
          latent_tbl$comp1_proteomics,
          method = "spearman",
          use = "complete.obs"
        ))

        latent_tbl_plot <- latent_tbl %>%
          left_join(
            merged %>%
              select(AnimalNum, any_of(c("Group", "StressGroup"))) %>%
              distinct(AnimalNum, .keep_all = TRUE),
            by = "AnimalNum"
          ) %>%
          {
            if (!"Group" %in% names(.)) .$Group <- NA_character_
            if (!"StressGroup" %in% names(.)) .$StressGroup <- NA_character_
            .
          } %>%
          mutate(
            PlotGroup = coalesce(
              as.character(.data[["Group"]]),
              as.character(.data[["StressGroup"]]),
              "All"
            )
          )

        p_pls <- latent_tbl_plot %>%
          ggplot(aes(comp1_behavior, comp1_proteomics, fill = PlotGroup)) +
          geom_point(
            shape = 21,
            size = 2.2,
            alpha = 0.9,
            color = "black",
            stroke = 0.25
          ) +
          geom_smooth(
            aes(group = 1),
            method = "lm",
            se = TRUE,
            linewidth = 0.4,
            color = "black",
            fill = "grey80"
          ) +
          labs(
            title = "Exploratory latent behavior–proteomics alignment",
            subtitle = paste0(
              subset_label,
              " neuropil; n = ",
              n_distinct(latent_tbl$AnimalNum),
              "; Spearman rho = ",
              round(rho_pls, 2)
            ),
            x = "Behavior latent component 1",
            y = "Proteomics latent component 1"
          ) +
          scale_fill_manual(values = mmm_group_colors, drop = FALSE) +
          theme_classic(base_size = 7) +
          theme(
            plot.title = element_text(size = 8.5, face = "bold"),
            plot.subtitle = element_text(size = 6.5),
            axis.title = element_text(size = 7),
            axis.text = element_text(size = 6),
            legend.position = "top",
            legend.title = element_blank()
          )

        save_plot_svg_pdf(
          p_pls,
          file.path(subset_dir, "figures", paste0("latent_behavior_proteomics_relationship_", subset_label)),
          width = 85,
          height = 75
        )
      }
    } else {
      warning("Skipping PLS for ", subset_label, ": insufficient numeric behavior/proteomics features.")
    }
  } else {
    message("Skipping high-dimensional PLS for subset: ", subset_label, " (requires mixOmics and n >= 12).")
  }

  message("Completed subset: ", subset_label)

  list(
    latent = latent_axis_associations,
    cross = cross_tbl
  )
}

# ------------------------------------------------
# 5) Load behavior once
# ------------------------------------------------

behavior_features <- build_behavior_feature_matrix()
beh_raw <- behavior_features$matrix

message(
  "Loaded integrated behavior feature matrix: ",
  nrow(beh_raw),
  " animals x ",
  ncol(beh_raw) - 3,
  " behavior features."
)

# ------------------------------------------------
# 6) Run once per proteomics file
# ------------------------------------------------

for (proteomics_file in proteomics_files) {

  proteomics_file_base <- tools::file_path_sans_ext(basename(proteomics_file))

  proteomics_label <- stringr::str_remove(
    proteomics_file_base,
    "^module_scores_"
  )

  output_dir <- file.path(
    base_output_dir,
    paste0("proteomics_integration_", proteomics_label)
  )

  message("\n========================================")
  message("Running behavior-proteomics integration")
  message("Proteomics file: ", proteomics_file)
  message("Proteomics label: ", proteomics_label)
  message("Output dir: ", output_dir)
  message("========================================\n")

  ensure_dir(output_dir)
  ensure_dir(file.path(output_dir, "tables"))
  ensure_dir(file.path(output_dir, "figures"))
  analysis_output_dirs(output_dir)

  write_output_manifest(
    output_dir,
    script_name = "12_behavior_proteomics_integration.R",
    analysis_name = paste(
      "sex-aware behavior-proteomics systems alignment:",
      proteomics_label
    ),
    primary_tables = c(
      "tables/integrated_behavior_feature_matrix.csv",
      "tables/integrated_behavior_feature_sources.csv",
      "tables/all_latent_axis_associations.csv",
      "tables/all_crossmodal_behavior_proteomics_correlations.csv",
      "tables/all_curated_behavior_proteomics_models.csv"
    ),
    primary_figures = c(
      "sex-specific subfolders/figures/behavior_proteomics_axis_heatmap_<subset>.svg",
      "sex-specific subfolders/figures/strongest_axis_behavior_proteomics_relationship_<subset>.svg",
      "sex-specific subfolders/figures/behavior_proteomics_heatmap_<subset>.svg",
      "sex-specific subfolders/figures/latent_behavior_proteomics_relationship_<subset>.svg"
    ),
    notes = c(
      paste("Proteomics input:", basename(proteomics_file)),
      paste("Proteomics label:", proteomics_label),
      "Behavior input is an integrated animal-level matrix built from scripts 03, 06, 07, 08, 09, 10, 11, 13, 15, 16 and 17 when available.",
      "If Sex is present, analyses are run by sex-specific matched subsets.",
      "Pooled analyses are only run if no usable Sex column exists.",
      "Use curated low-dimensional axes; feature-explosion correlations are retained only as exploratory supplements."
    )
  )

  write_table(
    beh_raw,
    file.path(output_dir, "tables", "integrated_behavior_feature_matrix.csv")
  )

  write_table(
    behavior_features$inventory,
    file.path(output_dir, "tables", "integrated_behavior_feature_sources.csv")
  )

  write_table(
    behavior_features$feature_long,
    file.path(output_dir, "tables", "integrated_behavior_features_long.csv")
  )

  if (!file.exists(proteomics_file)) {
    warning("Skipping missing proteomics file: ", proteomics_file)
    next
  }

  beh <- beh_raw
  prot <- read_behavior_table(proteomics_file)

  beh_animal_col <- first_existing_col(
    beh,
    c("AnimalNum", "AnimalID", "Animal", "MouseID", "Mouse", "ID"),
    TRUE,
    "behavior animal column"
  )

  prot_animal_col <- first_existing_col(
    prot,
    c("AnimalNum", "AnimalID", "Animal", "MouseID", "Mouse", "ID"),
    TRUE,
    "proteomics animal column"
  )

  beh <- beh %>%
    rename(AnimalNum = all_of(beh_animal_col)) %>%
    mutate(AnimalNum = normalize_animal_id(AnimalNum))

  prot <- prot %>%
    rename(AnimalNum = all_of(prot_animal_col)) %>%
    mutate(AnimalNum = normalize_animal_id(AnimalNum))

  matched_ids <- sort(intersect(beh$AnimalNum, prot$AnimalNum))

  cat("Behavior animals:", n_distinct(beh$AnimalNum), "\n")
  cat("Proteomics animals:", n_distinct(prot$AnimalNum), "\n")
  cat("Matched animals:", length(matched_ids), "\n")

  cat("\nMatched IDs:\n")
  print(matched_ids)

  cat("\nBehavior-only IDs:\n")
  print(sort(setdiff(unique(beh$AnimalNum), unique(prot$AnimalNum))))

  cat("\nProteomics-only IDs:\n")
  print(sort(setdiff(unique(prot$AnimalNum), unique(beh$AnimalNum))))

  if (length(matched_ids) < 3) {
    warning(
      "Skipping ",
      proteomics_label,
      ": too few matched animals after ID normalization: ",
      length(matched_ids)
    )
    next
  }

  # ------------------------------------------------
  # Sex-aware split
  # ------------------------------------------------

  beh_sex_col <- intersect(c("Sex", "sex"), names(beh))
  prot_sex_col <- intersect(c("Sex", "sex"), names(prot))

  if (length(beh_sex_col) > 0 && length(prot_sex_col) > 0) {

    beh <- beh %>%
      mutate(Sex = standardize_sex(.data[[beh_sex_col[1]]]))

    prot <- prot %>%
      mutate(Sex = standardize_sex(.data[[prot_sex_col[1]]]))

    common_sexes <- intersect(
      sort(unique(beh$Sex[!is.na(beh$Sex) & beh$Sex != ""])),
      sort(unique(prot$Sex[!is.na(prot$Sex) & prot$Sex != ""]))
    )

    message("Detected common sexes: ", paste(common_sexes, collapse = ", "))

    if (length(common_sexes) == 0) {
      warning("No common sex labels detected. Running pooled analysis instead.")

      results <- list(
        pooled = run_crossmodal_analysis(
          beh = beh,
          prot = prot,
          subset_label = "pooled",
          output_dir = output_dir,
          proteomics_label = proteomics_label,
          behavior_axis_specs = behavior_axis_specs,
          proteomic_axis_specs = proteomic_axis_specs
        )
      )

    } else {

      results <- purrr::map(
        common_sexes,
        function(sex_now) {

          beh_sub <- beh %>%
            filter(Sex == sex_now)

          prot_sub <- prot %>%
            filter(Sex == sex_now)

          run_crossmodal_analysis(
            beh = beh_sub,
            prot = prot_sub,
            subset_label = sex_now,
            output_dir = output_dir,
            proteomics_label = proteomics_label,
            behavior_axis_specs = behavior_axis_specs,
            proteomic_axis_specs = proteomic_axis_specs
          )
        }
      )

      names(results) <- common_sexes
    }

  } else {

    message("Sex column not detected in both behavior and proteomics. Running pooled analysis.")

    results <- list(
      pooled = run_crossmodal_analysis(
        beh = beh,
        prot = prot,
        subset_label = "pooled",
        output_dir = output_dir,
        proteomics_label = proteomics_label,
        behavior_axis_specs = behavior_axis_specs,
        proteomic_axis_specs = proteomic_axis_specs
      )
    )
  }

  # ------------------------------------------------
  # Combine subset-level results
  # ------------------------------------------------

  valid_results <- compact(results)

  all_latent <- purrr::map_dfr(valid_results, "latent")
  all_cross <- purrr::map_dfr(valid_results, "cross")

  if (nrow(all_latent) > 0) {
    write_table(
      all_latent,
      file.path(output_dir, "tables", "all_latent_axis_associations.csv")
    )

    all_curated <- all_latent %>%
      mutate(
        BiologicalModel = paste(BehaviorAxis, "aligned_with", ProteomicAxis),
        ClaimType = "associative",
        AllowedInterpretation = paste(
          "Behavioral axis covaries with a curated hippocampal proteomic module axis.",
          "Proteomics scope:",
          ProteomicsLabel,
          "Subset:",
          Subset
        ),
        ReviewerRisk = if_else(n < 12, "high_small_n", "medium_exploratory"),
        StableForMainText = FALSE
      )

    write_table(
      all_curated,
      file.path(output_dir, "tables", "all_curated_behavior_proteomics_models.csv")
    )
  }

  if (nrow(all_cross) > 0) {
    write_table(
      all_cross,
      file.path(output_dir, "tables", "all_crossmodal_behavior_proteomics_correlations.csv")
    )
  }

  message("Behavior-proteomics integration complete for: ", proteomics_label)
}

message("All behavior-proteomics integrations complete.")
