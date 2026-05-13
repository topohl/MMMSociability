# ================================================================
# Temporal Instability Analysis
# MMMSociability
# ================================================================
# Goal:
#   Quantify whether stress changes the temporal structure of
#   behavior beyond mean-level effects.
#
# Primary readout:
#   Movement. RMSSD/ACF1 have the cleanest interpretation for
#   locomotor temporal instability.
#
# Secondary exploratory readouts:
#   Entropy and Proximity. These are retained as robustness /
#   complementary behavioral-dynamics outputs, but should not be
#   interpreted as the primary temporal-instability phenotype.
#
# Input expectation:
#   Run Analysis/03_build_multiscale_behavior_metrics.R first.
#
# Recommended scale:
#   1-5 min bins. Phase-level data are too coarse for RMSSD/ACF metrics.
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(readr)
  library(zoo)
})

source("C:/Users/topohl/Documents/GitHub/MMMSociability/Functions/behavioral_dynamics_helpers.R")
source("C:/Users/topohl/Documents/GitHub/MMMSociability/Functions/behavioral_dynamics_stats_helpers.R")

# ------------------------------------------------
# USER INPUT
# ------------------------------------------------

bin_level <- "5min_based"
input_file <- file.path("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/analysis_ready/03_derived_metrics", bin_level, "all_behavior_metrics.csv")
output_dir <- file.path("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/analysis_ready/06_behavioral_dynamics/temporal_instability", bin_level)

# Movement is the primary metric for temporal instability.
# Entropy and proximity are retained as secondary exploratory readouts.
primary_metric <- "Movement"
secondary_metrics <- c("Entropy", "Proximity")
metrics_to_analyze <- c(primary_metric, secondary_metrics)

# Prefer normalized proximity for temporal instability because raw contact
# seconds scale with bin size. Falls back to Proximity if ProximityFraction is absent.
proximity_col <- "ProximityFraction"

group_colors <- c(
  "CON" = "#3d3b6e",
  "SUS" = "#e63947",
  "RES" = "#C6C3BB"
)
group_levels <- c("CON", "RES", "SUS")
pairwise_contrasts <- c("RES-CON", "SUS-CON", "SUS-RES")
instability_outcomes <- c("mean", "cv", "fano", "rmssd", "acf1")
metric_labels <- c(
  mean = "Mean",
  cv = "CV",
  fano = "Fano",
  rmssd = "RMSSD",
  acf1 = "ACF1",
  delta_mean = "Delta mean",
  delta_cv = "Delta CV",
  delta_fano = "Delta Fano",
  delta_rmssd = "Delta RMSSD",
  delta_acf1 = "Delta ACF1"
)

parse_cage_change_index <- function(x) {
  idx <- suppressWarnings(as.integer(stringr::str_extract(as.character(x), "\\d+")))
  ifelse(is.na(idx), dense_rank(as.character(x)) - 1L, idx)
}

add_temporal_context <- function(dat) {
  dat %>%
    mutate(
      Group = factor(as.character(Group), levels = group_levels),
      CageChangeIndex = parse_cage_change_index(CageChange)
    ) %>%
    group_by(AnimalNum, CageChange) %>%
    arrange(TimeIndex, .by_group = TRUE) %>%
    mutate(
      PhaseBlock = cumsum(as.character(Phase) != lag(as.character(Phase), default = first(as.character(Phase)))) + 1L
    ) %>%
    group_by(AnimalNum, CageChange, PhaseBlock) %>%
    mutate(
      TimeFromPhaseStart = TimeIndex - min(TimeIndex, na.rm = TRUE)
    ) %>%
    ungroup()
}

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

add_stat_columns <- function(tbl) {
  if (nrow(tbl) == 0) return(tbl)
  tbl %>%
    mutate(
      is_pairwise_contrast = contrast %in% pairwise_contrasts,
      p.sign = sig_from_p(p.adjust_bh_family),
      Significant = p.adjust_bh_family < 0.05,
      ReportingP = p.adjust_bh_family,
      ReportingSignificance = p.sign,
      ReportingCorrection = "BH across all pairwise group contrasts"
    ) %>%
    select(-any_of(c("is_primary_contrast", "p.adjust_bh_primary", "p.sign_primary")))
}

write_stats_package <- function(dat, analysis_name, value_cols, by_cols, summary_group_cols) {
  summary_tbl <- make_dynamics_group_summary(dat, value_cols = value_cols, group_cols = summary_group_cols)
  contrast_tbl <- make_dynamics_group_contrasts(dat, value_cols = value_cols, by_cols = by_cols) %>%
    add_stat_columns()

  write_table(summary_tbl, file.path(output_dir, "stats_tables", paste0(analysis_name, "_group_summary.csv")))
  write_table(contrast_tbl, file.path(output_dir, "stats_tables", paste0(analysis_name, "_group_contrasts.csv")))

  invisible(list(summary = summary_tbl, contrasts = contrast_tbl))
}

format_p_label <- function(p) {
  case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "p<.001",
    TRUE ~ paste0("p=", sub("^0", "", formatC(p, format = "f", digits = 3)))
  )
}

make_sig_labels <- function(data, contrast_tbl, outcome, facet_cols, value_col = NULL, pairwise_only = TRUE) {
  if (is.null(value_col)) value_col <- outcome
  empty_labels <- tibble(label = character(), y_pos = numeric(), x_pos = numeric())
  if (length(facet_cols) > 0) {
    for (nm in facet_cols) empty_labels[[nm]] <- data[[nm]][0]
  }
  if (nrow(contrast_tbl) == 0 || !value_col %in% names(data)) return(empty_labels)

  label_tbl <- contrast_tbl %>%
    filter(
      Outcome == outcome,
      if (pairwise_only) contrast %in% pairwise_contrasts else TRUE,
      !is.na(ReportingP),
      ReportingP < 0.10
    ) %>%
    mutate(label = paste0(contrast, " ", format_p_label(ReportingP))) %>%
    group_by(across(all_of(facet_cols))) %>%
    summarise(label = paste(label, collapse = "\n"), .groups = "drop")

  if (nrow(label_tbl) == 0) return(empty_labels)

  y_tbl <- data %>%
    group_by(across(all_of(facet_cols))) %>%
    summarise(
      y_pos = max(.data[[value_col]], na.rm = TRUE),
      y_min = min(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y_range = if_else(is.finite(y_pos - y_min) & y_pos > y_min, y_pos - y_min, abs(y_pos)),
      y_pos = y_pos + if_else(is.finite(y_range) & y_range > 0, 0.12 * y_range, 0.1)
    )

  label_tbl %>%
    left_join(y_tbl, by = facet_cols) %>%
    mutate(x_pos = 2)
}

make_pairwise_brackets <- function(data,
                                   contrast_tbl,
                                   outcome,
                                   facet_cols,
                                   value_col = NULL,
                                   p_threshold = 0.05,
                                   include_trends = FALSE) {
  if (is.null(value_col)) value_col <- outcome
  empty_brackets <- tibble(
    contrast = character(), x_start = numeric(), x_end = numeric(),
    y = numeric(), y_tip = numeric(), label = character()
  )
  if (length(facet_cols) > 0) {
    for (nm in facet_cols) empty_brackets[[nm]] <- data[[nm]][0]
  }
  if (nrow(contrast_tbl) == 0 || !value_col %in% names(data)) return(empty_brackets)

  threshold <- if (include_trends) 0.10 else p_threshold

  brackets <- contrast_tbl %>%
    filter(
      Outcome == outcome,
      contrast %in% pairwise_contrasts,
      !is.na(ReportingP),
      ReportingP < threshold
    ) %>%
    mutate(
      x_start = match(group_ref, group_levels),
      x_end = match(group_comp, group_levels),
      label = format_p_label(ReportingP)
    ) %>%
    filter(is.finite(x_start), is.finite(x_end))

  if (nrow(brackets) == 0) return(empty_brackets)

  y_tbl <- data %>%
    group_by(across(all_of(facet_cols))) %>%
    summarise(
      y_max = max(.data[[value_col]], na.rm = TRUE),
      y_min = min(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
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
    mutate(
      y = y_max + bracket_rank * y_step,
      y_tip = y - y_tip_size
    ) %>%
    select(any_of(facet_cols), contrast, x_start, x_end, y, y_tip, label, ReportingP, ReportingSignificance)
}

add_pairwise_brackets <- function(plot, bracket_tbl, text_size = 2.0) {
  if (nrow(bracket_tbl) == 0) return(plot)
  plot +
    geom_segment(
      data = bracket_tbl,
      aes(x = x_start, xend = x_end, y = y, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.25
    ) +
    geom_segment(
      data = bracket_tbl,
      aes(x = x_start, xend = x_start, y = y, yend = y_tip),
      inherit.aes = FALSE,
      linewidth = 0.25
    ) +
    geom_segment(
      data = bracket_tbl,
      aes(x = x_end, xend = x_end, y = y, yend = y_tip),
      inherit.aes = FALSE,
      linewidth = 0.25
    ) +
    geom_text(
      data = bracket_tbl,
      aes(x = (x_start + x_end) / 2, y = y, label = label),
      inherit.aes = FALSE,
      vjust = -0.25,
      size = text_size
    )
}

summarise_mean_ci <- function(dat, value_col, group_cols) {
  dat %>%
    filter(is.finite(.data[[value_col]])) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      n_animals = n_distinct(AnimalNum),
      n_obs = n(),
      mean = mean(.data[[value_col]], na.rm = TRUE),
      sd = sd(.data[[value_col]], na.rm = TRUE),
      sem = sd / sqrt(n_obs),
      ci95 = 1.96 * sem,
      ymin = mean - ci95,
      ymax = mean + ci95,
      .groups = "drop"
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

readout_tier_table <- function() {
  tibble(
    Metric = c("Movement", "Movement", "Movement", "Entropy", "Proximity", "Movement", "Movement"),
    Outcome = c("rmssd", "mean", "acf1", "rmssd", "rmssd", "cv", "fano"),
    ReadoutTier = c("Primary", "Secondary", "Secondary", "Secondary", "Secondary", "Exploratory", "Exploratory"),
    BehavioralReadout = c(
      "Locomotor temporal instability",
      "Overall locomotor output",
      "Behavioral inertia or persistence",
      "Instability in spatial exploration",
      "Instability in social contact",
      "Relative dispersion of behavior",
      "Overdispersion of behavior"
    ),
    Interpretation = c(
      "Best single readout for rapidly varying locomotor organization in long-term RFID tracking.",
      "Separates instability from simple hyperactivity or hypoactivity.",
      "Captures whether behavior is sticky and persistent versus rapidly changing.",
      "Tests whether animals destabilize their spatial-use strategy, not just movement amount.",
      "Tests whether social contact becomes temporally unstable.",
      "Useful robustness readout, but less temporally specific than RMSSD.",
      "Useful robustness readout, but strongly tied to mean-level count structure."
    )
  )
}

priority_rank <- function(tier) {
  recode(tier, Primary = 1L, Secondary = 2L, Exploratory = 3L, .default = 4L)
}

add_readout_context <- function(tbl, outcome_col = "Outcome") {
  readouts <- readout_tier_table()
  tbl %>%
    mutate(OutcomeJoin = as.character(.data[[outcome_col]]) %>% stringr::str_remove("^delta_")) %>%
    left_join(readouts, by = c("Metric", "OutcomeJoin" = "Outcome")) %>%
    mutate(
      ReadoutTier = coalesce(ReadoutTier, "Exploratory"),
      BehavioralReadout = coalesce(BehavioralReadout, paste(Metric, .data[[outcome_col]])),
      Interpretation = coalesce(Interpretation, "Exploratory temporal-dynamics readout.")
    ) %>%
    select(-OutcomeJoin)
}

make_finding_priority_table <- function(first_active_contrasts, overall_contrasts, delta_contrasts) {
  bind_rows(
    first_active_contrasts %>%
      mutate(AnalysisWindow = "Acute first active phase", TimeScale = "Acute"),
    overall_contrasts %>%
      mutate(AnalysisWindow = "Within cage change and phase", TimeScale = "Absolute"),
    delta_contrasts %>%
      mutate(AnalysisWindow = "Change from first cage change", TimeScale = "Longitudinal delta")
  ) %>%
    filter(status == "tested", contrast %in% pairwise_contrasts, !is.na(ReportingP)) %>%
    add_readout_context() %>%
    mutate(
      Direction = case_when(
        estimate > 0 ~ paste0(group_comp, " > ", group_ref),
        estimate < 0 ~ paste0(group_comp, " < ", group_ref),
        TRUE ~ "No direction"
      ),
      EvidenceClass = case_when(
        ReportingP < 0.05 & abs(cohen_d) >= 0.8 ~ "Strong effect",
        ReportingP < 0.05 ~ "FDR-supported",
        ReportingP < 0.10 ~ "Trend",
        abs(cohen_d) >= 0.8 ~ "Large effect, uncertain",
        TRUE ~ "No clear evidence"
      ),
      PriorityScore = case_when(
        is.finite(ReportingP) & ReportingP > 0 ~ -log10(ReportingP) * abs(cohen_d),
        TRUE ~ 0
      ),
      PriorityRank = priority_rank(ReadoutTier)
    ) %>%
    arrange(PriorityRank, desc(EvidenceClass %in% c("Strong effect", "FDR-supported")), desc(PriorityScore)) %>%
    select(
      AnalysisWindow, TimeScale, ReadoutTier, BehavioralReadout, Metric, Outcome,
      Sex, Phase, CageChange, contrast, Direction, estimate, conf.low, conf.high,
      cohen_d, ReportingP, ReportingSignificance, EvidenceClass, PriorityScore,
      Interpretation, everything()
    )
}

make_behavioral_pattern_table <- function(first_active_contrasts, overall_contrasts, delta_contrasts) {
  focal_readouts <- tibble(
    Metric = c("Movement", "Movement", "Movement", "Entropy", "Proximity"),
    Outcome = c("rmssd", "mean", "acf1", "rmssd", "rmssd")
  )

  acute_tbl <- first_active_contrasts %>%
    filter(status == "tested", contrast %in% pairwise_contrasts) %>%
    mutate(OutcomeBase = as.character(Outcome)) %>%
    semi_join(focal_readouts, by = c("Metric", "OutcomeBase" = "Outcome")) %>%
    mutate(
      Direction = case_when(
        estimate > 0 ~ paste0(group_comp, " > ", group_ref),
        estimate < 0 ~ paste0(group_comp, " < ", group_ref),
        TRUE ~ "No direction"
      )
    ) %>%
    group_by(Metric, OutcomeBase, Sex, contrast) %>%
    summarise(
      acute_min_p = min(ReportingP, na.rm = TRUE),
      acute_max_abs_d = max(abs(cohen_d), na.rm = TRUE),
      acute_direction = Direction[which.max(abs(cohen_d))],
      acute_supported = any(ReportingP < 0.05, na.rm = TRUE),
      .groups = "drop"
    )

  overall_tbl <- overall_contrasts %>%
    filter(status == "tested", contrast %in% pairwise_contrasts) %>%
    mutate(
      OutcomeBase = as.character(Outcome),
      CageChangeIndex = parse_cage_change_index(CageChange)
    ) %>%
    semi_join(focal_readouts, by = c("Metric", "OutcomeBase" = "Outcome")) %>%
    group_by(Metric, OutcomeBase, Sex, Phase, contrast) %>%
    summarise(
      first_cage_supported = any(CageChangeIndex == baseline_cage_change & ReportingP < 0.05, na.rm = TRUE),
      later_cage_supported = any(CageChangeIndex != baseline_cage_change & ReportingP < 0.05, na.rm = TRUE),
      n_supported_cage_changes = sum(ReportingP < 0.05, na.rm = TRUE),
      min_overall_p = min(ReportingP, na.rm = TRUE),
      max_overall_abs_d = max(abs(cohen_d), na.rm = TRUE),
      .groups = "drop"
    )

  delta_tbl <- delta_contrasts %>%
    filter(status == "tested", contrast %in% pairwise_contrasts) %>%
    mutate(
      OutcomeBase = stringr::str_remove(as.character(Outcome), "^delta_"),
      CageChangeIndex = parse_cage_change_index(CageChange)
    ) %>%
    semi_join(focal_readouts, by = c("Metric", "OutcomeBase" = "Outcome")) %>%
    group_by(Metric, OutcomeBase, Sex, Phase, contrast) %>%
    summarise(
      delta_supported = any(ReportingP < 0.05, na.rm = TRUE),
      delta_min_p = min(ReportingP, na.rm = TRUE),
      delta_max_abs_d = max(abs(cohen_d), na.rm = TRUE),
      delta_mean_estimate = mean(estimate[ReportingP < 0.10], na.rm = TRUE),
      .groups = "drop"
    )

  overall_tbl %>%
    left_join(acute_tbl, by = c("Metric", "OutcomeBase", "Sex", "contrast")) %>%
    left_join(delta_tbl, by = c("Metric", "OutcomeBase", "Sex", "Phase", "contrast")) %>%
    add_readout_context(outcome_col = "OutcomeBase") %>%
    mutate(
      across(c(acute_supported, first_cage_supported, later_cage_supported, delta_supported), ~replace_na(.x, FALSE)),
      Pattern = case_when(
        acute_supported & !later_cage_supported & !delta_supported ~ "Acute-only response",
        first_cage_supported & later_cage_supported ~ "Persistent group difference",
        !first_cage_supported & later_cage_supported ~ "Emergent later difference",
        delta_supported & delta_mean_estimate > 0 ~ "Sensitizing divergence from first cage change",
        delta_supported & delta_mean_estimate < 0 ~ "Habituating/converging change from first cage change",
        acute_supported & delta_supported ~ "Acute response with longitudinal change",
        TRUE ~ "No clear prioritized pattern"
      ),
      PhaseSpecificity = case_when(
        Phase == "Active" & (first_cage_supported | later_cage_supported | delta_supported) ~ "Active-phase evidence",
        Phase == "Inactive" & (first_cage_supported | later_cage_supported | delta_supported) ~ "Inactive-phase evidence",
        TRUE ~ "No phase-specific evidence"
      ),
      MinP = pmin(acute_min_p, min_overall_p, delta_min_p, na.rm = TRUE),
      MaxAbsD = pmax(acute_max_abs_d, max_overall_abs_d, delta_max_abs_d, na.rm = TRUE),
      PriorityScore = if_else(is.finite(MinP) & MinP > 0, -log10(MinP) * MaxAbsD, 0),
      PriorityRank = priority_rank(ReadoutTier)
    ) %>%
    arrange(PriorityRank, desc(Pattern != "No clear prioritized pattern"), desc(PriorityScore)) %>%
    select(
      ReadoutTier, BehavioralReadout, Metric, Outcome = OutcomeBase, Sex, Phase,
      contrast, Pattern, PhaseSpecificity, acute_supported, first_cage_supported,
      later_cage_supported, delta_supported, n_supported_cage_changes,
      MinP, MaxAbsD, PriorityScore, Interpretation, everything()
    )
}

# ------------------------------------------------
# LOAD + STANDARDIZE
# ------------------------------------------------

raw_dat <- read_behavior_table(input_file)
if (!proximity_col %in% names(raw_dat)) proximity_col <- "Proximity"

behav <- standardize_behavior_columns(raw_dat, proximity_col = proximity_col) %>%
  add_temporal_context()

ensure_dir(output_dir)
ensure_dir(file.path(output_dir, "tables"))
ensure_dir(file.path(output_dir, "stats_tables"))
ensure_dir(file.path(output_dir, "figures"))
ensure_dir(file.path(output_dir, "figures", "primary_movement"))
ensure_dir(file.path(output_dir, "figures", "supplementary_all_metrics"))
ensure_dir(file.path(output_dir, "figures", "first_active_first_cage_change"))
ensure_dir(file.path(output_dir, "figures", "cage_change_trajectories"))
ensure_dir(file.path(output_dir, "figures", "cage_change_deltas"))
ensure_dir(file.path(output_dir, "figures", "publication"))
ensure_dir(file.path(output_dir, "figures", "publication", "panels"))
ensure_dir(file.path(output_dir, "figures", "publication", "overview"))

# ------------------------------------------------
# LONG FORMAT
# ------------------------------------------------

long_dat <- behav %>%
  pivot_longer(
    cols = all_of(metrics_to_analyze),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    MetricRole = if_else(Metric == primary_metric, "Primary", "Secondary"),
    MetricRole = factor(MetricRole, levels = c("Primary", "Secondary")),
    Metric = factor(Metric, levels = metrics_to_analyze)
  )

movement_dat <- long_dat %>%
  filter(Metric == primary_metric)

# ------------------------------------------------
# INSTABILITY METRICS
# ------------------------------------------------

instability_tbl <- long_dat %>%
  group_by(Group, Sex, Phase, CageChange, CageChangeIndex, AnimalNum, Metric, MetricRole) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  summarise(calc_instability_metrics(Value), .groups = "drop") %>%
  mutate(BinLevel = bin_level, ProximityInput = proximity_col)

movement_instability_tbl <- instability_tbl %>%
  filter(Metric == primary_metric)

write_table(instability_tbl, file.path(output_dir, "tables", "temporal_instability_metrics_per_animal_all_metrics.csv"))
write_table(movement_instability_tbl, file.path(output_dir, "tables", "temporal_instability_metrics_per_animal_primary_movement.csv"))

overall_stats <- write_stats_package(
  instability_tbl,
  analysis_name = "overall_cage_change_phase_instability",
  value_cols = instability_outcomes,
  by_cols = c("BinLevel", "ProximityInput", "Metric", "Phase", "CageChange", "Sex"),
  summary_group_cols = c("BinLevel", "ProximityInput", "Metric", "Phase", "CageChange", "Sex", "Group")
)

baseline_cage_change <- min(instability_tbl$CageChangeIndex, na.rm = TRUE)

cage_change_delta_tbl <- instability_tbl %>%
  filter(is.finite(CageChangeIndex)) %>%
  select(
    BinLevel, ProximityInput, Group, Sex, Phase, AnimalNum, Metric, MetricRole,
    CageChange, CageChangeIndex, all_of(instability_outcomes)
  ) %>%
  left_join(
    instability_tbl %>%
      filter(CageChangeIndex == baseline_cage_change) %>%
      select(
        AnimalNum, Sex, Phase, Metric,
        baseline_CageChange = CageChange,
        all_of(instability_outcomes)
      ) %>%
      rename_with(~paste0("baseline_", .x), all_of(instability_outcomes)),
    by = c("AnimalNum", "Sex", "Phase", "Metric")
  ) %>%
  mutate(
    delta_mean = mean - baseline_mean,
    delta_cv = cv - baseline_cv,
    delta_fano = fano - baseline_fano,
    delta_rmssd = rmssd - baseline_rmssd,
    delta_acf1 = acf1 - baseline_acf1,
    CageChangeDelta = paste0(as.character(CageChange), " - ", as.character(baseline_CageChange))
  ) %>%
  filter(CageChangeIndex != baseline_cage_change)

movement_cage_change_delta_tbl <- cage_change_delta_tbl %>%
  filter(Metric == primary_metric)

write_table(cage_change_delta_tbl, file.path(output_dir, "tables", "cage_change_delta_instability_all_metrics.csv"))
write_table(movement_cage_change_delta_tbl, file.path(output_dir, "tables", "cage_change_delta_instability_primary_movement.csv"))

delta_stats <- write_stats_package(
  cage_change_delta_tbl,
  analysis_name = "delta_from_first_cage_change_instability",
  value_cols = paste0("delta_", instability_outcomes),
  by_cols = c("BinLevel", "ProximityInput", "Metric", "Phase", "CageChange", "Sex"),
  summary_group_cols = c("BinLevel", "ProximityInput", "Metric", "Phase", "CageChange", "Sex", "Group")
)

first_active_first_cage_change_long <- long_dat %>%
  filter(CageChangeIndex == baseline_cage_change, as.character(Phase) == "Active") %>%
  group_by(AnimalNum, CageChange, Metric) %>%
  mutate(FirstActiveBlock = min(PhaseBlock, na.rm = TRUE)) %>%
  filter(PhaseBlock == FirstActiveBlock) %>%
  group_by(AnimalNum, CageChange, Metric) %>%
  mutate(TimeFromFirstActiveStart = TimeIndex - min(TimeIndex, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Window = "First active phase of first cage change")

first_active_first_cage_change_instability_tbl <- first_active_first_cage_change_long %>%
  mutate(BinLevel = bin_level, ProximityInput = proximity_col) %>%
  group_by(BinLevel, ProximityInput, Window, Group, Sex, Phase, CageChange, CageChangeIndex, AnimalNum, Metric, MetricRole) %>%
  arrange(TimeFromFirstActiveStart, .by_group = TRUE) %>%
  summarise(calc_instability_metrics(Value), .groups = "drop")

movement_first_active_first_cage_change_instability_tbl <- first_active_first_cage_change_instability_tbl %>%
  filter(Metric == primary_metric)

write_table(first_active_first_cage_change_long, file.path(output_dir, "tables", "first_active_first_cage_change_timecourse_all_metrics.csv"))
write_table(first_active_first_cage_change_instability_tbl, file.path(output_dir, "tables", "first_active_first_cage_change_instability_all_metrics.csv"))
write_table(movement_first_active_first_cage_change_instability_tbl, file.path(output_dir, "tables", "first_active_first_cage_change_instability_primary_movement.csv"))

first_active_stats <- write_stats_package(
  first_active_first_cage_change_instability_tbl,
  analysis_name = "first_active_first_cage_change_instability",
  value_cols = instability_outcomes,
  by_cols = c("BinLevel", "ProximityInput", "Metric", "Sex"),
  summary_group_cols = c("BinLevel", "ProximityInput", "Metric", "Sex", "Group")
)

# ------------------------------------------------
# GROUP SUMMARY
# ------------------------------------------------

group_summary <- instability_tbl %>%
  group_by(BinLevel, ProximityInput, Group, Sex, Phase, CageChange, CageChangeIndex, Metric, MetricRole) %>%
  summarise(
    across(
      c(mean, cv, fano, rmssd, acf1),
      list(mean = ~mean(.x, na.rm = TRUE),
           sem = ~sd(.x, na.rm = TRUE) / sqrt(sum(is.finite(.x)))),
      .names = "{.col}_{.fn}"
    ),
    n_animals = n_distinct(AnimalNum),
    .groups = "drop"
  )

movement_group_summary <- group_summary %>%
  filter(Metric == primary_metric)

cage_change_levels <- instability_tbl %>%
  distinct(CageChange, CageChangeIndex) %>%
  arrange(CageChangeIndex, CageChange) %>%
  pull(CageChange) %>%
  as.character()

instability_tbl <- instability_tbl %>%
  mutate(CageChange = factor(as.character(CageChange), levels = cage_change_levels))
movement_instability_tbl <- movement_instability_tbl %>%
  mutate(CageChange = factor(as.character(CageChange), levels = cage_change_levels))
group_summary <- group_summary %>%
  mutate(CageChange = factor(as.character(CageChange), levels = cage_change_levels))
movement_group_summary <- movement_group_summary %>%
  mutate(CageChange = factor(as.character(CageChange), levels = cage_change_levels))
cage_change_delta_tbl <- cage_change_delta_tbl %>%
  mutate(CageChange = factor(as.character(CageChange), levels = cage_change_levels))
movement_cage_change_delta_tbl <- movement_cage_change_delta_tbl %>%
  mutate(CageChange = factor(as.character(CageChange), levels = cage_change_levels))

write_table(group_summary, file.path(output_dir, "tables", "temporal_instability_group_summary_all_metrics.csv"))
write_table(movement_group_summary, file.path(output_dir, "tables", "temporal_instability_group_summary_primary_movement.csv"))

# ------------------------------------------------
# PLOTS
# ------------------------------------------------

plot_metrics <- c("rmssd", "cv", "fano", "acf1")

for (metric_name in plot_metrics) {
  p_primary <- movement_instability_tbl %>%
    ggplot(aes(x = Group, y = .data[[metric_name]], fill = Group)) +
    geom_violin(alpha = 0.5, linewidth = 0.2, trim = FALSE) +
    geom_jitter(width = 0.1, size = 0.8, alpha = 0.7) +
    facet_grid(Sex ~ Phase, scales = "free_y") +
    labs(
      title = paste0(toupper(metric_name), " movement instability"),
      subtitle = paste0("Primary temporal-dynamics readout; bin level: ", bin_level),
      y = metric_name,
      x = NULL
    ) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    make_nature_theme()

  save_plot_svg_pdf(
    p_primary,
    file.path(output_dir, "figures", "primary_movement", paste0(metric_name, "_movement_group_comparison")),
    width = 180,
    height = 120
  )

  p_all <- instability_tbl %>%
    ggplot(aes(x = Group, y = .data[[metric_name]], fill = Group)) +
    geom_violin(alpha = 0.5, linewidth = 0.2, trim = FALSE) +
    geom_jitter(width = 0.1, size = 0.8, alpha = 0.7) +
    facet_grid(Metric ~ Phase, scales = "free_y") +
    labs(
      title = paste0(toupper(metric_name), " behavioral instability"),
      subtitle = paste0("Supplementary all-metric output; bin level: ", bin_level, "; proximity input: ", proximity_col),
      y = metric_name,
      x = NULL
    ) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    make_nature_theme()

  save_plot_svg_pdf(
    p_all,
    file.path(output_dir, "figures", "supplementary_all_metrics", paste0(metric_name, "_all_metrics_group_comparison")),
    width = 180,
    height = 120
  )
}

# ------------------------------------------------
# TIME-RESOLVED VOLATILITY
# ------------------------------------------------

rolling_tbl <- long_dat %>%
  group_by(Group, Sex, Phase, CageChange, AnimalNum, Metric, MetricRole) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  mutate(
    rolling_sd = zoo::rollapply(Value, width = 3, FUN = sd, align = "right", fill = NA, na.rm = TRUE),
    rolling_rmssd = zoo::rollapply(Value, width = 4, FUN = calc_rmssd, align = "right", fill = NA)
  ) %>%
  ungroup() %>%
  mutate(BinLevel = bin_level, ProximityInput = proximity_col)

movement_rolling_tbl <- rolling_tbl %>%
  filter(Metric == primary_metric)

write_table(rolling_tbl, file.path(output_dir, "tables", "rolling_instability_metrics_all_metrics.csv"))
write_table(movement_rolling_tbl, file.path(output_dir, "tables", "rolling_instability_metrics_primary_movement.csv"))

p_roll_primary <- movement_rolling_tbl %>%
  group_by(Group, Sex, Phase, CageChange, TimeFromPhaseStart) %>%
  summarise(
    mean_rmssd = mean(rolling_rmssd, na.rm = TRUE),
    sem_rmssd = sd(rolling_rmssd, na.rm = TRUE) / sqrt(sum(is.finite(rolling_rmssd))),
    .groups = "drop"
  ) %>%
  ggplot(aes(TimeFromPhaseStart, mean_rmssd, colour = Group, fill = Group)) +
  geom_line(linewidth = 0.5) +
  geom_ribbon(aes(ymin = mean_rmssd - sem_rmssd, ymax = mean_rmssd + sem_rmssd), alpha = 0.2, linewidth = 0) +
  facet_grid(Sex + Phase ~ CageChange, scales = "free_y") +
  labs(
    title = "Phase-aligned movement instability within cage changes",
    subtitle = paste0("Primary temporal-dynamics readout; x-axis resets within each phase block; bin level: ", bin_level),
    y = "Rolling movement RMSSD",
    x = "Bins from phase onset"
  ) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  make_nature_theme()

save_plot_svg_pdf(
  p_roll_primary,
  file.path(output_dir, "figures", "primary_movement", "rolling_rmssd_movement_timecourse"),
  width = 210,
  height = 160
)

p_roll_all <- rolling_tbl %>%
  group_by(Group, Metric, Phase, CageChange, TimeFromPhaseStart) %>%
  summarise(
    mean_rmssd = mean(rolling_rmssd, na.rm = TRUE),
    sem_rmssd = sd(rolling_rmssd, na.rm = TRUE) / sqrt(sum(is.finite(rolling_rmssd))),
    .groups = "drop"
  ) %>%
  ggplot(aes(TimeFromPhaseStart, mean_rmssd, colour = Group, fill = Group)) +
  geom_line(linewidth = 0.5) +
  geom_ribbon(aes(ymin = mean_rmssd - sem_rmssd, ymax = mean_rmssd + sem_rmssd), alpha = 0.2, linewidth = 0) +
  facet_grid(Metric + Phase ~ CageChange, scales = "free_y") +
  labs(
    title = "Phase-aligned behavioral instability within cage changes",
    subtitle = paste0("Supplementary all-metric output; x-axis resets within each phase block; bin level: ", bin_level),
    y = "Rolling RMSSD",
    x = "Bins from phase onset"
  ) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  make_nature_theme()

save_plot_svg_pdf(
  p_roll_all,
  file.path(output_dir, "figures", "supplementary_all_metrics", "rolling_rmssd_all_metrics_timecourse"),
  width = 220,
  height = 190
)

# ------------------------------------------------
# FIRST ACTIVE PHASE OF FIRST CAGE CHANGE
# ------------------------------------------------

first_active_plot_metrics <- c("rmssd", "cv", "acf1")

for (metric_name in first_active_plot_metrics) {
  stat_brackets <- make_pairwise_brackets(
    movement_first_active_first_cage_change_instability_tbl,
    first_active_stats$contrasts,
    outcome = metric_name,
    facet_cols = c("Sex"),
    value_col = metric_name
  )

  p_first_active <- movement_first_active_first_cage_change_instability_tbl %>%
    ggplot(aes(x = Group, y = .data[[metric_name]], fill = Group)) +
    geom_violin(alpha = 0.55, linewidth = 0.25, trim = FALSE, colour = "grey25") +
    geom_jitter(width = 0.08, size = 0.9, alpha = 0.75, shape = 21, colour = "grey20", stroke = 0.15, show.legend = FALSE) +
    facet_grid(. ~ Sex, scales = "free_y") +
    labs(
      title = paste0(toupper(metric_name), " during first active phase"),
      subtitle = paste0("First active phase of first cage change; bars show BH-adjusted pairwise p-values"),
      y = metric_name,
      x = NULL
    ) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    make_nature_theme()

  p_first_active <- add_pairwise_brackets(p_first_active, stat_brackets)

  save_plot_svg_pdf(
    p_first_active,
    file.path(output_dir, "figures", "first_active_first_cage_change", paste0(metric_name, "_movement_first_active_first_cage_change")),
    width = 150,
    height = 90
  )
}

first_active_rolling_tbl <- first_active_first_cage_change_long %>%
  group_by(Group, Sex, AnimalNum, Metric, MetricRole) %>%
  arrange(TimeFromFirstActiveStart, .by_group = TRUE) %>%
  mutate(
    rolling_sd = zoo::rollapply(Value, width = 3, FUN = sd, align = "right", fill = NA, na.rm = TRUE),
    rolling_rmssd = zoo::rollapply(Value, width = 4, FUN = calc_rmssd, align = "right", fill = NA)
  ) %>%
  ungroup() %>%
  mutate(BinLevel = bin_level, ProximityInput = proximity_col)

write_table(first_active_rolling_tbl, file.path(output_dir, "tables", "first_active_first_cage_change_rolling_instability_all_metrics.csv"))

p_first_active_roll <- first_active_rolling_tbl %>%
  filter(Metric == primary_metric) %>%
  group_by(Group, Sex, TimeFromFirstActiveStart) %>%
  summarise(
    mean_rmssd = mean(rolling_rmssd, na.rm = TRUE),
    sem_rmssd = sd(rolling_rmssd, na.rm = TRUE) / sqrt(sum(is.finite(rolling_rmssd))),
    .groups = "drop"
  ) %>%
  ggplot(aes(TimeFromFirstActiveStart, mean_rmssd, colour = Group, fill = Group)) +
  geom_line(linewidth = 0.55) +
  geom_ribbon(aes(ymin = mean_rmssd - sem_rmssd, ymax = mean_rmssd + sem_rmssd), alpha = 0.18, linewidth = 0) +
  facet_grid(. ~ Sex, scales = "free_y") +
  labs(
    title = "Movement instability within first active phase",
    subtitle = paste0("Rolling RMSSD aligned to active-phase onset; bin level: ", bin_level),
    y = "Rolling movement RMSSD",
    x = "Bins from active-phase onset"
  ) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  make_nature_theme()

save_plot_svg_pdf(
  p_first_active_roll,
  file.path(output_dir, "figures", "first_active_first_cage_change", "rolling_rmssd_movement_first_active_first_cage_change"),
  width = 150,
  height = 90
)

# ------------------------------------------------
# CAGE-CHANGE TRAJECTORIES AND DELTAS
# ------------------------------------------------

cage_change_sig_labels <- overall_stats$contrasts %>%
  filter(
    Metric == primary_metric,
    Outcome == "rmssd",
    contrast %in% pairwise_contrasts,
    !is.na(ReportingP),
    ReportingP < 0.10
  ) %>%
  mutate(label = paste0(contrast, " ", format_p_label(ReportingP))) %>%
  group_by(Sex, Phase, CageChange) %>%
  summarise(label = paste(label, collapse = "\n"), .groups = "drop") %>%
  mutate(CageChange = factor(as.character(CageChange), levels = cage_change_levels)) %>%
  left_join(
    movement_instability_tbl %>%
      group_by(Sex, Phase, CageChange) %>%
      summarise(
        y_pos = max(rmssd, na.rm = TRUE),
        y_min = min(rmssd, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        y_range = if_else(is.finite(y_pos - y_min) & y_pos > y_min, y_pos - y_min, abs(y_pos)),
        y_pos = y_pos + if_else(is.finite(y_range) & y_range > 0, 0.15 * y_range, 0.1)
      ),
    by = c("Sex", "Phase", "CageChange")
  )

p_cage_rmssd <- movement_group_summary %>%
  ggplot(aes(CageChange, rmssd_mean, colour = Group, fill = Group, group = Group)) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = rmssd_mean - rmssd_sem, ymax = rmssd_mean + rmssd_sem), width = 0.08, linewidth = 0.25) +
  geom_text(
    data = cage_change_sig_labels,
    aes(x = CageChange, y = y_pos, label = label),
    inherit.aes = FALSE,
    size = 2.2,
    lineheight = 0.85
  ) +
  facet_grid(Sex ~ Phase, scales = "free_y") +
  labs(
    title = "Movement instability across cage changes",
    subtitle = "Mean +/- SEM per animal; labels show BH-adjusted pairwise p-values",
    y = "Movement RMSSD",
    x = "Cage change"
  ) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  make_nature_theme()

save_plot_svg_pdf(
  p_cage_rmssd,
  file.path(output_dir, "figures", "cage_change_trajectories", "rmssd_movement_across_cage_changes"),
  width = 180,
  height = 120
)

movement_rmssd_stat_brackets <- make_pairwise_brackets(
  movement_instability_tbl,
  overall_stats$contrasts,
  outcome = "rmssd",
  facet_cols = c("Sex", "Phase", "CageChange"),
  value_col = "rmssd"
) %>%
  mutate(CageChange = factor(as.character(CageChange), levels = cage_change_levels))

p_movement_rmssd_pairwise <- movement_instability_tbl %>%
  ggplot(aes(Group, rmssd, fill = Group)) +
  geom_violin(alpha = 0.50, linewidth = 0.22, trim = FALSE, colour = "grey25") +
  geom_jitter(width = 0.08, size = 0.7, alpha = 0.65, shape = 21, colour = "grey20", stroke = 0.12, show.legend = FALSE) +
  facet_grid(Sex + Phase ~ CageChange, scales = "free_y") +
  labs(
    title = "Movement instability group differences across cage changes",
    subtitle = "Bars show BH-adjusted p-values for all significant pairwise group contrasts",
    y = "Movement RMSSD",
    x = NULL
  ) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  make_nature_theme()

p_movement_rmssd_pairwise <- add_pairwise_brackets(p_movement_rmssd_pairwise, movement_rmssd_stat_brackets, text_size = 1.7)

save_plot_svg_pdf(
  p_movement_rmssd_pairwise,
  file.path(output_dir, "figures", "primary_movement", "rmssd_movement_pairwise_bars_by_cage_change"),
  width = 210,
  height = 150
)

delta_stat_brackets <- make_pairwise_brackets(
  movement_cage_change_delta_tbl,
  delta_stats$contrasts,
  outcome = "delta_rmssd",
  facet_cols = c("Sex", "Phase", "CageChange"),
  value_col = "delta_rmssd"
) %>%
  mutate(CageChange = factor(as.character(CageChange), levels = cage_change_levels))

p_delta_rmssd <- movement_cage_change_delta_tbl %>%
  ggplot(aes(Group, delta_rmssd, fill = Group)) +
  geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey35") +
  geom_violin(alpha = 0.55, linewidth = 0.25, trim = FALSE, colour = "grey25") +
  geom_jitter(width = 0.08, size = 0.8, alpha = 0.7, shape = 21, colour = "grey20", stroke = 0.12, show.legend = FALSE) +
  facet_grid(Sex + Phase ~ CageChange, scales = "free_y") +
  labs(
    title = "Change in movement instability from first cage change",
    subtitle = "Delta RMSSD relative to first cage change; bars show BH-adjusted pairwise p-values",
    y = "Delta movement RMSSD",
    x = NULL
  ) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  make_nature_theme()

p_delta_rmssd <- add_pairwise_brackets(p_delta_rmssd, delta_stat_brackets, text_size = 1.8)

save_plot_svg_pdf(
  p_delta_rmssd,
  file.path(output_dir, "figures", "cage_change_deltas", "delta_rmssd_movement_from_first_cage_change"),
  width = 190,
  height = 130
)

# ------------------------------------------------
# PUBLICATION-STYLE SUMMARY FIGURES
# ------------------------------------------------

publication_panel_dir <- file.path(output_dir, "figures", "publication", "panels")
publication_overview_dir <- file.path(output_dir, "figures", "publication", "overview")

stress_stage_colors <- c(
  "Acute first exposure" = "#f3e8e5",
  "Early adaptation" = "#eef1f7",
  "Repeated exposure" = "#f4f4f0"
)

cage_change_order_tbl <- tibble(
  CageChange = factor(cage_change_levels, levels = cage_change_levels),
  CageChangeOrder = seq_along(cage_change_levels),
  CageChangeIndex = parse_cage_change_index(cage_change_levels)
)

followup_orders <- cage_change_order_tbl %>%
  filter(CageChangeIndex != baseline_cage_change) %>%
  pull(CageChangeOrder)

first_followup_order <- if (any(is.finite(followup_orders))) min(followup_orders, na.rm = TRUE) else NA_real_

stress_timing_tbl <- cage_change_order_tbl %>%
  mutate(
    StressModelWindow = case_when(
      CageChangeIndex == baseline_cage_change ~ "Acute first exposure",
      CageChangeOrder == first_followup_order ~ "Early adaptation",
      TRUE ~ "Repeated exposure"
    ),
    BiologicalTiming = case_when(
      StressModelWindow == "Acute first exposure" ~ "Initial response to the first cage-change challenge",
      StressModelWindow == "Early adaptation" ~ "First follow-up cage change after the initial challenge",
      TRUE ~ "Later repeated cage-change challenges"
    )
  )

stress_window_rect_tbl <- stress_timing_tbl %>%
  group_by(StressModelWindow) %>%
  summarise(
    xmin = min(CageChangeOrder) - 0.5,
    xmax = max(CageChangeOrder) + 0.5,
    .groups = "drop"
  ) %>%
  mutate(StressModelWindow = factor(StressModelWindow, levels = names(stress_stage_colors)))

stress_phase_tbl <- tibble(
  Phase = c("Active", "Inactive"),
  BiologicalTiming = c(
    "Active phase: high-activity/dark-cycle window, expected to reveal locomotor disruption",
    "Inactive phase: low-activity/light-cycle window, expected to reveal rest-state disruption"
  )
)

stress_comparison_tbl <- tibble(
  contrast = pairwise_contrasts,
  BiologicalComparison = c(
    "Resilient vs control",
    "Susceptible vs control",
    "Susceptible vs resilient"
  ),
  InterpretationFocus = c(
    "Stress exposure without susceptible phenotype",
    "Susceptible phenotype relative to non-stressed controls",
    "Phenotype specificity within stressed animals"
  )
)

write_table(stress_timing_tbl, file.path(output_dir, "tables", "publication_stress_model_timing_map.csv"))
write_table(stress_phase_tbl, file.path(output_dir, "tables", "publication_phase_timing_map.csv"))
write_table(stress_comparison_tbl, file.path(output_dir, "tables", "publication_stress_model_comparisons.csv"))

readout_definitions <- readout_tier_table()
write_table(readout_definitions, file.path(output_dir, "tables", "behavioral_readout_definitions.csv"))

finding_priority_tbl <- make_finding_priority_table(
  first_active_stats$contrasts,
  overall_stats$contrasts,
  delta_stats$contrasts
)

behavioral_pattern_tbl <- make_behavioral_pattern_table(
  first_active_stats$contrasts,
  overall_stats$contrasts,
  delta_stats$contrasts
)

priority_findings_shortlist <- finding_priority_tbl %>%
  filter(
    ReadoutTier %in% c("Primary", "Secondary"),
    EvidenceClass %in% c("Strong effect", "FDR-supported", "Trend", "Large effect, uncertain")
  ) %>%
  group_by(ReadoutTier, BehavioralReadout, Metric, Outcome, AnalysisWindow, contrast, Sex) %>%
  slice_max(order_by = PriorityScore, n = 2, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(priority_rank(ReadoutTier), desc(EvidenceClass %in% c("Strong effect", "FDR-supported")), desc(PriorityScore)) %>%
  slice_head(n = 50)

key_behavioral_patterns <- behavioral_pattern_tbl %>%
  filter(
    ReadoutTier %in% c("Primary", "Secondary"),
    Pattern != "No clear prioritized pattern"
  ) %>%
  arrange(priority_rank(ReadoutTier), desc(PriorityScore)) %>%
  select(
    ReadoutTier, BehavioralReadout, Metric, Outcome, Sex, Phase, contrast,
    Pattern, PhaseSpecificity, acute_supported, first_cage_supported,
    later_cage_supported, delta_supported, n_supported_cage_changes,
    MinP, MaxAbsD, PriorityScore, Interpretation
  )

write_table(finding_priority_tbl, file.path(output_dir, "tables", "behavioral_finding_priority_all_tests.csv"))
write_table(priority_findings_shortlist, file.path(output_dir, "tables", "behavioral_finding_priority_shortlist.csv"))
write_table(behavioral_pattern_tbl, file.path(output_dir, "tables", "behavioral_pattern_classification.csv"))
write_table(key_behavioral_patterns, file.path(output_dir, "tables", "behavioral_key_patterns.csv"))

priority_plot_tbl <- priority_findings_shortlist %>%
  mutate(
    PlotLabel = paste(AnalysisWindow, Sex, coalesce(as.character(Phase), "all phases"), coalesce(as.character(CageChange), "first cage"), sep = " | "),
    PlotLabel = factor(PlotLabel, levels = rev(unique(PlotLabel[order(ReadoutTier, BehavioralReadout, PriorityScore)]))),
    ReadoutTier = factor(ReadoutTier, levels = c("Primary", "Secondary", "Exploratory")),
    EvidenceClass = factor(EvidenceClass, levels = c("Strong effect", "FDR-supported", "Trend", "Large effect, uncertain", "No clear evidence"))
  )

if (nrow(priority_plot_tbl) > 0) {
  p_priority_readouts <- ggplot(priority_plot_tbl, aes(PriorityScore, PlotLabel, colour = contrast, shape = EvidenceClass)) +
    geom_point(size = 1.8) +
    facet_grid(ReadoutTier + BehavioralReadout ~ ., scales = "free_y", space = "free_y") +
    labs(
      title = "Prioritized behavioral readouts",
      subtitle = "Ranked by adjusted evidence and effect size; shortlist excludes unsupported tests",
      x = "Priority score",
      y = NULL,
      colour = NULL,
      shape = NULL
    ) +
    scale_colour_manual(values = c("RES-CON" = "#8f8a82", "SUS-CON" = group_colors[["SUS"]], "SUS-RES" = "#5f5a54"), drop = FALSE) +
    scale_shape_manual(values = c("Strong effect" = 16, "FDR-supported" = 16, "Trend" = 1, "Large effect, uncertain" = 2, "No clear evidence" = 4), drop = FALSE) +
    make_publication_theme(base_size = 6) +
    theme(panel.grid.major.y = element_blank())

  save_plot_svg_pdf(
    p_priority_readouts,
    file.path(publication_overview_dir, "behavioral_readout_priority_shortlist"),
    width = 190,
    height = 140
  )
}

pattern_plot_tbl <- behavioral_pattern_tbl %>%
  filter(ReadoutTier %in% c("Primary", "Secondary")) %>%
  mutate(
    Pattern = factor(
      Pattern,
      levels = c(
        "Persistent group difference",
        "Emergent later difference",
        "Sensitizing divergence from first cage change",
        "Habituating/converging change from first cage change",
        "Acute response with longitudinal change",
        "Acute-only response",
        "No clear prioritized pattern"
      )
    ),
    ReadoutLabel = paste(ReadoutTier, BehavioralReadout, sep = ": "),
    ReadoutLabel = factor(ReadoutLabel, levels = rev(unique(ReadoutLabel[order(priority_rank(ReadoutTier), BehavioralReadout)])))
  )

p_behavioral_patterns <- ggplot(pattern_plot_tbl, aes(Phase, ReadoutLabel, fill = Pattern)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = if_else(Pattern == "No clear prioritized pattern", "", sprintf("%.2f", MaxAbsD))), size = 2.0) +
  facet_grid(Sex ~ contrast) +
  labs(
    title = "Behavioral pattern classification",
    subtitle = "Tiles summarize acute, persistent, and longitudinal-change evidence; numbers show max |d|",
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  scale_fill_manual(
    values = c(
      "Persistent group difference" = "#e63947",
      "Emergent later difference" = "#f08a94",
      "Sensitizing divergence from first cage change" = "#b71d2a",
      "Habituating/converging change from first cage change" = "#3d3b6e",
      "Acute response with longitudinal change" = "#8f3f71",
      "Acute-only response" = "#C6C3BB",
      "No clear prioritized pattern" = "grey92"
    ),
    drop = FALSE
  ) +
  make_publication_theme(base_size = 6) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(7, "mm")
  )

save_plot_svg_pdf(
  p_behavioral_patterns,
  file.path(publication_overview_dir, "behavioral_pattern_classification"),
  width = 190,
  height = 125
)

timeline_tbl <- stress_timing_tbl %>%
  mutate(
    xmin = CageChangeOrder - 0.42,
    xmax = CageChangeOrder + 0.42,
    ymin = 0,
    ymax = 1,
    Fill = StressModelWindow
  )

timeline_label_tbl <- stress_window_rect_tbl %>%
  mutate(
    x = (xmin + xmax) / 2,
    y = 1.12,
    label = StressModelWindow
  )

phase_schematic_tbl <- tibble(
  Segment = c("Inactive", "Active"),
  xmin = c(0.05, 0.52),
  xmax = c(0.48, 0.95),
  ymin = c(1.28, 1.28),
  ymax = c(1.58, 1.58),
  Fill = c("Inactive", "Active")
)

rmssd_trace_tbl <- tibble(
  Time = seq_len(16),
  Movement = c(1, 1, 2, 1, 7, 8, 2, 1, 1, 5, 7, 6, 1, 1, 2, 1)
) %>%
  mutate(
    Diff = c(NA, diff(Movement)),
    SquaredDiff = Diff^2
  )

p_pub_schematic <- ggplot() +
  geom_rect(
    data = timeline_tbl,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Fill),
    colour = "grey25",
    linewidth = 0.25
  ) +
  geom_text(
    data = timeline_tbl,
    aes(x = CageChangeOrder, y = 0.48, label = CageChange),
    size = 2.5
  ) +
  geom_text(
    data = timeline_label_tbl,
    aes(x = x, y = y, label = label),
    size = 2.25
  ) +
  geom_rect(
    data = phase_schematic_tbl,
    aes(xmin = xmin * length(cage_change_levels) + 0.5,
        xmax = xmax * length(cage_change_levels) + 0.5,
        ymin = ymin + 0.22, ymax = ymax + 0.22, fill = Fill),
    colour = "grey25",
    linewidth = 0.2
  ) +
  geom_text(
    data = phase_schematic_tbl,
    aes(x = ((xmin + xmax) / 2) * length(cage_change_levels) + 0.5, y = 1.65, label = Segment),
    size = 2.3
  ) +
  geom_line(
    data = rmssd_trace_tbl,
    aes(x = scales::rescale(Time, to = c(0.65, length(cage_change_levels) + 0.35)),
        y = scales::rescale(Movement, to = c(1.88, 2.55))),
    linewidth = 0.35,
    colour = "grey15"
  ) +
  annotate("text", x = 0.55, y = 2.68, hjust = 0, size = 2.4, label = "RMSSD summarizes step-to-step volatility") +
  annotate("text", x = 0.55, y = 1.92, hjust = 0, size = 2.4, label = "Active and inactive phases test state-dependent disruption") +
  annotate("text", x = 0.55, y = 1.32, hjust = 0, size = 2.4, label = "Primary comparisons: SUS-CON, RES-CON, SUS-RES") +
  scale_fill_manual(values = c(stress_stage_colors, "Inactive" = "#d9d9d9", "Active" = "#ffffff")) +
  coord_cartesian(xlim = c(0.45, length(cage_change_levels) + 0.55), ylim = c(-0.05, 2.85), clip = "off") +
  labs(title = "Stress-model timing for temporal instability", x = NULL, y = NULL) +
  theme_void(base_size = 7) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0, size = 8),
    plot.margin = margin(4, 8, 4, 8)
  )

save_plot_svg_pdf(
  p_pub_schematic,
  file.path(publication_panel_dir, "panel_a_timeline_metric_schematic"),
  width = 180,
  height = 70
)

pub_first_active_summary <- first_active_rolling_tbl %>%
  filter(Metric == primary_metric) %>%
  summarise_mean_ci("rolling_rmssd", c("Group", "Sex", "TimeFromFirstActiveStart")) %>%
  filter(is.finite(mean), is.finite(ymin), is.finite(ymax))

write_table(pub_first_active_summary, file.path(output_dir, "tables", "publication_first_active_rolling_rmssd_summary.csv"))

p_pub_first_active <- ggplot(pub_first_active_summary, aes(TimeFromFirstActiveStart, mean, colour = Group, fill = Group)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.16, linewidth = 0) +
  geom_line(linewidth = 0.65) +
  facet_grid(. ~ Sex, scales = "free_y") +
  labs(
    title = "Acute active-phase movement instability",
    subtitle = "First cage-change active phase; mean rolling RMSSD +/- 95% CI",
    x = "Bins from active-phase onset",
    y = "Rolling movement RMSSD"
  ) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  make_publication_theme()

save_plot_svg_pdf(
  p_pub_first_active,
  file.path(publication_panel_dir, "panel_b_first_active_rolling_rmssd"),
  width = 160,
  height = 75
)

pub_movement_instability <- movement_instability_tbl %>%
  left_join(
    stress_timing_tbl %>%
      select(CageChange, CageChangeOrder, StressModelWindow),
    by = "CageChange"
  )

pub_cage_summary <- summarise_mean_ci(
  pub_movement_instability,
  "rmssd",
  c("Group", "Sex", "Phase", "CageChange", "CageChangeIndex", "CageChangeOrder", "StressModelWindow")
) %>%
  mutate(CageChange = factor(as.character(CageChange), levels = cage_change_levels))

write_table(pub_cage_summary, file.path(output_dir, "tables", "publication_cage_change_rmssd_summary.csv"))

p_pub_cage_trajectory <- ggplot() +
  geom_vline(
    data = stress_window_rect_tbl %>% filter(xmin > 0.5),
    aes(xintercept = xmin),
    inherit.aes = FALSE,
    linewidth = 0.2,
    linetype = "dashed",
    colour = "grey70"
  ) +
  geom_text(
    data = timeline_label_tbl,
    aes(x = x, y = Inf, label = label),
    inherit.aes = FALSE,
    vjust = 1.25,
    size = 2.1,
    colour = "grey25"
  ) +
  geom_line(
    data = pub_movement_instability,
    aes(CageChangeOrder, rmssd, group = interaction(AnimalNum, Group), colour = Group),
    alpha = 0.16,
    linewidth = 0.25
  ) +
  geom_line(
    data = pub_cage_summary,
    aes(CageChangeOrder, mean, group = Group, colour = Group),
    linewidth = 0.7
  ) +
  geom_errorbar(
    data = pub_cage_summary,
    aes(CageChangeOrder, ymin = ymin, ymax = ymax, colour = Group),
    width = 0.10,
    linewidth = 0.28
  ) +
  geom_point(
    data = pub_cage_summary,
    aes(CageChangeOrder, mean, fill = Group),
    shape = 21,
    colour = "grey20",
    stroke = 0.25,
    size = 1.6
  ) +
  facet_grid(Sex ~ Phase, scales = "free_y") +
  labs(
    title = "Movement instability across stress-model windows",
    subtitle = "Thin lines show animals; bold lines and error bars show group mean +/- 95% CI",
    x = "Cage change",
    y = "Movement RMSSD"
  ) +
  scale_x_continuous(breaks = stress_timing_tbl$CageChangeOrder, labels = as.character(stress_timing_tbl$CageChange)) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  make_publication_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

save_plot_svg_pdf(
  p_pub_cage_trajectory,
  file.path(publication_panel_dir, "panel_c_cage_change_rmssd_trajectory"),
  width = 180,
  height = 115
)

pub_movement_delta <- movement_cage_change_delta_tbl %>%
  left_join(
    stress_timing_tbl %>%
      select(CageChange, CageChangeOrder, StressModelWindow),
    by = "CageChange"
  )

pub_delta_summary <- summarise_mean_ci(
  pub_movement_delta,
  "delta_rmssd",
  c("Group", "Sex", "Phase", "CageChange", "CageChangeIndex", "CageChangeOrder", "StressModelWindow")
) %>%
  mutate(CageChange = factor(as.character(CageChange), levels = cage_change_levels))

write_table(pub_delta_summary, file.path(output_dir, "tables", "publication_delta_rmssd_summary.csv"))

p_pub_delta <- ggplot(pub_delta_summary, aes(CageChangeOrder, mean, colour = Group, group = Group)) +
  geom_vline(
    data = stress_window_rect_tbl %>% filter(xmin > min(pub_delta_summary$CageChangeOrder, na.rm = TRUE) - 0.5),
    aes(xintercept = xmin),
    inherit.aes = FALSE,
    linewidth = 0.2,
    linetype = "dashed",
    colour = "grey70"
  ) +
  geom_text(
    data = timeline_label_tbl %>% filter(x >= min(pub_delta_summary$CageChangeOrder, na.rm = TRUE)),
    aes(x = x, y = Inf, label = label),
    inherit.aes = FALSE,
    vjust = 1.25,
    size = 2.1,
    colour = "grey25"
  ) +
  geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey35") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.08, linewidth = 0.28) +
  geom_line(linewidth = 0.55) +
  geom_point(aes(fill = Group), size = 1.7, shape = 21, stroke = 0.25) +
  facet_grid(Sex ~ Phase, scales = "free_y") +
  labs(
    title = "Adaptation of movement instability after the first exposure",
    subtitle = "Delta RMSSD from first cage change; mean +/- 95% CI",
    x = "Cage change",
    y = "Delta movement RMSSD"
  ) +
  scale_x_continuous(breaks = stress_timing_tbl$CageChangeOrder, labels = as.character(stress_timing_tbl$CageChange)) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  make_publication_theme() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

save_plot_svg_pdf(
  p_pub_delta,
  file.path(publication_panel_dir, "panel_d_delta_rmssd_from_first_cage_change"),
  width = 180,
  height = 115
)

pub_forest_tbl <- bind_rows(
  overall_stats$contrasts %>%
    filter(Metric == primary_metric, Outcome == "rmssd", contrast %in% pairwise_contrasts) %>%
    mutate(Analysis = "Absolute instability"),
  delta_stats$contrasts %>%
    filter(Metric == primary_metric, Outcome == "delta_rmssd", contrast %in% pairwise_contrasts) %>%
    mutate(Analysis = "Adaptation from first exposure")
) %>%
  filter(status == "tested") %>%
  left_join(
    stress_timing_tbl %>%
      select(CageChange, CageChangeOrder, StressModelWindow),
    by = "CageChange"
  ) %>%
  left_join(stress_comparison_tbl, by = "contrast") %>%
  mutate(
    OutcomeLabel = recode(Outcome, rmssd = "RMSSD", delta_rmssd = "Delta RMSSD"),
    PanelLabel = paste(StressModelWindow, Sex, Phase, as.character(CageChange), sep = " | "),
    PanelLabel = factor(PanelLabel, levels = rev(unique(PanelLabel[order(CageChangeOrder, Phase, Sex)]))),
    contrast = factor(contrast, levels = pairwise_contrasts),
    BiologicalComparison = factor(BiologicalComparison, levels = stress_comparison_tbl$BiologicalComparison),
    SignificantLabel = if_else(ReportingP < 0.05, "FDR < 0.05", "n.s.")
  )

write_table(pub_forest_tbl, file.path(output_dir, "tables", "publication_pairwise_contrast_forest_data.csv"))

p_pub_forest <- ggplot(pub_forest_tbl, aes(x = estimate, y = PanelLabel, colour = contrast, shape = SignificantLabel)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey45") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", width = 0, linewidth = 0.28, alpha = 0.9) +
  geom_point(size = 1.6) +
  facet_grid(Analysis ~ BiologicalComparison, scales = "free_y", space = "free_y") +
  labs(
    title = "Stress-model contrasts for movement instability",
    subtitle = "Welch mean differences with 95% confidence intervals; rows encode timing, sex, phase, and cage change",
    x = "Mean difference",
    y = NULL,
    shape = NULL
  ) +
  scale_colour_manual(values = c("RES-CON" = "#8f8a82", "SUS-CON" = group_colors[["SUS"]], "SUS-RES" = "#5f5a54"), drop = FALSE) +
  scale_shape_manual(values = c("FDR < 0.05" = 16, "n.s." = 1), drop = FALSE) +
  make_publication_theme() +
  theme(panel.grid.major.y = element_blank())

save_plot_svg_pdf(
  p_pub_forest,
  file.path(publication_panel_dir, "panel_e_pairwise_contrast_forest"),
  width = 190,
  height = 130
)

pub_heat_tbl <- overall_stats$contrasts %>%
  filter(contrast %in% pairwise_contrasts, status == "tested") %>%
  left_join(
    stress_timing_tbl %>%
      select(CageChange, CageChangeOrder, StressModelWindow),
    by = "CageChange"
  ) %>%
  left_join(stress_comparison_tbl, by = "contrast") %>%
  mutate(
    Outcome = factor(Outcome, levels = instability_outcomes, labels = metric_labels[instability_outcomes]),
    Metric = factor(as.character(Metric), levels = metrics_to_analyze),
    CagePhase = paste(as.character(CageChange), Phase, sep = "\n"),
    CagePhase = factor(CagePhase, levels = unique(CagePhase[order(CageChangeOrder, Phase)])),
    contrast = factor(contrast, levels = pairwise_contrasts),
    BiologicalComparison = factor(BiologicalComparison, levels = stress_comparison_tbl$BiologicalComparison),
    SigLabel = if_else(ReportingP < 0.05, "*", "")
  )

write_table(pub_heat_tbl, file.path(output_dir, "tables", "publication_instability_effect_size_heatmap_data.csv"))

p_pub_heat <- ggplot(pub_heat_tbl, aes(CagePhase, Outcome, fill = cohen_d)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = SigLabel), size = 2.4, colour = "black") +
  facet_grid(Metric + Sex ~ BiologicalComparison, scales = "free_x", space = "free_x") +
  labs(
    title = "Temporal instability effect-size map across stress-model comparisons",
    subtitle = "Cohen's d for pairwise contrasts; columns are cage change and phase, asterisk marks BH-adjusted FDR < 0.05",
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
  p_pub_heat,
  file.path(publication_overview_dir, "effect_size_heatmap_all_metrics"),
  width = 190,
  height = 160
)

pub_delta_heat_tbl <- delta_stats$contrasts %>%
  filter(contrast %in% pairwise_contrasts, status == "tested") %>%
  left_join(
    stress_timing_tbl %>%
      select(CageChange, CageChangeOrder, StressModelWindow),
    by = "CageChange"
  ) %>%
  left_join(stress_comparison_tbl, by = "contrast") %>%
  mutate(
    Outcome = factor(Outcome, levels = paste0("delta_", instability_outcomes), labels = metric_labels[paste0("delta_", instability_outcomes)]),
    Metric = factor(as.character(Metric), levels = metrics_to_analyze),
    CagePhase = paste(as.character(CageChange), Phase, sep = "\n"),
    CagePhase = factor(CagePhase, levels = unique(CagePhase[order(CageChangeOrder, Phase)])),
    contrast = factor(contrast, levels = pairwise_contrasts),
    BiologicalComparison = factor(BiologicalComparison, levels = stress_comparison_tbl$BiologicalComparison),
    SigLabel = if_else(ReportingP < 0.05, "*", "")
  )

write_table(pub_delta_heat_tbl, file.path(output_dir, "tables", "publication_delta_effect_size_heatmap_data.csv"))

p_pub_delta_heat <- ggplot(pub_delta_heat_tbl, aes(CagePhase, Outcome, fill = cohen_d)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = SigLabel), size = 2.4, colour = "black") +
  facet_grid(Metric + Sex ~ BiologicalComparison, scales = "free_x", space = "free_x") +
  labs(
    title = "Adaptation effect-size map across stress-model comparisons",
    subtitle = "Cohen's d for delta metrics from first cage change; asterisk marks BH-adjusted FDR < 0.05",
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
  p_pub_delta_heat,
  file.path(publication_overview_dir, "delta_effect_size_heatmap_all_metrics"),
  width = 190,
  height = 160
)

message("Temporal instability analysis complete. Primary readout: Movement. Secondary exploratory readouts: Entropy and Proximity.")
