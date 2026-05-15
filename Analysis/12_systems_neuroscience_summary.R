# ================================================================
# Integrated Systems Neuroscience Summary
# MMMSociability
# ================================================================
# Goal:
#   Compress high-dimensional RFID/home-cage behavioral output into a
#   publication-facing systems neuroscience summary for stress resilience
#   and susceptibility in male and female mice.
#
# Conceptual framing:
#   This script does not replace the specialized scripts. It integrates them.
#   It is intended as the final graphical/statistical layer after running:
#     03_build_multiscale_behavior_metrics.R
#     06_burstiness_temporal_instability.R
#     07_behavioral_state_space.R
#     08_early_prediction_models.R
#     09_dynamic_social_networks.R
#     10_hmm_behavioral_states.R
#     11_gamm_trajectory_features.R
#
# Output:
#   - animal-level multiscale systems feature matrix
#   - feature dictionary and QC tables
#   - CON/RES/SUS group summaries and contrasts by sex
#   - PCA/UMAP behavioral state-space plots
#   - systems effect-size heatmaps
#   - sex-specific feature correlation networks
#   - early-prediction summary panels if outcome data are available
#   - a compact publication dashboard figure
#
# Design principles:
#   - one animal = one row in the integrated feature matrix
#   - all features are explicitly tagged by domain, time scale, and source
#   - statistics emphasize effect sizes, uncertainty, and FDR control
#   - plots are SVG/PDF export-ready and minimally styled
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

# ------------------------------------------------
# USER CONFIGURATION
# ------------------------------------------------

# Project root used in your existing scripts. Change only if needed.
project_root <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID"
repo_root <- "C:/Users/topohl/Documents/GitHub/MMMSociability"

# Main scale for the publication summary. Use 5min_based for the ACF/entropy story;
# use 10min_based for lower noise and prediction sensitivity.
primary_bin_level <- "5min_based"
sensitivity_bin_levels <- c("10sec_based", "5min_based", "10min_based", "30min_based")

# Main output root.
output_dir <- file.path(project_root, "analysis_ready/12_systems_neuroscience_summary", primary_bin_level)

# Optional endpoint file for physiology/behavioral burden/proteomics module data.
# Expected: one row per animal, with an AnimalNum-like column and endpoint columns.
# If NULL, the script tries to use endpoint columns already present in the integrated data.
endpoint_file <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/SIS_Analysis/E9_Behavior_Data.xlsx"
endpoint_sheet <- "zScore"
endpoint_cols <- c(
  "CombZ", "stress_z_score", "NOR", "sucrose_pref", "weight_dev", "delta_cort", "adrenal_weight", "spleen_weight",
  "SucrosePreference", "Corticosterone", "DeltaCorticosterone"
)
primary_outcome <- "CombZ"
primary_outcome_label <- "CombZ (lower = worse depressive-like endpoint; higher = more resilient-like endpoint)"
primary_outcome_direction <- "lower_worse"

# RES/SUS labels were derived from post-paradigm CombZ. Keep group contrasts
# descriptive, and restrict prediction claims to pre-endpoint behavioral features.
res_sus_derived_from_outcome <- TRUE
prospective_prediction_context <- "first_active_12h"

# Optional proteomics module score file. Same expectation: one row per animal.
# Useful columns could include RNP_module, Mito_module, Translation_module, etc.
proteomics_module_file <- NULL

# If TRUE, limits some plots to features that are easier to explain in a paper.
publication_focus_only <- TRUE

# Minimum animals per group for group contrasts.
min_n_per_group <- 2

# Plot colors. CON/RES/SUS are ordered deliberately.
group_levels <- c("CON", "RES", "SUS")
group_colors <- c("CON" = "#3d3b6e", "RES" = "#C6C3BB", "SUS" = "#e63947")
sex_levels <- c("Female", "Male")

# ------------------------------------------------
# SOURCE HELPERS
# ------------------------------------------------

source_if_exists <- function(path) {
  if (file.exists(path)) source(path)
}

source_if_exists(file.path(repo_root, "Functions/behavioral_dynamics_helpers.R"))
source_if_exists(file.path(repo_root, "Functions/behavioral_dynamics_stats_helpers.R"))

if (!exists("ensure_dir")) {
  ensure_dir <- function(path) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    invisible(path)
  }
}

if (!exists("write_table")) {
  write_table <- function(x, path) {
    ensure_dir(dirname(path))
    readr::write_csv(x, path)
    invisible(path)
  }
}

if (!exists("safe_name")) {
  safe_name <- function(x) {
    x %>% as.character() %>% str_replace_all("[^A-Za-z0-9]+", "_") %>% str_replace_all("^_|_$", "") %>% str_to_lower()
  }
}

if (TRUE) {
  make_nature_theme <- function(base_size = 7) {
    theme_classic(base_size = base_size, base_family = "Arial") +
      theme(
        axis.line = element_line(linewidth = 0.28, colour = "black"),
        axis.ticks = element_line(linewidth = 0.24, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.title = element_text(colour = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", colour = "black"),
        legend.title = element_blank(),
        legend.position = "top",
        legend.key.height = unit(3.5, "mm"),
        legend.key.width = unit(6, "mm"),
        plot.title = element_text(face = "bold", hjust = 0, margin = margin(b = 2)),
        plot.subtitle = element_text(hjust = 0, colour = "grey25", margin = margin(b = 3)),
        plot.caption = element_text(hjust = 0, colour = "grey30", size = rel(0.85)),
        plot.margin = margin(4, 4, 4, 4),
        panel.spacing = unit(1.1, "lines")
      )
  }
}

if (TRUE) {
  save_plot_svg_pdf <- function(plot, filename_base, width = 85, height = 65, units = "mm") {
    ensure_dir(dirname(filename_base))
    ggplot2::ggsave(paste0(filename_base, ".svg"), plot, width = width, height = height, units = units)
    pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else "pdf"
    ggplot2::ggsave(paste0(filename_base, ".pdf"), plot, width = width, height = height, units = units, device = pdf_device)
    ggplot2::ggsave(paste0(filename_base, ".png"), plot, width = width, height = height, units = units, dpi = 600)
    invisible(filename_base)
  }
}

if (!exists("first_existing_col")) {
  first_existing_col <- function(dat, candidates, required = TRUE, label = "column") {
    hit <- candidates[candidates %in% names(dat)][1]
    if (is.na(hit) && required) {
      stop("Could not find ", label, ". Tried: ", paste(candidates, collapse = ", "), call. = FALSE)
    }
    hit
  }
}

if (!exists("cohens_d_pooled")) {
  cohens_d_pooled <- function(x, y) {
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
    if (length(x) < 2 || length(y) < 2) return(NA_real_)
    sp <- sqrt(((length(x) - 1) * var(x) + (length(y) - 1) * var(y)) / (length(x) + length(y) - 2))
    if (!is.finite(sp) || sp == 0) return(NA_real_)
    (mean(y) - mean(x)) / sp
  }
}

# ------------------------------------------------
# GENERAL HELPERS
# ------------------------------------------------

ensure_dir(output_dir)
ensure_dir(file.path(output_dir, "tables"))
ensure_dir(file.path(output_dir, "stats_tables"))
ensure_dir(file.path(output_dir, "figures"))
ensure_dir(file.path(output_dir, "figures/publication_panels"))
ensure_dir(file.path(output_dir, "figures/qc"))

read_any_table <- function(path, sheet = NULL) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") return(readr::read_csv(path, show_col_types = FALSE))
  if (ext %in% c("tsv", "txt")) return(readr::read_tsv(path, show_col_types = FALSE))
  if (ext == "rds") return(readRDS(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("Install readxl to read Excel files: ", path)
    if (is.null(sheet)) return(readxl::read_excel(path))
    return(readxl::read_excel(path, sheet = sheet))
  }
  NULL
}

first_existing_path <- function(candidates) {
  candidates <- candidates[!is.na(candidates)]
  hit <- candidates[file.exists(candidates)][1]
  if (length(hit) == 0 || is.na(hit)) NA_character_ else hit
}

parse_cage_change_index <- function(x) {
  idx <- suppressWarnings(as.integer(str_extract(as.character(x), "\\d+")))
  fallback <- match(as.character(x), sort(unique(as.character(x))))
  ifelse(is.finite(idx), idx, fallback)
}

get_first_cage_change <- function(x) {
  ux <- unique(as.character(x))
  idx <- parse_cage_change_index(ux)
  ux[which.min(ifelse(is.finite(idx), idx, Inf))]
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

safe_scale <- function(x) {
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

safe_cor_p <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE)$p.value)
}

sem <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  sd(x) / sqrt(length(x))
}

mean_ci <- function(x, conf = 0.95) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(c(low = NA_real_, high = NA_real_))
  se <- sd(x) / sqrt(length(x))
  crit <- qt(1 - (1 - conf) / 2, df = length(x) - 1)
  c(low = mean(x) - crit * se, high = mean(x) + crit * se)
}

hedges_g <- function(x, y) {
  d <- cohens_d_pooled(x, y)
  n_x <- sum(is.finite(x))
  n_y <- sum(is.finite(y))
  df <- n_x + n_y - 2
  if (!is.finite(d) || df <= 1) return(NA_real_)
  correction <- 1 - 3 / (4 * df - 1)
  d * correction
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

p_label <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "p<0.001",
    TRUE ~ paste0("p=", formatC(p, format = "f", digits = 3))
  )
}

q_label <- function(q) {
  case_when(
    is.na(q) ~ "q=NA",
    q < 0.001 ~ "q<0.001",
    TRUE ~ paste0("q=", formatC(q, format = "f", digits = 3))
  )
}

rho_label <- function(rho) {
  ifelse(is.finite(rho), paste0("rho=", formatC(rho, format = "f", digits = 2)), "rho=NA")
}

sig_label <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    p < 0.10 ~ ".",
    TRUE ~ ""
  )
}

standardize_id_columns <- function(dat) {
  normalize_animal_id <- function(x) {
    x %>%
      as.character() %>%
      str_trim() %>%
      str_replace_all("\\s+", "") %>%
      str_to_upper()
  }

  animal_col <- first_existing_col(dat, c("AnimalNum", "Animal", "AnimalID", "MouseID", "Mouse", "ID", "RFID", "animal_id"), TRUE, "animal ID")
  group_col <- first_existing_col(dat, c("Group", "Phenotype", "Condition", "Treatment", "StressGroup"), FALSE, "group")
  sex_col <- first_existing_col(dat, c("Sex", "sex"), FALSE, "sex")

  out <- dat %>% mutate(AnimalNum = normalize_animal_id(.data[[animal_col]]))
  out$Group <- if (!is.na(group_col)) as.character(dat[[group_col]]) else NA_character_
  out$Sex <- if (!is.na(sex_col)) as.character(dat[[sex_col]]) else NA_character_
  out %>%
    mutate(
      Group = factor(Group, levels = unique(c(group_levels, sort(unique(na.omit(Group)))))),
      Sex = factor(Sex, levels = unique(c(sex_levels, sort(unique(na.omit(Sex))))))
    )
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

parse_system_feature <- function(feature) {
  parts <- str_split_fixed(feature, "__", 6)
  tibble(
    feature = feature,
    Source = parts[, 1],
    Domain = parts[, 2],
    Scale = parts[, 3],
    Metric = parts[, 4],
    Statistic = parts[, 5],
    Context = parts[, 6]
  )
}

# ------------------------------------------------
# INPUT PATH REGISTRY
# ------------------------------------------------

paths <- tibble(
  Source = c(
    "multiscale_metrics",
    "burstiness_instability",
    "state_space",
    "early_prediction",
    "dynamic_networks",
    "hmm_states",
    "gamm_trajectory"
  ),
  Path = c(
    file.path(project_root, "analysis_ready/03_derived_metrics", primary_bin_level, "all_behavior_metrics.csv"),
    file.path(project_root, "analysis_ready/06_behavioral_dynamics/burstiness", primary_bin_level, "tables"),
    file.path(project_root, "analysis_ready/06_behavioral_dynamics/state_space", primary_bin_level, "tables"),
    file.path(project_root, "analysis_ready/06_behavioral_dynamics/early_prediction", primary_bin_level, "tables"),
    file.path(project_root, "analysis_ready/06_behavioral_dynamics/dynamic_social_networks", primary_bin_level, "tables"),
    file.path(project_root, "analysis_ready/06_behavioral_dynamics/hmm_states", primary_bin_level, "tables"),
    file.path(project_root, "analysis_ready/06_behavioral_dynamics/gamm_trajectory_features", primary_bin_level, "tables")
  )
)

write_table(paths %>% mutate(Exists = file.exists(Path)), file.path(output_dir, "tables/input_path_registry.csv"))

# ------------------------------------------------
# LOAD BASE MULTISCALE METRICS
# ------------------------------------------------

base_file <- paths$Path[paths$Source == "multiscale_metrics"]
base_raw <- read_any_table(base_file)
if (is.null(base_raw)) {
  stop("Missing base metrics file. Run Analysis/03_build_multiscale_behavior_metrics.R first. Expected: ", base_file, call. = FALSE)
}

base <- standardize_id_columns(base_raw)

phase_col <- first_existing_col(base, c("Phase", "phase"), FALSE, "phase")
cage_col <- first_existing_col(base, c("CageChange", "CC", "CageChangeNum", "Regrouping"), FALSE, "cage-change")
time_col <- first_existing_col(base, c("TimeIndex", "BinStart", "HalfHourElapsed", "HalfHourWithinCC0", "Time", "DateTime"), TRUE, "time")

base <- base %>%
  mutate(
    Phase = if (!is.na(phase_col)) as.character(.data[[phase_col]]) else "All",
    CageChange = if (!is.na(cage_col)) as.character(.data[[cage_col]]) else "All",
    CageChangeIndex = parse_cage_change_index(CageChange),
    PhaseClass = case_when(
      str_detect(str_to_lower(Phase), "active|dark|night") ~ "Active",
      str_detect(str_to_lower(Phase), "inactive|light|day") ~ "Inactive",
      TRUE ~ Phase
    ),
    TimeIndex = .data[[time_col]],
    Movement = safe_numeric(.data[[first_existing_col(., c("Movement", "movement", "Distance", "Activity"), TRUE, "movement")]]),
    Entropy = safe_numeric(.data[[first_existing_col(., c("Entropy", "entropy", "ShannonEntropy", "PositionEntropy"), TRUE, "entropy")]]),
    Proximity = safe_numeric(.data[[first_existing_col(., c("ProximityFraction", "Proximity", "proximity", "MeanProximity"), TRUE, "proximity")]])
  ) %>%
  arrange(AnimalNum, CageChange, Phase, TimeIndex)

# ------------------------------------------------
# BUILD CORE ANIMAL-LEVEL FEATURE MATRIX
# ------------------------------------------------

calc_feature_set <- function(dat, context_name) {
  metric_cols <- intersect(c("Movement", "Entropy", "Proximity"), names(dat))
  if (length(metric_cols) == 0 || nrow(dat) == 0) return(tibble())

  dat %>%
    select(AnimalNum, Group, Sex, all_of(metric_cols)) %>%
    pivot_longer(all_of(metric_cols), names_to = "Metric", values_to = "Value") %>%
    group_by(AnimalNum, Group, Sex, Metric) %>%
    summarise(
      mean = mean(Value, na.rm = TRUE),
      sd = sd(Value, na.rm = TRUE),
      median = median(Value, na.rm = TRUE),
      q25 = if (all(!is.finite(Value))) NA_real_ else quantile(Value, 0.25, na.rm = TRUE, names = FALSE),
      q75 = if (all(!is.finite(Value))) NA_real_ else quantile(Value, 0.75, na.rm = TRUE, names = FALSE),
      p95 = if (all(!is.finite(Value))) NA_real_ else quantile(Value, 0.95, na.rm = TRUE, names = FALSE),
      max = if (all(!is.finite(Value))) NA_real_ else max(Value, na.rm = TRUE),
      rmssd = if (exists("calc_rmssd")) calc_rmssd(Value) else sqrt(mean(diff(Value[is.finite(Value)])^2, na.rm = TRUE)),
      acf1 = if (exists("calc_acf1")) calc_acf1(Value) else safe_cor(Value[-length(Value)], Value[-1], "pearson"),
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(mean, sd, median, q25, q75, p95, max, rmssd, acf1), names_to = "Statistic", values_to = "FeatureValue") %>%
    mutate(
      feature = pmap_chr(
        list(Source = "raw", Domain = "behavior", Scale = primary_bin_level, Metric = Metric, Statistic = Statistic),
        ~ make_feature_name(..1, ..2, ..3, ..4, ..5, context_name)
      )
    ) %>%
    select(AnimalNum, Group, Sex, feature, FeatureValue)
}

base_global_features <- calc_feature_set(base, "all")
base_phase_features <- base %>%
  group_split(Phase, .keep = TRUE) %>%
  map_dfr(~ calc_feature_set(.x, paste0("phase_", unique(.x$Phase)[1])))

first_cage_change <- get_first_cage_change(base$CageChange)
early_window_bins <- case_when(
  primary_bin_level == "10sec_based" ~ 12 * 60 * 6,
  primary_bin_level == "5min_based" ~ 12 * 60 / 5,
  primary_bin_level == "10min_based" ~ 12 * 60 / 10,
  primary_bin_level == "30min_based" ~ 12 * 60 / 30,
  TRUE ~ 144
)

first_active <- base %>%
  filter(
    as.character(CageChange) == first_cage_change,
    str_detect(str_to_lower(Phase), "active|dark|night")
  ) %>%
  group_by(AnimalNum, Phase) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  mutate(local_bin = row_number()) %>%
  filter(local_bin <= early_window_bins) %>%
  ungroup()

base_first_active_features <- calc_feature_set(first_active, prospective_prediction_context)

core_feature_long <- bind_rows(base_global_features, base_phase_features, base_first_active_features) %>%
  filter(is.finite(FeatureValue))

core_feature_wide <- core_feature_long %>%
  group_by(AnimalNum, Group, Sex, feature) %>%
  summarise(FeatureValue = mean(FeatureValue, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = feature, values_from = FeatureValue)

scaled_feature <- function(dat, feature_name) {
  if (!feature_name %in% names(dat)) return(rep(NA_real_, nrow(dat)))
  safe_scale(safe_numeric(dat[[feature_name]]))
}

# Add interpretable composite axes from early active window.
early_movement_mean <- make_feature_name("raw", "behavior", primary_bin_level, "Movement", "mean", prospective_prediction_context)
early_proximity_mean <- make_feature_name("raw", "behavior", primary_bin_level, "Proximity", "mean", prospective_prediction_context)
early_movement_rmssd <- make_feature_name("raw", "behavior", primary_bin_level, "Movement", "rmssd", prospective_prediction_context)
early_entropy_rmssd <- make_feature_name("raw", "behavior", primary_bin_level, "Entropy", "rmssd", prospective_prediction_context)
early_proximity_rmssd <- make_feature_name("raw", "behavior", primary_bin_level, "Proximity", "rmssd", prospective_prediction_context)
early_movement_acf1 <- make_feature_name("raw", "behavior", primary_bin_level, "Movement", "acf1", prospective_prediction_context)
early_entropy_acf1 <- make_feature_name("raw", "behavior", primary_bin_level, "Entropy", "acf1", prospective_prediction_context)
early_proximity_acf1 <- make_feature_name("raw", "behavior", primary_bin_level, "Proximity", "acf1", prospective_prediction_context)

core_feature_wide <- core_feature_wide %>%
  mutate(
    systems__composite__early__social_withdrawal__z__first_active_12h =
      scaled_feature(pick(everything()), early_movement_mean) -
      scaled_feature(pick(everything()), early_proximity_mean),
    systems__composite__early__behavioral_instability__z__first_active_12h =
      rowMeans(cbind(
        scaled_feature(pick(everything()), early_movement_rmssd),
        scaled_feature(pick(everything()), early_entropy_rmssd),
        scaled_feature(pick(everything()), early_proximity_rmssd)
      ), na.rm = TRUE),
    systems__composite__early__behavioral_inertia__z__first_active_12h =
      rowMeans(cbind(
        scaled_feature(pick(everything()), early_movement_acf1),
        scaled_feature(pick(everything()), early_entropy_acf1),
        scaled_feature(pick(everything()), early_proximity_acf1)
      ), na.rm = TRUE)
  )

# ------------------------------------------------
# OPTIONAL MODULE OUTPUTS FROM SPECIALIZED SCRIPTS
# ------------------------------------------------

load_optional_animal_table <- function(path, source_label, domain_label, scale_label = primary_bin_level) {
  dat <- read_any_table(path)
  if (is.null(dat) || nrow(dat) == 0) return(tibble())
  dat <- standardize_id_columns(dat)

  scale_vec <- if ("BinLevel" %in% names(dat)) as.character(dat$BinLevel) else rep(scale_label, nrow(dat))

  descriptor_cols <- intersect(
    c("Metric", "StateLabel", "State", "Phase", "CageChange", "CageChangeIndex"),
    names(dat)
  )

  numeric_cols <- dat %>%
    select(where(is.numeric)) %>%
    names()
  numeric_cols <- setdiff(
    numeric_cols,
    c(
      "AnimalNum", "TimeIndex", "BinSizeSec", "CageChangeIndex", "State",
      "n", "n_bins", "n_bins_state", "n_transitions"
    )
  )
  if (length(numeric_cols) == 0) return(tibble())

  out <- dat %>%
    select(AnimalNum, Group, Sex, any_of(descriptor_cols), all_of(numeric_cols)) %>%
    mutate(ScaleName = scale_vec) %>%
    mutate(across(any_of(descriptor_cols), ~ safe_name(as.character(.x)))) %>%
    pivot_longer(all_of(numeric_cols), names_to = "MetricName", values_to = "FeatureValue")

  if (length(descriptor_cols) > 0) {
    out <- out %>%
      tidyr::unite("ContextTag", all_of(descriptor_cols), sep = "_", remove = FALSE, na.rm = TRUE) %>%
      group_by(AnimalNum, Group, Sex, ScaleName, MetricName, ContextTag) %>%
      summarise(FeatureValue = mean(FeatureValue, na.rm = TRUE), .groups = "drop")
  } else {
    out <- out %>%
      mutate(ContextTag = "script_output") %>%
      group_by(AnimalNum, Group, Sex, ScaleName, MetricName, ContextTag) %>%
      summarise(FeatureValue = mean(FeatureValue, na.rm = TRUE), .groups = "drop")
  }

  out %>%
    filter(is.finite(FeatureValue)) %>%
    mutate(feature = make_feature_name(source_label, domain_label, ScaleName, MetricName, "mean", ContextTag)) %>%
    select(AnimalNum, Group, Sex, feature, FeatureValue)
}

optional_files <- tibble(
  source_label = c("burstiness", "state_space", "state_space", "early_prediction", "dynamic_network", "hmm", "gamm"),
  domain_label = c("temporal_dynamics", "latent_space", "behavioral_state", "prediction", "social_network", "behavioral_state", "trajectory"),
  path = c(
    first_existing_path(file.path(project_root, "analysis_ready/06_behavioral_dynamics/burstiness", sensitivity_bin_levels, "tables/burstiness_instability_features.csv")),
    first_existing_path(file.path(project_root, "analysis_ready/06_behavioral_dynamics/state_space", sensitivity_bin_levels, "tables/state_diversity_metrics.csv")),
    first_existing_path(file.path(project_root, "analysis_ready/06_behavioral_dynamics/state_space", sensitivity_bin_levels, "tables/state_switching_metrics.csv")),
    first_existing_path(file.path(project_root, "analysis_ready/06_behavioral_dynamics/early_prediction", sensitivity_bin_levels, "tables/early_behavior_features.csv")),
    first_existing_path(file.path(project_root, "analysis_ready/06_behavioral_dynamics/social_networks", sensitivity_bin_levels, "tables/animal_level_social_dynamics.csv")),
    first_existing_path(file.path(project_root, "analysis_ready/06_behavioral_dynamics/hmm_states", sensitivity_bin_levels, "tables/hmm_state_occupancy.csv")),
    first_existing_path(c(
      file.path(project_root, "analysis_ready/06_behavioral_dynamics/gamm_trajectory_features", sensitivity_bin_levels, "tables/gamm_trajectory_features.csv"),
      file.path(project_root, "analysis_ready/06_behavioral_dynamics/gamm_features", sensitivity_bin_levels, "tables/combined_gamm_features.csv")
    ))
  )
)

optional_long <- pmap_dfr(optional_files, function(source_label, domain_label, path) {
  load_optional_animal_table(path, source_label, domain_label)
})

optional_wide <- if (nrow(optional_long) > 0) {
  optional_long %>%
    group_by(AnimalNum, Group, Sex, feature) %>%
    summarise(FeatureValue = mean(FeatureValue, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = feature, values_from = FeatureValue)
} else {
  tibble(AnimalNum = character(), Group = factor(levels = group_levels), Sex = factor(levels = sex_levels))
}

systems_features <- full_join(core_feature_wide, optional_wide, by = c("AnimalNum", "Group", "Sex")) %>%
  mutate(
    Group = factor(as.character(Group), levels = group_levels),
    Sex = factor(as.character(Sex), levels = unique(c(sex_levels, sort(unique(as.character(Sex))))))
  )

# ------------------------------------------------
# OPTIONAL ENDPOINTS AND PROTEOMICS MODULES
# ------------------------------------------------

endpoint_dat <- read_any_table(endpoint_file, sheet = endpoint_sheet)
if (!is.null(endpoint_dat)) {
  endpoint_dat <- standardize_id_columns(endpoint_dat)
  endpoint_keep <- intersect(endpoint_cols, names(endpoint_dat))
  endpoint_dat <- endpoint_dat %>% select(AnimalNum, all_of(endpoint_keep))
  systems_features <- left_join(systems_features, endpoint_dat, by = "AnimalNum")
}

proteomics_dat <- read_any_table(proteomics_module_file)
if (!is.null(proteomics_dat)) {
  proteomics_dat <- standardize_id_columns(proteomics_dat)
  proteomics_numeric <- proteomics_dat %>% select(where(is.numeric)) %>% names()
  proteomics_numeric <- setdiff(proteomics_numeric, "AnimalNum")
  proteomics_dat <- proteomics_dat %>%
    select(AnimalNum, all_of(proteomics_numeric)) %>%
    rename_with(~ paste0("proteomics_module__", safe_name(.x)), all_of(proteomics_numeric))
  systems_features <- left_join(systems_features, proteomics_dat, by = "AnimalNum")
}

feature_cols <- systems_features %>% select(where(is.numeric)) %>% names()
feature_cols <- setdiff(feature_cols, c("AnimalNum", endpoint_cols))

feature_dictionary <- parse_system_feature(feature_cols) %>%
  mutate(
    Interpretation = case_when(
      str_detect(feature, "movement") & str_detect(feature, "mean") ~ "Mean locomotor output",
      str_detect(feature, "entropy") & str_detect(feature, "mean") ~ "Mean spatial dispersion / positional entropy",
      str_detect(feature, "proximity") & str_detect(feature, "mean") ~ "Mean social proximity",
      str_detect(feature, "rmssd") ~ "Short-timescale temporal instability",
      str_detect(feature, "acf1") ~ "Temporal inertia / persistence",
      str_detect(feature, "social_withdrawal") ~ "High movement with low proximity; exploratory-social decoupling",
      str_detect(feature, "behavioral_instability") ~ "Composite local volatility across movement, entropy and proximity",
      str_detect(feature, "behavioral_inertia") ~ "Composite autocorrelation/persistence across behavioral channels",
      TRUE ~ "Integrated behavioral feature"
    )
  )

feature_qc <- map_dfr(feature_cols, function(fc) {
  x <- systems_features[[fc]]
  tibble(
    feature = fc,
    n_animals = nrow(systems_features),
    n_finite = sum(is.finite(x)),
    missing_fraction = mean(!is.finite(x)),
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    zero_variance = sd(x, na.rm = TRUE) == 0
  )
}) %>%
  left_join(feature_dictionary, by = "feature") %>%
  arrange(Source, Domain, Scale, Metric, Statistic)

usable_features <- feature_qc %>%
  filter(n_finite >= 4, missing_fraction <= 0.50, !zero_variance) %>%
  pull(feature)

prospective_prediction_features <- feature_qc %>%
  filter(
    feature %in% usable_features,
    (
      Source == "raw" & Context == prospective_prediction_context
    ) |
      (
        Source == "systems" & Context == prospective_prediction_context
      )
  ) %>%
  pull(feature)

feature_sets <- tibble(
  FeatureSet = c("descriptive_full_experiment", "prospective_prediction"),
  n_features = c(length(usable_features), length(prospective_prediction_features)),
  IntendedUse = c(
    "Descriptive phenotype architecture after RES/SUS labels are known",
    "Prediction/association with post-paradigm endpoints using pre-endpoint early behavior only"
  ),
  IncludesCombZDerivedGroups = c(res_sus_derived_from_outcome, FALSE),
  Notes = c(
    "RES/SUS group contrasts are descriptive because groups were derived from CombZ.",
    paste0("Uses raw first active 12 h after ", first_cage_change, " plus early composites; excludes GAMM/HMM/full-experiment and mixed-phase early-prediction module features.")
  )
)

write_table(systems_features, file.path(output_dir, "tables/systems_animal_feature_matrix.csv"))
write_table(feature_dictionary, file.path(output_dir, "tables/systems_feature_dictionary.csv"))
write_table(feature_qc, file.path(output_dir, "tables/systems_feature_qc.csv"))
write_table(feature_sets, file.path(output_dir, "tables/systems_feature_sets.csv"))

# ------------------------------------------------
# GROUP SUMMARY AND CONTRASTS
# ------------------------------------------------

systems_long <- systems_features %>%
  select(AnimalNum, Group, Sex, all_of(usable_features)) %>%
  pivot_longer(all_of(usable_features), names_to = "feature", values_to = "Value") %>%
  filter(is.finite(Value)) %>%
  left_join(feature_dictionary, by = "feature")

group_summary <- systems_long %>%
  group_by(Sex, Group, feature, Source, Domain, Scale, Metric, Statistic, Context) %>%
  summarise(
    n = n_distinct(AnimalNum),
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    sem = sem(Value),
    ci95_low = mean_ci(Value)["low"],
    ci95_high = mean_ci(Value)["high"],
    median = median(Value, na.rm = TRUE),
    q25 = quantile(Value, 0.25, na.rm = TRUE, names = FALSE),
    q75 = quantile(Value, 0.75, na.rm = TRUE, names = FALSE),
    iqr = IQR(Value, na.rm = TRUE),
    ReportingN = paste0("n=", n, " animals"),
    SummaryStatistic = "mean +/- 95% CI; median and IQR also reported",
    .groups = "drop"
  )

contrast_pairs <- list(c("CON", "RES"), c("CON", "SUS"), c("RES", "SUS"))

group_contrasts <- systems_long %>%
  filter(!is.na(Group), as.character(Group) %in% group_levels) %>%
  group_by(Sex, feature, Source, Domain, Scale, Metric, Statistic, Context) %>%
  group_modify(~{
    d <- .x
    map_dfr(contrast_pairs, function(pair) {
      ref <- pair[1]
      comp <- pair[2]
      x <- d$Value[as.character(d$Group) == ref]
      y <- d$Value[as.character(d$Group) == comp]
      n_ref <- sum(is.finite(x))
      n_comp <- sum(is.finite(y))
      if (n_ref < min_n_per_group || n_comp < min_n_per_group) {
        return(tibble(contrast = paste0(comp, "-", ref), n_ref = n_ref, n_comp = n_comp,
                      mean_ref = NA_real_, mean_comp = NA_real_, median_ref = NA_real_, median_comp = NA_real_,
                      estimate = NA_real_, estimate_ci_low = NA_real_, estimate_ci_high = NA_real_,
                      p.value = NA_real_, wilcox_p = NA_real_, cohen_d = NA_real_, hedges_g = NA_real_,
                      test_method = "Welch two-sample t-test; Wilcoxon rank-sum also reported",
                      status = "low_n"))
      }
      tt <- try(t.test(y, x), silent = TRUE)
      wt <- try(wilcox.test(y, x, exact = FALSE), silent = TRUE)
      ci <- if (inherits(tt, "try-error") || is.null(tt$conf.int)) c(NA_real_, NA_real_) else as.numeric(tt$conf.int)
      tibble(
        contrast = paste0(comp, "-", ref),
        n_ref = n_ref,
        n_comp = n_comp,
        mean_ref = mean(x, na.rm = TRUE),
        mean_comp = mean(y, na.rm = TRUE),
        median_ref = median(x, na.rm = TRUE),
        median_comp = median(y, na.rm = TRUE),
        estimate = mean(y, na.rm = TRUE) - mean(x, na.rm = TRUE),
        estimate_ci_low = ci[1],
        estimate_ci_high = ci[2],
        p.value = if (inherits(tt, "try-error")) NA_real_ else tt$p.value,
        wilcox_p = if (inherits(wt, "try-error")) NA_real_ else wt$p.value,
        cohen_d = cohens_d_pooled(x, y),
        hedges_g = hedges_g(x, y),
        test_method = "Welch two-sample t-test; Wilcoxon rank-sum also reported",
        status = "tested"
      )
    })
  }) %>%
  ungroup() %>%
  group_by(Sex, contrast) %>%
  mutate(p_fdr_sex_contrast = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  group_by(Sex, Source, Domain, Scale, Context, contrast) %>%
  mutate(
    p_fdr_family = p.adjust(p.value, method = "BH"),
    p_fdr = p_fdr_family,
    sig = sig_label(p_fdr),
    ReportingCorrection = paste0("BH FDR within Sex x Source x Domain x Scale x Context x contrast; broad Sex x contrast FDR also provided"),
    evidence = case_when(
      is.na(p_fdr) ~ "not_tested",
      p_fdr < 0.05 & abs(hedges_g) >= 0.8 ~ "large_FDR_supported",
      p_fdr < 0.05 ~ "FDR_supported",
      p_fdr < 0.10 ~ "trend",
      abs(hedges_g) >= 0.8 ~ "large_effect_uncertain",
      TRUE ~ "weak_or_no_evidence"
    )
  ) %>%
  ungroup() %>%
  arrange(Sex, contrast, p_fdr, desc(abs(cohen_d)))

write_table(group_summary, file.path(output_dir, "stats_tables/systems_group_summary.csv"))
write_table(group_contrasts, file.path(output_dir, "stats_tables/systems_group_contrasts.csv"))

# ------------------------------------------------
# DIMENSIONALITY REDUCTION: PCA + OPTIONAL UMAP
# ------------------------------------------------

x <- systems_features %>% select(all_of(usable_features))
x_imp <- x %>% mutate(across(everything(), ~ replace_na(.x, median(.x, na.rm = TRUE))))
x_mat <- as.matrix(x_imp)
x_scaled <- scale(x_mat)
bad_scaled_cols <- !is.finite(colSums(x_scaled))
if (any(bad_scaled_cols)) x_scaled[, bad_scaled_cols] <- 0

pca <- prcomp(x_scaled, center = FALSE, scale. = FALSE)
var_exp <- pca$sdev^2 / sum(pca$sdev^2)

pca_scores <- systems_features %>%
  select(AnimalNum, Group, Sex, any_of(endpoint_cols)) %>%
  bind_cols(as_tibble(pca$x[, seq_len(min(6, ncol(pca$x))), drop = FALSE]))

pca_loadings <- as_tibble(pca$rotation[, seq_len(min(6, ncol(pca$rotation))), drop = FALSE], rownames = "feature") %>%
  left_join(feature_dictionary, by = "feature")

pca_variance <- tibble(
  PC = paste0("PC", seq_along(var_exp)),
  variance_explained = var_exp,
  cumulative_variance_explained = cumsum(var_exp)
)

write_table(pca_scores, file.path(output_dir, "tables/systems_pca_scores.csv"))
write_table(pca_loadings, file.path(output_dir, "tables/systems_pca_loadings.csv"))
write_table(pca_variance, file.path(output_dir, "tables/systems_pca_variance.csv"))

pca_plot_tbl <- pca_scores %>%
  add_count(Sex, name = "SexN") %>%
  mutate(SexPanel = paste0(Sex, " (n=", SexN, ")"))

p_pca <- pca_plot_tbl %>%
  ggplot(aes(PC1, PC2, colour = Group, fill = Group)) +
  geom_hline(yintercept = 0, linewidth = 0.15, colour = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.15, colour = "grey85") +
  stat_ellipse(aes(group = Group), linewidth = 0.35, alpha = 0.45, show.legend = FALSE) +
  geom_point(size = 2.15, alpha = 0.85, stroke = 0.25) +
  facet_wrap(~ SexPanel, nrow = 1) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  labs(
    title = "Integrated behavioral state space",
    subtitle = paste0("PCA on ", length(usable_features), " descriptive features; n=", nrow(pca_scores), " animals"),
    x = "Systems PC1",
    y = "Systems PC2",
    caption = paste0("PC1 ", round(100 * var_exp[1], 1), "%; PC2 ", round(100 * var_exp[2], 1), "% variance")
  ) +
  make_nature_theme(base_size = 7)

save_plot_svg_pdf(p_pca, file.path(output_dir, "figures/publication_panels/Fig_systems_state_space_PCA"), width = 135, height = 75)

if (requireNamespace("uwot", quietly = TRUE) && nrow(x_scaled) >= 8) {
  set.seed(1)
  n_neighbors <- max(3, min(10, floor(nrow(x_scaled) / 2)))
  um <- uwot::umap(x_scaled, n_neighbors = n_neighbors, min_dist = 0.25, metric = "euclidean")
  umap_scores <- systems_features %>%
    select(AnimalNum, Group, Sex, any_of(endpoint_cols)) %>%
    bind_cols(tibble(UMAP1 = um[, 1], UMAP2 = um[, 2]))
  write_table(umap_scores, file.path(output_dir, "tables/systems_umap_scores.csv"))

  umap_plot_tbl <- umap_scores %>%
    add_count(Sex, name = "SexN") %>%
    mutate(SexPanel = paste0(Sex, " (n=", SexN, ")"))

  p_umap <- umap_plot_tbl %>%
    ggplot(aes(UMAP1, UMAP2, colour = Group, fill = Group)) +
    geom_point(size = 2.2, alpha = 0.88, stroke = 0.25) +
    stat_ellipse(aes(group = Group), linewidth = 0.35, alpha = 0.45, show.legend = FALSE) +
    facet_wrap(~ SexPanel, nrow = 1) +
    scale_colour_manual(values = group_colors, drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    labs(
      title = "Nonlinear behavioral manifold",
      subtitle = paste0("UMAP on integrated animal-level systems features; n=", nrow(umap_scores), "; n_neighbors=", n_neighbors),
      x = "UMAP1",
      y = "UMAP2"
    ) +
    make_nature_theme(base_size = 7)

  save_plot_svg_pdf(p_umap, file.path(output_dir, "figures/publication_panels/Fig_systems_state_space_UMAP"), width = 135, height = 75)
}

# ------------------------------------------------
# TIME-RESOLVED LATENT TRAJECTORIES AND INSTABILITY
# ------------------------------------------------

epoch_metric_summary <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(tibble(mean = NA_real_, sd = NA_real_, rmssd = NA_real_, acf1 = NA_real_))
  }
  tibble(
    mean = mean(x),
    sd = if (length(x) >= 2) sd(x) else NA_real_,
    rmssd = if (length(x) >= 3) sqrt(mean(diff(x)^2, na.rm = TRUE)) else NA_real_,
    acf1 = if (length(x) >= 4) safe_cor(x[-length(x)], x[-1], "pearson") else NA_real_
  )
}

latent_epoch_features <- base %>%
  filter(is.finite(CageChangeIndex), PhaseClass %in% c("Active", "Inactive")) %>%
  arrange(AnimalNum, CageChangeIndex, PhaseClass, TimeIndex) %>%
  group_by(AnimalNum, Group, Sex, CageChange, CageChangeIndex, PhaseClass) %>%
  summarise(
    Movement_mean = mean(Movement, na.rm = TRUE),
    Movement_sd = sd(Movement, na.rm = TRUE),
    Movement_rmssd = if (sum(is.finite(Movement)) >= 3) sqrt(mean(diff(Movement[is.finite(Movement)])^2, na.rm = TRUE)) else NA_real_,
    Movement_acf1 = if (sum(is.finite(Movement)) >= 4) safe_cor(Movement[is.finite(Movement)][-sum(is.finite(Movement))], Movement[is.finite(Movement)][-1], "pearson") else NA_real_,
    Entropy_mean = mean(Entropy, na.rm = TRUE),
    Entropy_sd = sd(Entropy, na.rm = TRUE),
    Entropy_rmssd = if (sum(is.finite(Entropy)) >= 3) sqrt(mean(diff(Entropy[is.finite(Entropy)])^2, na.rm = TRUE)) else NA_real_,
    Entropy_acf1 = if (sum(is.finite(Entropy)) >= 4) safe_cor(Entropy[is.finite(Entropy)][-sum(is.finite(Entropy))], Entropy[is.finite(Entropy)][-1], "pearson") else NA_real_,
    Proximity_mean = mean(Proximity, na.rm = TRUE),
    Proximity_sd = sd(Proximity, na.rm = TRUE),
    Proximity_rmssd = if (sum(is.finite(Proximity)) >= 3) sqrt(mean(diff(Proximity[is.finite(Proximity)])^2, na.rm = TRUE)) else NA_real_,
    Proximity_acf1 = if (sum(is.finite(Proximity)) >= 4) safe_cor(Proximity[is.finite(Proximity)][-sum(is.finite(Proximity))], Proximity[is.finite(Proximity)][-1], "pearson") else NA_real_,
    n_bins = n(),
    .groups = "drop"
  ) %>%
  filter(n_bins >= 4)

latent_feature_cols <- latent_epoch_features %>%
  select(where(is.numeric)) %>%
  names()
latent_feature_cols <- setdiff(latent_feature_cols, c("CageChangeIndex", "n_bins"))

latent_embedding_tbl <- tibble()
latent_trajectory_summary <- tibble()
latent_instability_by_animal <- tibble()
p_latent_traj <- NULL
p_umap_traj <- NULL
p_phate_traj <- NULL
p_instability <- NULL

if (nrow(latent_epoch_features) >= 20 && length(latent_feature_cols) >= 3) {
  latent_mat <- latent_epoch_features %>%
    select(all_of(latent_feature_cols)) %>%
    mutate(across(everything(), ~ replace_na(.x, median(.x, na.rm = TRUE)))) %>%
    as.matrix()
  latent_mat_scaled <- scale(latent_mat)
  latent_mat_scaled[!is.finite(latent_mat_scaled)] <- 0

  temporal_pca <- prcomp(latent_mat_scaled, center = FALSE, scale. = FALSE)
  temporal_var_exp <- temporal_pca$sdev^2 / sum(temporal_pca$sdev^2)
  latent_embedding_tbl <- latent_epoch_features %>%
    bind_cols(tibble(
      TemporalPC1 = temporal_pca$x[, 1],
      TemporalPC2 = temporal_pca$x[, 2]
    ))

  if (requireNamespace("uwot", quietly = TRUE) && nrow(latent_mat_scaled) >= 30) {
    set.seed(2)
    temporal_umap <- uwot::umap(
      latent_mat_scaled,
      n_neighbors = max(5, min(30, floor(nrow(latent_mat_scaled) / 5))),
      min_dist = 0.20,
      metric = "euclidean"
    )
    latent_embedding_tbl <- latent_embedding_tbl %>%
      mutate(UMAP1 = temporal_umap[, 1], UMAP2 = temporal_umap[, 2])
  }
  if (!all(c("UMAP1", "UMAP2") %in% names(latent_embedding_tbl))) {
    latent_embedding_tbl <- latent_embedding_tbl %>%
      mutate(UMAP1 = NA_real_, UMAP2 = NA_real_)
  }
  if (requireNamespace("phateR", quietly = TRUE) && nrow(latent_mat_scaled) >= 30) {
    set.seed(3)
    temporal_phate <- try(
      phateR::phate(
        latent_mat_scaled,
        ndim = 2,
        knn = max(5, min(20, floor(nrow(latent_mat_scaled) / 8))),
        npca = min(20, ncol(latent_mat_scaled)),
        mds.solver = "smacof",
        verbose = 0,
        seed = 3,
        n.jobs = 1
      ),
      silent = TRUE
    )
    if (!inherits(temporal_phate, "try-error") && !is.null(temporal_phate$embedding)) {
      phate_embedding <- as.matrix(temporal_phate$embedding)
      latent_embedding_tbl <- latent_embedding_tbl %>%
        mutate(PHATE1 = phate_embedding[, 1], PHATE2 = phate_embedding[, 2])
    }
  }
  if (!all(c("PHATE1", "PHATE2") %in% names(latent_embedding_tbl))) {
    latent_embedding_tbl <- latent_embedding_tbl %>%
      mutate(PHATE1 = NA_real_, PHATE2 = NA_real_)
  }

  latent_trajectory_summary <- latent_embedding_tbl %>%
    group_by(Sex, Group, PhaseClass, CageChangeIndex) %>%
    summarise(
      n_animals = n_distinct(AnimalNum),
      TemporalPC1 = mean(TemporalPC1, na.rm = TRUE),
      TemporalPC2 = mean(TemporalPC2, na.rm = TRUE),
      UMAP1 = mean(UMAP1, na.rm = TRUE),
      UMAP2 = mean(UMAP2, na.rm = TRUE),
      PHATE1 = mean(PHATE1, na.rm = TRUE),
      PHATE2 = mean(PHATE2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Sex, Group, PhaseClass, CageChangeIndex)

  latent_instability_by_animal <- latent_embedding_tbl %>%
    arrange(AnimalNum, PhaseClass, CageChangeIndex) %>%
    group_by(AnimalNum, Group, Sex, PhaseClass) %>%
    mutate(
      step_distance = sqrt((TemporalPC1 - lag(TemporalPC1))^2 + (TemporalPC2 - lag(TemporalPC2))^2)
    ) %>%
    summarise(
      n_epochs = n(),
      latent_path_length = sum(step_distance, na.rm = TRUE),
      latent_net_displacement = sqrt((last(TemporalPC1) - first(TemporalPC1))^2 + (last(TemporalPC2) - first(TemporalPC2))^2),
      latent_roughness = latent_path_length / pmax(latent_net_displacement, 1e-6),
      latent_variance = var(TemporalPC1, na.rm = TRUE) + var(TemporalPC2, na.rm = TRUE),
      mean_step_distance = mean(step_distance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(systems_features %>% select(AnimalNum, any_of(endpoint_cols)), by = "AnimalNum")

  latent_instability_stats <- latent_instability_by_animal %>%
    filter(primary_outcome %in% names(.)) %>%
    group_by(Sex, PhaseClass) %>%
    summarise(
      n = sum(is.finite(.data[[primary_outcome]]) & is.finite(latent_roughness)),
      roughness_rho = safe_cor(latent_roughness, .data[[primary_outcome]], "spearman"),
      roughness_p = safe_cor_p(latent_roughness, .data[[primary_outcome]], "spearman"),
      path_length_rho = safe_cor(latent_path_length, .data[[primary_outcome]], "spearman"),
      path_length_p = safe_cor_p(latent_path_length, .data[[primary_outcome]], "spearman"),
      variance_rho = safe_cor(latent_variance, .data[[primary_outcome]], "spearman"),
      variance_p = safe_cor_p(latent_variance, .data[[primary_outcome]], "spearman"),
      .groups = "drop"
    ) %>%
    mutate(
      roughness_q = p.adjust(roughness_p, method = "BH"),
      path_length_q = p.adjust(path_length_p, method = "BH"),
      variance_q = p.adjust(variance_p, method = "BH")
    )

  write_table(latent_embedding_tbl, file.path(output_dir, "tables/systems_temporal_latent_epoch_embeddings.csv"))
  write_table(latent_trajectory_summary, file.path(output_dir, "tables/systems_temporal_latent_group_trajectories.csv"))
  write_table(latent_instability_by_animal, file.path(output_dir, "tables/systems_latent_instability_by_animal.csv"))
  write_table(latent_instability_stats, file.path(output_dir, "stats_tables/systems_latent_instability_combz_associations.csv"))

  p_latent_traj <- latent_trajectory_summary %>%
    ggplot(aes(TemporalPC1, TemporalPC2, colour = Group, group = Group)) +
    geom_path(arrow = arrow(length = unit(1.8, "mm"), type = "closed"), linewidth = 0.42, alpha = 0.85) +
    geom_point(aes(size = n_animals), alpha = 0.92) +
    geom_text(aes(label = CageChangeIndex), size = 1.7, vjust = -0.75, show.legend = FALSE) +
    facet_grid(Sex ~ PhaseClass) +
    scale_colour_manual(values = group_colors, drop = FALSE) +
    scale_size_continuous(range = c(1.2, 2.8), guide = "none") +
    labs(
      title = "Temporal evolution through behavioral state space",
      subtitle = paste0("Group mean trajectories across cage changes; temporal PC1/PC2 explain ", round(100 * temporal_var_exp[1], 1), "%/", round(100 * temporal_var_exp[2], 1), "%"),
      x = "Temporal latent PC1",
      y = "Temporal latent PC2",
      caption = "Numbers mark cage-change order; arrows show progression through the social-instability paradigm."
    ) +
    make_nature_theme(base_size = 6)

  save_plot_svg_pdf(p_latent_traj, file.path(output_dir, "figures/publication_panels/Fig_systems_temporal_latent_trajectories"), width = 170, height = 115)

  if (all(c("UMAP1", "UMAP2") %in% names(latent_trajectory_summary))) {
    p_umap_traj <- latent_trajectory_summary %>%
      ggplot(aes(UMAP1, UMAP2, colour = Group, group = Group)) +
      geom_path(arrow = arrow(length = unit(1.8, "mm"), type = "closed"), linewidth = 0.42, alpha = 0.85) +
      geom_point(aes(size = n_animals), alpha = 0.92) +
      geom_text(aes(label = CageChangeIndex), size = 1.7, vjust = -0.75, show.legend = FALSE) +
      facet_grid(Sex ~ PhaseClass) +
      scale_colour_manual(values = group_colors, drop = FALSE) +
      scale_size_continuous(range = c(1.2, 2.8), guide = "none") +
      labs(
        title = "Nonlinear temporal behavioral manifold",
        subtitle = "UMAP of animal-by-phase-by-cage-change dynamics",
        x = "Temporal UMAP1",
        y = "Temporal UMAP2",
        caption = "Numbers mark cage-change order; use as a nonlinear companion to the PCA trajectory."
      ) +
      make_nature_theme(base_size = 6)

    save_plot_svg_pdf(p_umap_traj, file.path(output_dir, "figures/publication_panels/Fig_systems_temporal_umap_trajectories"), width = 170, height = 115)
  }

  if (all(c("PHATE1", "PHATE2") %in% names(latent_trajectory_summary)) && any(is.finite(latent_trajectory_summary$PHATE1))) {
    p_phate_traj <- latent_trajectory_summary %>%
      ggplot(aes(PHATE1, PHATE2, colour = Group, group = Group)) +
      geom_path(arrow = arrow(length = unit(1.8, "mm"), type = "closed"), linewidth = 0.42, alpha = 0.85) +
      geom_point(aes(size = n_animals), alpha = 0.92) +
      geom_text(aes(label = CageChangeIndex), size = 1.7, vjust = -0.75, show.legend = FALSE) +
      facet_grid(Sex ~ PhaseClass) +
      scale_colour_manual(values = group_colors, drop = FALSE) +
      scale_size_continuous(range = c(1.2, 2.8), guide = "none") +
      labs(
        title = "PHATE temporal behavioral manifold",
        subtitle = "Diffusion-geometric embedding of animal-by-phase-by-cage-change dynamics",
        x = "Temporal PHATE1",
        y = "Temporal PHATE2",
        caption = "Numbers mark cage-change order; PHATE emphasizes gradual dynamical progression and attractor-like structure."
      ) +
      make_nature_theme(base_size = 6)

    save_plot_svg_pdf(p_phate_traj, file.path(output_dir, "figures/publication_panels/Fig_systems_temporal_phate_trajectories"), width = 170, height = 115)
  }

  p_instability <- latent_instability_by_animal %>%
    filter(is.finite(latent_roughness)) %>%
    ggplot(aes(Group, latent_roughness, colour = Group, fill = Group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.28, linewidth = 0.28) +
    geom_point(position = position_jitter(width = 0.11, height = 0), size = 1.55, alpha = 0.78, stroke = 0.20) +
    facet_grid(PhaseClass ~ Sex, scales = "free_y") +
    scale_colour_manual(values = group_colors, drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    labs(
      title = "Latent trajectory instability",
      subtitle = "Path roughness through temporal state space by animal",
      x = NULL,
      y = "Trajectory roughness"
    ) +
    make_nature_theme(base_size = 6)

  save_plot_svg_pdf(p_instability, file.path(output_dir, "figures/publication_panels/Fig_systems_latent_trajectory_instability"), width = 135, height = 95)
}

# ------------------------------------------------
# EFFECT-SIZE HEATMAPS
# ------------------------------------------------

focus_domains <- c("behavior", "composite", "temporal_dynamics", "latent_space", "behavioral_state", "trajectory")
heat_tbl <- group_contrasts %>%
  filter(status == "tested", Domain %in% focus_domains) %>%
  mutate(
    effect_size = hedges_g,
    DisplayFeature = paste(Metric, Statistic, Context, sep = " | "),
    DisplayFeature = str_replace_all(DisplayFeature, "_", " "),
    DisplayFeature = str_trunc(DisplayFeature, width = 48),
    contrast = factor(contrast, levels = c("RES-CON", "SUS-CON", "SUS-RES"))
  )

if (publication_focus_only) {
  keep_features <- heat_tbl %>%
    group_by(feature) %>%
    summarise(
      max_abs_d = if (all(is.na(effect_size))) NA_real_ else max(abs(effect_size), na.rm = TRUE),
      min_fdr = if (all(is.na(p_fdr))) NA_real_ else min(p_fdr, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(is.finite(max_abs_d)) %>%
    arrange(desc(max_abs_d), min_fdr) %>%
    slice_head(n = 40) %>%
    pull(feature)
  heat_tbl <- heat_tbl %>% filter(feature %in% keep_features)
}

p_heat <- heat_tbl %>%
  mutate(DisplayFeature = factor(DisplayFeature, levels = rev(unique(DisplayFeature)))) %>%
  ggplot(aes(contrast, DisplayFeature, fill = effect_size)) +
  geom_tile(colour = "white", linewidth = 0.20) +
  geom_text(aes(label = sig), size = 2.1) +
  facet_grid(Sex ~ Domain, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#3d3b6e", mid = "white", high = "#e63947", midpoint = 0, na.value = "grey90") +
  labs(
    title = "Systems-level group-difference map",
    subtitle = "Hedges g; symbols denote BH FDR within prespecified feature families",
    x = NULL,
    y = NULL,
    fill = "Hedges g",
    caption = "Symbols: * q<0.05, ** q<0.01, *** q<0.001; exact statistics in systems_group_contrasts.csv. RES/SUS labels are derived from post-paradigm CombZ."
  ) +
  make_nature_theme(base_size = 6) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "right",
    panel.grid = element_blank()
  )

save_plot_svg_pdf(p_heat, file.path(output_dir, "figures/publication_panels/Fig_systems_effect_size_heatmap"), width = 190, height = 165)

# ------------------------------------------------
# FEATURE NETWORKS: CORRELATION STRUCTURE BY SEX
# ------------------------------------------------

network_features <- group_contrasts %>%
  filter(status == "tested", !is.na(hedges_g)) %>%
  group_by(feature) %>%
  summarise(max_abs_d = if (all(is.na(hedges_g))) NA_real_ else max(abs(hedges_g), na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(max_abs_d)) %>%
  arrange(desc(max_abs_d)) %>%
  slice_head(n = 25) %>%
  pull(feature)

network_long <- map_dfr(levels(droplevels(systems_features$Sex)), function(sx) {
  dat <- systems_features %>% filter(as.character(Sex) == sx)
  if (nrow(dat) < 5) return(tibble())
  feats <- intersect(network_features, names(dat))
  if (length(feats) < 3) return(tibble())
  mat <- dat %>% select(all_of(feats)) %>% mutate(across(everything(), ~ replace_na(.x, median(.x, na.rm = TRUE)))) %>% as.matrix()
  cmat <- suppressWarnings(cor(mat, method = "spearman", use = "pairwise.complete.obs"))
  idx <- which(upper.tri(cmat), arr.ind = TRUE)
  tibble(
    Sex = sx,
    Feature1 = colnames(cmat)[idx[, 1]],
    Feature2 = colnames(cmat)[idx[, 2]],
    rho = cmat[idx]
  ) %>%
    filter(is.finite(rho), abs(rho) >= 0.50) %>%
    left_join(feature_dictionary %>% select(Feature1 = feature, Metric1 = Metric, Domain1 = Domain), by = "Feature1") %>%
    left_join(feature_dictionary %>% select(Feature2 = feature, Metric2 = Metric, Domain2 = Domain), by = "Feature2")
})

write_table(network_long, file.path(output_dir, "tables/systems_feature_correlation_network_edges.csv"))

if (nrow(network_long) > 0) {
  network_plot_tbl <- network_long %>%
    mutate(
      Pair = paste(str_trunc(Metric1, 18), str_trunc(Metric2, 18), sep = " <-> "),
      Pair = factor(Pair, levels = unique(Pair[order(abs(rho), decreasing = TRUE)]))
    ) %>%
    group_by(Sex) %>%
    slice_max(abs(rho), n = 20, with_ties = FALSE) %>%
    ungroup()

  p_net <- network_plot_tbl %>%
    ggplot(aes(rho, Pair, fill = rho)) +
    geom_col(width = 0.72) +
    facet_grid(Sex ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = "#3d3b6e", mid = "white", high = "#e63947", midpoint = 0) +
    labs(
      title = "Sex-specific systems coupling",
      subtitle = "Top feature-feature Spearman correlations among strongest group-discriminating features",
      x = "Spearman rho",
      y = NULL,
      fill = "rho"
    ) +
    make_nature_theme(base_size = 6) +
    theme(legend.position = "right")

  save_plot_svg_pdf(p_net, file.path(output_dir, "figures/publication_panels/Fig_systems_feature_coupling_network_summary"), width = 150, height = 120)
}

# ------------------------------------------------
# OUTCOME INTEGRATION AND PREDICTION SUMMARY
# ------------------------------------------------

available_outcomes <- intersect(endpoint_cols, names(systems_features))
available_outcomes <- available_outcomes[map_lgl(available_outcomes, ~ any(is.finite(safe_numeric(systems_features[[.x]]))))]

if (length(available_outcomes) > 0) {
  calc_outcome_assoc <- function(feature_set, feature_set_name) {
    map_dfr(available_outcomes, function(outcome) {
      y <- safe_numeric(systems_features[[outcome]])
      map_dfr(feature_set, function(fc) {
        n_pair <- sum(is.finite(y) & is.finite(systems_features[[fc]]))
        rho <- safe_cor(systems_features[[fc]], y, "spearman")
        rho_ci <- cor_ci_approx(rho, n_pair)
        pear <- safe_cor(systems_features[[fc]], y, "pearson")
        pear_ci <- cor_ci_approx(pear, n_pair)
        tibble(
          outcome = outcome,
          FeatureSet = feature_set_name,
          feature = fc,
          n = n_pair,
          spearman_rho = rho,
          spearman_ci_low = rho_ci["low"],
          spearman_ci_high = rho_ci["high"],
          spearman_p = safe_cor_p(systems_features[[fc]], y, "spearman"),
          pearson_r = pear,
          pearson_ci_low = pear_ci["low"],
          pearson_ci_high = pear_ci["high"],
          pearson_p = safe_cor_p(systems_features[[fc]], y, "pearson")
        )
      })
    }) %>%
      group_by(outcome, FeatureSet) %>%
      mutate(
        spearman_fdr = p.adjust(spearman_p, method = "BH"),
        sig = sig_label(spearman_fdr),
        ReportingCorrection = "BH FDR within outcome x feature set",
        TestMethod = "Spearman rank correlation; Pearson correlation also reported",
        ConfidenceInterval = "Approximate Fisher-z 95% CI for correlation coefficients",
        Direction = case_when(
          outcome == primary_outcome & primary_outcome_direction == "lower_worse" & spearman_rho > 0 ~ "higher_feature_less_worse",
          outcome == primary_outcome & primary_outcome_direction == "lower_worse" & spearman_rho < 0 ~ "higher_feature_more_worse",
          TRUE ~ "feature_outcome_association"
        )
      ) %>%
      ungroup() %>%
      left_join(feature_dictionary, by = "feature") %>%
      arrange(outcome, FeatureSet, spearman_fdr, desc(abs(spearman_rho)))
  }

  outcome_assoc <- calc_outcome_assoc(usable_features, "descriptive_full_experiment")
  prospective_outcome_assoc <- calc_outcome_assoc(prospective_prediction_features, "prospective_prediction")

  write_table(outcome_assoc, file.path(output_dir, "stats_tables/systems_outcome_associations.csv"))
  write_table(prospective_outcome_assoc, file.path(output_dir, "stats_tables/systems_prospective_outcome_associations.csv"))

  outcome_to_plot <- if (primary_outcome %in% available_outcomes) primary_outcome else available_outcomes[1]
  top_outcome_features <- outcome_assoc %>%
    filter(outcome == outcome_to_plot, !is.na(spearman_rho)) %>%
    arrange(spearman_fdr, desc(abs(spearman_rho))) %>%
    slice_head(n = 25)

  p_out_heat <- top_outcome_features %>%
    mutate(DisplayFeature = paste(Metric, Statistic, Context, sep = " | "),
           DisplayFeature = str_replace_all(DisplayFeature, "_", " "),
           DisplayFeature = factor(str_trunc(DisplayFeature, 52), levels = rev(str_trunc(DisplayFeature, 52))),
           StatLabel = paste(rho_label(spearman_rho), q_label(spearman_fdr), sep = "\n"),
           LabelX = if_else(spearman_rho >= 0, pmin(spearman_rho + 0.04, 0.96), pmax(spearman_rho - 0.04, -0.96)),
           LabelHJust = if_else(spearman_rho >= 0, 0, 1)) %>%
    ggplot(aes(spearman_rho, DisplayFeature, fill = spearman_rho)) +
    geom_col(width = 0.72) +
    geom_text(aes(x = LabelX, label = StatLabel, hjust = LabelHJust), size = 1.75, lineheight = 0.86) +
    scale_fill_gradient2(low = "#3d3b6e", mid = "white", high = "#e63947", midpoint = 0) +
    labs(
      title = paste0("Full-experiment systems features linked to ", outcome_to_plot),
      subtitle = paste0("Exploratory associations; labels show Spearman rho and BH q"),
      x = "Spearman rho",
      y = NULL,
      fill = "rho",
      caption = paste0("For ", primary_outcome, ", positive rho means less-worse/resilient-like endpoint; exact statistics in systems_outcome_associations.csv.")
    ) +
    coord_cartesian(xlim = c(-1, 1), clip = "off") +
    make_nature_theme(base_size = 6) +
    theme(legend.position = "right", plot.margin = margin(5.5, 18, 5.5, 5.5))

  save_plot_svg_pdf(p_out_heat, file.path(output_dir, "figures/publication_panels/Fig_systems_outcome_association_rank"), width = 135, height = 115)

  top_prospective_features <- prospective_outcome_assoc %>%
    filter(outcome == outcome_to_plot, !is.na(spearman_rho)) %>%
    arrange(spearman_fdr, desc(abs(spearman_rho))) %>%
    slice_head(n = 25)

  if (nrow(top_prospective_features) > 0) {
    p_prosp_heat <- top_prospective_features %>%
      mutate(DisplayFeature = paste(Metric, Statistic, Context, sep = " | "),
             DisplayFeature = str_replace_all(DisplayFeature, "_", " "),
             DisplayFeature = factor(str_trunc(DisplayFeature, 52), levels = rev(str_trunc(DisplayFeature, 52))),
             StatLabel = paste(rho_label(spearman_rho), q_label(spearman_fdr), sep = "\n"),
             LabelX = if_else(spearman_rho >= 0, pmin(spearman_rho + 0.04, 0.96), pmax(spearman_rho - 0.04, -0.96)),
             LabelHJust = if_else(spearman_rho >= 0, 0, 1)) %>%
      ggplot(aes(spearman_rho, DisplayFeature, fill = spearman_rho)) +
      geom_col(width = 0.72) +
      geom_text(aes(x = LabelX, label = StatLabel, hjust = LabelHJust), size = 1.75, lineheight = 0.86) +
      scale_fill_gradient2(low = "#3d3b6e", mid = "white", high = "#e63947", midpoint = 0) +
      labs(
        title = paste0("Prospective early behavior linked to ", outcome_to_plot),
        subtitle = paste0("First active 12 h after ", first_cage_change, "; labels show Spearman rho and BH q"),
        x = "Spearman rho",
        y = NULL,
        fill = "rho",
        caption = paste0("For ", primary_outcome, ", negative rho means higher early feature values predict worse later endpoint.")
      ) +
      coord_cartesian(xlim = c(-1, 1), clip = "off") +
      make_nature_theme(base_size = 6) +
      theme(legend.position = "right", plot.margin = margin(5.5, 18, 5.5, 5.5))

    save_plot_svg_pdf(p_prosp_heat, file.path(output_dir, "figures/publication_panels/Fig_systems_prospective_outcome_association_rank"), width = 135, height = 115)

    top_prospective_scatter_features <- top_prospective_features %>%
      slice_head(n = 6) %>%
      pull(feature)

    prospective_scatter_tbl <- systems_features %>%
      select(AnimalNum, Group, Sex, all_of(outcome_to_plot), all_of(top_prospective_scatter_features)) %>%
      rename(outcome = all_of(outcome_to_plot)) %>%
      pivot_longer(all_of(top_prospective_scatter_features), names_to = "feature", values_to = "FeatureValue") %>%
      filter(is.finite(outcome), is.finite(FeatureValue)) %>%
      left_join(
        top_prospective_features %>%
          select(feature, Metric, Statistic, Context, spearman_rho, spearman_fdr, Direction),
        by = "feature"
      ) %>%
      mutate(
        DisplayFeature = paste(Metric, Statistic, Context, sep = " | "),
        DisplayFeature = str_replace_all(DisplayFeature, "_", " "),
        DisplayFeature = factor(str_trunc(DisplayFeature, 46), levels = unique(str_trunc(DisplayFeature, 46)))
      )

    prospective_scatter_stats <- prospective_scatter_tbl %>%
      group_by(Sex, DisplayFeature) %>%
      summarise(
        n = n_distinct(AnimalNum),
        spearman_rho_sex = safe_cor(FeatureValue, outcome, "spearman"),
        spearman_p_sex = safe_cor_p(FeatureValue, outcome, "spearman"),
        .groups = "drop"
      ) %>%
      group_by(Sex) %>%
      mutate(
        spearman_q_sex = p.adjust(spearman_p_sex, method = "BH"),
        StatLabel = paste0(rho_label(spearman_rho_sex), "\n", q_label(spearman_q_sex), "\nn=", n)
      ) %>%
      ungroup()

    write_table(prospective_scatter_stats, file.path(output_dir, "stats_tables/systems_prospective_feature_scatter_sex_stats.csv"))

    p_prosp_scatter <- prospective_scatter_tbl %>%
      ggplot(aes(FeatureValue, outcome, colour = Group, fill = Group)) +
      geom_point(size = 1.7, alpha = 0.82, stroke = 0.22) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 0.38, colour = "grey20", fill = "grey70", alpha = 0.18) +
      geom_text(
        data = prospective_scatter_stats,
        aes(x = -Inf, y = Inf, label = StatLabel),
        inherit.aes = FALSE,
        hjust = -0.05,
        vjust = 1.10,
        size = 1.85,
        lineheight = 0.86,
        colour = "grey15"
      ) +
      facet_grid(Sex ~ DisplayFeature, scales = "free_x") +
      scale_colour_manual(values = group_colors, drop = FALSE) +
      scale_fill_manual(values = group_colors, drop = FALSE) +
      labs(
        title = paste0("Early behavioral predictors of ", outcome_to_plot),
        subtitle = "Top prospective first-active-window associations; panel labels show sex-stratified rho, BH q and n",
        x = "Feature value",
        y = paste0(outcome_to_plot, " (lower = worse)"),
        caption = "Lines show ordinary least-squares trend for visualization; panel statistics use Spearman correlation with BH FDR within sex across displayed features."
      ) +
      make_nature_theme(base_size = 6)

    save_plot_svg_pdf(p_prosp_scatter, file.path(output_dir, "figures/publication_panels/Fig_systems_prospective_feature_scatter"), width = 190, height = 118)
  }

  if (requireNamespace("glmnet", quietly = TRUE)) {
    model_dat <- systems_features %>%
      filter(is.finite(safe_numeric(.data[[outcome_to_plot]])))
    model_features <- prospective_prediction_features[prospective_prediction_features %in% names(model_dat)]
    model_feature_qc <- map_dfr(model_features, function(fc) {
      x <- safe_numeric(model_dat[[fc]])
      tibble(
        feature = fc,
        n_finite = sum(is.finite(x)),
        missing_fraction = mean(!is.finite(x)),
        zero_variance = sd(x, na.rm = TRUE) == 0
      )
    })
    model_features <- model_feature_qc %>%
      filter(n_finite >= 8, missing_fraction <= 0.50, !zero_variance) %>%
      pull(feature)

    x_model_raw <- model_dat %>%
      select(all_of(model_features)) %>%
      mutate(across(everything(), safe_numeric)) %>%
      as.matrix()
    y_model <- safe_numeric(model_dat[[outcome_to_plot]])

    if (nrow(x_model_raw) >= 8 && ncol(x_model_raw) >= 2) {
      pred <- rep(NA_real_, length(y_model))
      coef_tbl <- tibble()
      set.seed(1)
      for (i in seq_along(y_model)) {
        train_idx <- setdiff(seq_along(y_model), i)
        x_train <- x_model_raw[train_idx, , drop = FALSE]
        x_test <- x_model_raw[i, , drop = FALSE]
        train_medians <- apply(x_train, 2, median, na.rm = TRUE)
        train_medians[!is.finite(train_medians)] <- 0
        for (j in seq_along(train_medians)) {
          x_train[!is.finite(x_train[, j]), j] <- train_medians[j]
          x_test[!is.finite(x_test[, j]), j] <- train_medians[j]
        }
        fit <- glmnet::cv.glmnet(
          x_train, y_model[train_idx],
          alpha = 0.5,
          standardize = TRUE,
          nfolds = min(5, length(train_idx))
        )
        pred[i] <- as.numeric(predict(fit, x_test, s = "lambda.min"))
        cc <- as.matrix(coef(fit, s = "lambda.min"))[, 1]
        coef_tbl <- bind_rows(coef_tbl, tibble(held_out = model_dat$AnimalNum[i], feature = names(cc), coefficient = as.numeric(cc)))
      }

      prediction_tbl <- model_dat %>%
        transmute(AnimalNum, Group, Sex, outcome = y_model, predicted = pred, residual = outcome - predicted)
      baseline <- mean(prediction_tbl$outcome, na.rm = TRUE)
      loo_r <- safe_cor(prediction_tbl$outcome, prediction_tbl$predicted, "pearson")
      loo_r_ci <- cor_ci_approx(loo_r, nrow(prediction_tbl))
      loo_rho <- safe_cor(prediction_tbl$outcome, prediction_tbl$predicted, "spearman")
      loo_rho_ci <- cor_ci_approx(loo_rho, nrow(prediction_tbl))
      prediction_perf <- tibble(
        outcome = outcome_to_plot,
        FeatureSet = "prospective_prediction",
        n = nrow(prediction_tbl),
        n_features = length(model_features),
        loo_pearson_r = loo_r,
        loo_pearson_ci_low = loo_r_ci["low"],
        loo_pearson_ci_high = loo_r_ci["high"],
        loo_spearman_rho = loo_rho,
        loo_spearman_ci_low = loo_rho_ci["low"],
        loo_spearman_ci_high = loo_rho_ci["high"],
        loo_rmse = sqrt(mean((prediction_tbl$outcome - prediction_tbl$predicted)^2, na.rm = TRUE)),
        baseline_rmse = sqrt(mean((prediction_tbl$outcome - baseline)^2, na.rm = TRUE)),
        cv_r2_vs_mean = 1 - sum((prediction_tbl$outcome - prediction_tbl$predicted)^2, na.rm = TRUE) /
          sum((prediction_tbl$outcome - baseline)^2, na.rm = TRUE),
        Model = "Elastic net, leave-one-animal-out cross-validation",
        Preprocessing = "Missing feature values imputed from training-fold medians only"
      )

      set.seed(2)
      pearson_perm <- abs(replicate(1000, safe_cor(prediction_tbl$outcome, sample(prediction_tbl$predicted), "pearson")))
      spearman_perm <- abs(replicate(1000, safe_cor(prediction_tbl$outcome, sample(prediction_tbl$predicted), "spearman")))
      prediction_permutation <- tibble(
        n_permutations = 1000,
        pearson_r_permutation_p = (sum(pearson_perm >= abs(loo_r), na.rm = TRUE) + 1) / (sum(is.finite(pearson_perm)) + 1),
        spearman_rho_permutation_p = (sum(spearman_perm >= abs(loo_rho), na.rm = TRUE) + 1) / (sum(is.finite(spearman_perm)) + 1),
        NullModel = "Permutation of cross-validated predictions relative to observed endpoint; model not refit"
      )
      prediction_perf <- bind_cols(prediction_perf, prediction_permutation)

      prediction_auc_tbl <- tibble()
      if (requireNamespace("pROC", quietly = TRUE)) {
        prediction_auc_tbl <- prediction_tbl %>%
          filter(Group %in% c("RES", "SUS"), is.finite(predicted)) %>%
          group_by(Sex) %>%
          group_modify(~ {
            if (n_distinct(.x$Group) < 2) return(tibble(n = nrow(.x), auc_res_sus = NA_real_))
            roc_obj <- try(pROC::roc(response = .x$Group, predictor = -.x$predicted, levels = c("RES", "SUS"), quiet = TRUE), silent = TRUE)
            tibble(n = nrow(.x), auc_res_sus = if (inherits(roc_obj, "try-error")) NA_real_ else as.numeric(pROC::auc(roc_obj)))
          }) %>%
          ungroup() %>%
          bind_rows(
            prediction_tbl %>%
              filter(Group %in% c("RES", "SUS"), is.finite(predicted)) %>%
              { if (n_distinct(.$Group) < 2) tibble(Sex = "Pooled", n = nrow(.), auc_res_sus = NA_real_) else {
                roc_obj <- try(pROC::roc(response = .$Group, predictor = -.$predicted, levels = c("RES", "SUS"), quiet = TRUE), silent = TRUE)
                tibble(Sex = "Pooled", n = nrow(.), auc_res_sus = if (inherits(roc_obj, "try-error")) NA_real_ else as.numeric(pROC::auc(roc_obj)))
              }}
          ) %>%
          mutate(
            Outcome = "RES vs SUS classification from cross-validated predicted CombZ",
            ScoreDirection = "Higher score = lower predicted CombZ = more susceptible-like"
          )
      }

      coef_summary <- coef_tbl %>%
        filter(feature != "(Intercept)") %>%
        group_by(feature) %>%
        summarise(nonzero_frequency = mean(coefficient != 0), median_coefficient = median(coefficient), .groups = "drop") %>%
        left_join(feature_dictionary, by = "feature") %>%
        arrange(desc(nonzero_frequency), desc(abs(median_coefficient)))

      write_table(prediction_tbl, file.path(output_dir, "tables/systems_prospective_outcome_loo_predictions.csv"))
      write_table(prediction_perf, file.path(output_dir, "tables/systems_prospective_outcome_loo_performance.csv"))
      write_table(coef_summary, file.path(output_dir, "tables/systems_prospective_outcome_elastic_net_feature_stability.csv"))
      write_table(prediction_auc_tbl, file.path(output_dir, "tables/systems_prospective_outcome_res_sus_auc.csv"))
      write_table(prediction_tbl, file.path(output_dir, "tables/systems_outcome_loo_predictions.csv"))
      write_table(prediction_perf, file.path(output_dir, "tables/systems_outcome_loo_performance.csv"))
      write_table(coef_summary, file.path(output_dir, "tables/systems_outcome_elastic_net_feature_stability.csv"))
      write_table(prediction_auc_tbl, file.path(output_dir, "tables/systems_outcome_res_sus_auc.csv"))

      prediction_sex_stats <- prediction_tbl %>%
        group_by(Sex) %>%
        summarise(
          n = n_distinct(AnimalNum),
          pearson_r = safe_cor(outcome, predicted, "pearson"),
          spearman_rho = safe_cor(outcome, predicted, "spearman"),
          rmse = sqrt(mean((outcome - predicted)^2, na.rm = TRUE)),
          .groups = "drop"
        ) %>%
        mutate(
          StatLabel = paste0(
            "r=", formatC(pearson_r, format = "f", digits = 2),
            "\n", rho_label(spearman_rho),
            "\nRMSE=", formatC(rmse, format = "f", digits = 2),
            "\nn=", n
          )
        )

      write_table(prediction_sex_stats, file.path(output_dir, "tables/systems_prospective_outcome_loo_performance_by_sex.csv"))

      p_pred <- prediction_tbl %>%
        ggplot(aes(outcome, predicted, colour = Group, fill = Group)) +
        geom_abline(slope = 1, intercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey45") +
        geom_point(size = 2.2, alpha = 0.9) +
        geom_smooth(method = "lm", se = TRUE, linewidth = 0.45, colour = "grey25") +
        geom_text(
          data = prediction_sex_stats,
          aes(x = -Inf, y = Inf, label = StatLabel),
          inherit.aes = FALSE,
          hjust = -0.06,
          vjust = 1.10,
          size = 2.05,
          lineheight = 0.86,
          colour = "grey15"
        ) +
        facet_wrap(~ Sex, nrow = 1) +
        scale_colour_manual(values = group_colors, drop = FALSE) +
        scale_fill_manual(values = group_colors, drop = FALSE) +
        labs(
          title = paste0("Prospective early behavior predicts ", outcome_to_plot),
          subtitle = paste0("LOO elastic-net on ", length(model_features), " early features: pooled r=", round(prediction_perf$loo_pearson_r, 2), "; pooled CV R2=", round(prediction_perf$cv_r2_vs_mean, 2)),
          x = paste0("Observed ", outcome_to_plot, " (lower = worse)"),
          y = paste0("Predicted ", outcome_to_plot),
          caption = "Panel labels report sex-stratified prediction performance; model fitting was leave-one-animal-out across all animals."
        ) +
        make_nature_theme(base_size = 7)

      save_plot_svg_pdf(p_pred, file.path(output_dir, "figures/publication_panels/Fig_systems_prospective_crossvalidated_prediction"), width = 135, height = 75)
    }
  }
}

# ------------------------------------------------
# PUBLICATION DASHBOARD
# ------------------------------------------------

# A compact, paper-facing composite figure. If patchwork is not installed,
# individual panels above remain the primary outputs.
if (requireNamespace("patchwork", quietly = TRUE)) {
  dashboard_feature_label <- function(metric, statistic, context, width = 42) {
    label <- paste(
      str_replace_all(metric, "_", " "),
      str_replace_all(statistic, "_", " "),
      str_replace_all(context, "_", " "),
      sep = " | "
    )
    label <- str_replace(label, "mean \\| first active 12h", "early mean")
    label <- str_replace(label, "sd \\| first active 12h", "early variability")
    label <- str_replace(label, "p95 \\| first active 12h", "early high values")
    label <- str_replace(label, "q75 \\| first active 12h", "early upper quartile")
    label <- str_replace(label, "acf1 \\| first active 12h", "early autocorrelation")
    label <- str_replace(label, "rmssd \\| first active 12h", "early instability")
    str_trunc(label, width = width)
  }

  top_loadings <- bind_rows(
    pca_loadings %>% filter(PC1 > 0) %>% arrange(desc(PC1)) %>% slice_head(n = 6),
    pca_loadings %>% filter(PC1 < 0) %>% arrange(PC1) %>% slice_head(n = 6)
  )
  if (nrow(top_loadings) < 8) {
    top_loadings <- pca_loadings %>%
      arrange(desc(abs(PC1))) %>%
      slice_head(n = 12)
  }
  top_loadings <- top_loadings %>%
    mutate(
      DisplayFeature = dashboard_feature_label(Metric, Statistic, Context),
      DisplayFeature = factor(DisplayFeature, levels = rev(unique(DisplayFeature)))
    )

  p_load <- top_loadings %>%
    ggplot(aes(PC1, DisplayFeature, fill = PC1)) +
    geom_vline(xintercept = 0, linewidth = 0.20, colour = "grey55") +
    geom_col(width = 0.72) +
    scale_fill_gradient2(low = "#3d3b6e", mid = "white", high = "#e63947", midpoint = 0) +
    labs(
      title = "B. Dominant systems axis",
      subtitle = paste0("Top signed PC1 loadings; PC1 explains ", round(100 * var_exp[1], 1), "%"),
      x = "PC1 loading",
      y = NULL
    ) +
    make_nature_theme(base_size = 6) +
    theme(legend.position = "none")

  pca_outcome_tbl <- if (exists("outcome_to_plot") && outcome_to_plot %in% names(pca_scores)) {
    pca_scores %>%
      select(AnimalNum, Group, Sex, PC1, all_of(outcome_to_plot)) %>%
      rename(outcome = all_of(outcome_to_plot)) %>%
      filter(is.finite(PC1), is.finite(outcome))
  } else {
    tibble()
  }

  pca_outcome_stats <- if (nrow(pca_outcome_tbl) > 0) {
    pca_outcome_tbl %>%
      group_by(Sex) %>%
      summarise(
        n = n_distinct(AnimalNum),
        spearman_rho = safe_cor(PC1, outcome, "spearman"),
        spearman_p = safe_cor_p(PC1, outcome, "spearman"),
        .groups = "drop"
      ) %>%
      mutate(
        spearman_q = p.adjust(spearman_p, method = "BH"),
        StatLabel = paste0(rho_label(spearman_rho), "\n", q_label(spearman_q), "\nn=", n)
      )
  } else {
    tibble()
  }
  if (nrow(pca_outcome_stats) > 0) {
    write_table(pca_outcome_stats, file.path(output_dir, "stats_tables/systems_pc1_combz_association_by_sex.csv"))
  }

  p_pc1_outcome <- if (nrow(pca_outcome_tbl) > 0) {
    pca_outcome_tbl %>%
      ggplot(aes(PC1, outcome, colour = Group, fill = Group)) +
      geom_point(size = 1.8, alpha = 0.86, stroke = 0.22) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 0.35, colour = "grey25", fill = "grey70", alpha = 0.18) +
      geom_text(
        data = pca_outcome_stats,
        aes(x = -Inf, y = Inf, label = StatLabel),
        inherit.aes = FALSE,
        hjust = -0.06,
        vjust = 1.12,
        size = 1.9,
        lineheight = 0.86,
        colour = "grey15"
      ) +
      facet_wrap(~ Sex, nrow = 1) +
      scale_colour_manual(values = group_colors, drop = FALSE) +
      scale_fill_manual(values = group_colors, drop = FALSE) +
      labs(
        title = paste0("C. Integrated state axis vs ", outcome_to_plot),
        subtitle = "Sex-stratified association with post-paradigm endpoint",
        x = "Systems PC1",
        y = paste0(outcome_to_plot, " (lower = worse)")
      ) +
      make_nature_theme(base_size = 6)
  } else {
    ggplot() +
      annotate("text", x = 0, y = 0, label = "Endpoint association unavailable", size = 3) +
      theme_void()
  }

  p_pred_dash <- if (exists("prediction_tbl") && exists("prediction_sex_stats") && nrow(prediction_tbl) > 0) {
    prediction_tbl %>%
      ggplot(aes(outcome, predicted, colour = Group, fill = Group)) +
      geom_abline(slope = 1, intercept = 0, linewidth = 0.22, linetype = "dashed", colour = "grey45") +
      geom_point(size = 1.8, alpha = 0.88, stroke = 0.22) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 0.35, colour = "grey25", fill = "grey70", alpha = 0.18) +
      geom_text(
        data = prediction_sex_stats,
        aes(x = -Inf, y = Inf, label = StatLabel),
        inherit.aes = FALSE,
        hjust = -0.06,
        vjust = 1.12,
        size = 1.9,
        lineheight = 0.86,
        colour = "grey15"
      ) +
      facet_wrap(~ Sex, nrow = 1) +
      scale_colour_manual(values = group_colors, drop = FALSE) +
      scale_fill_manual(values = group_colors, drop = FALSE) +
      labs(
        title = "D. Prospective endpoint prediction",
        subtitle = paste0("First active 12 h; pooled r=", round(prediction_perf$loo_pearson_r, 2), ", CV R2=", round(prediction_perf$cv_r2_vs_mean, 2)),
        x = paste0("Observed ", outcome_to_plot),
        y = paste0("Predicted ", outcome_to_plot)
      ) +
      make_nature_theme(base_size = 6)
  } else {
    ggplot() +
      annotate("text", x = 0, y = 0, label = "Prospective prediction unavailable", size = 3) +
      theme_void()
  }

  group_n_panel <- systems_features %>%
    distinct(AnimalNum, Sex, Group) %>%
    count(Sex, Group, name = "n") %>%
    mutate(Group = factor(Group, levels = group_levels))
  write_table(group_n_panel, file.path(output_dir, "tables/systems_dashboard_group_n_by_sex.csv"))

  p_heat_small <- heat_tbl %>%
    group_by(Sex, contrast) %>%
    slice_max(abs(hedges_g), n = 6, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(DisplayFeature = dashboard_feature_label(Metric, Statistic, Context, width = 34),
           DisplayFeature = factor(DisplayFeature, levels = rev(unique(DisplayFeature)))) %>%
    ggplot(aes(contrast, DisplayFeature, fill = hedges_g)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(aes(label = sig), size = 1.8) +
    facet_grid(Sex ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = "#3d3b6e", mid = "white", high = "#e63947", midpoint = 0, na.value = "grey90") +
    labs(
      title = "F. Sex-stratified phenotype contrasts",
      subtitle = "Top Hedges g effects per sex and contrast; exact values in group-contrast table",
      x = NULL,
      y = NULL,
      fill = "g"
    ) +
    make_nature_theme(base_size = 5.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

  p_latent_dash <- if (inherits(p_latent_traj, "ggplot")) {
    p_latent_traj +
      labs(
        title = "A. Temporal behavioral trajectories",
        subtitle = "Group mean paths through latent state space across cage changes"
      ) +
      theme(legend.position = "top")
  } else {
    p_pca +
      labs(
        title = "A. Integrated behavioral state space",
        subtitle = paste0("PCA on ", length(usable_features), " animal-level features; points are animals")
      )
  }

  p_nonlinear_dash <- if (inherits(p_phate_traj, "ggplot")) {
    p_phate_traj +
      labs(
        title = "B. PHATE temporal manifold",
        subtitle = "Diffusion-geometric view of behavioral progression"
      ) +
      theme(legend.position = "top")
  } else if (inherits(p_umap_traj, "ggplot")) {
    p_umap_traj +
      labs(
        title = "B. Nonlinear temporal manifold",
        subtitle = "UMAP companion view of the same animal-epoch dynamics"
      ) +
      theme(legend.position = "top")
  } else {
    p_load
  }

  p_instability_dash <- if (inherits(p_instability, "ggplot")) {
    p_instability +
      labs(
        title = "E. Latent trajectory instability",
        subtitle = "Animal-level path roughness through behavioral state space"
      ) +
      theme(legend.position = "top")
  } else {
    ggplot() +
      annotate("text", x = 0, y = 0, label = "Latent instability unavailable", size = 3) +
      theme_void()
  }

  dashboard <- ((p_latent_dash | p_nonlinear_dash) / (p_pc1_outcome | p_pred_dash) / (p_instability_dash | p_heat_small)) +
    patchwork::plot_layout(heights = c(1.05, 0.95, 1.20)) +
    patchwork::plot_annotation(
      title = "Multiscale behavioral signatures of stress susceptibility",
      subtitle = paste0("Stress phenotypes as temporal dynamics across movement, entropy, proximity and state transitions; lower ", primary_outcome, " = worse endpoint"),
      caption = paste0(
        "CON/RES/SUS colours use original palette. RES/SUS are post-paradigm CombZ-derived labels; ",
        "prediction uses only first-active-window behavioral features. Heatmap symbols denote BH FDR: * q<0.05, ** q<0.01, *** q<0.001."
      )
    )

  save_plot_svg_pdf(dashboard, file.path(output_dir, "figures/Fig_integrated_systems_dashboard"), width = 210, height = 250)
}

# ------------------------------------------------
# TEXT SUMMARY FOR MANUSCRIPT / LAB MEETING
# ------------------------------------------------

strong_findings <- group_contrasts %>%
  filter(evidence %in% c("large_FDR_supported", "FDR_supported", "large_effect_uncertain")) %>%
  arrange(match(evidence, c("large_FDR_supported", "FDR_supported", "large_effect_uncertain")), desc(abs(hedges_g))) %>%
  select(
    Sex, contrast, Source, Domain, Scale, Metric, Statistic, Context,
    n_ref, n_comp, mean_ref, mean_comp, estimate, estimate_ci_low, estimate_ci_high,
    cohen_d, hedges_g, p.value, wilcox_p, p_fdr, p_fdr_sex_contrast, evidence,
    test_method, ReportingCorrection
  ) %>%
  slice_head(n = 50)

systems_summary <- tibble(
  Item = c(
    "Animals included",
    "Usable descriptive features",
    "Prospective prediction features",
    "Primary bin level",
    "Prospective window",
    "Primary outcome direction",
    "RES/SUS labels",
    "Main feature domains",
    "Interpretation"
  ),
  Value = c(
    as.character(n_distinct(systems_features$AnimalNum)),
    as.character(length(usable_features)),
    as.character(length(prospective_prediction_features)),
    primary_bin_level,
    paste0("First active 12 h after ", first_cage_change),
    primary_outcome_label,
    "Derived from post-paradigm CombZ; group contrasts are descriptive phenotype contrasts, not independent endpoint validation.",
    paste(sort(unique(feature_dictionary$Domain)), collapse = "; "),
    "Use PCA/UMAP and effect-size heatmaps for descriptive system-level organization. Use prospective prediction outputs for pre-endpoint prediction claims."
  )
)

write_table(strong_findings, file.path(output_dir, "tables/systems_top_findings_for_reporting.csv"))
write_table(systems_summary, file.path(output_dir, "tables/systems_analysis_summary.csv"))

stats_reporting_guide <- tibble(
  Output = c(
    "systems_group_summary.csv",
    "systems_group_contrasts.csv",
    "systems_outcome_associations.csv",
    "systems_prospective_outcome_associations.csv",
    "systems_prospective_outcome_loo_performance.csv"
  ),
  PrimaryStatistic = c(
    "Mean, 95% CI, median, IQR, animal n",
    "Welch mean difference, 95% CI, Hedges g, Cohen d, Wilcoxon p",
    "Spearman rho with approximate 95% CI; Pearson r also reported",
    "Spearman rho with approximate 95% CI for pre-endpoint early features only",
    "Leave-one-animal-out elastic-net r, rho, RMSE, CV R2 with approximate correlation CIs"
  ),
  MultipleTesting = c(
    "Not applicable",
    "BH FDR within Sex x Source x Domain x Scale x Context x contrast; broad Sex x contrast FDR provided",
    "BH FDR within outcome x descriptive feature set",
    "BH FDR within outcome x prospective feature set",
    "Cross-validation, no per-feature multiplicity claim"
  ),
  ManuscriptUse = c(
    "Figure legends and Supplementary tables",
    "Primary descriptive phenotype effect-size table",
    "Exploratory full-experiment endpoint associations",
    "Prospective endpoint association claims",
    "Prospective prediction performance claims"
  )
)

write_table(stats_reporting_guide, file.path(output_dir, "tables/systems_stats_reporting_guide.csv"))

message("Integrated systems neuroscience summary complete.")
message("Output: ", output_dir)
message("Primary figure candidates:")
message("  - figures/publication_panels/Fig_systems_state_space_PCA.svg")
message("  - figures/publication_panels/Fig_systems_effect_size_heatmap.svg")
message("  - figures/Fig_integrated_systems_dashboard.svg")
