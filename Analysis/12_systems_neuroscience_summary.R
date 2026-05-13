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
endpoint_file <- NULL
endpoint_cols <- c("CombZ", "stress_z_score", "SucrosePreference", "Corticosterone", "DeltaCorticosterone")
primary_outcome <- "CombZ"

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

if (!exists("make_nature_theme")) {
  make_nature_theme <- function(base_size = 7) {
    theme_classic(base_size = base_size) +
      theme(
        axis.line = element_line(linewidth = 0.25, colour = "black"),
        axis.ticks = element_line(linewidth = 0.2, colour = "black"),
        axis.text = element_text(colour = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_text(face = "bold", hjust = 0),
        plot.subtitle = element_text(hjust = 0)
      )
  }
}

if (!exists("save_plot_svg_pdf")) {
  save_plot_svg_pdf <- function(plot, filename_base, width = 85, height = 65, units = "mm") {
    ensure_dir(dirname(filename_base))
    ggplot2::ggsave(paste0(filename_base, ".svg"), plot, width = width, height = height, units = units)
    ggplot2::ggsave(paste0(filename_base, ".pdf"), plot, width = width, height = height, units = units)
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

read_any_table <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") return(readr::read_csv(path, show_col_types = FALSE))
  if (ext %in% c("tsv", "txt")) return(readr::read_tsv(path, show_col_types = FALSE))
  if (ext == "rds") return(readRDS(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("Install readxl to read Excel files: ", path)
    return(readxl::read_excel(path))
  }
  NULL
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

p_label <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "p<0.001",
    TRUE ~ paste0("p=", formatC(p, format = "f", digits = 3))
  )
}

sig_label <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    p < 0.10 ~ "·",
    TRUE ~ ""
  )
}

standardize_id_columns <- function(dat) {
  animal_col <- first_existing_col(dat, c("AnimalNum", "Animal", "AnimalID", "MouseID", "Mouse", "ID", "RFID", "animal_id"), TRUE, "animal ID")
  group_col <- first_existing_col(dat, c("Group", "Phenotype", "Condition", "Treatment", "StressGroup"), FALSE, "group")
  sex_col <- first_existing_col(dat, c("Sex", "sex"), FALSE, "sex")

  out <- dat %>% mutate(AnimalNum = .data[[animal_col]])
  out$Group <- if (!is.na(group_col)) as.character(dat[[group_col]]) else NA_character_
  out$Sex <- if (!is.na(sex_col)) as.character(dat[[sex_col]]) else NA_character_
  out %>%
    mutate(
      Group = factor(Group, levels = unique(c(group_levels, sort(unique(na.omit(Group)))))),
      Sex = factor(Sex, levels = unique(c(sex_levels, sort(unique(na.omit(Sex))))))
    )
}

make_feature_name <- function(source, domain, scale, metric, statistic, context = "global") {
  paste(safe_name(c(source, domain, scale, metric, statistic, context)), collapse = "__")
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
      q25 = quantile(Value, 0.25, na.rm = TRUE, names = FALSE),
      q75 = quantile(Value, 0.75, na.rm = TRUE, names = FALSE),
      p95 = quantile(Value, 0.95, na.rm = TRUE, names = FALSE),
      max = max(Value, na.rm = TRUE),
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

first_active <- base %>%
  filter(str_detect(str_to_lower(Phase), "active|dark|night")) %>%
  group_by(AnimalNum, CageChange, Phase) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  mutate(local_bin = row_number()) %>%
  filter(local_bin <= ifelse(primary_bin_level == "5min_based", 144, 72)) %>%
  ungroup()

base_first_active_features <- calc_feature_set(first_active, "first_active_12h")

core_feature_long <- bind_rows(base_global_features, base_phase_features, base_first_active_features) %>%
  filter(is.finite(FeatureValue))

core_feature_wide <- core_feature_long %>%
  group_by(AnimalNum, Group, Sex, feature) %>%
  summarise(FeatureValue = mean(FeatureValue, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = feature, values_from = FeatureValue)

# Add interpretable composite axes from early active window.
core_feature_wide <- core_feature_wide %>%
  mutate(
    systems__composite__early__social_withdrawal__z__first_active_12h =
      safe_scale(raw__behavior__5min_based__movement__mean__first_active_12h %||% NA_real_) -
      safe_scale(raw__behavior__5min_based__proximity__mean__first_active_12h %||% NA_real_),
    systems__composite__early__behavioral_instability__z__first_active_12h =
      rowMeans(cbind(
        safe_scale(raw__behavior__5min_based__movement__rmssd__first_active_12h %||% NA_real_),
        safe_scale(raw__behavior__5min_based__entropy__rmssd__first_active_12h %||% NA_real_),
        safe_scale(raw__behavior__5min_based__proximity__rmssd__first_active_12h %||% NA_real_)
      ), na.rm = TRUE),
    systems__composite__early__behavioral_inertia__z__first_active_12h =
      rowMeans(cbind(
        safe_scale(raw__behavior__5min_based__movement__acf1__first_active_12h %||% NA_real_),
        safe_scale(raw__behavior__5min_based__entropy__acf1__first_active_12h %||% NA_real_),
        safe_scale(raw__behavior__5min_based__proximity__acf1__first_active_12h %||% NA_real_)
      ), na.rm = TRUE)
  )

# ------------------------------------------------
# OPTIONAL MODULE OUTPUTS FROM SPECIALIZED SCRIPTS
# ------------------------------------------------

load_optional_animal_table <- function(path, source_label, domain_label, scale_label = primary_bin_level) {
  dat <- read_any_table(path)
  if (is.null(dat) || nrow(dat) == 0) return(tibble())
  dat <- standardize_id_columns(dat)

  numeric_cols <- dat %>%
    select(where(is.numeric)) %>%
    names()
  numeric_cols <- setdiff(numeric_cols, c("AnimalNum", "TimeIndex", "BinSizeSec"))
  if (length(numeric_cols) == 0) return(tibble())

  dat %>%
    group_by(AnimalNum, Group, Sex) %>%
    summarise(across(all_of(numeric_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_longer(all_of(numeric_cols), names_to = "Metric", values_to = "FeatureValue") %>%
    filter(is.finite(FeatureValue)) %>%
    mutate(feature = make_feature_name(source_label, domain_label, scale_label, Metric, "mean", "script_output")) %>%
    select(AnimalNum, Group, Sex, feature, FeatureValue)
}

optional_files <- tibble(
  source_label = c("burstiness", "state_space", "early_prediction", "dynamic_network", "hmm", "gamm"),
  domain_label = c("temporal_dynamics", "latent_space", "prediction", "social_network", "behavioral_state", "trajectory"),
  path = c(
    file.path(paths$Path[paths$Source == "burstiness_instability"], "burstiness_instability_features.csv"),
    file.path(paths$Path[paths$Source == "state_space"], "behavioral_state_space_features.csv"),
    file.path(paths$Path[paths$Source == "early_prediction"], "early_behavior_features.csv"),
    file.path(paths$Path[paths$Source == "dynamic_networks"], "dynamic_social_network_features.csv"),
    file.path(paths$Path[paths$Source == "hmm_states"], "hmm_animal_state_features.csv"),
    file.path(paths$Path[paths$Source == "gamm_trajectory"], "gamm_trajectory_features.csv")
  )
)

optional_long <- pmap_dfr(optional_files, function(source_label, domain_label, path) {
  load_optional_animal_table(path, source_label, domain_label)
})

optional_wide <- optional_long %>%
  group_by(AnimalNum, Group, Sex, feature) %>%
  summarise(FeatureValue = mean(FeatureValue, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = feature, values_from = FeatureValue)

systems_features <- full_join(core_feature_wide, optional_wide, by = c("AnimalNum", "Group", "Sex")) %>%
  mutate(
    Group = factor(as.character(Group), levels = group_levels),
    Sex = factor(as.character(Sex), levels = unique(c(sex_levels, sort(unique(as.character(Sex))))))
  )

# ------------------------------------------------
# OPTIONAL ENDPOINTS AND PROTEOMICS MODULES
# ------------------------------------------------

endpoint_dat <- read_any_table(endpoint_file)
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

write_table(systems_features, file.path(output_dir, "tables/systems_animal_feature_matrix.csv"))
write_table(feature_dictionary, file.path(output_dir, "tables/systems_feature_dictionary.csv"))
write_table(feature_qc, file.path(output_dir, "tables/systems_feature_qc.csv"))

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
    sem = sd / sqrt(n),
    median = median(Value, na.rm = TRUE),
    q25 = quantile(Value, 0.25, na.rm = TRUE, names = FALSE),
    q75 = quantile(Value, 0.75, na.rm = TRUE, names = FALSE),
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
                      estimate = NA_real_, p.value = NA_real_, cohen_d = NA_real_, status = "low_n"))
      }
      tt <- try(t.test(y, x), silent = TRUE)
      tibble(
        contrast = paste0(comp, "-", ref),
        n_ref = n_ref,
        n_comp = n_comp,
        estimate = mean(y, na.rm = TRUE) - mean(x, na.rm = TRUE),
        p.value = if (inherits(tt, "try-error")) NA_real_ else tt$p.value,
        cohen_d = cohens_d_pooled(x, y),
        status = "tested"
      )
    })
  }) %>%
  ungroup() %>%
  group_by(Sex, Source, Domain, Scale, Context, contrast) %>%
  mutate(
    p_fdr = p.adjust(p.value, method = "BH"),
    sig = sig_label(p_fdr),
    evidence = case_when(
      is.na(p_fdr) ~ "not_tested",
      p_fdr < 0.05 & abs(cohen_d) >= 0.8 ~ "large_FDR_supported",
      p_fdr < 0.05 ~ "FDR_supported",
      p_fdr < 0.10 ~ "trend",
      abs(cohen_d) >= 0.8 ~ "large_effect_uncertain",
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
x_scaled[, !is.finite(colSums(x_scaled)), drop = FALSE] <- 0

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

p_pca <- pca_scores %>%
  ggplot(aes(PC1, PC2, colour = Group, fill = Group, shape = Sex)) +
  geom_hline(yintercept = 0, linewidth = 0.15, colour = "grey85") +
  geom_vline(xintercept = 0, linewidth = 0.15, colour = "grey85") +
  stat_ellipse(aes(group = interaction(Group, Sex)), linewidth = 0.35, alpha = 0.45, show.legend = FALSE) +
  geom_point(size = 2.2, alpha = 0.85, stroke = 0.25) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  labs(
    title = "Integrated behavioral state space",
    subtitle = paste0("PC1 ", round(100 * var_exp[1], 1), "%; PC2 ", round(100 * var_exp[2], 1), "% variance"),
    x = "Systems PC1",
    y = "Systems PC2"
  ) +
  make_nature_theme(base_size = 7)

save_plot_svg_pdf(p_pca, file.path(output_dir, "figures/publication_panels/Fig_systems_state_space_PCA"), width = 90, height = 75)

if (requireNamespace("uwot", quietly = TRUE) && nrow(x_scaled) >= 8) {
  set.seed(1)
  n_neighbors <- max(3, min(10, floor(nrow(x_scaled) / 2)))
  um <- uwot::umap(x_scaled, n_neighbors = n_neighbors, min_dist = 0.25, metric = "euclidean")
  umap_scores <- systems_features %>%
    select(AnimalNum, Group, Sex, any_of(endpoint_cols)) %>%
    bind_cols(tibble(UMAP1 = um[, 1], UMAP2 = um[, 2]))
  write_table(umap_scores, file.path(output_dir, "tables/systems_umap_scores.csv"))

  p_umap <- umap_scores %>%
    ggplot(aes(UMAP1, UMAP2, colour = Group, fill = Group, shape = Sex)) +
    geom_point(size = 2.3, alpha = 0.88, stroke = 0.25) +
    stat_ellipse(aes(group = interaction(Group, Sex)), linewidth = 0.35, alpha = 0.45, show.legend = FALSE) +
    scale_colour_manual(values = group_colors, drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    labs(
      title = "Nonlinear behavioral manifold",
      subtitle = paste0("UMAP on integrated animal-level systems features; n_neighbors=", n_neighbors),
      x = "UMAP1",
      y = "UMAP2"
    ) +
    make_nature_theme(base_size = 7)

  save_plot_svg_pdf(p_umap, file.path(output_dir, "figures/publication_panels/Fig_systems_state_space_UMAP"), width = 90, height = 75)
}

# ------------------------------------------------
# EFFECT-SIZE HEATMAPS
# ------------------------------------------------

focus_domains <- c("behavior", "composite", "temporal_dynamics", "latent_space", "behavioral_state", "trajectory")
heat_tbl <- group_contrasts %>%
  filter(status == "tested", Domain %in% focus_domains) %>%
  mutate(
    DisplayFeature = paste(Metric, Statistic, Context, sep = " | "),
    DisplayFeature = str_replace_all(DisplayFeature, "_", " "),
    DisplayFeature = str_trunc(DisplayFeature, width = 48),
    contrast = factor(contrast, levels = c("RES-CON", "SUS-CON", "SUS-RES"))
  )

if (publication_focus_only) {
  keep_features <- heat_tbl %>%
    group_by(feature) %>%
    summarise(max_abs_d = max(abs(cohen_d), na.rm = TRUE), min_fdr = min(p_fdr, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_abs_d), min_fdr) %>%
    slice_head(n = 40) %>%
    pull(feature)
  heat_tbl <- heat_tbl %>% filter(feature %in% keep_features)
}

p_heat <- heat_tbl %>%
  mutate(DisplayFeature = factor(DisplayFeature, levels = rev(unique(DisplayFeature)))) %>%
  ggplot(aes(contrast, DisplayFeature, fill = cohen_d)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = sig), size = 2.1) +
  facet_grid(Sex ~ Domain, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#3d3b6e", mid = "white", high = "#e63947", midpoint = 0, na.value = "grey90") +
  labs(
    title = "Systems-level group-difference map",
    subtitle = "Cohen's d; symbols denote FDR within sex/source/domain/context/contrast",
    x = NULL,
    y = NULL,
    fill = "Cohen's d"
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
  filter(status == "tested", !is.na(cohen_d)) %>%
  group_by(feature) %>%
  summarise(max_abs_d = max(abs(cohen_d), na.rm = TRUE), .groups = "drop") %>%
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
      Pair = paste(str_trunc(Metric1, 18), str_trunc(Metric2, 18), sep = " ↔ "),
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
  outcome_assoc <- map_dfr(available_outcomes, function(outcome) {
    y <- safe_numeric(systems_features[[outcome]])
    map_dfr(usable_features, function(fc) {
      tibble(
        outcome = outcome,
        feature = fc,
        n = sum(is.finite(y) & is.finite(systems_features[[fc]])),
        spearman_rho = safe_cor(systems_features[[fc]], y, "spearman"),
        spearman_p = safe_cor_p(systems_features[[fc]], y, "spearman"),
        pearson_r = safe_cor(systems_features[[fc]], y, "pearson"),
        pearson_p = safe_cor_p(systems_features[[fc]], y, "pearson")
      )
    })
  }) %>%
    group_by(outcome) %>%
    mutate(
      spearman_fdr = p.adjust(spearman_p, method = "BH"),
      sig = sig_label(spearman_fdr)
    ) %>%
    ungroup() %>%
    left_join(feature_dictionary, by = "feature") %>%
    arrange(outcome, spearman_fdr, desc(abs(spearman_rho)))

  write_table(outcome_assoc, file.path(output_dir, "stats_tables/systems_outcome_associations.csv"))

  outcome_to_plot <- if (primary_outcome %in% available_outcomes) primary_outcome else available_outcomes[1]
  top_outcome_features <- outcome_assoc %>%
    filter(outcome == outcome_to_plot, !is.na(spearman_rho)) %>%
    arrange(spearman_fdr, desc(abs(spearman_rho))) %>%
    slice_head(n = 25)

  p_out_heat <- top_outcome_features %>%
    mutate(DisplayFeature = paste(Metric, Statistic, Context, sep = " | "),
           DisplayFeature = str_replace_all(DisplayFeature, "_", " "),
           DisplayFeature = factor(str_trunc(DisplayFeature, 52), levels = rev(str_trunc(DisplayFeature, 52)))) %>%
    ggplot(aes(spearman_rho, DisplayFeature, fill = spearman_rho)) +
    geom_col(width = 0.72) +
    geom_text(aes(label = sig), hjust = -0.15, size = 2.1) +
    scale_fill_gradient2(low = "#3d3b6e", mid = "white", high = "#e63947", midpoint = 0) +
    labs(
      title = paste0("Early systems features linked to ", outcome_to_plot),
      subtitle = "Top Spearman associations; symbols denote FDR across integrated features",
      x = "Spearman rho",
      y = NULL,
      fill = "rho"
    ) +
    coord_cartesian(xlim = c(-1, 1), clip = "off") +
    make_nature_theme(base_size = 6) +
    theme(legend.position = "right", plot.margin = margin(5.5, 18, 5.5, 5.5))

  save_plot_svg_pdf(p_out_heat, file.path(output_dir, "figures/publication_panels/Fig_systems_outcome_association_rank"), width = 135, height = 115)

  if (requireNamespace("glmnet", quietly = TRUE)) {
    model_dat <- systems_features %>%
      filter(is.finite(safe_numeric(.data[[outcome_to_plot]])))
    model_features <- usable_features[usable_features %in% names(model_dat)]
    x_model <- model_dat %>%
      select(all_of(model_features)) %>%
      mutate(across(everything(), ~ replace_na(.x, median(.x, na.rm = TRUE)))) %>%
      as.matrix()
    y_model <- safe_numeric(model_dat[[outcome_to_plot]])

    if (nrow(x_model) >= 8 && ncol(x_model) >= 2) {
      pred <- rep(NA_real_, length(y_model))
      coef_tbl <- tibble()
      set.seed(1)
      for (i in seq_along(y_model)) {
        train_idx <- setdiff(seq_along(y_model), i)
        fit <- glmnet::cv.glmnet(
          x_model[train_idx, , drop = FALSE], y_model[train_idx],
          alpha = 0.5,
          standardize = TRUE,
          nfolds = min(5, length(train_idx))
        )
        pred[i] <- as.numeric(predict(fit, x_model[i, , drop = FALSE], s = "lambda.min"))
        cc <- as.matrix(coef(fit, s = "lambda.min"))[, 1]
        coef_tbl <- bind_rows(coef_tbl, tibble(held_out = model_dat$AnimalNum[i], feature = names(cc), coefficient = as.numeric(cc)))
      }

      prediction_tbl <- model_dat %>%
        transmute(AnimalNum, Group, Sex, outcome = y_model, predicted = pred, residual = outcome - predicted)
      baseline <- mean(prediction_tbl$outcome, na.rm = TRUE)
      prediction_perf <- tibble(
        outcome = outcome_to_plot,
        n = nrow(prediction_tbl),
        loo_pearson_r = safe_cor(prediction_tbl$outcome, prediction_tbl$predicted, "pearson"),
        loo_spearman_rho = safe_cor(prediction_tbl$outcome, prediction_tbl$predicted, "spearman"),
        loo_rmse = sqrt(mean((prediction_tbl$outcome - prediction_tbl$predicted)^2, na.rm = TRUE)),
        baseline_rmse = sqrt(mean((prediction_tbl$outcome - baseline)^2, na.rm = TRUE)),
        cv_r2_vs_mean = 1 - sum((prediction_tbl$outcome - prediction_tbl$predicted)^2, na.rm = TRUE) /
          sum((prediction_tbl$outcome - baseline)^2, na.rm = TRUE)
      )
      coef_summary <- coef_tbl %>%
        filter(feature != "(Intercept)") %>%
        group_by(feature) %>%
        summarise(nonzero_frequency = mean(coefficient != 0), median_coefficient = median(coefficient), .groups = "drop") %>%
        left_join(feature_dictionary, by = "feature") %>%
        arrange(desc(nonzero_frequency), desc(abs(median_coefficient)))

      write_table(prediction_tbl, file.path(output_dir, "tables/systems_outcome_loo_predictions.csv"))
      write_table(prediction_perf, file.path(output_dir, "tables/systems_outcome_loo_performance.csv"))
      write_table(coef_summary, file.path(output_dir, "tables/systems_outcome_elastic_net_feature_stability.csv"))

      p_pred <- prediction_tbl %>%
        ggplot(aes(outcome, predicted, colour = Group, fill = Group, shape = Sex)) +
        geom_abline(slope = 1, intercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey45") +
        geom_point(size = 2.2, alpha = 0.9) +
        geom_smooth(method = "lm", se = TRUE, linewidth = 0.45, colour = "grey25") +
        scale_colour_manual(values = group_colors, drop = FALSE) +
        scale_fill_manual(values = group_colors, drop = FALSE) +
        labs(
          title = paste0("Integrated early behavior predicts ", outcome_to_plot),
          subtitle = paste0("LOO elastic-net: r=", round(prediction_perf$loo_pearson_r, 2), "; CV R2=", round(prediction_perf$cv_r2_vs_mean, 2)),
          x = paste0("Observed ", outcome_to_plot),
          y = paste0("Predicted ", outcome_to_plot)
        ) +
        make_nature_theme(base_size = 7)

      save_plot_svg_pdf(p_pred, file.path(output_dir, "figures/publication_panels/Fig_systems_crossvalidated_prediction"), width = 90, height = 75)
    }
  }
}

# ------------------------------------------------
# PUBLICATION DASHBOARD
# ------------------------------------------------

# A compact, paper-facing composite figure. If patchwork is not installed,
# individual panels above remain the primary outputs.
if (requireNamespace("patchwork", quietly = TRUE)) {
  top_loadings <- pca_loadings %>%
    arrange(desc(abs(PC1))) %>%
    slice_head(n = 12) %>%
    mutate(DisplayFeature = paste(Metric, Statistic, Context, sep = " | "),
           DisplayFeature = str_replace_all(DisplayFeature, "_", " "),
           DisplayFeature = factor(str_trunc(DisplayFeature, 45), levels = rev(str_trunc(DisplayFeature, 45))))

  p_load <- top_loadings %>%
    ggplot(aes(PC1, DisplayFeature, fill = PC1)) +
    geom_col(width = 0.72) +
    scale_fill_gradient2(low = "#3d3b6e", mid = "white", high = "#e63947", midpoint = 0) +
    labs(title = "Dominant systems axis", x = "PC1 loading", y = NULL) +
    make_nature_theme(base_size = 6) +
    theme(legend.position = "none")

  p_heat_small <- heat_tbl %>%
    group_by(Sex, contrast) %>%
    slice_max(abs(cohen_d), n = 8, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(DisplayFeature = paste(Metric, Statistic, Context, sep = " | "),
           DisplayFeature = str_replace_all(DisplayFeature, "_", " "),
           DisplayFeature = factor(str_trunc(DisplayFeature, 40), levels = rev(unique(str_trunc(DisplayFeature, 40))))) %>%
    ggplot(aes(contrast, DisplayFeature, fill = cohen_d)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(aes(label = sig), size = 1.8) +
    facet_grid(Sex ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = "#3d3b6e", mid = "white", high = "#e63947", midpoint = 0, na.value = "grey90") +
    labs(title = "Group effects", x = NULL, y = NULL, fill = "d") +
    make_nature_theme(base_size = 5.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

  dashboard <- (p_pca | p_load) / p_heat_small +
    patchwork::plot_annotation(
      title = "Multiscale behavioral systems architecture of stress adaptation",
      subtitle = "Animal-level integration of movement, entropy, social proximity, temporal instability, latent state-space and optional prediction modules"
    )

  save_plot_svg_pdf(dashboard, file.path(output_dir, "figures/Fig_integrated_systems_dashboard"), width = 180, height = 180)
}

# ------------------------------------------------
# TEXT SUMMARY FOR MANUSCRIPT / LAB MEETING
# ------------------------------------------------

strong_findings <- group_contrasts %>%
  filter(evidence %in% c("large_FDR_supported", "FDR_supported", "large_effect_uncertain")) %>%
  arrange(match(evidence, c("large_FDR_supported", "FDR_supported", "large_effect_uncertain")), desc(abs(cohen_d))) %>%
  select(Sex, contrast, Source, Domain, Metric, Statistic, Context, cohen_d, p.value, p_fdr, evidence) %>%
  slice_head(n = 50)

systems_summary <- tibble(
  Item = c(
    "Animals included",
    "Usable integrated features",
    "Primary bin level",
    "Main feature domains",
    "Interpretation"
  ),
  Value = c(
    as.character(n_distinct(systems_features$AnimalNum)),
    as.character(length(usable_features)),
    primary_bin_level,
    paste(sort(unique(feature_dictionary$Domain)), collapse = "; "),
    "Use PCA/UMAP for system-level organization, effect-size heatmaps for interpretable sex-by-group contrasts, and outcome associations/prediction only as endpoint-linked evidence."
  )
)

write_table(strong_findings, file.path(output_dir, "tables/systems_top_findings_for_reporting.csv"))
write_table(systems_summary, file.path(output_dir, "tables/systems_analysis_summary.csv"))

message("Integrated systems neuroscience summary complete.")
message("Output: ", output_dir)
message("Primary figure candidates:")
message("  - figures/publication_panels/Fig_systems_state_space_PCA.svg")
message("  - figures/publication_panels/Fig_systems_effect_size_heatmap.svg")
message("  - figures/Fig_integrated_systems_dashboard.svg")
