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
#     08b_early_prediction_model_ladder.R
#     09_dynamic_social_networks.R
#     10_hmm_behavioral_states.R
#     11_gamm_trajectory_features.R
#     13_nonlinear_systems_dynamics.R
#     14_nextgen_behavioral_phenotyping.R
#
# Output:
#   - animal-level multiscale systems feature matrix
#   - feature dictionary and QC tables
#   - CON/RES/SUS group summaries and contrasts by sex
#   - PCA/UMAP behavioral state-space plots
#   - systems effect-size heatmaps
#   - sex-specific feature correlation networks
#   - early-prediction summary panels if outcome data are available
#   - module-level prediction ladder with incremental CV-R2
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
optional_import_bin_levels <- unique(c(sensitivity_bin_levels, "1min_based"))

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
n_prediction_permutations <- 1000
n_prediction_bootstrap <- 1000

# Optional proteomics module score file. Same expectation: one row per animal.
# Useful columns could include RNP_module, Mito_module, Translation_module, etc.
proteomics_module_file <- NULL
preprocessed_position_dir <- file.path(project_root, "MMMSociability/preprocessed_data")

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
source_if_exists(file.path(repo_root, "Functions/duration_normalization_helpers.R"))

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
ensure_dir(file.path(output_dir, "tables/qc"))
ensure_dir(file.path(output_dir, "stats_tables"))
ensure_dir(file.path(output_dir, "figures"))
ensure_dir(file.path(output_dir, "figures/publication_panels"))
ensure_dir(file.path(output_dir, "figures/supplementary"))
ensure_dir(file.path(output_dir, "figures/exploratory"))
ensure_dir(file.path(output_dir, "figures/qc"))
output_dirs <- analysis_output_dirs(output_dir)
write_output_manifest(
  output_dir,
  script_name = "12_systems_neuroscience_summary.R",
  analysis_name = "integrated systems neuroscience summary",
  primary_tables = c(
    "tables/systems_animal_feature_matrix.csv",
    "tables/systems_feature_dictionary.csv",
    "tables/systems_computation_integration_audit.csv",
    "tables/systems_robustness_audit.csv",
    "tables/qc_chip_loss_flags.csv",
    "tables/systems_batch_system_cage_audit.csv",
    "tables/systems_primary_feature_set.csv",
    "tables/systems_primary_claims_table.csv",
    "tables/systems_sleep_like_inactivity_summary.csv",
    "tables/systems_endpoint_leakage_audit.csv",
    "tables/systems_leave_one_context_out_robustness.csv",
    "tables/systems_bin_size_sensitivity_summary.csv",
    "tables/systems_behavior_proteomics_bridge.csv",
    "stats_tables/systems_group_sex_interaction_models.csv",
    "tables/qc/chip_loss_downstream_feature_audit.csv",
    "tables/qc/group_balance_by_batch_system.csv",
    "tables/systems_module_scorecards.csv",
    "tables/systems_named_biological_scores.csv",
    "tables/systems_claim_hierarchy.csv",
    "tables/systems_primary_findings_summary.csv",
    "tables/systems_robustness_summary.csv",
    "tables/systems_interpretation_guide.csv",
    "tables/systems_visualization_guide.csv",
    "stats_tables/systems_group_contrasts.csv",
    "tables/duration_sensitivity_audit.csv"
  ),
  primary_figures = c(
    "figures/publication_panels/Fig_systems_state_space_PCA.svg",
    "figures/publication_panels/Fig_systems_effect_size_heatmap.svg",
    "figures/publication_panels/Fig_systems_module_scorecard.svg",
    "figures/publication_panels/Fig_systems_named_biological_scores.svg",
    "figures/publication_panels/Fig_sleep_like_inactivity_by_group_sex.svg",
    "figures/publication_panels/Fig_group_sex_interaction_effects.svg",
    "figures/publication_panels/Fig_first_active_movement_trajectory_by_group_sex.svg",
    "figures/publication_panels/Fig_first_active_entropy_acf1_or_instability.svg",
    "figures/publication_panels/Fig_behavior_proteomics_bridge.svg",
    "figures/qc/Fig_chip_loss_dropout_timeline.svg",
    "figures/qc/Fig_chip_loss_movement_proximity_diagnostics.svg",
    "figures/qc/Fig_batch_system_feature_bias.svg",
    "figures/qc/Fig_group_balance_by_batch_system.svg",
    "figures/qc/Fig_primary_feature_robustness.svg",
    "figures/publication_panels/Fig_systems_prediction_ladder.svg"
  ),
  notes = c("This is the publication-facing integration layer; use systems_visualization_guide.csv to choose panels.")
)

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

collapse_to_animal_rows <- function(dat) {
  if (is.null(dat) || nrow(dat) == 0 || !"AnimalNum" %in% names(dat)) return(dat)
  dat %>%
    group_by(AnimalNum) %>%
    summarise(
      across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
      across(
        where(~ !is.numeric(.x)),
        ~ {
          vals <- .x[!is.na(.x) & as.character(.x) != ""]
          if (length(vals) == 0) NA else vals[1]
        }
      ),
      .groups = "drop"
    ) %>%
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.x), NA_real_, .x)))
}

first_optional_col <- function(dat, candidates) {
  hit <- candidates[candidates %in% names(dat)][1]
  if (length(hit) == 0 || is.na(hit)) NA_character_ else hit
}

safe_model_p <- function(formula, dat, term_pattern) {
  fit <- try(lm(formula, data = dat), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  cf <- try(coef(summary(fit)), silent = TRUE)
  if (inherits(cf, "try-error") || is.null(rownames(cf))) return(NA_real_)
  rows <- str_detect(rownames(cf), term_pattern)
  if (!any(rows)) return(NA_real_)
  min(cf[rows, "Pr(>|t|)"], na.rm = TRUE)
}

bootstrap_cor_ci <- function(x, y, method = "spearman", n_boot = 500, seed = 1) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 6 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(c(low = NA_real_, high = NA_real_))
  set.seed(seed)
  xv <- x[ok]
  yv <- y[ok]
  vals <- replicate(n_boot, {
    idx <- sample(seq_along(xv), replace = TRUE)
    suppressWarnings(cor(xv[idx], yv[idx], method = method))
  })
  quantile(vals[is.finite(vals)], c(0.025, 0.975), na.rm = TRUE, names = FALSE) %>%
    setNames(c("low", "high"))
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

feature_entropy <- function(p) {
  p <- p[is.finite(p) & p > 0]
  if (length(p) == 0) return(NA_real_)
  -sum(p * log(p))
}

dominant_state_by_profile <- function(state_summary, prefer = c("inactive", "explore", "burst", "social")) {
  prefer <- match.arg(prefer)
  if (is.null(state_summary) || nrow(state_summary) == 0 || !"State" %in% names(state_summary)) return(NA_character_)
  prof <- state_summary %>%
    mutate(
      Movement_z = safe_numeric(.data[["Movement_z"]]),
      Entropy_z = safe_numeric(.data[["Entropy_z"]]),
      Proximity_z = safe_numeric(.data[["Proximity_z"]]),
      state_score = if (prefer == "inactive") {
        -Movement_z - Entropy_z
      } else if (prefer == "explore") {
        Entropy_z + Movement_z - Proximity_z
      } else if (prefer == "burst") {
        Movement_z
      } else if (prefer == "social") {
        Proximity_z
      } else {
        NA_real_
      }
    ) %>%
    filter(is.finite(state_score))
  if (nrow(prof) == 0) return(NA_character_)
  as.character(prof$State[which.max(prof$state_score)][1])
}

feature_module_from_parts <- function(source, domain, metric, statistic, context, feature) {
  key <- str_to_lower(paste(source, domain, metric, statistic, context, feature, sep = " "))
  case_when(
    str_detect(key, "biological_scores|rigidity_score|flexibility_score|withdrawal_score|fragmentation_score|recovery_score|adaptation_index|resilience_score") ~ "Predictive systems integration",
    str_detect(key, "prediction|integrated_systems_phenotype") ~ "Predictive systems integration",
    str_detect(key, "attractor|recurrence|manifold|multiscale|mse|sample_entropy|permutation_entropy|early_warning|criticalslowing|flickering|energy_landscape|complexity") ~ "Nonlinear systems dynamics",
    str_detect(key, "gamm|trajectory|auc|peak|trough|dynamic_range|time_to_peak|rebound_slope|half_life|decay|curvature|asymmetry|recovery") ~ "Trajectory geometry",
    str_detect(key, "dynamic_network|social_network|dyadic|centrality|degree|strength|betweenness|closeness|fragmentation|modularity|clustering|components|contact_entropy|isolation|engagement") ~ "Social topology",
    str_detect(key, "hmm|latent_state|state_occupancy|dwell|transition|switch_rate|state_bias|rigidity|flexibility|social_state_fraction|burst_state_fraction") ~ "Latent-state organization",
    str_detect(key, "rmssd|acf1|fano|cv|instability|inertia|persistence|burstiness") ~ "Temporal instability",
    str_detect(key, "movement.*mean|entropy.*mean|proximity.*mean|mean_locomotor|mean_spatial|mean_social") ~ "Magnitude",
    TRUE ~ "Other interpretable feature"
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
    "gamm_trajectory",
    "nonlinear_systems_dynamics",
    "adaptation_kinetics",
    "sleep_like_inactivity",
    "phase_organization",
    "behavior_proteomics",
    "nextgen_selective"
  ),
  Path = c(
    file.path(project_root, "analysis_ready/03_derived_metrics", primary_bin_level, "all_behavior_metrics.csv"),
    file.path(project_root, "analysis_ready/06_behavioral_dynamics/burstiness", primary_bin_level, "tables"),
    file.path(project_root, "analysis_ready/06_behavioral_dynamics/state_space", primary_bin_level, "tables"),
    file.path(project_root, "analysis_ready/06_behavioral_dynamics/early_prediction", primary_bin_level, "tables"),
    file.path(project_root, "analysis_ready/06_behavioral_dynamics/social_networks", primary_bin_level, "tables"),
    file.path(project_root, "analysis_ready/06_behavioral_dynamics/hmm_states", primary_bin_level, "tables"),
    file.path(project_root, "analysis_ready/06_behavioral_dynamics/gamm_trajectory_features", primary_bin_level, "tables"),
    file.path(project_root, "analysis_ready/13_nonlinear_systems_dynamics", primary_bin_level, "derived_data"),
    file.path(project_root, "analysis_ready/15_behavioral_adaptation_kinetics", primary_bin_level, "tables"),
    file.path(project_root, "analysis_ready/16_sleep_like_inactivity_metrics", primary_bin_level, "tables"),
    file.path(project_root, "analysis_ready/17_ethological_phase_organization", primary_bin_level, "tables"),
    file.path(project_root, "analysis_ready/12_behavior_proteomics_integration", "tables"),
    file.path(project_root, "analysis_ready/14_nextgen_behavioral_phenotyping", primary_bin_level, "tables")
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

epoch_duration_qc <- if (exists("write_epoch_duration_qc")) {
  write_epoch_duration_qc(base, output_dir, metric_source = "12_systems_neuroscience_summary", bin_size_sec = infer_bin_size_sec(base))
} else {
  tibble()
}

# ------------------------------------------------
# RAW/PREPROCESSED RFID CHIP-LOSS / DEAD-TAG QC
# ------------------------------------------------

read_preprocessed_position_files <- function(position_dir, max_files = Inf) {
  if (is.null(position_dir) || !dir.exists(position_dir)) return(tibble())
  files <- list.files(position_dir, pattern = "preprocessed.*\\.csv$", full.names = TRUE, ignore.case = TRUE)
  if (length(files) == 0) return(tibble())
  files <- files[seq_len(min(length(files), max_files))]
  map_dfr(files, function(path) {
    dat <- try(readr::read_csv(path, show_col_types = FALSE, progress = FALSE), silent = TRUE)
    if (inherits(dat, "try-error") || nrow(dat) == 0) return(tibble())
    dat %>% mutate(SourceFile = basename(path))
  })
}

raw_position_qc <- read_preprocessed_position_files(preprocessed_position_dir)

qc_chip_loss_flags <- if (nrow(raw_position_qc) > 0) {
  raw_id_col <- first_existing_col(raw_position_qc, c("AnimalNum", "AnimalID", "Animal", "MouseID", "RFID"), TRUE, "raw animal ID")
  raw_time_col <- first_existing_col(raw_position_qc, c("DateTime", "Timestamp", "Time", "BinStart"), TRUE, "raw timestamp")
  raw_position_col <- first_existing_col(raw_position_qc, c("PositionID", "Position", "Antenna", "AntennaID"), FALSE, "raw position")
  raw_batch_col <- first_existing_col(raw_position_qc, c("Batch", "batch"), FALSE, "batch")
  raw_system_col <- first_existing_col(raw_position_qc, c("System", "system", "Cage", "CageID"), FALSE, "system")
  raw_cage_col <- first_existing_col(raw_position_qc, c("CageChange", "CC", "Regrouping"), FALSE, "cage-change")
  raw_phase_col <- first_existing_col(raw_position_qc, c("Phase", "phase"), FALSE, "phase")

  raw_norm <- raw_position_qc %>%
    transmute(
      AnimalNum = as.character(.data[[raw_id_col]]) %>% str_trim() %>% str_replace_all("\\s+", "") %>% str_to_upper(),
      Time = suppressWarnings(as.POSIXct(.data[[raw_time_col]], tz = "UTC")),
      PositionID = if (!is.na(raw_position_col)) as.character(.data[[raw_position_col]]) else NA_character_,
      Batch = if (!is.na(raw_batch_col)) as.character(.data[[raw_batch_col]]) else NA_character_,
      System = if (!is.na(raw_system_col)) as.character(.data[[raw_system_col]]) else NA_character_,
      CageChange = if (!is.na(raw_cage_col)) as.character(.data[[raw_cage_col]]) else NA_character_,
      Phase = if (!is.na(raw_phase_col)) as.character(.data[[raw_phase_col]]) else NA_character_
    ) %>%
    filter(!is.na(AnimalNum), !is.na(Time)) %>%
    arrange(AnimalNum, Time)

  raw_epoch <- raw_norm %>%
    group_by(AnimalNum, Batch, System, CageChange, Phase) %>%
    arrange(Time, .by_group = TRUE) %>%
    summarise(
      first_time = min(Time, na.rm = TRUE),
      last_time = max(Time, na.rm = TRUE),
      n_reads = n(),
      n_positions = n_distinct(PositionID[!is.na(PositionID)]),
      position_switch_rate = if_else(n() >= 3, mean(PositionID != lag(PositionID), na.rm = TRUE), NA_real_),
      longest_gap_hours = {
        gaps <- as.numeric(diff(Time), units = "hours")
        if (length(gaps) == 0) NA_real_ else max(gaps, na.rm = TRUE)
      },
      .groups = "drop"
    ) %>%
    mutate(
      expected_reads = pmax(1, ceiling(as.numeric(difftime(last_time, first_time, units = "mins")) / pmax(infer_bin_size_sec(base) / 60, 1)) + 1),
      raw_observed_fraction = pmin(n_reads / expected_reads, 1)
    )

  base_epoch <- base %>%
    group_by(AnimalNum, Group, Sex, CageChange, Phase) %>%
    summarise(
      movement_epoch_mean = mean(Movement, na.rm = TRUE),
      proximity_epoch_mean = mean(Proximity, na.rm = TRUE),
      observed_bins_base = n(),
      .groups = "drop"
    )

  raw_epoch %>%
    left_join(base_epoch, by = c("AnimalNum", "CageChange", "Phase")) %>%
    group_by(AnimalNum) %>%
    arrange(first_time, .by_group = TRUE) %>%
    mutate(
      baseline_movement = median(movement_epoch_mean[row_number() <= ceiling(n() / 2)], na.rm = TRUE),
      baseline_proximity = median(proximity_epoch_mean[row_number() <= ceiling(n() / 2)], na.rm = TRUE),
      sudden_sustained_loss = raw_observed_fraction < 0.55,
      near_zero_movement_after_normal = is.finite(baseline_movement) & baseline_movement > 0.25 & movement_epoch_mean <= pmax(0.05, 0.15 * baseline_movement),
      near_zero_proximity_after_normal = is.finite(baseline_proximity) & baseline_proximity > 0.02 & proximity_epoch_mean <= pmax(0.005, 0.15 * baseline_proximity),
      long_flatline_immobile = n_positions <= 1 & position_switch_rate <= 0.01 & movement_epoch_mean <= 0.05,
      expected_observed_mismatch = raw_observed_fraction < 0.70,
      suspected_epoch = sudden_sustained_loss |
        long_flatline_immobile |
        (expected_observed_mismatch & (near_zero_movement_after_normal | near_zero_proximity_after_normal)),
      first_suspected_dropout_time = if (any(suspected_epoch %in% TRUE)) min(first_time[suspected_epoch %in% TRUE], na.rm = TRUE) else as.POSIXct(NA),
      last_valid_time = if_else(!is.na(first_suspected_dropout_time), first_suspected_dropout_time, last_time),
      qc_epoch_class = case_when(
        raw_observed_fraction < 0.10 | observed_bins_base < 4 ~ "insufficient_data",
        suspected_epoch & first_time >= first_suspected_dropout_time ~ "exclude_after_dropout",
        suspected_epoch ~ "suspected_chip_loss",
        raw_observed_fraction < 0.50 ~ "partial_epoch_usable",
        TRUE ~ "usable"
      ),
      movement_after_dropout = if_else(!is.na(first_suspected_dropout_time) & first_time >= first_suspected_dropout_time, movement_epoch_mean, NA_real_),
      proximity_after_dropout = if_else(!is.na(first_suspected_dropout_time) & first_time >= first_suspected_dropout_time, proximity_epoch_mean, NA_real_),
      recommended_action = case_when(
        qc_epoch_class == "exclude_after_dropout" ~ "exclude post-dropout epoch; keep earlier valid epochs",
        qc_epoch_class == "suspected_chip_loss" ~ "manual review; retain only if raw trace is biologically plausible",
        qc_epoch_class == "insufficient_data" ~ "exclude epoch from primary analyses",
        qc_epoch_class == "partial_epoch_usable" ~ "use with duration/missingness sensitivity",
        TRUE ~ "usable"
      )
    ) %>%
    ungroup() %>%
    transmute(
      AnimalNum, Group, Sex, Batch, System, CageChange, Phase,
      first_suspected_dropout_time, last_valid_time,
      observed_fraction = raw_observed_fraction,
      movement_after_dropout, proximity_after_dropout,
      recommended_action,
      qc_epoch_class,
      sudden_sustained_loss, near_zero_movement_after_normal,
      near_zero_proximity_after_normal, long_flatline_immobile,
      expected_observed_mismatch
    )
} else {
  tibble(
    AnimalNum = character(), Group = character(), Sex = character(), Batch = character(), System = character(),
    CageChange = character(), Phase = character(), first_suspected_dropout_time = as.POSIXct(character()),
    last_valid_time = as.POSIXct(character()), observed_fraction = numeric(),
    movement_after_dropout = numeric(), proximity_after_dropout = numeric(),
    recommended_action = character(), qc_epoch_class = character()
  )
}

write_table(qc_chip_loss_flags, file.path(output_dir, "tables/qc_chip_loss_flags.csv"))

chip_loss_epoch_classes <- qc_chip_loss_flags %>%
  select(AnimalNum, CageChange, Phase, qc_epoch_class, recommended_action)

if (nrow(chip_loss_epoch_classes) > 0) {
  base <- base %>%
    left_join(chip_loss_epoch_classes, by = c("AnimalNum", "CageChange", "Phase")) %>%
    mutate(qc_epoch_class = coalesce(qc_epoch_class, "usable")) %>%
    filter(!qc_epoch_class %in% c("exclude_after_dropout", "insufficient_data"))
}

if (nrow(qc_chip_loss_flags) > 0) {
  chip_timeline <- qc_chip_loss_flags %>%
    mutate(Epoch = paste(CageChange, Phase, sep = " | "), AnimalNum = factor(AnimalNum))
  p_chip_timeline <- chip_timeline %>%
    ggplot(aes(Epoch, AnimalNum, fill = qc_epoch_class)) +
    geom_tile(colour = "white", linewidth = 0.12) +
    scale_fill_manual(values = c(
      usable = "#5aa469", partial_epoch_usable = "#d9b44a", suspected_chip_loss = "#d95f02",
      exclude_after_dropout = "#b2182b", insufficient_data = "#777777"
    ), drop = FALSE) +
    labs(title = "RFID chip-loss/dropout QC timeline", x = NULL, y = "Animal", fill = "QC class") +
    make_nature_theme(base_size = 6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  save_plot_svg_pdf(p_chip_timeline, file.path(output_dir, "figures/qc/Fig_chip_loss_dropout_timeline"), width = 175, height = 120)

  p_chip_diag <- qc_chip_loss_flags %>%
    mutate(qc_epoch_class = factor(qc_epoch_class, levels = c("usable", "partial_epoch_usable", "suspected_chip_loss", "exclude_after_dropout", "insufficient_data"))) %>%
    ggplot(aes(movement_after_dropout, proximity_after_dropout, colour = qc_epoch_class)) +
    geom_point(alpha = 0.75, size = 1.7) +
    facet_grid(Sex ~ Group, scales = "free") +
    scale_colour_manual(values = c(
      usable = "#5aa469", partial_epoch_usable = "#d9b44a", suspected_chip_loss = "#d95f02",
      exclude_after_dropout = "#b2182b", insufficient_data = "#777777"
    ), drop = FALSE) +
    labs(title = "Movement/proximity diagnostics after suspected dropout", x = "Movement after dropout", y = "Proximity after dropout", colour = "QC class") +
    make_nature_theme(base_size = 6)
  save_plot_svg_pdf(p_chip_diag, file.path(output_dir, "figures/qc/Fig_chip_loss_movement_proximity_diagnostics"), width = 145, height = 95)
}

# ------------------------------------------------
# BUILD CORE ANIMAL-LEVEL FEATURE MATRIX
# ------------------------------------------------

calc_feature_set <- function(dat, context_name) {
  metric_cols <- intersect(c("Movement", "Entropy", "Proximity"), names(dat))
  if (length(metric_cols) == 0 || nrow(dat) == 0) {
    return(tibble(
      AnimalNum = character(),
      Group = factor(levels = group_levels),
      Sex = factor(levels = sex_levels),
      feature = character(),
      FeatureValue = numeric()
    ))
  }

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
    PhaseClass == "Active" | str_detect(str_to_lower(Phase), "active|dark|night")
  ) %>%
  group_by(AnimalNum, Phase) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  mutate(local_bin = row_number()) %>%
  filter(local_bin <= early_window_bins) %>%
  ungroup()

if (nrow(first_active) == 0) {
  first_active <- base %>%
    filter(as.character(CageChange) == first_cage_change) %>%
    group_by(AnimalNum, Phase) %>%
    arrange(TimeIndex, .by_group = TRUE) %>%
    mutate(local_bin = row_number()) %>%
    filter(local_bin <= early_window_bins) %>%
    ungroup()
}

if (nrow(first_active) > 0) {
  first_active_plot_tbl <- first_active %>%
    mutate(
      local_hour = (local_bin - 1) * infer_bin_size_sec(base) / 3600,
      AnimalNum = factor(AnimalNum)
    )
  first_active_summary_tbl <- first_active_plot_tbl %>%
    group_by(Group, Sex, local_hour) %>%
    summarise(
      mean_movement = mean(Movement, na.rm = TRUE),
      movement_ci_low = mean_ci(Movement)["low"],
      movement_ci_high = mean_ci(Movement)["high"],
      mean_entropy = mean(Entropy, na.rm = TRUE),
      entropy_ci_low = mean_ci(Entropy)["low"],
      entropy_ci_high = mean_ci(Entropy)["high"],
      .groups = "drop"
    )

  p_first_active_movement <- ggplot(first_active_plot_tbl, aes(local_hour, Movement, group = AnimalNum, colour = Group)) +
    geom_line(alpha = 0.16, linewidth = 0.18) +
    geom_ribbon(
      data = first_active_summary_tbl,
      aes(x = local_hour, ymin = movement_ci_low, ymax = movement_ci_high, y = mean_movement, fill = Group, group = Group),
      inherit.aes = FALSE,
      alpha = 0.16,
      colour = NA
    ) +
    geom_line(data = first_active_summary_tbl, aes(y = mean_movement, group = Group), linewidth = 0.55) +
    geom_vline(xintercept = 12, linetype = "dashed", linewidth = 0.2, colour = "grey45") +
    facet_grid(Sex ~ Group) +
    scale_colour_manual(values = group_colors, drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    labs(title = "First active-phase movement trajectory", subtitle = paste0("Raw binned movement during first 12 h after ", first_cage_change), x = "Hours from active-phase onset", y = "Movement") +
    make_nature_theme(base_size = 6) +
    theme(legend.position = "none")
  save_plot_svg_pdf(p_first_active_movement, file.path(output_dir, "figures/publication_panels/Fig_first_active_movement_trajectory_by_group_sex"), width = 175, height = 115)

  p_first_active_entropy <- ggplot(first_active_plot_tbl, aes(local_hour, Entropy, group = AnimalNum, colour = Group)) +
    geom_line(alpha = 0.16, linewidth = 0.18) +
    geom_ribbon(
      data = first_active_summary_tbl,
      aes(x = local_hour, ymin = entropy_ci_low, ymax = entropy_ci_high, y = mean_entropy, fill = Group, group = Group),
      inherit.aes = FALSE,
      alpha = 0.16,
      colour = NA
    ) +
    geom_line(data = first_active_summary_tbl, aes(y = mean_entropy, group = Group), linewidth = 0.55) +
    geom_vline(xintercept = 12, linetype = "dashed", linewidth = 0.2, colour = "grey45") +
    facet_grid(Sex ~ Group) +
    scale_colour_manual(values = group_colors, drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    labs(title = "First active-phase entropy trajectory", subtitle = "Raw entropy trajectory; entropy ACF1/RMSSD summaries are used for primary claims", x = "Hours from active-phase onset", y = "Entropy") +
    make_nature_theme(base_size = 6) +
    theme(legend.position = "none")
  save_plot_svg_pdf(p_first_active_entropy, file.path(output_dir, "figures/publication_panels/Fig_first_active_entropy_acf1_or_instability"), width = 175, height = 115)
}

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

load_hmm_system_features <- function(scale_label = primary_bin_level) {
  hmm_dir <- file.path(project_root, "analysis_ready/06_behavioral_dynamics/hmm_states", scale_label, "tables")
  occ <- read_any_table(file.path(hmm_dir, "hmm_state_occupancy.csv"))
  dwell <- read_any_table(file.path(hmm_dir, "hmm_state_dwell_times.csv"))
  trans <- read_any_table(file.path(hmm_dir, "hmm_transition_probabilities.csv"))
  state_summary <- read_any_table(file.path(hmm_dir, "hmm_state_summary.csv"))
  if (is.null(occ) && is.null(dwell) && is.null(trans)) return(tibble())

  inactive_state <- dominant_state_by_profile(state_summary, "inactive")
  explore_state <- dominant_state_by_profile(state_summary, "explore")
  burst_state <- dominant_state_by_profile(state_summary, "burst")
  social_state <- dominant_state_by_profile(state_summary, "social")

  occ_features <- if (!is.null(occ) && nrow(occ) > 0) {
    occ %>%
      standardize_id_columns() %>%
      mutate(State = as.character(State), frac_time = safe_numeric(frac_time)) %>%
      group_by(AnimalNum, Group, Sex) %>%
      summarise(
        state_occupancy_entropy = feature_entropy(frac_time / sum(frac_time, na.rm = TRUE)),
        social_state_fraction = sum(frac_time[State == social_state], na.rm = TRUE),
        burst_state_fraction = sum(frac_time[State == burst_state], na.rm = TRUE),
        inactive_state_fraction = sum(frac_time[State == inactive_state], na.rm = TRUE),
        exploratory_state_fraction = sum(frac_time[State == explore_state], na.rm = TRUE),
        .groups = "drop"
      )
  } else tibble()

  dwell_features <- if (!is.null(dwell) && nrow(dwell) > 0) {
    dwell %>%
      standardize_id_columns() %>%
      mutate(State = as.character(State)) %>%
      group_by(AnimalNum, Group, Sex) %>%
      summarise(
        mean_dwell_time_per_state = mean(safe_numeric(mean_dwell_bins), na.rm = TRUE),
        max_dwell_time_per_state = max(safe_numeric(max_dwell_bins), na.rm = TRUE),
        inactive_mean_dwell_time = mean(safe_numeric(mean_dwell_bins)[State == inactive_state], na.rm = TRUE),
        burst_mean_dwell_time = mean(safe_numeric(mean_dwell_bins)[State == burst_state], na.rm = TRUE),
        .groups = "drop"
      )
  } else tibble()

  trans_features <- if (!is.null(trans) && nrow(trans) > 0) {
    trans %>%
      standardize_id_columns() %>%
      mutate(
        State = as.character(State),
        NextState = as.character(NextState),
        TransitionProbability = safe_numeric(TransitionProbability)
      ) %>%
      group_by(AnimalNum, Group, Sex) %>%
      summarise(
        self_transition_probability = mean(TransitionProbability[State == NextState], na.rm = TRUE),
        transition_entropy = feature_entropy(TransitionProbability / sum(TransitionProbability, na.rm = TRUE)),
        inactive_to_explore_probability = mean(TransitionProbability[State == inactive_state & NextState == explore_state], na.rm = TRUE),
        inactive_to_burst_probability = mean(TransitionProbability[State == inactive_state & NextState == burst_state], na.rm = TRUE),
        number_of_state_transitions = sum(safe_numeric(Transitions), na.rm = TRUE),
        state_switch_rate = sum(safe_numeric(Transitions)[State != NextState], na.rm = TRUE) / sum(safe_numeric(Transitions), na.rm = TRUE),
        .groups = "drop"
      )
  } else tibble()

  hmm_parts <- list(occ_features, dwell_features, trans_features) %>% keep(~ nrow(.x) > 0)
  if (length(hmm_parts) == 0) return(tibble())

  hmm_wide <- reduce(hmm_parts, full_join, by = c("AnimalNum", "Group", "Sex"))
  for (nm in c(
    "social_state_fraction", "exploratory_state_fraction", "inactive_state_fraction",
    "self_transition_probability", "mean_dwell_time_per_state", "max_dwell_time_per_state",
    "transition_entropy", "state_switch_rate", "inactive_to_explore_probability"
  )) {
    if (!nm %in% names(hmm_wide)) hmm_wide[[nm]] <- NA_real_
  }

  hmm_wide %>%
    mutate(
      resilient_state_bias_index = safe_scale(social_state_fraction) + safe_scale(exploratory_state_fraction) - safe_scale(inactive_state_fraction),
      behavioral_rigidity_index = rowMeans(cbind(safe_scale(self_transition_probability), safe_scale(mean_dwell_time_per_state), safe_scale(max_dwell_time_per_state)), na.rm = TRUE),
      exploratory_flexibility_index = rowMeans(cbind(safe_scale(transition_entropy), safe_scale(state_switch_rate), safe_scale(inactive_to_explore_probability)), na.rm = TRUE)
    ) %>%
    pivot_longer(-c(AnimalNum, Group, Sex), names_to = "MetricName", values_to = "FeatureValue") %>%
    filter(is.finite(FeatureValue)) %>%
    mutate(feature = make_feature_name("hmm", "latent_state", scale_label, MetricName, "animal_level", "all")) %>%
    select(AnimalNum, Group, Sex, feature, FeatureValue)
}

load_gamm_shape_features <- function(scale_label = primary_bin_level) {
  path <- first_existing_path(c(
    file.path(project_root, "analysis_ready/06_behavioral_dynamics/gamm_trajectory_features", scale_label, "tables/gamm_trajectory_features.csv"),
    file.path(project_root, "analysis_ready/06_behavioral_dynamics/gamm_features", scale_label, "tables/combined_gamm_features.csv")
  ))
  dat <- read_any_table(path)
  if (is.null(dat) || nrow(dat) == 0) return(tibble())
  dat <- standardize_id_columns(dat) %>%
    mutate(Metric = if ("Metric" %in% names(.)) safe_name(Metric) else "behavior")

  for (nm in c("auc", "mean_pred", "peak", "trough", "dynamic_range", "time_to_peak", "rmssd_pred", "acf1_pred", "n_bins")) {
    if (!nm %in% names(dat)) dat[[nm]] <- NA_real_
  }

  dat <- dat %>%
    rowwise() %>%
    mutate(
      rebound_slope = ifelse(is.finite(dynamic_range) & is.finite(time_to_peak) & is.finite(n_bins) & n_bins > time_to_peak,
                             -dynamic_range / pmax(n_bins - time_to_peak, 1), NA_real_),
      adaptation_half_life = ifelse(is.finite(time_to_peak), time_to_peak / 2, NA_real_),
      early_peak_amplitude = ifelse(is.finite(peak) & is.finite(mean_pred), peak - mean_pred, NA_real_),
      late_phase_decay = ifelse(is.finite(peak) & is.finite(trough) & is.finite(n_bins), (peak - trough) / pmax(n_bins, 1), NA_real_),
      trajectory_curvature = ifelse(is.finite(dynamic_range) & is.finite(rmssd_pred), dynamic_range * rmssd_pred, NA_real_),
      trajectory_asymmetry = ifelse(is.finite(peak) & is.finite(trough) & is.finite(mean_pred), (peak - mean_pred) - (mean_pred - trough), NA_real_),
      cumulative_instability_auc = ifelse(is.finite(auc) & is.finite(rmssd_pred), abs(auc) * rmssd_pred, NA_real_),
      peak_to_baseline_recovery_time = ifelse(is.finite(n_bins) & is.finite(time_to_peak), pmax(n_bins - time_to_peak, 0), NA_real_)
    ) %>%
    ungroup()

  metric_cols <- intersect(c(
    "auc", "mean_pred", "peak", "trough", "dynamic_range", "time_to_peak", "rmssd_pred", "acf1_pred",
    "rebound_slope", "adaptation_half_life", "early_peak_amplitude", "late_phase_decay",
    "trajectory_curvature", "trajectory_asymmetry", "cumulative_instability_auc", "peak_to_baseline_recovery_time"
  ), names(dat))

  dat %>%
    select(AnimalNum, Group, Sex, Metric, any_of(c("Phase", "CageChange")), all_of(metric_cols)) %>%
    mutate(across(any_of(c("Phase", "CageChange")), ~ safe_name(as.character(.x)))) %>%
    pivot_longer(all_of(metric_cols), names_to = "MetricName", values_to = "FeatureValue") %>%
    unite("ContextTag", any_of(c("Metric", "Phase", "CageChange")), sep = "_", remove = FALSE, na.rm = TRUE) %>%
    filter(is.finite(FeatureValue)) %>%
    mutate(feature = make_feature_name("gamm", "trajectory_geometry", scale_label, MetricName, "animal_level", ContextTag)) %>%
    select(AnimalNum, Group, Sex, feature, FeatureValue)
}

load_nextgen_selective_features <- function(scale_label = primary_bin_level) {
  ng_dir <- file.path(project_root, "analysis_ready/14_nextgen_behavioral_phenotyping", scale_label, "tables")
  candidate_tables <- tibble(
    source_label = c("nextgen_complexity", "nextgen_early_warning", "nextgen_energy_landscape", "nextgen_coupling", "nextgen_integrated"),
    domain_label = c("nonlinear_dynamics", "nonlinear_dynamics", "nonlinear_dynamics", "social_topology", "systems_integration"),
    path = file.path(ng_dir, c(
      "multiscale_complexity_features.csv",
      "early_warning_summary_by_animal.csv",
      "behavioral_energy_landscape_summary_by_animal.csv",
      "animal_social_coupling_summary.csv",
      "nextgen_behavioral_phenotype_matrix.csv"
    ))
  )
  keep_regex <- "mse_auc|perm_entropy|CriticalSlowingIndex|FlickeringIndex|InstabilityRiseIndex|attractor|depth|energy|occupancy|coupling|ComplexityIndex|EnergyLandscapeIndex|NextGenPhenotypeIndex"
  pmap_dfr(candidate_tables, function(source_label, domain_label, path) {
    dat <- read_any_table(path)
    if (is.null(dat) || nrow(dat) == 0) return(tibble())
    dat <- standardize_id_columns(dat)
    metric_cols <- dat %>% select(where(is.numeric)) %>% names()
    metric_cols <- setdiff(metric_cols, c("AnimalNum", endpoint_cols))
    metric_cols <- metric_cols[str_detect(metric_cols, keep_regex)]
    if (length(metric_cols) == 0) return(tibble())
    dat %>%
      select(AnimalNum, Group, Sex, all_of(metric_cols)) %>%
      pivot_longer(all_of(metric_cols), names_to = "MetricName", values_to = "FeatureValue") %>%
      filter(is.finite(FeatureValue)) %>%
      mutate(feature = make_feature_name(source_label, domain_label, scale_label, MetricName, "animal_level", "selected")) %>%
      select(AnimalNum, Group, Sex, feature, FeatureValue)
  })
}

load_graph_period_features <- function(scale_label = primary_bin_level) {
  path <- file.path(project_root, "analysis_ready/06_behavioral_dynamics/social_networks", scale_label, "tables/dyadic_graph_period_summary.csv")
  dat <- read_any_table(path)
  if (is.null(dat) || nrow(dat) == 0) return(tibble())
  group_col <- first_existing_col(dat, c("Group", "Phenotype", "Condition", "Treatment", "StressGroup"), TRUE, "group")
  sex_col <- first_existing_col(dat, c("Sex", "sex"), TRUE, "sex")
  roster <- core_feature_wide %>% distinct(AnimalNum, Group, Sex)
  metric_cols <- intersect(c(
    "fragmentation_index", "mean_density", "mean_modularity", "mean_clustering",
    "mean_components", "mean_largest_component_fraction", "density_rmssd",
    "mean_edges", "mean_strength"
  ), names(dat))
  if (length(metric_cols) == 0) return(tibble())
  dat %>%
    mutate(
      Group = factor(as.character(.data[[group_col]]), levels = group_levels),
      Sex = factor(as.character(.data[[sex_col]]), levels = levels(roster$Sex))
    ) %>%
    select(Group, Sex, any_of(c("Phase", "CageChange", "System")), all_of(metric_cols)) %>%
    left_join(roster, by = c("Group", "Sex")) %>%
    filter(!is.na(AnimalNum)) %>%
    mutate(across(any_of(c("Phase", "CageChange", "System")), ~ safe_name(as.character(.x)))) %>%
    pivot_longer(all_of(metric_cols), names_to = "MetricName", values_to = "FeatureValue") %>%
    unite("ContextTag", any_of(c("Phase", "CageChange", "System")), sep = "_", remove = FALSE, na.rm = TRUE) %>%
    filter(is.finite(FeatureValue)) %>%
    mutate(feature = make_feature_name("dynamic_network", "social_network", scale_label, MetricName, "group_network", ContextTag)) %>%
    select(AnimalNum, Group, Sex, feature, FeatureValue)
}

optional_files <- tibble(
  source_label = c(
    "burstiness", "state_space", "state_space", "early_prediction", "dynamic_network", "dynamic_network",
    "nonlinear_systems", "adaptation_kinetics", "adaptation_kinetics", "sleep_like_inactivity",
    "phase_organization", "phase_organization", "phase_organization", "phase_organization", "phase_organization"
  ),
  domain_label = c(
    "temporal_dynamics", "latent_space", "behavioral_state", "prediction", "social_network", "social_network",
    "nonlinear_dynamics", "recovery_stabilization", "recovery_stabilization", "sleep_like_inactivity",
    "phase_organization", "phase_organization", "phase_organization", "phase_organization", "phase_organization"
  ),
  path = c(
    first_existing_path(c(
      file.path(project_root, "analysis_ready/06_behavioral_dynamics/burstiness", optional_import_bin_levels, "tables/temporal_instability_metrics_per_animal_all_metrics.csv"),
      file.path(project_root, "analysis_ready/06_behavioral_dynamics/burstiness", optional_import_bin_levels, "tables/burstiness_instability_features.csv")
    )),
    first_existing_path(file.path(project_root, "analysis_ready/06_behavioral_dynamics/state_space", optional_import_bin_levels, "tables/state_diversity_metrics.csv")),
    first_existing_path(file.path(project_root, "analysis_ready/06_behavioral_dynamics/state_space", optional_import_bin_levels, "tables/state_switching_metrics.csv")),
    first_existing_path(c(
      file.path(project_root, "analysis_ready/06_behavioral_dynamics/early_prediction", optional_import_bin_levels, "tables/early_behavior_features_wide.csv"),
      file.path(project_root, "analysis_ready/06_behavioral_dynamics/early_prediction", optional_import_bin_levels, "tables/early_behavior_features.csv"),
      file.path(project_root, "analysis_ready/06_behavioral_dynamics/early_prediction_model_ladder", optional_import_bin_levels, "tables/early_behavior_features_wide.csv")
    )),
    first_existing_path(file.path(project_root, "analysis_ready/06_behavioral_dynamics/social_networks", optional_import_bin_levels, "tables/animal_level_social_dynamics.csv")),
    first_existing_path(file.path(project_root, "analysis_ready/06_behavioral_dynamics/social_networks", optional_import_bin_levels, "tables/dyadic_node_summary.csv")),
    first_existing_path(file.path(project_root, "analysis_ready/13_nonlinear_systems_dynamics", optional_import_bin_levels, "derived_data/animal_level_nonlinear_feature_matrix.csv")),
    first_existing_path(file.path(project_root, "analysis_ready/15_behavioral_adaptation_kinetics", optional_import_bin_levels, "tables/adaptation_kinetics_features.csv")),
    first_existing_path(file.path(project_root, "analysis_ready/15_behavioral_adaptation_kinetics", optional_import_bin_levels, "tables/distance_to_control_trajectories.csv")),
    first_existing_path(file.path(project_root, "analysis_ready/16_sleep_like_inactivity_metrics", optional_import_bin_levels, "tables/sleep_like_inactivity_features.csv")),
    first_existing_path(file.path(project_root, "analysis_ready/17_ethological_phase_organization", optional_import_bin_levels, "tables/phase_contrast_features.csv")),
    first_existing_path(file.path(project_root, "analysis_ready/17_ethological_phase_organization", optional_import_bin_levels, "tables/phase_timing_features.csv")),
    first_existing_path(file.path(project_root, "analysis_ready/17_ethological_phase_organization", optional_import_bin_levels, "tables/phase_fragmentation_features.csv")),
    first_existing_path(file.path(project_root, "analysis_ready/17_ethological_phase_organization", optional_import_bin_levels, "tables/phase_recovery_kinetics.csv")),
    first_existing_path(file.path(project_root, "analysis_ready/17_ethological_phase_organization", optional_import_bin_levels, "tables/phase_predictability_features.csv"))
  )
)

optional_long <- pmap_dfr(optional_files, function(source_label, domain_label, path) {
  load_optional_animal_table(path, source_label, domain_label)
}) %>%
  bind_rows(
    map_dfr(optional_import_bin_levels, load_hmm_system_features),
    map_dfr(optional_import_bin_levels, load_gamm_shape_features),
    map_dfr(optional_import_bin_levels, load_graph_period_features),
    map_dfr(optional_import_bin_levels, load_nextgen_selective_features)
  )

optional_wide <- if (nrow(optional_long) > 0) {
  optional_long %>%
    group_by(AnimalNum, Group, Sex, feature) %>%
    summarise(FeatureValue = mean(FeatureValue, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = feature, values_from = FeatureValue)
} else {
  tibble(AnimalNum = character(), Group = factor(levels = group_levels), Sex = factor(levels = sex_levels))
}

systems_features <- full_join(core_feature_wide, optional_wide, by = c("AnimalNum", "Group", "Sex")) %>%
  collapse_to_animal_rows() %>%
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
  endpoint_dat <- endpoint_dat %>%
    select(AnimalNum, all_of(endpoint_keep)) %>%
    collapse_to_animal_rows()
  systems_features <- left_join(systems_features, endpoint_dat, by = "AnimalNum", relationship = "one-to-one")
}

proteomics_dat <- read_any_table(proteomics_module_file)
if (!is.null(proteomics_dat)) {
  proteomics_dat <- standardize_id_columns(proteomics_dat)
  proteomics_numeric <- proteomics_dat %>% select(where(is.numeric)) %>% names()
  proteomics_numeric <- setdiff(proteomics_numeric, "AnimalNum")
  proteomics_dat <- proteomics_dat %>%
    select(AnimalNum, all_of(proteomics_numeric)) %>%
    rename_with(~ paste0("proteomics_module__", safe_name(.x)), all_of(proteomics_numeric)) %>%
    collapse_to_animal_rows()
  systems_features <- left_join(systems_features, proteomics_dat, by = "AnimalNum", relationship = "one-to-one")
}

animal_context_metadata <- bind_rows(
  base %>%
    distinct(AnimalNum, Group, Sex) %>%
    mutate(Batch = NA_character_, System = NA_character_),
  qc_chip_loss_flags %>%
    distinct(AnimalNum, Group, Sex, Batch, System)
) %>%
  mutate(across(any_of(c("Group", "Sex", "Batch", "System")), as.character)) %>%
  group_by(AnimalNum) %>%
  summarise(
    Batch = {
      vals <- na.omit(Batch)
      if (length(vals) == 0) NA_character_ else names(sort(table(vals), decreasing = TRUE))[1]
    },
    System = {
      vals <- na.omit(System)
      if (length(vals) == 0) NA_character_ else names(sort(table(vals), decreasing = TRUE))[1]
    },
    .groups = "drop"
  )

chip_loss_animal_summary <- qc_chip_loss_flags %>%
  group_by(AnimalNum) %>%
  summarise(
    qc__chip_loss__animal__any_suspected_chip_loss__indicator__all = as.numeric(any(qc_epoch_class %in% c("suspected_chip_loss", "exclude_after_dropout"), na.rm = TRUE)),
    qc__chip_loss__animal__any_excluded_post_dropout__indicator__all = as.numeric(any(qc_epoch_class == "exclude_after_dropout", na.rm = TRUE)),
    qc__chip_loss__animal__min_observed_fraction__value__all = min(observed_fraction, na.rm = TRUE),
    qc__chip_loss__animal__n_affected_epochs__count__all = sum(qc_epoch_class %in% c("partial_epoch_usable", "suspected_chip_loss", "exclude_after_dropout", "insufficient_data"), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.infinite(.x) | is.nan(.x), NA_real_, .x)))

if (nrow(animal_context_metadata) > 0) {
  systems_features <- systems_features %>%
    left_join(animal_context_metadata, by = "AnimalNum")
}
if (nrow(chip_loss_animal_summary) > 0) {
  systems_features <- systems_features %>%
    left_join(chip_loss_animal_summary, by = "AnimalNum")
}

chip_loss_feature_audit <- tibble(
  qc_source = "raw/preprocessed animal-position files",
  n_animals_with_any_suspected_chip_loss = if (nrow(chip_loss_animal_summary) > 0) sum(chip_loss_animal_summary$qc__chip_loss__animal__any_suspected_chip_loss__indicator__all %in% 1, na.rm = TRUE) else 0L,
  n_animals_with_post_dropout_exclusion = if (nrow(chip_loss_animal_summary) > 0) sum(chip_loss_animal_summary$qc__chip_loss__animal__any_excluded_post_dropout__indicator__all %in% 1, na.rm = TRUE) else 0L,
  n_epochs_flagged = if (nrow(qc_chip_loss_flags) > 0) sum(!qc_chip_loss_flags$qc_epoch_class %in% "usable", na.rm = TRUE) else 0L,
  downstream_use = "post-dropout/insufficient epochs are removed before feature construction; animal-level QC indicators are retained as QC-only features"
)
write_table(chip_loss_feature_audit, file.path(output_dir, "tables/qc/chip_loss_downstream_feature_audit.csv"))

make_named_biological_scores <- function(dat) {
  candidate_features <- dat %>% select(where(is.numeric)) %>% names()
  candidate_features <- setdiff(candidate_features, c("AnimalNum", endpoint_cols))
  candidate_features <- candidate_features[
    !str_detect(candidate_features, "(__|^)n_bins($|__)|observed_bins|expected_bins|duration|count$|transitions$|switches$|auc$|path_length$")
  ]

  score_from_patterns <- function(positive_regex, negative_regex = NULL) {
    pos <- candidate_features[str_detect(candidate_features, regex(positive_regex, ignore_case = TRUE))]
    neg <- if (is.null(negative_regex)) character(0) else candidate_features[str_detect(candidate_features, regex(negative_regex, ignore_case = TRUE))]
    pos <- setdiff(pos, neg)
    parts <- list()
    if (length(pos) > 0) {
      parts$positive <- dat %>%
        select(all_of(pos)) %>%
        mutate(across(everything(), ~ safe_scale(safe_numeric(.x)))) %>%
        as.matrix()
    }
    if (length(neg) > 0) {
      parts$negative <- dat %>%
        select(all_of(neg)) %>%
        mutate(across(everything(), ~ safe_scale(safe_numeric(.x)))) %>%
        as.matrix()
      parts$negative <- -parts$negative
    }
    if (length(parts) == 0) return(rep(NA_real_, nrow(dat)))
    mat <- do.call(cbind, parts)
    out <- rowMeans(mat, na.rm = TRUE)
    out[!is.finite(out)] <- NA_real_
    out
  }

  rigidity <- score_from_patterns(
    "self_transition|mean_dwell|max_dwell|inactive_mean_dwell|behavioral_rigidity|acf1|inertia|persistence|attractor.*depth",
    "transition_entropy|state_switch_rate|switch.*per_hour|exploratory_flexibility"
  )
  flexibility <- score_from_patterns(
    "transition_entropy|state_switch_rate|switch.*per_hour|exploratory|entropy.*rmssd|movement.*rmssd|contact_entropy|partner_entropy|exploratory_flexibility",
    "self_transition|behavioral_rigidity|inactive_mean_dwell"
  )
  withdrawal <- score_from_patterns(
    "active_isolation|passive_isolation|social_withdrawal|top_partner_share|inactive_state_fraction",
    "proximity.*mean|social_engagement|contact_fraction|mean_contact|social_state_fraction|partner_entropy|partner_evenness"
  )
  fragmentation <- score_from_patterns(
    "fragmentation|modularity|components|active_isolation",
    "density|clustering|largest_component|social_engagement|mean_strength"
  )
  recovery <- score_from_patterns(
    "rebound_slope|recovery|late_phase_decay|trajectory_asymmetry|auc_per_hour",
    "time_to_peak|adaptation_half_life|acf1_pred|cumulative_instability_auc|trajectory_curvature|dynamic_range"
  )
  resilience_bias <- score_from_patterns(
    "resilient_state_bias|social_state_fraction|exploratory_state_fraction|social_engagement|partner_entropy|trajectory_recovery",
    "inactive_state_fraction|behavioral_rigidity|active_isolation|passive_isolation|fragmentation"
  )
  nonlinear_adaptability <- score_from_patterns(
    "complexity|mse|permutation_entropy|manifold_occupancy|coupling|entropy.*mean",
    "criticalslowing|flickering|attractor.*depth|energy.*depth"
  )
  sleep_like_fraction <- score_from_patterns(
    "sleep_like|inactivity_fraction|quiescence|inactive_fraction|prolonged_inactivity_fraction",
    "fragmentation|switch|bout_count|episodes"
  )
  sleep_bout_duration <- score_from_patterns(
    "mean_inactivity_bout|min|max_inactivity_bout|median_inactivity_bout|bout_duration",
    "fragmentation|switch|bout_count"
  )
  sleep_fragmentation <- score_from_patterns(
    "inactivity_bout_count|fragmentation|state_switch|transition|prolonged_inactivity_episodes",
    "fraction|duration|consolidation"
  )
  active_phase_quiescence <- score_from_patterns(
    "active.*inactivity|active.*quiescence|active.*inactive|phase.*active",
    "inactive_phase_consolidation"
  )
  inactive_phase_consolidation <- score_from_patterns(
    "inactive.*consolidation|inactive.*bout|inactive.*fraction|inactive.*quiescence|prolonged_inactivity_fraction",
    "fragmentation|switch"
  )

  named <- tibble(
    AnimalNum = dat$AnimalNum,
    !!make_feature_name("systems", "biological_scores", primary_bin_level, "behavioral_rigidity_score", "module_score", "all") := rigidity,
    !!make_feature_name("systems", "biological_scores", primary_bin_level, "exploratory_flexibility_score", "module_score", "all") := flexibility,
    !!make_feature_name("systems", "biological_scores", primary_bin_level, "social_withdrawal_score", "module_score", "all") := withdrawal,
    !!make_feature_name("systems", "biological_scores", primary_bin_level, "social_fragmentation_score", "module_score", "all") := fragmentation,
    !!make_feature_name("systems", "biological_scores", primary_bin_level, "trajectory_recovery_score", "module_score", "all") := recovery,
    !!make_feature_name("systems", "biological_scores", primary_bin_level, "resilient_state_bias_score", "module_score", "all") := resilience_bias,
    !!make_feature_name("systems", "biological_scores", primary_bin_level, "nonlinear_adaptability_score", "module_score", "all") := nonlinear_adaptability,
    !!make_feature_name("systems", "biological_scores", primary_bin_level, "sleep_like_fraction_score", "module_score", "all") := sleep_like_fraction,
    !!make_feature_name("systems", "biological_scores", primary_bin_level, "sleep_bout_duration_score", "module_score", "all") := sleep_bout_duration,
    !!make_feature_name("systems", "biological_scores", primary_bin_level, "sleep_fragmentation_score", "module_score", "all") := sleep_fragmentation,
    !!make_feature_name("systems", "biological_scores", primary_bin_level, "active_phase_quiescence_score", "module_score", "all") := active_phase_quiescence,
    !!make_feature_name("systems", "biological_scores", primary_bin_level, "inactive_phase_consolidation_score", "module_score", "all") := inactive_phase_consolidation
  ) %>%
    mutate(
      !!make_feature_name("systems", "biological_scores", primary_bin_level, "stress_adaptation_index", "module_score", "all") := rowMeans(
        cbind(
          safe_scale(.data[[make_feature_name("systems", "biological_scores", primary_bin_level, "exploratory_flexibility_score", "module_score", "all")]]),
          safe_scale(.data[[make_feature_name("systems", "biological_scores", primary_bin_level, "trajectory_recovery_score", "module_score", "all")]]),
          safe_scale(.data[[make_feature_name("systems", "biological_scores", primary_bin_level, "resilient_state_bias_score", "module_score", "all")]]),
          -safe_scale(.data[[make_feature_name("systems", "biological_scores", primary_bin_level, "behavioral_rigidity_score", "module_score", "all")]]),
          -safe_scale(.data[[make_feature_name("systems", "biological_scores", primary_bin_level, "social_withdrawal_score", "module_score", "all")]]),
          -safe_scale(.data[[make_feature_name("systems", "biological_scores", primary_bin_level, "social_fragmentation_score", "module_score", "all")]])
        ),
        na.rm = TRUE
      ),
      !!make_feature_name("systems", "biological_scores", primary_bin_level, "systems_resilience_score", "module_score", "all") := rowMeans(
        cbind(
          safe_scale(.data[[make_feature_name("systems", "biological_scores", primary_bin_level, "stress_adaptation_index", "module_score", "all")]]),
          safe_scale(.data[[make_feature_name("systems", "biological_scores", primary_bin_level, "nonlinear_adaptability_score", "module_score", "all")]])
        ),
        na.rm = TRUE
      )
    )

  named
}

named_biological_scores <- make_named_biological_scores(systems_features)
systems_features <- left_join(systems_features, named_biological_scores, by = "AnimalNum", relationship = "one-to-one")

feature_cols <- systems_features %>% select(where(is.numeric)) %>% names()
feature_cols <- setdiff(feature_cols, c("AnimalNum", endpoint_cols))

is_raw_duration_sensitive_feature <- function(feature, metric, statistic, context) {
  key <- str_to_lower(paste(feature, metric, statistic, context, sep = " "))
  key[is.na(key)] <- ""
  raw_signal <- str_detect(
    key,
    "transition|switch|count|n_bins|n_transitions|n_switches|path_length|roughness|auc|contact_bins|seconds|duration|dwell_bins|fragmentation_event"
  )
  normalized_signal <- str_detect(
    key,
    "per_hour|per_epoch|rate|fraction|probability|occupancy|mean|median|rmssd|acf1|entropy_rate|normalized|completeness"
  )
  isTRUE(raw_signal) && !isTRUE(normalized_signal)
}

source_script_from_source <- function(source) {
  case_when(
    source == "raw" ~ "03_build_multiscale_behavior_metrics.R",
    str_detect(source, "burstiness") ~ "06_burstiness_temporal_instability.R",
    str_detect(source, "state_space") ~ "07_behavioral_state_space.R",
    str_detect(source, "early_prediction|prediction") ~ "08b_early_prediction_model_ladder.R",
    str_detect(source, "dynamic_network") ~ "09_dynamic_social_networks.R",
    str_detect(source, "hmm") ~ "10_hmm_behavioral_states.R",
    str_detect(source, "gamm") ~ "11_gamm_trajectory_features.R",
    str_detect(source, "behavior_proteomics") ~ "12_behavior_proteomics_integration.R",
    str_detect(source, "nonlinear_systems") ~ "13_nonlinear_systems_dynamics.R",
    str_detect(source, "qc|chip_loss") ~ "00_tracking_qc_rfid_loss.R",
    str_detect(source, "adaptation_kinetics") ~ "15_behavioral_adaptation_kinetics.R",
    str_detect(source, "sleep_like_inactivity") ~ "16_sleep_like_inactivity_metrics.R",
    str_detect(source, "phase_organization") ~ "17_ethological_phase_organization.R",
    str_detect(source, "nextgen|nonlinear") ~ "14_nextgen_behavioral_phenotyping.R",
    str_detect(source, "systems") ~ "12_systems_neuroscience_summary.R",
    TRUE ~ "integrated_optional_output"
  )
}

biological_domain_from_module <- function(module, domain, metric, context) {
  key <- str_to_lower(paste(module, domain, metric, context))
  case_when(
    str_detect(key, "prediction|early") ~ "Early behavioral adaptation",
    str_detect(key, "temporal|rmssd|acf1|instability|persistence|burst") ~ "Temporal instability and organization",
    str_detect(key, "phase|active|inactive|sleep_like|quiescence") ~ "Ethological active/inactive phase organization",
    str_detect(key, "social|network|partner|proximity") ~ "Social reorganization after regrouping",
    str_detect(key, "proteomics|molecular") ~ "Behavioral-proteomic systems alignment",
    str_detect(key, "chip_loss|dropout|qc|batch|system") ~ "Quality control and confound audit",
    str_detect(key, "trajectory|recovery|stabilization|adaptation") ~ "Longitudinal adaptation and recovery",
    str_detect(key, "hmm|latent_state") ~ "Behavioral state organization",
    str_detect(key, "nonlinear|complexity|attractor|manifold") ~ "Exploratory nonlinear systems dynamics",
    TRUE ~ "Core behavioral phenotype"
  )
}

claim_type_from_feature <- function(module, source, feature) {
  key <- str_to_lower(paste(module, source, feature))
  case_when(
    str_detect(key, "prediction|combz|outcome") ~ "predictive",
    str_detect(key, "proteomics|molecular") ~ "associative",
    str_detect(key, "chip_loss|dropout|qc|batch|system") ~ "QC-only",
    str_detect(key, "nextgen|nonlinear|manifold|attractor|energy|recurrence") ~ "exploratory",
    str_detect(key, "hmm|state|network|trajectory|phase|inactivity") ~ "descriptive",
    TRUE ~ "descriptive"
  )
}

evidence_tier_from_feature <- function(module, biological_domain, source, feature) {
  key <- str_to_lower(paste(module, biological_domain, source, feature))
  case_when(
    str_detect(key, "first_active_12h|movement.*mean|movement.*rmssd|entropy.*acf1|prediction") ~ "Tier 1",
    str_detect(key, "hmm|latent|phase|social|trajectory|recovery|inactivity|instability") ~ "Tier 2",
    str_detect(key, "nextgen|nonlinear|manifold|attractor|energy|recurrence|complexity") ~ "Tier 3",
    TRUE ~ "Tier 2"
  )
}

feature_dictionary <- parse_system_feature(feature_cols) %>%
  mutate(
    Module = feature_module_from_parts(Source, Domain, Metric, Statistic, Context, feature),
    RawDurationSensitive = pmap_lgl(
      list(feature, Metric, Statistic, Context),
      is_raw_duration_sensitive_feature
    ),
    DurationUse = if_else(
      RawDurationSensitive,
      "descriptive_audit_only_use_normalized_counterpart_for_embedding_or_prediction",
      "duration_robust_or_normalized"
    ),
    Interpretation = case_when(
      str_detect(feature, "movement") & str_detect(feature, "mean") ~ "Mean locomotor output",
      str_detect(feature, "entropy") & str_detect(feature, "mean") ~ "Mean spatial dispersion / positional entropy",
      str_detect(feature, "proximity") & str_detect(feature, "mean") ~ "Mean social proximity",
      str_detect(feature, "rmssd") ~ "Short-timescale temporal instability",
      str_detect(feature, "acf1") ~ "Temporal inertia / persistence",
      str_detect(feature, "social_withdrawal") ~ "High movement with low proximity; exploratory-social decoupling",
      str_detect(feature, "behavioral_instability") ~ "Composite local volatility across movement, entropy and proximity",
      str_detect(feature, "behavioral_inertia") ~ "Composite autocorrelation/persistence across behavioral channels",
      Module == "Latent-state organization" ~ "Latent behavioral-state persistence, occupancy, dwell time or transition flexibility",
      Module == "Social topology" ~ "Social-network topology, centrality, fragmentation or contact organization",
      Module == "Trajectory geometry" ~ "Shape of smooth behavioral adaptation trajectories",
      Module == "Nonlinear systems dynamics" ~ "Exploratory nonlinear complexity, recurrence, attractor or early-warning dynamics",
      Module == "Predictive systems integration" ~ "Module-level or integrated predictive phenotype score",
      Domain == "biological_scores" ~ "Named biological composite score for systems-level interpretation",
      TRUE ~ "Integrated behavioral feature"
    ),
    Feature = feature,
    SourceScript = source_script_from_source(Source),
    BiologicalDomain = biological_domain_from_module(Module, Domain, Metric, Context),
    MathematicalDefinition = case_when(
      str_detect(Statistic, "mean|animal_level|module_score") ~ "Animal-level mean or module score computed from the source table.",
      str_detect(Statistic, "rmssd") ~ "Root mean square of successive differences across ordered time bins.",
      str_detect(Statistic, "acf1") ~ "Lag-1 autocorrelation across ordered time bins.",
      str_detect(Metric, "ratio") ~ "Ratio-style normalized phase or behavioral contrast.",
      str_detect(Metric, "per_hour|rate") ~ "Count or cumulative quantity divided by observed hours.",
      TRUE ~ "Source-defined animal-level summary; see SourceScript for exact preprocessing."
    ),
    TimeWindow = case_when(
      str_detect(Context, "first_active_12h|early") ~ "First active phase after first cage change, first 12 h",
      str_detect(Context, "phase") ~ "Phase-specific cage-change epoch",
      str_detect(Context, "all") ~ "All available cage-change epochs",
      TRUE ~ Context
    ),
    BinLevel = Scale,
    ClaimType = claim_type_from_feature(Module, Source, feature),
    EvidenceTier = evidence_tier_from_feature(Module, BiologicalDomain, Source, feature),
    AllowedInterpretation = Interpretation,
    ReviewerRisk = case_when(
      RawDurationSensitive ~ "high",
      Module %in% c("Nonlinear systems dynamics", "Trajectory geometry", "Latent-state organization") ~ "medium",
      TRUE ~ "low"
    ),
    DurationSensitive = RawDurationSensitive,
    StableForMainText = EvidenceTier %in% c("Tier 1", "Tier 2") & !RawDurationSensitive & ReviewerRisk != "high"
  ) %>%
  select(
    feature, Feature, Source, SourceScript, Domain, BiologicalDomain, Module,
    Metric, Statistic, Context, Scale, BinLevel, MathematicalDefinition,
    TimeWindow, EvidenceTier, ClaimType, AllowedInterpretation, ReviewerRisk,
    DurationSensitive, StableForMainText, RawDurationSensitive, DurationUse,
    Interpretation
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
  arrange(Module, Source, Domain, Scale, Metric, Statistic)

usable_features <- feature_qc %>%
  filter(n_finite >= 4, missing_fraction <= 0.50, !zero_variance) %>%
  pull(feature)

duration_robust_features <- feature_qc %>%
  filter(n_finite >= 4, missing_fraction <= 0.50, !zero_variance, !RawDurationSensitive) %>%
  pull(feature)

prospective_prediction_features <- feature_qc %>%
  filter(
    feature %in% duration_robust_features,
    (
      Source == "raw" & Context == prospective_prediction_context
    ) |
      (
        Source == "systems" & Context == prospective_prediction_context
      )
  ) %>%
  pull(feature)

module_feature_sets <- feature_qc %>%
  filter(feature %in% usable_features) %>%
  group_by(Module) %>%
  summarise(
    n_usable_features = n(),
    n_duration_robust_features = sum(!RawDurationSensitive, na.rm = TRUE),
    n_raw_duration_sensitive_features = sum(RawDurationSensitive, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(factor(Module, levels = c(
    "Magnitude", "Temporal instability", "Latent-state organization", "Social topology",
    "Trajectory geometry", "Nonlinear systems dynamics", "Predictive systems integration",
    "Other interpretable feature"
  )))

feature_sets <- tibble(
  FeatureSet = c("descriptive_full_experiment", "duration_robust_embedding_prediction", "prospective_prediction"),
  n_features = c(length(usable_features), length(duration_robust_features), length(prospective_prediction_features)),
  IntendedUse = c(
    "Descriptive phenotype architecture after RES/SUS labels are known",
    "PCA/UMAP and module-level predictive models using normalized or duration-robust features",
    "Prediction/association with post-paradigm endpoints using pre-endpoint early behavior only"
  ),
  IncludesCombZDerivedGroups = c(res_sus_derived_from_outcome, res_sus_derived_from_outcome, FALSE),
  Notes = c(
    "RES/SUS group contrasts are descriptive because groups were derived from CombZ.",
    "Excludes raw counts, cumulative AUC/path-length features, and other unnormalized duration-sensitive summaries when a normalized counterpart is available.",
    paste0("Uses raw first active 12 h after ", first_cage_change, " plus early composites; excludes GAMM/HMM/full-experiment and mixed-phase early-prediction module features.")
  )
)

pick_feature <- function(pattern, candidates = feature_qc$feature, require_robust = TRUE) {
  pool <- candidates[str_detect(candidates, regex(pattern, ignore_case = TRUE))]
  if (require_robust) pool <- intersect(pool, duration_robust_features)
  if (length(pool) == 0) return(NA_character_)
  pool[1]
}

primary_feature_registry <- tibble(
  PrimaryFeatureLabel = c(
    "First active 12 h Movement mean",
    "First active 12 h Movement RMSSD",
    "First active 12 h Entropy ACF1",
    "First active 12 h Proximity mean / social withdrawal",
    "Sleep-like inactivity fraction",
    "Sleep-like bout duration",
    "Sleep-like fragmentation",
    "HMM rigidity/flexibility",
    "Systems adaptation/resilience composite"
  ),
  feature = c(
    early_movement_mean,
    early_movement_rmssd,
    make_feature_name("raw", "behavior", primary_bin_level, "Entropy", "acf1", prospective_prediction_context),
    coalesce(
      early_proximity_mean,
      pick_feature("social_withdrawal_score|social_withdrawal")
    ),
    pick_feature("sleep_like_fraction_score|inactivity_fraction|quiescence"),
    pick_feature("sleep_bout_duration_score|mean_inactivity_bout|bout_duration"),
    pick_feature("sleep_fragmentation_score|fragmentation|inactivity_bout_count"),
    pick_feature("behavioral_rigidity_score|exploratory_flexibility_score|rigidity_index|flexibility_index|state_switch_rate"),
    pick_feature("stress_adaptation_index|systems_resilience_score")
  ),
  ClaimTier = c("primary", "primary", "primary", "primary", "primary", "primary", "primary", "secondary", "secondary"),
  BiologicalClaim = c(
    "Early locomotor output after first social instability exposure",
    "Early locomotor volatility/adaptation",
    "Early spatial entropy persistence/inertia",
    "Early social engagement or withdrawal",
    "Quiescence-like inactivity amount",
    "Quiescence-like bout architecture",
    "Quiescence-like fragmentation",
    "Latent-state rigidity/flexibility",
    "Integrated adaptation/resilience composite"
  )
) %>%
  filter(!is.na(feature), feature %in% feature_qc$feature) %>%
  distinct(feature, .keep_all = TRUE) %>%
  left_join(feature_dictionary %>% select(feature, SourceScript, Module, ReviewerRisk, DurationSensitive, ClaimType), by = "feature") %>%
  mutate(
    FeatureUseClass = ClaimTier,
    MainFigureEligible = ClaimTier %in% c("primary", "secondary") & !DurationSensitive %in% TRUE & ReviewerRisk != "high",
    LeakageClass = if_else(str_detect(feature, prospective_prediction_context), "safe_for_prediction", "descriptive_group_contrast")
  )

feature_dictionary <- feature_dictionary %>%
  left_join(primary_feature_registry %>% select(feature, FeatureUseClass, PrimaryFeatureLabel, BiologicalClaim, LeakageClass), by = "feature") %>%
  mutate(
    FeatureUseClass = coalesce(FeatureUseClass, case_when(
      Source == "qc" | BiologicalDomain == "Quality control and confound audit" ~ "QC-only",
      EvidenceTier == "Tier 3" ~ "exploratory",
      TRUE ~ "secondary"
    )),
    LeakageClass = coalesce(LeakageClass, case_when(
      str_detect(Context, prospective_prediction_context) & Source %in% c("raw", "systems") ~ "safe_for_prediction",
      ClaimType == "predictive" & !str_detect(Context, prospective_prediction_context) ~ "leaky_not_for_prediction",
      ClaimType == "associative" ~ "prospective_endpoint_association",
      TRUE ~ "descriptive_group_contrast"
    )),
    PrimaryFeatureLabel = coalesce(PrimaryFeatureLabel, Feature),
    BiologicalClaim = coalesce(BiologicalClaim, AllowedInterpretation)
  )

feature_qc <- feature_qc %>%
  select(-any_of(c("FeatureUseClass", "PrimaryFeatureLabel", "BiologicalClaim", "LeakageClass"))) %>%
  left_join(feature_dictionary %>% select(feature, FeatureUseClass, PrimaryFeatureLabel, BiologicalClaim, LeakageClass), by = "feature")

primary_claim_features <- primary_feature_registry %>% filter(ClaimTier == "primary") %>% pull(feature)
secondary_claim_features <- primary_feature_registry %>% filter(ClaimTier == "secondary") %>% pull(feature)
main_dashboard_features <- primary_feature_registry %>% filter(MainFigureEligible %in% TRUE) %>% pull(feature)

path_status <- function(paths) {
  paths <- paths[!is.na(paths)]
  if (length(paths) == 0) return(FALSE)
  any(file.exists(paths))
}

integration_audit_registry <- tibble(
  feature_module_name = c(
    "Core multiscale behavior",
    "Temporal burstiness and instability",
    "Behavioral state-space",
    "Early prediction features and model ladder",
    "Dynamic social networks",
    "HMM latent behavioral states",
    "GAMM trajectory geometry",
    "Nonlinear systems dynamics",
    "Next-generation nonlinear phenotyping",
    "Behavioral adaptation kinetics",
    "Sleep-like inactivity and quiescence",
    "Ethological active/inactive phase organization",
    "Behavior-proteomics integration",
    "Chip-loss and low-observation QC"
  ),
  SourceScript = c(
    "03_build_multiscale_behavior_metrics.R",
    "06_burstiness_temporal_instability.R",
    "07_behavioral_state_space.R",
    "08b_early_prediction_model_ladder.R",
    "09_dynamic_social_networks.R",
    "10_hmm_behavioral_states.R",
    "11_gamm_trajectory_features.R",
    "13_nonlinear_systems_dynamics.R",
    "14_nextgen_behavioral_phenotyping.R",
    "15_behavioral_adaptation_kinetics.R",
    "16_sleep_like_inactivity_metrics.R",
    "17_ethological_phase_organization.R",
    "12_behavior_proteomics_integration.R",
    "00_tracking_qc_rfid_loss.R"
  ),
  expected_outputs = I(list(
    file.path(project_root, "analysis_ready/03_derived_metrics", primary_bin_level, "all_behavior_metrics.csv"),
    file.path(project_root, "analysis_ready/06_behavioral_dynamics/burstiness", primary_bin_level, "tables/temporal_instability_metrics_per_animal_all_metrics.csv"),
    c(
      file.path(project_root, "analysis_ready/06_behavioral_dynamics/state_space", c("1min_based", optional_import_bin_levels), "tables/state_diversity_metrics.csv"),
      file.path(project_root, "analysis_ready/06_behavioral_dynamics/state_space", c("1min_based", optional_import_bin_levels), "tables/state_switching_metrics.csv")
    ),
    c(
      file.path(project_root, "analysis_ready/06_behavioral_dynamics/early_prediction", optional_import_bin_levels, "tables/early_behavior_features.csv"),
      file.path(project_root, "analysis_ready/06_behavioral_dynamics/early_prediction_model_ladder", optional_import_bin_levels, "tables/early_behavior_features_wide.csv")
    ),
    file.path(project_root, "analysis_ready/06_behavioral_dynamics/social_networks", primary_bin_level, "tables/animal_level_social_dynamics.csv"),
    file.path(project_root, "analysis_ready/06_behavioral_dynamics/hmm_states", optional_import_bin_levels, "tables/hmm_state_occupancy.csv"),
    c(
      file.path(project_root, "analysis_ready/06_behavioral_dynamics/gamm_trajectory_features", optional_import_bin_levels, "tables/combined_gamm_features.csv"),
      file.path(project_root, "analysis_ready/06_behavioral_dynamics/gamm_features", optional_import_bin_levels, "tables/combined_gamm_features.csv")
    ),
    file.path(project_root, "analysis_ready/13_nonlinear_systems_dynamics", primary_bin_level, "derived_data/animal_level_nonlinear_feature_matrix.csv"),
    file.path(project_root, "analysis_ready/14_nextgen_behavioral_phenotyping", primary_bin_level, "tables/nextgen_behavioral_phenotype_matrix.csv"),
    file.path(project_root, "analysis_ready/15_behavioral_adaptation_kinetics", primary_bin_level, "tables/adaptation_kinetics_features.csv"),
    file.path(project_root, "analysis_ready/16_sleep_like_inactivity_metrics", primary_bin_level, "tables/sleep_like_inactivity_features.csv"),
    file.path(project_root, "analysis_ready/17_ethological_phase_organization", primary_bin_level, "tables/phase_contrast_features.csv"),
    c(
      proteomics_module_file,
      file.path(project_root, "analysis_ready/12_behavior_proteomics_integration", "tables/behavior_proteomics_merged.csv")
    ),
    c(
      file.path(project_root, "analysis_ready/00_tracking_qc_rfid_loss", "tables/tracking_qc_by_animal.csv"),
      file.path(dirname(project_root), "raw_tracking_qc_rfid_loss", "tables/raw_tracking_qc_by_animal.csv")
    )
  )),
  expected_status = c(
    "implemented",
    "implemented",
    "implemented",
    "implemented",
    "implemented",
    "implemented",
    "implemented",
    "implemented",
    "implemented",
    "implemented",
    "implemented",
    "implemented",
    if_else(is.null(proteomics_module_file), "optional_not_requested", "implemented_if_table_valid"),
    "qc_available_if_upstream_run"
  ),
  reviewer_risk = c(
    "low", "medium", "medium", "medium", "medium", "medium", "medium",
    "high", "high", "medium", "medium", "medium", "medium", "high"
  ),
  suggested_figure_placement = c(
    "Main dashboard panel A/state-space plus module summary",
    "Main module summary; detailed temporal plots in supplement",
    "Supplement unless central state-space claim is emphasized",
    "Main dashboard prediction panel and robustness table",
    "Supplement; summarize as social topology module in main table",
    "Supplement; semantic state labels and transition differences only if stable",
    "Supplement; use trajectory/recovery module summary in main table",
    "Supplement/exploratory only",
    "Supplement/exploratory only",
    "Supplement; summarize recovery/stabilization module in main table",
    "Supplement; avoid EEG sleep claims",
    "Supplement; avoid circadian/sleep claims",
    "Supplement/associative only, and only when proteomics module table is provided",
    "QC/robustness table, not a biological figure"
  )
)

main_dashboard_source_scripts <- c(
  "03_build_multiscale_behavior_metrics.R",
  "06_burstiness_temporal_instability.R",
  "08b_early_prediction_model_ladder.R",
  "12_systems_neuroscience_summary.R"
)

feature_import_summary <- feature_dictionary %>%
  count(SourceScript, name = "n_imported_features") %>%
  mutate(imported_into_dashboard = n_imported_features > 0)

systems_computation_integration_audit <- integration_audit_registry %>%
  mutate(
    already_available = map_lgl(expected_outputs, path_status),
    expected_output_path = map_chr(expected_outputs, ~ paste(.x[!is.na(.x)], collapse = "; "))
  ) %>%
  left_join(feature_import_summary, by = "SourceScript") %>%
  mutate(
    n_imported_features = replace_na(n_imported_features, 0L),
    imported_into_dashboard = replace_na(imported_into_dashboard, FALSE),
    used_in_main_figure = SourceScript %in% main_dashboard_source_scripts & imported_into_dashboard,
    integration_status = case_when(
      expected_status == "optional_not_requested" ~ "optional_not_requested",
      !already_available ~ "missing_upstream_output",
      already_available & imported_into_dashboard ~ "available_and_imported",
      already_available & !imported_into_dashboard ~ "available_not_imported_or_qc_only",
      TRUE ~ "needs_review"
    ),
    reviewer_action = case_when(
      feature_module_name == "Behavior-proteomics integration" & is.null(proteomics_module_file) ~
        "Do not interpret behavior-proteomics; set proteomics_module_file to an animal-level module table to enable.",
      feature_module_name == "Chip-loss and low-observation QC" ~
        "Use as exclusion/sensitivity evidence; do not treat QC flags as biological phenotypes.",
      reviewer_risk == "high" ~
        "Keep out of main dashboard unless supported by robustness and independent validation.",
      TRUE ~ "Use source-specific statistics and robustness tables before main-text claims."
    )
  ) %>%
  select(
    feature_module_name, SourceScript, already_available, imported_into_dashboard,
    n_imported_features, used_in_main_figure, integration_status, reviewer_risk,
    suggested_figure_placement, reviewer_action, expected_output_path
  )

systems_robustness_audit <- tibble(
  robustness_check = c(
    "Permutation baseline",
    "Bootstrap/CI for prediction or correlations",
    "Cross-validation",
    "Duration sensitivity",
    "Chip-loss / detached-chip QC",
    "Low-observation epoch QC",
    "Proteomics module availability"
  ),
  dashboard_evidence = c(
    file.path(output_dir, "tables/systems_prediction_ladder_performance.csv"),
    file.path(output_dir, "tables/systems_prediction_ladder_performance.csv"),
    file.path(output_dir, "tables/systems_prediction_ladder_loo_predictions.csv"),
    file.path(output_dir, "tables/systems_prediction_ladder_performance_duration_sensitivity.csv"),
    file.path(project_root, "analysis_ready/00_tracking_qc_rfid_loss", "tables/tracking_qc_by_animal.csv"),
    file.path(output_dir, "tables/systems_duration_negative_control_by_animal.csv"),
    if (is.null(proteomics_module_file)) NA_character_ else proteomics_module_file
  ),
  upstream_evidence = c(
    first_existing_path(c(
      file.path(project_root, "analysis_ready/06_behavioral_dynamics/early_prediction_model_ladder", primary_bin_level, "tables/model_ladder_performance.csv"),
      file.path(project_root, "analysis_ready/06_behavioral_dynamics/early_prediction", primary_bin_level, "tables/model_ladder_performance.csv")
    )),
    first_existing_path(c(
      file.path(project_root, "analysis_ready/06_behavioral_dynamics/early_prediction_model_ladder", primary_bin_level, "tables/model_ladder_performance.csv"),
      file.path(project_root, "analysis_ready/06_behavioral_dynamics/early_prediction", primary_bin_level, "tables/model_ladder_performance.csv")
    )),
    first_existing_path(c(
      file.path(project_root, "analysis_ready/06_behavioral_dynamics/early_prediction_model_ladder", primary_bin_level, "tables/model_ladder_repeated_grouped_kfold_performance.csv"),
      file.path(project_root, "analysis_ready/06_behavioral_dynamics/early_prediction", primary_bin_level, "tables/model_ladder_repeated_grouped_kfold_performance.csv")
    )),
    first_existing_path(c(
      file.path(project_root, "analysis_ready/06_behavioral_dynamics/early_prediction_model_ladder", primary_bin_level, "tables/model_ladder_performance_duration_sensitivity.csv"),
      file.path(project_root, "analysis_ready/06_behavioral_dynamics/early_prediction", primary_bin_level, "tables/model_ladder_performance_duration_sensitivity.csv")
    )),
    file.path(project_root, "analysis_ready/00_tracking_qc_rfid_loss", "tables/tracking_qc_by_animal.csv"),
    file.path(output_dir, "tables/duration_sensitivity"),
    file.path(project_root, "analysis_ready/12_behavior_proteomics_integration", "tables/behavior_proteomics_merged.csv")
  ),
  interpretation = c(
    "Permutation p-values benchmark whether LOOCV performance exceeds label-shuffled baselines.",
    "Prediction tables contain bootstrap/Fisher-z intervals when available; group contrasts carry effect-size CIs.",
    "Main dashboard uses LOOCV; upstream 08b repeated grouped k-fold is an external robustness companion if present.",
    "Main prediction ladder is rerun after excluding animals with short-duration epochs when enough animals remain.",
    "RFID chip-loss QC is a data-validity guard, not a behavioral endpoint.",
    "Duration and completeness screens reduce false positives from low-observation epochs.",
    "Behavior-proteomics claims remain disabled unless an animal-level proteomics module table is configured."
  )
) %>%
  mutate(
    dashboard_available = !is.na(dashboard_evidence) & file.exists(dashboard_evidence),
    upstream_available = !is.na(upstream_evidence) & file.exists(upstream_evidence),
    status = case_when(
      dashboard_available ~ "dashboard_available",
      upstream_available ~ "upstream_available_not_dashboard_specific",
      robustness_check == "Proteomics module availability" & is.null(proteomics_module_file) ~ "optional_not_requested",
      TRUE ~ "missing_or_not_run"
    )
  )

write_table(systems_features, file.path(output_dir, "tables/systems_animal_feature_matrix.csv"))
write_table(feature_dictionary, file.path(output_dir, "tables/systems_feature_dictionary.csv"))
write_table(feature_qc, file.path(output_dir, "tables/systems_feature_qc.csv"))
write_table(feature_sets, file.path(output_dir, "tables/systems_feature_sets.csv"))
write_table(module_feature_sets, file.path(output_dir, "tables/systems_module_feature_inventory.csv"))
write_table(primary_feature_registry, file.path(output_dir, "tables/systems_primary_feature_set.csv"))
write_table(systems_computation_integration_audit, file.path(output_dir, "tables/systems_computation_integration_audit.csv"))
write_table(systems_robustness_audit, file.path(output_dir, "tables/systems_robustness_audit.csv"))

named_biological_scores_long <- systems_features %>%
  select(AnimalNum, Group, Sex, matches("^systems__biological_scores__")) %>%
  pivot_longer(matches("^systems__biological_scores__"), names_to = "feature", values_to = "Score") %>%
  left_join(feature_dictionary, by = "feature") %>%
  mutate(
    ScoreLabel = Metric %>%
      str_replace_all("_", " ") %>%
      str_to_sentence()
  )

write_table(named_biological_scores_long, file.path(output_dir, "tables/systems_named_biological_scores.csv"))

# ------------------------------------------------
# GROUP SUMMARY AND CONTRASTS
# ------------------------------------------------

systems_long <- systems_features %>%
  select(AnimalNum, Group, Sex, all_of(usable_features)) %>%
  pivot_longer(all_of(usable_features), names_to = "feature", values_to = "Value") %>%
  filter(is.finite(Value)) %>%
  left_join(feature_dictionary, by = "feature")

group_summary <- systems_long %>%
  group_by(Sex, Group, feature, Module, Source, Domain, Scale, Metric, Statistic, Context) %>%
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
  group_by(Sex, feature, Module, Source, Domain, Scale, Metric, Statistic, Context) %>%
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
# PRIMARY CLAIMS, CONFOUND AUDIT AND GROUP x SEX MODELS
# ------------------------------------------------

systems_primary_claims_table <- primary_feature_registry %>%
  left_join(
    group_contrasts %>%
      filter(feature %in% primary_feature_registry$feature) %>%
      group_by(feature) %>%
      summarise(
        strongest_abs_hedges_g = max(abs(hedges_g), na.rm = TRUE),
        min_group_q = min(p_fdr, na.rm = TRUE),
        strongest_contrast = contrast[which.max(abs(hedges_g))][1],
        .groups = "drop"
      ),
    by = "feature"
  ) %>%
  mutate(
    strongest_abs_hedges_g = ifelse(is.infinite(strongest_abs_hedges_g), NA_real_, strongest_abs_hedges_g),
    min_group_q = ifelse(is.infinite(min_group_q), NA_real_, min_group_q),
    ClaimStatus = case_when(
      ClaimTier == "primary" & !is.na(min_group_q) & min_group_q < 0.10 ~ "primary_descriptive_signal",
      ClaimTier == "primary" ~ "primary_feature_no_strong_group_signal",
      ClaimTier == "secondary" ~ "secondary_composite_or_mechanistic_context",
      TRUE ~ "exploratory"
    ),
    InterpretationRule = "RES/SUS group contrasts are descriptive because groups are post-paradigm CombZ-derived; prospective claims require pre-endpoint features only."
  )
write_table(systems_primary_claims_table, file.path(output_dir, "tables/systems_primary_claims_table.csv"))

primary_model_dat <- systems_features %>%
  select(AnimalNum, Group, Sex, any_of(c("Batch", "System")), all_of(unique(primary_feature_registry$feature))) %>%
  mutate(
    Group = factor(as.character(Group), levels = group_levels),
    Sex = factor(as.character(Sex), levels = sex_levels),
    Batch = factor(if_else(is.na(Batch) | Batch == "", "unknown_batch", as.character(Batch))),
    System = factor(if_else(is.na(System) | System == "", "unknown_system", as.character(System)))
  )

group_balance_by_batch_system <- systems_features %>%
  distinct(AnimalNum, Group, Sex, Batch, System) %>%
  mutate(
    Batch = if_else(is.na(Batch) | Batch == "", "unknown_batch", as.character(Batch)),
    System = if_else(is.na(System) | System == "", "unknown_system", as.character(System))
  ) %>%
  count(Batch, System, Sex, Group, name = "n")
write_table(group_balance_by_batch_system, file.path(output_dir, "tables/qc/group_balance_by_batch_system.csv"))

if (nrow(group_balance_by_batch_system) > 0) {
  p_group_balance <- group_balance_by_batch_system %>%
    ggplot(aes(Group, n, fill = Group)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.18) +
    facet_grid(Sex + Batch ~ System, scales = "free_y") +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    labs(title = "Group balance by batch and system", x = NULL, y = "Animals") +
    make_nature_theme(base_size = 6) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
  save_plot_svg_pdf(p_group_balance, file.path(output_dir, "figures/qc/Fig_group_balance_by_batch_system"), width = 170, height = 125)
}

batch_system_cage_audit <- map_dfr(unique(primary_feature_registry$feature), function(fc) {
  dat <- primary_model_dat %>% filter(is.finite(.data[[fc]]), !is.na(Group), !is.na(Sex))
  fd <- primary_feature_registry %>% filter(feature == fc) %>% slice_head(n = 1)
  if (nrow(dat) < 8 || n_distinct(dat$Group) < 2 || n_distinct(dat$Sex) < 1) {
    return(tibble(feature = fc, PrimaryFeatureLabel = fd$PrimaryFeatureLabel[1], n = nrow(dat), model_used = "not_estimable", group_p_uncovariate = NA_real_, group_p_covariate = NA_real_, batch_p = NA_real_, system_p = NA_real_, effect_size_without_covariates = NA_real_, effect_size_with_covariates = NA_real_, percent_change = NA_real_, interpretation_flag = "not_estimable"))
  }
  formula_uncov <- as.formula(paste0("`", fc, "` ~ Group * Sex"))
  use_system <- n_distinct(dat$System) > 1 && min(table(dat$System)) >= 2
  use_batch <- n_distinct(dat$Batch) > 1 && min(table(dat$Batch)) >= 2
  covars <- c(if (use_batch) "Batch", if (use_system) "System")
  formula_cov <- as.formula(paste0("`", fc, "` ~ Group * Sex", if (length(covars) > 0) paste0(" + ", paste(covars, collapse = " + ")) else ""))
  fit0 <- try(lm(formula_uncov, data = dat), silent = TRUE)
  fit1 <- try(lm(formula_cov, data = dat), silent = TRUE)
  eta_group <- function(fit) {
    if (inherits(fit, "try-error")) return(NA_real_)
    a <- try(anova(fit), silent = TRUE)
    if (inherits(a, "try-error") || !"Group" %in% rownames(a)) return(NA_real_)
    a["Group", "Sum Sq"] / sum(a[, "Sum Sq"], na.rm = TRUE)
  }
  eff0 <- eta_group(fit0)
  eff1 <- eta_group(fit1)
  tibble(
    feature = fc,
    PrimaryFeatureLabel = fd$PrimaryFeatureLabel[1],
    n = nrow(dat),
    model_used = deparse(formula_cov),
    group_p_uncovariate = safe_model_p(formula_uncov, dat, "^Group"),
    group_p_covariate = safe_model_p(formula_cov, dat, "^Group"),
    batch_p = if (use_batch) safe_model_p(formula_cov, dat, "^Batch") else NA_real_,
    system_p = if (use_system) safe_model_p(formula_cov, dat, "^System") else NA_real_,
    effect_size_without_covariates = eff0,
    effect_size_with_covariates = eff1,
    percent_change = 100 * (eff1 - eff0) / pmax(abs(eff0), 1e-9),
    interpretation_flag = case_when(
      nrow(dat) < 12 ~ "underpowered",
      is.finite(batch_p) & batch_p < 0.10 & abs(percent_change) > 30 ~ "possibly_batch-confounded",
      is.finite(system_p) & system_p < 0.10 & abs(percent_change) > 30 ~ "possibly_batch-confounded",
      is.finite(percent_change) & abs(percent_change) <= 30 ~ "robust",
      TRUE ~ "not_estimable"
    )
  )
}) %>%
  mutate(
    group_q_covariate = p.adjust(group_p_covariate, method = "BH"),
    batch_q = p.adjust(batch_p, method = "BH"),
    system_q = p.adjust(system_p, method = "BH")
  )
write_table(batch_system_cage_audit, file.path(output_dir, "tables/systems_batch_system_cage_audit.csv"))

if (nrow(batch_system_cage_audit) > 0) {
  p_batch_bias <- batch_system_cage_audit %>%
    mutate(
      PlotFeatureLabel = make.unique(str_trunc(PrimaryFeatureLabel, 42)),
      PlotFeatureLabel = factor(PlotFeatureLabel, levels = rev(PlotFeatureLabel))
    ) %>%
    ggplot(aes(percent_change, PlotFeatureLabel, fill = interpretation_flag)) +
    geom_vline(xintercept = c(-30, 30), linewidth = 0.2, linetype = "dashed", colour = "grey45") +
    geom_col(width = 0.65, colour = "white", linewidth = 0.15) +
    labs(title = "Batch/system sensitivity of primary feature group effects", x = "% change in group effect after covariates", y = NULL, fill = "Flag") +
    make_nature_theme(base_size = 6)
  save_plot_svg_pdf(p_batch_bias, file.path(output_dir, "figures/qc/Fig_batch_system_feature_bias"), width = 150, height = 95)
}

interaction_models <- map_dfr(unique(primary_feature_registry$feature), function(fc) {
  dat <- primary_model_dat %>% filter(is.finite(.data[[fc]]), !is.na(Group), !is.na(Sex))
  fd <- primary_feature_registry %>% filter(feature == fc) %>% slice_head(n = 1)
  if (nrow(dat) < 8 || n_distinct(dat$Group) < 2 || n_distinct(dat$Sex) < 2) {
    return(tibble(feature = fc, PrimaryFeatureLabel = fd$PrimaryFeatureLabel[1], term = NA_character_, estimate = NA_real_, p.value = NA_real_, model_component = "not_estimable", n = nrow(dat)))
  }
  use_system <- n_distinct(dat$System) > 1 && min(table(dat$System)) >= 2
  use_batch <- n_distinct(dat$Batch) > 1 && min(table(dat$Batch)) >= 2
  covars <- c(if (use_batch) "Batch", if (use_system) "System")
  form <- as.formula(paste0("`", fc, "` ~ Group * Sex", if (length(covars) > 0) paste0(" + ", paste(covars, collapse = " + ")) else ""))
  model_string <- paste(deparse(form), collapse = " ")
  fit <- try(lm(form, data = dat), silent = TRUE)
  if (inherits(fit, "try-error")) return(tibble(feature = fc, PrimaryFeatureLabel = fd$PrimaryFeatureLabel[1], term = NA_character_, estimate = NA_real_, p.value = NA_real_, model_component = "not_estimable", n = nrow(dat)))
  cf <- coef(summary(fit)) %>% as.data.frame() %>% rownames_to_column("term") %>% as_tibble()
  names(cf) <- str_replace_all(names(cf), fixed("Pr(>|t|)"), "p.value")
  names(cf) <- str_replace_all(names(cf), fixed("Estimate"), "estimate")
  coef_tbl <- cf %>%
    transmute(
      feature = fc,
      PrimaryFeatureLabel = fd$PrimaryFeatureLabel[1],
      term,
      estimate,
      p.value,
      model_component = case_when(
        str_detect(term, ":") ~ "Group_x_Sex_interaction",
        str_detect(term, "^Group") ~ "Group_main_or_contrast_term",
        str_detect(term, "^Sex") ~ "Sex_main_effect",
        TRUE ~ "covariate_or_intercept"
      ),
      n = nrow(dat),
      model = model_string
    )
  emm_tbl <- dat %>%
    group_by(Group, Sex) %>%
    summarise(estimated_marginal_mean = mean(.data[[fc]], na.rm = TRUE), sem = sem(.data[[fc]]), n_cell = n(), .groups = "drop") %>%
    mutate(
      feature = fc,
      PrimaryFeatureLabel = fd$PrimaryFeatureLabel[1],
      term = paste("EMM", Group, Sex, sep = "_"),
      estimate = estimated_marginal_mean,
      p.value = NA_real_,
      model_component = "estimated_marginal_mean",
      n = n_cell,
      model = model_string
    ) %>%
    select(feature, PrimaryFeatureLabel, term, estimate, p.value, model_component, n, model)
  sex_contrasts <- map_dfr(levels(dat$Sex), function(sx) {
    dd <- dat %>% filter(Sex == sx)
    map_dfr(contrast_pairs, function(pair) {
      ref <- pair[1]; comp <- pair[2]
      x <- dd[[fc]][as.character(dd$Group) == ref]
      y <- dd[[fc]][as.character(dd$Group) == comp]
      tibble(feature = fc, PrimaryFeatureLabel = fd$PrimaryFeatureLabel[1], term = paste0(comp, "-", ref, "_within_", sx), estimate = mean(y, na.rm = TRUE) - mean(x, na.rm = TRUE), p.value = tryCatch(t.test(y, x)$p.value, error = function(e) NA_real_), model_component = "within_sex_group_contrast", n = sum(is.finite(x)) + sum(is.finite(y)), model = model_string)
    })
  })
  bind_rows(coef_tbl, emm_tbl, sex_contrasts)
}) %>%
  mutate(
    p_holm = p.adjust(p.value, method = "holm"),
    p_fdr = p.adjust(p.value, method = "BH"),
    SexSpecificClaimRule = case_when(
      model_component == "Group_x_Sex_interaction" & !is.na(p_fdr) & p_fdr < 0.10 ~ "sex-specific claim supported by interaction",
      model_component == "within_sex_group_contrast" ~ "sex-stratified descriptive result unless interaction is supported",
      TRUE ~ "not a sex-specific claim"
    )
  )
write_table(interaction_models, file.path(output_dir, "stats_tables/systems_group_sex_interaction_models.csv"))

if (nrow(interaction_models) > 0) {
  p_interaction <- interaction_models %>%
    filter(model_component == "Group_x_Sex_interaction", !is.na(estimate)) %>%
    mutate(
      PlotFeatureLabel = make.unique(str_trunc(PrimaryFeatureLabel, 42)),
      PlotFeatureLabel = factor(PlotFeatureLabel, levels = rev(unique(PlotFeatureLabel)))
    ) %>%
    ggplot(aes(estimate, PlotFeatureLabel, fill = p_fdr < 0.10)) +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey55") +
    geom_col(width = 0.65, colour = "white", linewidth = 0.15) +
    labs(title = "Group x sex interaction screen", subtitle = "Sex-specific claims require supported interaction or descriptive labeling", x = "Interaction coefficient", y = NULL, fill = "FDR < 0.10") +
    make_nature_theme(base_size = 6)
  save_plot_svg_pdf(p_interaction, file.path(output_dir, "figures/publication_panels/Fig_group_sex_interaction_effects"), width = 145, height = 95)
}

# ------------------------------------------------
# INTERPRETABILITY LAYER: MODULE SCORECARDS + NAMED BIOLOGICAL SCORES
# ------------------------------------------------

module_scorecards_base <- feature_qc %>%
  filter(feature %in% usable_features) %>%
  group_by(Module) %>%
  summarise(
    n_features = n(),
    n_duration_robust = sum(!RawDurationSensitive, na.rm = TRUE),
    n_raw_duration_sensitive = sum(RawDurationSensitive, na.rm = TRUE),
    median_missing_fraction = median(missing_fraction, na.rm = TRUE),
    median_feature_sd = median(sd, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    group_contrasts %>%
      filter(status == "tested", is.finite(hedges_g)) %>%
      group_by(Module) %>%
      summarise(
        strongest_abs_hedges_g = max(abs(hedges_g), na.rm = TRUE),
        median_abs_hedges_g = median(abs(hedges_g), na.rm = TRUE),
        n_fdr_supported = sum(p_fdr < 0.05, na.rm = TRUE),
        strongest_contrast = contrast[which.max(abs(hedges_g))][1],
        strongest_sex = as.character(Sex[which.max(abs(hedges_g))][1]),
        .groups = "drop"
      ),
    by = "Module"
  ) %>%
  mutate(
    DurationRobustness = case_when(
      n_features == 0 ~ "no_features",
      n_raw_duration_sensitive == 0 ~ "all_duration_robust",
      n_duration_robust >= n_raw_duration_sensitive ~ "mixed_mostly_robust",
      TRUE ~ "duration_sensitive_review_needed"
    ),
    BiologicalReadout = case_when(
      Module == "Magnitude" ~ "Overall locomotor, spatial and social magnitude",
      Module == "Temporal instability" ~ "Local volatility, persistence and temporal organization",
      Module == "Latent-state organization" ~ "State occupancy, dwell time, switching and flexibility",
      Module == "Social topology" ~ "Engagement, isolation, centrality, fragmentation and partner structure",
      Module == "Trajectory geometry" ~ "Adaptation, recovery, peak timing and trajectory shape",
      Module == "Nonlinear systems dynamics" ~ "Complexity, attractor-like trapping, recurrence and early-warning dynamics",
      Module == "Predictive systems integration" ~ "Named biological composites and prediction-ready module scores",
      TRUE ~ "Other interpretable behavioral feature"
    )
  ) %>%
  arrange(factor(Module, levels = c(
    "Magnitude", "Temporal instability", "Latent-state organization", "Social topology",
    "Trajectory geometry", "Nonlinear systems dynamics", "Predictive systems integration",
    "Other interpretable feature"
  )))

write_table(module_scorecards_base, file.path(output_dir, "tables/systems_module_scorecards_base.csv"))
write_table(module_scorecards_base, file.path(output_dir, "tables/systems_module_scorecards.csv"))

feature_redundancy_tbl <- map_dfr(unique(feature_dictionary$Module), function(mod) {
  feats <- feature_dictionary %>%
    filter(Module == mod, feature %in% duration_robust_features) %>%
    pull(feature)
  feats <- intersect(feats, names(systems_features))
  if (length(feats) < 2) return(tibble())
  mat <- systems_features %>%
    select(all_of(feats)) %>%
    mutate(across(everything(), ~ replace_na(safe_numeric(.x), median(safe_numeric(.x), na.rm = TRUE)))) %>%
    as.matrix()
  cmat <- suppressWarnings(cor(mat, method = "spearman", use = "pairwise.complete.obs"))
  idx <- which(upper.tri(cmat), arr.ind = TRUE)
  tibble(
    Module = mod,
    Feature1 = colnames(cmat)[idx[, 1]],
    Feature2 = colnames(cmat)[idx[, 2]],
    spearman_rho = cmat[idx]
  ) %>%
    filter(is.finite(spearman_rho), abs(spearman_rho) >= 0.80) %>%
    mutate(
      redundancy_level = case_when(
        abs(spearman_rho) >= 0.95 ~ "near_duplicate",
        abs(spearman_rho) >= 0.90 ~ "very_high",
        TRUE ~ "high"
      )
    )
})

write_table(feature_redundancy_tbl, file.path(output_dir, "tables/systems_feature_redundancy_within_modules.csv"))

named_score_stats <- named_biological_scores_long %>%
  filter(is.finite(Score)) %>%
  group_by(Sex, ScoreLabel, Group) %>%
  summarise(
    n = n_distinct(AnimalNum),
    mean = mean(Score, na.rm = TRUE),
    ci_low = mean_ci(Score)["low"],
    ci_high = mean_ci(Score)["high"],
    median = median(Score, na.rm = TRUE),
    .groups = "drop"
  )

write_table(named_score_stats, file.path(output_dir, "stats_tables/systems_named_biological_score_group_summary.csv"))

if (nrow(named_biological_scores_long) > 0) {
  p_named_scores <- named_biological_scores_long %>%
    filter(is.finite(Score)) %>%
    mutate(ScoreLabel = factor(ScoreLabel, levels = unique(ScoreLabel))) %>%
    ggplot(aes(Group, Score, colour = Group, fill = Group)) +
    geom_violin(width = 0.82, alpha = 0.22, linewidth = 0.18, trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.70, linewidth = 0.20) +
    geom_point(position = position_jitter(width = 0.08, height = 0), size = 0.9, alpha = 0.75) +
    facet_grid(Sex ~ ScoreLabel, scales = "free_y") +
    scale_colour_manual(values = group_colors, drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    labs(
      title = "Named systems-level biological scores",
      subtitle = "Composite z-scores translate feature modules into interpretable behavioral constructs",
      x = NULL,
      y = "Composite score"
    ) +
    make_nature_theme(base_size = 6) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")

  save_plot_svg_pdf(p_named_scores, file.path(output_dir, "figures/publication_panels/Fig_systems_named_biological_scores"), width = 220, height = 125)

  if (requireNamespace("plotly", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
    ensure_dir(file.path(output_dir, "figures/interactive"))
    interactive_named <- named_biological_scores_long %>%
      filter(is.finite(Score)) %>%
      plotly::plot_ly(
        x = ~Score,
        y = ~ScoreLabel,
        color = ~Group,
        symbol = ~Sex,
        colors = group_colors,
        type = "scatter",
        mode = "markers",
        text = ~paste0("Animal: ", AnimalNum, "<br>Sex: ", Sex, "<br>Group: ", Group, "<br>Score: ", round(Score, 2)),
        hoverinfo = "text"
      ) %>%
      plotly::layout(
        title = "Interactive named biological scores",
        xaxis = list(title = "Composite score"),
        yaxis = list(title = "")
      )
    pandoc_available <- requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available()
    htmlwidgets::saveWidget(
      interactive_named,
      file.path(output_dir, "figures/interactive/systems_named_biological_scores.html"),
      selfcontained = pandoc_available
    )
  }
}

sleep_score_features <- named_biological_scores_long %>%
  filter(Metric %in% c(
    "sleep_like_fraction_score",
    "sleep_bout_duration_score",
    "sleep_fragmentation_score",
    "active_phase_quiescence_score",
    "inactive_phase_consolidation_score"
  ))

systems_sleep_like_inactivity_summary <- sleep_score_features %>%
  group_by(Sex, Group, Metric, ScoreLabel) %>%
  summarise(
    n = sum(is.finite(Score)),
    mean = mean(Score, na.rm = TRUE),
    ci_low = mean_ci(Score)["low"],
    ci_high = mean_ci(Score)["high"],
    median = median(Score, na.rm = TRUE),
    interpretation = "sleep-like inactivity / quiescence-like inactivity only; no EEG-validated sleep claim",
    .groups = "drop"
  )
write_table(systems_sleep_like_inactivity_summary, file.path(output_dir, "tables/systems_sleep_like_inactivity_summary.csv"))

if (nrow(sleep_score_features) > 0) {
  p_sleep_like <- sleep_score_features %>%
    filter(is.finite(Score)) %>%
    mutate(ScoreLabel = factor(ScoreLabel, levels = rev(unique(ScoreLabel)))) %>%
    ggplot(aes(Group, Score, colour = Group, fill = Group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.24, linewidth = 0.25) +
    geom_point(position = position_jitter(width = 0.10, height = 0), size = 1.25, alpha = 0.78, stroke = 0.18) +
    facet_grid(ScoreLabel ~ Sex, scales = "free_y") +
    scale_colour_manual(values = group_colors, drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    labs(
      title = "Sleep-like inactivity and quiescence-like organization",
      subtitle = "RFID-derived inactivity structure; not EEG-validated sleep",
      x = NULL,
      y = "Composite z-score"
    ) +
    make_nature_theme(base_size = 6) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
  save_plot_svg_pdf(p_sleep_like, file.path(output_dir, "figures/publication_panels/Fig_sleep_like_inactivity_by_group_sex"), width = 170, height = 125)
}

p_module_scorecard <- module_scorecards_base %>%
  mutate(Module = factor(Module, levels = rev(Module))) %>%
  ggplot(aes(strongest_abs_hedges_g, Module, fill = n_duration_robust)) +
  geom_col(width = 0.68, colour = "white", linewidth = 0.18) +
  geom_text(aes(label = paste0("n=", n_features, "\nrobust=", n_duration_robust)), hjust = -0.06, size = 1.85, lineheight = 0.85) +
  scale_fill_gradient(low = "#C6C3BB", high = "#3d3b6e", na.value = "grey85") +
  coord_cartesian(clip = "off") +
  labs(
    title = "Module scorecard",
    subtitle = "Strongest group effect per biological module with duration-robust feature count",
    x = "Strongest absolute Hedges g",
    y = NULL,
    fill = "Robust\nfeatures"
  ) +
  make_nature_theme(base_size = 6) +
  theme(plot.margin = margin(5.5, 18, 5.5, 5.5), legend.position = "right")

save_plot_svg_pdf(p_module_scorecard, file.path(output_dir, "figures/publication_panels/Fig_systems_module_scorecard"), width = 150, height = 88)

# ------------------------------------------------
# DIMENSIONALITY REDUCTION: PCA + OPTIONAL UMAP
# ------------------------------------------------

x <- systems_features %>% select(all_of(duration_robust_features))
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
    subtitle = paste0("PCA on ", length(duration_robust_features), " duration-robust descriptive features; n=", nrow(pca_scores), " animals"),
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
      subtitle = paste0("UMAP on duration-robust animal-level systems features; n=", nrow(umap_scores), "; n_neighbors=", n_neighbors),
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
  mutate(Phase = PhaseClass) %>%
  join_duration_qc(epoch_duration_qc, by_cols = c("AnimalNum", "Group", "Sex", "CageChange", "Phase")) %>%
  filter(n_bins >= 4, !short_epoch %in% TRUE)

latent_feature_cols <- latent_epoch_features %>%
  select(where(is.numeric)) %>%
  names()
latent_feature_cols <- setdiff(latent_feature_cols, c("CageChangeIndex", "n_bins"))
latent_feature_cols <- setdiff(latent_feature_cols, c("observed_bins", "expected_bins", "total_observation_duration_hours", "active_duration_hours", "inactive_duration_hours", "duration_completeness_fraction", "cage_change_duration_fraction"))
latent_feature_cols <- latent_feature_cols[
  map_lgl(latent_feature_cols, function(fc) {
    x <- safe_numeric(latent_epoch_features[[fc]])
    finite_x <- x[is.finite(x)]
    length(finite_x) >= 4 && isTRUE(sd(finite_x) > 0)
  })
]

latent_embedding_tbl <- tibble()
latent_trajectory_summary <- tibble()
latent_instability_by_animal <- tibble()
p_latent_traj <- NULL
p_umap_traj <- NULL
p_phate_traj <- NULL
p_instability <- NULL

if (nrow(latent_epoch_features) >= 20 && length(latent_feature_cols) >= 2) {
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
  phate_available <- requireNamespace("phateR", quietly = TRUE) &&
    requireNamespace("reticulate", quietly = TRUE) &&
    isTRUE(reticulate::py_module_available("phate"))
  if (phate_available && nrow(latent_mat_scaled) >= 30) {
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
      total_observation_duration_hours = sum(total_observation_duration_hours, na.rm = TRUE),
      latent_path_length_per_hour = latent_path_length / pmax(total_observation_duration_hours, 1e-9),
      latent_path_length_per_epoch = latent_path_length / pmax(n_epochs, 1),
      latent_roughness_normalized = latent_roughness / pmax(n_epochs, 1),
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
      roughness_rho = safe_cor(latent_roughness_normalized, .data[[primary_outcome]], "spearman"),
      roughness_p = safe_cor_p(latent_roughness_normalized, .data[[primary_outcome]], "spearman"),
      path_length_rho = safe_cor(latent_path_length_per_hour, .data[[primary_outcome]], "spearman"),
      path_length_p = safe_cor_p(latent_path_length_per_hour, .data[[primary_outcome]], "spearman"),
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
    filter(is.finite(latent_roughness_normalized)) %>%
    mutate(log_roughness = log10(latent_roughness_normalized + 1)) %>%
    ggplot(aes(Group, log_roughness, colour = Group, fill = Group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.28, linewidth = 0.28) +
    geom_point(position = position_jitter(width = 0.11, height = 0), size = 1.55, alpha = 0.78, stroke = 0.20) +
    facet_grid(PhaseClass ~ Sex, scales = "free_y") +
    scale_colour_manual(values = group_colors, drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    labs(
      title = "Latent trajectory instability",
      subtitle = "Duration-normalized path roughness through temporal state space by animal",
      x = NULL,
      y = "log10(normalized roughness + 1)"
    ) +
    make_nature_theme(base_size = 6)

  save_plot_svg_pdf(p_instability, file.path(output_dir, "figures/publication_panels/Fig_systems_latent_trajectory_instability"), width = 135, height = 95)
}

# ------------------------------------------------
# EFFECT-SIZE HEATMAPS
# ------------------------------------------------

focus_domains <- c(
  "behavior", "composite", "temporal_dynamics", "latent_space", "behavioral_state",
  "latent_state", "social_network", "trajectory", "trajectory_geometry",
  "nonlinear_dynamics", "systems_integration"
)
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
# MODERN INTERPRETABILITY VISUALS: HMM FLOW, SOCIAL MAPS, TRAJECTORY ADAPTATION
# ------------------------------------------------

hmm_dir <- file.path(project_root, "analysis_ready/06_behavioral_dynamics/hmm_states", primary_bin_level, "tables")
hmm_state_summary <- read_any_table(file.path(hmm_dir, "hmm_state_summary.csv"))
hmm_transition_prob <- read_any_table(file.path(hmm_dir, "hmm_transition_probabilities.csv"))
hmm_dwell <- read_any_table(file.path(hmm_dir, "hmm_state_dwell_times.csv"))
hmm_occupancy <- read_any_table(file.path(hmm_dir, "hmm_state_occupancy.csv"))

state_label_tbl <- if (!is.null(hmm_state_summary) && nrow(hmm_state_summary) > 0) {
  hmm_state_summary %>%
    mutate(
      State = as.character(State),
      Movement_z = safe_numeric(Movement_z),
      Entropy_z = safe_numeric(Entropy_z),
      Proximity_z = safe_numeric(Proximity_z),
      SemanticState = case_when(
        Movement_z <= median(Movement_z, na.rm = TRUE) & Entropy_z <= median(Entropy_z, na.rm = TRUE) ~ "inactive/low-exploration",
        Proximity_z >= quantile(Proximity_z, 0.67, na.rm = TRUE) ~ "social",
        Movement_z >= quantile(Movement_z, 0.67, na.rm = TRUE) ~ "burst/high-movement",
        Entropy_z >= quantile(Entropy_z, 0.67, na.rm = TRUE) ~ "exploratory",
        TRUE ~ "mixed"
      ),
      StateLabel = paste0("S", State, "\n", SemanticState)
    ) %>%
    distinct(State, StateLabel, SemanticState)
} else {
  tibble(State = character(), StateLabel = character(), SemanticState = character())
}

write_table(state_label_tbl, file.path(output_dir, "tables/systems_hmm_state_semantic_labels.csv"))

if (!is.null(hmm_transition_prob) && nrow(hmm_transition_prob) > 0) {
  hmm_transition_plot_tbl <- hmm_transition_prob %>%
    standardize_id_columns() %>%
    mutate(
      State = as.character(State),
      NextState = as.character(NextState),
      TransitionProbability = safe_numeric(TransitionProbability)
    ) %>%
    left_join(state_label_tbl %>% select(State, FromLabel = StateLabel), by = "State") %>%
    left_join(state_label_tbl %>% select(NextState = State, ToLabel = StateLabel), by = "NextState") %>%
    mutate(
      FromLabel = coalesce(FromLabel, paste0("S", State)),
      ToLabel = coalesce(ToLabel, paste0("S", NextState))
    ) %>%
    group_by(Sex, Group, FromLabel, ToLabel) %>%
    summarise(mean_probability = mean(TransitionProbability, na.rm = TRUE), .groups = "drop")

  transition_difference_wide <- hmm_transition_plot_tbl %>%
    filter(as.character(Group) %in% c("CON", "RES", "SUS")) %>%
    pivot_wider(names_from = Group, values_from = mean_probability)
  for (grp in c("CON", "RES", "SUS")) {
    if (!grp %in% names(transition_difference_wide)) transition_difference_wide[[grp]] <- NA_real_
  }
  transition_difference_tbl <- transition_difference_wide %>%
    mutate(
      `SUS - RES` = SUS - RES,
      `SUS - CON` = SUS - CON,
      `RES - CON` = RES - CON
    ) %>%
    pivot_longer(c(`SUS - RES`, `SUS - CON`, `RES - CON`), names_to = "Contrast", values_to = "DeltaProbability")

  write_table(hmm_transition_plot_tbl, file.path(output_dir, "tables/systems_hmm_transition_probability_summary.csv"))
  write_table(transition_difference_tbl, file.path(output_dir, "tables/systems_hmm_transition_probability_differences.csv"))

  if (nrow(transition_difference_tbl) > 0) {
    p_hmm_diff <- transition_difference_tbl %>%
      filter(is.finite(DeltaProbability)) %>%
      ggplot(aes(ToLabel, FromLabel, fill = DeltaProbability)) +
      geom_tile(colour = "white", linewidth = 0.22) +
      facet_grid(Sex ~ Contrast) +
      scale_fill_gradient2(low = "#3d3b6e", mid = "white", high = "#e63947", midpoint = 0, na.value = "grey90") +
      labs(
        title = "HMM transition architecture differences",
        subtitle = "Group differences in animal-level transition probabilities; positive values indicate higher probability in the first group",
        x = "To state",
        y = "From state",
        fill = "Delta P"
      ) +
      make_nature_theme(base_size = 6) +
      theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "right")

    save_plot_svg_pdf(p_hmm_diff, file.path(output_dir, "figures/publication_panels/Fig_systems_hmm_transition_difference"), width = 175, height = 110)
  }

  if (requireNamespace("ggalluvial", quietly = TRUE)) {
    hmm_flow_tbl <- hmm_transition_plot_tbl %>%
      group_by(Sex, Group) %>%
      slice_max(mean_probability, n = 8, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(Flow = paste(FromLabel, ToLabel, sep = " -> "))

    p_hmm_flow <- hmm_flow_tbl %>%
      ggplot(aes(axis1 = FromLabel, axis2 = ToLabel, y = mean_probability, fill = Group)) +
      ggalluvial::geom_alluvium(width = 0.18, alpha = 0.72) +
      ggalluvial::geom_stratum(width = 0.18, colour = "grey35", linewidth = 0.18, fill = "grey96") +
      ggplot2::geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 1.65, lineheight = 0.78) +
      facet_grid(Sex ~ Group) +
      scale_fill_manual(values = group_colors, drop = FALSE) +
      scale_x_discrete(limits = c("From", "To"), expand = c(0.08, 0.08)) +
      labs(
        title = "Dominant latent-state flows",
        subtitle = "Top HMM transitions by group and sex",
        x = NULL,
        y = "Mean transition probability"
      ) +
      make_nature_theme(base_size = 6) +
      theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())

    save_plot_svg_pdf(p_hmm_flow, file.path(output_dir, "figures/publication_panels/Fig_systems_hmm_state_flow_alluvial"), width = 180, height = 125)
  }
}

if (!is.null(hmm_dwell) && nrow(hmm_dwell) > 0) {
  dwell_plot_tbl <- hmm_dwell %>%
    standardize_id_columns() %>%
    mutate(State = as.character(State), mean_dwell_hours = safe_numeric(mean_dwell_hours)) %>%
    left_join(state_label_tbl %>% select(State, StateLabel), by = "State") %>%
    mutate(StateLabel = coalesce(StateLabel, paste0("S", State))) %>%
    filter(is.finite(mean_dwell_hours))

  if (requireNamespace("ggridges", quietly = TRUE) && nrow(dwell_plot_tbl) > 0) {
    p_dwell <- dwell_plot_tbl %>%
      ggplot(aes(mean_dwell_hours, StateLabel, fill = Group, colour = Group)) +
      ggridges::geom_density_ridges(alpha = 0.32, linewidth = 0.25, scale = 1.05, rel_min_height = 0.01) +
      facet_grid(Sex ~ Phase, scales = "free_y") +
      scale_colour_manual(values = group_colors, drop = FALSE) +
      scale_fill_manual(values = group_colors, drop = FALSE) +
      labs(
        title = "HMM dwell-time distributions",
        subtitle = "Long dwell times indicate state persistence or behavioral rigidity",
        x = "Mean dwell time (hours)",
        y = NULL
      ) +
      make_nature_theme(base_size = 6)

    save_plot_svg_pdf(p_dwell, file.path(output_dir, "figures/publication_panels/Fig_systems_hmm_dwell_time_ridges"), width = 170, height = 115)
  }
}

social_score_tbl <- named_biological_scores_long %>%
  select(AnimalNum, Group, Sex, Metric, Score) %>%
  filter(Metric %in% c("social_withdrawal_score", "social_fragmentation_score", "exploratory_flexibility_score")) %>%
  pivot_wider(names_from = Metric, values_from = Score)

if (nrow(social_score_tbl) > 0 && all(c("social_withdrawal_score", "social_fragmentation_score") %in% names(social_score_tbl))) {
  p_social_map <- social_score_tbl %>%
    mutate(flexibility_plot_size = replace_na(abs(exploratory_flexibility_score), 0)) %>%
    ggplot(aes(social_withdrawal_score, social_fragmentation_score, colour = Group, fill = Group)) +
    geom_hline(yintercept = 0, linewidth = 0.18, colour = "grey80") +
    geom_vline(xintercept = 0, linewidth = 0.18, colour = "grey80") +
    geom_point(aes(size = flexibility_plot_size), alpha = 0.82, stroke = 0.25) +
    stat_ellipse(aes(group = Group), linewidth = 0.35, alpha = 0.45, show.legend = FALSE) +
    facet_grid(. ~ Sex) +
    scale_colour_manual(values = group_colors, drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    scale_size_continuous(range = c(1.4, 4.0)) +
    labs(
      title = "Social phenotype map",
      subtitle = "Separates withdrawal from fragmented topology; point size reflects exploratory flexibility",
      x = "Social withdrawal score",
      y = "Social fragmentation score",
      size = "|Flexibility|"
    ) +
    make_nature_theme(base_size = 7) +
    theme(legend.position = "right")

  save_plot_svg_pdf(p_social_map, file.path(output_dir, "figures/publication_panels/Fig_systems_social_phenotype_map"), width = 140, height = 78)
}

trajectory_score_tbl <- named_biological_scores_long %>%
  select(AnimalNum, Group, Sex, Metric, Score) %>%
  filter(Metric %in% c("trajectory_recovery_score", "behavioral_rigidity_score", "stress_adaptation_index")) %>%
  pivot_wider(names_from = Metric, values_from = Score)

if (nrow(trajectory_score_tbl) > 0 && all(c("trajectory_recovery_score", "behavioral_rigidity_score") %in% names(trajectory_score_tbl))) {
  p_adaptation <- trajectory_score_tbl %>%
    mutate(adaptation_plot_size = replace_na(stress_adaptation_index, 0)) %>%
    ggplot(aes(behavioral_rigidity_score, trajectory_recovery_score, colour = Group, fill = Group)) +
    geom_hline(yintercept = 0, linewidth = 0.18, colour = "grey80") +
    geom_vline(xintercept = 0, linewidth = 0.18, colour = "grey80") +
    geom_point(aes(size = adaptation_plot_size), alpha = 0.84, stroke = 0.25) +
    geom_smooth(method = "lm", se = TRUE, colour = "grey25", fill = "grey70", linewidth = 0.35, alpha = 0.14) +
    facet_grid(. ~ Sex) +
    scale_colour_manual(values = group_colors, drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    scale_size_continuous(range = c(1.4, 4.2)) +
    labs(
      title = "Adaptation phase portrait",
      subtitle = "Recovery dynamics plotted against behavioral rigidity",
      x = "Behavioral rigidity score",
      y = "Trajectory recovery score",
      size = "Adaptation\nindex"
    ) +
    make_nature_theme(base_size = 7) +
    theme(legend.position = "right")

  save_plot_svg_pdf(p_adaptation, file.path(output_dir, "figures/publication_panels/Fig_systems_trajectory_adaptation_phase_portrait"), width = 140, height = 78)
}

# ------------------------------------------------
# MODULE-LEVEL PREDICTION LADDER HELPERS
# ------------------------------------------------

make_module_score <- function(dat, features, score_name) {
  features <- intersect(features, names(dat))
  if (length(features) == 0) return(tibble(AnimalNum = dat$AnimalNum, !!score_name := NA_real_))
  mat <- dat %>%
    select(all_of(features)) %>%
    mutate(across(everything(), safe_numeric)) %>%
    as.matrix()
  mat <- scale(mat)
  mat[!is.finite(mat)] <- NA_real_
  tibble(AnimalNum = dat$AnimalNum, !!score_name := rowMeans(mat, na.rm = TRUE))
}

loo_lm_module_predict <- function(dat, predictors, model_name) {
  predictors <- predictors[predictors %in% names(dat)]
  predictors <- predictors[map_lgl(predictors, ~ isTRUE(sd(dat[[.x]], na.rm = TRUE) > 0))]
  pred <- rep(NA_real_, nrow(dat))
  coef_tbl <- tibble()
  if (nrow(dat) < 6 || length(predictors) == 0) {
    return(list(predictions = dat %>% transmute(AnimalNum, Group, Sex, Model = model_name, outcome, predicted = NA_real_, residual = NA_real_),
                coefficients = coef_tbl))
  }

  for (i in seq_len(nrow(dat))) {
    train_idx <- setdiff(seq_len(nrow(dat)), i)
    train <- dat[train_idx, , drop = FALSE]
    test <- dat[i, , drop = FALSE]
    train_medians <- map_dbl(predictors, ~ median(train[[.x]], na.rm = TRUE))
    train_medians[!is.finite(train_medians)] <- 0
    for (j in seq_along(predictors)) {
      train[[predictors[j]]][!is.finite(train[[predictors[j]]])] <- train_medians[j]
      test[[predictors[j]]][!is.finite(test[[predictors[j]]])] <- train_medians[j]
    }
    form <- as.formula(paste("outcome ~", paste(predictors, collapse = " + ")))
    fit <- try(lm(form, data = train), silent = TRUE)
    if (inherits(fit, "try-error")) next
    pred[i] <- tryCatch(as.numeric(predict(fit, newdata = test)), error = function(e) NA_real_)
    cf <- coef(fit)
    coef_tbl <- bind_rows(coef_tbl, tibble(Model = model_name, held_out = dat$AnimalNum[i], term = names(cf), coefficient = as.numeric(cf)))
  }

  list(
    predictions = dat %>% transmute(AnimalNum, Group, Sex, Model = model_name, outcome, predicted = pred, residual = outcome - predicted),
    coefficients = coef_tbl
  )
}

prediction_performance_tbl <- function(pred_tbl, n_features, seed = 1) {
  ok <- pred_tbl %>% filter(is.finite(outcome), is.finite(predicted))
  if (nrow(ok) < 4) {
    return(tibble(n = nrow(ok), n_features = n_features, pearson_r = NA_real_, spearman_rho = NA_real_, rmse = NA_real_, mae = NA_real_, cv_r2 = NA_real_, permutation_p = NA_real_, pearson_ci_low = NA_real_, pearson_ci_high = NA_real_, spearman_ci_low = NA_real_, spearman_ci_high = NA_real_))
  }
  pear <- safe_cor(ok$outcome, ok$predicted, "pearson")
  spear <- safe_cor(ok$outcome, ok$predicted, "spearman")
  baseline <- mean(ok$outcome, na.rm = TRUE)
  set.seed(seed)
  perm <- replicate(n_prediction_permutations, abs(safe_cor(ok$outcome, sample(ok$predicted), "pearson")))
  set.seed(seed + 100)
  boot <- replicate(n_prediction_bootstrap, {
    idx <- sample(seq_len(nrow(ok)), replace = TRUE)
    c(pearson = safe_cor(ok$outcome[idx], ok$predicted[idx], "pearson"),
      spearman = safe_cor(ok$outcome[idx], ok$predicted[idx], "spearman"))
  })
  tibble(
    n = nrow(ok),
    n_features = n_features,
    pearson_r = pear,
    spearman_rho = spear,
    rmse = sqrt(mean((ok$outcome - ok$predicted)^2, na.rm = TRUE)),
    mae = mean(abs(ok$outcome - ok$predicted), na.rm = TRUE),
    cv_r2 = 1 - sum((ok$outcome - ok$predicted)^2, na.rm = TRUE) / sum((ok$outcome - baseline)^2, na.rm = TRUE),
    permutation_p = (sum(perm >= abs(pear), na.rm = TRUE) + 1) / (sum(is.finite(perm)) + 1),
    pearson_ci_low = quantile(boot["pearson", ], 0.025, na.rm = TRUE, names = FALSE),
    pearson_ci_high = quantile(boot["pearson", ], 0.975, na.rm = TRUE, names = FALSE),
    spearman_ci_low = quantile(boot["spearman", ], 0.025, na.rm = TRUE, names = FALSE),
    spearman_ci_high = quantile(boot["spearman", ], 0.975, na.rm = TRUE, names = FALSE)
  )
}

build_prediction_ladder <- function(outcome_to_plot, animal_filter = NULL, analysis_set = "full") {
  source_dat <- systems_features
  if (!is.null(animal_filter)) {
    source_dat <- source_dat %>% semi_join(tibble(AnimalNum = animal_filter), by = "AnimalNum")
  }

  fd <- feature_dictionary %>% filter(feature %in% duration_robust_features)
  early_raw <- fd %>% filter(Source == "raw", Context == prospective_prediction_context)
  movement_features <- early_raw %>% filter(Metric == "movement") %>% pull(feature)
  proximity_features <- early_raw %>% filter(Metric == "proximity") %>% pull(feature)
  entropy_acf_features <- early_raw %>% filter(Metric == "entropy", Statistic == "acf1") %>% pull(feature)
  instability_features <- fd %>% filter(Module == "Temporal instability", feature %in% prospective_prediction_features | Source %in% c("burstiness", "state_space")) %>% pull(feature)
  latent_features <- fd %>% filter(Module == "Latent-state organization") %>% pull(feature)
  social_features <- fd %>% filter(Module == "Social topology") %>% pull(feature)
  trajectory_features <- fd %>% filter(Module == "Trajectory geometry") %>% pull(feature)
  nonlinear_features <- fd %>% filter(Module == "Nonlinear systems dynamics") %>% pull(feature)

  score_tbl <- source_dat %>% select(AnimalNum, Group, Sex, all_of(outcome_to_plot)) %>% rename(outcome = all_of(outcome_to_plot))
  score_tbl <- reduce(list(
    make_module_score(source_dat, movement_features, "movement_score"),
    make_module_score(source_dat, proximity_features, "proximity_score"),
    make_module_score(source_dat, entropy_acf_features, "entropy_acf1_score"),
    make_module_score(source_dat, instability_features, "instability_score"),
    make_module_score(source_dat, latent_features, "latent_state_score"),
    make_module_score(source_dat, social_features, "social_topology_score"),
    make_module_score(source_dat, trajectory_features, "trajectory_geometry_score"),
    make_module_score(source_dat, nonlinear_features, "nonlinear_dynamics_score")
  ), left_join, by = "AnimalNum") %>%
    right_join(score_tbl, by = "AnimalNum") %>%
    select(AnimalNum, Group, Sex, outcome, everything()) %>%
    filter(is.finite(outcome))

  ladder <- tibble(
    ModelOrder = 1:7,
    Model = c(
      "movement-only",
      "movement + proximity",
      "movement + entropy_acf1",
      "movement + instability metrics",
      "movement + latent-state metrics",
      "movement + social-network metrics",
      "integrated systems model"
    ),
    Predictors = list(
      "movement_score",
      c("movement_score", "proximity_score"),
      c("movement_score", "proximity_score", "entropy_acf1_score"),
      c("movement_score", "proximity_score", "entropy_acf1_score", "instability_score"),
      c("movement_score", "proximity_score", "entropy_acf1_score", "instability_score", "latent_state_score"),
      c("movement_score", "proximity_score", "entropy_acf1_score", "instability_score", "latent_state_score", "social_topology_score"),
      c("movement_score", "proximity_score", "entropy_acf1_score", "instability_score", "latent_state_score", "social_topology_score", "trajectory_geometry_score", "nonlinear_dynamics_score")
    )
  )

  fits <- pmap(ladder, function(ModelOrder, Model, Predictors) {
    loo_lm_module_predict(score_tbl, Predictors, Model)
  })
  pred_tbl <- map_dfr(fits, "predictions")
  coef_tbl <- map_dfr(fits, "coefficients")
  perf_tbl <- map2_dfr(fits, seq_len(nrow(ladder)), function(fit, idx) {
    prediction_performance_tbl(fit$predictions, length(ladder$Predictors[[idx]]), seed = 100 + idx) %>%
      mutate(ModelOrder = ladder$ModelOrder[idx], Model = ladder$Model[idx], Predictors = paste(ladder$Predictors[[idx]], collapse = " + "))
  }) %>%
    arrange(ModelOrder) %>%
    mutate(
      delta_cv_r2_vs_previous = cv_r2 - lag(cv_r2),
      delta_cv_r2_vs_movement = cv_r2 - cv_r2[Model == "movement-only"][1],
      DurationAnalysisSet = analysis_set,
      DurationSensitivityStatus = "fit"
    )

  pred_tbl <- pred_tbl %>% mutate(DurationAnalysisSet = analysis_set)
  coef_tbl <- coef_tbl %>% mutate(DurationAnalysisSet = analysis_set)
  score_tbl <- score_tbl %>% mutate(DurationAnalysisSet = analysis_set)

  list(scores = score_tbl, predictions = pred_tbl, coefficients = coef_tbl, performance = perf_tbl)
}

# Module-score coupling uses module aggregates rather than hundreds of individual features.
module_score_matrix <- {
  fd <- feature_dictionary %>% filter(feature %in% duration_robust_features)
  module_features <- split(fd$feature, fd$Module)
  score_parts <- imap(module_features, function(features, module_name) {
    make_module_score(systems_features, features, paste0("module_score__", safe_name(module_name)))
  })
  if (length(score_parts) == 0) {
    systems_features %>% select(AnimalNum, Group, Sex)
  } else {
    reduce(score_parts, left_join, by = "AnimalNum") %>%
    left_join(systems_features %>% select(AnimalNum, Group, Sex), by = "AnimalNum") %>%
    relocate(AnimalNum, Group, Sex)
  }
}

write_table(module_score_matrix, file.path(output_dir, "tables/systems_module_scores_by_animal.csv"))

module_score_cols <- names(module_score_matrix)[str_detect(names(module_score_matrix), "^module_score__")]
module_coupling_tbl <- map_dfr(levels(droplevels(systems_features$Sex)), function(sx) {
  dat <- module_score_matrix %>% filter(as.character(Sex) == sx)
  if (nrow(dat) < 5 || length(module_score_cols) < 2) return(tibble())
  mat <- dat %>%
    select(all_of(module_score_cols)) %>%
    mutate(across(everything(), ~ replace_na(safe_numeric(.x), median(safe_numeric(.x), na.rm = TRUE)))) %>%
    as.matrix()
  cmat <- suppressWarnings(cor(mat, method = "spearman", use = "pairwise.complete.obs"))
  idx <- which(upper.tri(cmat), arr.ind = TRUE)
  tibble(
    Sex = sx,
    Module1 = str_remove(colnames(cmat)[idx[, 1]], "^module_score__") %>% str_replace_all("_", " ") %>% str_to_sentence(),
    Module2 = str_remove(colnames(cmat)[idx[, 2]], "^module_score__") %>% str_replace_all("_", " ") %>% str_to_sentence(),
    spearman_rho = cmat[idx]
  ) %>%
    filter(is.finite(spearman_rho))
})

write_table(module_coupling_tbl, file.path(output_dir, "tables/systems_module_coupling_edges.csv"))

if (nrow(module_coupling_tbl) > 0) {
  p_module_coupling <- module_coupling_tbl %>%
    filter(abs(spearman_rho) >= 0.25) %>%
    mutate(Pair = paste(Module1, Module2, sep = " <-> "),
           Pair = factor(Pair, levels = unique(Pair[order(abs(spearman_rho), decreasing = TRUE)]))) %>%
    ggplot(aes(spearman_rho, Pair, fill = spearman_rho)) +
    geom_col(width = 0.68) +
    facet_grid(Sex ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = "#3d3b6e", mid = "white", high = "#e63947", midpoint = 0) +
    labs(
      title = "Module-level coupling map",
      subtitle = "Spearman correlations among biologically interpretable module scores",
      x = "Spearman rho",
      y = NULL,
      fill = "rho"
    ) +
    make_nature_theme(base_size = 6) +
    theme(legend.position = "right")

  save_plot_svg_pdf(p_module_coupling, file.path(output_dir, "figures/publication_panels/Fig_systems_module_coupling_network"), width = 150, height = 108)
}

duration_negative_control_tbl <- if (nrow(epoch_duration_qc) > 0) {
  duration_qc_for_negative_control <- epoch_duration_qc
  if (!"short_regrouping_epoch" %in% names(duration_qc_for_negative_control)) {
    duration_qc_for_negative_control <- duration_qc_for_negative_control %>% mutate(short_regrouping_epoch = FALSE)
  }
  animal_duration_summary <- duration_qc_for_negative_control %>%
    group_by(AnimalNum) %>%
    summarise(
      mean_duration_completeness = mean(duration_completeness_fraction, na.rm = TRUE),
      min_duration_completeness = min(duration_completeness_fraction, na.rm = TRUE),
      fraction_short_epochs = mean(short_epoch %in% TRUE | cage_change_duration_class == "short" | short_regrouping_epoch %in% TRUE, na.rm = TRUE),
      total_observation_hours = sum(total_observation_duration_hours, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(across(c(min_duration_completeness, mean_duration_completeness), ~ if_else(is.infinite(.x), NA_real_, .x))) %>%
    left_join(systems_features %>% select(AnimalNum, Group, Sex, any_of(endpoint_cols)), by = "AnimalNum")

  duration_cols <- c("mean_duration_completeness", "min_duration_completeness", "fraction_short_epochs", "total_observation_hours")
  endpoint_available_for_qc <- intersect(endpoint_cols, names(animal_duration_summary))
  write_table(animal_duration_summary, file.path(output_dir, "tables/systems_duration_negative_control_by_animal.csv"))

  map_dfr(endpoint_available_for_qc, function(outcome) {
    y <- safe_numeric(animal_duration_summary[[outcome]])
    map_dfr(duration_cols, function(dc) {
      x <- safe_numeric(animal_duration_summary[[dc]])
      tibble(
        outcome = outcome,
        duration_metric = dc,
        n = sum(is.finite(x) & is.finite(y)),
        spearman_rho = safe_cor(x, y, "spearman"),
        spearman_p = safe_cor_p(x, y, "spearman"),
        reviewer_interpretation = "negative_control_duration_metric_should_not_drive_endpoint_association"
      )
    })
  }) %>%
    group_by(outcome) %>%
    mutate(spearman_fdr = p.adjust(spearman_p, method = "BH")) %>%
    ungroup()
} else {
  tibble()
}

write_table(duration_negative_control_tbl, file.path(output_dir, "stats_tables/systems_duration_negative_control_endpoint_correlations.csv"))

# ------------------------------------------------
# OUTCOME INTEGRATION AND PREDICTION SUMMARY
# ------------------------------------------------

available_outcomes <- intersect(endpoint_cols, names(systems_features))
available_outcomes <- available_outcomes[map_lgl(available_outcomes, ~ any(is.finite(safe_numeric(systems_features[[.x]]))))]

if (length(available_outcomes) > 0) {
  calc_outcome_assoc <- function(feature_set, feature_set_name) {
    empty_assoc <- tibble(
      outcome = character(),
      FeatureSet = character(),
      feature = character(),
      n = integer(),
      spearman_rho = numeric(),
      spearman_ci_low = numeric(),
      spearman_ci_high = numeric(),
      spearman_p = numeric(),
      pearson_r = numeric(),
      pearson_ci_low = numeric(),
      pearson_ci_high = numeric(),
      pearson_p = numeric(),
      spearman_fdr = numeric(),
      sig = character(),
      ReportingCorrection = character(),
      TestMethod = character(),
      ConfidenceInterval = character(),
      Direction = character()
    )
    if (length(feature_set) == 0) {
      return(empty_assoc %>% left_join(feature_dictionary, by = "feature"))
    }
    assoc <- map_dfr(available_outcomes, function(outcome) {
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
    })
    if (nrow(assoc) == 0) {
      return(empty_assoc %>% left_join(feature_dictionary, by = "feature"))
    }
    assoc %>%
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
           DisplayFeature = make.unique(str_trunc(DisplayFeature, 52)),
           DisplayFeature = factor(DisplayFeature, levels = rev(DisplayFeature)),
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
             DisplayFeature = make.unique(str_trunc(DisplayFeature, 52)),
             DisplayFeature = factor(DisplayFeature, levels = rev(DisplayFeature)),
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
        DisplayFeature = make.unique(str_trunc(DisplayFeature, 46)),
        DisplayFeature = factor(DisplayFeature, levels = unique(DisplayFeature))
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

  animal_duration_source <- epoch_duration_qc
  if (!"short_regrouping_epoch" %in% names(animal_duration_source)) {
    animal_duration_source <- animal_duration_source %>% mutate(short_regrouping_epoch = FALSE)
  }
  animal_duration_flags <- animal_duration_source %>%
    group_by(AnimalNum) %>%
    summarise(
      contains_short_duration_epoch = any(short_epoch %in% TRUE | cage_change_duration_class == "short" | short_regrouping_epoch %in% TRUE, na.rm = TRUE),
      .groups = "drop"
    )
  animals_without_short_duration <- systems_features %>%
    left_join(animal_duration_flags, by = "AnimalNum") %>%
    filter(!contains_short_duration_epoch %in% TRUE) %>%
    pull(AnimalNum)

  prediction_ladder <- build_prediction_ladder(outcome_to_plot, analysis_set = "full")
  prediction_ladder_no_short <- if (length(animals_without_short_duration) >= 8) {
    build_prediction_ladder(outcome_to_plot, animal_filter = animals_without_short_duration, analysis_set = "excluding_short_duration")
  } else {
    skipped_perf <- prediction_ladder$performance %>%
      mutate(
        across(any_of(c("n", "n_features", "pearson_r", "spearman_rho", "rmse", "mae", "cv_r2", "permutation_p",
                        "pearson_ci_low", "pearson_ci_high", "spearman_ci_low", "spearman_ci_high",
                        "delta_cv_r2_vs_previous", "delta_cv_r2_vs_movement")), ~ NA_real_),
        DurationAnalysisSet = "excluding_short_duration",
        DurationSensitivityStatus = "skipped_too_few_animals"
      )
    list(scores = tibble(), predictions = tibble(), coefficients = tibble(), performance = skipped_perf)
  }

  prediction_tbl_all_models_duration_sensitivity <- bind_rows(
    prediction_ladder$predictions,
    prediction_ladder_no_short$predictions
  )
  prediction_perf_duration_sensitivity <- bind_rows(
    prediction_ladder$performance,
    prediction_ladder_no_short$performance
  ) %>%
    group_by(Model) %>%
    mutate(
      full_cv_r2 = cv_r2[DurationAnalysisSet == "full"][1],
      full_pearson_r = pearson_r[DurationAnalysisSet == "full"][1],
      full_rmse = rmse[DurationAnalysisSet == "full"][1],
      delta_cv_r2_vs_full = cv_r2 - full_cv_r2,
      delta_pearson_r_vs_full = pearson_r - full_pearson_r,
      delta_rmse_vs_full = rmse - full_rmse
    ) %>%
    ungroup() %>%
    select(-full_cv_r2, -full_pearson_r, -full_rmse)

  prediction_tbl_all_models <- prediction_ladder$predictions
  prediction_tbl <- prediction_tbl_all_models %>% filter(Model == "integrated systems model")
  prediction_perf <- prediction_ladder$performance %>%
    mutate(
      outcome = outcome_to_plot,
      loo_pearson_r = pearson_r,
      loo_spearman_rho = spearman_rho,
      loo_rmse = rmse,
      cv_r2_vs_mean = cv_r2,
      pearson_r_permutation_p = permutation_p,
      ModelLeakageClass = "safe_for_prediction",
      ModelTemporalRule = "Uses pre-endpoint module scores built from first-active-window and duration-robust behavioral summaries; RES/SUS labels are not predictors.",
      ModelFamily = "Module-level linear model ladder, leave-one-animal-out cross-validation",
      Preprocessing = "Features are collapsed into biologically interpretable z-scored module scores before modeling"
    )
  if (!"delta_cv_r2_vs_previous" %in% names(prediction_perf)) {
    prediction_perf <- prediction_perf %>%
      arrange(ModelOrder) %>%
      mutate(delta_cv_r2_vs_previous = cv_r2 - lag(cv_r2))
  }
  if (!"delta_cv_r2_vs_movement" %in% names(prediction_perf)) {
    prediction_perf <- prediction_perf %>%
      arrange(ModelOrder) %>%
      mutate(delta_cv_r2_vs_movement = cv_r2 - cv_r2[Model == "movement-only"][1])
  }
  prediction_perf_integrated <- prediction_perf %>% filter(Model == "integrated systems model") %>% slice_head(n = 1)
  prediction_coef_tbl <- prediction_ladder$coefficients
  write_table(prediction_ladder$scores, file.path(output_dir, "tables/systems_prediction_module_scores.csv"))
  write_table(prediction_tbl_all_models, file.path(output_dir, "tables/systems_prediction_ladder_loo_predictions.csv"))
  write_table(prediction_perf, file.path(output_dir, "tables/systems_prediction_ladder_performance.csv"))
  write_table(prediction_coef_tbl, file.path(output_dir, "tables/systems_prediction_ladder_coefficients.csv"))
  write_table(bind_rows(prediction_ladder$scores, prediction_ladder_no_short$scores), file.path(output_dir, "tables/systems_prediction_module_scores_duration_sensitivity.csv"))
  write_table(prediction_tbl_all_models_duration_sensitivity, file.path(output_dir, "tables/systems_prediction_ladder_loo_predictions_duration_sensitivity.csv"))
  write_table(prediction_perf_duration_sensitivity, file.path(output_dir, "tables/systems_prediction_ladder_performance_duration_sensitivity.csv"))
  write_table(bind_rows(prediction_ladder$coefficients, prediction_ladder_no_short$coefficients), file.path(output_dir, "tables/systems_prediction_ladder_coefficients_duration_sensitivity.csv"))
  write_table(prediction_tbl, file.path(output_dir, "tables/systems_prospective_outcome_loo_predictions.csv"))
  write_table(prediction_perf_integrated, file.path(output_dir, "tables/systems_prospective_outcome_loo_performance.csv"))

  prediction_module_delta_tbl <- prediction_perf %>%
    transmute(
      ModelOrder,
      Model,
      AddedModule = case_when(
        Model == "movement-only" ~ "Magnitude",
        Model == "movement + proximity" ~ "Magnitude/social proximity",
        Model == "movement + entropy_acf1" ~ "Temporal persistence",
        Model == "movement + instability metrics" ~ "Temporal instability",
        Model == "movement + latent-state metrics" ~ "Latent-state organization",
        Model == "movement + social-network metrics" ~ "Social topology",
        Model == "integrated systems model" ~ "Trajectory + nonlinear integration",
        TRUE ~ Model
      ),
      cv_r2,
      delta_cv_r2_vs_previous,
      delta_cv_r2_vs_movement,
      pearson_r,
      rmse,
      permutation_p
    ) %>%
    mutate(delta_cv_r2_vs_previous = if_else(Model == "movement-only", cv_r2, delta_cv_r2_vs_previous))

  write_table(prediction_module_delta_tbl, file.path(output_dir, "tables/systems_prediction_module_delta_waterfall.csv"))

  p_delta <- prediction_module_delta_tbl %>%
    mutate(
      AddedModule = factor(AddedModule, levels = AddedModule),
      DeltaLabelVJust = if_else(delta_cv_r2_vs_previous >= 0, -0.2, 1.15)
    ) %>%
    ggplot(aes(AddedModule, delta_cv_r2_vs_previous, fill = delta_cv_r2_vs_previous)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_col(width = 0.68, colour = "white", linewidth = 0.18) +
    geom_text(aes(label = formatC(delta_cv_r2_vs_previous, format = "f", digits = 2), vjust = DeltaLabelVJust), size = 1.9) +
    scale_fill_gradient2(low = "#3d3b6e", mid = "white", high = "#e63947", midpoint = 0, na.value = "grey85") +
    labs(
      title = "Incremental predictive value by systems module",
      subtitle = paste0("Leave-one-animal-out prediction of ", outcome_to_plot, "; bars show delta CV-R2 versus previous ladder step"),
      x = NULL,
      y = "Delta CV-R2",
      fill = "Delta"
    ) +
    make_nature_theme(base_size = 6) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")

  save_plot_svg_pdf(p_delta, file.path(output_dir, "figures/publication_panels/Fig_systems_prediction_delta_waterfall"), width = 160, height = 86)

  p_duration_sensitivity <- prediction_perf_duration_sensitivity %>%
    filter(Model %in% c("movement-only", "movement + latent-state metrics", "movement + social-network metrics", "integrated systems model")) %>%
    mutate(Model = factor(Model, levels = unique(Model))) %>%
    ggplot(aes(Model, cv_r2, fill = DurationAnalysisSet)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_col(position = position_dodge(width = 0.70), width = 0.64, colour = "white", linewidth = 0.18) +
    labs(
      title = "Prediction robustness to short-duration epochs",
      subtitle = "Matched full-data versus excluding-short-duration model performance",
      x = NULL,
      y = "LOOCV R2",
      fill = "Analysis set"
    ) +
    scale_fill_manual(values = c(full = "#3d3b6e", excluding_short_duration = "#e63947"), drop = FALSE) +
    make_nature_theme(base_size = 6) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "top")

  save_plot_svg_pdf(p_duration_sensitivity, file.path(output_dir, "figures/publication_panels/Fig_systems_prediction_duration_sensitivity"), width = 150, height = 82)

  module_prediction_map <- tibble(
    Module = c("Magnitude", "Temporal instability", "Latent-state organization", "Social topology", "Trajectory geometry", "Nonlinear systems dynamics", "Predictive systems integration"),
    PredictionReadout = c(
      prediction_module_delta_tbl$cv_r2[prediction_module_delta_tbl$Model == "movement-only"][1],
      prediction_module_delta_tbl$delta_cv_r2_vs_previous[prediction_module_delta_tbl$Model == "movement + instability metrics"][1],
      prediction_module_delta_tbl$delta_cv_r2_vs_previous[prediction_module_delta_tbl$Model == "movement + latent-state metrics"][1],
      prediction_module_delta_tbl$delta_cv_r2_vs_previous[prediction_module_delta_tbl$Model == "movement + social-network metrics"][1],
      prediction_module_delta_tbl$delta_cv_r2_vs_previous[prediction_module_delta_tbl$Model == "integrated systems model"][1],
      prediction_module_delta_tbl$delta_cv_r2_vs_previous[prediction_module_delta_tbl$Model == "integrated systems model"][1],
      prediction_perf_integrated$cv_r2[1]
    )
  )

  module_scorecards <- module_scorecards_base %>%
    left_join(module_prediction_map, by = "Module") %>%
    mutate(
      PredictionInterpretation = case_when(
        !is.finite(PredictionReadout) ~ "not_available",
        PredictionReadout > 0.05 ~ "improves_prediction",
        PredictionReadout > 0 ~ "small_positive_increment",
        TRUE ~ "no_increment_or_overfit"
      )
    )

  write_table(module_scorecards, file.path(output_dir, "tables/systems_module_scorecards.csv"))

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

  p_ladder <- prediction_perf %>%
    mutate(Model = factor(Model, levels = Model)) %>%
    ggplot(aes(Model, cv_r2, fill = delta_cv_r2_vs_movement)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_col(width = 0.68, colour = "white", linewidth = 0.18) +
    geom_text(aes(label = paste0("r=", formatC(pearson_r, format = "f", digits = 2), "\nΔR2=", formatC(delta_cv_r2_vs_movement, format = "f", digits = 2))), size = 1.8, vjust = -0.15) +
    scale_fill_gradient2(low = "#3d3b6e", mid = "white", high = "#e63947", midpoint = 0, na.value = "grey85") +
    labs(
      title = paste0("Cross-validated prediction ladder for ", outcome_to_plot),
      subtitle = "Module scores test incremental value beyond early movement magnitude",
      x = NULL,
      y = "LOOCV R2 vs mean",
      fill = "ΔR2 vs movement",
      caption = "Permutation p-values and bootstrap CIs are exported in systems_prediction_ladder_performance.csv."
    ) +
    make_nature_theme(base_size = 6) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")
  save_plot_svg_pdf(p_ladder, file.path(output_dir, "figures/publication_panels/Fig_systems_prediction_ladder"), width = 165, height = 92)

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
      title = paste0("Integrated systems model predicts ", outcome_to_plot),
      subtitle = paste0("Module-level LOOCV: pooled r=", round(prediction_perf_integrated$pearson_r, 2), "; CV R2=", round(prediction_perf_integrated$cv_r2, 2)),
      x = paste0("Observed ", outcome_to_plot, " (lower = worse)"),
      y = paste0("Predicted ", outcome_to_plot)
    ) +
    make_nature_theme(base_size = 7)
  save_plot_svg_pdf(p_pred, file.path(output_dir, "figures/publication_panels/Fig_systems_prospective_crossvalidated_prediction"), width = 135, height = 75)

}

systems_endpoint_leakage_audit <- bind_rows(
  feature_dictionary %>%
    transmute(
      ItemType = "feature",
      Item = feature,
      SourceScript,
      FeatureUseClass,
      ClaimType,
      LeakageClass,
      TemporalRule = case_when(
        LeakageClass == "safe_for_prediction" ~ "pre-endpoint first-active-window behavioral feature",
        LeakageClass == "leaky_not_for_prediction" ~ "full-experiment or post-endpoint-relevant feature; do not use for prospective prediction",
        ClaimType == "associative" ~ "associative only unless temporal ordering is independently justified",
        TRUE ~ "descriptive only"
      ),
      ReviewerLabel = case_when(
        LeakageClass == "safe_for_prediction" ~ "safe_for_prediction",
        LeakageClass == "leaky_not_for_prediction" ~ "leaky_not_for_prediction",
        ClaimType == "associative" ~ "prospective_endpoint_association",
        TRUE ~ "descriptive_group_contrast"
      )
    ),
  tibble(
    ItemType = "group_contrast",
    Item = "RES/SUS/CON group contrasts",
    SourceScript = "12_systems_neuroscience_summary.R",
    FeatureUseClass = "primary_secondary_or_exploratory",
    ClaimType = "descriptive",
    LeakageClass = "descriptive_group_contrast",
    TemporalRule = "RES/SUS labels are derived from post-paradigm CombZ; group contrasts are descriptive phenotype contrasts.",
    ReviewerLabel = "descriptive_group_contrast"
  ),
  if (exists("prediction_perf")) {
    prediction_perf %>%
      transmute(
        ItemType = "prediction_model",
        Item = Model,
        SourceScript = "12_systems_neuroscience_summary.R",
        FeatureUseClass = "primary_or_secondary_module_scores",
        ClaimType = "predictive",
        LeakageClass = ModelLeakageClass,
        TemporalRule = ModelTemporalRule,
        ReviewerLabel = "safe_for_prediction"
      )
  } else tibble()
) %>%
  distinct()
write_table(systems_endpoint_leakage_audit, file.path(output_dir, "tables/systems_endpoint_leakage_audit.csv"))

# ------------------------------------------------
# OPTIONAL BEHAVIOR-PROTEOMICS BRIDGE
# ------------------------------------------------

proteomics_module_candidates <- c(
  "RNP_module", "RNA_processing_module", "Mito_module", "OXPHOS_module", "Translation_module",
  "Ribosomal_module", "Proteostasis_module", "Endolysosomal_module"
)
proteomics_feature_cols <- names(systems_features)[str_detect(names(systems_features), "^proteomics_module__")]
behavior_proteomics_axes <- unique(na.omit(c(
  early_movement_mean,
  make_feature_name("raw", "behavior", primary_bin_level, "Entropy", "acf1", prospective_prediction_context),
  pick_feature("Entropy.*rmssd|entropy.*rmssd", require_robust = TRUE),
  pick_feature("sleep_like_fraction_score|sleep_fragmentation_score|sleep_bout_duration_score", require_robust = TRUE),
  pick_feature("social_withdrawal_score|social_withdrawal", require_robust = TRUE),
  pick_feature("stress_adaptation_index|systems_resilience_score", require_robust = TRUE)
)))

systems_behavior_proteomics_bridge <- if (length(proteomics_feature_cols) > 0 && length(behavior_proteomics_axes) > 0) {
  map_dfr(c("pooled", as.character(sex_levels)), function(stratum) {
    dat <- systems_features
    if (stratum != "pooled") dat <- dat %>% filter(as.character(Sex) == stratum)
    map_dfr(behavior_proteomics_axes, function(bf) {
      map_dfr(proteomics_feature_cols, function(pf) {
        x <- safe_numeric(dat[[bf]])
        y <- safe_numeric(dat[[pf]])
        n_pair <- sum(is.finite(x) & is.finite(y))
        if (n_pair < 6) {
          return(tibble(Stratum = stratum, behavior_feature = bf, proteomics_module = pf, n = n_pair, spearman_rho = NA_real_, spearman_p = NA_real_, partial_rho_sex_adjusted = NA_real_, bootstrap_ci_low = NA_real_, bootstrap_ci_high = NA_real_))
        }
        ci <- bootstrap_cor_ci(x, y, "spearman", n_boot = 500, seed = 11)
        partial_rho <- NA_real_
        if (stratum == "pooled" && n_distinct(dat$Sex[is.finite(x) & is.finite(y)]) > 1) {
          ok <- is.finite(x) & is.finite(y) & !is.na(dat$Sex)
          rx <- residuals(lm(x[ok] ~ factor(dat$Sex[ok])))
          ry <- residuals(lm(y[ok] ~ factor(dat$Sex[ok])))
          partial_rho <- safe_cor(rx, ry, "spearman")
        }
        tibble(
          Stratum = stratum,
          behavior_feature = bf,
          proteomics_module = pf,
          n = n_pair,
          spearman_rho = safe_cor(x, y, "spearman"),
          spearman_p = safe_cor_p(x, y, "spearman"),
          partial_rho_sex_adjusted = partial_rho,
          bootstrap_ci_low = ci["low"],
          bootstrap_ci_high = ci["high"]
        )
      })
    })
  }) %>%
    left_join(feature_dictionary %>% select(behavior_feature = feature, BehaviorLabel = PrimaryFeatureLabel, FeatureUseClass), by = "behavior_feature") %>%
    mutate(
      proteomics_module = str_remove(proteomics_module, "^proteomics_module__"),
      spearman_fdr = p.adjust(spearman_p, method = "BH"),
      InterpretationLabel = case_when(
        n < 8 ~ "underpowered",
        FeatureUseClass %in% c("primary", "secondary") & !is.na(spearman_fdr) & spearman_fdr < 0.10 & abs(spearman_rho) >= 0.50 ~ "candidate mechanistic bridge",
        !is.na(spearman_fdr) & spearman_fdr < 0.10 ~ "associative molecular correlate",
        TRUE ~ "exploratory only"
      )
    )
} else {
  tibble(
    Stratum = "not_available",
    behavior_feature = NA_character_,
    proteomics_module = NA_character_,
    n = 0L,
    spearman_rho = NA_real_,
    spearman_p = NA_real_,
    partial_rho_sex_adjusted = NA_real_,
    bootstrap_ci_low = NA_real_,
    bootstrap_ci_high = NA_real_,
    BehaviorLabel = NA_character_,
    FeatureUseClass = NA_character_,
    spearman_fdr = NA_real_,
    InterpretationLabel = if (is.null(proteomics_module_file)) "proteomics_module_file_not_provided" else "proteomics_modules_not_detected"
  )
}
write_table(systems_behavior_proteomics_bridge, file.path(output_dir, "tables/systems_behavior_proteomics_bridge.csv"))

if (nrow(systems_behavior_proteomics_bridge) > 0 && any(is.finite(systems_behavior_proteomics_bridge$spearman_rho))) {
  p_prot_bridge <- systems_behavior_proteomics_bridge %>%
    filter(Stratum == "pooled", is.finite(spearman_rho)) %>%
    mutate(
      BehaviorLabel = factor(str_trunc(coalesce(BehaviorLabel, behavior_feature), 38), levels = rev(unique(str_trunc(coalesce(BehaviorLabel, behavior_feature), 38)))),
      proteomics_module = factor(str_replace_all(proteomics_module, "_", " "))
    ) %>%
    ggplot(aes(proteomics_module, BehaviorLabel, fill = spearman_rho)) +
    geom_tile(colour = "white", linewidth = 0.22) +
    geom_text(aes(label = sig_label(spearman_fdr)), size = 1.8) +
    scale_fill_gradient2(low = "#3d3b6e", mid = "white", high = "#e63947", midpoint = 0, na.value = "grey90") +
    labs(title = "Behavior-proteomics bridge", subtitle = "Associative molecular correlates of primary behavioral axes", x = NULL, y = NULL, fill = "rho") +
    make_nature_theme(base_size = 6) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "right")
  save_plot_svg_pdf(p_prot_bridge, file.path(output_dir, "figures/publication_panels/Fig_behavior_proteomics_bridge"), width = 160, height = 105)
}

# ------------------------------------------------
# LEAVE-ONE-CONTEXT AND BIN-SIZE ROBUSTNESS
# ------------------------------------------------

context_robustness_eval <- function(fc, context_col) {
  if (!context_col %in% names(systems_features)) return(tibble())
  dat <- systems_features %>%
    select(AnimalNum, Group, Sex, context = all_of(context_col), value = all_of(fc)) %>%
    filter(is.finite(value), !is.na(Group), !is.na(context), context != "")
  if (nrow(dat) < 8 || n_distinct(dat$context) < 2) return(tibble())
  full_eff <- group_contrasts %>%
    filter(feature == fc, contrast %in% c("SUS-CON", "RES-CON", "SUS-RES")) %>%
    arrange(desc(abs(hedges_g))) %>%
    slice_head(n = 1)
  map_dfr(unique(dat$context), function(level_out) {
    dd <- dat %>% filter(context != level_out)
    if (nrow(dd) < 6 || n_distinct(dd$Group) < 2) {
      return(tibble(feature = fc, ContextType = context_col, LevelLeftOut = as.character(level_out), n = nrow(dd), effect_after_leaveout = NA_real_, full_effect = ifelse(nrow(full_eff) > 0, full_eff$hedges_g[1], NA_real_), robustness_class = "not_estimable"))
    }
    pair <- if (nrow(full_eff) > 0) str_split(full_eff$contrast[1], "-", simplify = TRUE) else matrix(c("SUS", "CON"), nrow = 1)
    comp <- pair[1, 1]; ref <- pair[1, 2]
    eff <- hedges_g(dd$value[as.character(dd$Group) == ref], dd$value[as.character(dd$Group) == comp])
    full <- ifelse(nrow(full_eff) > 0, full_eff$hedges_g[1], NA_real_)
    tibble(
      feature = fc,
      ContextType = context_col,
      LevelLeftOut = as.character(level_out),
      n = nrow(dd),
      effect_after_leaveout = eff,
      full_effect = full,
      effect_percent_change = 100 * (eff - full) / pmax(abs(full), 1e-9),
      robustness_class = case_when(
        !is.finite(eff) | !is.finite(full) ~ "not_estimable",
        sign(eff) != sign(full) ~ paste0(str_to_lower(context_col), "_sensitive"),
        abs(effect_percent_change) > 50 ~ paste0(str_to_lower(context_col), "_sensitive"),
        nrow(dd) < 12 ~ "directionally stable but underpowered",
        TRUE ~ "robust"
      )
    )
  })
}

systems_leave_one_context_out_robustness <- map_dfr(unique(primary_feature_registry$feature), function(fc) {
  bind_rows(
    tibble(feature = fc, ContextType = "AnimalNum", LevelLeftOut = systems_features$AnimalNum, n = nrow(systems_features) - 1L) %>%
      mutate(effect_after_leaveout = NA_real_, full_effect = NA_real_, effect_percent_change = NA_real_, robustness_class = "leave-one-animal tracked by LOOCV where predictive"),
    context_robustness_eval(fc, "Batch"),
    context_robustness_eval(fc, "System")
  )
}) %>%
  left_join(primary_feature_registry %>% select(feature, PrimaryFeatureLabel, ClaimTier), by = "feature")
write_table(systems_leave_one_context_out_robustness, file.path(output_dir, "tables/systems_leave_one_context_out_robustness.csv"))

systems_bin_size_sensitivity_summary <- feature_dictionary %>%
  filter(feature %in% usable_features | FeatureUseClass %in% c("primary", "secondary")) %>%
  group_by(Metric, Statistic, Context, Source, Domain) %>%
  summarise(
    available_bin_levels = paste(sort(unique(BinLevel)), collapse = ";"),
    n_bin_levels = n_distinct(BinLevel),
    includes_primary_bin = primary_bin_level %in% BinLevel,
    any_primary_or_secondary = any(FeatureUseClass %in% c("primary", "secondary")),
    .groups = "drop"
  ) %>%
  mutate(
    robustness_class = case_when(
      n_bin_levels >= 3 & includes_primary_bin ~ "robust",
      n_bin_levels >= 2 ~ "directionally stable but underpowered",
      any_primary_or_secondary ~ "bin-size sensitive",
      TRUE ~ "not estimable"
    )
  )
write_table(systems_bin_size_sensitivity_summary, file.path(output_dir, "tables/systems_bin_size_sensitivity_summary.csv"))

primary_robustness_plot_tbl <- bind_rows(
  batch_system_cage_audit %>%
    transmute(PrimaryFeatureLabel, RobustnessAxis = "Batch/System covariates", robustness_class = interpretation_flag),
  systems_leave_one_context_out_robustness %>%
    filter(ContextType %in% c("Batch", "System")) %>%
    group_by(PrimaryFeatureLabel, RobustnessAxis = ContextType) %>%
    summarise(robustness_class = if_else(any(str_detect(robustness_class, "sensitive")), first(robustness_class[str_detect(robustness_class, "sensitive")]), "robust"), .groups = "drop"),
  systems_bin_size_sensitivity_summary %>%
    filter(any_primary_or_secondary %in% TRUE) %>%
    transmute(PrimaryFeatureLabel = paste(Metric, Statistic, Context, sep = " | "), RobustnessAxis = "Bin-size availability", robustness_class)
)
if (nrow(primary_robustness_plot_tbl) > 0) {
  p_primary_robust <- primary_robustness_plot_tbl %>%
    mutate(
      PrimaryFeatureLabel = factor(str_trunc(PrimaryFeatureLabel, 44), levels = rev(unique(str_trunc(PrimaryFeatureLabel, 44)))),
      robustness_class = factor(robustness_class)
    ) %>%
    ggplot(aes(RobustnessAxis, PrimaryFeatureLabel, fill = robustness_class)) +
    geom_tile(colour = "white", linewidth = 0.22) +
    labs(title = "Primary feature robustness summary", x = NULL, y = NULL, fill = "Classification") +
    make_nature_theme(base_size = 6) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")
  save_plot_svg_pdf(p_primary_robust, file.path(output_dir, "figures/qc/Fig_primary_feature_robustness"), width = 160, height = 105)
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
        subtitle = paste0(
          "Module ladder; pooled r=",
          round(if (exists("prediction_perf_integrated")) prediction_perf_integrated$loo_pearson_r else prediction_perf$loo_pearson_r[1], 2),
          ", CV R2=",
          round(if (exists("prediction_perf_integrated")) prediction_perf_integrated$cv_r2_vs_mean else prediction_perf$cv_r2_vs_mean[1], 2)
        ),
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

  module_effect_tbl <- heat_tbl %>%
    mutate(Module = if_else(is.na(Module), "Other interpretable feature", Module)) %>%
    group_by(Sex, contrast, Module) %>%
    summarise(
      n_features = n(),
      median_g = median(hedges_g, na.rm = TRUE),
      max_abs_g = hedges_g[which.max(abs(hedges_g))][1],
      min_q = min(p_fdr, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      module_effect = if_else(is.finite(median_g), median_g, max_abs_g),
      sig = sig_label(min_q),
      Module = factor(Module, levels = rev(c(
        "Magnitude", "Temporal instability", "Latent-state organization", "Social topology",
        "Trajectory geometry", "Nonlinear systems dynamics", "Predictive systems integration",
        "Other interpretable feature"
      )))
    ) %>%
    filter(is.finite(module_effect))

  write_table(module_effect_tbl, file.path(output_dir, "tables/systems_dashboard_module_effect_summary.csv"))

  p_heat_small <- module_effect_tbl %>%
    ggplot(aes(contrast, Module, fill = module_effect)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(aes(label = sig), size = 1.8) +
    facet_grid(Sex ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(low = "#3d3b6e", mid = "white", high = "#e63947", midpoint = 0, na.value = "grey90") +
    labs(
      title = "F. Multiscale phenotype modules",
      subtitle = "Median Hedges g across feature modules; detailed map exported separately",
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
        subtitle = paste0("PCA on ", length(duration_robust_features), " duration-robust animal-level features; points are animals")
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

  p_dashboard_trajectory <- if (exists("p_first_active_movement") && inherits(p_first_active_movement, "ggplot")) {
    p_first_active_movement +
      labs(title = "A. First active movement trajectory", subtitle = paste0("Raw binned behavior; first 12 h after ", first_cage_change)) +
      theme(legend.position = "none")
  } else {
    ggplot() + annotate("text", x = 0, y = 0, label = "First-active movement trajectory unavailable", size = 3) + theme_void()
  }

  primary_heat_tbl <- group_contrasts %>%
    filter(feature %in% c(primary_claim_features, secondary_claim_features), is.finite(hedges_g)) %>%
    left_join(primary_feature_registry %>% select(feature, PrimaryFeatureLabel, ClaimTier), by = "feature") %>%
    mutate(
      PrimaryFeatureLabel = factor(str_trunc(coalesce(PrimaryFeatureLabel, feature), 42), levels = rev(unique(str_trunc(coalesce(PrimaryFeatureLabel, feature), 42)))),
      contrast = factor(contrast, levels = c("RES-CON", "SUS-CON", "SUS-RES"))
    )
  p_dashboard_heat <- if (nrow(primary_heat_tbl) > 0) {
    primary_heat_tbl %>%
      ggplot(aes(contrast, PrimaryFeatureLabel, fill = hedges_g)) +
      geom_tile(colour = "white", linewidth = 0.22) +
      geom_text(aes(label = sig_label(p_fdr)), size = 1.75) +
      facet_grid(Sex ~ ., scales = "free_y", space = "free_y") +
      scale_fill_gradient2(low = "#3d3b6e", mid = "white", high = "#e63947", midpoint = 0, na.value = "grey90") +
      labs(title = "B. Primary feature effects", subtitle = "Hedges g; symbols denote BH FDR", x = NULL, y = NULL, fill = "g") +
      make_nature_theme(base_size = 5.5) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")
  } else {
    p_heat_small + labs(title = "B. Module effects", subtitle = "Primary-feature heatmap unavailable")
  }

  p_dashboard_sleep <- if (exists("p_sleep_like") && inherits(p_sleep_like, "ggplot")) {
    p_sleep_like + labs(title = "C. Sleep-like inactivity", subtitle = "Quiescence-like inactivity, not EEG sleep") + theme(legend.position = "none")
  } else {
    ggplot() + annotate("text", x = 0, y = 0, label = "Sleep-like inactivity panel unavailable", size = 3) + theme_void()
  }

  p_dashboard_interaction <- if (exists("p_interaction") && inherits(p_interaction, "ggplot")) {
    p_interaction + labs(title = "D. Group x sex interaction screen", subtitle = "Sex-specific claims require interaction support")
  } else {
    ggplot() + annotate("text", x = 0, y = 0, label = "Group x sex interaction model unavailable", size = 3) + theme_void()
  }

  p_dashboard_prediction <- p_pred_dash +
    labs(title = paste0("E. Prospective ", primary_outcome, " association"), subtitle = "Pre-endpoint features only")

  p_dashboard_proteomics <- if (exists("p_prot_bridge") && inherits(p_prot_bridge, "ggplot")) {
    p_prot_bridge + labs(title = "F. Behavior-proteomics bridge", subtitle = "Associative molecular correlate")
  } else {
    ggplot() + annotate("text", x = 0, y = 0, label = "F. Proteomics bridge unavailable\n(proteomics_module_file not provided)", size = 2.8) + theme_void()
  }

  dashboard <- ((p_dashboard_trajectory | p_dashboard_heat) / (p_dashboard_sleep | p_dashboard_interaction) / (p_dashboard_prediction | p_dashboard_proteomics)) +
    patchwork::plot_layout(heights = c(1.1, 0.95, 0.95)) +
    patchwork::plot_annotation(
      title = "Behavioral dynamics of stress susceptibility",
      subtitle = paste0("Primary claim dashboard: early trajectories, primary features, quiescence-like inactivity, interaction tests and pre-endpoint ", primary_outcome, " prediction"),
      caption = "RES/SUS are CombZ-derived labels; group contrasts are descriptive. Prediction uses pre-endpoint first-active-window features only. PCA/UMAP and high-dimensional nonlinear outputs are supplementary architecture views."
    )

  save_plot_svg_pdf(dashboard, file.path(output_dir, "figures/Fig_integrated_systems_dashboard"), width = 230, height = 245)
}

# ------------------------------------------------
# SIS-SPECIFIC BIOLOGICAL DASHBOARD
# ------------------------------------------------

# This layer reframes the integrated outputs around the biology of adolescent
# social instability stress (SIS): first perturbation response, repeated
# regrouping adaptation, phase-specific organization, social spatial structure,
# rest-like fragmentation and state architecture. It intentionally uses
# domain-level scores to avoid letting locomotion dominate the dashboard.

standardize_within_context <- function(dat, value_col, group_cols = c("Sex", "PhaseClass", "CageChangeIndex")) {
  dat %>%
    group_by(across(any_of(group_cols))) %>%
    mutate("{value_col}_z" := safe_scale(safe_numeric(.data[[value_col]]))) %>%
    ungroup()
}

score_mean <- function(dat, cols) {
  cols <- intersect(cols, names(dat))
  if (length(cols) == 0) return(rep(NA_real_, nrow(dat)))
  out <- rowMeans(as.matrix(dat[, cols, drop = FALSE]), na.rm = TRUE)
  out[!is.finite(out)] <- NA_real_
  out
}

safe_metric_summary <- function(x) {
  x <- x[is.finite(x)]
  tibble(
    mean = if (length(x) > 0) mean(x) else NA_real_,
    sd = if (length(x) > 1) sd(x) else NA_real_,
    rmssd = if (length(x) > 2) sqrt(mean(diff(x)^2, na.rm = TRUE)) else NA_real_,
    acf1 = if (length(x) > 3) safe_cor(x[-length(x)], x[-1], "pearson") else NA_real_
  )
}

extract_lmm_stats <- function(dat, response, model_label, include_cage_change = TRUE) {
  dat <- dat %>%
    filter(is.finite(.data[[response]]), !is.na(Group), !is.na(Sex), !is.na(PhaseClass)) %>%
    mutate(
      Group = factor(as.character(Group), levels = group_levels),
      Sex = factor(as.character(Sex), levels = sex_levels),
      PhaseClass = factor(as.character(PhaseClass), levels = c("Active", "Inactive")),
      CageChangeIndex = safe_numeric(CageChangeIndex)
    )
  if (nrow(dat) < 8 || n_distinct(dat$AnimalNum) < 4 || n_distinct(dat$Group) < 2) {
    return(tibble(model = model_label, domain = response, term = NA_character_, estimate = NA_real_, std.error = NA_real_, statistic = NA_real_, p.value = NA_real_, model_engine = "not_estimable", n = nrow(dat)))
  }

  fixed <- if (include_cage_change && n_distinct(dat$CageChangeIndex[is.finite(dat$CageChangeIndex)]) >= 2) {
    paste0(response, " ~ Group * Sex * PhaseClass + CageChangeIndex")
  } else {
    paste0(response, " ~ Group * Sex * PhaseClass")
  }
  random_form <- as.formula(paste0(fixed, " + (1 | AnimalNum)"))
  fixed_form <- as.formula(fixed)

  if (requireNamespace("lmerTest", quietly = TRUE)) {
    fit <- try(lmerTest::lmer(random_form, data = dat), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      cf <- coef(summary(fit)) %>%
        as.data.frame() %>%
        rownames_to_column("term") %>%
        as_tibble()
      return(cf %>%
        transmute(
          model = model_label,
          domain = response,
          term,
          estimate = Estimate,
          std.error = `Std. Error`,
          statistic = `t value`,
          p.value = `Pr(>|t|)`,
          model_engine = "lmerTest::lmer",
          n = nrow(dat)
        ))
    }
  }

  fit <- try(lm(fixed_form, data = dat), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(tibble(model = model_label, domain = response, term = NA_character_, estimate = NA_real_, std.error = NA_real_, statistic = NA_real_, p.value = NA_real_, model_engine = "not_estimable", n = nrow(dat)))
  }
  coef(summary(fit)) %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    as_tibble() %>%
    transmute(
      model = model_label,
      domain = response,
      term,
      estimate = Estimate,
      std.error = `Std. Error`,
      statistic = `t value`,
      p.value = `Pr(>|t|)`,
      model_engine = "lm_fallback",
      n = nrow(dat)
    )
}

sis_domain_interpretation <- tibble(
  Domain = c(
    "Early adaptation / prediction",
    "Active-phase adaptation/exploration",
    "Repeated adaptation / recovery dynamics",
    "Inactive-phase rest/circadian regulation",
    "Behavioral flexibility / predictability",
    "Social spatial organization",
    "Behavioral state architecture",
    "Behavioral volatility / fragmentation",
    "Psychomotor activation"
  ),
  BiologicalInterpretation = c(
    "Initial active-phase response to the first social regrouping; candidate vulnerability/resilience signal.",
    "Active-phase exploration, psychomotor engagement, social investigation and novelty adaptation.",
    "Habituation, sensitization or recovery across repeated post-regrouping active phases.",
    "Inactive-phase recovery, rest-like consolidation and circadian-context stability.",
    "Spatial diversity and predictability; high entropy can mean flexible exploration or disorganization depending on context.",
    "Social proximity and stability of spatial co-organization; not direct sociability.",
    "Occupancy, dwell time and switching among inferred behavioral states.",
    "Short-timescale volatility, fragmented behavior or unstable engagement.",
    "Locomotor output treated as one domain of adaptation, not the full phenotype."
  ),
  PrimaryRawContributors = c(
    "Movement, entropy, proximity, RMSSD and ACF1 during first 12 h active phase after first regrouping.",
    "Active-phase movement, entropy, proximity and low persistence during post-regrouping windows.",
    "Cage-change trajectories, early-window slopes, volatility decay and distance-to-control summaries.",
    "Inactive-phase inactivity fraction, bout duration, fragmentation, movement stability and active/inactive contrasts.",
    "Entropy mean, entropy RMSSD and entropy ACF1 in phase-specific epochs.",
    "Mean proximity, proximity RMSSD/ACF1 and social-network features when available.",
    "HMM/state-space occupancy, transition probability, dwell time and switch-rate outputs when available.",
    "RMSSD, burstiness, inactivity fragmentation and state-switching metrics.",
    "Movement mean and movement trajectory summaries."
  ),
  Directness = c("indirect/prospective association", "direct RFID behavioral summary", "direct longitudinal behavioral summary", "indirect rest-like RFID summary", "indirect spatial-organization summary", "indirect social-spatial proxy", "model-derived latent construct", "direct temporal-structure summary", "direct locomotor summary"),
  Caveat = c(
    "Associative unless externally validated; RES/SUS labels are CombZ-derived.",
    "Active phase should not be pooled with inactive phase without biological justification.",
    "Trajectory differences do not prove recovery mechanisms.",
    "Sleep-like/rest-like inactivity is not EEG-confirmed sleep.",
    "Entropy direction is context-dependent.",
    "Proximity is not equivalent to sociability or preference.",
    "HMM states are data-derived and require semantic caution.",
    "Volatility can reflect adaptive exploration or maladaptive fragmentation.",
    "Avoid reducing SIS biology to hypoactivity or hyperactivity."
  )
)
write_table(sis_domain_interpretation, file.path(output_dir, "tables/systems_sis_domain_interpretation_guide.csv"))

sleep_features_epoch <- read_any_table(file.path(project_root, "analysis_ready/16_sleep_like_inactivity_metrics", primary_bin_level, "tables/sleep_like_inactivity_features.csv"))
if (!is.null(sleep_features_epoch) && nrow(sleep_features_epoch) > 0) {
  sleep_features_epoch <- sleep_features_epoch %>%
    standardize_id_columns() %>%
    mutate(
      CageChange = as.character(CageChange),
      PhaseClass = if ("PhaseClass" %in% names(.)) as.character(PhaseClass) else as.character(Phase),
      CageChangeIndex = if ("CageChangeIndex" %in% names(.)) safe_numeric(CageChangeIndex) else parse_cage_change_index(CageChange)
    ) %>%
    select(
      AnimalNum, CageChange, CageChangeIndex, PhaseClass,
      any_of(c(
        "inactivity_fraction", "zero_like_inactivity_fraction", "mean_inactivity_bout_min",
        "median_inactivity_bout_min", "max_inactivity_bout_min", "inactivity_fragmentation",
        "active_inactive_transition_rate", "prolonged_inactivity_episodes_per_hour"
      ))
    )
} else {
  sleep_features_epoch <- tibble()
}

sis_epoch_raw <- base %>%
  filter(PhaseClass %in% c("Active", "Inactive"), is.finite(CageChangeIndex)) %>%
  group_by(AnimalNum, Group, Sex, CageChange, CageChangeIndex, PhaseClass) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  summarise(
    Movement_mean = mean(Movement, na.rm = TRUE),
    Movement_rmssd = if (sum(is.finite(Movement)) >= 3) sqrt(mean(diff(Movement[is.finite(Movement)])^2, na.rm = TRUE)) else NA_real_,
    Movement_acf1 = if (sum(is.finite(Movement)) >= 4) safe_cor(Movement[is.finite(Movement)][-sum(is.finite(Movement))], Movement[is.finite(Movement)][-1], "pearson") else NA_real_,
    Entropy_mean = mean(Entropy, na.rm = TRUE),
    Entropy_rmssd = if (sum(is.finite(Entropy)) >= 3) sqrt(mean(diff(Entropy[is.finite(Entropy)])^2, na.rm = TRUE)) else NA_real_,
    Entropy_acf1 = if (sum(is.finite(Entropy)) >= 4) safe_cor(Entropy[is.finite(Entropy)][-sum(is.finite(Entropy))], Entropy[is.finite(Entropy)][-1], "pearson") else NA_real_,
    Proximity_mean = mean(Proximity, na.rm = TRUE),
    Proximity_rmssd = if (sum(is.finite(Proximity)) >= 3) sqrt(mean(diff(Proximity[is.finite(Proximity)])^2, na.rm = TRUE)) else NA_real_,
    Proximity_acf1 = if (sum(is.finite(Proximity)) >= 4) safe_cor(Proximity[is.finite(Proximity)][-sum(is.finite(Proximity))], Proximity[is.finite(Proximity)][-1], "pearson") else NA_real_,
    n_bins = n(),
    .groups = "drop"
  ) %>%
  left_join(sleep_features_epoch, by = c("AnimalNum", "CageChange", "CageChangeIndex", "PhaseClass"))

sis_z_cols <- c(
  "Movement_mean", "Movement_rmssd", "Movement_acf1",
  "Entropy_mean", "Entropy_rmssd", "Entropy_acf1",
  "Proximity_mean", "Proximity_rmssd", "Proximity_acf1",
  "inactivity_fraction", "mean_inactivity_bout_min", "inactivity_fragmentation",
  "active_inactive_transition_rate", "prolonged_inactivity_episodes_per_hour"
)
sis_epoch_score_base <- reduce(
  intersect(sis_z_cols, names(sis_epoch_raw)),
  standardize_within_context,
  .init = sis_epoch_raw
)
for (needed_col in paste0(sis_z_cols, "_z")) {
  if (!needed_col %in% names(sis_epoch_score_base)) sis_epoch_score_base[[needed_col]] <- NA_real_
}

sis_epoch_scores <- sis_epoch_score_base %>%
  mutate(
    `Psychomotor activation` = Movement_mean_z,
    `Behavioral flexibility / predictability` = score_mean(pick(everything()), c("Entropy_mean_z", "Entropy_rmssd_z")) - coalesce(Entropy_acf1_z, 0),
    `Social spatial organization` = score_mean(pick(everything()), c("Proximity_mean_z", "Proximity_acf1_z")) - coalesce(Proximity_rmssd_z, 0),
    `Behavioral volatility / fragmentation` = score_mean(pick(everything()), c("Movement_rmssd_z", "Entropy_rmssd_z", "Proximity_rmssd_z", "inactivity_fragmentation_z", "active_inactive_transition_rate_z")),
    `Inactive-phase rest/circadian regulation` = if_else(
      PhaseClass == "Inactive",
      score_mean(pick(everything()), c("inactivity_fraction_z", "mean_inactivity_bout_min_z", "Movement_acf1_z")) -
        score_mean(pick(everything()), c("Movement_rmssd_z", "inactivity_fragmentation_z", "active_inactive_transition_rate_z")),
      NA_real_
    ),
    `Active-phase adaptation/exploration` = if_else(
      PhaseClass == "Active",
      score_mean(pick(everything()), c("Movement_mean_z", "Entropy_mean_z", "Proximity_mean_z")) -
        score_mean(pick(everything()), c("Movement_acf1_z", "Entropy_acf1_z")),
      NA_real_
    )
  ) %>%
  mutate(
    `Repeated adaptation / recovery dynamics` = if_else(
      PhaseClass == "Active",
      score_mean(pick(everything()), c("Active-phase adaptation/exploration", "Behavioral flexibility / predictability", "Social spatial organization")) -
        coalesce(`Behavioral volatility / fragmentation`, 0),
      NA_real_
    ),
    `Early adaptation / prediction` = if_else(
      PhaseClass == "Active" & CageChangeIndex == min(CageChangeIndex, na.rm = TRUE),
      `Active-phase adaptation/exploration`,
      NA_real_
    )
  )

hmm_epoch_scores <- if (!is.null(hmm_occupancy) && nrow(hmm_occupancy) > 0) {
  state_summary_for_arch <- state_label_tbl
  hmm_occupancy %>%
    standardize_id_columns() %>%
    mutate(
      PhaseClass = case_when(
        str_detect(str_to_lower(as.character(Phase)), "inactive|light|day") ~ "Inactive",
        str_detect(str_to_lower(as.character(Phase)), "active|dark|night") ~ "Active",
        TRUE ~ as.character(Phase)
      ),
      CageChange = as.character(CageChange),
      CageChangeIndex = parse_cage_change_index(CageChange),
      State = as.character(State),
      frac_time = safe_numeric(frac_time)
    ) %>%
    left_join(state_summary_for_arch %>% select(State, SemanticState), by = "State") %>%
    group_by(AnimalNum, CageChange, CageChangeIndex, PhaseClass) %>%
    summarise(
      state_occupancy_entropy = feature_entropy(frac_time / sum(frac_time, na.rm = TRUE)),
      inactive_state_fraction = sum(frac_time[SemanticState == "inactive/low-exploration"], na.rm = TRUE),
      social_state_fraction = sum(frac_time[SemanticState == "social"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    standardize_within_context("state_occupancy_entropy") %>%
    standardize_within_context("inactive_state_fraction") %>%
    standardize_within_context("social_state_fraction") %>%
    mutate(`Behavioral state architecture` = score_mean(pick(everything()), c("state_occupancy_entropy_z", "social_state_fraction_z")) - coalesce(inactive_state_fraction_z, 0)) %>%
    select(AnimalNum, CageChange, CageChangeIndex, PhaseClass, `Behavioral state architecture`)
} else {
  tibble()
}

if (nrow(hmm_epoch_scores) > 0) {
  sis_epoch_scores <- sis_epoch_scores %>%
    left_join(hmm_epoch_scores, by = c("AnimalNum", "CageChange", "CageChangeIndex", "PhaseClass"))
} else {
  sis_epoch_scores <- sis_epoch_scores %>%
    mutate(`Behavioral state architecture` = score_mean(pick(everything()), c("Behavioral flexibility / predictability", "Behavioral volatility / fragmentation")))
}

sis_domain_cols <- c(
  "Early adaptation / prediction",
  "Active-phase adaptation/exploration",
  "Repeated adaptation / recovery dynamics",
  "Inactive-phase rest/circadian regulation",
  "Behavioral flexibility / predictability",
  "Social spatial organization",
  "Behavioral state architecture",
  "Behavioral volatility / fragmentation",
  "Psychomotor activation"
)

sis_domain_scores <- sis_epoch_scores %>%
  select(AnimalNum, Group, Sex, CageChange, CageChangeIndex, PhaseClass, any_of(sis_domain_cols)) %>%
  pivot_longer(any_of(sis_domain_cols), names_to = "Domain", values_to = "DomainScore") %>%
  filter(is.finite(DomainScore)) %>%
  left_join(sis_domain_interpretation, by = "Domain")

write_table(sis_epoch_raw, file.path(output_dir, "tables/systems_sis_raw_phase_epoch_features.csv"))
write_table(sis_domain_scores, file.path(output_dir, "tables/systems_sis_domain_scores.csv"))

sis_active_inactive_contrasts <- sis_domain_scores %>%
  filter(PhaseClass %in% c("Active", "Inactive")) %>%
  group_by(AnimalNum, Group, Sex, CageChange, CageChangeIndex, PhaseClass, Domain, BiologicalInterpretation) %>%
  summarise(DomainScore = mean(DomainScore, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = PhaseClass, values_from = DomainScore)
if (!"Active" %in% names(sis_active_inactive_contrasts)) sis_active_inactive_contrasts$Active <- NA_real_
if (!"Inactive" %in% names(sis_active_inactive_contrasts)) sis_active_inactive_contrasts$Inactive <- NA_real_
sis_active_inactive_contrasts <- sis_active_inactive_contrasts %>%
  mutate(
    active_minus_inactive = Active - Inactive,
    active_inactive_abs_contrast = abs(active_minus_inactive)
  )
write_table(sis_active_inactive_contrasts, file.path(output_dir, "tables/systems_sis_active_inactive_domain_contrasts.csv"))

sis_model_stats <- map_dfr(sis_domain_cols, function(dom) {
  model_dat <- sis_domain_scores %>%
    filter(Domain == dom) %>%
    transmute(AnimalNum, Group, Sex, PhaseClass, CageChangeIndex, domain_score = DomainScore)
  extract_lmm_stats(model_dat, "domain_score", dom, include_cage_change = TRUE)
}) %>%
  group_by(model) %>%
  mutate(p_fdr = p.adjust(p.value, method = "BH")) %>%
  ungroup()
write_table(sis_model_stats, file.path(output_dir, "stats_tables/systems_sis_domain_mixed_model_stats.csv"))

sis_feature_redundancy_inputs <- sis_epoch_scores %>%
  select(AnimalNum, Group, Sex, CageChangeIndex, PhaseClass, any_of(sis_z_cols), any_of(paste0(sis_z_cols, "_z"))) %>%
  select(where(~ !all(is.na(.x))))
redundancy_numeric <- names(sis_feature_redundancy_inputs)[map_lgl(sis_feature_redundancy_inputs, is.numeric)]
redundancy_numeric <- setdiff(redundancy_numeric, "CageChangeIndex")
sis_feature_redundancy <- if (length(redundancy_numeric) >= 2) {
  red_mat <- sis_feature_redundancy_inputs %>%
    select(all_of(redundancy_numeric)) %>%
    mutate(across(everything(), ~ replace_na(safe_numeric(.x), median(safe_numeric(.x), na.rm = TRUE)))) %>%
    as.matrix()
  bad_redundancy_cols <- !is.finite(colSums(red_mat))
  if (any(bad_redundancy_cols)) red_mat[, bad_redundancy_cols] <- 0
  cmat <- suppressWarnings(cor(red_mat, method = "spearman", use = "pairwise.complete.obs"))
  idx <- which(upper.tri(cmat), arr.ind = TRUE)
  movement_ref <- grep("Movement_mean", colnames(cmat), value = TRUE)[1]
  tibble(
    Feature1 = colnames(cmat)[idx[, 1]],
    Feature2 = colnames(cmat)[idx[, 2]],
    spearman_rho = cmat[idx]
  ) %>%
    filter(is.finite(spearman_rho)) %>%
    mutate(
      redundancy_level = case_when(
        abs(spearman_rho) >= 0.95 ~ "near_duplicate",
        abs(spearman_rho) >= 0.85 ~ "highly_redundant",
        abs(spearman_rho) >= 0.70 ~ "moderately_redundant",
        TRUE ~ "not_highly_redundant"
      ),
      locomotion_driven = if (!is.na(movement_ref)) {
        abs(cmat[cbind(match(Feature1, colnames(cmat)), match(movement_ref, colnames(cmat)))]) >= 0.70 |
          abs(cmat[cbind(match(Feature2, colnames(cmat)), match(movement_ref, colnames(cmat)))]) >= 0.70
      } else {
        FALSE
      }
    ) %>%
    arrange(desc(abs(spearman_rho)))
} else {
  tibble()
}
write_table(sis_feature_redundancy, file.path(output_dir, "tables/systems_sis_feature_redundancy_diagnostics.csv"))

sis_locomotion_dominance <- sis_domain_scores %>%
  select(AnimalNum, Sex, Group, CageChangeIndex, PhaseClass, Domain, DomainScore) %>%
  left_join(sis_epoch_scores %>% select(AnimalNum, CageChangeIndex, PhaseClass, Movement_mean_z), by = c("AnimalNum", "CageChangeIndex", "PhaseClass")) %>%
  group_by(Domain, PhaseClass) %>%
  summarise(
    n = sum(is.finite(DomainScore) & is.finite(Movement_mean_z)),
    spearman_rho_with_movement = safe_cor(DomainScore, Movement_mean_z, "spearman"),
    spearman_p = safe_cor_p(DomainScore, Movement_mean_z, "spearman"),
    .groups = "drop"
  ) %>%
  group_by(PhaseClass) %>%
  mutate(
    spearman_fdr = p.adjust(spearman_p, method = "BH"),
    locomotion_dominance_flag = abs(spearman_rho_with_movement) >= 0.70
  ) %>%
  ungroup()
write_table(sis_locomotion_dominance, file.path(output_dir, "tables/systems_sis_locomotion_dominance_audit.csv"))

sis_first_active_raw <- first_active %>%
  group_by(AnimalNum, Group, Sex) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  summarise(
    Movement_mean = mean(Movement, na.rm = TRUE),
    Movement_rmssd = if (sum(is.finite(Movement)) >= 3) sqrt(mean(diff(Movement[is.finite(Movement)])^2, na.rm = TRUE)) else NA_real_,
    Movement_acf1 = if (sum(is.finite(Movement)) >= 4) safe_cor(Movement[is.finite(Movement)][-sum(is.finite(Movement))], Movement[is.finite(Movement)][-1], "pearson") else NA_real_,
    Entropy_mean = mean(Entropy, na.rm = TRUE),
    Entropy_rmssd = if (sum(is.finite(Entropy)) >= 3) sqrt(mean(diff(Entropy[is.finite(Entropy)])^2, na.rm = TRUE)) else NA_real_,
    Entropy_acf1 = if (sum(is.finite(Entropy)) >= 4) safe_cor(Entropy[is.finite(Entropy)][-sum(is.finite(Entropy))], Entropy[is.finite(Entropy)][-1], "pearson") else NA_real_,
    Proximity_mean = mean(Proximity, na.rm = TRUE),
    Proximity_rmssd = if (sum(is.finite(Proximity)) >= 3) sqrt(mean(diff(Proximity[is.finite(Proximity)])^2, na.rm = TRUE)) else NA_real_,
    Proximity_acf1 = if (sum(is.finite(Proximity)) >= 4) safe_cor(Proximity[is.finite(Proximity)][-sum(is.finite(Proximity))], Proximity[is.finite(Proximity)][-1], "pearson") else NA_real_,
    n_bins = n(),
    .groups = "drop"
  )
sis_first_active_score_base <- reduce(
  intersect(c(
    "Movement_mean", "Movement_rmssd", "Movement_acf1",
    "Entropy_mean", "Entropy_rmssd", "Entropy_acf1",
    "Proximity_mean", "Proximity_rmssd", "Proximity_acf1"
  ), names(sis_first_active_raw)),
  function(dat, value_col) standardize_within_context(dat, value_col, group_cols = "Sex"),
  .init = sis_first_active_raw
)
for (needed_col in paste0(c(
  "Movement_mean", "Movement_rmssd", "Movement_acf1",
  "Entropy_mean", "Entropy_rmssd", "Entropy_acf1",
  "Proximity_mean", "Proximity_rmssd", "Proximity_acf1"
), "_z")) {
  if (!needed_col %in% names(sis_first_active_score_base)) sis_first_active_score_base[[needed_col]] <- NA_real_
}
first_active_domain_scores <- sis_first_active_score_base %>%
  mutate(
    `Psychomotor activation` = Movement_mean_z,
    `Behavioral flexibility / predictability` = score_mean(pick(everything()), c("Entropy_mean_z", "Entropy_rmssd_z")) - coalesce(Entropy_acf1_z, 0),
    `Social spatial organization` = score_mean(pick(everything()), c("Proximity_mean_z", "Proximity_acf1_z")) - coalesce(Proximity_rmssd_z, 0),
    `Behavioral volatility / fragmentation` = score_mean(pick(everything()), c("Movement_rmssd_z", "Entropy_rmssd_z", "Proximity_rmssd_z")),
    `Active-phase adaptation/exploration` = score_mean(pick(everything()), c("Movement_mean_z", "Entropy_mean_z", "Proximity_mean_z")) -
      score_mean(pick(everything()), c("Movement_acf1_z", "Entropy_acf1_z")),
    `Early adaptation / prediction` = `Active-phase adaptation/exploration`
  ) %>%
  select(AnimalNum, Group, Sex, any_of(c("Psychomotor activation", "Behavioral flexibility / predictability", "Social spatial organization", "Behavioral volatility / fragmentation", "Active-phase adaptation/exploration", "Early adaptation / prediction"))) %>%
  pivot_longer(-c(AnimalNum, Group, Sex), names_to = "Domain", values_to = "EarlyDomainScore") %>%
  filter(is.finite(EarlyDomainScore)) %>%
  left_join(sis_domain_interpretation, by = "Domain") %>%
  left_join(systems_features %>% select(AnimalNum, any_of(endpoint_cols)), by = "AnimalNum")

first_active_prediction_table <- map_dfr(intersect(endpoint_cols, names(first_active_domain_scores)), function(outcome) {
  first_active_domain_scores %>%
    filter(is.finite(.data[[outcome]]), is.finite(EarlyDomainScore)) %>%
    group_by(Sex, Domain, BiologicalInterpretation, Caveat) %>%
    summarise(
      outcome = outcome,
      n = n_distinct(AnimalNum),
      spearman_rho = safe_cor(EarlyDomainScore, .data[[outcome]], "spearman"),
      spearman_p = safe_cor_p(EarlyDomainScore, .data[[outcome]], "spearman"),
      pearson_r = safe_cor(EarlyDomainScore, .data[[outcome]], "pearson"),
      .groups = "drop"
    )
}) %>%
  group_by(outcome, Sex) %>%
  mutate(
    spearman_fdr = p.adjust(spearman_p, method = "BH"),
    Interpretation = "Candidate early vulnerability/resilience signal from first 12 h ACTIVE phase after first regrouping; associative unless externally validated."
  ) %>%
  ungroup()
write_table(first_active_domain_scores, file.path(output_dir, "tables/systems_sis_first_active_12h_domain_scores.csv"))
write_table(first_active_prediction_table, file.path(output_dir, "stats_tables/systems_sis_first_active_12h_prediction_table.csv"))

sis_group_summary <- sis_domain_scores %>%
  group_by(Domain, PhaseClass, CageChangeIndex, Sex, Group) %>%
  summarise(
    n = n_distinct(AnimalNum),
    mean = mean(DomainScore, na.rm = TRUE),
    ci_low = mean_ci(DomainScore)["low"],
    ci_high = mean_ci(DomainScore)["high"],
    sem = sem(DomainScore),
    .groups = "drop"
  )
write_table(sis_group_summary, file.path(output_dir, "stats_tables/systems_sis_domain_group_summary.csv"))

domain_effect_summary <- sis_domain_scores %>%
  group_by(Domain, PhaseClass, Sex) %>%
  group_modify(~{
    d <- .x
    map_dfr(contrast_pairs, function(pair) {
      ref <- pair[1]; comp <- pair[2]
      x <- d$DomainScore[as.character(d$Group) == ref]
      y <- d$DomainScore[as.character(d$Group) == comp]
      tibble(
        contrast = paste0(comp, "-", ref),
        n_ref = sum(is.finite(x)),
        n_comp = sum(is.finite(y)),
        hedges_g = hedges_g(x, y),
        mean_difference = mean(y, na.rm = TRUE) - mean(x, na.rm = TRUE),
        p.value = tryCatch(t.test(y, x)$p.value, error = function(e) NA_real_)
      )
    })
  }) %>%
  ungroup() %>%
  group_by(PhaseClass, Sex) %>%
  mutate(p_fdr = p.adjust(p.value, method = "BH")) %>%
  ungroup()
write_table(domain_effect_summary, file.path(output_dir, "stats_tables/systems_sis_domain_effect_summary.csv"))

plot_domain_trajectory <- function(domain_name, phase = "Active", title, subtitle, y_lab = "Domain score") {
  plot_tbl <- sis_domain_scores %>%
    filter(Domain == domain_name, PhaseClass == phase) %>%
    group_by(AnimalNum, Group, Sex, CageChangeIndex) %>%
    summarise(DomainScore = mean(DomainScore, na.rm = TRUE), .groups = "drop")
  if (nrow(plot_tbl) == 0) return(ggplot() + annotate("text", x = 0, y = 0, label = paste(title, "unavailable"), size = 3) + theme_void())
  summary_tbl <- plot_tbl %>%
    group_by(Group, Sex, CageChangeIndex) %>%
    summarise(mean = mean(DomainScore, na.rm = TRUE), ci_low = mean_ci(DomainScore)["low"], ci_high = mean_ci(DomainScore)["high"], .groups = "drop")
  ggplot(plot_tbl, aes(CageChangeIndex, DomainScore, group = AnimalNum, colour = Group)) +
    geom_line(alpha = 0.16, linewidth = 0.18) +
    geom_ribbon(data = summary_tbl, aes(y = mean, ymin = ci_low, ymax = ci_high, fill = Group, group = Group), inherit.aes = FALSE, alpha = 0.16, colour = NA) +
    geom_line(data = summary_tbl, aes(y = mean, group = Group), linewidth = 0.58) +
    geom_point(data = summary_tbl, aes(y = mean, group = Group), size = 1.1) +
    facet_wrap(~ Sex, nrow = 1) +
    scale_colour_manual(values = group_colors, drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    labs(title = title, subtitle = subtitle, x = "Regrouping / cage-change order", y = y_lab) +
    make_nature_theme(base_size = 6) +
    theme(legend.position = "top")
}

early_panel_outcome <- if ("outcome" %in% names(first_active_prediction_table) && primary_outcome %in% first_active_prediction_table$outcome) {
  primary_outcome
} else if ("outcome" %in% names(first_active_prediction_table) && nrow(first_active_prediction_table) > 0) {
  first(first_active_prediction_table$outcome)
} else {
  NA_character_
}
early_panel_domain <- if (!is.na(early_panel_outcome) && "outcome" %in% names(first_active_prediction_table)) {
  first_active_prediction_table %>%
    filter(outcome == early_panel_outcome, is.finite(spearman_rho), Domain != "Psychomotor activation") %>%
    arrange(spearman_fdr, desc(abs(spearman_rho))) %>%
    slice_head(n = 1) %>%
    pull(Domain)
} else {
  character(0)
}
if (length(early_panel_domain) == 0 || is.na(early_panel_domain)) early_panel_domain <- "Early adaptation / prediction"
outcome_to_plot_sis <- if (primary_outcome %in% names(first_active_domain_scores)) primary_outcome else intersect(endpoint_cols, names(first_active_domain_scores))[1]

if (!is.na(outcome_to_plot_sis) && outcome_to_plot_sis %in% names(first_active_domain_scores)) {
  early_plot_tbl <- first_active_domain_scores %>%
    filter(Domain == early_panel_domain, is.finite(EarlyDomainScore), is.finite(.data[[outcome_to_plot_sis]])) %>%
    rename(outcome = all_of(outcome_to_plot_sis))
  early_stats_tbl <- early_plot_tbl %>%
    group_by(Sex) %>%
    summarise(
      n = n_distinct(AnimalNum),
      spearman_rho = safe_cor(EarlyDomainScore, outcome, "spearman"),
      spearman_p = safe_cor_p(EarlyDomainScore, outcome, "spearman"),
      .groups = "drop"
    ) %>%
    mutate(
      spearman_q = p.adjust(spearman_p, method = "BH"),
      StatLabel = paste0(rho_label(spearman_rho), "\n", q_label(spearman_q), "\nn=", n)
    )
} else {
  early_plot_tbl <- tibble()
  early_stats_tbl <- tibble()
  outcome_to_plot_sis <- primary_outcome
}

p_sis_early_prediction <- if (nrow(early_plot_tbl) > 0) {
  ggplot(early_plot_tbl, aes(EarlyDomainScore, outcome, colour = Group, fill = Group)) +
    geom_point(size = 1.9, alpha = 0.86, stroke = 0.22) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.34, colour = "grey25", fill = "grey70", alpha = 0.18) +
    geom_text(data = early_stats_tbl, aes(x = -Inf, y = Inf, label = StatLabel), inherit.aes = FALSE, hjust = -0.06, vjust = 1.12, size = 1.9, lineheight = 0.86, colour = "grey15") +
    facet_wrap(~ Sex, nrow = 1) +
    scale_colour_manual(values = group_colors, drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    labs(
      title = "A. Candidate early vulnerability/resilience signal",
      subtitle = paste0("First 12 h ACTIVE phase after ", first_cage_change, ": ", early_panel_domain),
      x = "Early active-phase domain score",
      y = paste0(outcome_to_plot_sis, " (lower = worse)")
    ) +
    make_nature_theme(base_size = 6)
} else {
  ggplot() + annotate("text", x = 0, y = 0, label = "First-active endpoint association unavailable", size = 3) + theme_void()
}
save_plot_svg_pdf(p_sis_early_prediction, file.path(output_dir, "figures/publication_panels/Fig_sis_early_prediction_first_active_12h"), width = 135, height = 78)

p_sis_repeated_adaptation <- plot_domain_trajectory(
  "Repeated adaptation / recovery dynamics",
  "Active",
  "B. Repeated active-phase adaptation",
  "First active phase after each regrouping; trajectories suggest habituation or sensitization",
  "Active adaptation/recovery score"
)
save_plot_svg_pdf(p_sis_repeated_adaptation, file.path(output_dir, "figures/publication_panels/Fig_sis_repeated_active_phase_adaptation"), width = 140, height = 82)

p_sis_phase_heatmap <- domain_effect_summary %>%
  filter(Domain %in% c("Active-phase adaptation/exploration", "Inactive-phase rest/circadian regulation", "Behavioral flexibility / predictability", "Social spatial organization", "Behavioral state architecture", "Behavioral volatility / fragmentation", "Psychomotor activation")) %>%
  mutate(
    Domain = factor(Domain, levels = rev(c("Active-phase adaptation/exploration", "Inactive-phase rest/circadian regulation", "Behavioral flexibility / predictability", "Social spatial organization", "Behavioral state architecture", "Behavioral volatility / fragmentation", "Psychomotor activation"))),
    contrast = factor(contrast, levels = c("RES-CON", "SUS-CON", "SUS-RES"))
  ) %>%
  ggplot(aes(contrast, Domain, fill = hedges_g)) +
  geom_tile(colour = "white", linewidth = 0.22) +
  geom_text(aes(label = sig_label(p_fdr)), size = 1.8) +
  facet_grid(Sex ~ PhaseClass) +
  scale_fill_gradient2(low = "#3d3b6e", mid = "white", high = "#e63947", midpoint = 0, na.value = "grey90") +
  labs(
    title = "C. Active and inactive phase organization",
    subtitle = "Domain-level Hedges g; active/inactive phases are distinct biological regimes",
    x = NULL,
    y = NULL,
    fill = "g"
  ) +
  make_nature_theme(base_size = 5.5) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")
save_plot_svg_pdf(p_sis_phase_heatmap, file.path(output_dir, "figures/publication_panels/Fig_sis_active_inactive_domain_heatmap"), width = 165, height = 102)

p_sis_flexibility <- plot_domain_trajectory(
  "Behavioral flexibility / predictability",
  "Active",
  "D. Behavioral flexibility across regroupings",
  "Entropy-derived flexibility/predictability, separated from psychomotor output",
  "Flexibility / predictability score"
)
save_plot_svg_pdf(p_sis_flexibility, file.path(output_dir, "figures/publication_panels/Fig_sis_behavioral_flexibility_trajectory"), width = 140, height = 82)

p_sis_social <- plot_domain_trajectory(
  "Social spatial organization",
  "Active",
  "E. Social spatial organization",
  "Proximity level and stability; proximity is not direct sociability",
  "Social spatial organization score"
)
save_plot_svg_pdf(p_sis_social, file.path(output_dir, "figures/publication_panels/Fig_sis_social_spatial_organization"), width = 140, height = 82)

p_sis_rest <- plot_domain_trajectory(
  "Inactive-phase rest/circadian regulation",
  "Inactive",
  "F. Sleep-like/rest-like inactivity architecture",
  "Inactive phase recovery and fragmentation; not EEG-confirmed sleep",
  "Rest-like organization score"
)

p_sis_state <- if (!is.null(hmm_occupancy) && nrow(hmm_occupancy) > 0) {
  hmm_occupancy %>%
    standardize_id_columns() %>%
    mutate(
      PhaseClass = case_when(
        str_detect(str_to_lower(as.character(Phase)), "inactive|light|day") ~ "Inactive",
        str_detect(str_to_lower(as.character(Phase)), "active|dark|night") ~ "Active",
        TRUE ~ as.character(Phase)
      ),
      State = as.character(State),
      frac_time = safe_numeric(frac_time)
    ) %>%
    left_join(state_label_tbl %>% select(State, SemanticState), by = "State") %>%
    group_by(Sex, Group, PhaseClass, SemanticState) %>%
    summarise(mean_fraction = mean(frac_time, na.rm = TRUE), .groups = "drop") %>%
    filter(!is.na(SemanticState), is.finite(mean_fraction)) %>%
    ggplot(aes(Group, mean_fraction, fill = Group, colour = Group)) +
    geom_col(width = 0.68, alpha = 0.72, colour = "white", linewidth = 0.18) +
    facet_grid(PhaseClass + Sex ~ SemanticState, scales = "free_y") +
    scale_colour_manual(values = group_colors, drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    labs(title = "F. Behavioral state architecture", subtitle = "HMM state occupancy by phase", x = NULL, y = "Mean occupancy") +
    make_nature_theme(base_size = 5.5) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "none")
} else {
  p_sis_rest
}
save_plot_svg_pdf(p_sis_state, file.path(output_dir, "figures/publication_panels/Fig_sis_rest_or_state_architecture"), width = 150, height = 105)

sis_domain_pca_features <- sis_domain_scores %>%
  group_by(AnimalNum, Group, Sex, Domain) %>%
  summarise(DomainScore = mean(DomainScore, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Domain, values_from = DomainScore)
sis_domain_pca_cols <- intersect(sis_domain_cols, names(sis_domain_pca_features))
if (length(sis_domain_pca_cols) >= 2 && nrow(sis_domain_pca_features) >= 4) {
  sis_domain_pca_mat <- sis_domain_pca_features %>%
    select(all_of(sis_domain_pca_cols)) %>%
    mutate(across(everything(), ~ replace_na(safe_numeric(.x), median(safe_numeric(.x), na.rm = TRUE)))) %>%
    as.matrix()
  sis_domain_pca_scaled <- scale(sis_domain_pca_mat)
  sis_domain_pca_scaled[!is.finite(sis_domain_pca_scaled)] <- 0
  sis_domain_pca <- prcomp(sis_domain_pca_scaled, center = FALSE, scale. = FALSE)
  sis_domain_pca_scores <- sis_domain_pca_features %>%
    select(AnimalNum, Group, Sex) %>%
    bind_cols(as_tibble(sis_domain_pca$x[, 1:2, drop = FALSE]))
  sis_domain_pca_loadings <- as_tibble(sis_domain_pca$rotation[, 1:2, drop = FALSE], rownames = "Domain") %>%
    mutate(abs_PC1 = abs(PC1), locomotion_axis_flag = Domain == "Psychomotor activation" & abs_PC1 == max(abs_PC1, na.rm = TRUE))
  write_table(sis_domain_pca_scores, file.path(output_dir, "tables/systems_sis_domain_pca_scores.csv"))
  write_table(sis_domain_pca_loadings, file.path(output_dir, "tables/systems_sis_domain_pca_loadings.csv"))
}

sis_figure_legend_draft <- tibble(
  Figure = "Fig_integrated_systems_dashboard",
  Panel = c("A", "B", "C", "D", "E", "F"),
  DraftLegend = c(
    paste0("Candidate early vulnerability/resilience signal from the first 12 h active phase after the first regrouping (", first_cage_change, "). Points are animals, colored by later CON/RES/SUS classification; lines show OLS fits for visualization and panel statistics report sex-stratified Spearman correlations with ", outcome_to_plot_sis, "."),
    "Repeated active-phase post-regrouping adaptation/recovery trajectories across cage changes. Thin lines show animals; thick lines and ribbons show group means and 95% confidence intervals.",
    "Domain-level active versus inactive phase organization. Tiles show Hedges g for CON, RES and SUS contrasts, emphasizing that active and inactive phases are separate biological regimes.",
    "Entropy-derived behavioral flexibility/predictability across regroupings, shown separately from psychomotor activation.",
    "Social spatial organization based on proximity level and stability. Proximity is interpreted as spatial co-organization, not direct sociability.",
    "Sleep-like/rest-like inactivity architecture or, when available, HMM behavioral state architecture. Inactivity is RFID-derived and not EEG-confirmed sleep."
  ),
  Caveat = c(
    "Associative early signal only; not a validated biomarker.",
    "Trajectory differences do not prove causal recovery.",
    "RES/SUS contrasts are descriptive because labels are post-endpoint derived.",
    "Entropy can reflect adaptive exploration or disorganization depending on phase and direction.",
    "Graph/proximity metrics are social-spatial proxies.",
    "Use sleep-like/rest-like wording; HMM states are data-derived."
  )
)
write_table(sis_figure_legend_draft, file.path(output_dir, "tables/systems_sis_dashboard_figure_legend_draft.csv"))

sis_dashboard_visualization_rows <- tibble(
  Figure = c(
    "Fig_integrated_systems_dashboard",
    "Fig_sis_early_prediction_first_active_12h",
    "Fig_sis_repeated_active_phase_adaptation",
    "Fig_sis_active_inactive_domain_heatmap",
    "Fig_sis_behavioral_flexibility_trajectory",
    "Fig_sis_social_spatial_organization",
    "Fig_sis_rest_or_state_architecture"
  ),
  PrimaryQuestion = c(
    "How does repeated social instability reorganize behavior across active/inactive contexts without reducing the story to locomotion?",
    "Does the first 12 h active response to the first social perturbation associate with later stress burden?",
    "Do repeated regroupings produce habituation, sensitization or recovery-like trajectories?",
    "Which biological domains differ by group within active versus inactive phases?",
    "Are flexibility/predictability trajectories separable from psychomotor output?",
    "How is social spatial organization altered after regroupings?",
    "Is inactive-phase rest-like organization or behavioral state architecture disrupted?"
  ),
  ManuscriptUse = c(
    "Main figure candidate",
    "Main panel A",
    "Main panel B",
    "Main panel C",
    "Main panel D",
    "Main panel E",
    "Main panel F"
  ),
  Caution = c(
    "Domain scores are interpretable summaries; raw-feature tables remain the audit trail.",
    "Candidate early predictor only; validation requires independent data or stronger cross-validation.",
    "Use descriptive trajectory language unless model statistics support interactions.",
    "Active/inactive effects may reflect different biological processes.",
    "Entropy direction is context-dependent.",
    "Proximity is not equivalent to sociability.",
    "Do not claim EEG sleep or fixed ethological HMM states."
  )
)

if (requireNamespace("patchwork", quietly = TRUE)) {
  sis_dashboard <- ((p_sis_early_prediction | p_sis_repeated_adaptation) /
    (p_sis_phase_heatmap | p_sis_flexibility) /
    (p_sis_social | p_sis_state)) +
    patchwork::plot_layout(heights = c(0.95, 1.05, 1.0), guides = "collect") +
    patchwork::plot_annotation(
      title = "Adolescent social instability stress reorganizes home-cage behavior across circadian contexts",
      subtitle = "Domain-level dashboard prioritizing first-perturbation response, repeated active-phase adaptation, active/inactive organization, social spatial structure and rest-like/state architecture",
      caption = "CON/RES/SUS group contrasts are descriptive because RES/SUS labels are post-paradigm CombZ-derived. The first-active 12 h window is a candidate early vulnerability/resilience signal, not a validated biomarker. Inactivity is sleep-like/rest-like RFID behavior, not EEG-defined sleep."
    )
  save_plot_svg_pdf(sis_dashboard, file.path(output_dir, "figures/Fig_integrated_systems_dashboard"), width = 230, height = 230)
  save_plot_svg_pdf(sis_dashboard, file.path(output_dir, "figures/publication_panels/Fig_sis_systems_neuroscience_dashboard"), width = 230, height = 230)
}

# ------------------------------------------------
# TEXT SUMMARY FOR MANUSCRIPT / LAB MEETING
# ------------------------------------------------

duration_methods_text <- tibble(
  Section = "Methods - duration normalization",
  Text = paste(
    "Observation duration was quantified for every animal x cage-change x phase epoch before downstream behavioral analysis.",
    "Epochs were classified as short when their observed bin count was substantially below the median epoch duration for the corresponding phase, and cage changes were classified without hard-coding cage-change labels.",
    "Per-bin means and medians were retained as duration-robust summaries, whereas cumulative quantities including AUCs, transition counts, switch counts, burst counts, path lengths and interaction counts were converted to per-hour or per-epoch rates.",
    "Temporal-stability metrics such as RMSSD, ACF1, transition entropy and HMM transition summaries were set to missing when too few bins were available for stable estimation.",
    "Major outputs include full-data results and duration-sensitivity tables excluding automatically detected short-duration epochs, so biological CC4 structure is preserved while duration artifacts are not interpreted as phenotypes."
  )
)

systems_visualization_guide <- tibble(
  Figure = c(
    "Fig_integrated_systems_dashboard",
    "Fig_systems_module_scorecard",
    "Fig_systems_named_biological_scores",
    "Fig_systems_hmm_transition_difference",
    "Fig_systems_hmm_state_flow_alluvial",
    "Fig_systems_social_phenotype_map",
    "Fig_systems_trajectory_adaptation_phase_portrait",
    "Fig_systems_prediction_delta_waterfall",
    "Fig_systems_prediction_duration_sensitivity",
    "Fig_systems_module_coupling_network",
    "systems_named_biological_scores.html"
  ),
  PrimaryQuestion = c(
    "What is the smallest coherent systems-neuroscience story supported by robust, integrated outputs?",
    "Which biological modules carry the strongest group effects and how duration-robust are they?",
    "How do animals distribute across named constructs such as rigidity, flexibility, withdrawal and recovery?",
    "Which latent-state transitions differ between SUS, RES and CON?",
    "What are the dominant HMM state flows by group and sex?",
    "Can social withdrawal be separated from fragmented/unstable social topology?",
    "Do recovery dynamics dissociate from behavioral rigidity?",
    "Which module adds predictive value beyond the previous model step?",
    "Does prediction survive exclusion of short-duration epochs?",
    "Which biological modules are coupled or dissociated within sex?",
    "Exploratory hoverable view of animal-level named scores"
  ),
  ManuscriptUse = c(
    "Main figure candidate",
    "Main dashboard subpanel or compact supplement",
    "Supplementary biological interpretation panel",
    "Supplementary HMM/latent-state panel",
    "Supplementary or talk figure if ggalluvial is installed",
    "Supplementary social topology interpretation panel",
    "Supplementary trajectory adaptation panel",
    "Main prediction ladder companion or reviewer supplement",
    "Reviewer robustness panel",
    "Supplementary systems architecture panel",
    "Lab meeting/exploration only"
  ),
  Caution = c(
    "Keeps exploratory nonlinear and state-flow views outside the main composite unless independently central to the claim.",
    "Summarizes strongest effects; use detailed contrast table for exact statistics.",
    "Composite scores are interpretable indices, not independent raw measurements.",
    "Depends on HMM state labeling; semantic labels are data-derived.",
    "Displays top flows only to reduce clutter.",
    "Uses named social composites; inspect source features for mechanistic detail.",
    "Phase portrait is descriptive; use trajectory-feature stats for inference.",
    "LOOCV estimates can be noisy with small n.",
    "Short-excluded model can lose power if many animals contain short epochs.",
    "Correlation is descriptive and does not imply causality.",
    "Interactive HTML is not intended as a manuscript figure."
  )
)

systems_visualization_guide <- bind_rows(
  if (exists("sis_dashboard_visualization_rows")) sis_dashboard_visualization_rows else tibble(),
  tibble(
    Figure = c(
      "Fig_chip_loss_dropout_timeline",
      "Fig_chip_loss_movement_proximity_diagnostics",
      "Fig_batch_system_feature_bias",
      "Fig_group_balance_by_batch_system",
      "Fig_group_sex_interaction_effects",
      "Fig_sleep_like_inactivity_by_group_sex",
      "Fig_first_active_movement_trajectory_by_group_sex",
      "Fig_first_active_entropy_acf1_or_instability",
      "Fig_primary_feature_robustness",
      "Fig_behavior_proteomics_bridge"
    ),
    PrimaryQuestion = c(
      "Which animal epochs show RFID dropout/dead-tag signatures?",
      "Do suspected dropout periods show implausible movement/social-contact collapse?",
      "Are primary feature effects sensitive to batch/system covariates?",
      "Are group and sex balanced across acquisition batch/system?",
      "Which primary features show Group x Sex interaction evidence?",
      "Do groups differ in RFID-derived sleep-like inactivity/quiescence organization?",
      "What is the raw first-active movement trajectory underlying the primary early-behavior claim?",
      "What is the raw first-active entropy trajectory underlying entropy persistence/instability summaries?",
      "Are primary features stable across context and bin-size checks?",
      "Do primary behavioral axes associate with optional proteomics module scores?"
    ),
    ManuscriptUse = c(
      "QC supplement",
      "QC supplement",
      "Reviewer robustness panel",
      "QC supplement",
      "Main or reviewer panel for sex-specific claims",
      "Main sleep-like inactivity panel",
      "Main trajectory panel",
      "Main or supplementary early-entropy panel",
      "Reviewer robustness panel",
      "Optional associative panel only when proteomics table is provided"
    ),
    Caution = c(
      "Flags are conservative and should trigger raw-trace review.",
      "Movement/proximity collapse can reflect biology or tag failure; use with timeline.",
      "Covariate models are underpowered when batch/system cells are sparse.",
      "Imbalance is a design/QC issue, not a phenotype.",
      "Sex-stratified claims need interaction support or descriptive labeling.",
      "Use sleep-like/quiescence-like wording only; no EEG sleep claim.",
      "Raw trajectories are descriptive; statistics are in primary feature tables.",
      "Entropy trajectory supports but does not replace ACF1/RMSSD summaries.",
      "Classifications are conservative reviewer-facing labels.",
      "Proteomics associations are molecular correlates, not causal mechanisms."
    )
  ),
  systems_visualization_guide
) %>%
  mutate(
    LeakageControlLabel = case_when(
      str_detect(Figure, "prediction|CombZ|prospective") ~ "safe_for_prediction_if_pre_endpoint_features_only",
      str_detect(Figure, "proteomics") ~ "associative_or_exploratory",
      str_detect(Figure, "chip|batch|robustness|balance") ~ "QC-only",
      TRUE ~ "descriptive_group_contrast"
    )
  )

systems_output_naming_conventions <- tibble(
  Convention = c(
    "Folder: tables/",
    "Folder: stats_tables/",
    "Folder: tables/qc/",
    "Folder: tables/duration_sensitivity/",
    "Folder: figures/publication_panels/",
    "Folder: figures/supplementary/",
    "Folder: figures/qc/",
    "Folder: figures/exploratory/",
    "Prefix: Fig_systems_",
    "Prefix: systems_",
    "Palette",
    "File formats"
  ),
  Meaning = c(
    "Analysis-ready data tables used by downstream scripts.",
    "Statistical summaries, contrasts, model performance and FDR-corrected results.",
    "Quality-control tables, missingness, duration checks and input audits.",
    "Full versus excluding-short-duration sensitivity outputs.",
    "Static manuscript-style panels exported as SVG, PDF and PNG.",
    "Detailed supportive figures that are too dense for the main figure.",
    "Diagnostic figures for model/data-quality checks.",
    "Discovery and interactive/exploratory figures; do not cite as primary inference.",
    "Publication-facing figure generated by the integrated systems script.",
    "Publication-facing table generated by the integrated systems script.",
    "CON = deep indigo, RES = warm grey, SUS = red; use the same palette across all panels.",
    "SVG/PDF are manuscript-preferred; PNG is exported for quick previews and slide decks."
  ),
  ReviewerBenefit = c(
    "Separates raw analysis products from statistical claims.",
    "Makes inferential outputs easy to find.",
    "Makes preprocessing and missingness assumptions auditable.",
    "Makes unequal-duration robustness explicit.",
    "Keeps main figures visually consistent.",
    "Prevents main figures from becoming overcrowded.",
    "Provides transparency without bloating manuscript figures.",
    "Keeps modern visualizations available but clearly labeled.",
    "Easy figure tracking in manuscripts and slide decks.",
    "Easy table tracking in supplements.",
    "Prevents group-color drift across analyses.",
    "Supports both publication production and fast review."
  )
)

duration_sensitivity_audit <- tibble::tribble(
  ~script, ~output_file, ~metric, ~duration_sensitive, ~reason, ~correction_applied, ~normalized_metric_name, ~cc4_exclusion_sensitivity_available, ~reviewer_risk,
  "03_build_multiscale_behavior_metrics.R", "all_behavior_metrics.csv", "Movement, MovementDistance, ProximitySeconds", TRUE, "Raw counts/seconds scale with observed duration.", "Added observation-duration QC and per-hour rate columns.", "MovementPerHour; MovementDistancePerHour; ProximitySecondsPerHour", TRUE, "high",
  "04_gamm_movement_proximity_phase_and_early_window.R", "gamm AUC tables", "GAMM trajectory AUC", TRUE, "GAMM AUC integrates over the available time grid and can scale with trajectory duration.", "Observation-duration QC exported and AUC tables include per-hour normalized values.", "AUC_per_hour; AUC_diff_per_hour", TRUE, "high",
  "05_build_dyadic_rfid_contacts.R", "dyadic_network_ready.csv", "same/adjacent position seconds", TRUE, "Dyadic contact seconds scale with observed duration and feed network edge weights.", "Added epoch duration QC, BinSizeSec, and per-hour contact-second rates.", "same_position_seconds_per_hour; adjacent_seconds_per_hour", TRUE, "high",
  "06_burstiness_temporal_instability.R", "temporal_instability_metrics_per_animal_all_metrics.csv", "RMSSD, ACF1, CV, Fano", TRUE, "Temporal estimates become unstable with short epochs.", "Epoch duration QC joined; helper safeguards retain minimum-bin NA rules; excluding-short-duration table exported.", "same metric, duration-tagged", TRUE, "medium",
  "07_behavioral_state_space.R", "state_switching_metrics.csv", "n_switches, n_transitions, switch_rate", TRUE, "Raw transition counts scale with number of bins.", "Counts converted to per-hour rates; transition metrics require minimum bins.", "n_switches_per_hour; n_transitions_per_hour", TRUE, "high",
  "07_behavioral_state_space.R", "state_diversity_metrics.csv", "state entropy", TRUE, "Entropy stability depends on epoch length.", "Duration QC joined and entropy_rate exported.", "entropy_rate", TRUE, "medium",
  "08_early_prediction_models.R", "early_behavior_features.csv", "early-window features", TRUE, "Prediction features can leak duration if QC columns enter the model.", "Duration QC exported; duration columns excluded from predictors; short-duration feature table exported.", "early_behavior_features_excluding_short_duration.csv", TRUE, "medium",
  "08b_early_prediction_model_ladder.R", "model_ladder_performance_duration_sensitivity.csv", "prediction ladder performance", TRUE, "LOOCV performance may change if short-duration early-window rows contribute different feature reliability.", "Full-data and excluding-short-duration prediction ladders exported with performance deltas.", "duration-sensitive prediction-performance comparison", TRUE, "medium",
  "09_dynamic_social_networks.R", "animal_level_social_dynamics.csv", "contact switches", TRUE, "Switch counts scale with observation duration.", "Contact switches converted to per-hour rates and duration-sensitivity table exported.", "contact_switch_count_per_hour", TRUE, "high",
  "09_dynamic_social_networks.R", "dyadic_pair_summary.csv", "contact bins", TRUE, "Dyadic contact counts scale with observed bins.", "Contact bins converted to per-hour rates.", "contact_bins_per_hour", FALSE, "high",
  "10_hmm_behavioral_states.R", "hmm_transition_counts.csv", "Transitions", TRUE, "Raw HMM transitions scale with sequence length.", "Minimum sequence length raised and transitions converted to per-hour rates.", "Transitions_per_hour", TRUE, "high",
  "10_hmm_behavioral_states.R", "hmm_state_dwell_times.csv", "dwell bins", TRUE, "Dwell time in bins depends on bin size and available duration.", "Dwell bins converted to hours and duration QC joined.", "mean_dwell_hours; median_dwell_hours; max_dwell_hours", TRUE, "medium",
  "11_gamm_trajectory_features.R", "combined_gamm_features.csv", "AUC", TRUE, "AUC is cumulative over trajectory duration.", "AUC set to NA for short trajectories and converted to per-hour rate.", "auc_per_hour", TRUE, "high",
  "15_behavioral_adaptation_kinetics.R", "adaptation_kinetics_features.csv", "stabilization time, volatility decay, recovery slope", TRUE, "Recovery/stabilization estimates depend on enough bins to estimate early and late trajectory structure.", "Epoch duration QC joined; short-duration sensitivity output exported.", "recovery_slope_per_hour; stabilization_time_hours; adaptation_half_life_hours", TRUE, "medium",
  "16_sleep_like_inactivity_metrics.R", "sleep_like_inactivity_features.csv", "bout counts and transition rates", TRUE, "Inactivity bout counts and prolonged episodes scale with observation duration.", "Bout/prolonged episode counts converted to rates and duration sensitivity output exported.", "inactivity_bout_count_per_hour; prolonged_inactivity_episodes_per_hour", TRUE, "medium",
  "17_ethological_phase_organization.R", "phase_contrast_features.csv", "active/inactive phase contrasts", TRUE, "Phase organization estimates require comparable active and inactive observation structure.", "Phase-level duration QC exported and full/excluding-short-duration feature tables written.", "active_minus_inactive_mean; active_inactive_ratio_mean", TRUE, "medium",
  "12_systems_neuroscience_summary.R", "systems_temporal_latent_epoch_embeddings.csv", "PCA/UMAP/PHATE epoch state space", TRUE, "Embedding can be biased by short/cumulative epoch summaries.", "Duration QC joined; short epochs excluded from latent trajectory embedding.", "duration-tagged normalized features", TRUE, "medium",
  "12_systems_neuroscience_summary.R", "systems_pca_scores.csv; systems_umap_scores.csv", "animal-level systems state space", TRUE, "Unnormalized cumulative features can dominate embedding axes when observation duration differs.", "PCA/UMAP use duration-robust feature pool excluding raw cumulative/count summaries.", "duration_robust_embedding_prediction feature set", TRUE, "medium",
  "12_systems_neuroscience_summary.R", "systems_latent_instability_by_animal.csv", "latent path length and roughness", TRUE, "Path length scales with number of epochs.", "Path length normalized per hour and per epoch; roughness normalized by epoch count.", "latent_path_length_per_hour; latent_path_length_per_epoch; latent_roughness_normalized", TRUE, "high",
  "12_systems_neuroscience_summary.R", "systems_prediction_ladder_performance_duration_sensitivity.csv", "integrated systems prediction ladder", TRUE, "Module scores can inherit duration artifacts if raw cumulative features are included.", "Prediction ladder uses duration-robust feature pool and exports full vs excluding-short-duration performance.", "duration-robust module scores", TRUE, "medium"
)

systems_claim_hierarchy <- tibble(
  EvidenceTier = c("Tier 1", "Tier 1", "Tier 1", "Tier 2", "Tier 2", "Tier 2", "Tier 2", "Tier 2", "Tier 3"),
  BiologicalDomain = c(
    "Early behavioral adaptation",
    "Temporal instability and organization",
    "Longitudinal adaptation and recovery",
    "Ethological active/inactive phase organization",
    "Social reorganization after regrouping",
    "Behavioral state organization",
    "Sleep-like inactivity/quiescence",
    "Behavioral-proteomic systems alignment",
    "Exploratory nonlinear systems dynamics"
  ),
  CentralClaim = c(
    "Early first-active-phase behavior after the first social instability exposure predicts later stress burden.",
    "Temporal organization contributes information beyond movement magnitude.",
    "Adaptive stabilization/recovery after regrouping is a biologically interpretable resilience axis.",
    "Stress phenotypes can alter active/inactive behavioral structure without claiming validated circadian or sleep disruption.",
    "Regrouping induces social reorganization that should be separated from simple proximity magnitude.",
    "HMM/state-space outputs decompose behavioral regimes and flexibility but do not define the primary claim.",
    "RFID inactivity structure supports rest-like/quiescence interpretations, not EEG sleep claims.",
    "Low-dimensional behavioral axes can be compared with curated hippocampal molecular adaptation modules.",
    "Complexity, attractor and manifold views are supplementary hypothesis-generating visualizations."
  ),
  PrimaryEvidence = c(
    "08b behavior-only prediction ladder; Movement_mean, Movement_rmssd, Entropy_acf1.",
    "06 temporal instability metrics; early prediction associations.",
    "11/15 trajectory and adaptation kinetics features.",
    "17 phase contrast, timing, fragmentation and recovery outputs.",
    "09 animal-level social dynamics plus dyadic threshold sensitivity if available.",
    "10 HMM occupancy, dwell and transition outputs.",
    "16 sleep-like inactivity features.",
    "12_behavior_proteomics_integration low-dimensional axis associations.",
    "14 next-generation behavioral phenotyping outputs."
  ),
  ClaimType = c("predictive", "predictive", "descriptive", "descriptive", "descriptive", "descriptive", "descriptive", "associative", "exploratory"),
  AllowedInterpretation = c(
    "Prospective behavioral adaptation dynamics are associated with later stress burden.",
    "Resilience/susceptibility may differ in temporal organization, not only mean locomotion.",
    "Animals differ in stabilization or convergence after social perturbation.",
    "Phase-specific behavioral organization is altered or preserved.",
    "Social engagement, partner stability and topology are separable constructs.",
    "Latent state persistence/flexibility provides mechanistic decomposition.",
    "Quiescence fragmentation is an RFID-derived rest-like behavioral phenotype.",
    "Behavioral adaptation axes align with molecular adaptation modules.",
    "Nonlinear systems signatures suggest candidate dynamics for follow-up."
  ),
  DisallowedInterpretation = c(
    "Behavior causes later stress burden.",
    "ACF1 or RMSSD alone proves resilience.",
    "Recovery is causal repair.",
    "Circadian disruption or sleep disruption without validation.",
    "Graph-theory social claims when dyadic identity is unavailable.",
    "HMM states are fixed ethological states without validation.",
    "EEG sleep architecture.",
    "Molecular mechanism is proven by correlation.",
    "Attractor/manifold metrics are validated biomarkers."
  )
)

systems_claim_hierarchy <- systems_claim_hierarchy %>%
  mutate(
    LeakageControlLabel = case_when(
      ClaimType == "predictive" ~ "prospective_endpoint_association_requires_pre_endpoint_features",
      ClaimType == "associative" ~ "associative_not_causal_or_prospective_without_temporal_justification",
      ClaimType == "exploratory" ~ "exploratory_supplement_only",
      TRUE ~ "descriptive_group_contrast"
    )
  )

systems_interpretation_guide <- feature_dictionary %>%
  distinct(
    Feature, SourceScript, BiologicalDomain, MathematicalDefinition, TimeWindow,
    BinLevel, EvidenceTier, ClaimType, AllowedInterpretation, ReviewerRisk,
    DurationSensitive, StableForMainText
  ) %>%
  arrange(EvidenceTier, BiologicalDomain, ReviewerRisk, Feature)

systems_primary_findings_summary <- group_contrasts %>%
  left_join(feature_dictionary %>% select(feature, Feature, BiologicalDomain, EvidenceTier, ClaimType, AllowedInterpretation, StableForMainText), by = "feature") %>%
  mutate(
    StableForMainText = StableForMainText %in% TRUE & !is.na(p_fdr) & abs(hedges_g) >= 0.50,
    ReportingPriority = case_when(
      EvidenceTier == "Tier 1" & StableForMainText ~ "main_text_candidate",
      EvidenceTier == "Tier 2" & StableForMainText ~ "supplement_or_mechanistic_panel",
      EvidenceTier == "Tier 3" ~ "exploratory_supplement",
      TRUE ~ "supporting_or_unstable"
    )
  ) %>%
  arrange(EvidenceTier, desc(StableForMainText), p_fdr, desc(abs(hedges_g))) %>%
  select(
    Feature, Source, BiologicalDomain, EvidenceTier, ClaimType, Sex, contrast,
    mean_ref, mean_comp, estimate, cohen_d, hedges_g, p.value, p_fdr,
    evidence, AllowedInterpretation, StableForMainText, ReportingPriority
  ) %>%
  slice_head(n = 100)

systems_robustness_summary <- {
  prediction_robust <- if (exists("prediction_perf_duration_sensitivity")) {
    prediction_perf_duration_sensitivity %>%
      select(any_of(c("Model", "DurationAnalysisSet", "cv_r2", "pearson_r", "delta_cv_r2_vs_full", "delta_pearson_r_vs_full", "delta_rmse_vs_full"))) %>%
      mutate(
        Analysis = "systems_prediction_ladder",
        Metric = Model,
        full_data_effect = if_else(DurationAnalysisSet == "full", cv_r2, NA_real_),
        excluding_short_duration_effect = if_else(DurationAnalysisSet == "excluding_short_duration", cv_r2, NA_real_),
        delta_estimate = delta_cv_r2_vs_full,
        delta_cohen_d = NA_real_,
        direction_stable = sign(cv_r2) == sign(cv_r2 - delta_cv_r2_vs_full),
        StableForMainText = direction_stable %in% TRUE & (is.na(delta_cohen_d) | abs(delta_cohen_d) < 0.30)
      ) %>%
      select(Analysis, Metric, DurationAnalysisSet, full_data_effect, excluding_short_duration_effect, delta_estimate, delta_cohen_d, direction_stable, StableForMainText)
  } else tibble()

  contrast_robust <- group_contrasts %>%
    left_join(feature_dictionary %>% select(feature, BiologicalDomain, EvidenceTier), by = "feature") %>%
    group_by(feature, Sex, contrast) %>%
    summarise(
      Analysis = first(BiologicalDomain),
      Metric = first(feature),
      full_data_effect = first(hedges_g),
      excluding_short_duration_effect = NA_real_,
      delta_estimate = NA_real_,
      delta_cohen_d = NA_real_,
      direction_stable = TRUE,
      StableForMainText = EvidenceTier[1] %in% c("Tier 1", "Tier 2") & abs(full_data_effect) >= 0.50,
      .groups = "drop"
    ) %>%
    select(Analysis, Metric, full_data_effect, excluding_short_duration_effect, delta_estimate, delta_cohen_d, direction_stable, StableForMainText)

  bind_rows(prediction_robust, contrast_robust)
}

write_table(duration_methods_text, file.path(output_dir, "tables/duration_normalization_methods_text.csv"))
write_table(systems_visualization_guide, file.path(output_dir, "tables/systems_visualization_guide.csv"))
write_table(systems_output_naming_conventions, file.path(output_dir, "tables/systems_output_naming_conventions.csv"))
write_table(duration_sensitivity_audit, file.path(output_dir, "tables/duration_sensitivity_audit.csv"))
write_table(systems_claim_hierarchy, file.path(output_dir, "tables/systems_claim_hierarchy.csv"))
write_table(systems_primary_findings_summary, file.path(output_dir, "tables/systems_primary_findings_summary.csv"))
write_table(systems_robustness_summary, file.path(output_dir, "tables/systems_robustness_summary.csv"))
write_table(systems_interpretation_guide, file.path(output_dir, "tables/systems_interpretation_guide.csv"))

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
    "Duration-robust embedding/prediction features",
    "Prospective prediction features",
    "Named biological scores",
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
    as.character(length(duration_robust_features)),
    as.character(length(prospective_prediction_features)),
    as.character(n_distinct(named_biological_scores_long$Metric)),
    primary_bin_level,
    paste0("First active 12 h after ", first_cage_change),
    primary_outcome_label,
    "Derived from post-paradigm CombZ; group contrasts are descriptive phenotype contrasts, not independent endpoint validation.",
    paste(sort(unique(feature_dictionary$Domain)), collapse = "; "),
    "Use PCA/UMAP for duration-robust descriptive system-level organization. Use module-level prediction ladder outputs and short-duration sensitivity tables for incremental predictive-value claims."
      )
    )

write_table(strong_findings, file.path(output_dir, "tables/systems_top_findings_for_reporting.csv"))
write_table(systems_summary, file.path(output_dir, "tables/systems_analysis_summary.csv"))

stats_reporting_guide <- tibble(
  Output = c(
    "systems_sis_domain_scores.csv",
    "systems_sis_first_active_12h_prediction_table.csv",
    "systems_sis_domain_mixed_model_stats.csv",
    "systems_sis_feature_redundancy_diagnostics.csv",
    "systems_sis_dashboard_figure_legend_draft.csv",
    "systems_group_summary.csv",
    "systems_group_contrasts.csv",
    "systems_outcome_associations.csv",
    "systems_prospective_outcome_associations.csv",
    "systems_prospective_outcome_loo_performance.csv",
    "systems_prediction_ladder_performance.csv",
    "systems_prediction_ladder_performance_duration_sensitivity.csv",
    "systems_module_scorecards.csv",
    "systems_named_biological_scores.csv",
    "systems_duration_negative_control_endpoint_correlations.csv",
    "systems_module_feature_inventory.csv"
  ),
  PrimaryStatistic = c(
    "Sex/phase/cage-change z-scored biologically interpretable SIS domain scores",
    "Sex-stratified Spearman and Pearson associations for first 12 h ACTIVE phase domain scores",
    "Group x Sex x Phase mixed-model coefficients with raw and FDR-adjusted p-values",
    "Spearman feature-feature redundancy and locomotion-dominance flags",
    "Publication-ready draft legend text and caveats for the SIS dashboard panels",
    "Mean, 95% CI, median, IQR, animal n",
    "Welch mean difference, 95% CI, Hedges g, Cohen d, Wilcoxon p",
    "Spearman rho with approximate 95% CI; Pearson r also reported",
    "Spearman rho with approximate 95% CI for pre-endpoint early features only",
    "Integrated module-level LOOCV r, rho, RMSE, MAE and CV R2",
    "LOOCV r, rho, RMSE, MAE, CV R2, permutation p, bootstrap CIs and delta CV R2",
    "Matched full-data and excluding-short-duration LOOCV performance with deltas versus full data",
    "Feature count, duration robustness, strongest group effect and prediction readout by module",
    "Named z-scored biological constructs such as rigidity, flexibility, social withdrawal and recovery",
    "Spearman association between duration/QC variables and endpoint variables",
    "Number of usable features per biologically interpretable module"
  ),
  MultipleTesting = c(
    "Not applicable; domain-score construction table",
    "BH FDR within outcome x sex across early active-domain scores",
    "BH FDR within each displayed domain model",
    "Diagnostic only; no inferential claim",
    "Not applicable",
    "Not applicable",
    "BH FDR within Sex x Source x Domain x Scale x Context x contrast; broad Sex x contrast FDR provided",
    "BH FDR within outcome x descriptive feature set",
    "BH FDR within outcome x prospective feature set",
    "Cross-validation, no per-feature multiplicity claim",
    "Cross-validation plus permutation test per ladder model",
    "Sensitivity analysis; compare effect direction and predictive-performance deltas",
    "Descriptive module-level summary; detailed tests in group contrast and prediction tables",
    "BH FDR available through the standard group-contrast/outcome-association outputs when scores enter the feature matrix",
    "BH FDR within endpoint; should be interpreted as a negative-control screen",
    "Not applicable"
  ),
  ManuscriptUse = c(
    "Primary dashboard source table",
    "Primary early-prediction panel statistics",
    "Primary active/inactive phase model statistics",
    "Reviewer-facing redundancy and locomotion-dominance audit",
    "Figure legend drafting",
    "Figure legends and Supplementary tables",
    "Primary descriptive phenotype effect-size table",
    "Exploratory full-experiment endpoint associations",
    "Prospective endpoint association claims",
    "Integrated systems prediction performance claims",
    "Incremental predictive value beyond movement magnitude",
    "Reviewer-facing robustness check for unequal observation duration",
    "Main/supplementary module interpretability table",
    "Main biological construct table and raincloud-style figure",
    "Reviewer-facing duration confound check",
    "Feature architecture and redundancy-control reporting"
  )
)

write_table(stats_reporting_guide, file.path(output_dir, "tables/systems_stats_reporting_guide.csv"))

message("Integrated systems neuroscience summary complete.")
message("Output: ", output_dir)
message("Primary figure candidates:")
message("  - figures/publication_panels/Fig_systems_state_space_PCA.svg")
message("  - figures/publication_panels/Fig_systems_effect_size_heatmap.svg")
message("  - figures/publication_panels/Fig_systems_module_scorecard.svg")
message("  - figures/publication_panels/Fig_systems_named_biological_scores.svg")
message("  - figures/publication_panels/Fig_systems_hmm_transition_difference.svg")
message("  - figures/publication_panels/Fig_systems_social_phenotype_map.svg")
message("  - figures/publication_panels/Fig_systems_trajectory_adaptation_phase_portrait.svg")
message("  - figures/publication_panels/Fig_systems_prediction_delta_waterfall.svg")
message("  - figures/publication_panels/Fig_systems_prediction_ladder.svg")
message("  - figures/Fig_integrated_systems_dashboard.svg")
message("  - figures/publication_panels/Fig_sis_systems_neuroscience_dashboard.svg")
message("  - figures/publication_panels/Fig_sis_early_prediction_first_active_12h.svg")
