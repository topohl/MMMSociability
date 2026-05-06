#' @title Reviewer-ready analysis of Shannon entropy in animal positioning
#'
#' @description
#' Full adjusted pipeline for Shannon entropy analysis of animal spatial positioning.
#' This version preserves the original workflow but reorganizes the inferential layer into
#' a more defensible framework:
#'   1. Raw/preprocessed entropy extraction, phase-based and half-hour based.
#'   2. Clean consecutive phase and half-hour alignment.
#'   3. Primary mixed-effects model for individual entropy.
#'   4. Secondary cage-level mixed model.
#'   5. Conservative time-resolved GAMM with residual autocorrelation diagnostics.
#'   6. Entropy-based exploration metrics at animal level.
#'   7. Animal-level transition-network metrics instead of pooled edgewise group tests.
#'   8. Publication-ready SVG figures and structured tables.
#'
#' @details
#' Important assumptions:
#'   - 8 spatial positions are used for Shannon entropy; maximum entropy = log2(8) = 3.
#'   - 24 half-hours correspond to one 12 h active or inactive phase.
#'   - If half-hour files already contain true ConsecActive/ConsecInactive/Phase columns,
#'     those are preferred over inferred phase labels.
#'   - Entropy-based exploration is defined as animalEntropy > 0. This is not identical to
#'     outside-home-cage exploration unless home-cage occupancy is explicitly computed.
#'
#' @date Revised version, May 2026
#' @authors Tobias Pohl, Anja Magister

# ===================================================================
# 0. PACKAGE LOADING
# ===================================================================
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}

# Use only packages needed for the pipeline.
# DHARMa and kableExtra are intentionally omitted here because they often
# create avoidable installation problems on locked-down Windows/R setups.
pacman::p_load(
  readr,
  dplyr,
  tidyr,
  tibble,
  purrr,
  stringr,
  lubridate,
  ggplot2,
  patchwork,
  scales,
  svglite,
  lmerTest,
  emmeans,
  broom,
  broom.mixed,
  mgcv,
  performance,
  igraph,
  ggraph,
  jsonlite
)

# ===================================================================
# 1. RUNTIME CONFIGURATION
# ===================================================================
load_existing_data <- TRUE
show_plots <- FALSE
save_plots <- TRUE
save_tables <- TRUE
analyze_by_halfhour <- TRUE
exclude_homecage <- TRUE

# Reviewer-oriented defaults
RUN_PRIMARY_LMM <- TRUE
RUN_CAGE_LMM <- FALSE
RUN_GAMM <- TRUE
RUN_CC1_ACTIVE_GAMM <- TRUE
RUN_EXPLORATION_MODELS <- TRUE
RUN_TRANSITION_METRICS <- TRUE

# CC4 filtering
filter_cc4_late_phases <- TRUE
cc4_max_active_phase <- 2
cc4_max_inactive_phase <- 2

# GAMM defaults
GAMM_K_GLOBAL <- 6L
GAMM_K_GROUP_DEVIATION <- 4L
GAMM_K_HALFHOUR_PHASE <- 5L
GAMM_K_CC1_ACTIVE <- 6L
GAMM_GAMMA <- 1.5
GAMM_RHO_MIN <- 0
GAMM_RHO_MAX <- 0.8
GAMM_ESTIMATE_RHO <- TRUE

# Entropy constants
n_positions <- 8
max_position_entropy <- log2(n_positions)
low_entropy_cutoff <- max_position_entropy / 3
moderate_entropy_cutoff <- 2 * max_position_entropy / 3

# Batches and cage changes
batches <- c("B1", "B2", "B3", "B4", "B5", "B6")
cageChanges <- c("CC1", "CC2", "CC3", "CC4")
female_batches <- c("B3", "B4", "B6")

# ===================================================================
# 2. DIRECTORY STRUCTURE
# ===================================================================
working_directory <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability"
saving_directory <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability/new_structure"

raw_data_dir <- file.path(working_directory, "raw_data")
preprocessed_data_dir <- file.path(working_directory, "preprocessed_data")

# -------------------------------------------------------------------
# Legacy folders retained for compatibility with previous outputs
# -------------------------------------------------------------------
tables_dir <- file.path(saving_directory, "tables")
plots_dir <- file.path(saving_directory, "plots")
logs_dir <- file.path(saving_directory, "logs")

combined_tables_dir <- file.path(tables_dir, "combined")
batch_cage_tables_dir <- file.path(tables_dir, "batch_cagechange")
exploration_metrics_dir <- file.path(tables_dir, "exploration_metrics")
metadata_dir <- file.path(combined_tables_dir, "metadata")

combined_plots_dir <- file.path(plots_dir, "combined")
batch_cage_plots_dir <- file.path(plots_dir, "batch_cagechange")
misc_plots_dir <- file.path(plots_dir, "miscellaneous")
transition_networks_dir <- file.path(combined_plots_dir, "transition_networks")

# -------------------------------------------------------------------
# Reviewer-ready output structure
# -------------------------------------------------------------------
review_dir <- file.path(saving_directory, "analysis_ready_reviewed")

# Clean and intermediate data
review_data_dir <- file.path(review_dir, "01_data")
review_raw_loaded_dir <- file.path(review_data_dir, "01_raw_loaded")
review_clean_data_dir <- file.path(review_data_dir, "02_clean")
review_derived_data_dir <- file.path(review_data_dir, "03_derived")

# Statistics and model outputs
review_stats_dir <- file.path(review_dir, "02_statistics")
review_descriptive_stats_dir <- file.path(review_stats_dir, "01_descriptive")
review_lmm_stats_dir <- file.path(review_stats_dir, "02_lmm")
review_gamm_stats_dir <- file.path(review_stats_dir, "03_gamm")
review_cc1_active_gamm_stats_dir <- file.path(review_gamm_stats_dir, "01_cc1_active_phase1")
review_emm_stats_dir <- file.path(review_stats_dir, "04_emmeans_contrasts")
review_effect_stats_dir <- file.path(review_stats_dir, "05_effect_sizes")
review_exploration_stats_dir <- file.path(review_stats_dir, "06_exploration")
review_transition_stats_dir <- file.path(review_stats_dir, "07_transition_metrics")

# Serialized model objects and diagnostics
review_model_dir <- file.path(review_dir, "03_models_rds")
review_diagnostics_dir <- file.path(review_dir, "04_diagnostics")
review_lmm_diagnostics_dir <- file.path(review_diagnostics_dir, "01_lmm")
review_gamm_diagnostics_dir <- file.path(review_diagnostics_dir, "02_gamm")
review_cc1_active_gamm_diagnostics_dir <- file.path(review_gamm_diagnostics_dir, "01_cc1_active_phase1")
review_plot_diagnostics_dir <- file.path(review_diagnostics_dir, "03_diagnostic_plots")

# Figures: explicit publication hierarchy
review_figures_dir <- file.path(review_dir, "05_figures")
review_main_figures_dir <- file.path(review_figures_dir, "01_main_figures")
review_main_panels_dir <- file.path(review_main_figures_dir, "single_panels")
review_main_multipanel_dir <- file.path(review_main_figures_dir, "multi_panel")
review_supp_figures_dir <- file.path(review_figures_dir, "02_supplementary_figures")
review_supp_panels_dir <- file.path(review_supp_figures_dir, "single_panels")
review_supp_multipanel_dir <- file.path(review_supp_figures_dir, "multi_panel")
review_exploratory_figures_dir <- file.path(review_figures_dir, "03_exploratory_quality_control")
review_gamm_figures_dir <- file.path(review_figures_dir, "04_gamm_timecourse")
review_cc1_active_gamm_figures_dir <- file.path(review_gamm_figures_dir, "01_cc1_active_phase1")
review_transition_figures_dir <- file.path(review_figures_dir, "05_transition_networks")

# Tables for manuscript/supplement
review_tables_dir <- file.path(review_dir, "06_tables_for_manuscript")
review_main_tables_dir <- file.path(review_tables_dir, "01_main_tables")
review_supp_tables_dir <- file.path(review_tables_dir, "02_supplementary_tables")
review_machine_tables_dir <- file.path(review_tables_dir, "03_machine_readable_full_outputs")

# Logs and metadata
review_logs_dir <- file.path(review_dir, "07_logs_metadata")
review_transition_dir <- file.path(review_dir, "08_transition_networks_full")

all_dirs <- c(
  tables_dir,
  plots_dir,
  logs_dir,
  combined_tables_dir,
  file.path(combined_tables_dir, "phase_based"),
  file.path(combined_tables_dir, "halfhour_based"),
  batch_cage_tables_dir,
  exploration_metrics_dir,
  metadata_dir,
  combined_plots_dir,
  batch_cage_plots_dir,
  misc_plots_dir,
  transition_networks_dir,
  review_dir,
  review_data_dir,
  review_raw_loaded_dir,
  review_clean_data_dir,
  review_derived_data_dir,
  review_stats_dir,
  review_descriptive_stats_dir,
  review_lmm_stats_dir,
  review_gamm_stats_dir,
  review_cc1_active_gamm_stats_dir,
  review_emm_stats_dir,
  review_effect_stats_dir,
  review_exploration_stats_dir,
  review_transition_stats_dir,
  review_model_dir,
  review_diagnostics_dir,
  review_lmm_diagnostics_dir,
  review_gamm_diagnostics_dir,
  review_cc1_active_gamm_diagnostics_dir,
  review_plot_diagnostics_dir,
  review_figures_dir,
  review_main_figures_dir,
  review_main_panels_dir,
  review_main_multipanel_dir,
  review_supp_figures_dir,
  review_supp_panels_dir,
  review_supp_multipanel_dir,
  review_exploratory_figures_dir,
  review_gamm_figures_dir,
  review_cc1_active_gamm_figures_dir,
  review_transition_figures_dir,
  review_tables_dir,
  review_main_tables_dir,
  review_supp_tables_dir,
  review_machine_tables_dir,
  review_logs_dir,
  review_transition_dir
)

for (dir_path in all_dirs) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

# ===================================================================
# 3. LOGGING
# ===================================================================
processing_start_time <- Sys.time()
log_file <- file.path(review_logs_dir, "reviewed_entropy_pipeline_log.txt")
if (file.exists(log_file)) file.remove(log_file)

log_message <- function(message_text, batch = NULL, cageChange = NULL, phase = NULL) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  context <- paste(Filter(Negate(is.null), list(batch, cageChange, phase)), collapse = " | ")
  if (nzchar(context)) {
    full_message <- paste0("[", timestamp, "] ", context, ": ", message_text)
  } else {
    full_message <- paste0("[", timestamp, "] ", message_text)
  }
  message(full_message)
  write(full_message, file = log_file, append = TRUE)
}

log_message("Starting reviewer-ready Shannon entropy pipeline")

# ===================================================================
# 4. PLOTTING AND OUTPUT UTILITIES
# ===================================================================
group_cols <- c(
  con = "#3d3b6e",
  res = "#c6c3bb",
  sus = "#e63947"
)

sex_cols <- c(
  male = "#6F6F91",
  female = "#C9B27C"
)

sem <- function(x) {
  stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
}

# Use base_family = "sans" instead of "Arial" for Windows-safe export.
# This avoids 'invalid font type' errors during SVG/PDF rendering.
theme_publication <- function(base_size = 7, base_family = "sans") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.28, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.28, colour = "black"),
      axis.ticks.length = grid::unit(1.4, "mm"),
      axis.text = ggplot2::element_text(size = base_size, colour = "black"),
      axis.title = ggplot2::element_text(size = base_size + 1, colour = "black"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = base_size + 1, colour = "black"),
      plot.title = ggplot2::element_text(size = base_size + 2, face = "plain", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = base_size, hjust = 0),
      legend.title = ggplot2::element_text(size = base_size),
      legend.text = ggplot2::element_text(size = base_size),
      legend.key.size = grid::unit(3, "mm"),
      panel.grid = ggplot2::element_blank(),
      plot.margin = grid::unit(c(2, 2, 2, 2), "mm")
    )
}

theme_nature_gamm <- function(base_size = 7, base_family = "sans") {
  theme_publication(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      legend.position = "top",
      legend.justification = "left",
      legend.box.spacing = grid::unit(0.5, "mm"),
      panel.spacing = grid::unit(2.4, "mm"),
      plot.title = ggplot2::element_text(size = base_size + 1, face = "plain", hjust = 0),
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_text(size = base_size - 1, colour = "grey25", hjust = 0),
      plot.caption.position = "plot"
    )
}

format_group_labels <- function(x) {
  factor(x, levels = c("con", "res", "sus"), labels = c("CON", "RES", "SUS"))
}

sanitize_filename <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "") %>%
    stringr::str_to_lower()
}

analysis_file <- function(prefix, descriptor, level = NULL, ext = "csv") {
  parts <- c(prefix, descriptor, level)
  parts <- parts[!is.na(parts) & nzchar(parts)]
  paste0(paste(sanitize_filename(parts), collapse = "__"), ".", ext)
}

write_review_csv <- function(data, directory, prefix, descriptor, level = NULL) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(directory, analysis_file(prefix, descriptor, level, ext = "csv"))
  readr::write_csv(data, path)
  invisible(path)
}

save_review_figure <- function(plot, directory, prefix, descriptor, level = NULL,
                               width = 90, height = 75, units = "mm", ext = "svg") {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  filename <- analysis_file(prefix, descriptor, level, ext = ext)
  path <- file.path(directory, filename)
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    units = units,
    dpi = 600,
    bg = "white"
  )
  invisible(path)
}

save_review_figure_dual <- function(plot, directory, prefix, descriptor, level = NULL,
                                    width = 90, height = 75, units = "mm") {
  svg_path <- save_review_figure(
    plot = plot,
    directory = directory,
    prefix = prefix,
    descriptor = descriptor,
    level = level,
    width = width,
    height = height,
    units = units,
    ext = "svg"
  )

  pdf_path <- save_review_figure(
    plot = plot,
    directory = directory,
    prefix = prefix,
    descriptor = descriptor,
    level = level,
    width = width,
    height = height,
    units = units,
    ext = "pdf"
  )

  invisible(c(svg = svg_path, pdf = pdf_path))
}

save_model_diagnostic_plots <- function(model, directory, prefix) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)

  aug <- broom.mixed::augment(model)

  p_resid_fitted <- ggplot(aug, aes(x = .fitted, y = .resid)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey45") +
    geom_point(size = 0.55, alpha = 0.35) +
    labs(x = "Fitted values", y = "Residuals", title = "Residuals vs fitted") +
    theme_publication()

  p_qq <- ggplot(aug, aes(sample = .resid)) +
    stat_qq(size = 0.55, alpha = 0.35) +
    stat_qq_line(linewidth = 0.25, colour = "black") +
    labs(x = "Theoretical quantiles", y = "Residual quantiles", title = "Residual Q-Q plot") +
    theme_publication()

  save_review_figure_dual(p_resid_fitted, directory, prefix, "diagnostic_residuals_vs_fitted", NULL, width = 80, height = 65)
  save_review_figure_dual(p_qq, directory, prefix, "diagnostic_qq_residuals", NULL, width = 80, height = 65)
}

cohens_d_from_summary <- function(mean_1, mean_2, sd_1, sd_2, n_1, n_2) {
  pooled_sd <- sqrt(((n_1 - 1) * sd_1^2 + (n_2 - 1) * sd_2^2) / (n_1 + n_2 - 2))
  (mean_1 - mean_2) / pooled_sd
}

make_pairwise_effect_sizes <- function(data, value_col, group_col = "Group", by_vars = c("Sex", "Phase")) {
  group_levels <- levels(data[[group_col]])
  comps <- combn(group_levels[!is.na(group_levels)], 2, simplify = FALSE)

  data %>%
    group_by(across(all_of(by_vars))) %>%
    group_modify(function(.x, .y) {
      purrr::map_dfr(comps, function(comp) {
        d1 <- .x %>% filter(.data[[group_col]] == comp[1]) %>% pull(.data[[value_col]])
        d2 <- .x %>% filter(.data[[group_col]] == comp[2]) %>% pull(.data[[value_col]])

        d1 <- d1[!is.na(d1)]
        d2 <- d2[!is.na(d2)]

        if (length(d1) < 2 || length(d2) < 2) {
          return(tibble(
            contrast = paste(comp[1], "-", comp[2]),
            n_1 = length(d1),
            n_2 = length(d2),
            mean_1 = mean(d1, na.rm = TRUE),
            mean_2 = mean(d2, na.rm = TRUE),
            mean_difference = mean(d1, na.rm = TRUE) - mean(d2, na.rm = TRUE),
            cohens_d = NA_real_
          ))
        }

        tibble(
          contrast = paste(comp[1], "-", comp[2]),
          n_1 = length(d1),
          n_2 = length(d2),
          mean_1 = mean(d1, na.rm = TRUE),
          mean_2 = mean(d2, na.rm = TRUE),
          sd_1 = sd(d1, na.rm = TRUE),
          sd_2 = sd(d2, na.rm = TRUE),
          mean_difference = mean(d1, na.rm = TRUE) - mean(d2, na.rm = TRUE),
          cohens_d = cohens_d_from_summary(
            mean(d1, na.rm = TRUE),
            mean(d2, na.rm = TRUE),
            sd(d1, na.rm = TRUE),
            sd(d2, na.rm = TRUE),
            length(d1),
            length(d2)
          )
        )
      })
    }) %>%
    ungroup() %>%
    mutate(across(where(is.numeric), ~ round(.x, 6)))
}

safe_factor <- function(x, levels) {
  factor(as.character(x), levels = levels)
}

trapezoid_auc <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 2) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

# ===================================================================
# 5. CUSTOM FUNCTIONS AND GROUP DEFINITIONS
# ===================================================================
custom_functions_path <- file.path(working_directory, "E9_SIS_AnimalPos-functions.R")
if (file.exists(custom_functions_path)) {
  source(custom_functions_path)
  log_message("Loaded custom AnimalPos functions")
} else {
  stop("Custom function file not found: ", custom_functions_path)
}

sus_animals <- readLines(file.path(raw_data_dir, "sus_animals.csv"))
con_animals <- readLines(file.path(raw_data_dir, "con_animals.csv"))

assign_group <- function(animal_id) {
  dplyr::case_when(
    animal_id %in% sus_animals ~ "sus",
    animal_id %in% con_animals ~ "con",
    TRUE ~ "res"
  )
}

batch_to_sex <- function(batch) {
  ifelse(batch %in% female_batches, "female", "male")
}

# ===================================================================
# 6. ENTROPY CALCULATION HELPERS
# ===================================================================
calc_entropy_safe <- function(prob_vec) {
  prob_vec <- prob_vec[!is.na(prob_vec) & !is.nan(prob_vec)]
  prob_vec <- prob_vec[prob_vec > 0]
  if (length(prob_vec) == 0 || sum(prob_vec) <= 0) return(NA_real_)
  prob_vec <- prob_vec / sum(prob_vec)
  -sum(prob_vec * log2(prob_vec))
}

init_cage_probability <- function() {
  replicate(8, c(NA, 0, 0, 0), simplify = FALSE) %>%
    purrr::imap(function(x, i) {
      x[1] <- i
      x
    })
}

init_animal_probability <- function(animal_ids) {
  tibble::tibble(
    AnimalID = rep(animal_ids, each = 8),
    Position = rep(1:8, length.out = length(animal_ids) * 8),
    Seconds = 0,
    SumPercentage = 0,
    Prob = 0
  )
}

compute_entropy_for_window <- function(data_window, animal_ids, system_complete) {
  animal_list <- list(
    animal_1 = list(name = "", time = "", position = 0),
    animal_2 = list(name = "", time = "", position = 0),
    animal_3 = list(name = "", time = "", position = 0),
    animal_4 = list(name = "", time = "", position = 0),
    data_temp = list(elapsed_seconds = 0, current_row = 0)
  )

  cage_position_probability <- init_cage_probability()
  animal_position_probability <- init_animal_probability(animal_ids)

  animal_list <- initialize_animal_positions(animal_ids, data_window, animal_list)
  initial_time <- animal_list[[1]][[2]]
  current_row <- 5
  total_rows <- nrow(data_window) + 1

  while (current_row < total_rows) {
    previous_animal_positions <- animal_list
    animal_list <- update_animal_list(
      animal_ids,
      animal_list,
      data_window,
      initial_time,
      current_row
    )

    elapsed_seconds <- animal_list[["data_temp"]][["elapsed_seconds"]]

    if (system_complete) {
      cage_position_probability <- update_cage_position_probability(
        previous_animal_positions,
        animal_list,
        cage_position_probability,
        elapsed_seconds
      )
    }

    animal_position_probability <- update_animal_position_probability(
      previous_animal_positions,
      animal_list,
      animal_position_probability,
      elapsed_seconds
    )

    current_row <- animal_list[["data_temp"]][["current_row"]]
    initial_time <- animal_list[[1]][[2]]
  }

  for (i in 1:8) {
    cage_position_probability[[i]][4] <- ifelse(
      cage_position_probability[[i]][2] > 0,
      cage_position_probability[[i]][3] / cage_position_probability[[i]][2],
      0
    )
  }

  animal_position_probability <- animal_position_probability %>%
    mutate(Prob = ifelse(Seconds > 0, SumPercentage / Seconds, 0))

  cage_prob_vec <- sapply(cage_position_probability, `[`, 4)
  cage_entropy <- if (system_complete) calc_entropy_safe(cage_prob_vec) else NA_real_

  animal_entropy_tbl <- animal_position_probability %>%
    group_by(AnimalID) %>%
    summarise(
      animalEntropy = calc_entropy_safe(Prob),
      .groups = "drop"
    )

  list(
    cage_prob = cage_position_probability,
    cage_entropy = cage_entropy,
    animal_entropy = animal_entropy_tbl,
    animal_position_probability = animal_position_probability
  )
}

# ===================================================================
# 7. CREATE BATCH/CAGE OUTPUT FOLDERS
# ===================================================================
create_batch_cage_folders <- function(batches, cageChanges) {
  for (batch in batches) {
    for (cageChange in cageChanges) {
      phase_tab <- file.path(batch_cage_tables_dir, batch, cageChange, "phase_based")
      halfhour_tab <- file.path(batch_cage_tables_dir, batch, cageChange, "halfhour_based")
      phase_plot <- file.path(batch_cage_plots_dir, batch, cageChange, "phase_based")
      halfhour_plot <- file.path(batch_cage_plots_dir, batch, cageChange, "halfhour_based")
      exploration_plot <- file.path(batch_cage_plots_dir, batch, cageChange, "exploration")
      for (d in c(phase_tab, halfhour_tab, phase_plot, halfhour_plot, exploration_plot)) {
        dir.create(d, recursive = TRUE, showWarnings = FALSE)
      }
    }
  }
}

create_batch_cage_folders(batches, cageChanges)

# ===================================================================
# 8. LOAD EXISTING DATA OR PROCESS RAW PREPROCESSED DATA
# ===================================================================
files_to_load <- list(
  cagePosProb = file.path(combined_tables_dir, "phase_based", "all_batch_CC_cagePosProb.csv"),
  cagePosEntropy = file.path(combined_tables_dir, "phase_based", "all_batch_CC_cagePosEntropy.csv"),
  animalPosEntropy = file.path(combined_tables_dir, "phase_based", "all_batch_CC_animalPosEntropy.csv"),
  cagePosEntropy_halfhour = file.path(combined_tables_dir, "halfhour_based", "all_batch_CC_cagePosEntropy_halfhour.csv"),
  animalPosEntropy_halfhour = file.path(combined_tables_dir, "halfhour_based", "all_batch_CC_animalPosEntropy_halfhour.csv")
)

if (load_existing_data && all(file.exists(unlist(files_to_load)))) {
  log_message("Loading existing processed entropy CSVs")

  all_cagePosProb <- readr::read_csv(files_to_load$cagePosProb, show_col_types = FALSE)
  all_cagePosEntropy <- readr::read_csv(files_to_load$cagePosEntropy, show_col_types = FALSE)
  all_animalPosEntropy <- readr::read_csv(files_to_load$animalPosEntropy, show_col_types = FALSE)
  all_cagePosEntropy_halfhour <- readr::read_csv(files_to_load$cagePosEntropy_halfhour, show_col_types = FALSE)
  all_animalPosEntropy_halfhour <- readr::read_csv(files_to_load$animalPosEntropy_halfhour, show_col_types = FALSE)

  log_message("Existing entropy CSVs loaded")

} else {
  log_message("Processing entropy from preprocessed position files")

  all_cagePosProb <- tibble(
    Batch = character(),
    System = character(),
    CageChange = character(),
    Phase = character(),
    Position = numeric(),
    Probability = numeric()
  )

  all_cagePosEntropy <- tibble(
    Batch = character(),
    Sex = character(),
    System = character(),
    CageChange = character(),
    Phase = character(),
    CageEntropy = numeric()
  )

  all_animalPosEntropy <- tibble(
    Batch = character(),
    Sex = character(),
    System = character(),
    CageChange = character(),
    Phase = character(),
    AnimalID = character(),
    animalEntropy = numeric(),
    Group = character()
  )

  all_cagePosEntropy_halfhour <- tibble(
    Batch = character(),
    Sex = character(),
    System = character(),
    CageChange = character(),
    HalfHour = numeric(),
    CageEntropy = numeric()
  )

  all_animalPosEntropy_halfhour <- tibble(
    Batch = character(),
    Sex = character(),
    System = character(),
    CageChange = character(),
    HalfHour = numeric(),
    AnimalID = character(),
    animalEntropy = numeric(),
    Group = character()
  )

  for (batch in batches) {
    sex <- batch_to_sex(batch)

    for (cageChange in cageChanges) {
      log_message("Processing batch/cage", batch = batch, cageChange = cageChange)

      filename <- paste0("E9_SIS_", batch, "_", cageChange, "_AnimalPos")
      csvFilePath <- file.path(preprocessed_data_dir, paste0(filename, "_preprocessed.csv"))

      if (!file.exists(csvFilePath)) {
        log_message(paste("File not found; skipping:", csvFilePath), batch, cageChange)
        next
      }

      data_preprocessed <- readr::read_csv(csvFilePath, show_col_types = FALSE) %>%
        as_tibble()

      if (!"Batch" %in% names(data_preprocessed)) data_preprocessed$Batch <- batch
      if (!"CageChange" %in% names(data_preprocessed)) data_preprocessed$CageChange <- cageChange

      unique_systems <- stringr::str_sort(unique(data_preprocessed$System))
      phases <- c("Active", "Inactive")
      active_phases <- sort(unique(data_preprocessed$ConsecActive))
      inactive_phases <- sort(unique(data_preprocessed$ConsecInactive))

      active_phases <- active_phases[!is.na(active_phases) & active_phases > 0]
      inactive_phases <- inactive_phases[!is.na(inactive_phases) & inactive_phases > 0]

      if (filter_cc4_late_phases && cageChange == "CC4") {
        active_phases <- active_phases[active_phases <= cc4_max_active_phase]
        inactive_phases <- inactive_phases[inactive_phases <= cc4_max_inactive_phase]
      }

      cagePosProb_list <- list()
      cagePosEntropy_list <- list()
      animalPosEntropy_list <- list()
      cageHalfhour_list <- list()
      animalHalfhour_list <- list()

      # ---------------------------------------------------------------
      # Phase-based entropy
      # ---------------------------------------------------------------
      for (system_id in unique_systems) {
        system_data <- data_preprocessed %>%
          filter(System == system_id) %>%
          as_tibble()

        animal_ids_raw <- unique(system_data$AnimalID)
        system_complete <- length(animal_ids_raw) >= 4
        animal_ids <- animal_ids_raw
        while (length(animal_ids) < 4) animal_ids <- append(animal_ids, NA)

        for (phase in phases) {
          phase_nums <- if (phase == "Active") active_phases else inactive_phases
          phase_code <- ifelse(phase == "Active", "A", "I")

          for (phaseN in phase_nums) {
            data_phase <- system_data %>%
              filter(
                ConsecActive == ifelse(phase == "Active", phaseN, 0),
                ConsecInactive == ifelse(phase == "Inactive", phaseN, 0)
              ) %>%
              as_tibble()

            if (nrow(data_phase) == 0) next

            entropy_result <- compute_entropy_for_window(data_phase, animal_ids, system_complete)
            phase_label <- paste0(phase_code, phaseN)

            if (system_complete) {
              cage_prob_vec <- sapply(entropy_result$cage_prob, `[`, 4)
              cagePosProb_list[[length(cagePosProb_list) + 1]] <- tibble(
                Batch = batch,
                System = system_id,
                CageChange = cageChange,
                Phase = phase_label,
                Position = 1:8,
                Probability = cage_prob_vec
              )

              cagePosEntropy_list[[length(cagePosEntropy_list) + 1]] <- tibble(
                Batch = batch,
                Sex = sex,
                System = system_id,
                CageChange = cageChange,
                Phase = phase_label,
                CageEntropy = entropy_result$cage_entropy
              )
            }

            animal_tbl <- entropy_result$animal_entropy %>%
              filter(!is.na(AnimalID)) %>%
              mutate(
                Batch = batch,
                Sex = sex,
                System = system_id,
                CageChange = cageChange,
                Phase = phase_label,
                Group = assign_group(AnimalID)
              ) %>%
              select(Batch, Sex, System, CageChange, Phase, AnimalID, animalEntropy, Group)

            animalPosEntropy_list[[length(animalPosEntropy_list) + 1]] <- animal_tbl
          }
        }
      }

      # ---------------------------------------------------------------
      # Half-hour entropy
      # ---------------------------------------------------------------
      if (analyze_by_halfhour) {
        for (system_id in unique_systems) {
          system_data <- data_preprocessed %>%
            filter(System == system_id) %>%
            as_tibble()

          animal_ids_raw <- unique(system_data$AnimalID)
          system_complete <- length(animal_ids_raw) >= 4
          animal_ids <- animal_ids_raw
          while (length(animal_ids) < 4) animal_ids <- append(animal_ids, NA)

          valid_halfhours <- sort(unique(system_data$HalfHoursElapsed))
          valid_halfhours <- valid_halfhours[!is.na(valid_halfhours)]

          if (filter_cc4_late_phases && cageChange == "CC4") {
            valid_halfhours <- system_data %>%
              filter(
                (ConsecActive > 0 & ConsecActive <= cc4_max_active_phase) |
                  (ConsecInactive > 0 & ConsecInactive <= cc4_max_inactive_phase)
              ) %>%
              pull(HalfHoursElapsed) %>%
              unique() %>%
              sort()
          }

          for (halfhour in valid_halfhours) {
            data_halfhour <- system_data %>%
              filter(HalfHoursElapsed == halfhour) %>%
              as_tibble()

            if (nrow(data_halfhour) == 0) next

            entropy_result <- compute_entropy_for_window(data_halfhour, animal_ids, system_complete)

            if (system_complete) {
              cageHalfhour_list[[length(cageHalfhour_list) + 1]] <- tibble(
                Batch = batch,
                Sex = sex,
                System = system_id,
                CageChange = cageChange,
                HalfHour = halfhour,
                CageEntropy = entropy_result$cage_entropy
              )
            }

            animal_tbl <- entropy_result$animal_entropy %>%
              filter(!is.na(AnimalID)) %>%
              mutate(
                Batch = batch,
                Sex = sex,
                System = system_id,
                CageChange = cageChange,
                HalfHour = halfhour,
                Group = assign_group(AnimalID)
              ) %>%
              select(Batch, Sex, System, CageChange, HalfHour, AnimalID, animalEntropy, Group)

            animalHalfhour_list[[length(animalHalfhour_list) + 1]] <- animal_tbl
          }
        }
      }

      cagePosProb <- bind_rows(cagePosProb_list)
      cagePosEntropy <- bind_rows(cagePosEntropy_list)
      animalPosEntropy <- bind_rows(animalPosEntropy_list)
      cagePosEntropy_halfhour <- bind_rows(cageHalfhour_list)
      animalPosEntropy_halfhour <- bind_rows(animalHalfhour_list)

      # Batch-level save
      current_tables_dir_phase <- file.path(batch_cage_tables_dir, batch, cageChange, "phase_based")
      current_tables_dir_halfhour <- file.path(batch_cage_tables_dir, batch, cageChange, "halfhour_based")

      if (save_tables) {
        readr::write_csv(cagePosProb, file.path(current_tables_dir_phase, paste0(batch, "_", cageChange, "_cagePosProb.csv")))
        readr::write_csv(cagePosEntropy, file.path(current_tables_dir_phase, paste0(batch, "_", cageChange, "_cagePosEntropy.csv")))
        readr::write_csv(animalPosEntropy, file.path(current_tables_dir_phase, paste0(batch, "_", cageChange, "_animalPosEntropy.csv")))

        if (analyze_by_halfhour) {
          readr::write_csv(cagePosEntropy_halfhour, file.path(current_tables_dir_halfhour, paste0(batch, "_", cageChange, "_cagePosEntropy_halfhour.csv")))
          readr::write_csv(animalPosEntropy_halfhour, file.path(current_tables_dir_halfhour, paste0(batch, "_", cageChange, "_animalPosEntropy_halfhour.csv")))
        }
      }

      all_cagePosProb <- bind_rows(all_cagePosProb, cagePosProb)
      all_cagePosEntropy <- bind_rows(all_cagePosEntropy, cagePosEntropy)
      all_animalPosEntropy <- bind_rows(all_animalPosEntropy, animalPosEntropy)
      all_cagePosEntropy_halfhour <- bind_rows(all_cagePosEntropy_halfhour, cagePosEntropy_halfhour)
      all_animalPosEntropy_halfhour <- bind_rows(all_animalPosEntropy_halfhour, animalPosEntropy_halfhour)

      log_message("Finished batch/cage", batch = batch, cageChange = cageChange)
    }
  }

  if (save_tables) {
    readr::write_csv(all_cagePosProb, files_to_load$cagePosProb)
    readr::write_csv(all_cagePosEntropy, files_to_load$cagePosEntropy)
    readr::write_csv(all_animalPosEntropy, files_to_load$animalPosEntropy)
    readr::write_csv(all_cagePosEntropy_halfhour, files_to_load$cagePosEntropy_halfhour)
    readr::write_csv(all_animalPosEntropy_halfhour, files_to_load$animalPosEntropy_halfhour)
  }
}

# ===================================================================
# 9. CLEAN PHASE-BASED ENTROPY TABLES
# ===================================================================
standardize_group_levels <- function(data) {
  data %>%
    mutate(
      Group = if ("Group" %in% names(.)) as.character(Group) else assign_group(AnimalID),
      Group = str_to_lower(str_trim(Group)),
      Sex = str_to_lower(str_trim(as.character(Sex))),
      CageChange = str_to_upper(str_trim(as.character(CageChange))),
      Batch = as.character(Batch),
      System = as.character(System),
      Group = factor(Group, levels = c("con", "res", "sus")),
      Sex = factor(Sex, levels = c("male", "female")),
      CageChange = factor(CageChange, levels = c("CC1", "CC2", "CC3", "CC4")),
      Batch = factor(Batch),
      System = factor(System)
    )
}

standardize_cage_levels <- function(data) {
  data %>%
    mutate(
      Sex = str_to_lower(str_trim(as.character(Sex))),
      CageChange = str_to_upper(str_trim(as.character(CageChange))),
      Batch = as.character(Batch),
      System = as.character(System),
      Sex = factor(Sex, levels = c("male", "female")),
      CageChange = factor(CageChange, levels = c("CC1", "CC2", "CC3", "CC4")),
      Batch = factor(Batch),
      System = factor(System)
    )
}

parse_phase_table <- function(data, entropy_col) {
  data %>%
    tidyr::separate(
      Phase,
      into = c("PhaseCode", "PhaseNumber"),
      sep = 1,
      remove = FALSE,
      convert = FALSE
    ) %>%
    mutate(
      PhaseNumber = suppressWarnings(as.numeric(PhaseNumber)),
      Phase = case_when(
        PhaseCode == "A" ~ "active",
        PhaseCode == "I" ~ "inactive",
        TRUE ~ NA_character_
      ),
      Phase = factor(Phase, levels = c("active", "inactive")),
      ConsecActive = ifelse(Phase == "active", PhaseNumber, 0),
      ConsecInactive = ifelse(Phase == "inactive", PhaseNumber - 1, 0),
      EntropyNorm = .data[[entropy_col]] / max_position_entropy
    )
}

animal_phase <- all_animalPosEntropy %>%
  standardize_group_levels() %>%
  parse_phase_table("animalEntropy") %>%
  filter(!is.na(animalEntropy), !is.na(Group), !is.na(Sex), !is.na(Phase))

cage_phase <- all_cagePosEntropy %>%
  standardize_cage_levels() %>%
  parse_phase_table("CageEntropy") %>%
  filter(!is.na(CageEntropy), !is.na(Sex), !is.na(Phase))

if (filter_cc4_late_phases) {
  animal_phase <- animal_phase %>%
    filter(!(CageChange == "CC4" & Phase == "active" & ConsecActive > cc4_max_active_phase)) %>%
    filter(!(CageChange == "CC4" & Phase == "inactive" & ConsecInactive > cc4_max_inactive_phase))

  cage_phase <- cage_phase %>%
    filter(!(CageChange == "CC4" & Phase == "active" & ConsecActive > cc4_max_active_phase)) %>%
    filter(!(CageChange == "CC4" & Phase == "inactive" & ConsecInactive > cc4_max_inactive_phase))
}

# Consecutive phase counts across cage changes
add_consecutive_phase_counts <- function(data) {
  out <- data
  max_consecAct <- 0
  max_consecInact <- 0

  out$ConsecActiveGlobal <- out$ConsecActive
  out$ConsecInactiveGlobal <- out$ConsecInactive

  for (change in levels(out$CageChange)) {
    if (change != "CC1") {
      out <- out %>%
        mutate(
          ConsecActiveGlobal = ifelse(CageChange == change & Phase == "active", ConsecActive + max_consecAct, ConsecActiveGlobal),
          ConsecInactiveGlobal = ifelse(CageChange == change & Phase == "inactive", ConsecInactive + max_consecInact, ConsecInactiveGlobal)
        )
    }

    current_act <- out %>%
      filter(CageChange == change) %>%
      pull(ConsecActiveGlobal)
    current_inact <- out %>%
      filter(CageChange == change) %>%
      pull(ConsecInactiveGlobal)

    if (length(current_act) > 0 && any(is.finite(current_act))) max_consecAct <- max(current_act, na.rm = TRUE)
    if (length(current_inact) > 0 && any(is.finite(current_inact))) max_consecInact <- max(current_inact, na.rm = TRUE)
  }

  out
}

animal_phase <- add_consecutive_phase_counts(animal_phase)
cage_phase <- add_consecutive_phase_counts(cage_phase)

# ===================================================================
# 10. CLEAN HALF-HOUR TABLES AND ALIGN PHASES
# ===================================================================
load_all_preprocessed_for_lookup <- function() {
  all_preprocessed <- list()

  for (batch in batches) {
    for (cageChange in cageChanges) {
      filename <- paste0("E9_SIS_", batch, "_", cageChange, "_AnimalPos")
      csvFilePath <- file.path(preprocessed_data_dir, paste0(filename, "_preprocessed.csv"))

      if (!file.exists(csvFilePath)) next

      temp_data <- readr::read_csv(csvFilePath, show_col_types = FALSE) %>%
        as_tibble() %>%
        mutate(
          Batch = if ("Batch" %in% colnames(.)) Batch else batch,
          CageChange = if ("CageChange" %in% colnames(.)) CageChange else cageChange
        )

      all_preprocessed[[paste(batch, cageChange, sep = "_")]] <- temp_data
    }
  }

  bind_rows(all_preprocessed)
}

make_phase_lookup <- function(preprocessed_data) {
  preprocessed_data %>%
    mutate(
      Batch = as.character(Batch),
      System = as.character(System),
      CageChange = str_to_upper(str_trim(as.character(CageChange))),
      HalfHour = HalfHoursElapsed
    ) %>%
    group_by(Batch, System, CageChange, HalfHour) %>%
    summarise(
      ConsecActive = if (all(is.na(ConsecActive))) NA_real_ else max(ConsecActive, na.rm = TRUE),
      ConsecInactive = if (all(is.na(ConsecInactive))) NA_real_ else max(ConsecInactive, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Phase = case_when(
        !is.na(ConsecActive) & ConsecActive > 0 ~ "active",
        !is.na(ConsecInactive) & ConsecInactive > 0 ~ "inactive",
        TRUE ~ NA_character_
      )
    )
}

preprocessed_lookup_data <- load_all_preprocessed_for_lookup()
phase_lookup <- if (nrow(preprocessed_lookup_data) > 0) make_phase_lookup(preprocessed_lookup_data) else NULL

make_halfhour_clean <- function(data, entropy_col, animal_level = TRUE) {
  cleaned <- data %>%
    mutate(
      Batch = as.character(Batch),
      System = as.character(System),
      CageChange = str_to_upper(str_trim(as.character(CageChange))),
      Sex = str_to_lower(str_trim(as.character(Sex))),
      HalfHour = if ("HalfHour" %in% names(.)) HalfHour else HalfHoursElapsed,
      Sex = factor(Sex, levels = c("male", "female")),
      CageChange = factor(CageChange, levels = c("CC1", "CC2", "CC3", "CC4")),
      Batch = factor(Batch),
      System = factor(System),
      EntropyNorm = .data[[entropy_col]] / max_position_entropy
    )

  if (animal_level) {
    cleaned <- cleaned %>%
      mutate(
        Group = if ("Group" %in% names(.)) as.character(Group) else assign_group(AnimalID),
        Group = factor(str_to_lower(str_trim(Group)), levels = c("con", "res", "sus"))
      )
  }

  if (!is.null(phase_lookup)) {
    cleaned <- cleaned %>%
      left_join(
        phase_lookup,
        by = c("Batch", "System", "CageChange", "HalfHour")
      )
  } else {
    cleaned <- cleaned %>%
      mutate(
        ConsecActive = NA_real_,
        ConsecInactive = NA_real_,
        Phase = NA_character_
      )
  }

  cleaned <- cleaned %>%
    group_by(Batch, System) %>%
    arrange(CageChange, HalfHour, .by_group = TRUE) %>%
    mutate(
      ConsecHalfHour = dense_rank(interaction(CageChange, HalfHour, drop = TRUE)),
      HalfHourWithinCC = HalfHour,
      HalfHourWithinCCOffset = as.integer(
        !any(HalfHourWithinCC == 0, na.rm = TRUE) &&
          is.finite(min(HalfHourWithinCC, na.rm = TRUE)) &&
          min(HalfHourWithinCC, na.rm = TRUE) == 1
      ),
      HalfHourWithinCC0 = HalfHourWithinCC - HalfHourWithinCCOffset,
      HalfHourWithinCC0 = as.integer(HalfHourWithinCC0),
      HalfHourInPhase = (HalfHourWithinCC0 %% 24) + 1,
      PhaseBlockInferred = floor(HalfHourWithinCC0 / 24) + 1,
      PhaseInferred = ifelse(PhaseBlockInferred %% 2 == 1, "active", "inactive"),
      Phase = ifelse(is.na(Phase), PhaseInferred, Phase),
      Phase = factor(Phase, levels = c("active", "inactive")),
      HalfHourInPhase_z = as.numeric(scale(HalfHourInPhase)),
      ConsecHalfHour_z = as.numeric(scale(ConsecHalfHour))
    ) %>%
    ungroup()

  if (filter_cc4_late_phases) {
    cleaned <- cleaned %>%
      filter(!(CageChange == "CC4" & Phase == "active" & !is.na(ConsecActive) & ConsecActive > cc4_max_active_phase)) %>%
      filter(!(CageChange == "CC4" & Phase == "inactive" & !is.na(ConsecInactive) & ConsecInactive > cc4_max_inactive_phase))
  }

  cleaned
}

animal_halfhour <- make_halfhour_clean(all_animalPosEntropy_halfhour, "animalEntropy", animal_level = TRUE) %>%
  filter(!is.na(animalEntropy), !is.na(Group), !is.na(Sex), !is.na(Phase), !is.na(ConsecHalfHour))

cage_halfhour <- make_halfhour_clean(all_cagePosEntropy_halfhour, "CageEntropy", animal_level = FALSE) %>%
  filter(!is.na(CageEntropy), !is.na(Sex), !is.na(Phase), !is.na(ConsecHalfHour))

if (save_tables) {
  write_review_csv(animal_phase, review_clean_data_dir, "clean", "animal_entropy", "phase")
  write_review_csv(cage_phase, review_clean_data_dir, "clean", "cage_entropy", "phase")
  write_review_csv(animal_halfhour, review_clean_data_dir, "clean", "animal_entropy", "halfhour")
  write_review_csv(cage_halfhour, review_clean_data_dir, "clean", "cage_entropy", "halfhour")
}

# ===================================================================
# 11. DESCRIPTIVE TABLES
# ===================================================================
make_entropy_summary_table <- function(data, entropy_col, grouping_vars, animal_col = NULL) {
  summary_tbl <- data %>%
    group_by(across(all_of(grouping_vars))) %>%
    summarise(
      n_observations = sum(!is.na(.data[[entropy_col]])),
      mean = mean(.data[[entropy_col]], na.rm = TRUE),
      sd = sd(.data[[entropy_col]], na.rm = TRUE),
      sem = sem(.data[[entropy_col]]),
      median = median(.data[[entropy_col]], na.rm = TRUE),
      q25 = quantile(.data[[entropy_col]], 0.25, na.rm = TRUE),
      q75 = quantile(.data[[entropy_col]], 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  if (!is.null(animal_col) && animal_col %in% names(data)) {
    animal_n <- data %>%
      group_by(across(all_of(grouping_vars))) %>%
      summarise(n_animals = n_distinct(.data[[animal_col]]), .groups = "drop")

    summary_tbl <- summary_tbl %>%
      left_join(animal_n, by = grouping_vars) %>%
      relocate(n_animals, .after = n_observations)
  }

  summary_tbl %>%
    mutate(
      mean_sem = sprintf("%.3f ± %.3f", mean, sem),
      median_iqr = sprintf("%.3f [%.3f, %.3f]", median, q25, q75)
    )
}

summary_animal_phase <- make_entropy_summary_table(
  animal_phase,
  entropy_col = "animalEntropy",
  grouping_vars = c("Group", "Sex", "Phase", "CageChange"),
  animal_col = "AnimalID"
)

summary_cage_phase <- make_entropy_summary_table(
  cage_phase,
  entropy_col = "CageEntropy",
  grouping_vars = c("Sex", "Phase", "CageChange")
)

summary_animal_halfhour <- animal_halfhour %>%
  group_by(Group, Sex, Phase, ConsecHalfHour) %>%
  summarise(
    n_observations = n(),
    n_animals = n_distinct(AnimalID),
    mean_entropy = mean(animalEntropy, na.rm = TRUE),
    sem_entropy = sem(animalEntropy),
    .groups = "drop"
  )

summary_cage_halfhour <- cage_halfhour %>%
  group_by(Sex, Phase, ConsecHalfHour) %>%
  summarise(
    n_observations = n(),
    n_systems = n_distinct(System),
    mean_entropy = mean(CageEntropy, na.rm = TRUE),
    sem_entropy = sem(CageEntropy),
    .groups = "drop"
  )

if (save_tables) {
  write_review_csv(summary_animal_phase, review_tables_dir, "summary", "animal_entropy", "phase")
  write_review_csv(summary_cage_phase, review_tables_dir, "summary", "cage_entropy", "phase")
  write_review_csv(summary_animal_halfhour, review_tables_dir, "summary", "animal_entropy", "halfhour")
  write_review_csv(summary_cage_halfhour, review_tables_dir, "summary", "cage_entropy", "halfhour")
}

# ===================================================================
# 12. PRIMARY MIXED MODEL: INDIVIDUAL ENTROPY
# ===================================================================
format_lmm_fixed_effects <- function(model, model_name) {
  broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
    mutate(
      model = model_name,
      across(where(is.numeric), ~ round(.x, 6))
    ) %>%
    select(model, everything())
}

format_type3_anova <- function(model, model_name) {
  as.data.frame(anova(model, type = 3)) %>%
    rownames_to_column("term") %>%
    mutate(
      model = model_name,
      across(where(is.numeric), ~ round(.x, 6))
    ) %>%
    select(model, everything())
}

primary_lmm <- NULL
primary_lmm_norm <- NULL

if (RUN_PRIMARY_LMM) {
  log_message("Fitting primary LMM for individual entropy")

  primary_lmm <- lmerTest::lmer(
    animalEntropy ~ Sex * Group * Phase + CageChange +
      (1 | Batch/System) +
      (1 | AnimalID),
    data = animal_phase,
    REML = TRUE
  )

  primary_lmm_norm <- lmerTest::lmer(
    EntropyNorm ~ Sex * Group * Phase + CageChange +
      (1 | Batch/System) +
      (1 | AnimalID),
    data = animal_phase,
    REML = TRUE
  )

  primary_lmm_anova <- format_type3_anova(primary_lmm, "primary_lmm_animal_entropy")
  primary_lmm_fixed <- format_lmm_fixed_effects(primary_lmm, "primary_lmm_animal_entropy")

  primary_emm_group <- emmeans::emmeans(primary_lmm, ~ Group | Sex * Phase)
  primary_group_contrasts <- emmeans::contrast(primary_emm_group, method = "pairwise", adjust = "holm") %>%
    as.data.frame(infer = TRUE) %>%
    mutate(across(where(is.numeric), ~ round(.x, 6)))

  primary_emm_sex <- emmeans::emmeans(primary_lmm, ~ Sex | Group * Phase)
  primary_sex_contrasts <- emmeans::contrast(primary_emm_sex, method = "pairwise", adjust = "holm") %>%
    as.data.frame(infer = TRUE) %>%
    mutate(across(where(is.numeric), ~ round(.x, 6)))

  primary_effect_sizes <- make_pairwise_effect_sizes(
    animal_phase,
    value_col = "animalEntropy",
    group_col = "Group",
    by_vars = c("Sex", "Phase")
  )

  if (save_tables) {
    write_review_csv(primary_lmm_anova, review_lmm_stats_dir, "primary_lmm", "type3_anova", "animal_entropy")
    write_review_csv(primary_lmm_fixed, review_lmm_stats_dir, "primary_lmm", "fixed_effects", "animal_entropy")
    write_review_csv(as.data.frame(primary_emm_group), review_emm_stats_dir, "primary_lmm", "emmeans_group_by_sex_phase", "animal_entropy")
    write_review_csv(primary_group_contrasts, review_emm_stats_dir, "primary_lmm", "group_contrasts_holm", "animal_entropy")
    write_review_csv(as.data.frame(primary_emm_sex), review_emm_stats_dir, "primary_lmm", "emmeans_sex_by_group_phase", "animal_entropy")
    write_review_csv(primary_sex_contrasts, review_emm_stats_dir, "primary_lmm", "sex_contrasts_holm", "animal_entropy")
    write_review_csv(primary_effect_sizes, review_effect_stats_dir, "effect_size", "cohens_d_pairwise", "animal_entropy")

    write_review_csv(primary_lmm_anova, review_main_tables_dir, "table", "primary_lmm_type3_anova", "animal_entropy")
    write_review_csv(primary_group_contrasts, review_main_tables_dir, "table", "primary_group_contrasts_holm", "animal_entropy")
    write_review_csv(primary_effect_sizes, review_supp_tables_dir, "table_s", "effect_sizes_pairwise", "animal_entropy")
  }

  saveRDS(primary_lmm, file.path(review_model_dir, "primary_lmm_animal_entropy.rds"))
  saveRDS(primary_lmm_norm, file.path(review_model_dir, "primary_lmm_animal_entropy_normalized.rds"))

  if (save_plots) {
    save_model_diagnostic_plots(
      primary_lmm,
      directory = review_lmm_diagnostics_dir,
      prefix = "primary_lmm_animal_entropy"
    )
  }

  sink(file.path(review_diagnostics_dir, "primary_lmm_animal_entropy_diagnostics.txt"))
  print(summary(primary_lmm))
  print(anova(primary_lmm, type = 3))
  print(performance::check_singularity(primary_lmm))
  print(performance::check_collinearity(primary_lmm))
  print(performance::icc(primary_lmm))
  print(performance::r2(primary_lmm))
  sink()
}

# ===================================================================
# 13. SECONDARY MIXED MODEL: CAGE ENTROPY
# ===================================================================
cage_lmm <- NULL

if (RUN_CAGE_LMM) {
  log_message("Fitting secondary LMM for cage entropy")

  cage_lmm <- lmerTest::lmer(
    CageEntropy ~ Sex * Phase + CageChange +
      (1 | Batch/System),
    data = cage_phase,
    REML = TRUE
  )

  cage_lmm_anova <- format_type3_anova(cage_lmm, "secondary_lmm_cage_entropy")
  cage_lmm_fixed <- format_lmm_fixed_effects(cage_lmm, "secondary_lmm_cage_entropy")

  cage_emm_phase <- emmeans::emmeans(cage_lmm, ~ Phase | Sex)
  cage_phase_contrasts <- emmeans::contrast(cage_emm_phase, method = "pairwise", adjust = "holm") %>%
    as.data.frame(infer = TRUE) %>%
    mutate(across(where(is.numeric), ~ round(.x, 6)))

  if (save_tables) {
    write_review_csv(cage_lmm_anova, review_tables_dir, "secondary_lmm", "type3_anova", "cage_entropy")
    write_review_csv(cage_lmm_fixed, review_tables_dir, "secondary_lmm", "fixed_effects", "cage_entropy")
    write_review_csv(as.data.frame(cage_emm_phase), review_tables_dir, "secondary_lmm", "emmeans_phase_by_sex", "cage_entropy")
    write_review_csv(cage_phase_contrasts, review_tables_dir, "secondary_lmm", "phase_contrasts_holm", "cage_entropy")
  }

  saveRDS(cage_lmm, file.path(review_model_dir, "secondary_lmm_cage_entropy.rds"))

  if (save_plots) {
    save_model_diagnostic_plots(
      cage_lmm,
      directory = review_lmm_diagnostics_dir,
      prefix = "secondary_lmm_cage_entropy"
    )
  }

  sink(file.path(review_diagnostics_dir, "secondary_lmm_cage_entropy_diagnostics.txt"))
  print(summary(cage_lmm))
  print(anova(cage_lmm, type = 3))
  print(performance::check_singularity(cage_lmm))
  print(performance::check_collinearity(cage_lmm))
  print(performance::icc(cage_lmm))
  print(performance::r2(cage_lmm))
  sink()
}

# ===================================================================
# 14. SECONDARY TIME-RESOLVED GAMM
# ===================================================================
estimate_rho_from_series <- function(data, response_col = "animalEntropy") {
  rho_tbl <- data %>%
    arrange(AnimalID, ConsecHalfHour) %>%
    group_by(AnimalID) %>%
    summarise(
      rho1 = {
        x <- .data[[response_col]]
        x <- x[!is.na(x)]
        if (length(x) > 3) {
          as.numeric(stats::acf(x, plot = FALSE, lag.max = 1)$acf[2])
        } else {
          NA_real_
        }
      },
      .groups = "drop"
    )

  rho_est <- median(rho_tbl$rho1, na.rm = TRUE)
  if (!is.finite(rho_est)) rho_est <- 0
  rho_est <- max(GAMM_RHO_MIN, min(GAMM_RHO_MAX, rho_est))
  rho_est
}

fit_gamm_by_sex <- function(data, sex_value) {
  data_sex <- data %>%
    filter(Sex == sex_value) %>%
    arrange(Batch, System, AnimalID, CageChange, ConsecHalfHour) %>%
    mutate(
      AnimalID = factor(AnimalID),
      System = factor(System),
      Batch = factor(Batch),
      Group = factor(Group, levels = c("con", "res", "sus")),
      Phase = factor(Phase, levels = c("active", "inactive")),
      AR_start = row_number() == 1L |
        AnimalID != lag(AnimalID) |
        System != lag(System) |
        CageChange != lag(CageChange)
    ) %>%
    filter(
      !is.na(animalEntropy),
      !is.na(Group),
      !is.na(Phase),
      !is.na(CageChange),
      is.finite(ConsecHalfHour_z),
      is.finite(HalfHourInPhase_z)
    ) %>%
    droplevels()

  log_message(paste0(
    "GAMM input for ", sex_value,
    ": rows = ", nrow(data_sex),
    "; groups = ", paste(levels(data_sex$Group), collapse = ","),
    "; phases = ", paste(levels(data_sex$Phase), collapse = ","),
    "; cage changes = ", paste(levels(data_sex$CageChange), collapse = ",")
  ))

  if (nrow(data_sex) < 50) {
    log_message(paste("Skipping GAMM for", sex_value, ": fewer than 50 rows"))
    return(NULL)
  }

  if (n_distinct(data_sex$Group) < 2) {
    log_message(paste("Skipping GAMM for", sex_value, ": fewer than two groups"))
    return(NULL)
  }

  rho_est <- if (GAMM_ESTIMATE_RHO) estimate_rho_from_series(data_sex) else 0
  log_message(paste0("GAMM rho for ", sex_value, " = ", round(rho_est, 3)))

  gamm_fit <- tryCatch(
    mgcv::bam(
      animalEntropy ~ Group * Phase + CageChange +
        s(ConsecHalfHour_z, k = GAMM_K_GLOBAL) +
        s(ConsecHalfHour_z, by = Group, k = GAMM_K_GROUP_DEVIATION) +
        s(HalfHourInPhase_z, k = GAMM_K_HALFHOUR_PHASE) +
        s(AnimalID, bs = "re") +
        s(System, bs = "re"),
      data = data_sex,
      method = "fREML",
      discrete = TRUE,
      gamma = GAMM_GAMMA,
      rho = rho_est,
      AR.start = data_sex$AR_start
    ),
    error = function(e) {
      log_message(paste("Skipping GAMM for", sex_value, ":", conditionMessage(e)))
      return(NULL)
    }
  )

  if (is.null(gamm_fit)) return(NULL)

  list(
    model = gamm_fit,
    data = data_sex,
    rho = rho_est
  )
}

gamm_results <- list()

if (RUN_GAMM && analyze_by_halfhour) {
  log_message("Fitting conservative sex-stratified GAMMs")

  for (sex_value in c("male", "female")) {
    gamm_results[[sex_value]] <- fit_gamm_by_sex(animal_halfhour, sex_value)
  }

  for (sex_value in names(gamm_results)) {
    res <- gamm_results[[sex_value]]
    if (is.null(res)) next

    gamm_fit <- res$model

    gamm_parametric <- as.data.frame(summary(gamm_fit)$p.table) %>%
      rownames_to_column("term") %>%
      mutate(
        Sex = sex_value,
        rho = res$rho,
        across(where(is.numeric), ~ round(.x, 6))
      )

    gamm_smooths <- as.data.frame(summary(gamm_fit)$s.table) %>%
      rownames_to_column("smooth_term") %>%
      mutate(
        Sex = sex_value,
        rho = res$rho,
        across(where(is.numeric), ~ round(.x, 6))
      )

    if (save_tables) {
      write_review_csv(gamm_parametric, review_tables_dir, "secondary_gamm", "parametric", sex_value)
      write_review_csv(gamm_smooths, review_tables_dir, "secondary_gamm", "smooths", sex_value)
    }

    saveRDS(gamm_fit, file.path(review_model_dir, paste0("secondary_gamm_animal_entropy_", sex_value, ".rds")))

    sink(file.path(review_diagnostics_dir, paste0("secondary_gamm_diagnostics_", sex_value, ".txt")))
    print(summary(gamm_fit))
    print(mgcv::gam.check(gamm_fit))
    print(mgcv::concurvity(gamm_fit, full = TRUE))
    sink()
  }
}

prepare_cc1_active_phase1_data <- function(data) {
  data_cc1 <- data %>%
    mutate(
      CageChange_chr = str_to_upper(str_trim(as.character(CageChange))),
      Phase_chr = str_to_lower(str_trim(as.character(Phase))),
      ConsecActive_num = suppressWarnings(as.numeric(as.character(ConsecActive))),
      HalfHourWithinCC_raw = suppressWarnings(as.numeric(as.character(HalfHourWithinCC))),
      HalfHourWithinCCOffset = as.integer(
        !any(HalfHourWithinCC_raw == 0, na.rm = TRUE) &&
          is.finite(min(HalfHourWithinCC_raw, na.rm = TRUE)) &&
          min(HalfHourWithinCC_raw, na.rm = TRUE) == 1
      ),
      HalfHourWithinCC0 = HalfHourWithinCC_raw - HalfHourWithinCCOffset,
      HalfHourWithinCC0 = as.integer(HalfHourWithinCC0)
    ) %>%
    filter(
      CageChange_chr == "CC1",
      Phase_chr == "active",
      !is.na(HalfHourWithinCC0),
      HalfHourWithinCC0 >= 0,
      HalfHourWithinCC0 <= 23
    )

  if (nrow(data_cc1) == 0) return(data_cc1)

  first_active_rows <- data_cc1 %>%
    filter(is.na(ConsecActive_num) | ConsecActive_num == 1)

  if (n_distinct(first_active_rows$HalfHourWithinCC0) >= 5) {
    out <- first_active_rows
  } else {
    log_message("CC1 active-phase-1 GAMM using HalfHourWithinCC0 0-23 window because ConsecActive == 1 filter did not retain enough bins")
    out <- data_cc1
  }

  log_message(paste0(
    "CC1 active-phase-1 candidate rows = ", nrow(out),
    "; half-hour bins = ", paste(sort(unique(out$HalfHourWithinCC0)), collapse = ","),
    "; groups = ", paste(sort(unique(as.character(out$Group))), collapse = ",")
  ))

  out %>%
    arrange(Batch, System, AnimalID, HalfHourWithinCC0) %>%
    mutate(
      CC1ActivePhase1Definition = "CC1 active phase, HalfHourWithinCC0 0-23; ConsecActive == 1 where available",
      AnimalID = factor(AnimalID),
      System = factor(System),
      Batch = factor(Batch),
      Sex = droplevels(factor(Sex, levels = c("male", "female"))),
      Group = droplevels(factor(Group, levels = c("con", "res", "sus"))),
      AR_start = row_number() == 1L |
        AnimalID != lag(AnimalID) |
        System != lag(System)
    )
}

fit_cc1_active_phase1_gamm_by_sex <- function(data, sex_value) {
  model_data <- prepare_cc1_active_phase1_data(data) %>%
    filter(Sex == sex_value) %>%
    droplevels()

  cc1_precheck <- data %>%
    mutate(
      CageChange_chr = str_to_upper(str_trim(as.character(CageChange))),
      Phase_chr = str_to_lower(str_trim(as.character(Phase))),
      ConsecActive_num = suppressWarnings(as.numeric(as.character(ConsecActive))),
      HalfHourWithinCC_raw = suppressWarnings(as.numeric(as.character(HalfHourWithinCC))),
      HalfHourWithinCCOffset = as.integer(
        !any(HalfHourWithinCC_raw == 0, na.rm = TRUE) &&
          is.finite(min(HalfHourWithinCC_raw, na.rm = TRUE)) &&
          min(HalfHourWithinCC_raw, na.rm = TRUE) == 1
      ),
      HalfHourWithinCC0 = as.integer(HalfHourWithinCC_raw - HalfHourWithinCCOffset),
      is_cc1_active_window = CageChange_chr == "CC1" &
        Phase_chr == "active" &
        !is.na(HalfHourWithinCC0) &
        HalfHourWithinCC0 >= 0 &
        HalfHourWithinCC0 <= 23,
      is_cc1_active_consec1_window = is_cc1_active_window &
        (is.na(ConsecActive_num) | ConsecActive_num == 1)
    ) %>%
    summarise(
      Sex = sex_value,
      n_cc1_active_window_rows = sum(is_cc1_active_window, na.rm = TRUE),
      n_cc1_active_consec1_rows = sum(is_cc1_active_consec1_window, na.rm = TRUE),
      n_model_rows = nrow(model_data),
      n_model_halfhour_bins = n_distinct(model_data$HalfHourWithinCC0),
      n_model_animals = n_distinct(model_data$AnimalID),
      n_model_groups = n_distinct(model_data$Group),
      model_halfhour_bins = paste(sort(unique(model_data$HalfHourWithinCC0)), collapse = ";"),
      .groups = "drop"
    )

  if (save_tables) {
    write_review_csv(cc1_precheck, review_cc1_active_gamm_stats_dir, "cc1_active_phase1", "model_data_precheck", sex_value)
  }

  if (nrow(model_data) < 50) {
    log_message(paste("Skipping CC1 active-phase-1 GAMM for", sex_value, ": fewer than 50 rows"))
    writeLines(
      c(paste("CC1 active-phase-1 GAMM skipped for", sex_value, ": fewer than 50 model rows."), capture.output(print(cc1_precheck))),
      file.path(review_cc1_active_gamm_diagnostics_dir, paste0("cc1_active_phase1_gamm_skipped_", sex_value, ".txt"))
    )
    return(NULL)
  }

  if (n_distinct(model_data$Group) < 2) {
    log_message(paste("Skipping CC1 active-phase-1 GAMM for", sex_value, ": fewer than two groups"))
    writeLines(
      c(paste("CC1 active-phase-1 GAMM skipped for", sex_value, ": fewer than two groups."), capture.output(print(cc1_precheck))),
      file.path(review_cc1_active_gamm_diagnostics_dir, paste0("cc1_active_phase1_gamm_skipped_", sex_value, ".txt"))
    )
    return(NULL)
  }

  n_timepoints <- n_distinct(model_data$HalfHourWithinCC0)
  if (n_timepoints < 5) {
    log_message(paste("Skipping CC1 active-phase-1 GAMM for", sex_value, ": fewer than five half-hour bins"))
    writeLines(
      c(paste("CC1 active-phase-1 GAMM skipped for", sex_value, ": fewer than five half-hour bins."), capture.output(print(cc1_precheck))),
      file.path(review_cc1_active_gamm_diagnostics_dir, paste0("cc1_active_phase1_gamm_skipped_", sex_value, ".txt"))
    )
    return(NULL)
  }

  cc1_k <- max(3L, min(GAMM_K_CC1_ACTIVE, n_timepoints - 1L))
  rho_est <- if (GAMM_ESTIMATE_RHO) estimate_rho_from_series(model_data) else 0
  log_message(paste0("CC1 active-phase-1 GAMM rho for ", sex_value, " = ", round(rho_est, 3)))

  model <- mgcv::bam(
    animalEntropy ~ Group +
      s(HalfHourWithinCC0, k = cc1_k) +
      s(HalfHourWithinCC0, by = Group, k = cc1_k) +
      s(AnimalID, bs = "re") +
      s(System, bs = "re"),
    data = model_data,
    method = "fREML",
    discrete = TRUE,
    gamma = GAMM_GAMMA,
    rho = rho_est,
    AR.start = model_data$AR_start
  )

  list(
    model = model,
    data = model_data,
    sex = sex_value,
    rho = rho_est,
    k = cc1_k
  )
}

make_cc1_active_prediction_grid <- function(cc1_result) {
  if (is.null(cc1_result)) return(tibble())

  model_data <- cc1_result$data
  model <- cc1_result$model

  ref_animal <- levels(model_data$AnimalID)[1]
  ref_system <- levels(model_data$System)[1]

  grid <- tidyr::expand_grid(
    HalfHourWithinCC0 = 0:23,
    Group = levels(model_data$Group)
  ) %>%
    filter(!is.na(Group)) %>%
    mutate(
      Group = factor(Group, levels = levels(model_data$Group)),
      Sex = factor(cc1_result$sex, levels = c("male", "female")),
      AnimalID = factor(ref_animal, levels = levels(model_data$AnimalID)),
      System = factor(ref_system, levels = levels(model_data$System))
    )

  pred <- tryCatch(
    predict(
      model,
      newdata = grid,
      se.fit = TRUE,
      exclude = c("s(AnimalID)", "s(System)")
    ),
    error = function(e) {
      log_message(paste("CC1 active-phase-1 prediction failed:", conditionMessage(e)))
      return(NULL)
    }
  )

  if (is.null(pred)) return(tibble())

  grid %>%
    mutate(
      fit = as.numeric(pred$fit),
      se = as.numeric(pred$se.fit),
      lower = fit - 1.96 * se,
      upper = fit + 1.96 * se,
      GroupLabel = format_group_labels(Group)
    )
}

cc1_active_gamm_results <- list()
cc1_active_gamm_result <- NULL
cc1_active_gamm_prediction_df <- tibble()
cc1_active_gamm_auc <- tibble()
cc1_active_gamm_auc_contrasts <- tibble()
cc1_active_observed_auc <- tibble()
cc1_active_observed_auc_anova <- tibble()
cc1_active_observed_auc_contrasts <- tibble()

if (RUN_CC1_ACTIVE_GAMM && analyze_by_halfhour) {
  log_message("Fitting focused sex-stratified GAMMs for CC1 active phase 1")

  for (sex_value in c("male", "female")) {
    cc1_active_gamm_results[[sex_value]] <- fit_cc1_active_phase1_gamm_by_sex(animal_halfhour, sex_value)
  }

  cc1_active_gamm_results <- cc1_active_gamm_results[!purrr::map_lgl(cc1_active_gamm_results, is.null)]
  cc1_active_gamm_result <- if (length(cc1_active_gamm_results) > 0) cc1_active_gamm_results[[1]] else NULL

  for (sex_value in names(cc1_active_gamm_results)) {
    cc1_active_res <- cc1_active_gamm_results[[sex_value]]
    cc1_active_model <- cc1_active_res$model
    cc1_active_data <- cc1_active_res$data

    cc1_active_parametric <- as.data.frame(summary(cc1_active_model)$p.table) %>%
      rownames_to_column("term") %>%
      mutate(
        model = "cc1_active_phase1_gamm",
        Sex = sex_value,
        rho = cc1_active_res$rho,
        k = cc1_active_res$k,
        across(where(is.numeric), ~ round(.x, 6))
      )

    cc1_active_smooths <- as.data.frame(summary(cc1_active_model)$s.table) %>%
      rownames_to_column("smooth_term") %>%
      mutate(
        model = "cc1_active_phase1_gamm",
        Sex = sex_value,
        rho = cc1_active_res$rho,
        k = cc1_active_res$k,
        across(where(is.numeric), ~ round(.x, 6))
      )

    cc1_active_prediction_sex <- make_cc1_active_prediction_grid(cc1_active_res)
    cc1_active_gamm_prediction_df <- bind_rows(cc1_active_gamm_prediction_df, cc1_active_prediction_sex)

    cc1_active_individual_output <- cc1_active_data %>%
      mutate(
        fitted_subject = as.numeric(predict(cc1_active_model, newdata = cc1_active_data)),
        fitted_population = as.numeric(predict(
          cc1_active_model,
          newdata = cc1_active_data,
          exclude = c("s(AnimalID)", "s(System)")
        )),
        residual_subject = animalEntropy - fitted_subject,
        residual_population = animalEntropy - fitted_population
      ) %>%
      select(
        AnimalID, Group, Sex, Batch, System, CageChange, Phase, ConsecActive,
        HalfHourWithinCC0, animalEntropy, fitted_subject, fitted_population,
        residual_subject, residual_population
      )

    cc1_active_animal_summary <- cc1_active_data %>%
      group_by(AnimalID, Group, Sex, Batch, System) %>%
      summarise(
        n_halfhours = n_distinct(HalfHourWithinCC0),
        n_observations = n(),
        animal_mean_entropy = mean(animalEntropy, na.rm = TRUE),
        animal_median_entropy = median(animalEntropy, na.rm = TRUE),
        animal_sd_entropy = sd(animalEntropy, na.rm = TRUE),
        animal_min_entropy = min(animalEntropy, na.rm = TRUE),
        animal_max_entropy = max(animalEntropy, na.rm = TRUE),
        .groups = "drop"
      )

    cc1_active_group_halfhour_summary <- cc1_active_data %>%
      group_by(Group, Sex, HalfHourWithinCC0) %>%
      summarise(
        n_observations = n(),
        n_animals = n_distinct(AnimalID),
        mean_entropy = mean(animalEntropy, na.rm = TRUE),
        sem_entropy = sem(animalEntropy),
        median_entropy = median(animalEntropy, na.rm = TRUE),
        .groups = "drop"
      )

    cc1_active_group_summary <- cc1_active_animal_summary %>%
      group_by(Group, Sex) %>%
      summarise(
        n_animals = n_distinct(AnimalID),
        mean_animal_entropy = mean(animal_mean_entropy, na.rm = TRUE),
        sem_animal_entropy = sem(animal_mean_entropy),
        median_animal_entropy = median(animal_mean_entropy, na.rm = TRUE),
        .groups = "drop"
      )

    if (save_tables) {
      write_review_csv(cc1_active_individual_output, review_cc1_active_gamm_stats_dir, "cc1_active_phase1", "individual_observed_fitted", sex_value)
      write_review_csv(cc1_active_animal_summary, review_cc1_active_gamm_stats_dir, "cc1_active_phase1", "animal_level_summary", sex_value)
      write_review_csv(cc1_active_group_halfhour_summary, review_cc1_active_gamm_stats_dir, "cc1_active_phase1", "group_halfhour_summary", sex_value)
      write_review_csv(cc1_active_group_summary, review_cc1_active_gamm_stats_dir, "cc1_active_phase1", "group_level_summary", sex_value)
      write_review_csv(cc1_active_parametric, review_cc1_active_gamm_stats_dir, "cc1_active_phase1_gamm", "parametric", sex_value)
      write_review_csv(cc1_active_smooths, review_cc1_active_gamm_stats_dir, "cc1_active_phase1_gamm", "smooths", sex_value)
      if (nrow(cc1_active_prediction_sex) > 0) {
        write_review_csv(cc1_active_prediction_sex, review_cc1_active_gamm_stats_dir, "cc1_active_phase1_gamm", "prediction_grid", sex_value)
      }
    }

    saveRDS(cc1_active_model, file.path(review_model_dir, paste0("cc1_active_phase1_gamm_animal_entropy_", sex_value, ".rds")))

    sink(file.path(review_cc1_active_gamm_diagnostics_dir, paste0("cc1_active_phase1_gamm_diagnostics_", sex_value, ".txt")))
    print(summary(cc1_active_model))
    print(mgcv::gam.check(cc1_active_model))
    print(mgcv::concurvity(cc1_active_model, full = TRUE))
    sink()
  }

  if (save_tables && nrow(cc1_active_gamm_prediction_df) > 0) {
    write_review_csv(cc1_active_gamm_prediction_df, review_cc1_active_gamm_stats_dir, "cc1_active_phase1_gamm", "prediction_grid", "all_sexes")
  }
}

if (nrow(cc1_active_gamm_prediction_df) > 0) {
  cc1_active_gamm_auc <- cc1_active_gamm_prediction_df %>%
    group_by(Sex, Group, GroupLabel) %>%
    summarise(
      gamm_predicted_auc = trapezoid_auc(HalfHourWithinCC0, fit),
      gamm_predicted_auc_lower = trapezoid_auc(HalfHourWithinCC0, lower),
      gamm_predicted_auc_upper = trapezoid_auc(HalfHourWithinCC0, upper),
      .groups = "drop"
    )

  cc1_active_gamm_auc_contrasts <- cc1_active_gamm_auc %>%
    group_by(Sex) %>%
    group_modify(function(.x, .y) {
      groups_present <- as.character(.x$Group)
      if (length(groups_present) < 2) {
        return(tibble(
          contrast = character(),
          group_1 = character(),
          group_2 = character(),
          gamm_predicted_auc_1 = numeric(),
          gamm_predicted_auc_2 = numeric(),
          gamm_predicted_auc_difference = numeric()
        ))
      }
      comps <- combn(groups_present, 2, simplify = FALSE)
      purrr::map_dfr(comps, function(comp) {
        auc_1 <- .x %>% filter(Group == comp[1]) %>% pull(gamm_predicted_auc)
        auc_2 <- .x %>% filter(Group == comp[2]) %>% pull(gamm_predicted_auc)
        tibble(
          contrast = paste(comp[1], "-", comp[2]),
          group_1 = comp[1],
          group_2 = comp[2],
          gamm_predicted_auc_1 = auc_1,
          gamm_predicted_auc_2 = auc_2,
          gamm_predicted_auc_difference = auc_1 - auc_2
        )
      })
    }) %>%
    ungroup() %>%
    mutate(across(where(is.numeric), ~ round(.x, 6)))

  if (save_tables) {
    write_review_csv(cc1_active_gamm_auc, review_cc1_active_gamm_stats_dir, "cc1_active_phase1_gamm", "predicted_auc", "all_sexes")
    write_review_csv(cc1_active_gamm_auc_contrasts, review_cc1_active_gamm_stats_dir, "cc1_active_phase1_gamm", "predicted_auc_contrasts", "all_sexes")
  }
}

if (length(cc1_active_gamm_results) > 0) {
  cc1_active_model_data_all <- bind_rows(purrr::map(cc1_active_gamm_results, "data"))

  cc1_active_observed_auc <- cc1_active_model_data_all %>%
    group_by(AnimalID, Group, Sex, Batch, System) %>%
    summarise(
      n_halfhours = n_distinct(HalfHourWithinCC0),
      observed_auc = trapezoid_auc(HalfHourWithinCC0, animalEntropy),
      mean_entropy = mean(animalEntropy, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(!is.na(observed_auc), n_halfhours >= 5)

  cc1_active_observed_auc_lm <- purrr::map(
    split(cc1_active_observed_auc, cc1_active_observed_auc$Sex),
    function(sex_df) {
      sex_df <- sex_df %>% droplevels()
      if (nrow(sex_df) < 6 || n_distinct(sex_df$Group) < 2) return(NULL)
      stats::lm(observed_auc ~ Group, data = sex_df)
    }
  )
  cc1_active_observed_auc_lm <- cc1_active_observed_auc_lm[!purrr::map_lgl(cc1_active_observed_auc_lm, is.null)]

  cc1_active_observed_auc_anova <- purrr::imap_dfr(
    cc1_active_observed_auc_lm,
    ~ broom::tidy(stats::anova(.x)) %>%
      mutate(Sex = .y, model = "cc1_active_phase1_observed_auc_lm") %>%
      relocate(model, Sex)
  ) %>%
    mutate(across(where(is.numeric), ~ round(.x, 6)))

  cc1_active_observed_auc_contrasts <- purrr::imap_dfr(
    cc1_active_observed_auc_lm,
    ~ emmeans::emmeans(.x, ~ Group) %>%
      emmeans::contrast(method = "pairwise", adjust = "holm") %>%
      as.data.frame(infer = TRUE) %>%
      mutate(Sex = .y, model = "cc1_active_phase1_observed_auc_lm") %>%
      relocate(model, Sex)
  ) %>%
    mutate(across(where(is.numeric), ~ round(.x, 6)))

  if (save_tables) {
    write_review_csv(cc1_active_observed_auc, review_cc1_active_gamm_stats_dir, "cc1_active_phase1", "observed_animal_auc", "all_sexes")
    write_review_csv(cc1_active_observed_auc_anova, review_cc1_active_gamm_stats_dir, "cc1_active_phase1", "observed_auc_lm_anova", "all_sexes")
    write_review_csv(cc1_active_observed_auc_contrasts, review_cc1_active_gamm_stats_dir, "cc1_active_phase1", "observed_auc_group_contrasts_holm", "all_sexes")
  }
}

# ===================================================================
# 15. ENTROPY-BASED EXPLORATION METRICS
# ===================================================================
exploration_phase <- animal_phase %>%
  mutate(
    explored_entropy_based = !is.na(animalEntropy) & animalEntropy > 0,
    exploration_entropy = ifelse(is.na(animalEntropy), 0, animalEntropy),
    exploration_category = case_when(
      is.na(animalEntropy) ~ "missing",
      animalEntropy == 0 ~ "minimal_1_position",
      animalEntropy < low_entropy_cutoff ~ "low_diversity",
      animalEntropy < moderate_entropy_cutoff ~ "moderate_diversity",
      TRUE ~ "high_diversity"
    ),
    exploration_category = factor(
      exploration_category,
      levels = c("minimal_1_position", "low_diversity", "moderate_diversity", "high_diversity", "missing")
    )
  )

exploration_animal_level <- exploration_phase %>%
  group_by(AnimalID, Group, Sex, Phase) %>%
  summarise(
    n_observations = n(),
    pct_entropy_based_exploration = 100 * mean(explored_entropy_based, na.rm = TRUE),
    mean_entropy = mean(animalEntropy, na.rm = TRUE),
    median_entropy = median(animalEntropy, na.rm = TRUE),
    .groups = "drop"
  )

exploration_lmm <- NULL

if (RUN_EXPLORATION_MODELS) {
  log_message("Fitting animal-level exploration model")

  exploration_lmm <- lmerTest::lmer(
    pct_entropy_based_exploration ~ Sex * Group * Phase +
      (1 | AnimalID),
    data = exploration_animal_level,
    REML = TRUE
  )

  exploration_anova <- format_type3_anova(exploration_lmm, "exploration_entropy_based_lmm")
  exploration_fixed <- format_lmm_fixed_effects(exploration_lmm, "exploration_entropy_based_lmm")
  exploration_emm <- emmeans::emmeans(exploration_lmm, ~ Group | Sex * Phase)
  exploration_contrasts <- emmeans::contrast(exploration_emm, method = "pairwise", adjust = "holm") %>%
    as.data.frame(infer = TRUE) %>%
    mutate(across(where(is.numeric), ~ round(.x, 6)))

  if (save_tables) {
    write_review_csv(exploration_phase, review_tables_dir, "exploration", "entropy_based_observations", "phase")
    write_review_csv(exploration_animal_level, review_tables_dir, "exploration", "entropy_based_animal_level", "phase")
    write_review_csv(exploration_anova, review_tables_dir, "exploration", "lmm_type3_anova", "phase")
    write_review_csv(exploration_fixed, review_tables_dir, "exploration", "lmm_fixed_effects", "phase")
    write_review_csv(as.data.frame(exploration_emm), review_tables_dir, "exploration", "emmeans", "phase")
    write_review_csv(exploration_contrasts, review_tables_dir, "exploration", "group_contrasts_holm", "phase")
  }

  saveRDS(exploration_lmm, file.path(review_model_dir, "exploration_entropy_based_lmm.rds"))
}

# ===================================================================
# 16. PUBLICATION-READY FIGURES
# ===================================================================
# Main manuscript emphasis:
#   A. Individual-level entropy: animal-level spatial diversity.
#   B. CON/RES/SUS group differences in individual entropy.
#   C. Time-resolved GAMM trajectories of individual entropy, styled consistently
#      with the trajectory_firstChangeActive reference figure.
#
# Cage/system-level entropy is intentionally not used for the main figures.
# System remains only as a random effect / blocking factor where needed.

# -------------------------------------------------------------------
# Trajectory-style plotting utilities matching the reference GAMM figure
# -------------------------------------------------------------------
theme_trajectory_reference <- function(base_size = 7, base_family = "sans") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.28, colour = "#1A1A1A"),
      axis.ticks = ggplot2::element_line(linewidth = 0.28, colour = "#1A1A1A"),
      axis.ticks.length = grid::unit(1.4, "mm"),
      axis.text = ggplot2::element_text(size = base_size, colour = "#1A1A1A"),
      axis.title = ggplot2::element_text(size = base_size + 1, colour = "#1A1A1A"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = base_size + 1, colour = "#1A1A1A"),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = base_size, colour = "#1A1A1A"),
      legend.text = ggplot2::element_text(size = base_size, colour = "#1A1A1A"),
      legend.key.size = grid::unit(3.0, "mm"),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = base_size + 2, face = "plain", hjust = 0, colour = "#1A1A1A"),
      plot.subtitle = ggplot2::element_text(size = base_size, hjust = 0, colour = "#1A1A1A"),
      plot.margin = grid::unit(c(2, 2, 2, 2), "mm")
    )
}

primary_emm_df <- if (exists("primary_emm_group")) as.data.frame(primary_emm_group) else NULL

fig_animal_phase <- animal_phase %>%
  ggplot(aes(x = Group, y = animalEntropy, colour = Group)) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0),
    size = 0.65,
    alpha = 0.25
  ) +
  {
    if (!is.null(primary_emm_df)) {
      geom_pointrange(
        data = primary_emm_df,
        aes(x = Group, y = emmean, ymin = lower.CL, ymax = upper.CL),
        inherit.aes = FALSE,
        colour = "black",
        linewidth = 0.28,
        size = 0.35
      )
    }
  } +
  facet_grid(Phase ~ Sex) +
  scale_colour_manual(values = group_cols, breaks = c("con", "res", "sus")) +
  labs(
    x = NULL,
    y = "Individual spatial entropy",
    title = "Individual spatial entropy"
  ) +
  theme_publication() +
  theme(legend.position = "none")

fig_animal_time <- summary_animal_halfhour %>%
  ggplot(aes(x = ConsecHalfHour, y = mean_entropy, colour = Group, fill = Group)) +
  geom_ribbon(
    aes(ymin = mean_entropy - sem_entropy, ymax = mean_entropy + sem_entropy),
    alpha = 0.16,
    colour = NA
  ) +
  geom_line(linewidth = 0.45) +
  facet_grid(Sex ~ Phase, scales = "free_x") +
  scale_colour_manual(values = group_cols, breaks = c("con", "res", "sus")) +
  scale_fill_manual(values = group_cols, breaks = c("con", "res", "sus")) +
  labs(
    x = "Consecutive half-hour",
    y = "Individual spatial entropy",
    title = "Time-resolved individual entropy"
  ) +
  theme_publication() +
  theme(legend.position = "top")

fig_exploration <- exploration_animal_level %>%
  ggplot(aes(x = Group, y = pct_entropy_based_exploration, colour = Group)) +
  geom_point(
    position = position_jitter(width = 0.12, height = 0),
    size = 0.9,
    alpha = 0.45
  ) +
  stat_summary(fun = mean, geom = "point", colour = "black", size = 1.2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", colour = "black", width = 0.15, linewidth = 0.25) +
  facet_grid(Phase ~ Sex) +
  scale_colour_manual(values = group_cols, breaks = c("con", "res", "sus")) +
  labs(
    x = NULL,
    y = "Entropy-based exploration (%)",
    title = "Entropy-based exploration frequency"
  ) +
  theme_publication() +
  theme(legend.position = "none")

# -------------------------------------------------------------------
# Additional publication-style panels
# -------------------------------------------------------------------
fig_animal_violin <- animal_phase %>%
  ggplot(aes(x = Group, y = animalEntropy, fill = Group, colour = Group)) +
  geom_violin(width = 0.82, alpha = 0.35, linewidth = 0.22, trim = FALSE) +
  geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.7, linewidth = 0.22, colour = "black") +
  geom_point(
    position = position_jitter(width = 0.08, height = 0),
    size = 0.45,
    alpha = 0.18
  ) +
  facet_grid(Phase ~ Sex) +
  scale_fill_manual(values = group_cols, breaks = c("con", "res", "sus")) +
  scale_colour_manual(values = group_cols, breaks = c("con", "res", "sus")) +
  labs(
    x = NULL,
    y = "Individual spatial entropy",
    title = "Distribution of individual entropy"
  ) +
  theme_publication() +
  theme(legend.position = "none")

fig_animal_effect_sizes <- if (exists("primary_effect_sizes")) {
  primary_effect_sizes %>%
    ggplot(aes(x = contrast, y = cohens_d)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_point(size = 1.1, colour = "black") +
    facet_grid(Phase ~ Sex) +
    labs(
      x = NULL,
      y = "Cohen's d",
      title = "Pairwise standardized effects"
    ) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
} else {
  NULL
}

fig_sex_difference <- animal_phase %>%
  group_by(AnimalID, Group, Sex, Phase) %>%
  summarise(mean_entropy = mean(animalEntropy, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = Sex, y = mean_entropy, colour = Sex)) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0),
    size = 0.75,
    alpha = 0.45
  ) +
  stat_summary(fun = mean, geom = "point", colour = "black", size = 1.2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", colour = "black", width = 0.14, linewidth = 0.25) +
  facet_grid(Phase ~ Group) +
  scale_colour_manual(values = sex_cols, breaks = c("male", "female")) +
  labs(
    x = NULL,
    y = "Mean entropy per animal",
    title = "Sex differences within stress-response groups"
  ) +
  theme_publication() +
  theme(legend.position = "none")

fig_halfhour_heatmap <- animal_halfhour %>%
  group_by(Group, Sex, Phase, ConsecHalfHour) %>%
  summarise(mean_entropy = mean(animalEntropy, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = ConsecHalfHour, y = Group, fill = mean_entropy)) +
  geom_tile(linewidth = 0) +
  facet_grid(Phase ~ Sex, scales = "free_x") +
  scale_fill_gradient(low = "white", high = "black", name = "Entropy") +
  labs(
    x = "Consecutive half-hour",
    y = NULL,
    title = "Group-level entropy heatmap"
  ) +
  theme_publication() +
  theme(legend.position = "right")

fig_entropy_phase_trajectory <- animal_phase %>%
  mutate(
    PhaseNumberGlobal = ifelse(Phase == "active", ConsecActiveGlobal, ConsecInactiveGlobal)
  ) %>%
  group_by(Group, Sex, Phase, PhaseNumberGlobal) %>%
  summarise(
    mean_entropy = mean(animalEntropy, na.rm = TRUE),
    sem_entropy = sem(animalEntropy),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = PhaseNumberGlobal, y = mean_entropy, colour = Group, fill = Group)) +
  geom_ribbon(aes(ymin = mean_entropy - sem_entropy, ymax = mean_entropy + sem_entropy), alpha = 0.14, colour = NA) +
  geom_line(linewidth = 0.42) +
  geom_point(size = 0.65) +
  facet_grid(Phase ~ Sex, scales = "free_x") +
  scale_colour_manual(values = group_cols, breaks = c("con", "res", "sus")) +
  scale_fill_manual(values = group_cols, breaks = c("con", "res", "sus")) +
  labs(
    x = "Consecutive phase number",
    y = "Individual spatial entropy",
    title = "Phase-resolved entropy trajectory"
  ) +
  theme_publication() +
  theme(legend.position = "top")

# -------------------------------------------------------------------
# Exploratory quality-control figures
# -------------------------------------------------------------------
fig_qc_raw_entropy_by_phase <- animal_phase %>%
  ggplot(aes(x = Phase, y = animalEntropy, colour = Group)) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.55),
    size = 0.55,
    alpha = 0.22
  ) +
  facet_grid(Sex ~ CageChange) +
  scale_colour_manual(values = group_cols, breaks = c("con", "res", "sus")) +
  labs(
    x = NULL,
    y = "Individual spatial entropy",
    title = "Raw phase-level entropy by cage change"
  ) +
  theme_publication() +
  theme(legend.position = "top")

fig_qc_animal_means <- animal_phase %>%
  group_by(AnimalID, Group, Sex, Phase) %>%
  summarise(mean_entropy = mean(animalEntropy, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = Group, y = mean_entropy, colour = Group)) +
  geom_point(
    position = position_jitter(width = 0.11, height = 0),
    size = 0.8,
    alpha = 0.55
  ) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.45, colour = "black", linewidth = 0.25) +
  facet_grid(Phase ~ Sex) +
  scale_colour_manual(values = group_cols, breaks = c("con", "res", "sus")) +
  labs(
    x = NULL,
    y = "Mean entropy per animal",
    title = "Animal-level entropy means"
  ) +
  theme_publication() +
  theme(legend.position = "none")

fig_qc_missingness <- bind_rows(
  animal_phase %>%
    count(Group, Sex, Phase, CageChange, name = "n_observations") %>%
    mutate(metric = "animal_phase_rows"),
  animal_halfhour %>%
    count(Group, Sex, Phase, CageChange, name = "n_observations") %>%
    mutate(metric = "animal_halfhour_rows")
) %>%
  ggplot(aes(x = CageChange, y = n_observations, fill = Group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  facet_grid(metric + Phase ~ Sex, scales = "free_y") +
  scale_fill_manual(values = group_cols, breaks = c("con", "res", "sus")) +
  labs(
    x = NULL,
    y = "Number of observations",
    title = "Observation counts by analysis table"
  ) +
  theme_publication() +
  theme(legend.position = "top")

fig_qc_entropy_density <- animal_phase %>%
  ggplot(aes(x = animalEntropy, colour = Group, fill = Group)) +
  geom_density(alpha = 0.16, linewidth = 0.35, adjust = 1.1) +
  facet_grid(Phase ~ Sex) +
  scale_colour_manual(values = group_cols, breaks = c("con", "res", "sus")) +
  scale_fill_manual(values = group_cols, breaks = c("con", "res", "sus")) +
  labs(
    x = "Individual spatial entropy",
    y = "Density",
    title = "Entropy distribution density"
  ) +
  theme_publication() +
  theme(legend.position = "top")

# -------------------------------------------------------------------
# GAMM prediction figures
# -------------------------------------------------------------------
make_gamm_prediction_grid <- function(gamm_result, sex_value, n_grid = 200) {
  if (is.null(gamm_result)) return(NULL)

  data_sex <- gamm_result$data %>%
    mutate(
      Group = droplevels(factor(Group, levels = c("con", "res", "sus"))),
      Phase = droplevels(factor(Phase, levels = c("active", "inactive"))),
      CageChange = droplevels(factor(as.character(CageChange))),
      AnimalID = droplevels(factor(as.character(AnimalID))),
      System = droplevels(factor(as.character(System)))
    )

  model <- gamm_result$model

  group_levels_present <- levels(data_sex$Group)
  phase_levels_present <- levels(data_sex$Phase)
  cage_levels_present <- levels(data_sex$CageChange)
  animal_levels_present <- levels(data_sex$AnimalID)
  system_levels_present <- levels(data_sex$System)

  group_levels_present <- group_levels_present[!is.na(group_levels_present)]
  phase_levels_present <- phase_levels_present[!is.na(phase_levels_present)]
  cage_levels_present <- cage_levels_present[!is.na(cage_levels_present)]
  animal_levels_present <- animal_levels_present[!is.na(animal_levels_present)]
  system_levels_present <- system_levels_present[!is.na(system_levels_present)]

  if (length(group_levels_present) < 1 || length(phase_levels_present) < 1 ||
      length(cage_levels_present) < 1 || length(animal_levels_present) < 1 ||
      length(system_levels_present) < 1) {
    log_message(paste("Skipping GAMM prediction grid for", sex_value, ": missing valid factor levels"))
    return(NULL)
  }

  most_common_cage <- names(sort(table(data_sex$CageChange), decreasing = TRUE))[1]
  ref_animal <- animal_levels_present[1]
  ref_system <- system_levels_present[1]

  z_seq <- seq(
    min(data_sex$ConsecHalfHour_z, na.rm = TRUE),
    max(data_sex$ConsecHalfHour_z, na.rm = TRUE),
    length.out = n_grid
  )

  # Use the empirical mean time per scaled z value for a readable secondary table.
  # The plot itself uses TimeScaled, matching the reference trajectory style.
  hh_center <- mean(data_sex$ConsecHalfHour, na.rm = TRUE)
  hh_scale <- stats::sd(data_sex$ConsecHalfHour, na.rm = TRUE)
  if (!is.finite(hh_scale) || hh_scale == 0) hh_scale <- 1

  grid <- tidyr::expand_grid(
    Group = group_levels_present,
    Phase = phase_levels_present,
    ConsecHalfHour_z = z_seq
  ) %>%
    mutate(
      Group = factor(Group, levels = levels(data_sex$Group)),
      Phase = factor(Phase, levels = levels(data_sex$Phase)),
      CageChange = factor(most_common_cage, levels = levels(data_sex$CageChange)),
      HalfHourInPhase_z = 0,
      AnimalID = factor(ref_animal, levels = levels(data_sex$AnimalID)),
      System = factor(ref_system, levels = levels(data_sex$System)),
      Sex = sex_value,
      TimeScaled = ConsecHalfHour_z,
      TimeOriginal = ConsecHalfHour_z * hh_scale + hh_center
    )

  if (anyNA(grid$CageChange) || anyNA(grid$AnimalID) || anyNA(grid$System)) {
    log_message(paste("Skipping GAMM prediction grid for", sex_value, ": prediction grid contains NA factor levels"))
    return(NULL)
  }

  pred <- tryCatch(
    predict(
      model,
      newdata = grid,
      se.fit = TRUE,
      exclude = c("s(AnimalID)", "s(System)")
    ),
    error = function(e) {
      log_message(paste("GAMM prediction failed for", sex_value, ":", conditionMessage(e)))
      return(NULL)
    }
  )

  if (is.null(pred)) return(NULL)

  out <- grid %>%
    mutate(
      fit = as.numeric(pred$fit),
      se = as.numeric(pred$se.fit),
      lower = fit - 1.96 * se,
      upper = fit + 1.96 * se,
      GroupLabel = format_group_labels(Group)
    )

  if (all(is.na(out$fit))) {
    log_message(paste("GAMM prediction for", sex_value, "returned only NA fitted values"))
  }

  out
}

gamm_prediction_df <- if (exists("gamm_results") && length(gamm_results) > 0) {
  pred_list <- purrr::imap(gamm_results, make_gamm_prediction_grid)
  pred_list <- pred_list[!purrr::map_lgl(pred_list, is.null)]
  if (length(pred_list) > 0) dplyr::bind_rows(pred_list) else tibble()
} else {
  tibble()
}

fig_gamm_prediction <- if (nrow(gamm_prediction_df) > 0) {
  gamm_prediction_df %>%
    ggplot(aes(x = TimeScaled, y = fit, colour = GroupLabel, fill = GroupLabel)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.16, colour = NA) +
    geom_line(linewidth = 0.55) +
    facet_grid(Phase ~ Sex) +
    scale_colour_manual(values = c(CON = group_cols[["con"]], RES = group_cols[["res"]], SUS = group_cols[["sus"]])) +
    scale_fill_manual(values = c(CON = group_cols[["con"]], RES = group_cols[["res"]], SUS = group_cols[["sus"]])) +
    labs(
      x = "Time (scaled)",
      y = "Predicted entropy",
      colour = "Group",
      fill = "Group",
      title = "GAMM trajectories by group",
      subtitle = "Individual-level Shannon entropy"
    ) +
    theme_trajectory_reference()
} else {
  NULL
}

fig_gamm_raw_overlay <- if (nrow(gamm_prediction_df) > 0) {
  animal_halfhour %>%
    mutate(GroupLabel = format_group_labels(Group)) %>%
    ggplot(aes(x = ConsecHalfHour_z, y = animalEntropy, colour = GroupLabel)) +
    geom_point(size = 0.30, alpha = 0.07) +
    geom_ribbon(
      data = gamm_prediction_df,
      aes(x = TimeScaled, y = fit, ymin = lower, ymax = upper, fill = GroupLabel),
      inherit.aes = FALSE,
      alpha = 0.16,
      colour = NA
    ) +
    geom_line(
      data = gamm_prediction_df,
      aes(x = TimeScaled, y = fit, colour = GroupLabel),
      inherit.aes = FALSE,
      linewidth = 0.55
    ) +
    facet_grid(Phase ~ Sex) +
    scale_colour_manual(values = c(CON = group_cols[["con"]], RES = group_cols[["res"]], SUS = group_cols[["sus"]])) +
    scale_fill_manual(values = c(CON = group_cols[["con"]], RES = group_cols[["res"]], SUS = group_cols[["sus"]])) +
    labs(
      x = "Time (scaled)",
      y = "Individual entropy",
      colour = "Group",
      fill = "Group",
      title = "GAMM trajectories over raw observations",
      subtitle = "Individual-level Shannon entropy"
    ) +
    theme_trajectory_reference()
} else {
  NULL
}

fig_gamm_group_level <- if (nrow(gamm_prediction_df) > 0) {
  # Same model predictions summarized visually as group-level trajectories.
  # This keeps the visual focus on group separation while preserving animal-level inference.
  gamm_prediction_df %>%
    ggplot(aes(x = TimeScaled, y = fit, colour = GroupLabel, fill = GroupLabel)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.16, colour = NA) +
    geom_line(linewidth = 0.65) +
    facet_wrap(~ Sex, nrow = 1) +
    scale_colour_manual(values = c(CON = group_cols[["con"]], RES = group_cols[["res"]], SUS = group_cols[["sus"]])) +
    scale_fill_manual(values = c(CON = group_cols[["con"]], RES = group_cols[["res"]], SUS = group_cols[["sus"]])) +
    labs(
      x = "Time (scaled)",
      y = "Predicted entropy",
      colour = "Group",
      fill = "Group",
      title = "Group-level entropy trajectories",
      subtitle = "GAMM estimates from individual entropy"
    ) +
    theme_trajectory_reference()
} else {
  NULL
}

cc1_active_gamm_caption <- if (exists("cc1_active_gamm_results") && length(cc1_active_gamm_results) > 0) {
  rho_text <- purrr::imap_chr(
    cc1_active_gamm_results,
    ~ paste0(.y, " rho = ", round(.x$rho, 3), ", k = ", .x$k)
  )
  paste0(
    "Sex-stratified GAMMs: animalEntropy ~ Group + s(HalfHourWithinCC0) + ",
    "s(HalfHourWithinCC0, by = Group) + s(AnimalID, bs = 're') + s(System, bs = 're'); ",
    "bam(fREML, discrete = TRUE, gamma = ",
    GAMM_GAMMA,
    ", AR.start by AnimalID/System; ",
    paste(rho_text, collapse = "; "),
    "."
  )
} else {
  NULL
}

fig_cc1_active_individual_gamm <- if (!is.null(cc1_active_gamm_result) && nrow(cc1_active_gamm_prediction_df) > 0) {
  bind_rows(purrr::map(cc1_active_gamm_results, "data")) %>%
    mutate(GroupLabel = format_group_labels(Group)) %>%
    ggplot(aes(x = HalfHourWithinCC0, y = animalEntropy, colour = GroupLabel)) +
    geom_line(aes(group = AnimalID), linewidth = 0.16, alpha = 0.14) +
    geom_point(size = 0.26, alpha = 0.12) +
    geom_ribbon(
      data = cc1_active_gamm_prediction_df,
      aes(x = HalfHourWithinCC0, y = fit, ymin = lower, ymax = upper, fill = GroupLabel),
      inherit.aes = FALSE,
      alpha = 0.14,
      colour = NA
    ) +
    geom_line(
      data = cc1_active_gamm_prediction_df,
      aes(x = HalfHourWithinCC0, y = fit, colour = GroupLabel),
      inherit.aes = FALSE,
      linewidth = 0.55
    ) +
    facet_wrap(~ Sex, nrow = 1) +
    scale_x_continuous(breaks = c(0, 6, 12, 18, 23), limits = c(0, 23)) +
    scale_y_continuous(limits = c(0, max_position_entropy), expand = expansion(mult = c(0.02, 0.04))) +
    scale_colour_manual(values = c(CON = group_cols[["con"]], RES = group_cols[["res"]], SUS = group_cols[["sus"]])) +
    scale_fill_manual(values = c(CON = group_cols[["con"]], RES = group_cols[["res"]], SUS = group_cols[["sus"]])) +
    labs(
      x = "Half-hour within CC1 active phase 1",
      y = "Individual spatial entropy",
      colour = "Group",
      fill = "Group",
      title = "CC1 active phase 1: individual trajectories",
      caption = cc1_active_gamm_caption
    ) +
    theme_nature_gamm()
} else {
  NULL
}

fig_cc1_active_group_gamm <- if (nrow(cc1_active_gamm_prediction_df) > 0) {
  cc1_active_gamm_prediction_df %>%
    ggplot(aes(x = HalfHourWithinCC0, y = fit, colour = GroupLabel, fill = GroupLabel)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.16, colour = NA) +
    geom_line(linewidth = 0.62) +
    facet_wrap(~ Sex, nrow = 1) +
    scale_x_continuous(breaks = c(0, 6, 12, 18, 23), limits = c(0, 23)) +
    scale_y_continuous(limits = c(0, max_position_entropy), expand = expansion(mult = c(0.02, 0.04))) +
    scale_colour_manual(values = c(CON = group_cols[["con"]], RES = group_cols[["res"]], SUS = group_cols[["sus"]])) +
    scale_fill_manual(values = c(CON = group_cols[["con"]], RES = group_cols[["res"]], SUS = group_cols[["sus"]])) +
    labs(
      x = "Half-hour within CC1 active phase 1",
      y = "Predicted entropy",
      colour = "Group",
      fill = "Group",
      title = "CC1 active phase 1: GAMM group estimates",
      caption = cc1_active_gamm_caption
    ) +
    theme_nature_gamm()
} else {
  NULL
}

fig_cc1_active_group_observed <- if (!is.null(cc1_active_gamm_result)) {
  bind_rows(purrr::map(cc1_active_gamm_results, "data")) %>%
    group_by(Group, Sex, HalfHourWithinCC0) %>%
    summarise(
      mean_entropy = mean(animalEntropy, na.rm = TRUE),
      sem_entropy = sem(animalEntropy),
      .groups = "drop"
    ) %>%
    mutate(GroupLabel = format_group_labels(Group)) %>%
    ggplot(aes(x = HalfHourWithinCC0, y = mean_entropy, colour = GroupLabel, fill = GroupLabel)) +
    geom_ribbon(aes(ymin = mean_entropy - sem_entropy, ymax = mean_entropy + sem_entropy), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.42) +
    geom_point(size = 0.6) +
    facet_wrap(~ Sex, nrow = 1) +
    scale_x_continuous(breaks = c(0, 6, 12, 18, 23), limits = c(0, 23)) +
    scale_y_continuous(limits = c(0, max_position_entropy), expand = expansion(mult = c(0.02, 0.04))) +
    scale_colour_manual(values = c(CON = group_cols[["con"]], RES = group_cols[["res"]], SUS = group_cols[["sus"]])) +
    scale_fill_manual(values = c(CON = group_cols[["con"]], RES = group_cols[["res"]], SUS = group_cols[["sus"]])) +
    labs(
      x = "Half-hour within CC1 active phase 1",
      y = "Mean individual entropy",
      colour = "Group",
      fill = "Group",
      title = "CC1 active phase 1: observed group means",
      caption = "Observed means are descriptive: points/lines show group mean individual entropy; ribbons show mean +/- SEM."
    ) +
    theme_nature_gamm()
} else {
  NULL
}

fig_cc1_active_auc <- if (nrow(cc1_active_gamm_auc) > 0) {
  auc_points <- if (nrow(cc1_active_observed_auc) > 0) {
    cc1_active_observed_auc %>%
      mutate(GroupLabel = format_group_labels(Group))
  } else {
    tibble()
  }

  cc1_active_gamm_auc %>%
    mutate(GroupLabel = format_group_labels(Group)) %>%
    ggplot(aes(x = GroupLabel, y = gamm_predicted_auc, colour = GroupLabel, fill = GroupLabel)) +
    {
      if (nrow(auc_points) > 0) {
        geom_point(
          data = auc_points,
          aes(x = GroupLabel, y = observed_auc, colour = GroupLabel),
          inherit.aes = FALSE,
          position = position_jitter(width = 0.08, height = 0),
          size = 0.75,
          alpha = 0.35
        )
      }
    } +
    geom_pointrange(
      aes(ymin = gamm_predicted_auc_lower, ymax = gamm_predicted_auc_upper),
      linewidth = 0.28,
      size = 0.45,
      fatten = 1.6
    ) +
    facet_wrap(~ Sex, nrow = 1) +
    scale_colour_manual(values = c(CON = group_cols[["con"]], RES = group_cols[["res"]], SUS = group_cols[["sus"]])) +
    scale_fill_manual(values = c(CON = group_cols[["con"]], RES = group_cols[["res"]], SUS = group_cols[["sus"]])) +
    labs(
      x = NULL,
      y = "AUC across CC1 active phase 1",
      colour = "Group",
      fill = "Group",
      title = "CC1 active phase 1: entropy AUC",
      caption = "Large points/ranges show GAMM-predicted AUC integrated over half-hours 0-23 using fitted values and pointwise 95% CI curves; small points show animal-level observed AUC."
    ) +
    theme_nature_gamm() +
    theme(legend.position = "none")
} else {
  NULL
}

fig_cc1_active_gamm_auc_contrasts <- if (nrow(cc1_active_gamm_auc_contrasts) > 0) {
  cc1_active_gamm_auc_contrasts %>%
    mutate(
      contrast = factor(contrast, levels = c("con - res", "con - sus", "res - sus"))
    ) %>%
    ggplot(aes(x = contrast, y = gamm_predicted_auc_difference)) +
    geom_hline(yintercept = 0, linewidth = 0.28, colour = "grey45") +
    geom_point(size = 1.15, colour = "black") +
    geom_segment(
      aes(x = contrast, xend = contrast, y = 0, yend = gamm_predicted_auc_difference),
      linewidth = 0.3,
      colour = "black"
    ) +
    facet_wrap(~ Sex, nrow = 1) +
    labs(
      x = NULL,
      y = "GAMM-predicted AUC difference",
      title = "CC1 active phase 1: GAMM AUC contrasts",
      caption = "Contrasts are pairwise differences between sex-stratified GAMM-predicted AUC values integrated over half-hours 0-23. Positive values indicate higher AUC for the first group in the contrast."
    ) +
    theme_nature_gamm() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 35, hjust = 1)
    )
} else {
  NULL
}

if (save_plots) {
  # Main figure panels: saved as SVG + PDF
  save_review_figure_dual(fig_animal_phase, review_main_panels_dir, "figure", "fig1a_individual_entropy_emm", "phase", width = 90, height = 75)
  save_review_figure_dual(fig_exploration, review_main_panels_dir, "figure", "fig1b_entropy_based_exploration", "phase", width = 90, height = 75)
  save_review_figure_dual(fig_animal_time, review_main_panels_dir, "figure", "fig2a_individual_entropy_timecourse", "halfhour", width = 180, height = 80)

  # Exploratory QC panels
  save_review_figure_dual(fig_qc_raw_entropy_by_phase, review_exploratory_figures_dir, "qc", "raw_entropy_by_phase_cagechange", NULL, width = 150, height = 110)
  save_review_figure_dual(fig_qc_animal_means, review_exploratory_figures_dir, "qc", "animal_level_entropy_means", NULL, width = 95, height = 80)
  save_review_figure_dual(fig_qc_missingness, review_exploratory_figures_dir, "qc", "observation_counts", NULL, width = 150, height = 130)
  save_review_figure_dual(fig_qc_entropy_density, review_exploratory_figures_dir, "qc", "entropy_density", NULL, width = 95, height = 80)

  # GAMM timecourse panels
  if (!is.null(fig_gamm_prediction)) {
    save_review_figure_dual(fig_gamm_prediction, review_gamm_figures_dir, "gamm", "predicted_entropy_trajectories", NULL, width = 150, height = 95)
  }

  if (!is.null(fig_gamm_raw_overlay)) {
    save_review_figure_dual(fig_gamm_raw_overlay, review_gamm_figures_dir, "gamm", "raw_overlay_predicted_trajectories", NULL, width = 150, height = 95)
  }

  if (!is.null(fig_gamm_group_level)) {
    save_review_figure_dual(fig_gamm_group_level, review_gamm_figures_dir, "gamm", "group_level_entropy_trajectories", NULL, width = 150, height = 75)
  }

  if (!is.null(fig_cc1_active_individual_gamm)) {
    save_review_figure_dual(fig_cc1_active_individual_gamm, review_cc1_active_gamm_figures_dir, "cc1_active_phase1_gamm", "individual_trajectories", "animal_entropy", width = 150, height = 72)
  }

  if (!is.null(fig_cc1_active_group_gamm)) {
    save_review_figure_dual(fig_cc1_active_group_gamm, review_cc1_active_gamm_figures_dir, "cc1_active_phase1_gamm", "group_estimates", "animal_entropy", width = 120, height = 68)
  }

  if (!is.null(fig_cc1_active_group_observed)) {
    save_review_figure_dual(fig_cc1_active_group_observed, review_cc1_active_gamm_figures_dir, "cc1_active_phase1", "observed_group_means", "animal_entropy", width = 120, height = 68)
  }

  if (!is.null(fig_cc1_active_auc)) {
    save_review_figure_dual(fig_cc1_active_auc, review_cc1_active_gamm_figures_dir, "cc1_active_phase1_gamm", "auc_contrasts", "animal_entropy", width = 110, height = 72)
  }

  if (!is.null(fig_cc1_active_gamm_auc_contrasts)) {
    save_review_figure_dual(fig_cc1_active_gamm_auc_contrasts, review_cc1_active_gamm_figures_dir, "cc1_active_phase1_gamm", "auc_pairwise_contrast_differences", "animal_entropy", width = 95, height = 68)
  }

  if (exists("gamm_prediction_df") && nrow(gamm_prediction_df) > 0 && save_tables) {
    write_review_csv(gamm_prediction_df, review_gamm_stats_dir, "gamm", "prediction_grid", "animal_entropy")
  }

  # Supplementary panels
  save_review_figure_dual(fig_animal_violin, review_supp_panels_dir, "figure_s", "individual_entropy_violin_raw", "phase", width = 90, height = 75)
  save_review_figure_dual(fig_sex_difference, review_supp_panels_dir, "figure_s", "sex_differences_by_group", "animal_level", width = 100, height = 85)
  save_review_figure_dual(fig_halfhour_heatmap, review_supp_panels_dir, "figure_s", "halfhour_entropy_heatmap", "group", width = 150, height = 70)
  save_review_figure_dual(fig_entropy_phase_trajectory, review_supp_panels_dir, "figure_s", "phase_resolved_entropy_trajectory", NULL, width = 150, height = 85)

  if (!is.null(fig_animal_effect_sizes)) {
    save_review_figure_dual(fig_animal_effect_sizes, review_supp_panels_dir, "figure_s", "pairwise_effect_sizes", "cohens_d", width = 95, height = 80)
  }

  # Main multipanel figures
  figure_1 <- (fig_animal_phase | fig_exploration) +
    patchwork::plot_annotation(tag_levels = "A")

  figure_2 <- fig_animal_time +
    patchwork::plot_annotation(tag_levels = "A")

  save_review_figure_dual(figure_1, review_main_multipanel_dir, "figure", "figure_1_individual_entropy_overview", NULL, width = 140, height = 80)
  save_review_figure_dual(figure_2, review_main_multipanel_dir, "figure", "figure_2_entropy_timecourse", NULL, width = 180, height = 80)

  # Supplementary multipanel figure
  figure_s1 <- (fig_animal_violin | fig_sex_difference) / (fig_entropy_phase_trajectory | fig_halfhour_heatmap) +
    patchwork::plot_annotation(tag_levels = "A")

  save_review_figure_dual(figure_s1, review_supp_multipanel_dir, "figure_s", "supplementary_entropy_distributions_and_time", NULL, width = 180, height = 150)
}

if (show_plots) {
  print(fig_animal_phase)
  if (exists("fig_cage_phase")) print(fig_cage_phase)
  print(fig_exploration)
  print(fig_animal_time)
  if (!is.null(fig_cc1_active_individual_gamm)) print(fig_cc1_active_individual_gamm)
  if (!is.null(fig_cc1_active_group_gamm)) print(fig_cc1_active_group_gamm)
}

# ===================================================================
# 17. TRANSITION NETWORK METRICS AT ANIMAL LEVEL
# ===================================================================
load_position_data_for_networks <- function() {
  all_preprocessed <- list()

  for (batch in batches) {
    for (cageChange in cageChanges) {
      filename <- paste0("E9_SIS_", batch, "_", cageChange, "_AnimalPos")
      csvFilePath <- file.path(preprocessed_data_dir, paste0(filename, "_preprocessed.csv"))

      if (!file.exists(csvFilePath)) next

      temp_data <- readr::read_csv(csvFilePath, show_col_types = FALSE) %>%
        as_tibble() %>%
        mutate(
          Batch = if ("Batch" %in% colnames(.)) Batch else batch,
          CageChange = if ("CageChange" %in% colnames(.)) CageChange else cageChange,
          Sex = batch_to_sex(Batch),
          Group = assign_group(AnimalID)
        )

      all_preprocessed[[paste(batch, cageChange, sep = "_")]] <- temp_data
    }
  }

  bind_rows(all_preprocessed)
}

detect_position_column <- function(data) {
  possible_pos_cols <- c("PositionID", "Position", "position", "Pos", "pos", "position_id", "Zone", "zone")
  pos_candidates <- intersect(possible_pos_cols, colnames(data))
  if (length(pos_candidates) == 0) return(NA_character_)
  if ("PositionID" %in% pos_candidates) return("PositionID")
  pos_candidates[1]
}

ensure_datetime_column <- function(data) {
  if ("DateTime" %in% names(data)) return(data)

  dt_candidates <- intersect(
    c("DateTime", "Date_Time", "Timestamp", "timestamp", "datetime", "Date", "Time"),
    names(data)
  )

  if (length(dt_candidates) > 0) {
    data$DateTime <- data[[dt_candidates[1]]]
  } else {
    data$DateTime <- seq_len(nrow(data))
  }

  data
}

compute_transition_metrics <- function(df, pos_col) {
  seq_positions <- df %>%
    arrange(DateTime) %>%
    pull(.data[[pos_col]])

  seq_positions <- seq_positions[!is.na(seq_positions)]

  if (length(seq_positions) < 2) {
    return(tibble(
      n_transitions = 0,
      n_unique_edges = NA_real_,
      transition_entropy = NA_real_,
      normalized_transition_entropy = NA_real_,
      self_transition_rate = NA_real_
    ))
  }

  edges <- tibble(
    from = as.character(head(seq_positions, -1)),
    to = as.character(tail(seq_positions, -1))
  )

  edge_counts <- edges %>%
    count(from, to, name = "weight")

  probs <- edge_counts$weight / sum(edge_counts$weight)
  transition_entropy <- -sum(probs * log2(probs))
  max_edge_entropy <- log2(nrow(edge_counts))

  tibble(
    n_transitions = nrow(edges),
    n_unique_edges = nrow(edge_counts),
    transition_entropy = transition_entropy,
    normalized_transition_entropy = ifelse(max_edge_entropy > 0, transition_entropy / max_edge_entropy, NA_real_),
    self_transition_rate = mean(edges$from == edges$to, na.rm = TRUE)
  )
}

plot_transition_network_descriptive <- function(df, pos_col, title_text) {
  seq_positions <- df %>%
    arrange(DateTime) %>%
    pull(.data[[pos_col]])

  seq_positions <- seq_positions[!is.na(seq_positions)]

  if (length(seq_positions) < 2) return(NULL)

  edges <- tibble(
    from = as.character(head(seq_positions, -1)),
    to = as.character(tail(seq_positions, -1))
  ) %>%
    count(from, to, name = "weight")

  g <- igraph::graph_from_data_frame(edges, directed = TRUE)
  igraph::E(g)$weight <- edges$weight

  ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_fan(
      aes(width = weight, alpha = weight),
      arrow = grid::arrow(length = grid::unit(2.2, "mm"), type = "closed"),
      colour = "grey30",
      show.legend = FALSE
    ) +
    ggraph::scale_edge_width(range = c(0.15, 2.4)) +
    ggraph::geom_node_point(size = 2.2, colour = "black") +
    ggraph::geom_node_text(aes(label = name), repel = TRUE, size = 2.2) +
    ggtitle(title_text) +
    theme_void(base_family = "sans") +
    theme(plot.title = element_text(size = 7, hjust = 0))
}

transition_metrics <- NULL
transition_lmm <- NULL

if (RUN_TRANSITION_METRICS) {
  log_message("Computing animal-level transition-network metrics")

  position_data <- load_position_data_for_networks()
  pos_col <- detect_position_column(position_data)

  if (nrow(position_data) == 0 || is.na(pos_col)) {
    log_message("Transition metrics skipped: no position data or no position column")
  } else {
    position_data <- position_data %>%
      ensure_datetime_column() %>%
      mutate(
        Group = factor(Group, levels = c("con", "res", "sus")),
        Sex = factor(Sex, levels = c("male", "female")),
        CageChange = factor(CageChange, levels = c("CC1", "CC2", "CC3", "CC4")),
        Phase = case_when(
          "ConsecActive" %in% names(.) & !is.na(ConsecActive) & ConsecActive > 0 ~ "active",
          "ConsecInactive" %in% names(.) & !is.na(ConsecInactive) & ConsecInactive > 0 ~ "inactive",
          TRUE ~ NA_character_
        ),
        Phase = factor(Phase, levels = c("active", "inactive"))
      )

    transition_metrics <- position_data %>%
      filter(!is.na(Phase), !is.na(Group), !is.na(Sex)) %>%
      group_by(AnimalID, Group, Sex, Batch, Phase) %>%
      group_modify(~ compute_transition_metrics(.x, pos_col)) %>%
      ungroup()

    transition_model_data <- transition_metrics %>%
      filter(!is.na(transition_entropy), n_transitions >= 10)

    if (nrow(transition_model_data) > 10 && n_distinct(transition_model_data$Group) > 1) {
      transition_lmm <- lmerTest::lmer(
        transition_entropy ~ Sex * Group * Phase +
          (1 | AnimalID),
        data = transition_model_data,
        REML = TRUE
      )

      transition_anova <- format_type3_anova(transition_lmm, "transition_entropy_lmm")
      transition_fixed <- format_lmm_fixed_effects(transition_lmm, "transition_entropy_lmm")
      transition_emm <- emmeans::emmeans(transition_lmm, ~ Group | Sex * Phase)
      transition_contrasts <- emmeans::contrast(transition_emm, method = "pairwise", adjust = "holm") %>%
        as.data.frame(infer = TRUE) %>%
        mutate(across(where(is.numeric), ~ round(.x, 6)))

      if (save_tables) {
        write_review_csv(transition_metrics, review_transition_dir, "transition", "animal_level_metrics", NULL)
        write_review_csv(transition_anova, review_transition_dir, "transition", "lmm_type3_anova", "entropy")
        write_review_csv(transition_fixed, review_transition_dir, "transition", "lmm_fixed_effects", "entropy")
        write_review_csv(as.data.frame(transition_emm), review_transition_dir, "transition", "emmeans", "entropy")
        write_review_csv(transition_contrasts, review_transition_dir, "transition", "group_contrasts_holm", "entropy")
      }

      saveRDS(transition_lmm, file.path(review_model_dir, "transition_entropy_lmm.rds"))

      fig_transition_entropy <- transition_model_data %>%
        ggplot(aes(x = Group, y = transition_entropy, colour = Group)) +
        geom_point(
          position = position_jitter(width = 0.12, height = 0),
          size = 0.9,
          alpha = 0.45
        ) +
        stat_summary(fun = mean, geom = "point", colour = "black", size = 1.2) +
        stat_summary(fun.data = mean_se, geom = "errorbar", colour = "black", width = 0.15, linewidth = 0.25) +
        facet_grid(Phase ~ Sex) +
        scale_colour_manual(values = group_cols, breaks = c("con", "res", "sus")) +
        labs(
          x = NULL,
          y = "Transition entropy",
          title = "Animal-level transition entropy"
        ) +
        theme_publication() +
        theme(legend.position = "none")

      if (save_plots) {
        save_review_figure_dual(fig_transition_entropy, review_transition_figures_dir, "figure", "transition_entropy_animal_level", NULL, width = 90, height = 75)
      }
    }

    # Descriptive network plots only, grouped by Sex x Group x Phase
    network_plot_dir <- file.path(review_transition_dir, "descriptive_network_plots")
    dir.create(network_plot_dir, recursive = TRUE, showWarnings = FALSE)

    for (sex_value in c("male", "female")) {
      for (group_value in c("con", "res", "sus")) {
        for (phase_value in c("active", "inactive")) {
          plot_df <- position_data %>%
            filter(Sex == sex_value, Group == group_value, Phase == phase_value)

          if (nrow(plot_df) < 2) next

          p_net <- plot_transition_network_descriptive(
            plot_df,
            pos_col,
            paste(sex_value, group_value, phase_value, sep = " | ")
          )

          if (!is.null(p_net) && save_plots) {
            file_name <- paste0("network_", sex_value, "_", group_value, "_", phase_value, ".svg")
            ggsave(
              file.path(network_plot_dir, file_name),
              p_net,
              width = 80,
              height = 80,
              units = "mm",
              dpi = 600,
              bg = "white"
            )
          }
        }
      }
    }
  }
}

# ===================================================================
# 18. OUTPUT AUDIT
# ===================================================================
audit_output_directories <- function() {
  expected_outputs <- tibble::tribble(
    ~output_area, ~directory, ~enabled,
    "clean_data", review_clean_data_dir, save_tables,
    "descriptive_tables", review_tables_dir, save_tables,
    "lmm_statistics", review_lmm_stats_dir, RUN_PRIMARY_LMM && save_tables,
    "gamm_statistics", review_gamm_stats_dir, RUN_GAMM && save_tables,
    "cc1_active_phase1_gamm_statistics", review_cc1_active_gamm_stats_dir, RUN_CC1_ACTIVE_GAMM && save_tables,
    "model_rds", review_model_dir, TRUE,
    "diagnostics", review_diagnostics_dir, TRUE,
    "cc1_active_phase1_gamm_diagnostics", review_cc1_active_gamm_diagnostics_dir, RUN_CC1_ACTIVE_GAMM,
    "main_figures", review_main_figures_dir, save_plots,
    "supplementary_figures", review_supp_figures_dir, save_plots,
    "gamm_figures", review_gamm_figures_dir, RUN_GAMM && save_plots,
    "cc1_active_phase1_gamm_figures", review_cc1_active_gamm_figures_dir, RUN_CC1_ACTIVE_GAMM && save_plots,
    "transition_outputs", review_transition_dir, RUN_TRANSITION_METRICS,
    "transition_figures", review_transition_figures_dir, RUN_TRANSITION_METRICS && save_plots
  )

  expected_outputs %>%
    mutate(
      directory_exists = dir.exists(directory),
      n_files_direct = purrr::map_int(directory, ~ if (dir.exists(.x)) length(list.files(.x, recursive = FALSE, all.files = FALSE, no.. = TRUE)) else 0L),
      n_files_recursive = purrr::map_int(directory, ~ if (dir.exists(.x)) length(list.files(.x, recursive = TRUE, all.files = FALSE, no.. = TRUE)) else 0L),
      audit_status = case_when(
        !enabled ~ "not_assessed_by_runtime_flag",
        !directory_exists ~ "missing_directory",
        n_files_recursive == 0 ~ "created_but_no_outputs_detected",
        TRUE ~ "outputs_detected"
      )
    )
}

output_audit <- audit_output_directories()

if (save_tables) {
  write_review_csv(output_audit, review_logs_dir, "audit", "output_directories_and_empty_folders")
}

empty_enabled_outputs <- output_audit %>%
  filter(enabled, audit_status != "outputs_detected")

if (nrow(empty_enabled_outputs) > 0) {
  log_message(paste(
    "Output audit flagged enabled directories without detected files:",
    paste(empty_enabled_outputs$output_area, collapse = ", ")
  ))
}

# ===================================================================
# 19. RUN MANIFEST AND SESSION INFO
# ===================================================================
processing_end_time <- Sys.time()

run_manifest <- tibble::tibble(
  field = c(
    "script",
    "run_time",
    "processing_start_time",
    "processing_end_time",
    "load_existing_data",
    "analyze_by_halfhour",
    "filter_cc4_late_phases",
    "cc4_max_active_phase",
    "cc4_max_inactive_phase",
    "max_position_entropy",
    "primary_model",
    "cage_model",
    "gamm_model",
    "cc1_active_phase1_gamm_model",
    "gamm_gamma",
    "gamm_k_cc1_active",
    "exploration_definition",
    "transition_network_inference",
    "n_animal_phase_rows",
    "n_cage_phase_rows",
    "n_animal_halfhour_rows",
    "n_cage_halfhour_rows",
    "working_directory",
    "saving_directory",
    "review_directory",
    "R_version"
  ),
  value = c(
    "reviewer_ready_shannon_entropy_pipeline.R",
    as.character(Sys.time()),
    as.character(processing_start_time),
    as.character(processing_end_time),
    as.character(load_existing_data),
    as.character(analyze_by_halfhour),
    as.character(filter_cc4_late_phases),
    as.character(cc4_max_active_phase),
    as.character(cc4_max_inactive_phase),
    as.character(max_position_entropy),
    "animalEntropy ~ Sex * Group * Phase + CageChange + (1 | Batch/System) + (1 | AnimalID)",
    "CageEntropy ~ Sex * Phase + CageChange + (1 | Batch/System)",
    "animalEntropy ~ Group * Phase + CageChange + s(global time) + s(group deviations) + s(within-phase time) + random effects",
    "Sex-stratified: animalEntropy ~ Group + s(HalfHourWithinCC0) + s(HalfHourWithinCC0, by = Group) + random effects; CC1, active, ConsecActive == 1, HalfHourWithinCC0 0-23",
    as.character(GAMM_GAMMA),
    as.character(GAMM_K_CC1_ACTIVE),
    "entropy-based exploration = animalEntropy > 0; not equivalent to outside-home-cage occupancy unless occupancy is explicitly computed",
    "animal-level transition metrics only; pooled edgewise tests intentionally removed",
    as.character(nrow(animal_phase)),
    as.character(nrow(cage_phase)),
    as.character(nrow(animal_halfhour)),
    as.character(nrow(cage_halfhour)),
    working_directory,
    saving_directory,
    review_dir,
    R.Version()$version.string
  )
)

if (save_tables) {
  write_review_csv(run_manifest, review_tables_dir, "manifest", "run_settings_and_counts")
}

metadata <- list(
  parameters = list(
    load_existing_data = load_existing_data,
    analyze_by_halfhour = analyze_by_halfhour,
    filter_cc4_late_phases = filter_cc4_late_phases,
    cc4_max_active_phase = cc4_max_active_phase,
    cc4_max_inactive_phase = cc4_max_inactive_phase,
    max_position_entropy = max_position_entropy,
    run_cc1_active_gamm = RUN_CC1_ACTIVE_GAMM,
    gamm_k_cc1_active = GAMM_K_CC1_ACTIVE,
    gamm_gamma = GAMM_GAMMA
  ),
  primary_model = "animalEntropy ~ Sex * Group * Phase + CageChange + (1 | Batch/System) + (1 | AnimalID)",
  cage_model = "CageEntropy ~ Sex * Phase + CageChange + (1 | Batch/System)",
  cc1_active_phase1_gamm_model = "Sex-stratified: animalEntropy ~ Group + s(HalfHourWithinCC0) + s(HalfHourWithinCC0, by = Group) + random effects",
  exploration_definition = "entropy-based exploration = animalEntropy > 0",
  transition_network_inference = "animal-level metrics only",
  start_time = as.character(processing_start_time),
  end_time = as.character(processing_end_time),
  R_version = R.Version()$version.string,
  session_info = capture.output(sessionInfo())
)

jsonlite::write_json(metadata, file.path(review_logs_dir, "reviewed_processing_metadata.json"), pretty = TRUE)
writeLines(metadata$session_info, file.path(review_diagnostics_dir, "session_info.txt"))

log_message("Saved metadata and session info")
log_message("Reviewer-ready entropy analysis complete")

message("=======================================================================")
message(" REVIEWER-READY ENTROPY ANALYSIS PIPELINE COMPLETE")
message("=======================================================================")
message("Clean data: ", review_clean_data_dir)
message("Tables: ", review_tables_dir)
message("Figures: ", review_figures_dir)
message("Models: ", review_model_dir)
message("Diagnostics: ", review_diagnostics_dir)
message("Transition outputs: ", review_transition_dir)
message("=======================================================================")
