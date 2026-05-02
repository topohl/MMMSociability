#' @title Analysis of Shannon Entropy in Animal Positioning
#'
#' @description
#' This script performs a structured analysis of Shannon entropy to
#' quantify spatial distribution and diversity of animal positions within
#' the environment, using both phase-based and half-hour–resolved measures.
#'
#' @details
#' The pipeline supports both raw data processing and loading of precomputed results.
#' Set load_existing_data = TRUE to skip processing and load existing CSVs.
#'
#' @date October 2025
#' @authors Tobias Pohl, Anja Magister

# ===================================================================
# Load required packages
# ===================================================================
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(readr, dplyr, lubridate, tibble, purrr, ggplot2, reshape2,
               scales, stringr, gridExtra, lmerTest, mgcv, emmeans, ggpubr, tidyr, GGally,
               igraph, ggraph, jsonlite, patchwork, broom, broom.mixed, rstatix, knitr,
               kableExtra)

# ===================================================================
# Runtime configuration
# ===================================================================
load_existing_data <- TRUE
show_plots <- FALSE
save_plots <- TRUE
save_tables <- TRUE
exclude_homecage <- TRUE
analyze_by_halfhour <- TRUE

# Autocorrelation parameter for GAMM residuals.
# rho = 0 disables AR1 correction; set >0 only after residual ACF diagnostics.
gamm_rho <- 0

filter_cc4_late_phases <- TRUE
cc4_max_active_phase <- 2
cc4_max_inactive_phase <- 2

# ===================================================================
# Directory structure
# ===================================================================
working_directory <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability"
saving_directory <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability/new_structure"

raw_data_dir <- file.path(working_directory, "raw_data")
preprocessed_data_dir <- file.path(working_directory, "preprocessed_data")

tables_dir <- file.path(saving_directory, "tables")
plots_dir <- file.path(saving_directory, "plots")
logs_dir <- file.path(saving_directory, "logs")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

combined_tables_dir <- file.path(tables_dir, "combined")
batch_cage_tables_dir <- file.path(tables_dir, "batch_cagechange")
exploration_metrics_dir <- file.path(tables_dir, "exploration_metrics")
metadata_dir <- file.path(combined_tables_dir, "metadata")

dir.create(combined_tables_dir, recursive = TRUE, showWarnings = FALSE)
# Ensure expected combined subfolders exist before saving/loading
dir.create(file.path(combined_tables_dir, "phase_based"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(combined_tables_dir, "halfhour_based"), recursive = TRUE, showWarnings = FALSE)

dir.create(batch_cage_tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(exploration_metrics_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

combined_plots_dir <- file.path(plots_dir, "combined")
batch_cage_plots_dir <- file.path(plots_dir, "batch_cagechange")
misc_plots_dir <- file.path(plots_dir, "miscellaneous")
transition_networks_dir <- file.path(combined_plots_dir, "transition_networks")

dir.create(combined_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(batch_cage_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(misc_plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(transition_networks_dir, recursive = TRUE, showWarnings = FALSE)

# ===================================================================
# Structured output directories
# ===================================================================
pub_dir <- file.path(saving_directory, "publication_ready")
pub_tables_dir <- file.path(pub_dir, "tables")
pub_figures_dir <- file.path(pub_dir, "figures")
pub_figure_panels_dir <- file.path(pub_figures_dir, "single_panels")
pub_multipanel_dir <- file.path(pub_figures_dir, "multi_panel")
pub_supplement_dir <- file.path(pub_dir, "supplementary")

for (dir_path in c(pub_dir, pub_tables_dir, pub_figures_dir,
                  pub_figure_panels_dir, pub_multipanel_dir,
                  pub_supplement_dir)) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

# ===================================================================
# Analysis-ready output structure
# ===================================================================
analysis_dir <- file.path(saving_directory, "analysis_ready")
analysis_raw_loaded_dir <- file.path(analysis_dir, "01_raw_loaded")
analysis_intermediate_phase_dir <- file.path(analysis_dir, "02_intermediate", "phase_based")
analysis_intermediate_halfhour_dir <- file.path(analysis_dir, "02_intermediate", "halfhour_based")
analysis_derived_phase_dir <- file.path(analysis_dir, "03_derived_metrics", "phase_based")
analysis_derived_halfhour_dir <- file.path(analysis_dir, "03_derived_metrics", "halfhour_based")
analysis_model_lmm_dir <- file.path(analysis_dir, "04_model_outputs", "lmm")
analysis_model_gamm_dir <- file.path(analysis_dir, "04_model_outputs", "gamm")
analysis_model_emm_dir <- file.path(analysis_dir, "04_model_outputs", "estimated_marginal_means")
analysis_figures_exploratory_dir <- file.path(analysis_dir, "05_figures", "exploratory")
analysis_figures_model_dir <- file.path(analysis_dir, "05_figures", "model_based")
analysis_figures_multipanel_dir <- file.path(analysis_dir, "05_figures", "multi_panel")
analysis_logs_dir <- file.path(analysis_dir, "06_logs")

for (dir_path in c(
  analysis_raw_loaded_dir,
  analysis_intermediate_phase_dir,
  analysis_intermediate_halfhour_dir,
  analysis_derived_phase_dir,
  analysis_derived_halfhour_dir,
  analysis_model_lmm_dir,
  analysis_model_gamm_dir,
  analysis_model_emm_dir,
  analysis_figures_exploratory_dir,
  analysis_figures_model_dir,
  analysis_figures_multipanel_dir,
  analysis_logs_dir
)) {
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

# ===================================================================
# Plotting and formatting utilities
# ===================================================================
group_cols <- c(
  con = "#457B9D",
  res = "#BFBFBF",
  sus = "#E63946"
)

sex_cols <- c(
  male = "#6F6F91",
  female = "#C9B27C"
)

sem <- function(x) {
  stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
}

save_pub_plot <- function(plot, filename, width = 85, height = 70, units = "mm") {
  ggplot2::ggsave(
    filename = file.path(pub_figure_panels_dir, filename),
    plot = plot,
    width = width,
    height = height,
    units = units,
    dpi = 600,
    bg = "white"
  )
}

save_pub_multipanel <- function(plot, filename, width = 180, height = 120, units = "mm") {
  ggplot2::ggsave(
    filename = file.path(pub_multipanel_dir, filename),
    plot = plot,
    width = width,
    height = height,
    units = units,
    dpi = 600,
    bg = "white"
  )
}

theme_publication <- function(base_size = 7, base_family = "Arial") {
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

# -------------------------------------------------------------------
# Analysis file naming and output helpers
# -------------------------------------------------------------------
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

write_analysis_csv <- function(data, directory, prefix, descriptor, level = NULL) {
  readr::write_csv(
    data,
    file.path(directory, analysis_file(prefix, descriptor, level, ext = "csv"))
  )
}

save_analysis_figure <- function(plot, directory, prefix, descriptor, level = NULL,
                               width = 90, height = 75, units = "mm", ext = "svg") {
  filename <- analysis_file(prefix, descriptor, level, ext = ext)
  ggplot2::ggsave(
    filename = file.path(directory, filename),
    plot = plot,
    width = width,
    height = height,
    units = units,
    dpi = 600,
    bg = "white"
  )
  invisible(filename)
}

make_entropy_summary_table <- function(data, entropy_col, grouping_vars) {
  data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) %>%
    dplyr::summarise(
      n = sum(!is.na(.data[[entropy_col]])),
      mean = mean(.data[[entropy_col]], na.rm = TRUE),
      sd = stats::sd(.data[[entropy_col]], na.rm = TRUE),
      sem = sem(.data[[entropy_col]]),
      median = stats::median(.data[[entropy_col]], na.rm = TRUE),
      q25 = stats::quantile(.data[[entropy_col]], 0.25, na.rm = TRUE),
      q75 = stats::quantile(.data[[entropy_col]], 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      mean_sem = sprintf("%.3f ± %.3f", mean, sem),
      median_iqr = sprintf("%.3f [%.3f, %.3f]", median, q25, q75)
    )
}

# -------------------------------------------------------------------
# Model helper functions for standardized outputs
# -------------------------------------------------------------------
standardize_group_levels <- function(data) {
  data %>%
    dplyr::mutate(
      Group = stringr::str_to_lower(stringr::str_trim(as.character(Group))),
      Sex = stringr::str_to_lower(stringr::str_trim(as.character(Sex))),
      Phase = stringr::str_to_lower(stringr::str_trim(as.character(Phase))),
      CageChange = stringr::str_to_upper(stringr::str_trim(as.character(CageChange))),
      Group = factor(Group, levels = c("con", "res", "sus")),
      Sex = factor(Sex, levels = c("male", "female")),
      Phase = factor(Phase, levels = c("active", "inactive")),
      CageChange = factor(CageChange, levels = c("CC1", "CC2", "CC3", "CC4"))
    )
}

standardize_cage_levels <- function(data) {
  data %>%
    dplyr::mutate(
      Sex = stringr::str_to_lower(stringr::str_trim(as.character(Sex))),
      Phase = stringr::str_to_lower(stringr::str_trim(as.character(Phase))),
      CageChange = stringr::str_to_upper(stringr::str_trim(as.character(CageChange))),
      Sex = factor(Sex, levels = c("male", "female")),
      Phase = factor(Phase, levels = c("active", "inactive")),
      CageChange = factor(CageChange, levels = c("CC1", "CC2", "CC3", "CC4"))
    )
}

format_lmm_fixed_effects <- function(model, model_name) {
  broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE) %>%
    dplyr::mutate(
      model = model_name,
      dplyr::across(where(is.numeric), ~ round(.x, 5))
    ) %>%
    dplyr::select(model, dplyr::everything())
}

format_type3_anova <- function(model, model_name) {
  as.data.frame(anova(model, type = 3)) %>%
    tibble::rownames_to_column("term") %>%
    dplyr::mutate(
      model = model_name,
      dplyr::across(where(is.numeric), ~ round(.x, 5))
    ) %>%
    dplyr::select(model, dplyr::everything())
}

format_emmeans_table <- function(emm_object, table_name) {
  as.data.frame(emm_object) %>%
    dplyr::mutate(
      table = table_name,
      dplyr::across(where(is.numeric), ~ round(.x, 5))
    ) %>%
    dplyr::select(table, dplyr::everything())
}

format_gamm_parametric <- function(model, model_name) {
  as.data.frame(summary(model)$p.table) %>%
    tibble::rownames_to_column("term") %>%
    dplyr::rename(
      estimate = Estimate,
      std_error = `Std. Error`,
      statistic = `t value`,
      p_value = `Pr(>|t|)`
    ) %>%
    dplyr::mutate(
      model = model_name,
      dplyr::across(where(is.numeric), ~ round(.x, 5))
    ) %>%
    dplyr::select(model, dplyr::everything())
}

format_gamm_smooths <- function(model, model_name) {
  as.data.frame(summary(model)$s.table) %>%
    tibble::rownames_to_column("smooth_term") %>%
    dplyr::mutate(
      model = model_name,
      dplyr::across(where(is.numeric), ~ round(.x, 5))
    ) %>%
    dplyr::select(model, dplyr::everything())
}

# -------------------------------------------------------------------
# Sex-stratified model fitting functions
# -------------------------------------------------------------------
fit_sex_specific_animal_lmm <- function(data, sex_value) {
  data_sex <- data %>% dplyr::filter(Sex == sex_value)

  if (nrow(data_sex) == 0) {
    message("Skipping animal entropy LMM for ", sex_value, ": no rows after filtering.")
    return(NULL)
  }

  if (dplyr::n_distinct(data_sex$Group) < 2) {
    message("Skipping animal entropy LMM for ", sex_value, ": fewer than two groups.")
    return(NULL)
  }

  lmerTest::lmer(
    animalEntropy ~ Group * Phase + CageChange + (1 | AnimalID) + (1 | System),
    data = data_sex
  )
}

fit_sex_specific_cage_lmm <- function(data, sex_value) {
  data_sex <- data %>% dplyr::filter(Sex == sex_value)

  if (nrow(data_sex) == 0) {
    message("Skipping cage entropy LMM for ", sex_value, ": no rows after filtering.")
    return(NULL)
  }

  if (dplyr::n_distinct(data_sex$Phase) < 2) {
    message("Skipping cage entropy LMM for ", sex_value, ": fewer than two phases.")
    return(NULL)
  }

  lmerTest::lmer(
    CageEntropy ~ Phase + CageChange + (1 | System),
    data = data_sex
  )
}

fit_sex_specific_animal_gamm <- function(data, sex_value, rho = gamm_rho) {
  data_sex <- data %>%
    dplyr::filter(Sex == sex_value) %>%
    dplyr::arrange(System, AnimalID, CageChange, Phase, ConsecHalfHour) %>%
    dplyr::mutate(
      AnimalID = factor(AnimalID),
      System = factor(System),
      GroupPhase = interaction(Group, Phase, drop = TRUE),
      AR_start = dplyr::row_number() == 1L |
        AnimalID != dplyr::lag(AnimalID) |
        System != dplyr::lag(System) |
        CageChange != dplyr::lag(CageChange) |
        Phase != dplyr::lag(Phase),
      ConsecHalfHour_z = as.numeric(scale(ConsecHalfHour))
    )

  if (nrow(data_sex) == 0) {
    message("Skipping animal entropy GAMM for ", sex_value, ": no rows after filtering.")
    return(NULL)
  }

  if (dplyr::n_distinct(data_sex$Group) < 2) {
    message("Skipping animal entropy GAMM for ", sex_value, ": fewer than two groups.")
    return(NULL)
  }

  if (dplyr::n_distinct(data_sex$GroupPhase) < 2) {
    message("Skipping animal entropy GAMM for ", sex_value, ": fewer than two group-phase combinations.")
    return(NULL)
  }

  mgcv::bam(
    animalEntropy ~ Group * Phase + CageChange +
      s(ConsecHalfHour_z, by = GroupPhase, k = 8) +
      s(AnimalID, bs = "re") +
      s(System, bs = "re"),
    data = data_sex,
    method = "fREML",
    discrete = TRUE,
    rho = rho,
    AR.start = data_sex$AR_start
  )
}

create_batch_cage_folders <- function(batches, cageChanges) {
  for (batch in batches) {
    for (cageChange in cageChanges) {
      phase_tab <- file.path(batch_cage_tables_dir, batch, cageChange, "phase_based")
      halfhour_tab <- file.path(batch_cage_tables_dir, batch, cageChange, "halfhour_based")
      dir.create(phase_tab, recursive = TRUE, showWarnings = FALSE)
      dir.create(halfhour_tab, recursive = TRUE, showWarnings = FALSE)
      phase_plot <- file.path(batch_cage_plots_dir, batch, cageChange, "phase_based")
      halfhour_plot <- file.path(batch_cage_plots_dir, batch, cageChange, "halfhour_based")
      exploration_plot <- file.path(batch_cage_plots_dir, batch, cageChange, "exploration")
      dir.create(phase_plot, recursive = TRUE, showWarnings = FALSE)
      dir.create(halfhour_plot, recursive = TRUE, showWarnings = FALSE)
      dir.create(exploration_plot, recursive = TRUE, showWarnings = FALSE)
    }
  }
}

batches <- c("B1", "B2", "B3", "B4", "B5", "B6")
cageChanges <- c("CC1", "CC2", "CC3", "CC4")
create_batch_cage_folders(batches, cageChanges)

# ===================================================================
# Initialize logging
# ===================================================================
processing_start_time <- Sys.time()
log_file <- file.path(logs_dir, "processing_log.txt")
if (file.exists(log_file)) file.remove(log_file)

log_message <- function(message_text, batch = NULL, cageChange = NULL, phase = NULL) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  context <- paste(Filter(Negate(is.null), list(batch, cageChange, phase)), collapse = " | ")
  full_message <- paste0("[", timestamp, "] ", context, ": ", message_text)
  message(full_message)
  write(full_message, file = log_file, append = TRUE)
}
log_message("Starting Shannon entropy analysis pipeline")

# ===================================================================
# Load custom functions and animal group definitions
# ===================================================================
source(file.path(working_directory, "E9_SIS_AnimalPos-functions.R"))

sus_animals <- readLines(file.path(raw_data_dir, "sus_animals.csv"))
con_animals <- readLines(file.path(raw_data_dir, "con_animals.csv"))

# ===================================================================
# Option 1: Load preprocessed data
# ===================================================================
if (load_existing_data) {
  message("=======================================================================")
  message("LOADING EXISTING PROCESSED DATA")
  message("=======================================================================")
  
  files_to_load <- list(
    cagePosProb = file.path(combined_tables_dir, "phase_based", "all_batch_CC_cagePosProb.csv"),
    cagePosEntropy = file.path(combined_tables_dir, "phase_based", "all_batch_CC_cagePosEntropy.csv"),
    animalPosEntropy = file.path(combined_tables_dir, "phase_based", "all_batch_CC_animalPosEntropy.csv")
  )
  
  if (analyze_by_halfhour) {
    files_to_load$cagePosEntropy_halfhour <- file.path(combined_tables_dir, "halfhour_based", "all_batch_CC_cagePosEntropy_halfhour.csv")
    files_to_load$animalPosEntropy_halfhour <- file.path(combined_tables_dir, "halfhour_based", "all_batch_CC_animalPosEntropy_halfhour.csv")
  }
  
  all_files_exist <- all(file.exists(unlist(files_to_load)))
  
  if (all_files_exist) {
    message("✓ All required files found. Loading data...")
    
    all_cagePosProb <- read_csv(files_to_load$cagePosProb, show_col_types = FALSE)
    message("   ✓ Loaded: all_batch_CC_cagePosProb.csv")
    
    all_cagePosEntropy <- read_csv(files_to_load$cagePosEntropy, show_col_types = FALSE)
    message("   ✓ Loaded: all_batch_CC_cagePosEntropy.csv")
    
    all_animalPosEntropy <- read_csv(files_to_load$animalPosEntropy, show_col_types = FALSE)
    message("   ✓ Loaded: all_batch_CC_animalPosEntropy.csv")
    
    if (analyze_by_halfhour) {
      all_cagePosEntropy_halfhour <- read_csv(files_to_load$cagePosEntropy_halfhour, show_col_types = FALSE)
      message("   ✓ Loaded: all_batch_CC_cagePosEntropy_halfhour.csv")
    
      all_animalPosEntropy_halfhour <- read_csv(files_to_load$animalPosEntropy_halfhour, show_col_types = FALSE)
      message("   ✓ Loaded: all_batch_CC_animalPosEntropy_halfhour.csv")
    }
    
    message("=======================================================================")
    message("DATA LOADING COMPLETE")
    message("=======================================================================")
  
  } else {
    message("✗ Some files not found. Switching to processing mode...")
    message("   Missing files:")
    for (f in unlist(files_to_load)[!file.exists(unlist(files_to_load))]) {
      message("   - ", basename(f))
    }
    load_existing_data <- FALSE
  }
}

# ===================================================================
# Option 2: Process raw data
# ===================================================================
if (!load_existing_data) {
  message("=======================================================================")
  message("STARTING DATA PROCESSING")
  message("=======================================================================")
  
  # Initialize combined empty tibbles for cumulative results
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
    animalEntropy = numeric()
  )
  
  if (analyze_by_halfhour) {
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
      animalEntropy = numeric()
    )
  }
  
  for (batch in batches) {
  sex <- ifelse(batch %in% c("B3", "B4", "B6"), "female", "male")
  
  for (cageChange in cageChanges) {
    message("=======================================================================")
    message("PROCESSING: ", batch, " ", cageChange)
    message("=======================================================================")
    
    # Define paths for tables and plots
    current_tables_dir_phase <- file.path(batch_cage_tables_dir, batch, cageChange, "phase_based")
    current_tables_dir_halfhour <- file.path(batch_cage_tables_dir, batch, cageChange, "halfhour_based")
    current_plots_dir_phase <- file.path(plots_dir, "batch_cagechange", batch, cageChange, "phase_based")
    current_plots_dir_halfhour <- file.path(plots_dir, "batch_cagechange", batch, cageChange, "halfhour_based")
    dir.create(current_tables_dir_phase, recursive = TRUE, showWarnings = FALSE)
    dir.create(current_tables_dir_halfhour, recursive = TRUE, showWarnings = FALSE)
    dir.create(current_plots_dir_phase, recursive = TRUE, showWarnings = FALSE)
    dir.create(current_plots_dir_halfhour, recursive = TRUE, showWarnings = FALSE)
    
    # Initialize empty tibbles with columns defined
    cagePosProb <- tibble(
      Batch = character(),
      System = character(),
      CageChange = character(),
      Phase = character(),
      Position = integer(),
      Probability = numeric()
    )
    
    cagePosEntropy <- tibble(
      Batch = character(),
      Sex = character(),
      System = character(),
      CageChange = character(),
      Phase = character(),
      CageEntropy = numeric()
    )
    
    animalPosEntropy <- tibble(
      Batch = character(),
      Sex = character(),
      System = character(),
      CageChange = character(),
      Phase = character(),
      AnimalID = character(),
      animalEntropy = numeric()
    )
    
    if (analyze_by_halfhour) {
      cagePosEntropy_halfhour <- tibble(
        Batch = character(),
        Sex = character(),
        System = character(),
        CageChange = character(),
        HalfHour = integer(),
        CageEntropy = numeric()
      )
      
      animalPosEntropy_halfhour <- tibble(
        Batch = character(),
        Sex = character(),
        System = character(),
        CageChange = character(),
        HalfHour = integer(),
        AnimalID = character(),
        animalEntropy = numeric()
      )
    }
    
    # Load preprocessed data
    filename <- paste0("E9_SIS_", batch, "_", cageChange, "_AnimalPos")
    csvFilePath <- file.path(preprocessed_data_dir, paste0(filename, "_preprocessed.csv"))
    if (!file.exists(csvFilePath)) {
      message("✗ File not found: ", csvFilePath)
      message("   Skipping ", batch, " ", cageChange)
      next
    }
    
    data_preprocessed <- read_delim(csvFilePath, delim = ",", show_col_types = FALSE) %>% as_tibble()
    message("✓ Loaded preprocessed data: ", nrow(data_preprocessed), " rows")
    
    unique_systems <- str_sort(unique(data_preprocessed$System))
    phases <- c("Active", "Inactive")
    active_phases <- unique(data_preprocessed$ConsecActive)
    inactive_phases <- unique(data_preprocessed$ConsecInactive)
    
    # Filter late phases in CC4
    if (cageChange == "CC4") {
      active_phases <- active_phases[active_phases <= cc4_max_active_phase]
      inactive_phases <- inactive_phases[inactive_phases <= cc4_max_inactive_phase]
    }
    
    # Phase-based entropy and occupancy probability calculation
    for (system_id in unique_systems) {
      system_data <- data_preprocessed %>% filter(System == system_id) %>% as_tibble()
      animal_ids <- unique(system_data$AnimalID)
      system_complete <- length(animal_ids) >= 4
      while (length(animal_ids) < 4) animal_ids <- append(animal_ids, NA)
      
      for (phase in phases) {
        phase_nums <- if (phase == "Active") active_phases else inactive_phases
        
        for (phaseN in phase_nums) {
          data_phase <- system_data %>% 
            filter(ConsecActive == ifelse(phase == "Active", phaseN, 0)) %>% 
            filter(ConsecInactive == ifelse(phase == "Inactive", phaseN, 0)) %>% 
            as_tibble()
          
          if (nrow(data_phase) == 0) next
          
          animal_list <- list(
            animal_1 = list(name = "", time = "", position = 0),
            animal_2 = list(name = "", time = "", position = 0),
            animal_3 = list(name = "", time = "", position = 0),
            animal_4 = list(name = "", time = "", position = 0),
            data_temp = list(elapsed_seconds = 0, current_row = 0)
          )
          
          cage_position_probability <- list(
            c(1, 0, 0, 0), c(2, 0, 0, 0), c(3, 0, 0, 0), c(4, 0, 0, 0),
            c(5, 0, 0, 0), c(6, 0, 0, 0), c(7, 0, 0, 0), c(8, 0, 0, 0)
          )
          
          animal_position_probability <- tibble(AnimalID = rep(animal_ids, each = 8),
                                                Position = rep(1:8, length.out = 32),
                                                Seconds = 0,
                                                SumPercentage = 0,
                                                Prob = 0)
          
          animal_list <- initialize_animal_positions(animal_ids, data_phase, animal_list)
          initial_time <- animal_list[[1]][[2]]
          current_row <- 5
          total_rows <- nrow(data_phase) + 1
          
          while (current_row != total_rows && current_row < total_rows) {
            previous_animal_positions <- animal_list
            animal_list <- update_animal_list(animal_ids, animal_list, data_phase, initial_time, current_row)
            elapsed_seconds <- animal_list[["data_temp"]][["elapsed_seconds"]]
            
            if (system_complete) {
              cage_position_probability <- update_cage_position_probability(previous_animal_positions, animal_list, cage_position_probability, elapsed_seconds)
            }
            
            animal_position_probability <- update_animal_position_probability(previous_animal_positions, animal_list, animal_position_probability, elapsed_seconds)
            current_row <- animal_list[["data_temp"]][["current_row"]]
            initial_time <- animal_list[[1]][[2]]
          }
          
          for (i in 1:8) {
            if (cage_position_probability[[i]][[2]] > 0) {
              cage_position_probability[[i]][[4]] <- cage_position_probability[[i]][[3]] / cage_position_probability[[i]][[2]]
            } else {
              cage_position_probability[[i]][[4]] <- 0
            }
          }
          
          animal_position_probability <- animal_position_probability %>% 
            mutate(Prob = ifelse(Seconds > 0, SumPercentage / Seconds, 0))
          
          if (system_complete) {
            for (i in 1:8) {
              p <- paste0(substr(phase, 1, 1), phaseN)
              cagePosProb <- cagePosProb %>% 
                add_row(Batch = batch, System = system_id, CageChange = cageChange,
                        Phase = p, Position = i, Probability = cage_position_probability[[i]][[4]])
            }
            
            cage_prob_vec <- sapply(cage_position_probability, function(x) x[4])
            cage_prob_vec <- cage_prob_vec[!is.na(cage_prob_vec) & !is.nan(cage_prob_vec)]
            if (length(cage_prob_vec) > 0 && sum(cage_prob_vec) > 0) {
              cage_prob_vec <- cage_prob_vec / sum(cage_prob_vec)
              cage_shannon_entropy <- calc_shannon_entropy(cage_prob_vec)
            } else {
              cage_shannon_entropy <- NA
            }
            
            cagePosEntropy <- cagePosEntropy %>% 
              add_row(Batch = batch, Sex = sex, System = system_id, CageChange = cageChange,
                      Phase = p, CageEntropy = cage_shannon_entropy)
          }
          
          for (animal in animal_ids) {
            if (is.na(animal)) next
            
            animal_prob_vec <- animal_position_probability %>% 
              filter(AnimalID == animal) %>% pull(Prob)
            
            animal_prob_vec <- animal_prob_vec[!is.na(animal_prob_vec) & !is.nan(animal_prob_vec)]
            
            if (length(animal_prob_vec) > 0 && sum(animal_prob_vec) > 0) {
              animal_prob_vec <- animal_prob_vec / sum(animal_prob_vec)
              animal_shannon_entropy <- calc_shannon_entropy(animal_prob_vec)
            } else {
              animal_shannon_entropy <- NA
            }
            
            animalPosEntropy <- animalPosEntropy %>% add_row(Batch=batch, Sex=sex, System=system_id,
                                                            CageChange=cageChange, Phase=p, AnimalID=animal, animalEntropy=animal_shannon_entropy)
          }
        }
      }
    }
    
    # -------------------------------------------------------------------
    # Half-hour processing (CC4 filtering applied prior to iteration)
    # -------------------------------------------------------------------
    if (analyze_by_halfhour) {
      for (system_id in unique_systems) {
        system_data <- data_preprocessed %>% 
          filter(System == system_id) %>% 
          as_tibble()

        animal_ids <- unique(system_data$AnimalID)
        system_complete <- length(animal_ids) >= 4
        while (length(animal_ids) < 4) animal_ids <- append(animal_ids, NA)

        # ---- CC4-consistent half-hour filtering ----
        valid_halfhours <- system_data$HalfHoursElapsed

        if (cageChange == "CC4") {
          valid_halfhours <- system_data %>%
            filter(
              (ConsecActive > 0 & ConsecActive <= cc4_max_active_phase) |
                (ConsecInactive > 0 & ConsecInactive <= cc4_max_inactive_phase)
            ) %>%
            pull(HalfHoursElapsed) %>%
            unique()
        }

        for (halfhour in sort(valid_halfhours)) {

          data_halfhour <- system_data %>%
            filter(HalfHoursElapsed == halfhour) %>%
            as_tibble()

          if (nrow(data_halfhour) == 0) next

          animal_list <- list(
            animal_1 = list(name = "", time = "", position = 0),
            animal_2 = list(name = "", time = "", position = 0),
            animal_3 = list(name = "", time = "", position = 0),
            animal_4 = list(name = "", time = "", position = 0),
            data_temp = list(elapsed_seconds = 0, current_row = 0)
          )

          cage_position_probability <- replicate(8, c(NA, 0, 0, 0), simplify = FALSE)
          for (i in 1:8) cage_position_probability[[i]][1] <- i

          animal_position_probability <- tibble(
            AnimalID = rep(animal_ids, each = 8),
            Position = rep(1:8, length.out = 32),
            Seconds = 0,
            SumPercentage = 0,
            Prob = 0
          )

          animal_list <- initialize_animal_positions(animal_ids, data_halfhour, animal_list)
          initial_time <- animal_list[[1]][[2]]
          current_row <- 5
          total_rows <- nrow(data_halfhour) + 1

          while (current_row < total_rows) {
            previous_positions <- animal_list
            animal_list <- update_animal_list(
              animal_ids, animal_list, data_halfhour, initial_time, current_row
            )

            elapsed_seconds <- animal_list$data_temp$elapsed_seconds

            if (system_complete) {
              cage_position_probability <- update_cage_position_probability(
                previous_positions, animal_list,
                cage_position_probability, elapsed_seconds
              )
            }

            animal_position_probability <- update_animal_position_probability(
              previous_positions, animal_list,
              animal_position_probability, elapsed_seconds
            )

            current_row <- animal_list$data_temp$current_row
            initial_time <- animal_list[[1]][[2]]
          }

          for (i in 1:8) {
            cage_position_probability[[i]][4] <-
              ifelse(cage_position_probability[[i]][2] > 0,
                    cage_position_probability[[i]][3] /
                      cage_position_probability[[i]][2], 0)
          }

          animal_position_probability <- animal_position_probability %>%
            mutate(Prob = ifelse(Seconds > 0, SumPercentage / Seconds, 0))

          if (system_complete) {
            cage_prob_vec <- sapply(cage_position_probability, `[`, 4)
            cage_prob_vec <- cage_prob_vec[!is.na(cage_prob_vec) & !is.nan(cage_prob_vec)]

            cage_entropy <- if (sum(cage_prob_vec) > 0) {
              calc_shannon_entropy(cage_prob_vec / sum(cage_prob_vec))
            } else NA

            cagePosEntropy_halfhour <- cagePosEntropy_halfhour %>%
              add_row(
                Batch = batch, Sex = sex, System = system_id,
                CageChange = cageChange, HalfHour = halfhour,
                CageEntropy = cage_entropy
              )
          }

          for (animal in animal_ids) {
            if (is.na(animal)) next

            prob_vec <- animal_position_probability %>%
              filter(AnimalID == animal) %>%
              pull(Prob)

            entropy <- if (sum(prob_vec, na.rm = TRUE) > 0) {
              calc_shannon_entropy(prob_vec / sum(prob_vec))
            } else NA

            animalPosEntropy_halfhour <- animalPosEntropy_halfhour %>%
              add_row(
                Batch = batch, Sex = sex, System = system_id,
                CageChange = cageChange, HalfHour = halfhour,
                AnimalID = animal, animalEntropy = entropy
              )
          }
        }
      }
    }
    
    # Assign experimental group labels
    animalPosEntropy <- animalPosEntropy %>%
      mutate(Group = ifelse(AnimalID %in% sus_animals, "sus",
                            ifelse(AnimalID %in% con_animals, "con", "res")))
    if (analyze_by_halfhour) {
      animalPosEntropy_halfhour <- animalPosEntropy_halfhour %>%
        mutate(Group = ifelse(AnimalID %in% sus_animals, "sus",
                              ifelse(AnimalID %in% con_animals, "con", "res")))
    }
    
    # Generate phase-based entropy plots per system
for (system_id in unique_systems) {
  for (phase in phases) {
    # Insert diagnostics here:
    message(sprintf("Batch %s Cage %s: animalPosEntropy rows: %d", batch, cageChange, nrow(animalPosEntropy)))
    message(sprintf("Systems: %s", paste(unique(animalPosEntropy$System), collapse = ", ")))
    message(sprintf("Phases: %s", paste(unique(animalPosEntropy$Phase), collapse = ", ")))
    phase_code <- ifelse(phase == "Active", "A", "I")
    plot_folder <- file.path(saving_directory, "plots", "batch_cagechange", batch, cageChange, "phase_based")
    dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)
    
    plot_data <- animalPosEntropy %>% filter(System == system_id, Phase == phase_code)

    message(sprintf("Batch %s Cage %s: system=%s phase=%s plot_data rows: %d",
                    batch, cageChange, system_id, phase_code, nrow(plot_data)))
    
    if (nrow(plot_data) == 0) next
    
    p_phase_plot <- ggplot(plot_data, aes(x = Group, y = animalEntropy, color = Group)) +
      geom_jitter(width = 0.2, size = 3, alpha = 0.7) +
      stat_summary(fun = mean, geom = "point", shape = 18, size = 5, color = "black") +
      scale_color_manual(values = c(sus = "#E63946", res = "grey60", con = "#457B9D")) +
      labs(title = paste("Animal Entropy -", batch, cageChange, system_id, phase),
           x = "Group", y = "Shannon Entropy") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 18),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12),
            legend.position = "none")
    
    ggsave(filename = file.path(plot_folder, paste0("animal_entropy_", system_id, "_", phase_code, ".svg")), 
           plot = p_phase_plot, width = 10, height = 7, dpi = 300)
  }
}

    # Generate half-hour entropy plots per system
# -------------------------------------------------------------------
# Plot half-hour animal entropy (FIXED: uses half-hour data)
# -------------------------------------------------------------------
if (analyze_by_halfhour) {
  for (system_id in unique_systems) {

    plot_folder <- file.path(
      saving_directory, "plots",
      "batch_cagechange", batch, cageChange, "halfhour_based"
    )
    dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)

    halfhours <- sort(unique(
      animalPosEntropy_halfhour %>%
        filter(System == system_id, CageChange == cageChange) %>%
        pull(HalfHour)
    ))

    for (halfhour in halfhours) {

      plot_data <- animalPosEntropy_halfhour %>%
        filter(
          System == system_id,
          CageChange == cageChange,
          HalfHour == halfhour,
          !is.na(animalEntropy),
          !is.na(Group)
        )

      if (nrow(plot_data) == 0) next

      p <- ggplot(plot_data, aes(x = Group, y = animalEntropy, color = Group)) +
        geom_jitter(width = 0.2, size = 3, alpha = 0.7) +
        stat_summary(fun = mean, geom = "point",
                     shape = 18, size = 5, color = "black") +
        scale_color_manual(
          values = c(sus = "#E63946", res = "grey60", con = "#457B9D")
        ) +
        labs(
          title = paste(
            "Animal Entropy Over Time (Half-Hour)",
            batch, cageChange, system_id,
            paste("HH", halfhour)
          ),
          x = "Group",
          y = "Shannon Entropy"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 18),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.position = "none"
        )

      ggsave(
        filename = file.path(
          plot_folder,
          paste0("animal_entropy_", system_id, "_HH", halfhour, ".svg")
        ),
        plot = p, width = 10, height = 7, dpi = 300
      )
    }
  }
}


    # Append batch-level results to combined datasets
    all_cagePosProb <- bind_rows(all_cagePosProb, cagePosProb)
    all_cagePosEntropy <- bind_rows(all_cagePosEntropy, cagePosEntropy)
    all_animalPosEntropy <- bind_rows(all_animalPosEntropy, animalPosEntropy)
    if (analyze_by_halfhour) {
      all_cagePosEntropy_halfhour <- bind_rows(all_cagePosEntropy_halfhour, cagePosEntropy_halfhour)
      all_animalPosEntropy_halfhour <- bind_rows(all_animalPosEntropy_halfhour, animalPosEntropy_halfhour)
    }

    # Save batch-level tables
    write.csv(cagePosProb, file=file.path(current_tables_dir_phase, paste0(batch, "_", cageChange, "_cagePosProb.csv")), row.names=FALSE)
    write.csv(cagePosEntropy, file=file.path(current_tables_dir_phase, paste0(batch, "_", cageChange, "_cagePosEntropy.csv")), row.names=FALSE)
    write.csv(animalPosEntropy, file=file.path(current_tables_dir_phase, paste0(batch, "_", cageChange, "_animalPosEntropy.csv")), row.names=FALSE)
    if (analyze_by_halfhour) {
      write.csv(cagePosEntropy_halfhour, file=file.path(current_tables_dir_halfhour, paste0(batch, "_", cageChange, "_cagePosEntropy_halfhour.csv")), row.names=FALSE)
      write.csv(animalPosEntropy_halfhour, file=file.path(current_tables_dir_halfhour, paste0(batch, "_", cageChange, "_animalPosEntropy_halfhour.csv")), row.names=FALSE)
    }
    message("Finished ", batch, " ", cageChange)
  }
}
  
  # Save combined results after all batches/cage changes processed
  if (save_tables == TRUE) {
    message("=======================================================================")
    message("SAVING COMBINED RESULTS")
    message("=======================================================================")
    
    write.csv(all_cagePosProb, file = file.path(combined_tables_dir, "phase_based", "all_batch_CC_cagePosProb.csv"), row.names = FALSE)
    message("   ✓ Saved: all_batch_CC_cagePosProb.csv")
    
    write.csv(all_cagePosEntropy, file = file.path(combined_tables_dir, "phase_based", "all_batch_CC_cagePosEntropy.csv"), row.names = FALSE)
    message("   ✓ Saved: all_batch_CC_cagePosEntropy.csv")
    
    write.csv(all_animalPosEntropy, file = file.path(combined_tables_dir, "phase_based", "all_batch_CC_animalPosEntropy.csv"), row.names = FALSE)
    message("   ✓ Saved: all_batch_CC_animalPosEntropy.csv")
    
    if(analyze_by_halfhour) {
      write.csv(all_cagePosEntropy_halfhour, file = file.path(combined_tables_dir, "halfhour_based", "all_batch_CC_cagePosEntropy_halfhour.csv"), row.names = FALSE)
      message("   ✓ Saved: all_batch_CC_cagePosEntropy_halfhour.csv")
      
      write.csv(all_animalPosEntropy_halfhour, file = file.path(combined_tables_dir, "halfhour_based", "all_batch_CC_animalPosEntropy_halfhour.csv"), row.names = FALSE)
      message("   ✓ Saved: all_batch_CC_animalPosEntropy_halfhour.csv")
    }
    message("      Location: ", combined_tables_dir)
  }
  
  message("=======================================================================")
  message("DATA PROCESSING COMPLETE")
  message("=======================================================================")
}

# ===================================================================
# Consecutive entropy processing (phase-based)
# ===================================================================
message("=======================================================================")
message("PREPROCESSING CONSECUTIVE ENTROPY DATA (PHASE-BASED)")
message("=======================================================================")

consec_cage_entropy <- all_cagePosEntropy
consec_animal_entropy <- all_animalPosEntropy

entropy_tibble_names <- c("consec_cage_entropy", "consec_animal_entropy")

for (i in seq_along(entropy_tibble_names)) {
  name <- entropy_tibble_names[i]
  entropy_list <- if (name == "consec_cage_entropy") consec_cage_entropy else consec_animal_entropy
  
  entropy_list[c("Phase", "Consec")] <- str_split_fixed(entropy_list$Phase, '', 2)
  
  entropy_list <- entropy_list %>%
    mutate(ConsecActive = ifelse(Phase == "A", as.numeric(Consec), 0)) %>%
    mutate(ConsecInactive = ifelse(Phase == "I", as.numeric(Consec) - 1, 0)) %>%
    mutate(Phase = ifelse(Phase == "A", "active", "inactive"))
  
  if (filter_cc4_late_phases) {
    entropy_list <- entropy_list %>%
      filter(!(CageChange == "CC4" & ConsecActive > cc4_max_active_phase)) %>%
      filter(!(CageChange == "CC4" & ConsecInactive > cc4_max_inactive_phase))
    message("✓ Filtered CC4 phases: keeping ≤ A", cc4_max_active_phase, " and I", cc4_max_inactive_phase)
  }
  
  if (name == "consec_animal_entropy") {
    entropy_list <- entropy_list %>%
      mutate(Group = ifelse(AnimalID %in% sus_animals, "sus",
                            ifelse(AnimalID %in% con_animals, "con", "res")))
  }
  
  max_consecAct <- 0
  max_consecInact <- 0

  for (change in c("CC1", "CC2", "CC3", "CC4")) {
    if (change != "CC1") {
      entropy_list <- entropy_list %>%
        mutate(ConsecActive = ifelse(CageChange == change & Phase == "active", ConsecActive + max_consecAct, ConsecActive)) %>%
        mutate(ConsecInactive = ifelse(CageChange == change & Phase == "inactive", ConsecInactive + max_consecInact, ConsecInactive))
    }
    max_consecAct <- entropy_list %>% filter(CageChange == change) %>% pull(ConsecActive) %>% unique() %>% max()
    max_consecInact <- entropy_list %>% filter(CageChange == change) %>% pull(ConsecInactive) %>% unique() %>% max()
  }
  
  if (name == "consec_animal_entropy") {
    entropy_list <- entropy_list[c("CageChange", "Batch", "System", "AnimalID", "Sex", "Group", "Phase", "ConsecActive", "ConsecInactive", "animalEntropy")]
    consec_animal_entropy <- entropy_list
  } else {
    entropy_list <- entropy_list[c("CageChange", "Batch", "System", "Sex", "Phase", "ConsecActive", "ConsecInactive", "CageEntropy")]
    consec_cage_entropy <- entropy_list
  }
}

# ===================================================================
# Consecutive entropy processing (half-hour resolution with phase alignment)
# ===================================================================
if (analyze_by_halfhour) {
  message("=======================================================================")
  message("PREPROCESSING CONSECUTIVE ENTROPY DATA (HALF-HOUR) WITH PHASE COUNTS")
  message("=======================================================================")

  consec_cage_entropy_halfhour <- all_cagePosEntropy_halfhour
  consec_animal_entropy_halfhour <- all_animalPosEntropy_halfhour

  # Ensure we have preprocessed data to extract per-halfhour phase counts
  data_for_filter <- NULL
  if (exists("data_preprocessed") && !is.null(data_preprocessed) && nrow(data_preprocessed) > 0) {
    data_for_filter <- data_preprocessed
  } else {
    all_preprocessed <- tibble()
    for (batch in batches) {
      for (cageChange in cageChanges) {
        filename <- paste0("E9_SIS_", batch, "_", cageChange, "_AnimalPos")
        csvFilePath <- file.path(preprocessed_data_dir, paste0(filename, "_preprocessed.csv"))
        if (file.exists(csvFilePath)) {
          temp_data <- read_delim(csvFilePath, delim = ",", show_col_types = FALSE) %>%
            as_tibble() %>%
            mutate(
              Batch = if ("Batch" %in% colnames(.)) Batch else batch,
              CageChange = if ("CageChange" %in% colnames(.)) CageChange else cageChange
            )
          all_preprocessed <- bind_rows(all_preprocessed, temp_data)
        }
      }
    }
    if (nrow(all_preprocessed) > 0) {
      data_for_filter <- all_preprocessed
      message("   ✓ Loaded combined preprocessed data with ", nrow(data_for_filter), " rows for phase lookup")
    } else {
      message("   Warning: No preprocessed CSVs found; phase counts cannot be joined to half-hour tables.")
    }
  }

  # Build phase lookup: per Batch/System/CageChange/HalfHour get ConsecActive & ConsecInactive
  phase_lookup <- NULL
  if (!is.null(data_for_filter)) {
    phase_lookup <- data_for_filter %>%
      mutate(HalfHour = HalfHoursElapsed) %>%
      group_by(Batch, System, CageChange, HalfHour) %>%
      summarise(
        ConsecActive = if (all(is.na(ConsecActive))) NA_real_ else max(ConsecActive, na.rm = TRUE),
        ConsecInactive = if (all(is.na(ConsecInactive))) NA_real_ else max(ConsecInactive, na.rm = TRUE),
        .groups = "drop"
      )
  }

  # Add ConsecHalfHour and (if available) phase counts to entropy tables
  for (name in c("consec_cage_entropy_halfhour", "consec_animal_entropy_halfhour")) {
    max_halfhour_offset <- 0
    entropy_tbl <- if (name == "consec_cage_entropy_halfhour") consec_cage_entropy_halfhour else consec_animal_entropy_halfhour

    # add Group for animal table
    if (name == "consec_animal_entropy_halfhour") {
      entropy_tbl <- entropy_tbl %>%
        mutate(Group = ifelse(AnimalID %in% sus_animals, "sus",
                              ifelse(AnimalID %in% con_animals, "con", "res")))
    }

    # ensure HalfHour column name matches
    if (!"HalfHour" %in% colnames(entropy_tbl) && "HalfHoursElapsed" %in% colnames(entropy_tbl)) {
      entropy_tbl <- entropy_tbl %>% rename(HalfHour = HalfHoursElapsed)
    }

    entropy_tbl <- entropy_tbl %>% mutate(ConsecHalfHour = HalfHour)

    # perform cumulative offset across cage changes so ConsecHalfHour is continuous
    for (change in c("CC1", "CC2", "CC3", "CC4")) {
      if (change != "CC1") {
        entropy_tbl <- entropy_tbl %>%
          mutate(ConsecHalfHour = ifelse(CageChange == change, HalfHour + max_halfhour_offset, ConsecHalfHour))
      }
      current_max <- entropy_tbl %>% filter(CageChange == change) %>% pull(ConsecHalfHour) %>% max(na.rm = TRUE)
      if (is.finite(current_max)) max_halfhour_offset <- current_max
    }

    # join phase counts if available
    if (!is.null(phase_lookup)) {
      entropy_tbl <- entropy_tbl %>%
        left_join(phase_lookup, by = c("Batch", "System", "CageChange", "HalfHour"))
    } else {
      entropy_tbl <- entropy_tbl %>%
        mutate(ConsecActive = NA_real_, ConsecInactive = NA_real_)
    }

    # infer Phase for the half-hour (active if ConsecActive>0 else inactive)
    entropy_tbl <- entropy_tbl %>%
      mutate(Phase = ifelse(!is.na(ConsecActive) & ConsecActive > 0, "active",
                            ifelse(!is.na(ConsecInactive) & ConsecInactive > 0, "inactive", NA_character_)))

    # select / reorder columns
    if (name == "consec_animal_entropy_halfhour") {
      entropy_tbl <- entropy_tbl %>%
        select(CageChange, Batch, System, AnimalID, Sex, Group, HalfHour, ConsecHalfHour, ConsecActive, ConsecInactive, Phase, animalEntropy)
      consec_animal_entropy_halfhour <- entropy_tbl
    } else {
      entropy_tbl <- entropy_tbl %>%
        select(CageChange, Batch, System, Sex, HalfHour, ConsecHalfHour, ConsecActive, ConsecInactive, Phase, CageEntropy)
      consec_cage_entropy_halfhour <- entropy_tbl
    }
  }

  # Apply CC4 filtering based on phase counts if requested
  if (filter_cc4_late_phases && !is.null(phase_lookup)) {
    message("=======================================================================")
    message("FILTERING CC4 HALF-HOURS USING IMPORTED PHASE COUNTS")
    message("=======================================================================")

    consec_animal_entropy_halfhour <- consec_animal_entropy_halfhour %>%
      filter(!(CageChange == "CC4" & Phase == "active" & !is.na(ConsecActive) & ConsecActive > cc4_max_active_phase)) %>%
      filter(!(CageChange == "CC4" & Phase == "inactive" & !is.na(ConsecInactive) & ConsecInactive > cc4_max_inactive_phase))

    consec_cage_entropy_halfhour <- consec_cage_entropy_halfhour %>%
      filter(!(CageChange == "CC4" & Phase == "active" & !is.na(ConsecActive) & ConsecActive > cc4_max_active_phase)) %>%
      filter(!(CageChange == "CC4" & Phase == "inactive" & !is.na(ConsecInactive) & ConsecInactive > cc4_max_inactive_phase))

    message("✓ CC4 half-hours filtered by imported phase counts")
    message("  Remaining animal half-hour rows: ", nrow(consec_animal_entropy_halfhour))
    message("  Remaining cage half-hour rows: ", nrow(consec_cage_entropy_halfhour))
  } else if (filter_cc4_late_phases) {
    message("   Warning: phase lookup not available; CC4 half-hour filtering by phase skipped.")
  }
}

# ===================================================================
# Save consecutive entropy tables
# ===================================================================
if (save_tables == TRUE) {
  message("Saving consecutive entropy tables...")
  
  write.csv(consec_cage_entropy, 
            file = file.path(combined_tables_dir, "all_batches_all_cageChanges_consec_cage_entropy.csv"), 
            row.names = FALSE)
  message("   ✓ Saved: all_batches_all_cageChanges_consec_cage_entropy.csv")
  
  write.csv(consec_animal_entropy,
            file = file.path(combined_tables_dir, "all_batches_all_cageChanges_consec_animal_entropy.csv"),
            row.names = FALSE)
  message("   ✓ Saved: all_batches_all_cageChanges_consec_animal_entropy.csv")
  
  if (analyze_by_halfhour) {
    write.csv(consec_cage_entropy_halfhour,
              file = file.path(combined_tables_dir, "all_batches_all_cageChanges_consec_cage_entropy_halfhour.csv"),
              row.names = FALSE)
    message("   ✓ Saved: all_batches_all_cageChanges_consec_cage_entropy_halfhour.csv")
    
    write.csv(consec_animal_entropy_halfhour,
              file = file.path(combined_tables_dir, "all_batches_all_cageChanges_consec_animal_entropy_halfhour.csv"),
              row.names = FALSE)
    message("   ✓ Saved: all_batches_all_cageChanges_consec_animal_entropy_halfhour.csv")
  }
  message("    Location: ", combined_tables_dir)
  
  # Also save to analysis_ready structure for clean output
  message("\n  Saving to analysis_ready structure...")
  write_analysis_csv(
    consec_animal_entropy,
    analysis_derived_phase_dir,
    prefix = "derived",
    descriptor = "animal_entropy",
    level = "phase"
  )
  message("   ✓ Saved: derived__animal_entropy__phase.csv")
  
  write_analysis_csv(
    consec_cage_entropy,
    analysis_derived_phase_dir,
    prefix = "derived",
    descriptor = "cage_entropy",
    level = "phase"
  )
  message("   ✓ Saved: derived__cage_entropy__phase.csv")
  
  if (analyze_by_halfhour) {
    write_analysis_csv(
      consec_animal_entropy_halfhour,
      analysis_derived_halfhour_dir,
      prefix = "derived",
      descriptor = "animal_entropy",
      level = "halfhour"
    )
    message("   ✓ Saved: derived__animal_entropy__halfhour.csv")
    
    write_analysis_csv(
      consec_cage_entropy_halfhour,
      analysis_derived_halfhour_dir,
      prefix = "derived",
      descriptor = "cage_entropy",
      level = "halfhour"
    )
    message("   ✓ Saved: derived__cage_entropy__halfhour.csv")
  }
  message("    Location: ", analysis_derived_phase_dir, " / ", analysis_derived_halfhour_dir)
}

# ===================================================================
# Exploration metric computation
# ===================================================================
message("=======================================================================")
message("PROCESSING EXPLORATION METRICS")
message("=======================================================================")

max_position_entropy <- log2(8)
low_entropy_cutoff <- max_position_entropy / 3
moderate_entropy_cutoff <- 2 * max_position_entropy / 3

consec_animal_entropy <- consec_animal_entropy %>%
  mutate(
    explored = !is.na(animalEntropy),
    exploration_entropy = ifelse(is.na(animalEntropy), 0, animalEntropy),
    exploration_category = case_when(
      is.na(animalEntropy) ~ "No exploration",
      animalEntropy == 0 ~ "Minimal (1 position)",
      animalEntropy < low_entropy_cutoff ~ "Low diversity",
      animalEntropy < moderate_entropy_cutoff ~ "Moderate diversity",
      TRUE ~ "High diversity"
    ),
    exploration_category = factor(exploration_category, 
                                  levels = c("No exploration", "Minimal (1 position)", 
                                             "Low diversity", "Moderate diversity", "High diversity"))
  )

if (analyze_by_halfhour) {
  consec_animal_entropy_halfhour <- consec_animal_entropy_halfhour %>%
    mutate(
      explored = !is.na(animalEntropy),
      exploration_entropy = ifelse(is.na(animalEntropy), 0, animalEntropy),
      exploration_category = case_when(
        is.na(animalEntropy) ~ "No exploration",
        animalEntropy == 0 ~ "Minimal (1 position)",
        animalEntropy < low_entropy_cutoff ~ "Low diversity",
        animalEntropy < moderate_entropy_cutoff ~ "Moderate diversity",
        TRUE ~ "High diversity"
      ),
      exploration_category = factor(exploration_category, 
                                    levels = c("No exploration", "Minimal (1 position)", 
                                               "Low diversity", "Moderate diversity", "High diversity"))
    )
}

exploration_summary_phase <- consec_animal_entropy %>%
  group_by(Group, Sex, CageChange, Phase) %>%
  summarise(
    n_observations = n(),
    n_explored = sum(explored),
    pct_exploration = round(100 * n_explored / n_observations, 2),
    mean_entropy_overall = round(mean(exploration_entropy, na.rm = TRUE), 3),
    mean_entropy_when_exploring = round(mean(animalEntropy, na.rm = TRUE), 3),
    sd_entropy_when_exploring = round(sd(animalEntropy, na.rm = TRUE), 3),
    .groups = "drop"
  )

message("Phase-based exploration summary:")
print(exploration_summary_phase)

if (analyze_by_halfhour) {
  exploration_summary_halfhour <- consec_animal_entropy_halfhour %>%
    group_by(Group, Sex, CageChange) %>%
    summarise(
      n_periods = n(),
      n_explored = sum(explored),
      pct_exploration = round(100 * n_explored / n_periods, 2),
      mean_entropy_overall = round(mean(exploration_entropy, na.rm = TRUE), 3),
      mean_entropy_when_exploring = round(mean(animalEntropy, na.rm = TRUE), 3),
      sd_entropy_when_exploring = round(sd(animalEntropy, na.rm = TRUE), 3),
      .groups = 'drop'
    )
  
  message("Half-hour exploration summary:")
  print(exploration_summary_halfhour)
}

animal_exploration_profile <- consec_animal_entropy %>%
  group_by(AnimalID, Group, Sex, Batch) %>%
  summarise(
    total_observations = n(),
    times_explored = sum(explored),
    pct_exploration = round(100 * times_explored / total_observations, 2),
    mean_entropy_overall = round(mean(exploration_entropy, na.rm = TRUE), 3),
    mean_entropy_when_exploring = round(mean(animalEntropy, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>% arrange(desc(pct_exploration))

# ===================================================================
# Save exploration metrics
# ===================================================================
if (save_tables == TRUE) {
  message("Saving enhanced exploration datasets...")
  
  write.csv(consec_animal_entropy, 
            file = file.path(exploration_metrics_dir, "consec_animal_entropy_with_exploration_metrics.csv"),
            row.names = FALSE)
  message("   ✓ Saved: consec_animal_entropy_with_exploration_metrics.csv")
  
  if(analyze_by_halfhour) {
    write.csv(consec_animal_entropy_halfhour, 
              file = file.path(exploration_metrics_dir, "consec_animal_entropy_halfhour_with_exploration_metrics.csv"),
              row.names = FALSE)
    message("   ✓ Saved: consec_animal_entropy_halfhour_with_exploration_metrics.csv")
  }
  
  write.csv(exploration_summary_phase,
            file = file.path(exploration_metrics_dir, "exploration_summary_by_phase.csv"),
            row.names = FALSE)
  message("   ✓ Saved: exploration_summary_by_phase.csv")
  
  if(analyze_by_halfhour) {
    write.csv(exploration_summary_halfhour,
              file = file.path(exploration_metrics_dir, "exploration_summary_by_halfhour.csv"),
              row.names = FALSE)
    message("   ✓ Saved: exploration_summary_by_halfhour.csv")
  }
  
  write.csv(animal_exploration_profile,
            file = file.path(exploration_metrics_dir, "animal_exploration_profiles.csv"),
            row.names = FALSE)
  message("   ✓ Saved: animal_exploration_profiles.csv")
  
  # Also save to analysis_ready structure
  message("\n  Saving exploration summaries to analysis_ready structure...")
  write_analysis_csv(
    exploration_summary_phase,
    analysis_derived_phase_dir,
    prefix = "summary",
    descriptor = "animal_exploration_descriptive_statistics",
    level = "phase"
  )
  message("   ✓ Saved: summary__animal_exploration_descriptive_statistics__phase.csv")
  
  if(analyze_by_halfhour) {
    write_analysis_csv(
      exploration_summary_halfhour,
      analysis_derived_halfhour_dir,
      prefix = "summary",
      descriptor = "animal_exploration_descriptive_statistics",
      level = "halfhour"
    )
    message("   ✓ Saved: summary__animal_exploration_descriptive_statistics__halfhour.csv")
  }
}

# ===================================================================
# Summary tables
# ===================================================================
message("Generating publication-ready summary tables...")

pub_animal_entropy_phase <- make_entropy_summary_table(
  consec_animal_entropy,
  entropy_col = "animalEntropy",
  grouping_vars = c("Group", "Sex", "Phase", "CageChange")
)

pub_cage_entropy_phase <- make_entropy_summary_table(
  consec_cage_entropy,
  entropy_col = "CageEntropy",
  grouping_vars = c("Sex", "Phase", "CageChange")
)

pub_exploration_frequency <- consec_animal_entropy %>%
  dplyr::group_by(Group, Sex, Phase, CageChange) %>%
  dplyr::summarise(
    n_observations = dplyr::n(),
    n_explored = sum(explored, na.rm = TRUE),
    pct_explored = 100 * n_explored / n_observations,
    .groups = "drop"
  ) %>%
  dplyr::mutate(pct_explored = round(pct_explored, 2))

if (analyze_by_halfhour) {
  pub_animal_entropy_halfhour <- make_entropy_summary_table(
    consec_animal_entropy_halfhour,
    entropy_col = "animalEntropy",
    grouping_vars = c("Group", "Sex", "Phase", "CageChange", "ConsecHalfHour")
  )

  pub_cage_entropy_halfhour <- make_entropy_summary_table(
    consec_cage_entropy_halfhour,
    entropy_col = "CageEntropy",
    grouping_vars = c("Sex", "Phase", "CageChange", "ConsecHalfHour")
  )
}

if (save_tables) {
  readr::write_csv(pub_animal_entropy_phase,
                   file.path(pub_tables_dir, "Table_1_animal_entropy_phase_summary.csv"))
  readr::write_csv(pub_cage_entropy_phase,
                   file.path(pub_tables_dir, "Table_2_cage_entropy_phase_summary.csv"))
  readr::write_csv(pub_exploration_frequency,
                   file.path(pub_tables_dir, "Table_3_exploration_frequency_summary.csv"))


  if (analyze_by_halfhour) {
    readr::write_csv(pub_animal_entropy_halfhour,
                     file.path(pub_tables_dir, "Table_S1_animal_entropy_halfhour_summary.csv"))
    readr::write_csv(pub_cage_entropy_halfhour,
                     file.path(pub_tables_dir, "Table_S2_cage_entropy_halfhour_summary.csv"))
  }
}

# ===================================================================
# Model output tables
# ===================================================================
message("Generating publication-ready sex-specific model tables...")

animal_model_data <- consec_animal_entropy %>%
  dplyr::filter(!is.na(animalEntropy), !is.na(Group), !is.na(Sex), !is.na(Phase), !is.na(CageChange)) %>%
  standardize_group_levels()

cage_model_data <- consec_cage_entropy %>%
  dplyr::filter(!is.na(CageEntropy), !is.na(Sex), !is.na(Phase), !is.na(CageChange)) %>%
  standardize_cage_levels()

publication_model_tables <- list()
publication_emmeans_tables <- list()
publication_gamm_tables <- list()

for (sex_value in c("male", "female")) {

  animal_entropy_lmm <- fit_sex_specific_animal_lmm(animal_model_data, sex_value)
  if (!is.null(animal_entropy_lmm)) {
    animal_entropy_emm <- emmeans::emmeans(animal_entropy_lmm, ~ Group | Phase)
    animal_entropy_group_contrasts <- emmeans::contrast(
      animal_entropy_emm, method = "pairwise", adjust = "holm"
    )

    emm_df <- as.data.frame(animal_entropy_emm) %>%
      dplyr::mutate(Sex = sex_value, dplyr::across(where(is.numeric), ~ round(.x, 5)))

    eff_df <- as.data.frame(animal_entropy_group_contrasts, infer = TRUE) %>%
      dplyr::mutate(
        Sex = sex_value,
        contrast = stringr::str_replace_all(contrast, " ", ""),
        dplyr::across(where(is.numeric), ~ round(.x, 5))
      )

    publication_model_tables[[paste0("animal_entropy_lmm_fixed_", sex_value)]] <-
      format_lmm_fixed_effects(animal_entropy_lmm, paste0("animal_entropy_lmm_", sex_value))

    publication_model_tables[[paste0("animal_entropy_lmm_anova_", sex_value)]] <-
      format_type3_anova(animal_entropy_lmm, paste0("animal_entropy_lmm_", sex_value))

    publication_emmeans_tables[[paste0("animal_entropy_emmeans_", sex_value)]] <- emm_df
    publication_emmeans_tables[[paste0("animal_entropy_contrasts_", sex_value)]] <- eff_df
  }

  cage_entropy_lmm <- fit_sex_specific_cage_lmm(cage_model_data, sex_value)
  if (!is.null(cage_entropy_lmm)) {
    cage_entropy_emm <- emmeans::emmeans(cage_entropy_lmm, ~ Phase)
    cage_entropy_contrasts <- emmeans::contrast(
      cage_entropy_emm, method = "pairwise", adjust = "holm"
    )

    publication_model_tables[[paste0("cage_entropy_lmm_fixed_", sex_value)]] <-
      format_lmm_fixed_effects(cage_entropy_lmm, paste0("cage_entropy_lmm_", sex_value))

    publication_model_tables[[paste0("cage_entropy_lmm_anova_", sex_value)]] <-
      format_type3_anova(cage_entropy_lmm, paste0("cage_entropy_lmm_", sex_value))

    publication_emmeans_tables[[paste0("cage_entropy_emmeans_", sex_value)]] <-
      format_emmeans_table(cage_entropy_emm, paste0("cage_entropy_emmeans_", sex_value))

    publication_emmeans_tables[[paste0("cage_entropy_contrasts_", sex_value)]] <-
      format_emmeans_table(cage_entropy_contrasts, paste0("cage_entropy_contrasts_", sex_value))
  }

  if (analyze_by_halfhour) {
    animal_gamm_data <- consec_animal_entropy_halfhour %>%
      dplyr::filter(!is.na(animalEntropy), !is.na(Group), !is.na(Sex), !is.na(Phase), !is.na(ConsecHalfHour)) %>%
      standardize_group_levels()

    gamm_fit <- fit_sex_specific_animal_gamm(animal_gamm_data, sex_value)
    if (!is.null(gamm_fit)) {
      publication_gamm_tables[[paste0("animal_entropy_gamm_param_", sex_value)]] <-
        format_gamm_parametric(gamm_fit, paste0("animal_entropy_gamm_", sex_value))

      publication_gamm_tables[[paste0("animal_entropy_gamm_smooth_", sex_value)]] <-
        format_gamm_smooths(gamm_fit, paste0("animal_entropy_gamm_", sex_value))
    }
  }
}

# Combine EMM tables for plotting
animal_entropy_emm_df <- purrr::map_dfr(
  publication_emmeans_tables[stringr::str_detect(names(publication_emmeans_tables), "animal_entropy_emmeans_")],
  dplyr::bind_rows
)

animal_entropy_effects <- purrr::map_dfr(
  publication_emmeans_tables[stringr::str_detect(names(publication_emmeans_tables), "animal_entropy_contrasts_")],
  dplyr::bind_rows
)

# Save model/EMM/GAMM tables for publication if requested
if (save_tables) {
  purrr::iwalk(publication_model_tables, ~ readr::write_csv(
    .x,
    file.path(pub_tables_dir, paste0("Model_", .y, ".csv"))
  ))

  purrr::iwalk(publication_emmeans_tables, ~ readr::write_csv(
    .x,
    file.path(pub_tables_dir, paste0("EMM_", .y, ".csv"))
  ))

  purrr::iwalk(publication_gamm_tables, ~ readr::write_csv(
    .x,
    file.path(pub_tables_dir, paste0("GAMM_", .y, ".csv"))
  ))
  
  # Also save to analysis_ready structure with consistent naming
  message("\nSaving model outputs to analysis_ready structure...")
  purrr::iwalk(publication_model_tables, ~ write_analysis_csv(
    .x,
    analysis_model_lmm_dir,
    prefix = "model",
    descriptor = .y
  ))
  message("   ✓ Saved LMM model tables")
  
  purrr::iwalk(publication_emmeans_tables, ~ write_analysis_csv(
    .x,
    analysis_model_emm_dir,
    prefix = "emm",
    descriptor = .y
  ))
  message("   ✓ Saved estimated marginal means tables")
  
  purrr::iwalk(publication_gamm_tables, ~ write_analysis_csv(
    .x,
    analysis_model_gamm_dir,
    prefix = "gamm",
    descriptor = .y
  ))
  message("   ✓ Saved GAMM model tables")
}

# ===================================================================
# Run manifest
# ===================================================================
run_manifest <- tibble::tibble(
  field = c(
    "script",
    "run_time",
    "load_existing_data",
    "exclude_homecage",
    "analyze_by_halfhour",
    "gamm_rho",
    "max_position_entropy",
    "low_entropy_cutoff",
    "moderate_entropy_cutoff",
    "filter_cc4_late_phases",
    "cc4_max_active_phase",
    "cc4_max_inactive_phase",
    "n_animal_phase_rows",
    "n_cage_phase_rows",
    "n_animal_halfhour_rows",
    "n_cage_halfhour_rows",
    "working_directory",
    "saving_directory"
  ),
  value = c(
    "E9_SIS_AnimalPos-analyzing-shannon v.1.0.1.r",
    as.character(Sys.time()),
    as.character(load_existing_data),
    as.character(exclude_homecage),
    as.character(analyze_by_halfhour),
    as.character(gamm_rho),
    as.character(max_position_entropy),
    as.character(low_entropy_cutoff),
    as.character(moderate_entropy_cutoff),
    as.character(filter_cc4_late_phases),
    as.character(cc4_max_active_phase),
    as.character(cc4_max_inactive_phase),
    as.character(nrow(consec_animal_entropy)),
    as.character(nrow(consec_cage_entropy)),
    ifelse(analyze_by_halfhour, as.character(nrow(consec_animal_entropy_halfhour)), NA_character_),
    ifelse(analyze_by_halfhour, as.character(nrow(consec_cage_entropy_halfhour)), NA_character_),
    working_directory,
    saving_directory
  )
)

if (save_tables) {
  readr::write_csv(run_manifest, file.path(pub_tables_dir, "Run_manifest.csv"))
  
  # Also save to analysis_ready structure
  write_analysis_csv(
    run_manifest,
    analysis_logs_dir,
    prefix = "manifest",
    descriptor = "run_settings_and_counts"
  )
}

# ===================================================================
# Generate exploratory plots
# These plots are retained for data inspection. Publication figures are
# generated separately using model-based estimates.
# ===================================================================
if (save_plots == TRUE || show_plots == TRUE) {
  message("=======================================================================")
  message("GENERATING ENTROPY PLOTS")
  message("=======================================================================")
  
# Phase-based animal entropy plots grouped by various factors
columns_to_group <- list(
    c("Sex", "AnimalID", "Group"),
    c("Sex", "AnimalID", "CageChange", "Group"),
    c("Sex", "AnimalID", "Phase", "Group")
)
x_axis <- list("Group", "CageChange", "Group")
phasecount <- 1
plots <- list()

  
  for (i in seq_along(columns_to_group)) {
    group_vars <- columns_to_group[[i]]
    x_var <- x_axis[[i]]
    data_to_plot <- consec_animal_entropy
    
    if ("Phase" %in% group_vars) {
      data_to_plot <- filter(data_to_plot, Phase == ifelse(phasecount == 1, "active", "inactive"))
      phasecount <- phasecount + 1
    }
    
    result <- data_to_plot %>%
      group_by(across(all_of(group_vars))) %>%
      summarise(Mean_Entropy = mean(animalEntropy, na.rm = TRUE), .groups = 'drop')
    
    p <- ggplot(data = result, aes(x = !!sym(x_var), y = Mean_Entropy, color = Group, group = Group)) +
      geom_jitter(aes(fill = Group), size = 4, alpha = 0.7, shape = 16,
                  position = position_dodge(width = ifelse("CageChange" %in% group_vars, 0.75, 0))) +
      labs(title = "Animal Entropy", subtitle = paste("Grouped by:", paste(group_vars, collapse = ", ")),
           x = x_var, y = "Mean Entropy") +
      scale_color_manual(values = c("sus" = "#E63946", "res" = "grey60", "con" = "#457B9D")) +
      facet_grid(ifelse("Phase" %in% group_vars, "Phase~Sex", "~Sex")) +
      stat_summary(fun.min = function(z) {quantile(z, 0.25)},
                   fun.max = function(z) {quantile(z, 0.75)},
                   fun = median, color = "black", size = 0.8, shape = 16,
                   position = position_dodge(width = ifelse("CageChange" %in% group_vars, 0.75, 0))) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 20),
            plot.subtitle = element_text(hjust = 0.5, size = 16),
            legend.key.size = unit(3, "lines"),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 20),
            axis.text = element_text(size = 20),
            axis.title = element_text(size = 22),
            strip.text = element_text(size = 20))
    
    plots[[i]] <- p
    
    if (save_plots) {
      ggsave(filename = file.path(combined_plots_dir, paste0("animal_entropy_phase_plot_", i, ".svg")),
             plot = p, width = 12, height = 8, dpi = 300)
    }
  }
  
  message("   ✓ Saved ", length(plots), " phase-based animal entropy plots")
  
  # Phase-based cage entropy plots active/inactive
  cagePosEntropy_act <- consec_cage_entropy %>% filter(Phase == "active")
  cage_ent_plot_act <- ggplot(data = cagePosEntropy_act,
                              aes(x = ConsecActive, y = CageEntropy, color = System)) +
    geom_jitter(aes(fill = System), size = 4, alpha = 0.7, width = 0.2, shape = 16) +
    scale_y_continuous("Cage Entropy") +
    scale_x_continuous("Active Phase Number") +
    labs(title = "Cage Entropy - Active Phases") +
    stat_summary(fun.min = function(z) {quantile(z, 0.25)},
                 fun.max = function(z) {quantile(z, 0.75)},
                 fun = median, color = "black", size = 0.8, shape = 16) +
    facet_grid(Sex ~ CageChange) +
    theme_minimal() +
    theme(title = element_text(size = 20),
          legend.key.size = unit(3, "lines"),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 22),
          strip.text = element_text(size = 20))
  
  cagePosEntropy_inact <- consec_cage_entropy %>% filter(Phase == "inactive")
  cage_ent_plot_inact <- ggplot(data = cagePosEntropy_inact,
                                aes(x = ConsecInactive, y = CageEntropy, color = System)) +
    geom_jitter(aes(fill = System), size = 4, alpha = 0.7, width = 0.2, shape = 16) +
    scale_y_continuous("Cage Entropy") +
    scale_x_continuous("Inactive Phase Number") +
    labs(title = "Cage Entropy - Inactive Phases") +
    stat_summary(fun.min = function(z) {quantile(z, 0.25)},
                 fun.max = function(z) {quantile(z, 0.75)},
                 fun = median, color = "black", size = 0.8, shape = 16) +
    facet_grid(Sex ~ CageChange) +
    theme_minimal() +
    theme(title = element_text(size = 20),
          legend.key.size = unit(3, "lines"),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 22),
          strip.text = element_text(size = 20))
  
  if (save_plots) {
    ggsave(filename = file.path(combined_plots_dir, "cage_entropy_active.svg"),
           plot = cage_ent_plot_act, width = 14, height = 10, dpi = 300)
    ggsave(filename = file.path(combined_plots_dir, "cage_entropy_inactive.svg"),
           plot = cage_ent_plot_inact, width = 14, height = 10, dpi = 300)
    message("   ✓ Saved cage entropy plots (phase-based)")
  }
  
  # Half-hour plots with consecutive half-hours and per cagechange
  if (analyze_by_halfhour) {
    message("Generating half-hour entropy plots with consecutive numbering (boundaries after each 4th active phase)...")

    # Prefer using phase counts per half-hour to place boundaries after every 4th active phase
    ph_tbl <- consec_animal_entropy_halfhour

    cage_boundaries <- c()
    boundary_labels <- data.frame(x = numeric(0), label = character(0), stringsAsFactors = FALSE)

    if (!is.null(ph_tbl) && nrow(ph_tbl) > 0 && "ConsecHalfHour" %in% colnames(ph_tbl) && "ConsecActive" %in% colnames(ph_tbl)) {
      for (cc in unique(ph_tbl$CageChange)) {
        cc_tbl <- ph_tbl %>% filter(CageChange == cc & Phase == "active" & !is.na(ConsecActive))
        if (nrow(cc_tbl) == 0) next
        max_act <- max(cc_tbl$ConsecActive, na.rm = TRUE)
        if (!is.finite(max_act) || max_act < 4) next
        targets <- seq(4, max_act, by = 4)
        for (t in targets) {
          hh_val <- cc_tbl %>% filter(ConsecActive == t) %>% pull(ConsecHalfHour)
          if (length(hh_val) == 0) next
          hh_max <- suppressWarnings(max(hh_val, na.rm = TRUE))
          if (is.finite(hh_max)) {
            cage_boundaries <- c(cage_boundaries, hh_max)
            boundary_labels <- bind_rows(boundary_labels,
                                         data.frame(x = hh_max,
                                                    label = paste0(cc, " after A", t),
                                                    stringsAsFactors = FALSE))
          }
        }
      }
      if (length(cage_boundaries) > 0) {
        cage_boundaries <- sort(unique(cage_boundaries))
      } else {
        cage_boundaries <- numeric(0)
      }
    } else {
      # Fallback: previous cumulative-max approach if phase counts unavailable
      message("   Warning: phase counts not available for refined boundaries; using fallback cumulative approach.")
      cage_boundaries <- consec_animal_entropy_halfhour %>%
        group_by(CageChange) %>%
        summarise(max_consec = max(ConsecHalfHour, na.rm = TRUE)) %>%
        pull(max_consec) %>%
        cumsum()
      if (length(cage_boundaries) > 0) cage_boundaries <- cage_boundaries[-length(cage_boundaries)]
      if (length(cage_boundaries) > 0) {
        boundary_labels <- data.frame(
          x = cage_boundaries,
          label = paste("CC", 1:length(cage_boundaries), "→", 2:(length(cage_boundaries) + 1)),
          stringsAsFactors = FALSE
        )
      }
    }
    
    animal_halfhour_plot_consec <- ggplot(data = consec_animal_entropy_halfhour,
                                          aes(x = ConsecHalfHour, y = animalEntropy, color = Group)) +
      geom_point(size = 2, alpha = 0.3, na.rm = TRUE) +
      geom_smooth(aes(group = Group), method = "loess", se = TRUE, span = 0.3, na.rm = TRUE) +
      labs(title = "Animal Entropy Over Continuous Time",
           subtitle = "Consecutive half-hour periods across all cage changes",
           x = "Consecutive Half-Hour Period",
           y = "Entropy") +
      scale_color_manual(values = c("sus" = "#E63946", "res" = "grey60", "con" = "#457B9D")) +
      facet_grid(Sex ~ .) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 20),
            plot.subtitle = element_text(hjust = 0.5, size = 16),
            legend.key.size = unit(3, "lines"),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 20),
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 18),
            strip.text = element_text(size = 16))
    
    if (length(cage_boundaries) > 0) {
      animal_halfhour_plot_consec <- animal_halfhour_plot_consec +
        geom_vline(xintercept = cage_boundaries, linetype = "dashed", color = "black", alpha = 0.5) +
        geom_text(data = boundary_labels, aes(x = x, y = Inf, label = label),
                  inherit.aes = FALSE, vjust = 1.5, size = 3, angle = 90, color = "black")
    }
    
    cage_halfhour_plot_consec <- ggplot(data = consec_cage_entropy_halfhour,
                                        aes(x = ConsecHalfHour, y = CageEntropy, color = System)) +
      geom_point(size = 2, alpha = 0.3, na.rm = TRUE) +
      geom_smooth(aes(group = System), method = "loess", se = TRUE, span = 0.3, na.rm = TRUE) +
      labs(title = "Cage Entropy Over Continuous Time",
           subtitle = "Consecutive half-hour periods across all cage changes",
           x = "Consecutive Half-Hour Period",
           y = "Cage Entropy") +
      facet_grid(Sex ~ .) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 20),
            plot.subtitle = element_text(hjust = 0.5, size = 16),
            legend.key.size = unit(3, "lines"),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 20),
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 18),
            strip.text = element_text(size = 16))
    
    if (length(cage_boundaries) > 0) {
      cage_halfhour_plot_consec <- cage_halfhour_plot_consec +
        geom_vline(xintercept = cage_boundaries, linetype = "dashed", color = "black", alpha = 0.5) +
        geom_text(data = boundary_labels, aes(x = x, y = Inf, label = label),
                  inherit.aes = FALSE, vjust = 1.5, size = 3, angle = 90, color = "black")
    }
    
    animal_halfhour_plot_bycage <- ggplot(data = consec_animal_entropy_halfhour,
                                         aes(x = HalfHour, y = animalEntropy, color = Group)) +
      geom_point(size = 2, alpha = 0.5, na.rm = TRUE) +
      geom_smooth(aes(group = Group), method = "loess", se = TRUE, span = 0.5, na.rm = TRUE) +
      labs(title = "Animal Entropy Over Time (Half-Hour Intervals)",
           subtitle = "Separated by cage change",
           x = "Half-Hour Period (within cage change)",
           y = "Entropy") +
      scale_color_manual(values = c("sus" = "#E63946", "res" = "grey60", "con" = "#457B9D")) +
      facet_grid(Sex ~ CageChange) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 20),
            plot.subtitle = element_text(hjust = 0.5, size = 16),
            legend.key.size = unit(3, "lines"),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 20),
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 18),
            strip.text = element_text(size = 16))
    
    cage_halfhour_plot_bycage <- ggplot(data = consec_cage_entropy_halfhour,
                                       aes(x = HalfHour, y = CageEntropy, color = System)) +
      geom_point(size = 2, alpha = 0.5, na.rm = TRUE) +
      geom_smooth(aes(group = System), method = "loess", se = TRUE, span = 0.5, na.rm = TRUE) +
      labs(title = "Cage Entropy Over Time (Half-Hour Intervals)",
           subtitle = "Separated by cage change",
           x = "Half-Hour Period (within cage change)",
           y = "Cage Entropy") +
      facet_grid(Sex ~ CageChange) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 20),
            plot.subtitle = element_text(hjust = 0.5, size = 16),
            legend.key.size = unit(3, "lines"),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 20),
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 18),
            strip.text = element_text(size = 16))
    
    if (save_plots) {
      ggsave(filename = file.path(combined_plots_dir, "animal_entropy_halfhour_consecutive.svg"),
             plot = animal_halfhour_plot_consec, width = 14, height = 10, dpi = 300)
      ggsave(filename = file.path(combined_plots_dir, "cage_entropy_halfhour_consecutive.svg"),
             plot = cage_halfhour_plot_consec, width = 14, height = 10, dpi = 300)
      ggsave(filename = file.path(combined_plots_dir, "animal_entropy_halfhour_by_cagechange.svg"),
             plot = animal_halfhour_plot_bycage, width = 16, height = 10, dpi = 300)
      ggsave(filename = file.path(combined_plots_dir, "cage_entropy_halfhour_by_cagechange.svg"),
             plot = cage_halfhour_plot_bycage, width = 16, height = 10, dpi = 300)
      message("   ✓ Saved 4 half-hour entropy plots")
    }
  }
  
  if (show_plots) {
    for (p in plots) { print(p) }
    print(cage_ent_plot_act)
    print(cage_ent_plot_inact)
    if (analyze_by_halfhour) {
      print(animal_halfhour_plot_consec)
      print(cage_halfhour_plot_consec)
      print(animal_halfhour_plot_bycage)
      print(cage_halfhour_plot_bycage)
    }
  }
}

# ===================================================================
# Additional exploratory and multivariate visualizations
# ===================================================================
message("Generating additional exploration and multivariate plots...")

# 1. Exploration Frequency Plot (% of time animals left home cage)
exploration_freq_plot <- exploration_summary_phase %>%
  ggplot(aes(x = CageChange, y = pct_exploration, color = Group, group = Group)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  facet_grid(Phase ~ Sex) +
  scale_color_manual(values = c("sus" = "#E63946",
                                "res" = "grey60",
                                "con" = "#457B9D")) +
  labs(title = "Exploration Frequency by Group",
       y = "% of time outside home cage",
       x = "Cage Change") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 14))

if (save_plots) {
  ggsave(filename = file.path(exploration_metrics_dir, "exploration_frequency_by_group.svg"),
         exploration_freq_plot, width = 12, height = 8, dpi = 300)
}

# 2. Exploration Diversity Plot (entropy only when animals explored)
explored_only <- consec_animal_entropy %>% filter(explored == TRUE)

exploration_diversity_plot <- ggplot(explored_only, aes(x = Group, y = animalEntropy, fill = Group)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.alpha = 0.4) +
  facet_grid(Phase ~ Sex) +
  stat_compare_means(comparisons = list(c("sus", "con"), c("sus", "res"), c("con", "res")),
                     label = "p.signif", method = "t.test") +
  scale_fill_manual(values = c("sus" = "#E63946", "res" = "grey60", "con" = "#457B9D")) +
  labs(title = "Exploration Diversity (During Exploratory Periods)",
       y = "Shannon Entropy",
       x = "Group") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18))

if (save_plots) {
  ggsave(filename = file.path(exploration_metrics_dir, "exploration_diversity_when_exploring.svg"),
         exploration_diversity_plot, width = 12, height = 8, dpi = 300)
}

# 3. Exploration Category Distribution (stacked bar of categories per group)
exploration_cat_dist <- consec_animal_entropy %>%
  count(Group, Sex, Phase, exploration_category) %>%
  group_by(Group, Sex, Phase) %>%
  mutate(pct = 100 * n / sum(n))

exploration_category_plot <- ggplot(exploration_cat_dist, aes(x = Group, y = pct, fill = exploration_category)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(Phase ~ Sex) +
  scale_fill_manual(values = c(
    "No exploration" = "#d62828",
    "Minimal (1 position)" = "#f77f00",
    "Low diversity" = "#fcbf49",
    "Moderate diversity" = "#06a77d",
    "High diversity" = "#003049"
  ), name = "Exploration Type") +
  labs(title = "Distribution of Exploration Categories",
       y = "Percentage (%)",
       x = "Group") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        legend.key.size = unit(2, "lines"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

if (save_plots) {
  ggsave(filename = file.path(exploration_metrics_dir, "exploration_category_distribution.svg"),
         exploration_category_plot, width = 12, height = 8, dpi = 300)
}

# ===================================================================
# PUBLICATION-READY FIGURES
# ===================================================================
message("Generating publication-ready figures...")

animal_entropy_panel <- animal_model_data %>%
  ggplot2::ggplot(ggplot2::aes(x = Group, y = animalEntropy, colour = Group)) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0.12, height = 0),
    size = 0.65,
    alpha = 0.28
  ) +
  ggplot2::geom_pointrange(
    data = animal_entropy_emm_df,
    ggplot2::aes(x = Group, y = emmean, ymin = lower.CL, ymax = upper.CL),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.28,
    size = 0.35
  ) +
  ggplot2::facet_grid(Phase ~ Sex) +
  ggplot2::scale_colour_manual(values = group_cols, breaks = c("con", "res", "sus")) +
  ggplot2::labs(x = NULL, y = "Individual entropy", title = "Individual spatial entropy") +
  theme_publication() +
  ggplot2::theme(legend.position = "none")

cage_entropy_panel <- consec_cage_entropy %>%
  dplyr::filter(!is.na(CageEntropy), !is.na(Sex), !is.na(Phase)) %>%
  ggplot2::ggplot(ggplot2::aes(x = CageChange, y = CageEntropy, group = System, colour = System)) +
  ggplot2::geom_line(linewidth = 0.25, alpha = 0.45) +
  ggplot2::geom_point(size = 0.8, alpha = 0.55) +
  ggplot2::stat_summary(ggplot2::aes(group = 1), fun = mean, geom = "line",
                        colour = "black", linewidth = 0.55) +
  ggplot2::facet_grid(Phase ~ Sex) +
  ggplot2::labs(x = NULL, y = "Cage entropy", title = "Group spatial dispersion") +
  theme_publication() +
  ggplot2::theme(legend.position = "none")

exploration_frequency_panel <- pub_exploration_frequency %>%
  ggplot2::ggplot(ggplot2::aes(x = CageChange, y = pct_explored,
                               colour = Group, group = Group)) +
  ggplot2::geom_line(linewidth = 0.45) +
  ggplot2::geom_point(size = 1.1) +
  ggplot2::facet_grid(Phase ~ Sex) +
  ggplot2::scale_colour_manual(values = group_cols, breaks = c("con", "res", "sus")) +
  ggplot2::labs(x = NULL, y = "Exploration (%)", title = "Exploration frequency") +
  theme_publication() +
  ggplot2::theme(legend.position = "top")

if (analyze_by_halfhour) {
  animal_entropy_time_panel <- consec_animal_entropy_halfhour %>%
    dplyr::filter(!is.na(animalEntropy), !is.na(Group), !is.na(ConsecHalfHour)) %>%
    dplyr::group_by(Group, Sex, ConsecHalfHour) %>%
    dplyr::summarise(
      mean_entropy = mean(animalEntropy, na.rm = TRUE),
      sem_entropy = sem(animalEntropy),
      .groups = "drop"
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = ConsecHalfHour, y = mean_entropy,
                                 colour = Group, fill = Group)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = mean_entropy - sem_entropy,
                                      ymax = mean_entropy + sem_entropy),
                         alpha = 0.16, colour = NA) +
    ggplot2::geom_line(linewidth = 0.45) +
    ggplot2::facet_grid(Sex ~ .) +
    ggplot2::scale_colour_manual(values = group_cols, breaks = c("con", "res", "sus")) +
    ggplot2::scale_fill_manual(values = group_cols, breaks = c("con", "res", "sus")) +
    ggplot2::labs(x = "Consecutive half-hour", y = "Individual entropy",
                  title = "Individual entropy over time") +
    theme_publication() +
    ggplot2::theme(legend.position = "top")
}

if (save_plots) {
  save_pub_plot(animal_entropy_panel, "Fig1A_individual_entropy_phase.svg", width = 90, height = 75)
  save_pub_plot(cage_entropy_panel, "Fig1B_cage_entropy_phase.svg", width = 90, height = 75)
  save_pub_plot(exploration_frequency_panel, "Fig1C_exploration_frequency.svg", width = 90, height = 75)

  if (analyze_by_halfhour) {
    save_pub_plot(animal_entropy_time_panel, "Fig2A_individual_entropy_timecourse.svg", width = 180, height = 80)
  }

  figure_1 <- (animal_entropy_panel | cage_entropy_panel | exploration_frequency_panel) +
    patchwork::plot_annotation(tag_levels = "A")

  save_pub_multipanel(figure_1, "Figure_1_entropy_exploration_overview.svg", width = 180, height = 80)

  if (analyze_by_halfhour) {
    figure_2 <- animal_entropy_time_panel +
      patchwork::plot_annotation(tag_levels = "A")

    save_pub_multipanel(figure_2, "Figure_2_entropy_timecourse.svg", width = 180, height = 80)
  }
  
  # Also save to analysis_ready structure
  message("\nSaving figures to analysis_ready structure...")
  save_analysis_figure(
    animal_entropy_panel,
    analysis_figures_model_dir,
    prefix = "figure",
    descriptor = "fig1a_individual_entropy_phase",
    width = 90, height = 75
  )
  message("   ✓ Saved: figure__fig1a_individual_entropy_phase.svg")
  
  save_analysis_figure(
    cage_entropy_panel,
    analysis_figures_model_dir,
    prefix = "figure",
    descriptor = "fig1b_cage_entropy_phase",
    width = 90, height = 75
  )
  message("   ✓ Saved: figure__fig1b_cage_entropy_phase.svg")
  
  save_analysis_figure(
    exploration_frequency_panel,
    analysis_figures_model_dir,
    prefix = "figure",
    descriptor = "fig1c_exploration_frequency",
    width = 90, height = 75
  )
  message("   ✓ Saved: figure__fig1c_exploration_frequency.svg")
  
  if (analyze_by_halfhour) {
    save_analysis_figure(
      animal_entropy_time_panel,
      analysis_figures_model_dir,
      prefix = "figure",
      descriptor = "fig2a_individual_entropy_timecourse",
      width = 180, height = 80
    )
    message("   ✓ Saved: figure__fig2a_individual_entropy_timecourse.svg")
  }
  
  save_analysis_figure(
    figure_1,
    analysis_figures_multipanel_dir,
    prefix = "figure",
    descriptor = "figure_1_entropy_exploration_overview",
    width = 180, height = 80
  )
  message("   ✓ Saved: figure__figure_1_entropy_exploration_overview.svg")
  
  if (analyze_by_halfhour) {
    save_analysis_figure(
      figure_2,
      analysis_figures_multipanel_dir,
      prefix = "figure",
      descriptor = "figure_2_entropy_timecourse",
      width = 180, height = 80
    )
    message("   ✓ Saved: figure__figure_2_entropy_timecourse.svg")
  }
}

# ===================================================================
# STATISTICAL ANALYSIS: SEX x GROUP (2x3 DESIGN)
# ===================================================================

# Define specific sub-folders for Sex x Group analysis output
sex_stats_dir <- file.path(exploration_metrics_dir, "sex_comparisons")
if (save_tables && !dir.exists(sex_stats_dir)) dir.create(sex_stats_dir, recursive = TRUE)

# We check if a corresponding plots folder exists in the main plots structure, otherwise we create it
# As 'exploration_metrics_dir' is inside 'tables', we should put plots in 'plots/exploration_metrics/sex_comparisons'
sex_plots_dir <- file.path(plots_dir, "exploration_metrics", "sex_comparisons")
if (save_plots && !dir.exists(sex_plots_dir)) dir.create(sex_plots_dir, recursive = TRUE)

# 1. 2x3 ANOVA for Exploration Diversity (Entropy when exploring)
# -------------------------------------------------------------------
message("Running 2x3 ANOVA (Sex x Group) for Exploration Diversity...")

anova_results <- list()
posthoc_results_list <- list()

for (phase_val in c("active", "inactive")) {
  # Filter for explore-only events to analyze diversity quality
  model_data <- consec_animal_entropy %>%
    filter(Phase == phase_val, explored == TRUE)
  
  if (nrow(model_data) > 10) { # Ensure enough data points
    # Two-way ANOVA with Interaction
    aov_model <- aov(animalEntropy ~ Sex * Group, data = model_data)
    
    # Store summary
    tidy_aov <- broom::tidy(aov_model)
    tidy_aov$Phase <- phase_val
    anova_results[[phase_val]] <- tidy_aov
    
    message(paste0("   Phase: ", phase_val))
    print(summary(aov_model))
    
    # ADEQUATE POST-HOC: Estimated Marginal Means with Tukey Adjustment
    # Using emmeans for robust pairwise comparisons and Type I error control
    
    # Load emmeans Model
    emm_model <- emmeans(aov_model, ~ Sex * Group)
    
    # A. Pairwise comparisons: Groups within each Sex
    pairs_group_by_sex <- pairs(emm_model, by = "Sex", adjust = "tukey")
    message("   --- Pairwise Comparisons: Group within Sex (Tukey) ---")
    print(pairs_group_by_sex)
    
    # B. Pairwise comparisons: Sex within each Group
    pairs_sex_by_group <- pairs(emm_model, by = "Group", adjust = "tukey")
    message("   --- Pairwise Comparisons: Sex within Group (Tukey) ---")
    print(pairs_sex_by_group)
    
    # Save posthoc stats
    res_g_s <- as.data.frame(pairs_group_by_sex)
    res_g_s$ComparisonType <- "Group_within_Sex"
    res_g_s$Phase <- phase_val
    
    res_s_g <- as.data.frame(pairs_sex_by_group)
    res_s_g$ComparisonType <- "Sex_within_Group"
    res_s_g$Phase <- phase_val
    
    posthoc_results_list[[paste0(phase_val, "_Gs")]] <- res_g_s
    posthoc_results_list[[paste0(phase_val, "_Sg")]] <- res_s_g
  }
}

if (length(anova_results) > 0) {
  all_anova_results <- bind_rows(anova_results)
  if (save_tables) {
    write.csv(all_anova_results, file = file.path(sex_stats_dir, "SexGroup_Diversity_ANOVA.csv"), row.names = FALSE)
  }
}

if (length(posthoc_results_list) > 0) {
  all_posthoc_results <- bind_rows(posthoc_results_list)
  if (save_tables) {
    write.csv(all_posthoc_results, file = file.path(sex_stats_dir, "SexGroup_Diversity_PostHoc.csv"), row.names = FALSE)
  }
}

# 2. Categorical Analysis: Sex differences WITHIN Groups
# -------------------------------------------------------------------
message("Running Chi-Square analysis for Sex differences within each Group (Sex impacts on Category)...")

sex_group_cat_stats <- list()

for (phase_val in c("active", "inactive")) {
  for (group_val in unique(consec_animal_entropy$Group)) {
    
    # Create contingency table: Sex vs Exploration Category for specific Group & Phase
    test_data <- consec_animal_entropy %>%
      filter(Phase == phase_val, Group == group_val) %>%
      select(Sex, exploration_category) %>%
      table()
    
    # Only run if we have data for both sexes
    if (nrow(test_data) == 2 && ncol(test_data) > 1) {
      
      expected_counts <- tryCatch(chisq.test(test_data)$expected, error = function(e) matrix(0,1,1))
      min_expected <- min(expected_counts)
      
      # Choose test based on sample size
      if (min_expected >= 5) {
        test_res <- chisq.test(test_data)
        method <- "Chi-square"
        pval <- test_res$p.value
        stat <- test_res$statistic
      } else {
        test_res <- fisher.test(test_data, simulate.p.value = TRUE)
        method <- "Fisher's Exact"
        pval <- test_res$p.value
        stat <- NA
      }
      
      sex_group_cat_stats[[paste(phase_val, group_val, sep="_")]] <- tibble(
        Phase = phase_val,
        Group = group_val,
        Test_Method = method,
        Statistic = as.numeric(stat),
        P_Value = as.numeric(pval),
        Significant = pval < 0.05
      )
    }
  }
}

if (length(sex_group_cat_stats) > 0) {
  sex_group_stats_df <- bind_rows(sex_group_cat_stats)
  print(sex_group_stats_df)
  if (save_tables) {
    write.csv(sex_group_stats_df, file = file.path(sex_stats_dir, "SexGroup_Category_ChiSq.csv"), row.names = FALSE)
  }
}

# 2b. Exploration Diversity - Sex Comparison (within Group) PLOT
exploration_diversity_sex_plot <- ggplot(explored_only, aes(x = Group, y = animalEntropy, color = Sex, group = Sex)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, position = position_dodge(width = 0.8), color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(width = 0.8), color = "black") +
  facet_grid(Phase ~ .) +
  stat_compare_means(aes(group = Sex), label = "p.signif", method = "t.test", hide.ns = FALSE) +
  scale_color_manual(values = c("male" = "#8585ab", "female" = "#d7c39d")) +
  labs(title = "Exploration Diversity: Sex Differences within Groups",
       subtitle = "Shannon Entropy during exploratory periods",
       y = "Shannon Entropy",
       x = "Group") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 14))

if (save_plots) {
  ggsave(filename = file.path(sex_plots_dir, "SexGroup_Diversity_Plot.svg"),
         exploration_diversity_sex_plot, width = 8, height = 8, dpi = 300)
}

# ===================================================================
# 3b. Exploration Diversity - Averaged per Animal (Sex Comparison)
# ===================================================================
# Create an averaged summary per animal to avoid pseudo-replication
explored_averaged <- explored_only %>%
  group_by(AnimalID, Group, Sex, Phase) %>%
  summarise(mean_animalEntropy = mean(animalEntropy, na.rm = TRUE), .groups = "drop")

# -------------------------------------------------------------------
# Statistical Analysis: 2x3 ANOVA & Tukey Post-Hoc for Averaged Data
# -------------------------------------------------------------------
message("Running 2x3 ANOVA (Sex x Group) for Averaged Exploration Diversity...")

anova_avg_results <- list()
posthoc_avg_results <- list()

for (phase_val in c("active", "inactive")) {
  # Filter data for specific phase
  model_data <- explored_averaged %>% filter(Phase == phase_val)
  
  # Ensure sufficient data points
  if (nrow(model_data) > 6) {
    # 1. Run ANOVA
    aov_model <- aov(mean_animalEntropy ~ Sex * Group, data = model_data)
    
    # Store ANOVA results
    tidy_res <- broom::tidy(aov_model)
    tidy_res$Phase <- phase_val
    anova_avg_results[[phase_val]] <- tidy_res
    
    # Print summary to console
    message(paste0("   Phase (Avg): ", phase_val))
    print(summary(aov_model))
    
    # 2. Tukey Post-Hoc Tests via emmeans
    emm_model <- emmeans(aov_model, ~ Sex * Group)
    
    # A. Pairwise: Sex within Group
    pairs_sex_within_group <- pairs(emm_model, by = "Group", adjust = "tukey")
    res_s_g <- as.data.frame(pairs_sex_within_group)
    res_s_g$ComparisonType <- "Sex_within_Group"
    res_s_g$Phase <- phase_val
    res_s_g$PostHoc_Test <- "Tukey"
    res_s_g$Adjustment_Method <- "tukey"
    
    # B. Pairwise: Group within Sex
    pairs_group_within_sex <- pairs(emm_model, by = "Sex", adjust = "tukey")
    res_g_s <- as.data.frame(pairs_group_within_sex)
    res_g_s$ComparisonType <- "Group_within_Sex"
    res_g_s$Phase <- phase_val
    res_g_s$PostHoc_Test <- "Tukey"
    res_g_s$Adjustment_Method <- "tukey"
    
    # Rename columns for clarity in CSV if needed (emmeans usually gives p.value)
    # Ensure p-value info is explicit
    # emmeans output usually has columns: contrast, <factors>, estimate, SE, df, t.ratio, p.value
    
    posthoc_avg_results[[paste0(phase_val, "_Sg")]] <- res_s_g
    posthoc_avg_results[[paste0(phase_val, "_Gs")]] <- res_g_s
  }
}

# Save ANOVA and Post-Hoc Tables
if (length(anova_avg_results) > 0 && save_tables) {
  all_anova_avg <- bind_rows(anova_avg_results)
  write.csv(all_anova_avg, 
            file = file.path(sex_stats_dir, "SexGroup_Diversity_Avg_ANOVA.csv"), 
            row.names = FALSE)
  message("   ✓ Saved averaged ANOVA results")
}

if (length(posthoc_avg_results) > 0 && save_tables) {
  all_posthoc_avg <- bind_rows(posthoc_avg_results)
  write.csv(all_posthoc_avg, 
            file = file.path(sex_stats_dir, "SexGroup_Diversity_Avg_PostHoc.csv"), 
            row.names = FALSE)
  message("   ✓ Saved averaged Post-hoc results")
}


# -------------------------------------------------------------------
# Plotting Averaged Data
# -------------------------------------------------------------------

# Prepare annotation dataframe from previously computed Post-Hoc results to ensure plot matches tables
stat_test_df <- tibble()
if (exists("posthoc_avg_results") && length(posthoc_avg_results) > 0) {
  
  # Bind all results
  all_ph <- bind_rows(posthoc_avg_results)
  
  # Filter for Sex differences within Groups (Sg)
  if ("ComparisonType" %in% colnames(all_ph)) {
    stat_test_df <- all_ph %>%
      filter(ComparisonType == "Sex_within_Group")
    
    if (nrow(stat_test_df) > 0) {
      # Prepare for ggpubr: needs group1, group2, y.position
      stat_test_df <- stat_test_df %>%
        mutate(
          Group = as.character(Group),
          group1 = "male",
          group2 = "female",
          p.signif = case_when(
            p.value < 0.001 ~ "***",
            p.value < 0.01  ~ "**",
            p.value < 0.05  ~ "*",
            TRUE            ~ "ns"
          )
        )
      
      # Determine Y-position for brackets (max data value + buffer)
      y_max_values <- explored_averaged %>%
        group_by(Group, Phase) %>%
        summarise(max_val = max(mean_animalEntropy, na.rm = TRUE), .groups = "drop")
      
      stat_test_df <- stat_test_df %>%
        left_join(y_max_values, by = c("Group", "Phase")) %>%
        mutate(y.position = max_val + 0.3) # Shift bracket slightly higher
    }
  }
}

group_cols <- c("con" = "#477c9e", "res" = "#c6c3bb", "sus" = "#e63947")

exploration_diversity_sex_plot_avg <- ggplot(explored_averaged, aes(x = Group, y = mean_animalEntropy)) +
  geom_jitter(aes(color = Group, shape = Sex, group = Sex), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
              size = 4, alpha = 0.7) +
  # Changed from point (shape 18) to crossbar or errorbar with just the y to represent mean as a horizontal line
  stat_summary(aes(group = Sex), fun = mean, geom = "crossbar", width = 0.6, 
               position = position_dodge(width = 0.8), color = "black", show.legend = FALSE, fatten = 1, size = 0.4) +
  stat_summary(aes(group = Sex), fun.data = mean_se, geom = "errorbar", width = 0.2, 
               position = position_dodge(width = 0.8), color = "black", show.legend = FALSE) +
  facet_grid(Phase ~ ., scales = "free_y") +
  scale_color_manual(values = group_cols) +
  scale_shape_manual(values = c("female" = 16, "male" = 1)) +
  labs(title = "Exploration Diversity: Sex Differences within Groups",
       subtitle = "Mean Shannon Entropy per animal (averaged across explored phases)",
       y = "Mean Shannon Entropy",
       x = "Group") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Add statistical brackets: Use pre-computed ANOVA/Tukey if available, else fallback to t-test
if (nrow(stat_test_df) > 0) {
  exploration_diversity_sex_plot_avg <- exploration_diversity_sex_plot_avg +
    stat_pvalue_manual(
      stat_test_df, 
      x = "Group",
      label = "p.signif", 
      tip.length = 0.01,
      hide.ns = FALSE
    )
} else {
  exploration_diversity_sex_plot_avg <- exploration_diversity_sex_plot_avg +
    geom_pwc(
      aes(group = Sex), 
      method = "t_test", 
      label = "p.signif",
      bracket.nudge.y = 0.2, 
      hide.ns = FALSE
    )
}

if (save_plots) {
  ggsave(filename = file.path(sex_plots_dir, "SexGroup_Diversity_Avg_Plot.svg"),
         exploration_diversity_sex_plot_avg, width = 5, height = 5, dpi = 300)
}

# 3b. Exploration Category Distribution - Sex Comparison (within Group) PLOT & STATS
# -------------------------------------------------------------------
# Statistical Analysis: 2x3 ANOVA (Parametric) or Rank-Based (Non-Parametric)
# -------------------------------------------------------------------
# We treat 'pct' of exploration categories as the dependent variable.
# We analyze the proportion of "meaningful exploration" (100% - "No exploration").

message("Running Sex x Group Analysis for Proportion of Meaningful Exploration...")
message("   Step 1: Checking Normal Distribution of Residuals (Shapiro-Wilk)")
message("   Step 2: Selecting Parametric (ANOVA) or Non-Parametric (Wilcoxon/Kruskal) tests")

# First, calculate the per-animal percentage of time spent in meaningful exploration categories
# to avoid pseudo-replication on the summary table.
meaningful_exploration_data <- consec_animal_entropy %>%
  mutate(is_meaningful = exploration_category != "No exploration") %>%
  group_by(AnimalID, Group, Sex, Phase) %>%
  summarise(
    total_obs = n(),
    meaningful_obs = sum(is_meaningful),
    pct_meaningful = 100 * meaningful_obs / total_obs,
    .groups = "drop"
  )

anova_cat_results <- list()
posthoc_cat_results <- list()
normality_log <- list()

for (phase_val in c("active", "inactive")) {
  # Filter data for specific phase
  model_data <- meaningful_exploration_data %>% filter(Phase == phase_val)
  
  if (nrow(model_data) > 6 && n_distinct(model_data$Sex) > 1 && n_distinct(model_data$Group) > 1) {
    
    # Run initial ANOVA model to check residuals
    aov_model <- aov(pct_meaningful ~ Sex * Group, data = model_data)
    
    # Check Normality
    # Use residuals from the model
    resid_shapiro <- tryCatch(shapiro.test(residuals(aov_model)), error = function(e) list(p.value = 1))
    
    is_normal <- resid_shapiro$p.value > 0.05
    test_method <- ifelse(is_normal, "Parametric (ANOVA + Tukey)", "Non-Parametric (Split Wilcoxon/Kruskal)")
    
    normality_log[[phase_val]] <- tibble(
      Phase = phase_val,
      Shapiro_P = resid_shapiro$p.value,
      Is_Normal = is_normal,
      Selected_Method = test_method
    )
    
    message(sprintf("   Phase: %s | Normality p=%.4f | Method: %s", phase_val, resid_shapiro$p.value, test_method))
    
    if (is_normal) {
      # --- PARAMETRIC ---
      
      # 1. Main ANOVA
      tidy_res <- broom::tidy(aov_model)
      tidy_res$Phase <- phase_val
      tidy_res$Test_Type <- "ANOVA"
      anova_cat_results[[phase_val]] <- tidy_res
      
      print(summary(aov_model))
      
      # 2. Post-Hoc (Tukey)
      emm_model <- emmeans(aov_model, ~ Sex * Group)
      
      # Sex within Group
      pairs_sex_within_group <- pairs(emm_model, by = "Group", adjust = "tukey")
      res_s_g <- as.data.frame(pairs_sex_within_group)
      res_s_g$ComparisonType <- "Sex_within_Group"
      res_s_g$Phase <- phase_val
      res_s_g$Test_Type <- "Tukey"
      
      # Group within Sex
      pairs_group_within_sex <- pairs(emm_model, by = "Sex", adjust = "tukey")
      res_g_s <- as.data.frame(pairs_group_within_sex)
      res_g_s$ComparisonType <- "Group_within_Sex"
      res_g_s$Phase <- phase_val
      res_g_s$Test_Type <- "Tukey"
      
      posthoc_cat_results[[paste0(phase_val, "_Sg")]] <- res_s_g
      posthoc_cat_results[[paste0(phase_val, "_Gs")]] <- res_g_s
      
    } else {
      # --- NON-PARAMETRIC ---
      # Assumptions violated. Perform stratified non-parametric tests.
      
      # 1. Sex within Group (Mann-Whitney U / Wilcoxon Rank Sum)
      # Equivalent to the plot comparison
      s_g_list <- list()
      for(g in unique(model_data$Group)) {
        sub_d <- model_data %>% filter(Group == g)
        if(n_distinct(sub_d$Sex) == 2) {
          wt <- wilcox.test(pct_meaningful ~ Sex, data = sub_d, exact = FALSE)
          s_g_list[[g]] <- tibble(
            Group = g,
            contrast = "male - female", # Direction label
            estimate = NA,
            p.value = wt$p.value,
            ComparisonType = "Sex_within_Group",
            Phase = phase_val,
            Test_Type = "Wilcoxon"
          )
        }
      }
      if(length(s_g_list) > 0) posthoc_cat_results[[paste0(phase_val, "_Sg")]] <- bind_rows(s_g_list)
      
      # 2. Group within Sex (Kruskal-Wallis + Pairwise Wilcoxon with BH adjust)
      g_s_list <- list()
      for(s in unique(model_data$Sex)) {
        sub_d <- model_data %>% filter(Sex == s)
        if(n_distinct(sub_d$Group) > 1) {
          # Pairwise Wilcoxon
          pwt <- pairwise.wilcox.test(sub_d$pct_meaningful, sub_d$Group, p.adjust.method = "BH")
          
          # Convert triangular p-value matrix to table
          if (!is.null(pwt$p.value)) {
            p_mat <- as.table(pwt$p.value)
            p_df <- as.data.frame(p_mat) %>% 
              na.omit() %>%
              rename(group1 = Var1, group2 = Var2, p.value = Freq) %>%
              mutate(
                Sex = s,
                contrast = paste(group1, "-", group2),
                estimate = NA,
                ComparisonType = "Group_within_Sex",
                Phase = phase_val,
                Test_Type = "WilcoxBH"
              )
            g_s_list[[s]] <- p_df
          }
        }
      }
      if(length(g_s_list) > 0) posthoc_cat_results[[paste0(phase_val, "_Gs")]] <- bind_rows(g_s_list)
      
      # We do not produce a global ANOVA table for non-parametric split path
      # but we can log the decision
      anova_cat_results[[phase_val]] <- tibble(
        term = c("Normality_Check", "Method"),
        statistic = c(resid_shapiro$statistic, NA),
        p.value = c(resid_shapiro$p.value, NA),
        Phase = phase_val,
        Test_Type = "Assumption Check"
      )
    }
  }
}

# Save Stats Tables
if (save_tables) {
  if (length(normality_log) > 0) {
    write.csv(bind_rows(normality_log),
              file = file.path(sex_stats_dir, "SexGroup_MeaningfulExpl_Normality.csv"),
              row.names = FALSE)
  }
  if (length(anova_cat_results) > 0) {
    write.csv(bind_rows(anova_cat_results), 
              file = file.path(sex_stats_dir, "SexGroup_MeaningfulExpl_Main.csv"), 
              row.names = FALSE)
  }
  if (length(posthoc_cat_results) > 0) {
    write.csv(bind_rows(posthoc_cat_results), 
              file = file.path(sex_stats_dir, "SexGroup_MeaningfulExpl_PostHoc.csv"), 
              row.names = FALSE)
  }
}

# Prepare annotation for plot (Sex differences within Group)
stat_cat_df <- tibble()
if (length(posthoc_cat_results) > 0) {
  all_cat_ph <- bind_rows(posthoc_cat_results)
  
  if ("ComparisonType" %in% colnames(all_cat_ph)) {
    stat_cat_df <- all_cat_ph %>%
      filter(ComparisonType == "Sex_within_Group") %>%
      mutate(
        Group = as.character(Group),
        group1 = "male", # mapping for ggpubr/plots usually needs x-axis checks
        group2 = "female",
        p.signif = case_when(
          p.value < 0.001 ~ "***",
          p.value < 0.01  ~ "**",
          p.value < 0.05  ~ "*",
          TRUE            ~ "ns"
        )
      ) %>%
      filter(p.signif != "ns") # Keep only significant
  }
}

# Plot Generation
exploration_category_sex_plot <- ggplot(exploration_cat_dist, aes(x = Sex, y = pct, fill = exploration_category)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(Phase ~ Group) +
  scale_fill_manual(values = c(
    "No exploration" = "#d62828",
    "Minimal (1 position)" = "#f77f00",
    "Low diversity" = "#fcbf49",
    "Moderate diversity" = "#06a77d",
    "High diversity" = "#003049"
  ), name = "Exploration Type") +
  labs(title = "Exploration Categories: Sex Differences within Groups",
       subtitle = "Significant differences in meaningful exploration (Test depends on Normality)",
       y = "Percentage (%)",
       x = "Sex") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_text(size = 12),
        strip.text = element_text(size = 14),
        legend.key.size = unit(2, "lines"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

# Add statistical annotations
if (nrow(stat_cat_df) > 0) {
  # Add y.position for brackets (above 100%)
  stat_cat_df$y.position <- 105
  stat_cat_df$xmin <- 1 # male on factor axis (usually alphabetic order by default in R factor, check levels if custom)
  stat_cat_df$xmax <- 2 # female
  
  # Ensure Sex order is alphabetic or matches plot
  # If ggplot default, female=1, male=2 (f comes before m). 
  # Let's verify levels of Sex in data frame data_preprocessed?
  # "female", "male". -> female=1, male=2.
  stat_cat_df$xmin <- 1 # female
  stat_cat_df$xmax <- 2 # male
  
  exploration_category_sex_plot <- exploration_category_sex_plot +
    geom_signif(
      data = stat_cat_df,
      aes(xmin = xmin, xmax = xmax, annotations = p.signif, y_position = y.position, group = NULL),
      manual = TRUE,
      inherit.aes = FALSE,
      tip_length = 0.01,
      vjust = 0.5
    ) +
    scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, 25))
}

if (save_plots) {
  ggsave(filename = file.path(sex_plots_dir, "SexGroup_Category_Dist_Plot.svg"),
         exploration_category_sex_plot, width = 8, height = 8, dpi = 300)
}

# ===================================================================
# STATISTICAL TESTS FOR EXPLORATION CATEGORIES
# ===================================================================
message("Running statistical tests on exploration categories...")

cat_stats_results <- list()

for (sex_val in c("male", "female")) {
  for (phase_val in c("active", "inactive")) {
    test_data <- consec_animal_entropy %>%
      filter(Sex == sex_val, Phase == phase_val) %>%
      select(Group, exploration_category) %>%
      table()
    
    if (sum(test_data) > 0 && nrow(test_data) > 1 && ncol(test_data) > 1) {
      expected_counts <- chisq.test(test_data)$expected
      min_expected <- min(expected_counts)
      
      if (min_expected >= 5) {
        chi_test <- chisq.test(test_data)
        cat_stats_results[[paste(sex_val, phase_val, sep = "_")]] <- list(
          sex = sex_val,
          phase = phase_val,
          chi_square = as.numeric(chi_test$statistic),
          p_value = as.numeric(chi_test$p.value),
          df = as.numeric(chi_test$parameter),
          test_type = "Chi-square test",
          min_expected_count = min_expected,
          significant = as.numeric(chi_test$p.value) < 0.05
        )
        
        message(sprintf("   %s %s: χ² = %.3f, df = %d, p = %.4f (Significant: %s)",
                        sex_val, phase_val, chi_test$statistic, chi_test$parameter,
                        chi_test$p.value, ifelse(chi_test$p.value < 0.05, "Yes", "No")))
      } else {
        fisher_test <- fisher.test(test_data, simulate.p.value = TRUE)
        cat_stats_results[[paste(sex_val, phase_val, sep = "_")]] <- list(
          sex = sex_val,
          phase = phase_val,
          chi_square = NA,
          p_value = as.numeric(fisher_test$p.value),
          df = NA,
          test_type = "Fisher's exact test",
          min_expected_count = min_expected,
          significant = as.numeric(fisher_test$p.value) < 0.05
        )
        
        message(sprintf("   %s %s: Fisher's exact test p = %.4f (min expected = %.2f, Significant: %s)",
                        sex_val, phase_val, fisher_test$p.value, min_expected,
                        ifelse(fisher_test$p.value < 0.05, "Yes", "No")))
      }
    } else {
      message(sprintf("   %s %s: Insufficient data for test", sex_val, phase_val))
    }
  }
}

# Adjust p-values for multiple testing using Benjamini-Hochberg procedure
p_values <- sapply(cat_stats_results, function(x) x$p_value)
adjusted_p_values <- p.adjust(p_values, method = "BH")
for (i in seq_along(cat_stats_results)) {
  cat_stats_results[[i]]$adjusted_p_value <- adjusted_p_values[i]
}

# Pairwise comparisons for proportions of each category between groups
pairwise_results <- list()

for (sex_val in c("male", "female")) {
  for (phase_val in c("active", "inactive")) {
    for (cat in levels(consec_animal_entropy$exploration_category)) {
      test_data <- consec_animal_entropy %>%
        filter(Sex == sex_val, Phase == phase_val) %>%
        mutate(is_category = exploration_category == cat)
      
      if (nrow(test_data) > 0) {
        group_props <- test_data %>%
          group_by(Group) %>%
          summarise(n_total = n(), n_category = sum(is_category), prop = mean(is_category), .groups = 'drop')
        
        groups <- unique(test_data$Group)
        comparisons <- combn(groups, 2, simplify = FALSE)
        
        for (comp in comparisons) {
          g1_data <- filter(group_props, Group == comp[1])
          g2_data <- filter(group_props, Group == comp[2])
          if (nrow(g1_data) > 0 && nrow(g2_data) > 0) {
            tryCatch({
              prop_test <- prop.test(c(g1_data$n_category, g2_data$n_category), c(g1_data$n_total, g2_data$n_total))
              key <- paste(sex_val, phase_val, cat, paste(comp, collapse = "_vs_"), sep = "_")
              pairwise_results[[key]] <- list(
                sex = sex_val,
                phase = phase_val,
                category = cat,
                comparison = paste(comp, collapse = " vs "),
                prop1 = g1_data$prop,
                prop2 = g2_data$prop,
                chi_square = as.numeric(prop_test$statistic),
                p_value = as.numeric(prop_test$p.value),
                significant = as.numeric(prop_test$p.value) < 0.05
              )
            }, error = function(e) {
              message(sprintf("   Warning: Could not compute proportion test for %s %s %s: %s",
                              sex_val, phase_val, cat, e$message))
            })
          }
        }
      }
    }
  }
}

# Convert results to data.frames and save
if (length(cat_stats_results) > 0) {
  cat_stats_df <- bind_rows(cat_stats_results)
} else {
  cat_stats_df <- tibble(sex = character(), phase = character(),
                        chi_square = numeric(), p_value = numeric(),
                        adjusted_p_value = numeric(), df = numeric(),
                        test_type = character(), min_expected_count = numeric(),
                        significant = logical())
}
if (length(pairwise_results) > 0) {
  pairwise_stats_df <- bind_rows(pairwise_results)
} else {
  pairwise_stats_df <- tibble(sex = character(), phase = character(),
                             category = character(), comparison = character(),
                             prop1 = numeric(), prop2 = numeric(),
                             chi_square = numeric(), p_value = numeric(),
                             significant = logical())
}

if (save_tables) {
  write.csv(cat_stats_df, file = file.path(exploration_metrics_dir, "exploration_category_chisquare_tests.csv"), row.names = FALSE)
  write.csv(pairwise_stats_df, file = file.path(exploration_metrics_dir, "exploration_category_pairwise_tests.csv"), row.names = FALSE)
  message("   ✓ Saved exploration category statistical tests")
}

# 4. Individual Animal Exploration Profiles (frequency vs mean entropy when exploring)
individual_profile_plot <- ggplot(animal_exploration_profile, aes(x = pct_exploration, y = mean_entropy_when_exploring, color = Group)) +
  geom_point(size = 3, alpha = 0.6) +
  facet_wrap(~ Sex) +
  scale_color_manual(values = c("sus" = "#E63946", "res" = "grey60", "con" = "#457B9D")) +
  labs(title = "Individual Animal Exploration Profiles",
       subtitle = "Exploration Frequency vs Diversity",
       x = "Exploration Frequency (%)",
       y = "Mean Entropy When Exploring") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        legend.position = "bottom")

if (save_plots) {
  ggsave(filename = file.path(exploration_metrics_dir, "animal_exploration_profiles.svg"),
         individual_profile_plot, width = 12, height = 10, dpi = 300)
}

# ===================================================================
# TRANSITION NETWORK PLOTS PER ANIMAL
# ===================================================================
message("Generating transition network plots for each animal...")

# Provide a fallback implementation if plot_transition_network is not defined
if (!exists("plot_transition_network") || !is.function(plot_transition_network)) {
  plot_transition_network <- function(seq_positions, animal_id = NULL) {
    if (length(seq_positions) < 2) stop("Need at least two positions to build transitions.")
    edges <- data.frame(from = head(seq_positions, -1), to = tail(seq_positions, -1), stringsAsFactors = FALSE)
    edges <- edges %>% group_by(from, to) %>% summarise(weight = n(), .groups = "drop")
    g <- graph_from_data_frame(edges, directed = TRUE)
    E(g)$weight <- edges$weight
    V(g)$label <- V(g)$name
    node_degree <- degree(g, mode = "all")
    p <- tryCatch({
      ggraph(g, layout = "fr") +
        geom_edge_fan(aes(width = weight, alpha = weight), arrow = arrow(length = unit(3, "mm"), type = "closed"), show.legend = TRUE) +
        geom_node_point(aes(size = node_degree), color = "steelblue") +
        geom_node_text(aes(label = label), repel = TRUE, size = 3) +
        scale_edge_width(range = c(0.2, 3)) +
        scale_size(range = c(3, 10)) +
        theme_void() +
        ggtitle(paste("Transition network", ifelse(is.null(animal_id), "", paste0("-", animal_id))))
    }, error = function(e) stop("ggraph plotting failed: ", e$message))
    return(p)
  }
  message("✓ Defined fallback plot_transition_network()")
}

# Load preprocessed data if needed
if (!exists("data_preprocessed") || is.null(data_preprocessed) || nrow(data_preprocessed) == 0) {
  message("Reloading preprocessed data for transition network analysis...")
  all_preprocessed <- tibble()
  for (batch in batches) {
    for (cageChange in cageChanges) {
      filename <- paste0("E9_SIS_", batch, "_", cageChange, "_AnimalPos")
      csvFilePath <- file.path(preprocessed_data_dir, paste0(filename, "_preprocessed.csv"))
      if (file.exists(csvFilePath)) {
        temp_data <- read_delim(csvFilePath, delim = ",", show_col_types = FALSE) %>% as_tibble()
        all_preprocessed <- bind_rows(all_preprocessed, temp_data)
      }
    }
  }
  data_preprocessed <- all_preprocessed
  message("   ✓ Combined ", nrow(data_preprocessed), " rows from all batches/cage changes")
}

# Ensure there is a DateTime column to allow arranging; if missing, try common alternatives or create a fallback ordering column
if (!"DateTime" %in% colnames(data_preprocessed)) {
  dt_candidates <- intersect(c("DateTime", "Date_Time", "Timestamp", "timestamp", "datetime", "Date", "Time"), colnames(data_preprocessed))
  if (length(dt_candidates) > 0) {
    chosen_dt <- dt_candidates[1]
    data_preprocessed$DateTime <- data_preprocessed[[chosen_dt]]
    message("   ✓ Mapped existing column '", chosen_dt, "' to DateTime for ordering")
  } else {
    data_preprocessed$DateTime <- seq_len(nrow(data_preprocessed))
    message("   Warning: No Date/Time column found; created synthetic DateTime ordering column")
  }
}

animals <- unique(data_preprocessed$AnimalID)

# Directory for transition networks
trans_net_dir <- transition_networks_dir
dir.create(trans_net_dir, recursive = TRUE, showWarnings = FALSE)

# Detect position column
possible_pos_cols <- c("PositionID", "Position", "position", "Pos", "pos", "position_id", "Zone", "zone")
pos_candidates <- intersect(possible_pos_cols, colnames(data_preprocessed))

if (length(pos_candidates) == 0) {
  message("Warning: No Position column found in preprocessed data; skipping transition network generation.")
} else {
  pos_col <- if ("PositionID" %in% pos_candidates) "PositionID" else pos_candidates[1]
  message("Using position column: ", pos_col)

  processed_count <- 0L

  for (animal in animals) {
    seq_df <- data_preprocessed %>% filter(AnimalID == animal) %>% arrange(.data$DateTime)
    if (nrow(seq_df) == 0) { message("   Skipping animal ", animal, ": no rows"); next }
    if (!pos_col %in% colnames(seq_df)) { message("   Skipping animal ", animal, ": position column '", pos_col, "' missing"); next }

    seq_positions <- seq_df[[pos_col]]
    seq_positions <- seq_positions[!is.na(seq_positions)]

    if (length(seq_positions) < 2) { message("   Skipping animal ", animal, ": not enough positions (", length(seq_positions), ")"); next }

    p <- tryCatch(plot_transition_network(seq_positions, animal), error = function(e) { message("   Error plotting for ", animal, ": ", e$message); NULL })

    if (!is.null(p)) {
      fname_base <- gsub("[^A-Za-z0-9_\\-]", "_", as.character(animal))
      svg_path <- file.path(trans_net_dir, paste0("transition_network_", fname_base, ".svg"))
      tryCatch(ggsave(svg_path, plot = p, width = 10, height = 10, dpi = 300),
               error = function(e) message("   Failed to save plot for ", animal, ": ", e$message))
      if (show_plots) print(p)
    }

    edges <- data.frame(from = head(seq_positions, -1), to = tail(seq_positions, -1), stringsAsFactors = FALSE)
    g <- tryCatch(graph_from_data_frame(edges, directed = TRUE), error = function(e) { message("   Graph build error for ", animal, ": ", e$message); NULL })
    if (!is.null(g)) {
      rds_path <- file.path(trans_net_dir, paste0("network_", gsub("[^A-Za-z0-9_\\-]", "_", as.character(animal)), ".rds"))
      tryCatch(saveRDS(g, file = rds_path), error = function(e) message("   Failed to save RDS for ", animal, ": ", e$message))
    }

    processed_count <- processed_count + 1L
  }

  message("✓ Transition network plots and graphs saved for ", processed_count, " / ", length(animals), " animals.")
}

# Ensure Group and Sex columns exist in the combined preprocessed data
if (!"Group" %in% colnames(data_preprocessed)) {
  data_preprocessed <- data_preprocessed %>%
    mutate(Group = ifelse(AnimalID %in% sus_animals, "sus",
                          ifelse(AnimalID %in% con_animals, "con", "res")))
}
if (!"Sex" %in% colnames(data_preprocessed)) {
  data_preprocessed <- data_preprocessed %>%
    mutate(Sex = ifelse(Batch %in% c("B3", "B4", "B6"), "female", "male"))
}

# Ensure Phase exists (infer from ConsecActive/ConsecInactive if necessary)
if (!"Phase" %in% colnames(data_preprocessed)) {
  data_preprocessed <- data_preprocessed %>%
    mutate(
      Phase = case_when(
        !is.na(ConsecActive) & ConsecActive > 0 ~ "Active",
        !is.na(ConsecInactive) & ConsecInactive > 0 ~ "Inactive",
        TRUE ~ NA_character_
      )
    )
}

# Per-group/sex combined networks (overall, and split by phase)
for (phase_val in c(NA, "Active", "Inactive")) {
  phase_label <- if (is.na(phase_val)) "overall" else phase_val
  phase_dir <- file.path(trans_net_dir, paste0("phase_", phase_label))
  dir.create(phase_dir, recursive = TRUE, showWarnings = FALSE)

  message("Generating combined transition network plots by group and sex (phase: ", phase_label, ")...")

  for (sex_val in c("male", "female")) {
    for (group_val in c("sus", "res", "con")) {
      seq_df <- data_preprocessed %>%
        filter(Sex == sex_val, Group == group_val)

      if (!is.na(phase_val)) seq_df <- seq_df %>% filter(Phase == phase_val)

      seq_df <- seq_df %>% arrange(.data$DateTime) %>% select(AnimalID, !!sym(pos_col)) %>% filter(!is.na(!!sym(pos_col)))

      if (nrow(seq_df) == 0) { message("   Skipping group/sex/phase combination: ", sex_val, "/", group_val, "/", phase_label); next }

      combined_positions <- unlist(seq_df[[pos_col]])
      combined_positions <- combined_positions[!is.na(combined_positions)]

      if (length(combined_positions) < 2) { message("   Skipping group/sex/phase combination: ", sex_val, "/", group_val, "/", phase_label, " - not enough positions"); next }

      p <- tryCatch(plot_transition_network(combined_positions, paste0(sex_val, "_", group_val, "_", phase_label)), error = function(e) { message("   Error plotting for ", sex_val, "/", group_val, " (", phase_label, "): ", e$message); NULL })

      if (!is.null(p)) {
        fname_base <- paste0(sex_val, "_", group_val, "_", phase_label)
        svg_path <- file.path(phase_dir, paste0("transition_network_", fname_base, ".svg"))
        tryCatch(ggsave(svg_path, plot = p, width = 10, height = 10, dpi = 300),
                 error = function(e) message("   Failed to save plot for ", sex_val, "/", group_val, " (", phase_label, "): ", e$message))
        if (show_plots) print(p)
      }
    }

    # Merged network across groups for this sex, edges colored by Group (optionally filtered by phase)
    message("Generating merged transition network for sex: ", sex_val, " (phase: ", phase_label, ")")
    sex_seq_df <- data_preprocessed %>%
      filter(Sex == sex_val)
    if (!is.na(phase_val)) sex_seq_df <- sex_seq_df %>% filter(Phase == phase_val)

    sex_seq_df <- sex_seq_df %>% arrange(.data$DateTime) %>% select(AnimalID, Group, DateTime, !!sym(pos_col)) %>% filter(!is.na(!!sym(pos_col)))

    if (nrow(sex_seq_df) == 0) {
      message("   No data for sex: ", sex_val, " (phase: ", phase_label, ") — skipping merged network.")
    } else {
      sex_edges <- sex_seq_df %>%
        group_by(AnimalID) %>%
        arrange(.data$DateTime) %>%
        mutate(next_pos = lead(!!sym(pos_col))) %>%
        ungroup() %>%
        filter(!is.na(!!sym(pos_col)) & !is.na(next_pos)) %>%
        transmute(from = as.character(!!sym(pos_col)), to = as.character(next_pos), Group = as.character(Group))

      if (nrow(sex_edges) == 0) {
        message("   Not enough transitions to build merged network for sex: ", sex_val, " (phase: ", phase_label, ")")
      } else {
        edges_agg_by_group <- sex_edges %>%
          group_by(from, to, Group) %>%
          summarise(weight = n(), .groups = "drop")

        g_merged <- tryCatch({
          graph_from_data_frame(edges_agg_by_group, directed = TRUE)
        }, error = function(e) { message("   Graph build error for merged ", sex_val, " (", phase_label, "): ", e$message); NULL })

        if (!is.null(g_merged)) {
          E(g_merged)$weight <- edges_agg_by_group$weight
          E(g_merged)$group <- edges_agg_by_group$Group

          p_sex <- tryCatch({
            ggraph(g_merged, layout = "fr") +
              geom_edge_fan(aes(width = weight, color = group, alpha = weight),
                            arrow = arrow(length = unit(3, "mm"), type = "closed"),
                            show.legend = TRUE) +
              scale_edge_width(range = c(0.2, 3)) +
              scale_edge_colour_manual(values = c(sus = "#E63946", res = "grey60", con = "#457B9D")) +
              geom_node_point(aes(size = degree(g_merged, mode = "all")), color = "steelblue") +
              geom_node_text(aes(label = name), repel = TRUE, size = 3) +
              theme_void() +
              guides(alpha = "none") +
              ggtitle(paste0("Merged transition network (sex=", sex_val, ", phase=", phase_label, ")"))
          }, error = function(e) { message("   Error plotting merged network for ", sex_val, " (", phase_label, "): ", e$message); NULL })

          if (!is.null(p_sex)) {
            fname_sex <- paste0(sex_val, "_merged_by_group_", phase_label)
            svg_path_sex <- file.path(phase_dir, paste0("transition_network_", fname_sex, ".svg"))
            tryCatch(ggsave(svg_path_sex, plot = p_sex, width = 10, height = 10, dpi = 300),
                     error = function(e) message("   Failed to save merged plot for ", sex_val, " (", phase_label, "): ", e$message))
            if (show_plots) print(p_sex)
          }

          rds_path_sex <- file.path(phase_dir, paste0("network_", paste0(sex_val, "_merged_by_group_", phase_label), ".rds"))
          tryCatch(saveRDS(g_merged, file = rds_path_sex), error = function(e) message("   Failed to save merged RDS for ", sex_val, " (", phase_label, "): ", e$message))
        }
      }
    }
  } # end sex loop
} # end phase loop
message("✓ Combined transition network plots saved.")


# -------------------------------------------------------------------
# STATISTICAL COMPARISONS OF TRANSITION NETWORKS
# a) Within sex between groups (edgewise tests + network entropy)
# b) Between sex for group "con"
# c) Add comparison plots for male vs female group "con"
# -------------------------------------------------------------------

# helper: aggregate transitions for a filtered df
aggregate_transitions <- function(df, pos_col) {
  df %>%
    arrange(.data$DateTime) %>%
    group_by(AnimalID) %>%
    mutate(next_pos = lead(!!sym(pos_col))) %>%
    ungroup() %>%
    filter(!is.na(!!sym(pos_col)) & !is.na(next_pos)) %>%
    transmute(from = as.character(!!sym(pos_col)), to = as.character(next_pos))
}

# helper: build edge counts table
edge_counts_table <- function(trans_df) {
  trans_df %>% group_by(from, to) %>% summarise(weight = n(), .groups = "drop")
}

# helper: entropy of transition probability distribution (works for arbitrary edge counts)
transition_entropy <- function(edge_counts) {
  if (nrow(edge_counts) == 0) return(NA_real_)
  probs <- edge_counts$weight / sum(edge_counts$weight)
  probs <- probs[!is.na(probs) & probs > 0]
  if (length(probs) == 0) return(NA_real_)
  # compute Shannon entropy in bits without relying on calc_shannon_entropy (which may expect a fixed-length vector)
  -sum(probs * log2(probs))
}

# Helper: edgewise contingency tests between two groups (within same sex)
edgewise_comparisons <- function(edgesA, edgesB, out_prefix) {
  totalA <- sum(edgesA$weight)
  totalB <- sum(edgesB$weight)
  all_edges <- full_join(edgesA, edgesB, by = c("from", "to"), suffix = c(".A", ".B")) %>%
    replace_na(list(weight.A = 0, weight.B = 0))
  res <- all_edges %>% rowwise() %>% mutate({
    mat <- matrix(c(weight.A, totalA - weight.A, weight.B, totalB - weight.B), nrow = 2)
    test <- tryCatch({
      if (min(chisq.test(mat)$expected, na.rm = TRUE) >= 5) {
        t <- chisq.test(mat)
        list(test_type = "chisq", p_value = as.numeric(t$p.value), statistic = as.numeric(t$statistic))
      } else {
        t <- fisher.test(mat, simulate.p.value = TRUE)
        list(test_type = "fisher", p_value = as.numeric(t$p.value), statistic = NA_real_)
      }
    }, error = function(e) list(test_type = "error", p_value = NA_real_, statistic = NA_real_))
    tibble(test_type = test$test_type, p_value = test$p_value, statistic = test$statistic)
  }) %>% ungroup()
  res_adj <- res %>% mutate(p_adj = p.adjust(p_value, method = "BH"))
  write.csv(bind_cols(all_edges, res_adj), file = file.path(trans_net_dir, paste0(out_prefix, "_edgewise_tests.csv")), row.names = FALSE)
  return(bind_cols(all_edges, res_adj))
}

# Helper: permutation test comparing entropies via resampling transitions
perm_test_entropy <- function(edgesA, edgesB, n_perm = 1000, seed = 42) {
  set.seed(seed)
  vecA <- rep(paste0(edgesA$from, "->", edgesA$to), edgesA$weight)
  vecB <- rep(paste0(edgesB$from, "->", edgesB$to), edgesB$weight)
  nA <- length(vecA); nB <- length(vecB)
  pool <- c(vecA, vecB)
  obs_entropy_diff <- transition_entropy(edge_counts_table(tibble(from = sub("->.*", "", vecA), to = sub(".*->", "", vecA)))) -
    transition_entropy(edge_counts_table(tibble(from = sub("->.*", "", vecB), to = sub(".*->", "", vecB))))
  if (length(pool) < 2) return(list(obs = obs_entropy_diff, p_value = NA_real_))
  perm_diffs <- replicate(n_perm, {
    samp <- sample(pool, length(pool), replace = FALSE)
    a <- samp[1:nA]; b <- samp[(nA + 1):(nA + nB)]
    ea <- edge_counts_table(tibble(from = sub("->.*", "", a), to = sub(".*->", "", a)))
    eb <- edge_counts_table(tibble(from = sub("->.*", "", b), to = sub(".*->", "", b)))
    transition_entropy(ea) - transition_entropy(eb)
  })
  pval <- mean(abs(perm_diffs) >= abs(obs_entropy_diff))
  list(obs = obs_entropy_diff, p_value = pval, perm_diffs = perm_diffs)
}

analysis_dir <- file.path(trans_net_dir, "analysis")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

# a) Within sex, compare groups pairwise
groups <- c("sus", "res", "con")
for (sex_val in c("male", "female")) {
  message("Running within-sex group comparisons for sex = ", sex_val)
  sex_df <- data_preprocessed %>% filter(Sex == sex_val)
  if (nrow(sex_df) == 0) next
  trans_all <- aggregate_transitions(sex_df, pos_col)
  # build per-group edge counts
  per_group_edges <- sex_df %>%
    arrange(.data$DateTime) %>%
    group_by(AnimalID, Group) %>%
    mutate(next_pos = lead(!!sym(pos_col))) %>%
    ungroup() %>%
    filter(!is.na(!!sym(pos_col)) & !is.na(next_pos)) %>%
    transmute(from = as.character(!!sym(pos_col)), to = as.character(next_pos), Group = Group)
  edges_by_group <- per_group_edges %>% group_by(Group, from, to) %>% summarise(weight = n(), .groups = "drop")
  # pairwise comparisons
  combs <- combn(groups, 2, simplify = FALSE)
  entropy_summary <- list()
  for (cmp in combs) {
    g1 <- cmp[1]; g2 <- cmp[2]
    e1 <- edges_by_group %>% filter(Group == g1) %>% select(from, to, weight)
    e2 <- edges_by_group %>% filter(Group == g2) %>% select(from, to, weight)
    e1_tbl <- edge_counts_table(e1)
    e2_tbl <- edge_counts_table(e2)
    edgewise_res <- edgewise_comparisons(e1_tbl, e2_tbl, paste0("within_sex_", sex_val, "_", g1, "_vs_", g2))
    # entropy permutation test
    perm_res <- perm_test_entropy(e1_tbl, e2_tbl, n_perm = 1000)
    entropy_summary[[paste0(sex_val, "_", g1, "_vs_", g2)]] <- list(groupA = g1, groupB = g2,
                                                                   entropy_diff = perm_res$obs, p_value = perm_res$p_value)
    message("   ", sex_val, " ", g1, " vs ", g2, ": entropy diff = ", signif(perm_res$obs, 4), ", p = ", signif(perm_res$p_value, 4))
  }
  entropy_df <- bind_rows(lapply(names(entropy_summary), function(k) tibble(comparison = k, !!!entropy_summary[[k]])))
  write.csv(entropy_df, file = file.path(analysis_dir, paste0("within_sex_", sex_val, "_entropy_comparisons.csv")), row.names = FALSE)
}

# b) Between sex for group "con"
message("Comparing male vs female for group 'con'")
con_df <- data_preprocessed %>% filter(Group == "con")
if (nrow(con_df) > 0) {
  male_df <- con_df %>% filter(Sex == "male")
  female_df <- con_df %>% filter(Sex == "female")
  e_male <- aggregate_transitions(male_df, pos_col) %>% edge_counts_table()
  e_female <- aggregate_transitions(female_df, pos_col) %>% edge_counts_table()
  edgewise_consex <- edgewise_comparisons(e_male, e_female, "con_male_vs_female")
  perm_consex <- perm_test_entropy(e_male, e_female, n_perm = 2000)
  write.csv(tibble(metric = c("entropy_diff_obs", "entropy_perm_p"), value = c(perm_consex$obs, perm_consex$p_value)),
            file = file.path(analysis_dir, "con_male_vs_female_entropy_test.csv"), row.names = FALSE)
  message("   con male vs female: entropy diff = ", signif(perm_consex$obs,4), ", p = ", signif(perm_consex$p_value,4))
} else {
  message("   No 'con' data available for male/female comparison.")
}

# c) Transition plot comparison between male and female group "con" animals
message("Generating differential transition plot for 'con' males vs females")
if (nrow(con_df) > 0 && nrow(male_df) > 0 && nrow(female_df) > 0) {
  male_edges <- aggregate_transitions(male_df, pos_col) %>% edge_counts_table() %>% mutate(sex = "male")
  female_edges <- aggregate_transitions(female_df, pos_col) %>% edge_counts_table() %>% mutate(sex = "female")
  all_edges <- full_join(male_edges, female_edges, by = c("from", "to"), suffix = c(".male", ".female")) %>%
    replace_na(list(weight.male = 0, weight.female = 0)) %>%
    mutate(diff = weight.male - weight.female,
           diff_norm = diff / pmax(1, weight.male + weight.female))
  # create graph with diff attribute
  g_diff <- graph_from_data_frame(all_edges %>% select(from, to, weight.male, weight.female, diff, diff_norm), directed = TRUE)
  E(g_diff)$weight_male <- all_edges$weight.male
  E(g_diff)$weight_female <- all_edges$weight.female
  E(g_diff)$diff <- all_edges$diff
  # plot: color edges by diff (male>female red, female>male blue)
  p_diff <- tryCatch({
    ggraph(g_diff, layout = "fr") +
      geom_edge_fan(aes(width = abs(diff), color = diff), arrow = arrow(length = unit(3, "mm"), type = "closed"), show.legend = TRUE) +
      scale_edge_colour_gradient2(low = "#457B9D", mid = "grey80", high = "#E63946", midpoint = 0, name = "male - female") +
      scale_edge_width(range = c(0.2, 4)) +
      geom_node_point(aes(size = degree(g_diff, mode = "all")), color = "black") +
      geom_node_text(aes(label = name), repel = TRUE, size = 3) +
      theme_void() +
      ggtitle("Differential transitions: con males vs con females")
  }, error = function(e) { message("   Error creating differential plot: ", e$message); NULL })
  if (!is.null(p_diff)) {
    svg_path_diff <- file.path(trans_net_dir, "con_male_vs_female_difference.svg")
    tryCatch(ggsave(svg_path_diff, plot = p_diff, width = 10, height = 10, dpi = 300),
             error = function(e) message("   Failed to save differential plot: ", e$message))
    if (show_plots) print(p_diff)
  }
  # also save combined male and female network plots side-by-side
  pm <- tryCatch(plot_transition_network(rep(male_edges$from, male_edges$weight.male) %>% c(rep(male_edges$to, male_edges$weight.male)), "con_male"), error = function(e) NULL)
  pf <- tryCatch(plot_transition_network(rep(female_edges$from, female_edges$weight.female) %>% c(rep(female_edges$to, female_edges$weight.female)), "con_female"), error = function(e) NULL)
  if (!is.null(pm)) {
    tryCatch(ggsave(file.path(trans_net_dir, "con_male_network.svg"), plot = pm, width = 10, height = 10, dpi = 300), error = function(e) NULL)
  }
  if (!is.null(pf)) {
    tryCatch(ggsave(file.path(trans_net_dir, "con_female_network.svg"), plot = pf, width = 10, height = 10, dpi = 300), error = function(e) NULL)
  }
  write.csv(all_edges, file = file.path(analysis_dir, "con_male_vs_female_edge_counts.csv"), row.names = FALSE)
  message("   Differential plots and tables saved.")
} else {
  message("   Insufficient con male/female data to create comparison plots.")
}

message("✓ Transition network statistical comparisons complete. Results in: ", analysis_dir)


# ===================================================================
# ENTROPY DYNAMICS OVER TIME (Per animal)
# ===================================================================
message("Plotting entropy velocity over time...")
entropy_dynamics <- consec_animal_entropy_halfhour %>%
  group_by(AnimalID) %>%
  arrange(HalfHour) %>%
  mutate(entropy_change = animalEntropy - lag(animalEntropy),
         entropy_velocity = abs(entropy_change))

p_entropy <- ggplot(entropy_dynamics, aes(x = HalfHour, y = entropy_velocity, group = AnimalID, color = AnimalID)) +
  geom_line(alpha = 0.5) +
  labs(title = "Entropy Velocity over Time", x = "Half-Hour Period", y = "Change in Entropy") +
  theme_minimal() +
  theme(legend.position = "none")

# Save or show
if (save_plots) {
  ggsave(file.path(combined_plots_dir, "entropy_velocity_over_time.svg"), p_entropy, width=12, height=6, dpi=300)
}
if (show_plots) print(p_entropy)

# ===================================================================
# SAVE METADATA AND SESSION INFO
# ===================================================================
processing_end_time <- Sys.time()

metadata <- list(
  parameters = list(
    load_existing_data = load_existing_data,
    filter_cc4_late_phases = filter_cc4_late_phases,
    cc4_max_active_phase = cc4_max_active_phase,
    cc4_max_inactive_phase = cc4_max_inactive_phase,
    analyze_by_halfhour = analyze_by_halfhour,
    exclude_homecage = exclude_homecage
  ),
  processing_summary = list(
    batches_processed = batches,
    cage_changes_processed = cageChanges,
    total_rows_preprocessed = if(exists("data_preprocessed")) nrow(data_preprocessed) else NA
  ),
  start_time = as.character(processing_start_time),
  end_time = as.character(processing_end_time),
  git_commit = ifelse(nzchar(Sys.getenv("GIT_COMMIT")), Sys.getenv("GIT_COMMIT"), NA),
  R_version = R.Version()$version.string,
  session_info = capture.output(sessionInfo())
)

# Save JSON metadata
metadata_path <- file.path(metadata_dir, "processing_metadata.json")
write_json(metadata, metadata_path, pretty=TRUE)

# Save session info text
session_info_path <- file.path(metadata_dir, "session_info.txt")
writeLines(metadata$session_info, session_info_path)

# Log completion
log_message("Saved metadata and session info.")
message("=======================================================================")
message(" ENTIRE ENTROPY ANALYSIS PIPELINE COMPLETE! ")
message("=======================================================================")
message("=======================================================================")
message("PUBLICATION-READY OUTPUT COMPLETE")
message("=======================================================================")
message("Publication tables: ", pub_tables_dir)
message("Single-panel figures: ", pub_figure_panels_dir)
message("Multi-panel figures: ", pub_multipanel_dir)