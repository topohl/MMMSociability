#' @title Analysis of Shannon Entropy in Animal Positioning
#'
#' @description
#' This script conducts a comprehensive analysis of Shannon entropy to
#' evaluate the distribution and diversity of animal positions within a
#' specified environment, analyzing both by phase and by half-hour periods.
#'
#' @details
#' This script can either process raw data or load previously processed results.
#' Set load_existing_data = TRUE to skip processing and load existing CSVs.
#'
#' @date October 2025
#' @authors Tobias Pohl, Anja Magister

# ===================================================================
# LOAD REQUIRED PACKAGES
# ===================================================================
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(readr, dplyr, lubridate, tibble, purrr, ggplot2, reshape2,
               scales, stringr, gridExtra, lme4, emmeans, ggpubr, tidyr, GGally,
               igraph, ggraph, jsonlite)

# ===================================================================
# RUN MODE CONFIGURATION
# ===================================================================
load_existing_data <- TRUE
show_plots <- FALSE
save_plots <- TRUE
save_tables <- TRUE
exclude_homecage <- TRUE
analyze_by_halfhour <- TRUE

filter_cc4_late_phases <- TRUE
cc4_max_active_phase <- 2
cc4_max_inactive_phase <- 2

# ===================================================================
# PATHS WITH OPTIMIZED FOLDER STRUCTURE
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
# INITIALIZE LOGGING
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
# LOAD CUSTOM FUNCTIONS AND ANIMAL LISTS
# ===================================================================
source(file.path(working_directory, "E9_SIS_AnimalPos-functions.R"))

sus_animals <- readLines(file.path(raw_data_dir, "sus_animals.csv"))
con_animals <- readLines(file.path(raw_data_dir, "con_animals.csv"))

# ===================================================================
# OPTION 1: LOAD EXISTING PROCESSED DATA
# ===================================================================
if (load_existing_data) {
  message("=======================================================================")
  message("LOADING EXISTING PROCESSED DATA")
  message("=======================================================================")
  
  files_to_load <- list(
    cagePosProb = file.path(combined_tables_dir, "phase_based", "all_batches_all_cageChanges_cagePosProb.csv"),
    cagePosEntropy = file.path(combined_tables_dir, "phase_based", "all_batches_all_cageChanges_cagePosEntropy.csv"),
    animalPosEntropy = file.path(combined_tables_dir, "phase_based", "all_batches_all_cageChanges_animalPosEntropy.csv")
  )
  
  if (analyze_by_halfhour) {
    files_to_load$cagePosEntropy_halfhour <- file.path(combined_tables_dir, "halfhour_based", "all_batches_all_cageChanges_cagePosEntropy_halfhour.csv")
    files_to_load$animalPosEntropy_halfhour <- file.path(combined_tables_dir, "halfhour_based", "all_batches_all_cageChanges_animalPosEntropy_halfhour.csv")
  }
  
  all_files_exist <- all(file.exists(unlist(files_to_load)))
  
  if (all_files_exist) {
    message("✓ All required files found. Loading data...")
    
    all_cagePosProb <- read_csv(files_to_load$cagePosProb, show_col_types = FALSE)
    message("   ✓ Loaded: all_batches_all_cageChanges_cagePosProb.csv")
    
    all_cagePosEntropy <- read_csv(files_to_load$cagePosEntropy, show_col_types = FALSE)
    message("   ✓ Loaded: all_batches_all_cageChanges_cagePosEntropy.csv")
    
    all_animalPosEntropy <- read_csv(files_to_load$animalPosEntropy, show_col_types = FALSE)
    message("   ✓ Loaded: all_batches_all_cageChanges_animalPosEntropy.csv")
    
    if (analyze_by_halfhour) {
      all_cagePosEntropy_halfhour <- read_csv(files_to_load$cagePosEntropy_halfhour, show_col_types = FALSE)
      message("   ✓ Loaded: all_batches_all_cageChanges_cagePosEntropy_halfhour.csv")
    
      all_animalPosEntropy_halfhour <- read_csv(files_to_load$animalPosEntropy_halfhour, show_col_types = FALSE)
      message("   ✓ Loaded: all_batches_all_cageChanges_animalPosEntropy_halfhour.csv")
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
# OPTION 2: PROCESS DATA FROM SCRATCH
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
    
    # Loop: Phase-based entropy and probability calculations
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
    
    # Half-hour processing
    if (analyze_by_halfhour) {
      for (system_id in unique_systems) {
        system_data <- data_preprocessed %>% filter(System == system_id) %>% as_tibble()
        animal_ids <- unique(system_data$AnimalID)
        system_complete <- length(animal_ids) >= 4
        while (length(animal_ids) < 4) animal_ids <- append(animal_ids, NA)
        
        for (halfhour in sort(unique(system_data$HalfHoursElapsed))) {
          data_halfhour <- system_data %>% filter(HalfHoursElapsed == halfhour) %>% as_tibble()
          if (nrow(data_halfhour) == 0) next
          
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
          
          animal_position_probability <- tibble(AnimalID=rep(animal_ids, each=8),
                                                Position=rep(1:8, length.out=32),
                                                Seconds=0, SumPercentage=0, Prob=0)
          
          animal_list <- initialize_animal_positions(animal_ids, data_halfhour, animal_list)
          initial_time <- animal_list[[1]][[2]]
          current_row <- 5
          total_rows <- nrow(data_halfhour) + 1
          
          while (current_row != total_rows && current_row < total_rows) {
            previous_positions <- animal_list
            animal_list <- update_animal_list(animal_ids, animal_list, data_halfhour, initial_time, current_row)
            elapsed_seconds <- animal_list[["data_temp"]][["elapsed_seconds"]]
            if (system_complete) {
              cage_position_probability <- update_cage_position_probability(previous_positions, animal_list, cage_position_probability, elapsed_seconds)
            }
            animal_position_probability <- update_animal_position_probability(previous_positions, animal_list, animal_position_probability, elapsed_seconds)
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
            cage_prob_vec <- sapply(cage_position_probability, function(x) x[4])
            cage_prob_vec <- cage_prob_vec[!is.na(cage_prob_vec) & !is.nan(cage_prob_vec)]
            cage_shannon_entropy <- NA
            if (length(cage_prob_vec) > 0 && sum(cage_prob_vec) > 0) {
              cage_prob_vec <- cage_prob_vec / sum(cage_prob_vec)
              cage_shannon_entropy <- calc_shannon_entropy(cage_prob_vec)
            }
            
            cagePosEntropy_halfhour <- cagePosEntropy_halfhour %>%
              add_row(Batch=batch, Sex=sex, System=system_id, CageChange=cageChange,
                      HalfHour=halfhour, CageEntropy=cage_shannon_entropy)
          }
          
          for (animal in animal_ids) {
            if (is.na(animal)) next
            
            animal_prob_vec <- animal_position_probability %>%
              filter(AnimalID == animal) %>%
              pull(Prob)
            
            animal_prob_vec <- animal_prob_vec[!is.na(animal_prob_vec) & !is.nan(animal_prob_vec)]
            
            animal_shannon_entropy <- NA
            if (length(animal_prob_vec) > 0 && sum(animal_prob_vec) > 0) {
              animal_prob_vec <- animal_prob_vec / sum(animal_prob_vec)
              animal_shannon_entropy <- calc_shannon_entropy(animal_prob_vec)
            }
            
            animalPosEntropy_halfhour <- animalPosEntropy_halfhour %>%
              add_row(Batch=batch, Sex=sex, System=system_id, CageChange=cageChange,
                      HalfHour=halfhour, AnimalID=animal, animalEntropy=animal_shannon_entropy)
          }
        }
      }
    }
    
    # Assign Group columns
    animalPosEntropy <- animalPosEntropy %>%
      mutate(Group = ifelse(AnimalID %in% sus_animals, "sus",
                            ifelse(AnimalID %in% con_animals, "con", "res")))
    if (analyze_by_halfhour) {
      animalPosEntropy_halfhour <- animalPosEntropy_halfhour %>%
        mutate(Group = ifelse(AnimalID %in% sus_animals, "sus",
                              ifelse(AnimalID %in% con_animals, "con", "res")))
    }
    
    # Plot phase-based entropy per system and phase
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

    
    # Plot half-hour entropy per system and halfhour period
if (analyze_by_halfhour) {
  for (system_id in unique_systems) {
    plot_folder <- file.path(saving_directory, "plots", "batch_cagechange", batch, cageChange, "halfhour_based")
    dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)
    
    halfhours <- sort(unique(animalPosEntropy_halfhour$HalfHour))
    
    for (halfhour in halfhours) {
      plot_data <- animalPosEntropy %>% filter(System == system_id, Phase == phase_code) %>%
  filter(!is.na(animalEntropy), !is.na(Group))  # Remove rows with missing essential data

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

full_plot_path <- file.path(plot_folder, paste0("animal_entropy_", system_id, "_", phase_code, ".svg"))
message("Trying to save plot at: ", full_plot_path)
print(p_phase_plot)

tryCatch({
  ggsave(filename = full_plot_path, plot = p_phase_plot, width = 10, height = 7, dpi = 300)
  message("Plot saved successfully.")
}, error = function(e) {
  message("Error saving plot: ", e$message)
})

    }
  }
}

    
    # Save to CSV
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
    
    write.csv(all_cagePosProb, file = file.path(combined_tables_dir, "phase_based", "all_batches_all_cageChanges_cagePosProb.csv"), row.names = FALSE)
    message("   ✓ Saved: all_batches_all_cageChanges_cagePosProb.csv")
    
    write.csv(all_cagePosEntropy, file = file.path(combined_tables_dir, "phase_based", "all_batches_all_cageChanges_cagePosEntropy.csv"), row.names = FALSE)
    message("   ✓ Saved: all_batches_all_cageChanges_cagePosEntropy.csv")
    
    write.csv(all_animalPosEntropy, file = file.path(combined_tables_dir, "phase_based", "all_batches_all_cageChanges_animalPosEntropy.csv"), row.names = FALSE)
    message("   ✓ Saved: all_batches_all_cageChanges_animalPosEntropy.csv")
    
    if(analyze_by_halfhour) {
      write.csv(all_cagePosEntropy_halfhour, file = file.path(combined_tables_dir, "halfhour_based", "all_batches_all_cageChanges_cagePosEntropy_halfhour.csv"), row.names = FALSE)
      message("   ✓ Saved: all_batches_all_cageChanges_cagePosEntropy_halfhour.csv")
      
      write.csv(all_animalPosEntropy_halfhour, file = file.path(combined_tables_dir, "halfhour_based", "all_batches_all_cageChanges_animalPosEntropy_halfhour.csv"), row.names = FALSE)
      message("   ✓ Saved: all_batches_all_cageChanges_animalPosEntropy_halfhour.csv")
    }
    message("      Location: ", combined_tables_dir)
  }
  
  message("=======================================================================")
  message("DATA PROCESSING COMPLETE")
  message("=======================================================================")
}

# ===================================================================
# CONSECUTIVE ENTROPY PROCESSING (PHASE-BASED)
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
# CONSECUTIVE ENTROPY PROCESSING (HALF-HOUR)
# ===================================================================
if (analyze_by_halfhour) {
  message("=======================================================================")
  message("PREPROCESSING CONSECUTIVE ENTROPY DATA (HALF-HOUR)")
  message("=======================================================================")
  
  consec_cage_entropy_halfhour <- all_cagePosEntropy_halfhour
  consec_animal_entropy_halfhour <- all_animalPosEntropy_halfhour
  
  entropy_halfhour_names <- c("consec_cage_entropy_halfhour", "consec_animal_entropy_halfhour")
  max_halfhour <- 0
  
  for (i in seq_along(entropy_halfhour_names)) {
    name <- entropy_halfhour_names[i]
    entropy_list <- if (name == "consec_cage_entropy_halfhour") consec_cage_entropy_halfhour else consec_animal_entropy_halfhour
    
    if (name == "consec_animal_entropy_halfhour") {
      entropy_list <- entropy_list %>%
        mutate(Group = ifelse(AnimalID %in% sus_animals, "sus",
                              ifelse(AnimalID %in% con_animals, "con", "res")))
    }
    entropy_list <- entropy_list %>% mutate(ConsecHalfHour = HalfHour)
    
    for (change in c("CC1", "CC2", "CC3", "CC4")) {
      if (change != "CC1") {
        entropy_list <- entropy_list %>%
          mutate(ConsecHalfHour = ifelse(CageChange == change, HalfHour + max_halfhour, ConsecHalfHour))
      }
      max_halfhour <- entropy_list %>% filter(CageChange == change) %>% pull(ConsecHalfHour) %>% max(na.rm = TRUE)
    }
    
    if (name == "consec_animal_entropy_halfhour") {
      entropy_list <- entropy_list[c("CageChange", "Batch", "System", "AnimalID", "Sex", "Group", "HalfHour", "ConsecHalfHour", "animalEntropy")]
      consec_animal_entropy_halfhour <- entropy_list
    } else {
      entropy_list <- entropy_list[c("CageChange", "Batch", "System", "Sex", "HalfHour", "ConsecHalfHour", "CageEntropy")]
      consec_cage_entropy_halfhour <- entropy_list
    }
  }
}

if (analyze_by_halfhour && filter_cc4_late_phases) {
  message("=======================================================================")
  message("FILTERING CC4 HALF-HOURS BASED ON PHASE INFORMATION")
  message("=======================================================================")
  
  if (!exists("data_preprocessed") || is.null(data_preprocessed) || nrow(data_preprocessed) == 0) {
    message("   Warning: Cannot filter CC4 half-hours by phase - preprocessed data not available")
    message("   Run with load_existing_data = FALSE to enable or use cc4_halfhour_cutoff instead")
  } else {
    cc4_valid_halfhours <- data_preprocessed %>%
      filter(CageChange == "CC4") %>%
      filter(ConsecActive <= cc4_max_active_phase | ConsecInactive <= cc4_max_inactive_phase) %>%
      select(Batch, System, CageChange, HalfHour = HalfHoursElapsed) %>%
      distinct() %>%
      mutate(keep_halfhour = TRUE)
    
    consec_animal_entropy_halfhour <- consec_animal_entropy_halfhour %>%
      left_join(cc4_valid_halfhours, by = c("Batch", "System", "CageChange", "HalfHour")) %>%
      filter(CageChange != "CC4" | !is.na(keep_halfhour)) %>%
      select(-keep_halfhour)
    
    consec_cage_entropy_halfhour <- consec_cage_entropy_halfhour %>%
      left_join(cc4_valid_halfhours, by = c("Batch", "System", "CageChange", "HalfHour")) %>%
      filter(CageChange != "CC4" | !is.na(keep_halfhour)) %>%
      select(-keep_halfhour)
    
    message("✓ Filtered CC4 half-hours based on phase information")
    message("  Remaining half-hour observations: ", nrow(consec_animal_entropy_halfhour))
  }
}

# ===================================================================
# SAVE CONSECUTIVE ENTROPY TABLES
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
}

# ===================================================================
# EXPLORATION METRICS CALCULATION
# ===================================================================
message("=======================================================================")
message("PROCESSING EXPLORATION METRICS")
message("=======================================================================")

consec_animal_entropy <- consec_animal_entropy %>%
  mutate(
    explored = !is.na(animalEntropy),
    exploration_entropy = ifelse(is.na(animalEntropy), 0, animalEntropy),
    exploration_category = case_when(
      is.na(animalEntropy) ~ "No exploration",
      animalEntropy == 0 ~ "Minimal (1 position)",
      animalEntropy < 1.5 ~ "Low diversity",
      animalEntropy < 2.5 ~ "Moderate diversity",
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
        animalEntropy < 1.5 ~ "Low diversity",
        animalEntropy < 2.5 ~ "Moderate diversity",
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
# SAVE EXPLORATION METRICS TABLES
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
}

# ===================================================================
# GENERATE PLOTS
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
    message("Generating half-hour entropy plots with consecutive numbering...")
    cage_boundaries <- consec_animal_entropy_halfhour %>%
      group_by(CageChange) %>%
      summarise(max_consec = max(ConsecHalfHour, na.rm = TRUE)) %>%
      pull(max_consec) %>%
      cumsum()
    if (length(cage_boundaries) > 0) {
      cage_boundaries <- cage_boundaries[-length(cage_boundaries)]
    }
    if (length(cage_boundaries) > 0) {
      boundary_labels <- data.frame(
        x = cage_boundaries,
        label = paste("CC", 1:length(cage_boundaries), "→", 2:(length(cage_boundaries) + 1))
      )
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
# ADDITIONAL EXPLORATION & ADVANCED PLOTS
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

animals <- unique(data_preprocessed$AnimalID)

# Directory for transition networks
trans_net_dir <- transition_networks_dir
# Generate and save transition networks for each animal
for (animal in animals) {
  seq_positions <- data_preprocessed %>%
    filter(AnimalID == animal) %>%
    arrange(DateTime) %>%
    pull(Position)

  p <- plot_transition_network(seq_positions, animal)

  # Save plot
  filename <- file.path(trans_net_dir, paste0("transition_network_", animal, ".svg"))
  ggsave(filename, plot = p, width = 10, height = 10, dpi = 300)

  # Optionally show plot
  if (show_plots) print(p)

  # Save graph object as RDS
  edges <- data.frame(from = head(seq_positions, -1), to = tail(seq_positions, -1))
  g <- graph_from_data_frame(edges, directed = TRUE)
  saveRDS(g, file = file.path(trans_net_dir, paste0("network_", animal, ".rds")))
}

message("✓ Transition network plots and graphs saved for ", length(animals), " animals.")

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
