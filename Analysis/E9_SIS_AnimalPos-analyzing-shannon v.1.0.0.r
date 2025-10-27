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
               scales, stringr, gridExtra)

# ===================================================================
# RUN MODE CONFIGURATION
# ===================================================================
load_existing_data <- TRUE
show_plots <- FALSE
save_plots <- TRUE
save_tables <- TRUE
exclude_homecage <- FALSE
analyze_by_halfhour <- TRUE

# Filter late CC4 phases (after grid installation)
filter_cc4_late_phases <- TRUE  # Set FALSE to keep all CC4 phases
cc4_max_active_phase <- 2       # Keep only A1, A2
cc4_max_inactive_phase <- 2     # Keep only I2

# ===================================================================
# PATHS
# ===================================================================
# Alternative local directory
working_directory <- "D:/MMMSociability"
saving_directory <- "D:/MMMSociability"

# Define base directories for entropy analysis
#plots_base <- "/plots/Entropy"
#tables_base <- "/tables/Entropy"
plots_base <- paste0("/plots/", ifelse(exclude_homecage, "noHomeCage", "withHomeCage"))
tables_base <- paste0("/tables/", ifelse(exclude_homecage, "noHomeCage", "withHomeCage"))

# ===================================================================
# LOAD CUSTOM FUNCTIONS AND ANIMAL LISTS
# ===================================================================
source("C:/Users/Tobias Pohl/Documents/GitHub/MMMSociability/Functions/E9_SIS_AnimalPos-functions.R")

# Load lists of susceptible and control animals for analysis
sus_animals <- readLines(paste0(working_directory, "/raw_data/sus_animals.csv"))
con_animals <- readLines(paste0(working_directory, "/raw_data/con_animals.csv"))

# ===================================================================
# OPTION 1: LOAD EXISTING PROCESSED DATA
# ===================================================================
if (load_existing_data) {
  message("=======================================================================")
  message("LOADING EXISTING PROCESSED DATA")
  message("=======================================================================")
  
  combined_tables_dir <- paste0(tables_base, "/combined")
  
  # Define files to load
  files_to_load <- list(
    cagePosProb = paste0(saving_directory, combined_tables_dir, "/all_batches_all_cageChanges_cagePosProb.csv"),
    cagePosEntropy = paste0(saving_directory, combined_tables_dir, "/all_batches_all_cageChanges_cagePosEntropy.csv"),
    animalPosEntropy = paste0(saving_directory, combined_tables_dir, "/all_batches_all_cageChanges_animalPosEntropy.csv")
  )
  
  if (analyze_by_halfhour) {
    files_to_load$cagePosEntropy_halfhour <- paste0(saving_directory, combined_tables_dir, 
                                                     "/all_batches_all_cageChanges_cagePosEntropy_halfhour.csv")
    files_to_load$animalPosEntropy_halfhour <- paste0(saving_directory, combined_tables_dir, 
                                                       "/all_batches_all_cageChanges_animalPosEntropy_halfhour.csv")
  }
  
  # Check if all files exist
  all_files_exist <- all(file.exists(unlist(files_to_load)))
  
  if (all_files_exist) {
    message("✓ All required files found. Loading data...")
    
    # Load phase-based data
    all_cagePosProb <- read_csv(files_to_load$cagePosProb, show_col_types = FALSE)
    message("   ✓ Loaded: all_batches_all_cageChanges_cagePosProb.csv")
    
    all_cagePosEntropy <- read_csv(files_to_load$cagePosEntropy, show_col_types = FALSE)
    message("   ✓ Loaded: all_batches_all_cageChanges_cagePosEntropy.csv")
    
    all_animalPosEntropy <- read_csv(files_to_load$animalPosEntropy, show_col_types = FALSE)
    message("   ✓ Loaded: all_batches_all_cageChanges_animalPosEntropy.csv")
    
    # Load half-hour data
    if (analyze_by_halfhour) {
      all_cagePosEntropy_halfhour <- read_csv(files_to_load$cagePosEntropy_halfhour, show_col_types = FALSE)
      message("   ✓ Loaded: all_batches_all_cageChanges_cagePosEntropy_halfhour.csv")
      
      all_animalPosEntropy_halfhour <- read_csv(files_to_load$animalPosEntropy_halfhour, show_col_types = FALSE)
      message("   ✓ Loaded: all_batches_all_cageChanges_animalPosEntropy_halfhour.csv")
    }
    
    message("=======================================================================")
    message("DATA LOADING COMPLETE")
    message("=======================================================================")
    message("")
    
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
  
  # Initialize global result tibbles
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
  
  # Initialize half-hour tibbles
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
  
  # ===================================================================
  # MAIN ANALYSIS LOOP
  # ===================================================================
  batches <- c("B1", "B2", "B3", "B4", "B5", "B6")
  cageChanges <- c("CC1", "CC2", "CC3", "CC4")
  
  for (batch in batches) {
    
    # Define sex based on batch
    sex <- ifelse(batch %in% c("B3", "B4", "B6"), "female", "male")
    
    for (cageChange in cageChanges) {
      
      message("=======================================================================")
      message("PROCESSING: ", batch, " ", cageChange)
      message("=======================================================================")
      
      # Create folder structure
      current_plots_dir <- paste0(plots_base, "/", batch, "/", cageChange)
      current_tables_dir <- paste0(tables_base, "/", batch, "/", cageChange)
      
      plots_dir_phase <- paste0(current_plots_dir, "/by_phase")
      tables_dir_phase <- paste0(current_tables_dir, "/by_phase")
      tables_dir_metadata <- paste0(current_tables_dir, "/metadata")
      
      dir.create(paste0(saving_directory, plots_dir_phase), recursive = TRUE, showWarnings = FALSE)
      dir.create(paste0(saving_directory, tables_dir_phase), recursive = TRUE, showWarnings = FALSE)
      dir.create(paste0(saving_directory, tables_dir_metadata), recursive = TRUE, showWarnings = FALSE)
      
      if (analyze_by_halfhour) {
        plots_dir_halfhour <- paste0(current_plots_dir, "/by_halfhour")
        tables_dir_halfhour <- paste0(current_tables_dir, "/by_halfhour")
        
        dir.create(paste0(saving_directory, plots_dir_halfhour), recursive = TRUE, showWarnings = FALSE)
        dir.create(paste0(saving_directory, tables_dir_halfhour), recursive = TRUE, showWarnings = FALSE)
      }
      
      message("✓ Created directory structure")
      
      # Initialize local result tibbles
      cagePosProb <- tibble(
        Batch = character(), System = character(), CageChange = character(),
        Phase = character(), Position = numeric(), Probability = numeric()
      )
      
      cagePosEntropy <- tibble(
        Batch = character(), Sex = character(), System = character(),
        CageChange = character(), Phase = character(), CageEntropy = numeric()
      )
      
      animalPosEntropy <- tibble(
        Batch = character(), Sex = character(), System = character(),
        CageChange = character(), Phase = character(), AnimalID = character(),
        animalEntropy = numeric()
      )
      
      if (analyze_by_halfhour) {
        cagePosEntropy_halfhour <- tibble(
          Batch = character(), Sex = character(), System = character(),
          CageChange = character(), HalfHour = numeric(), CageEntropy = numeric()
        )
        
        animalPosEntropy_halfhour <- tibble(
          Batch = character(), Sex = character(), System = character(),
          CageChange = character(), HalfHour = numeric(), AnimalID = character(),
          animalEntropy = numeric()
        )
      }
      
      # Read preprocessed data
      filename <- paste0("E9_SIS_", batch, "_", cageChange, "_AnimalPos")
      csvFilePath <- paste0(working_directory, "/preprocessed_data/", filename, "_preprocessed.csv")
      
      if (!file.exists(csvFilePath)) {
        message("✗ File not found: ", csvFilePath)
        message("   Skipping ", batch, " ", cageChange)
        next
      }
      
      data_preprocessed <- read_delim(csvFilePath, delim = ",", show_col_types = FALSE) %>%
        as_tibble()
      
      message("✓ Loaded preprocessed data: ", nrow(data_preprocessed), " rows")
      
      # Define systems and phases
      unique_systems <- unique(data_preprocessed$System)
      unique_systems <- str_sort(unique_systems)
      
      phases <- c("Active", "Inactive")
      
      active_phase_count <- unique(data_preprocessed$ConsecActive)
      inactive_phase_count <- unique(data_preprocessed$ConsecInactive)
      
      active_phase_count <- active_phase_count[!active_phase_count %in% 0]
      inactive_phase_count <- inactive_phase_count[!inactive_phase_count %in% c(0, 1, max(inactive_phase_count))]
      
      if (cageChange == "CC4") {
        active_phase_count <- active_phase_count[active_phase_count <= 2]
        inactive_phase_count <- inactive_phase_count[inactive_phase_count <= 2]
      }
      
      active_phase_count <- sort(active_phase_count)
      inactive_phase_count <- sort(inactive_phase_count)
      
      if (analyze_by_halfhour) {
        halfhour_periods <- sort(unique(data_preprocessed$HalfHoursElapsed))
      }
      
      unique_animals <- unique(data_preprocessed$AnimalID)
      
      message("✓ Found ", length(unique_systems), " systems, ", 
              length(active_phase_count), " active phases, ",
              length(inactive_phase_count), " inactive phases")
      
      # ===================================================================
      # PHASE-BASED ANALYSIS
      # ===================================================================
      message("Starting PHASE-BASED entropy analysis...")
      
      for (system_id in unique_systems) {
        
        message("   Processing system: ", system_id)
        
        system_data <- data_preprocessed %>%
          filter(System == system_id) %>%
          as_tibble()
        
        animal_ids <- unique(system_data$AnimalID)
        system_complete <- ifelse(length(animal_ids) < 4, FALSE, TRUE)
        
        while (length(animal_ids) < 4) {
          animal_ids <- append(animal_ids, NA)
        }
        
        for (phase in phases) {
          
          phase_numbers <- if (phase == "Active") active_phase_count else inactive_phase_count
          
          for (phase_number in phase_numbers) {
            
            message("      ", phase, " phase ", phase_number)
            
            data_system_phase <- system_data %>%
              filter(ConsecActive == ifelse(phase == "Active", phase_number, 0)) %>%
              filter(ConsecInactive == ifelse(phase == "Inactive", phase_number, 0)) %>%
              as_tibble()
            
            if (nrow(data_system_phase) == 0) {
              message("         No data found. Skipping...")
              next
            }
            
            # Initialize lists
            animal_list <- list(
              "animal_1" = list(name = "", time = "", position = 0),
              "animal_2" = list(name = "", time = "", position = 0),
              "animal_3" = list(name = "", time = "", position = 0),
              "animal_4" = list(name = "", time = "", position = 0),
              "data_temp" = list(elapsed_seconds = 0, current_row = 0)
            )
            
            cage_position_probability <- list(
              c(1, 0, 0, 0), c(2, 0, 0, 0), c(3, 0, 0, 0), c(4, 0, 0, 0),
              c(5, 0, 0, 0), c(6, 0, 0, 0), c(7, 0, 0, 0), c(8, 0, 0, 0)
            )
            
            animal_position_probability <- tibble(
              AnimalID = rep(animal_ids, each = 8),
              Position = rep(1:8, length.out = 32),
              Seconds = 0,
              SumPercentage = 0,
              Prob = 0
            )
            
            # Calculate probabilities
            animal_list <- initialize_animal_positions(animal_ids, data_system_phase, animal_list)
            
            initial_time <- animal_list[[1]][[2]]
            current_row <- 5
            elapsed_seconds <- 0
            total_rows <- nrow(data_system_phase) + 1
            
            while (current_row != total_rows && current_row < total_rows) {
              
              previous_animal_positions <- animal_list
              
              animal_list <- update_animal_list(animal_ids, animal_list, data_system_phase, 
                                                initial_time, current_row)
              
              elapsed_seconds <- animal_list[["data_temp"]][["elapsed_seconds"]]
              
              if (system_complete) {
                cage_position_probability <- update_cage_position_probability(
                  previous_animal_positions, animal_list, cage_position_probability, elapsed_seconds
                )
              }
              
              animal_position_probability <- update_animal_position_probability(
                previous_animal_positions, animal_list, animal_position_probability, elapsed_seconds
              )
              
              current_row <- animal_list[["data_temp"]][["current_row"]]
              initial_time <- animal_list[[1]][[2]]
            }
            
            # Calculate final probabilities
            for (i in 1:8) {
              if (cage_position_probability[[i]][[2]] > 0) {
                cage_position_probability[[i]][[4]] <- cage_position_probability[[i]][[3]] / 
                                                        cage_position_probability[[i]][[2]]
              } else {
                cage_position_probability[[i]][[4]] <- 0
              }
            }
            
            animal_position_probability <- animal_position_probability %>%
              mutate(Prob = ifelse(Seconds > 0, SumPercentage / Seconds, 0))
            
            # Store cage position probabilities
            if (system_complete) {
              for (i in 1:8) {
                p <- paste0(substr(phase, 1, 1), phase_number)
                cagePosProb <- cagePosProb %>%
                  add_row(Batch = batch, System = system_id, CageChange = cageChange,
                          Phase = p, Position = i,
                          Probability = cage_position_probability[[i]][[4]])
              }
            }
            
            # Calculate Shannon entropy
            if (system_complete) {
              cage_prob_vec <- sapply(cage_position_probability, function(x) x[4])
              cage_prob_vec <- cage_prob_vec[!is.na(cage_prob_vec) & !is.nan(cage_prob_vec)]
              
              if (length(cage_prob_vec) > 0 && sum(cage_prob_vec) > 0) {
                cage_prob_vec <- cage_prob_vec / sum(cage_prob_vec)
                cage_shannon_entropy <- calc_shannon_entropy(cage_prob_vec)
              } else {
                cage_shannon_entropy <- NA
              }
              
              cagePosEntropy <- cagePosEntropy %>%
                add_row(Batch = batch, Sex = sex, System = system_id,
                        CageChange = cageChange,
                        Phase = paste0(substr(phase, 1, 1), phase_number),
                        CageEntropy = cage_shannon_entropy)
            }
            
            # Animal entropy
            for (animal in animal_ids) {
              if (is.na(animal)) next
              
              animal_prob_vec <- animal_position_probability %>%
                filter(AnimalID == animal) %>%
                pull(Prob)
              
              animal_prob_vec <- animal_prob_vec[!is.na(animal_prob_vec) & !is.nan(animal_prob_vec)]
              
              if (length(animal_prob_vec) > 0 && sum(animal_prob_vec) > 0) {
                animal_prob_vec <- animal_prob_vec / sum(animal_prob_vec)
                animal_shannon_entropy <- calc_shannon_entropy(animal_prob_vec)
              } else {
                animal_shannon_entropy <- NA
              }
              
              animalPosEntropy <- animalPosEntropy %>%
                add_row(Batch = batch, Sex = sex, System = system_id,
                        CageChange = cageChange,
                        Phase = paste0(substr(phase, 1, 1), phase_number),
                        AnimalID = animal,
                        animalEntropy = animal_shannon_entropy)
            }
          }
        }
      }
      
      # ===================================================================
      # HALF-HOUR ANALYSIS
      # ===================================================================
      if (analyze_by_halfhour) {
        message("Starting HALF-HOUR entropy analysis...")
        
        for (system_id in unique_systems) {
          
          message("   Processing system: ", system_id)
          
          system_data <- data_preprocessed %>%
            filter(System == system_id) %>%
            as_tibble()
          
          animal_ids <- unique(system_data$AnimalID)
          system_complete <- ifelse(length(animal_ids) < 4, FALSE, TRUE)
          
          while (length(animal_ids) < 4) {
            animal_ids <- append(animal_ids, NA)
          }
          
          for (halfhour_period in halfhour_periods) {
            
            message("      Half-hour period: ", halfhour_period)
            
            data_system_halfhour <- system_data %>%
              filter(HalfHoursElapsed == halfhour_period) %>%
              as_tibble()
            
            if (nrow(data_system_halfhour) == 0) {
              message("         No data found. Skipping...")
              next
            }
            
            # Initialize lists
            animal_list <- list(
              "animal_1" = list(name = "", time = "", position = 0),
              "animal_2" = list(name = "", time = "", position = 0),
              "animal_3" = list(name = "", time = "", position = 0),
              "animal_4" = list(name = "", time = "", position = 0),
              "data_temp" = list(elapsed_seconds = 0, current_row = 0)
            )
            
            cage_position_probability <- list(
              c(1, 0, 0, 0), c(2, 0, 0, 0), c(3, 0, 0, 0), c(4, 0, 0, 0),
              c(5, 0, 0, 0), c(6, 0, 0, 0), c(7, 0, 0, 0), c(8, 0, 0, 0)
            )
            
            animal_position_probability <- tibble(
              AnimalID = rep(animal_ids, each = 8),
              Position = rep(1:8, length.out = 32),
              Seconds = 0,
              SumPercentage = 0,
              Prob = 0
            )
            
            # Calculate probabilities
            animal_list <- initialize_animal_positions(animal_ids, data_system_halfhour, animal_list)
            
            initial_time <- animal_list[[1]][[2]]
            current_row <- 5
            elapsed_seconds <- 0
            total_rows <- nrow(data_system_halfhour) + 1
            
            while (current_row != total_rows && current_row < total_rows) {
              
              previous_animal_positions <- animal_list
              
              animal_list <- update_animal_list(animal_ids, animal_list, data_system_halfhour, 
                                                initial_time, current_row)
              
              elapsed_seconds <- animal_list[["data_temp"]][["elapsed_seconds"]]
              
              if (system_complete) {
                cage_position_probability <- update_cage_position_probability(
                  previous_animal_positions, animal_list, cage_position_probability, elapsed_seconds
                )
              }
              
              animal_position_probability <- update_animal_position_probability(
                previous_animal_positions, animal_list, animal_position_probability, elapsed_seconds
              )
              
              current_row <- animal_list[["data_temp"]][["current_row"]]
              initial_time <- animal_list[[1]][[2]]
            }
            
            # Calculate final probabilities
            for (i in 1:8) {
              if (cage_position_probability[[i]][[2]] > 0) {
                cage_position_probability[[i]][[4]] <- cage_position_probability[[i]][[3]] / 
                                                        cage_position_probability[[i]][[2]]
              } else {
                cage_position_probability[[i]][[4]] <- 0
              }
            }
            
            animal_position_probability <- animal_position_probability %>%
              mutate(Prob = ifelse(Seconds > 0, SumPercentage / Seconds, 0))
            
            # Calculate Shannon entropy
            if (system_complete) {
              cage_prob_vec <- sapply(cage_position_probability, function(x) x[4])
              cage_prob_vec <- cage_prob_vec[!is.na(cage_prob_vec) & !is.nan(cage_prob_vec)]
              
              if (length(cage_prob_vec) > 0 && sum(cage_prob_vec) > 0) {
                cage_prob_vec <- cage_prob_vec / sum(cage_prob_vec)
                cage_shannon_entropy <- calc_shannon_entropy(cage_prob_vec)
              } else {
                cage_shannon_entropy <- NA
              }
              
              cagePosEntropy_halfhour <- cagePosEntropy_halfhour %>%
                add_row(Batch = batch, Sex = sex, System = system_id,
                        CageChange = cageChange, HalfHour = halfhour_period,
                        CageEntropy = cage_shannon_entropy)
            }
            
            # Animal entropy
            for (animal in animal_ids) {
              if (is.na(animal)) next
              
              animal_prob_vec <- animal_position_probability %>%
                filter(AnimalID == animal) %>%
                pull(Prob)
              
              animal_prob_vec <- animal_prob_vec[!is.na(animal_prob_vec) & !is.nan(animal_prob_vec)]
              
              if (length(animal_prob_vec) > 0 && sum(animal_prob_vec) > 0) {
                animal_prob_vec <- animal_prob_vec / sum(animal_prob_vec)
                animal_shannon_entropy <- calc_shannon_entropy(animal_prob_vec)
              } else {
                animal_shannon_entropy <- NA
              }
              
              animalPosEntropy_halfhour <- animalPosEntropy_halfhour %>%
                add_row(Batch = batch, Sex = sex, System = system_id,
                        CageChange = cageChange, HalfHour = halfhour_period,
                        AnimalID = animal, animalEntropy = animal_shannon_entropy)
            }
          }
        }
      }
      
      # ===================================================================
      # SAVE BATCH/CAGECHANGE RESULTS
      # ===================================================================
      if (save_tables == TRUE) {
        message("Saving results for ", batch, " ", cageChange)
        
        write.csv(cagePosProb, 
                  file = paste0(saving_directory, tables_dir_phase, "/", 
                               batch, "_", cageChange, "_cagePosProb.csv"), 
                  row.names = FALSE)
        
        write.csv(cagePosEntropy, 
                  file = paste0(saving_directory, tables_dir_phase, "/", 
                               batch, "_", cageChange, "_cagePosEntropy.csv"), 
                  row.names = FALSE)
        
        write.csv(animalPosEntropy, 
                  file = paste0(saving_directory, tables_dir_phase, "/", 
                               batch, "_", cageChange, "_animalPosEntropy.csv"), 
                  row.names = FALSE)
        
        if (analyze_by_halfhour) {
          write.csv(cagePosEntropy_halfhour, 
                    file = paste0(saving_directory, tables_dir_halfhour, "/", 
                                 batch, "_", cageChange, "_cagePosEntropy_halfhour.csv"), 
                    row.names = FALSE)
          
          write.csv(animalPosEntropy_halfhour, 
                    file = paste0(saving_directory, tables_dir_halfhour, "/", 
                                 batch, "_", cageChange, "_animalPosEntropy_halfhour.csv"), 
                    row.names = FALSE)
        }
      }
      
      # Accumulate in global tibbles
      all_cagePosProb <- bind_rows(all_cagePosProb, cagePosProb)
      all_cagePosEntropy <- bind_rows(all_cagePosEntropy, cagePosEntropy)
      all_animalPosEntropy <- bind_rows(all_animalPosEntropy, animalPosEntropy)
      
      if (analyze_by_halfhour) {
        all_cagePosEntropy_halfhour <- bind_rows(all_cagePosEntropy_halfhour, cagePosEntropy_halfhour)
        all_animalPosEntropy_halfhour <- bind_rows(all_animalPosEntropy_halfhour, animalPosEntropy_halfhour)
      }
      
      message("FINISHED ", batch, " ", cageChange)
      message("")
    }
  }
  
  # ===================================================================
  # SAVE COMBINED RESULTS
  # ===================================================================
  if (save_tables == TRUE) {
    message("=======================================================================")
    message("SAVING COMBINED RESULTS")
    message("=======================================================================")
    
    combined_tables_dir <- paste0(tables_base, "/combined")
    dir.create(paste0(saving_directory, combined_tables_dir), recursive = TRUE, showWarnings = FALSE)
    
    write.csv(all_cagePosProb, 
              file = paste0(saving_directory, combined_tables_dir, 
                           "/all_batches_all_cageChanges_cagePosProb.csv"), 
              row.names = FALSE)
    message("   ✓ Saved: all_batches_all_cageChanges_cagePosProb.csv")
    
    write.csv(all_cagePosEntropy, 
              file = paste0(saving_directory, combined_tables_dir, 
                           "/all_batches_all_cageChanges_cagePosEntropy.csv"), 
              row.names = FALSE)
    message("   ✓ Saved: all_batches_all_cageChanges_cagePosEntropy.csv")
    
    write.csv(all_animalPosEntropy, 
              file = paste0(saving_directory, combined_tables_dir, 
                           "/all_batches_all_cageChanges_animalPosEntropy.csv"), 
              row.names = FALSE)
    message("   ✓ Saved: all_batches_all_cageChanges_animalPosEntropy.csv")
    
    if (analyze_by_halfhour) {
      write.csv(all_cagePosEntropy_halfhour, 
                file = paste0(saving_directory, combined_tables_dir, 
                             "/all_batches_all_cageChanges_cagePosEntropy_halfhour.csv"), 
                row.names = FALSE)
      message("   ✓ Saved: all_batches_all_cageChanges_cagePosEntropy_halfhour.csv")
      
      write.csv(all_animalPosEntropy_halfhour, 
                file = paste0(saving_directory, combined_tables_dir, 
                             "/all_batches_all_cageChanges_animalPosEntropy_halfhour.csv"), 
                row.names = FALSE)
      message("   ✓ Saved: all_batches_all_cageChanges_animalPosEntropy_halfhour.csv")
    }
    
    message("      Location: ", combined_tables_dir)
  }
  
  message("=======================================================================")
  message("DATA PROCESSING COMPLETE")
  message("=======================================================================")
}

# ===================================================================
# CONSECUTIVE ENTROPY PROCESSING (ALWAYS RUNS)
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
  
  # Apply CC4 filter if enabled
  if (filter_cc4_late_phases) {
    entropy_list <- entropy_list %>%
      filter(!(CageChange == "CC4" & ConsecActive > cc4_max_active_phase)) %>%
      filter(!(CageChange == "CC4" & ConsecInactive > cc4_max_inactive_phase))
    
    message("✓ Filtered CC4: keeping only phases ≤ A", cc4_max_active_phase, " and I", cc4_max_inactive_phase)
  }
  
  if (name == "consec_animal_entropy") {
    entropy_list <- entropy_list %>%
      mutate(Group = ifelse(AnimalID %in% sus_animals, "sus", 
                           ifelse(AnimalID %in% con_animals, "con", "res")))
  }
  
  for (change in c("CC1", "CC2", "CC3", "CC4")) {
    if (change != "CC1") {
      entropy_list <- entropy_list %>%
        mutate(ConsecActive = ifelse(CageChange == change & Phase == "active", 
                                     ConsecActive + max_consecAct, ConsecActive)) %>%
        mutate(ConsecInactive = ifelse(CageChange == change & Phase == "inactive", 
                                       ConsecInactive + max_consecInact, ConsecInactive))
    }
    
    max_consecAct <- entropy_list %>%
      filter(CageChange == change) %>%
      pull(ConsecActive) %>%
      unique() %>%
      max()
    
    max_consecInact <- entropy_list %>%
      filter(CageChange == change) %>%
      pull(ConsecInactive) %>%
      unique() %>%
      max()
  }
  
  if (name == "consec_animal_entropy") {
    entropy_list <- entropy_list[c("CageChange", "Batch", "System", "AnimalID",
                                   "Sex", "Group", "Phase", "ConsecActive",
                                   "ConsecInactive", "animalEntropy")]
  } else {
    entropy_list <- entropy_list[c("CageChange", "Batch", "System", "Sex",
                                   "Phase", "ConsecActive", "ConsecInactive",
                                   "CageEntropy")]
  }
  
  if (name == "consec_cage_entropy") {
    consec_cage_entropy <- entropy_list
  } else {
    consec_animal_entropy <- entropy_list
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
  
  for (i in seq_along(entropy_halfhour_names)) {
    name <- entropy_halfhour_names[i]
    entropy_list <- if (name == "consec_cage_entropy_halfhour") {
      consec_cage_entropy_halfhour
    } else {
      consec_animal_entropy_halfhour
    }
    
    if (name == "consec_animal_entropy_halfhour") {
      entropy_list <- entropy_list %>%
        mutate(Group = ifelse(AnimalID %in% sus_animals, "sus", 
                             ifelse(AnimalID %in% con_animals, "con", "res")))
    }
    
    entropy_list <- entropy_list %>%
      mutate(ConsecHalfHour = HalfHour)
    
    max_halfhour <- 0
    
    for (change in c("CC1", "CC2", "CC3", "CC4")) {
      if (change != "CC1") {
        entropy_list <- entropy_list %>%
          mutate(ConsecHalfHour = ifelse(CageChange == change, 
                                         HalfHour + max_halfhour, 
                                         ConsecHalfHour))
      }
      
      max_halfhour <- entropy_list %>%
        filter(CageChange == change) %>%
        pull(ConsecHalfHour) %>%
        max(na.rm = TRUE)
    }
    
    if (name == "consec_animal_entropy_halfhour") {
      entropy_list <- entropy_list[c("CageChange", "Batch", "System", "AnimalID",
                                     "Sex", "Group", "HalfHour", "ConsecHalfHour",
                                     "animalEntropy")]
    } else {
      entropy_list <- entropy_list[c("CageChange", "Batch", "System", "Sex",
                                     "HalfHour", "ConsecHalfHour", "CageEntropy")]
    }
    
    if (name == "consec_cage_entropy_halfhour") {
      consec_cage_entropy_halfhour <- entropy_list
    } else {
      consec_animal_entropy_halfhour <- entropy_list
    }
  }
}

if (analyze_by_halfhour && filter_cc4_late_phases) {
  message("=======================================================================")
  message("FILTERING CC4 HALF-HOURS BASED ON PHASE INFORMATION")
  message("=======================================================================")
  
  # Check if we have preprocessed data (only available if we ran processing, not if we loaded)
  if (!exists("data_preprocessed")) {
    message("   Warning: Cannot filter CC4 half-hours by phase - preprocessed data not available")
    message("   Run with load_existing_data = FALSE to enable this feature, or use cc4_halfhour_cutoff instead")
  } else {
    
    # Create lookup from preprocessed data
    cc4_valid_halfhours <- data_preprocessed %>%
      filter(CageChange == "CC4") %>%
      filter(ConsecActive <= cc4_max_active_phase | ConsecInactive <= cc4_max_inactive_phase) %>%
      select(Batch, System, CageChange, HalfHour = HalfHoursElapsed) %>%
      distinct() %>%
      mutate(keep_halfhour = TRUE)
    
    # Filter animal entropy
    consec_animal_entropy_halfhour <- consec_animal_entropy_halfhour %>%
      left_join(cc4_valid_halfhours,
                by = c("Batch", "System", "CageChange", "HalfHour")) %>%
      filter(CageChange != "CC4" | !is.na(keep_halfhour)) %>%
      select(-keep_halfhour)
    
    # Filter cage entropy  
    consec_cage_entropy_halfhour <- consec_cage_entropy_halfhour %>%
      left_join(cc4_valid_halfhours,
                by = c("Batch", "System", "CageChange", "HalfHour")) %>%
      filter(CageChange != "CC4" | !is.na(keep_halfhour)) %>%
      select(-keep_halfhour)
    
    message("✓ Filtered CC4 half-hours based on phase information")
    message("   Remaining half-hour observations: ", nrow(consec_animal_entropy_halfhour))
  }
}

# ===================================================================
# SAVE CONSECUTIVE ENTROPY TABLES
# ===================================================================
if (save_tables == TRUE) {
  message("Saving consecutive entropy tables...")
  
  combined_tables_dir <- paste0(tables_base, "/combined")
  
  write.csv(consec_cage_entropy, 
            file = paste0(saving_directory, combined_tables_dir, 
                         "/all_batches_all_cageChanges_consec_cage_entropy.csv"), 
            row.names = FALSE)
  message("   ✓ Saved: all_batches_all_cageChanges_consec_cage_entropy.csv")
  
  write.csv(consec_animal_entropy, 
            file = paste0(saving_directory, combined_tables_dir, 
                         "/all_batches_all_cageChanges_consec_animal_entropy.csv"), 
            row.names = FALSE)
  message("   ✓ Saved: all_batches_all_cageChanges_consec_animal_entropy.csv")
  
  if (analyze_by_halfhour) {
    write.csv(consec_cage_entropy_halfhour, 
              file = paste0(saving_directory, combined_tables_dir, 
                           "/all_batches_all_cageChanges_consec_cage_entropy_halfhour.csv"), 
              row.names = FALSE)
    message("   ✓ Saved: all_batches_all_cageChanges_consec_cage_entropy_halfhour.csv")
    
    write.csv(consec_animal_entropy_halfhour, 
              file = paste0(saving_directory, combined_tables_dir, 
                           "/all_batches_all_cageChanges_consec_animal_entropy_halfhour.csv"), 
              row.names = FALSE)
    message("   ✓ Saved: all_batches_all_cageChanges_consec_animal_entropy_halfhour.csv")
  }
  
  message("      Location: ", combined_tables_dir)
}

# ===================================================================
# EXPLORATION METRICS (OPTION 2)
# ===================================================================
message("=======================================================================")
message("PROCESSING EXPLORATION METRICS (OPTION 2)")
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

# Calculate summary statistics
exploration_summary_phase <- consec_animal_entropy %>%
  group_by(Group, Sex, CageChange, Phase) %>%
  summarise(
    n_observations = n(),
    n_explored = sum(explored),
    pct_exploration = round(100 * n_explored / n_observations, 2),
    mean_entropy_overall = round(mean(exploration_entropy, na.rm = TRUE), 3),
    mean_entropy_when_exploring = round(mean(animalEntropy, na.rm = TRUE), 3),
    sd_entropy_when_exploring = round(sd(animalEntropy, na.rm = TRUE), 3),
    .groups = 'drop'
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
    .groups = 'drop'
  ) %>%
  arrange(desc(pct_exploration))

# Save enhanced datasets
if (save_tables == TRUE) {
  message("Saving enhanced exploration datasets...")
  
  combined_tables_dir <- paste0(tables_base, "/combined")
  
  write.csv(consec_animal_entropy, 
            file = paste0(saving_directory, combined_tables_dir, 
                         "/consec_animal_entropy_with_exploration_metrics.csv"), 
            row.names = FALSE)
  message("   ✓ Saved: consec_animal_entropy_with_exploration_metrics.csv")
  
  if (analyze_by_halfhour) {
    write.csv(consec_animal_entropy_halfhour, 
              file = paste0(saving_directory, combined_tables_dir, 
                           "/consec_animal_entropy_halfhour_with_exploration_metrics.csv"), 
              row.names = FALSE)
    message("   ✓ Saved: consec_animal_entropy_halfhour_with_exploration_metrics.csv")
  }
  
  write.csv(exploration_summary_phase, 
            file = paste0(saving_directory, combined_tables_dir, 
                         "/exploration_summary_by_phase.csv"), 
            row.names = FALSE)
  message("   ✓ Saved: exploration_summary_by_phase.csv")
  
  if (analyze_by_halfhour) {
    write.csv(exploration_summary_halfhour, 
              file = paste0(saving_directory, combined_tables_dir, 
                           "/exploration_summary_by_halfhour.csv"), 
              row.names = FALSE)
    message("   ✓ Saved: exploration_summary_by_halfhour.csv")
  }
  
  write.csv(animal_exploration_profile, 
            file = paste0(saving_directory, combined_tables_dir, 
                         "/animal_exploration_profiles.csv"), 
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
  
  combined_plots_dir <- paste0(plots_base, "/combined")
  dir.create(paste0(saving_directory, combined_plots_dir), recursive = TRUE, showWarnings = FALSE)
  
  # Phase-based animal entropy plots
  message("Generating phase-based animal entropy plots...")
  
  columns_to_group <- list(
    c("Sex", "AnimalID", "Group"),
    c("Sex", "AnimalID", "CageChange", "Group"),
    c("Sex", "AnimalID", "Phase", "Group"),
    c("Sex", "AnimalID", "Phase", "Group")
  )
  x_axis <- list("Group", "CageChange", "Group", "Group")
  phasecount <- 1
  plots <- list()
  
  for (i in seq_along(columns_to_group)) {
    group_vars <- columns_to_group[[i]]
    x_var <- x_axis[[i]]
    data <- consec_animal_entropy
    
    if ("Phase" %in% group_vars) {
      data <- filter(data, Phase == ifelse(phasecount == 1, "active", "inactive"))
      phasecount <- phasecount + 1
    }
    
    result <- data %>%
      group_by(across(all_of(group_vars))) %>%
      summarise(Mean_Entropy = mean(animalEntropy, na.rm = TRUE), .groups = 'drop')
    
    p <- ggplot(data = result, aes(x = !!sym(x_var),
                                   y = Mean_Entropy,
                                   color = Group,
                                   group = Group)) +
      geom_jitter(aes(fill = Group),
                 size = 4,
                 alpha = 0.7,
                 shape = 16,
                 position = position_dodge(width = ifelse("CageChange" %in% group_vars, 0.75, 0))) +
      labs(title = "Animal Entropy",
           subtitle = paste("Grouped by:", paste(group_vars, collapse = ", ")),
           x = x_var,
           y = "Mean Entropy") +
      scale_color_manual(values = c("sus" = "#E63946",
                                   "res" = "grey60",
                                   "con" = "#457B9D")) +
      facet_grid(ifelse("Phase" %in% group_vars, "Phase~Sex", "~Sex")) +
      stat_summary(fun.min = function(z) {quantile(z, 0.25)},
                  fun.max = function(z) {quantile(z, 0.75)},
                  fun = median,
                  color = "black",
                  size = 0.8,
                  shape = 16,
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
      ggsave(filename = paste0(saving_directory, combined_plots_dir, 
                              "/animal_entropy_phase_plot_", i, ".svg"),
             plot = p,
             width = 12,
             height = 8,
             dpi = 300)
    }
  }
  
  message("   ✓ Saved ", length(plots), " phase-based animal entropy plots")
  
  # Cage entropy plots
  message("Generating phase-based cage entropy plots...")
  
  cagePosEntropy_act <- consec_cage_entropy %>%
    filter(Phase == "active")
  
  cage_ent_plot_act <- ggplot(data = cagePosEntropy_act, 
                              aes(x = ConsecActive, y = CageEntropy, color = System)) +
    geom_jitter(aes(fill = System), size = 4, alpha = 0.7, width = 0.2, shape = 16) +
    scale_y_continuous("Cage Entropy") +
    scale_x_continuous("Active Phase Number") +
    labs(title = "Cage Entropy - Active Phases") +
    stat_summary(fun.min = function(z) {quantile(z, 0.25)},
                fun.max = function(z) {quantile(z, 0.75)},
                fun = median,
                color = "black",
                size = 0.8,
                shape = 16) +
    facet_grid(Sex ~ CageChange) +
    theme_minimal() +
    theme(title = element_text(size = 20),
          legend.key.size = unit(3, "lines"),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 22),
          strip.text = element_text(size = 20))
  
  cagePosEntropy_inact <- consec_cage_entropy %>%
    filter(Phase == "inactive")
  
  cage_ent_plot_inact <- ggplot(data = cagePosEntropy_inact, 
                                aes(x = ConsecInactive, y = CageEntropy, color = System)) +
    geom_jitter(aes(fill = System), size = 4, alpha = 0.7, width = 0.2, shape = 16) +
    scale_y_continuous("Cage Entropy") +
    scale_x_continuous("Inactive Phase Number") +
    labs(title = "Cage Entropy - Inactive Phases") +
    stat_summary(fun.min = function(z) {quantile(z, 0.25)},
                fun.max = function(z) {quantile(z, 0.75)},
                fun = median,
                color = "black",
                size = 0.8,
                shape = 16) +
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
    ggsave(filename = paste0(saving_directory, combined_plots_dir, 
                            "/cage_entropy_active.svg"),
           plot = cage_ent_plot_act,
           width = 14,
           height = 10,
           dpi = 300)
    
    ggsave(filename = paste0(saving_directory, combined_plots_dir, 
                            "/cage_entropy_inactive.svg"),
           plot = cage_ent_plot_inact,
           width = 14,
           height = 10,
           dpi = 300)
    
    message("   ✓ Saved cage entropy plots (phase-based)")
  }
  
  # Half-hour plots
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
      scale_color_manual(values = c("sus" = "#E63946",
                                   "res" = "grey60",
                                   "con" = "#457B9D")) +
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
        geom_text(data = boundary_labels, 
                 aes(x = x, y = Inf, label = label),
                 inherit.aes = FALSE,
                 vjust = 1.5, size = 3, angle = 90, color = "black")
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
        geom_text(data = boundary_labels, 
                 aes(x = x, y = Inf, label = label),
                 inherit.aes = FALSE,
                 vjust = 1.5, size = 3, angle = 90, color = "black")
    }
    
    animal_halfhour_plot_bycage <- ggplot(data = consec_animal_entropy_halfhour,
                                          aes(x = HalfHour, y = animalEntropy, color = Group)) +
      geom_point(size = 2, alpha = 0.5, na.rm = TRUE) +
      geom_smooth(aes(group = Group), method = "loess", se = TRUE, span = 0.5, na.rm = TRUE) +
      labs(title = "Animal Entropy Over Time (Half-Hour Intervals)",
           subtitle = "Separated by cage change",
           x = "Half-Hour Period (within cage change)",
           y = "Entropy") +
      scale_color_manual(values = c("sus" = "#E63946",
                                   "res" = "grey60",
                                   "con" = "#457B9D")) +
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
      ggsave(filename = paste0(saving_directory, combined_plots_dir, 
                              "/animal_entropy_halfhour_consecutive.svg"),
             plot = animal_halfhour_plot_consec,
             width = 14,
             height = 10,
             dpi = 300)
      
      ggsave(filename = paste0(saving_directory, combined_plots_dir, 
                              "/cage_entropy_halfhour_consecutive.svg"),
             plot = cage_halfhour_plot_consec,
             width = 14,
             height = 10,
             dpi = 300)
      
      ggsave(filename = paste0(saving_directory, combined_plots_dir, 
                              "/animal_entropy_halfhour_by_cagechange.svg"),
             plot = animal_halfhour_plot_bycage,
             width = 16,
             height = 10,
             dpi = 300)
      
      ggsave(filename = paste0(saving_directory, combined_plots_dir, 
                              "/cage_entropy_halfhour_by_cagechange.svg"),
             plot = cage_halfhour_plot_bycage,
             width = 16,
             height = 10,
             dpi = 300)
      
      message("   ✓ Saved 4 half-hour entropy plots")
    }
    
    if (show_plots) {
      print(animal_halfhour_plot_consec)
      print(cage_halfhour_plot_consec)
      print(animal_halfhour_plot_bycage)
      print(cage_halfhour_plot_bycage)
    }
  }
  
  if (show_plots) {
    for (p in plots) {
      print(p)
    }
    print(cage_ent_plot_act)
    print(cage_ent_plot_inact)
  }
  
  message("      Location: ", combined_plots_dir)
}

message("=======================================================================")
message("ENTROPY ANALYSIS COMPLETE!")
message("=======================================================================")
