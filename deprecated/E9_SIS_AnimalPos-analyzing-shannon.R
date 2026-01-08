#' @title Analysis of Shannon Entropy in Animal Positioning
#'
#' @description
#' This script conducts a comprehensive analysis of Shannon entropy to
#' evaluate the distribution and diversity of animal positions within a
#' specified environment, analyzing both by phase and by half-hour periods.
#'
#' @details
#' This script is an integral component of the MMMSociability project,
#' aimed at studying sociability patterns among animals. Results are now
#' organized by batch and cage change for better data management.
#'
#' @date October 2025
#' @authors Tobias Pohl, Anja Magister

# Load required packages using pacman
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(readr, dplyr, lubridate, tibble, purrr, ggplot2, reshape2,
               scales, stringr, gridExtra)

# ===================================================================
# CUSTOMIZABLE VARIABLES
# ===================================================================
show_plots <- FALSE
save_plots <- TRUE
save_tables <- TRUE
exclude_homecage <- TRUE  # Set to TRUE to exclude home cage positions
analyze_by_halfhour <- TRUE  # Set to TRUE to also analyze by half-hour periods

# ===================================================================
# PATHS
# ===================================================================
# Choose working directory (comment/uncomment as needed)
#working_directory <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability"
#saving_directory <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability"

# Alternative local directory
working_directory <- "D:/MMMSociability"
saving_directory <- "D:/MMMSociability"

# Define base directories for entropy analysis
#plots_base <- "/plots/Entropy"
#tables_base <- "/tables/Entropy"
plots_base <- paste0("/plots/", ifelse(exclude_homecage, "noHomeCage", "withHomeCage"))
tables_base <- paste0("/tables/", ifelse(exclude_homecage, "noHomeCage", "withHomeCage"))

# ===================================================================
# LOAD CUSTOM FUNCTIONS
# ===================================================================
source("C:/Users/Tobias Pohl/Documents/GitHub/MMMSociability/Functions/E9_SIS_AnimalPos-functions.R")

# Load lists of susceptible and control animals for analysis
sus_animals <- readLines(paste0(working_directory, "/raw_data/sus_animals.csv"))
con_animals <- readLines(paste0(working_directory, "/raw_data/con_animals.csv"))

# ===================================================================
# INITIALIZE GLOBAL RESULT TIBBLES (PHASE-BASED)
# ===================================================================
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

# ===================================================================
# INITIALIZE GLOBAL RESULT TIBBLES (HALF-HOUR)
# ===================================================================
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
    
    # ===================================================================
    # CREATE FOLDER STRUCTURE FOR CURRENT BATCH AND CAGE CHANGE
    # ===================================================================
    current_plots_dir <- paste0(plots_base, "/", batch, "/", cageChange)
    current_tables_dir <- paste0(tables_base, "/", batch, "/", cageChange)
    
    # Create subdirectories
    plots_dir_phase <- paste0(current_plots_dir, "/by_phase")
    tables_dir_phase <- paste0(current_tables_dir, "/by_phase")
    tables_dir_metadata <- paste0(current_tables_dir, "/metadata")
    
    # Create all directories
    dir.create(paste0(saving_directory, plots_dir_phase), recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0(saving_directory, tables_dir_phase), recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0(saving_directory, tables_dir_metadata), recursive = TRUE, showWarnings = FALSE)
    
    # Create half-hour directories if needed
    if (analyze_by_halfhour) {
      plots_dir_halfhour <- paste0(current_plots_dir, "/by_halfhour")
      tables_dir_halfhour <- paste0(current_tables_dir, "/by_halfhour")
      
      dir.create(paste0(saving_directory, plots_dir_halfhour), recursive = TRUE, showWarnings = FALSE)
      dir.create(paste0(saving_directory, tables_dir_halfhour), recursive = TRUE, showWarnings = FALSE)
    }
    
    message("✓ Created directory structure for ", batch, " ", cageChange)
    
    # ===================================================================
    # INITIALIZE LOCAL RESULT TIBBLES FOR THIS BATCH/CAGECHANGE
    # ===================================================================
    # Phase-based tibbles
    cagePosProb <- tibble(
      Batch = character(),
      System = character(),
      CageChange = character(),
      Phase = character(),
      Position = numeric(),
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
    
    # Half-hour tibbles
    if (analyze_by_halfhour) {
      cagePosEntropy_halfhour <- tibble(
        Batch = character(),
        Sex = character(),
        System = character(),
        CageChange = character(),
        HalfHour = numeric(),
        CageEntropy = numeric()
      )
      
      animalPosEntropy_halfhour <- tibble(
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
    # READ PREPROCESSED DATA
    # ===================================================================
    filename <- paste0("E9_SIS_", batch, "_", cageChange, "_AnimalPos")
    csvFilePath <- paste0(working_directory, "/preprocessed_data_test/", filename, "_preprocessed.csv")
    
    # Check if file exists
    if (!file.exists(csvFilePath)) {
      message("✗ File not found: ", csvFilePath)
      message("   Skipping ", batch, " ", cageChange)
      next
    }
    
    data_preprocessed <- read_delim(csvFilePath, delim = ",", show_col_types = FALSE) %>%
      as_tibble()
    
    message("✓ Loaded preprocessed data: ", nrow(data_preprocessed), " rows")
    
    # ===================================================================
    # DEFINE SYSTEMS AND PHASES
    # ===================================================================
    unique_systems <- unique(data_preprocessed$System)
    unique_systems <- str_sort(unique_systems)
    
    phases <- c("Active", "Inactive")
    
    # Identify phase numbers
    active_phase_count <- unique(data_preprocessed$ConsecActive)
    inactive_phase_count <- unique(data_preprocessed$ConsecInactive)
    
    # Remove zeros and filter incomplete phases
    active_phase_count <- active_phase_count[!active_phase_count %in% 0]
    inactive_phase_count <- inactive_phase_count[!inactive_phase_count %in% c(0, 1, max(inactive_phase_count))]
    
    # Remove phases > 2 from CC4
    if (cageChange == "CC4") {
      active_phase_count <- active_phase_count[active_phase_count <= 2]
      inactive_phase_count <- inactive_phase_count[inactive_phase_count <= 2]
    }
    
    # Sort to ensure correct order
    active_phase_count <- sort(active_phase_count)
    inactive_phase_count <- sort(inactive_phase_count)
    
    # Get unique half-hour periods
    if (analyze_by_halfhour) {
      halfhour_periods <- sort(unique(data_preprocessed$HalfHoursElapsed))
      message("✓ Found ", length(halfhour_periods), " half-hour periods")
    }
    
    # Define unique animals
    unique_animals <- unique(data_preprocessed$AnimalID)
    
    message("✓ Found ", length(unique_systems), " systems, ", 
            length(active_phase_count), " active phases, ",
            length(inactive_phase_count), " inactive phases")
    
    # ===================================================================
    # ANALYSIS BY PHASE: ITERATE THROUGH SYSTEMS AND PHASES
    # ===================================================================
    message("Starting PHASE-BASED entropy analysis...")
    
    for (system_id in unique_systems) {
      
      message("   Processing system: ", system_id)
      
      # Filter data for current system
      system_data <- data_preprocessed %>%
        filter(System == system_id) %>%
        as_tibble()
      
      # Define animal IDs for current system
      animal_ids <- unique(system_data$AnimalID)
      system_complete <- ifelse(length(animal_ids) < 4, FALSE, TRUE)
      
      # Fill vector with NAs if system is incomplete
      while (length(animal_ids) < 4) {
        animal_ids <- append(animal_ids, NA)
      }
      
      # Loop through phases
      for (phase in phases) {
        
        # Select phase numbers
        phase_numbers <- if (phase == "Active") active_phase_count else inactive_phase_count
        
        # Loop through phase numbers
        for (phase_number in phase_numbers) {
          
          message("      ", phase, " phase ", phase_number)
          
          # Filter data for current phase
          data_system_phase <- system_data %>%
            filter(ConsecActive == ifelse(phase == "Active", phase_number, 0)) %>%
            filter(ConsecInactive == ifelse(phase == "Inactive", phase_number, 0)) %>%
            as_tibble()
          
          # Skip if no data
          if (nrow(data_system_phase) == 0) {
            message("         No data found. Skipping...")
            next
          }
          
          # ===================================================================
          # INITIALIZE LISTS
          # ===================================================================
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
          
          # ===================================================================
          # ANALYSIS: CALCULATE PROBABILITIES
          # ===================================================================
          animal_list <- initialize_animal_positions(animal_ids, data_system_phase, animal_list)
          
          initial_time <- animal_list[[1]][[2]]
          current_row <- 5
          elapsed_seconds <- 0
          total_rows <- nrow(data_system_phase) + 1
          
          # Process each row
          while (current_row != total_rows && current_row < total_rows) {
            
            previous_animal_positions <- animal_list
            
            animal_list <- update_animal_list(animal_ids, animal_list, data_system_phase, 
                                              initial_time, current_row)
            
            elapsed_seconds <- animal_list[["data_temp"]][["elapsed_seconds"]]
            
            # Update probabilities
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
          
          # ===================================================================
          # CALCULATE FINAL PROBABILITIES (PHASE)
          # ===================================================================
          # Calculate the probability for each cage position
          for (i in 1:8) {
            # Handle division by zero
            if (cage_position_probability[[i]][[2]] > 0) {
              cage_position_probability[[i]][[4]] <- cage_position_probability[[i]][[3]] / 
                                                      cage_position_probability[[i]][[2]]
            } else {
              cage_position_probability[[i]][[4]] <- 0
            }
          }
          
          # Calculate the probability for each animal position, handling division by zero
          animal_position_probability <- animal_position_probability %>%
            mutate(Prob = ifelse(Seconds > 0, SumPercentage / Seconds, 0))
          
          # Print the results
          message("Results: ")
          message("Animal Position Probability")
          print(animal_position_probability)
          
          # ===================================================================
          # STORE RESULTS IN TIBBLES (PHASE)
          # ===================================================================
          if (system_complete) {
            message("Entering data from this phase into the total result tibble")
            # Enter information into the result tibble for each phase and system
            for (i in 1:8) { # For each position in the cage
              p <- paste0(substr(phase, 1, 1), phase_number) # The current phase
              cagePosProb <- cagePosProb %>%
                add_row(Batch = batch,
                        System = system_id,
                        CageChange = cageChange,
                        Phase = p,
                        Position = i,
                        Probability = cage_position_probability[[i]][[4]])
            }
          }
          
          # ===================================================================
          # SHANNON ENTROPY (PHASE)
          # ===================================================================
          if (system_complete) {
            # Cage entropy - extract DIRECTLY from cage_position_probability list
            cage_prob_vec <- sapply(cage_position_probability, function(x) x[4])
            
            # Remove NA/NaN values
            cage_prob_vec <- cage_prob_vec[!is.na(cage_prob_vec) & !is.nan(cage_prob_vec)]
            
            # Only calculate if we have valid probabilities
            if (length(cage_prob_vec) > 0 && sum(cage_prob_vec) > 0) {
              # Normalize probabilities
              cage_prob_vec <- cage_prob_vec / sum(cage_prob_vec)
              cage_shannon_entropy <- calc_shannon_entropy(cage_prob_vec)
            } else {
              cage_shannon_entropy <- NA
            }
            
            # Add new row with Shannon entropy into result tibble
            cagePosEntropy <- cagePosEntropy %>%
              add_row(Batch = batch,
                      Sex = sex,
                      System = system_id,
                      CageChange = cageChange,
                      Phase = paste0(substr(phase, 1, 1), phase_number),
                      CageEntropy = cage_shannon_entropy)
          }
          
          ## Animals
          # Calculate probabilities and entropies for each animal
          for (animal in animal_ids) {
            
            # Skip if animal is not tracked (incomplete system)
            if (is.na(animal)) {next}
            
            # Extract probability vector for the current animal
            animal_prob_vec <- animal_position_probability %>%
              filter(AnimalID == animal) %>%
              pull(Prob)
            
            # Remove NA/NaN values
            animal_prob_vec <- animal_prob_vec[!is.na(animal_prob_vec) & !is.nan(animal_prob_vec)]
            
            # Only calculate if we have valid probabilities
            if (length(animal_prob_vec) > 0 && sum(animal_prob_vec) > 0) {
              # Normalize probabilities
              animal_prob_vec <- animal_prob_vec / sum(animal_prob_vec)
              animal_shannon_entropy <- calc_shannon_entropy(animal_prob_vec)
            } else {
              animal_shannon_entropy <- NA
            }
            
            # Add the calculated entropy to the result tibble
            animalPosEntropy <- animalPosEntropy %>%
              add_row(Batch = batch,
                      Sex = sex,
                      System = system_id,
                      CageChange = cageChange,
                      Phase = paste0(substr(phase, 1, 1), phase_number),
                      AnimalID = animal,
                      animalEntropy = animal_shannon_entropy)
          }
        }
      }
    }
    
    # ===================================================================
    # ANALYSIS BY HALF-HOUR: ITERATE THROUGH SYSTEMS AND PERIODS
    # ===================================================================
    if (analyze_by_halfhour) {
      message("Starting HALF-HOUR entropy analysis...")
      
      for (system_id in unique_systems) {
        
        message("   Processing system: ", system_id)
        
        # Filter data for current system
        system_data <- data_preprocessed %>%
          filter(System == system_id) %>%
          as_tibble()
        
        # Define animal IDs for current system
        animal_ids <- unique(system_data$AnimalID)
        system_complete <- ifelse(length(animal_ids) < 4, FALSE, TRUE)
        
        # Fill vector with NAs if system is incomplete
        while (length(animal_ids) < 4) {
          animal_ids <- append(animal_ids, NA)
        }
        
        # Loop through half-hour periods
        for (halfhour_period in halfhour_periods) {
          
          message("      Half-hour period: ", halfhour_period)
          
          # Filter data for current half-hour period
          data_system_halfhour <- system_data %>%
            filter(HalfHoursElapsed == halfhour_period) %>%
            as_tibble()
          
          # Skip if no data
          if (nrow(data_system_halfhour) == 0) {
            message("         No data found. Skipping...")
            next
          }
          
          # ===================================================================
          # INITIALIZE LISTS
          # ===================================================================
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
          
          # ===================================================================
          # ANALYSIS: CALCULATE PROBABILITIES
          # ===================================================================
          animal_list <- initialize_animal_positions(animal_ids, data_system_halfhour, animal_list)
          
          initial_time <- animal_list[[1]][[2]]
          current_row <- 5
          elapsed_seconds <- 0
          total_rows <- nrow(data_system_halfhour) + 1
          
          # Process each row
          while (current_row != total_rows && current_row < total_rows) {
            
            previous_animal_positions <- animal_list
            
            animal_list <- update_animal_list(animal_ids, animal_list, data_system_halfhour, 
                                              initial_time, current_row)
            
            elapsed_seconds <- animal_list[["data_temp"]][["elapsed_seconds"]]
            
            # Update probabilities
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
          
          # ===================================================================
          # CALCULATE FINAL PROBABILITIES (HALF-HOUR)
          # ===================================================================
          for (i in 1:8) {
            # Handle division by zero - if no time was recorded, probability is 0
            if (cage_position_probability[[i]][[2]] > 0) {
              cage_position_probability[[i]][[4]] <- cage_position_probability[[i]][[3]] / 
                                                      cage_position_probability[[i]][[2]]
            } else {
              cage_position_probability[[i]][[4]] <- 0
            }
          }
          
          # Calculate animal probabilities, handling division by zero
          animal_position_probability <- animal_position_probability %>%
            mutate(Prob = ifelse(Seconds > 0, SumPercentage / Seconds, 0))
          
          # ===================================================================
          # CALCULATE SHANNON ENTROPY FOR HALF-HOUR
          # ===================================================================
          if (system_complete) {
            # Cage entropy
            cage_prob_vec <- sapply(cage_position_probability, function(x) x[4])
            
            # Remove NA/NaN values and filter out zeros for valid entropy calculation
            cage_prob_vec <- cage_prob_vec[!is.na(cage_prob_vec) & !is.nan(cage_prob_vec)]
            
            # Only calculate entropy if we have valid probabilities
            if (length(cage_prob_vec) > 0 && sum(cage_prob_vec) > 0) {
              # Normalize probabilities to sum to 1
              cage_prob_vec <- cage_prob_vec / sum(cage_prob_vec)
              cage_shannon_entropy <- calc_shannon_entropy(cage_prob_vec)
            } else {
              cage_shannon_entropy <- NA  # No valid data for this period
            }
            
            cagePosEntropy_halfhour <- cagePosEntropy_halfhour %>%
              add_row(Batch = batch, Sex = sex, System = system_id,
                      CageChange = cageChange,
                      HalfHour = halfhour_period,
                      CageEntropy = cage_shannon_entropy)
          }
          
          # Animal entropy
          for (animal in animal_ids) {
            if (is.na(animal)) next
            
            animal_prob_vec <- animal_position_probability %>%
              filter(AnimalID == animal) %>%
              pull(Prob)
            
            # Remove NA/NaN values and filter out zeros
            animal_prob_vec <- animal_prob_vec[!is.na(animal_prob_vec) & !is.nan(animal_prob_vec)]
            
            # Only calculate entropy if we have valid probabilities
            if (length(animal_prob_vec) > 0 && sum(animal_prob_vec) > 0) {
              # Normalize probabilities to sum to 1
              animal_prob_vec <- animal_prob_vec / sum(animal_prob_vec)
              animal_shannon_entropy <- calc_shannon_entropy(animal_prob_vec)
            } else {
              animal_shannon_entropy <- NA  # No valid data for this animal in this period
            }
            
            animalPosEntropy_halfhour <- animalPosEntropy_halfhour %>%
              add_row(Batch = batch, Sex = sex, System = system_id,
                      CageChange = cageChange,
                      HalfHour = halfhour_period,
                      AnimalID = animal,
                      animalEntropy = animal_shannon_entropy)
          }

        }
      }
    }
    
    # ===================================================================
    # SAVE RESULTS FOR THIS BATCH/CAGECHANGE
    # ===================================================================
    if (save_tables == TRUE) {
      message("=======================================================================")
      message("SAVING RESULTS FOR ", batch, " ", cageChange)
      message("=======================================================================")
      
      # Save phase-based results
      message("1. Saving phase-based analysis tables...")
      
      write.csv(cagePosProb, 
                file = paste0(saving_directory, tables_dir_phase, "/", 
                             batch, "_", cageChange, "_cagePosProb.csv"), 
                row.names = FALSE)
      message("   ✓ Saved: ", batch, "_", cageChange, "_cagePosProb.csv")
      
      write.csv(cagePosEntropy, 
                file = paste0(saving_directory, tables_dir_phase, "/", 
                             batch, "_", cageChange, "_cagePosEntropy.csv"), 
                row.names = FALSE)
      message("   ✓ Saved: ", batch, "_", cageChange, "_cagePosEntropy.csv")
      
      write.csv(animalPosEntropy, 
                file = paste0(saving_directory, tables_dir_phase, "/", 
                             batch, "_", cageChange, "_animalPosEntropy.csv"), 
                row.names = FALSE)
      message("   ✓ Saved: ", batch, "_", cageChange, "_animalPosEntropy.csv")
      message("      Location: ", tables_dir_phase)
      
      # Save half-hour results
      if (analyze_by_halfhour) {
        message("2. Saving half-hour analysis tables...")
        
        write.csv(cagePosEntropy_halfhour, 
                  file = paste0(saving_directory, tables_dir_halfhour, "/", 
                               batch, "_", cageChange, "_cagePosEntropy_halfhour.csv"), 
                  row.names = FALSE)
        message("   ✓ Saved: ", batch, "_", cageChange, "_cagePosEntropy_halfhour.csv")
        
        write.csv(animalPosEntropy_halfhour, 
                  file = paste0(saving_directory, tables_dir_halfhour, "/", 
                               batch, "_", cageChange, "_animalPosEntropy_halfhour.csv"), 
                  row.names = FALSE)
        message("   ✓ Saved: ", batch, "_", cageChange, "_animalPosEntropy_halfhour.csv")
        message("      Location: ", tables_dir_halfhour)
      }
    }
    
    # ===================================================================
    # ACCUMULATE IN GLOBAL TIBBLES
    # ===================================================================
    all_cagePosProb <- bind_rows(all_cagePosProb, cagePosProb)
    all_cagePosEntropy <- bind_rows(all_cagePosEntropy, cagePosEntropy)
    all_animalPosEntropy <- bind_rows(all_animalPosEntropy, animalPosEntropy)
    
    if (analyze_by_halfhour) {
      all_cagePosEntropy_halfhour <- bind_rows(all_cagePosEntropy_halfhour, cagePosEntropy_halfhour)
      all_animalPosEntropy_halfhour <- bind_rows(all_animalPosEntropy_halfhour, animalPosEntropy_halfhour)
    }
    
    message("=======================================================================")
    message("FINISHED ", batch, " ", cageChange)
    message("=======================================================================")
    message("")
  }
}

# ===================================================================
# SAVE COMBINED RESULTS (ALL BATCHES AND CAGE CHANGES)
# ===================================================================
if (save_tables == TRUE) {
  message("=======================================================================")
  message("SAVING COMBINED RESULTS (ALL BATCHES/CAGE CHANGES)")
  message("=======================================================================")
  
  # Create main tables directory
  combined_tables_dir <- paste0(tables_base, "/combined")
  dir.create(paste0(saving_directory, combined_tables_dir), recursive = TRUE, showWarnings = FALSE)
  
  # Save combined phase-based results
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
  
  # Save combined half-hour results
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

# ===================================================================
# DATA PREPROCESSING FOR CONSECUTIVE ENTROPY (PHASE-BASED)
# ===================================================================
message("=======================================================================")
message("PREPROCESSING CONSECUTIVE ENTROPY DATA (PHASE-BASED)")
message("=======================================================================")

consec_cage_entropy <- all_cagePosEntropy
consec_animal_entropy <- all_animalPosEntropy

# Process both tibbles
entropy_tibble_names <- c("consec_cage_entropy", "consec_animal_entropy")

for (i in seq_along(entropy_tibble_names)) {
  name <- entropy_tibble_names[i]
  entropy_list <- if (name == "consec_cage_entropy") consec_cage_entropy else consec_animal_entropy
  
  # Split Phase column
  entropy_list[c("Phase", "Consec")] <- str_split_fixed(entropy_list$Phase, '', 2)
  
  # Update consecutive columns
  entropy_list <- entropy_list %>%
    filter(Batch != "B6") %>%
    mutate(ConsecActive = ifelse(Phase == "A", as.numeric(Consec), 0)) %>%
    mutate(ConsecInactive = ifelse(Phase == "I", as.numeric(Consec) - 1, 0)) %>%
    mutate(Phase = ifelse(Phase == "A", "active", "inactive"))
  
  # Add Group column for animal entropy
  if (name == "consec_animal_entropy") {
    entropy_list <- entropy_list %>%
      mutate(Group = ifelse(AnimalID %in% sus_animals, "sus", 
                           ifelse(AnimalID %in% con_animals, "con", "res")))
  }
  
  # Adjust consecutive numbering across cage changes
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
  
  # Reorder columns
  if (name == "consec_animal_entropy") {
    entropy_list <- entropy_list[c("CageChange", "Batch", "System", "AnimalID",
                                   "Sex", "Group", "Phase", "ConsecActive",
                                   "ConsecInactive", "animalEntropy")]
  } else {
    entropy_list <- entropy_list[c("CageChange", "Batch", "System", "Sex",
                                   "Phase", "ConsecActive", "ConsecInactive",
                                   "CageEntropy")]
  }
  
  # Update original tibbles
  if (name == "consec_cage_entropy") {
    consec_cage_entropy <- entropy_list
  } else {
    consec_animal_entropy <- entropy_list
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
  message("      Location: ", combined_tables_dir)
}

# ===================================================================
# DATA PREPROCESSING FOR CONSECUTIVE ENTROPY (HALF-HOUR)
# ===================================================================
if (analyze_by_halfhour) {
  message("=======================================================================")
  message("PREPROCESSING CONSECUTIVE ENTROPY DATA (HALF-HOUR)")
  message("=======================================================================")
  
  consec_cage_entropy_halfhour <- all_cagePosEntropy_halfhour
  consec_animal_entropy_halfhour <- all_animalPosEntropy_halfhour
  
  # Process both half-hour tibbles
  entropy_halfhour_names <- c("consec_cage_entropy_halfhour", "consec_animal_entropy_halfhour")
  
  for (i in seq_along(entropy_halfhour_names)) {
    name <- entropy_halfhour_names[i]
    entropy_list <- if (name == "consec_cage_entropy_halfhour") {
      consec_cage_entropy_halfhour
    } else {
      consec_animal_entropy_halfhour
    }
    
    # Add Group column for animal entropy
    if (name == "consec_animal_entropy_halfhour") {
      entropy_list <- entropy_list %>%
        mutate(Group = ifelse(AnimalID %in% sus_animals, "sus", 
                             ifelse(AnimalID %in% con_animals, "con", "res")))
    }
    
    # Filter out B6
    entropy_list <- entropy_list %>%
      filter(Batch != "B6")
    
    # Create consecutive half-hour numbering across cage changes
    # Initialize ConsecHalfHour column
    entropy_list <- entropy_list %>%
      mutate(ConsecHalfHour = HalfHour)
    
    # Adjust consecutive numbering across cage changes
    max_halfhour <- 0
    
    for (change in c("CC1", "CC2", "CC3", "CC4")) {
      if (change != "CC1") {
        # Add the maximum from previous cage change
        entropy_list <- entropy_list %>%
          mutate(ConsecHalfHour = ifelse(CageChange == change, 
                                         HalfHour + max_halfhour, 
                                         ConsecHalfHour))
      }
      
      # Get the maximum half-hour value for this cage change
      max_halfhour <- entropy_list %>%
        filter(CageChange == change) %>%
        pull(ConsecHalfHour) %>%
        max(na.rm = TRUE)
    }
    
    # Reorder columns
    if (name == "consec_animal_entropy_halfhour") {
      entropy_list <- entropy_list[c("CageChange", "Batch", "System", "AnimalID",
                                     "Sex", "Group", "HalfHour", "ConsecHalfHour",
                                     "animalEntropy")]
    } else {
      entropy_list <- entropy_list[c("CageChange", "Batch", "System", "Sex",
                                     "HalfHour", "ConsecHalfHour", "CageEntropy")]
    }
    
    # Update original tibbles
    if (name == "consec_cage_entropy_halfhour") {
      consec_cage_entropy_halfhour <- entropy_list
    } else {
      consec_animal_entropy_halfhour <- entropy_list
    }
  }
  
  # ===================================================================
  # SAVE CONSECUTIVE HALF-HOUR ENTROPY TABLES
  # ===================================================================
  if (save_tables == TRUE) {
    message("Saving consecutive half-hour entropy tables...")
    
    combined_tables_dir <- paste0(tables_base, "/combined")
    
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
    message("      Location: ", combined_tables_dir)
  }
}

# ===================================================================
# OPTION 2: CREATE EXPLORATION METRICS (HOME CAGE EXCLUSION)
# ===================================================================
message("=======================================================================")
message("PROCESSING EXPLORATION METRICS (OPTION 2)")
message("=======================================================================")

# ===================================================================
# 1. ENHANCE PHASE-BASED DATA
# ===================================================================
message("Processing phase-based exploration metrics...")

consec_animal_entropy <- consec_animal_entropy %>%
  mutate(
    # Binary: Did animal explore outside home cage?
    explored = !is.na(animalEntropy),
    
    # Replace NA with 0 for analysis (0 = no exploration)
    exploration_entropy = ifelse(is.na(animalEntropy), 0, animalEntropy),
    
    # Categorize exploration behavior
    exploration_category = case_when(
      is.na(animalEntropy) ~ "No exploration",
      animalEntropy == 0 ~ "Minimal (1 position)",
      animalEntropy < 1.5 ~ "Low diversity",
      animalEntropy < 2.5 ~ "Moderate diversity",
      TRUE ~ "High diversity"
    ),
    
    # Make category a factor with proper ordering
    exploration_category = factor(exploration_category, 
                                  levels = c("No exploration", "Minimal (1 position)", 
                                           "Low diversity", "Moderate diversity", "High diversity"))
  )

# ===================================================================
# 2. ENHANCE HALF-HOUR DATA
# ===================================================================
if (analyze_by_halfhour) {
  message("Processing half-hour exploration metrics...")
  
  consec_animal_entropy_halfhour <- consec_animal_entropy_halfhour %>%
    mutate(
      # Binary: Did animal explore outside home cage?
      explored = !is.na(animalEntropy),
      
      # Replace NA with 0 for analysis
      exploration_entropy = ifelse(is.na(animalEntropy), 0, animalEntropy),
      
      # Categorize exploration behavior
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

# ===================================================================
# 3. CALCULATE EXPLORATION SUMMARY STATISTICS
# ===================================================================
message("Calculating exploration summary statistics...")

# Phase-based exploration summary
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

# Half-hour exploration summary
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

# Individual animal exploration profiles
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

message("\nTop 10 most exploratory animals:")
print(head(animal_exploration_profile, 10))

message("\nBottom 10 least exploratory animals:")
print(tail(animal_exploration_profile, 10))

# ===================================================================
# 4. SAVE ENHANCED DATASETS
# ===================================================================
if (save_tables == TRUE) {
  message("Saving enhanced exploration datasets...")
  
  combined_tables_dir <- paste0(tables_base, "/combined")
  
  # Save enhanced phase-based data
  write.csv(consec_animal_entropy, 
            file = paste0(saving_directory, combined_tables_dir, 
                         "/consec_animal_entropy_with_exploration_metrics.csv"), 
            row.names = FALSE)
  message("   ✓ Saved: consec_animal_entropy_with_exploration_metrics.csv")
  
  # Save enhanced half-hour data
  if (analyze_by_halfhour) {
    write.csv(consec_animal_entropy_halfhour, 
              file = paste0(saving_directory, combined_tables_dir, 
                           "/consec_animal_entropy_halfhour_with_exploration_metrics.csv"), 
              row.names = FALSE)
    message("   ✓ Saved: consec_animal_entropy_halfhour_with_exploration_metrics.csv")
  }
  
  # Save summary statistics
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
  
  # Save individual animal profiles
  write.csv(animal_exploration_profile, 
            file = paste0(saving_directory, combined_tables_dir, 
                         "/animal_exploration_profiles.csv"), 
            row.names = FALSE)
  message("   ✓ Saved: animal_exploration_profiles.csv")
  
  message("      Location: ", combined_tables_dir)
}

# ===================================================================
# 5. GENERATE EXPLORATION PLOTS
# ===================================================================
if (save_plots == TRUE || show_plots == TRUE) {
  message("=======================================================================")
  message("GENERATING EXPLORATION PLOTS")
  message("=======================================================================")
  
  combined_plots_dir <- paste0(plots_base, "/combined")
  
  # ---------------------------------------------------------------
  # PLOT 1: Exploration Frequency by Group
  # ---------------------------------------------------------------
  exploration_freq_plot <- ggplot(exploration_summary_phase, 
                                  aes(x = CageChange, y = pct_exploration, 
                                      color = Group, group = Group)) +
    geom_point(size = 4) +
    geom_line(linewidth = 1) +
    facet_grid(Phase ~ Sex) +
    scale_color_manual(values = c("sus" = "#E63946",
                                 "res" = "grey60",
                                 "con" = "#457B9D")) +
    labs(title = "Exploration Frequency by Group",
         subtitle = "% of observations where animals left home cage",
         x = "Cage Change",
         y = "% Exploration") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.key.size = unit(3, "lines"),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 16),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 14))
  
  # ---------------------------------------------------------------
  # PLOT 2: Exploration Diversity (when exploring)
  # ---------------------------------------------------------------
  data_explored_only <- consec_animal_entropy %>%
    filter(explored == TRUE)
  
  exploration_diversity_plot <- ggplot(data_explored_only, 
                                       aes(x = Group, y = animalEntropy, fill = Group)) +
    geom_violin(alpha = 0.6) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.alpha = 0.3) +
    facet_grid(Phase ~ Sex) +
    scale_fill_manual(values = c("sus" = "#E63946",
                                "res" = "grey60",
                                "con" = "#457B9D")) +
    labs(title = "Exploration Diversity (During Exploratory Bouts)",
         subtitle = "Entropy when animals left home cage",
         x = "Group",
         y = "Shannon Entropy") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.position = "none",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 14))
  
  # ---------------------------------------------------------------
  # PLOT 3: Exploration Category Distribution
  # ---------------------------------------------------------------
  exploration_category_plot <- consec_animal_entropy %>%
    count(Group, Sex, Phase, exploration_category) %>%
    group_by(Group, Sex, Phase) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    ggplot(aes(x = Group, y = pct, fill = exploration_category)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_grid(Phase ~ Sex) +
    scale_fill_manual(values = c("No exploration" = "#d62828",
                                "Minimal (1 position)" = "#f77f00",
                                "Low diversity" = "#fcbf49",
                                "Moderate diversity" = "#06a77d",
                                "High diversity" = "#003049"),
                     name = "Exploration Type") +
    labs(title = "Distribution of Exploration Categories",
         x = "Group",
         y = "Percentage (%)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.key.size = unit(2, "lines"),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 14))
  
  # ---------------------------------------------------------------
  # PLOT 4: Individual Animal Exploration Profiles
  # ---------------------------------------------------------------
  animal_profile_plot <- ggplot(animal_exploration_profile, 
                                aes(x = pct_exploration, 
                                    y = mean_entropy_when_exploring,
                                    color = Group)) +
    geom_point(size = 3, alpha = 0.6) +
    facet_wrap(~ Sex) +
    scale_color_manual(values = c("sus" = "#E63946",
                                 "res" = "grey60",
                                 "con" = "#457B9D")) +
    labs(title = "Individual Animal Exploration Profiles",
         subtitle = "Each point represents one animal",
         x = "Exploration Frequency (%)",
         y = "Mean Entropy When Exploring") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          plot.subtitle = element_text(hjust = 0.5, size = 16),
          legend.key.size = unit(3, "lines"),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          strip.text = element_text(size = 14))
  
  # ---------------------------------------------------------------
  # PLOT 5: Time-course of exploration (half-hour)
  # ---------------------------------------------------------------
  if (analyze_by_halfhour) {
    exploration_timecourse_plot <- consec_animal_entropy_halfhour %>%
      group_by(Group, Sex, ConsecHalfHour) %>%
      summarise(
        pct_explored = 100 * mean(explored),
        mean_entropy = mean(exploration_entropy, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      ggplot(aes(x = ConsecHalfHour, y = pct_explored, color = Group)) +
      geom_line(linewidth = 1) +
      geom_smooth(se = TRUE, alpha = 0.2) +
      facet_wrap(~ Sex, ncol = 1) +
      scale_color_manual(values = c("sus" = "#E63946",
                                   "res" = "grey60",
                                   "con" = "#457B9D")) +
      labs(title = "Exploration Frequency Over Time",
           subtitle = "% of animals exploring at each half-hour period",
           x = "Consecutive Half-Hour Period",
           y = "% Exploring") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 20),
            plot.subtitle = element_text(hjust = 0.5, size = 16),
            legend.key.size = unit(3, "lines"),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 14),
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            strip.text = element_text(size = 14))
  }
  
  # ---------------------------------------------------------------
  # SAVE PLOTS
  # ---------------------------------------------------------------
  if (save_plots) {
    ggsave(filename = paste0(saving_directory, combined_plots_dir, 
                            "/exploration_frequency_by_group.svg"),
           plot = exploration_freq_plot,
           width = 12,
           height = 10,
           dpi = 300)
    
    ggsave(filename = paste0(saving_directory, combined_plots_dir, 
                            "/exploration_diversity_when_exploring.svg"),
           plot = exploration_diversity_plot,
           width = 12,
           height = 10,
           dpi = 300)
    
    ggsave(filename = paste0(saving_directory, combined_plots_dir, 
                            "/exploration_category_distribution.svg"),
           plot = exploration_category_plot,
           width = 12,
           height = 10,
           dpi = 300)
    
    ggsave(filename = paste0(saving_directory, combined_plots_dir, 
                            "/animal_exploration_profiles.svg"),
           plot = animal_profile_plot,
           width = 12,
           height = 8,
           dpi = 300)
    
    if (analyze_by_halfhour) {
      ggsave(filename = paste0(saving_directory, combined_plots_dir, 
                              "/exploration_timecourse.svg"),
             plot = exploration_timecourse_plot,
             width = 14,
             height = 10,
             dpi = 300)
      
      message("   ✓ Saved 5 exploration plots")
    } else {
      message("   ✓ Saved 4 exploration plots")
    }
    
    message("      Location: ", combined_plots_dir)
  }
  
  # ---------------------------------------------------------------
  # SHOW PLOTS
  # ---------------------------------------------------------------
  if (show_plots) {
    print(exploration_freq_plot)
    print(exploration_diversity_plot)
    print(exploration_category_plot)
    print(animal_profile_plot)
    if (analyze_by_halfhour) {
      print(exploration_timecourse_plot)
    }
  }
}

message("=======================================================================")
message("EXPLORATION METRICS COMPLETE!")
message("=======================================================================")


# ===================================================================
# GENERATE PLOTS
# ===================================================================
if (save_plots == TRUE || show_plots == TRUE) {
  message("=======================================================================")
  message("GENERATING ENTROPY PLOTS")
  message("=======================================================================")
  
  # Create plots directory
  combined_plots_dir <- paste0(plots_base, "/combined")
  dir.create(paste0(saving_directory, combined_plots_dir), recursive = TRUE, showWarnings = FALSE)
  
  # ===================================================================
  # ANIMAL ENTROPY PLOTS (PHASE-BASED)
  # ===================================================================
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
    
    # Filter by phase if needed
    if ("Phase" %in% group_vars) {
      data <- filter(data, Phase == ifelse(phasecount == 1, "active", "inactive"))
      phasecount <- phasecount + 1
    }
    
    # Group and calculate mean
    result <- data %>%
      group_by(across(all_of(group_vars))) %>%
      summarise(Mean_Entropy = mean(animalEntropy, na.rm = TRUE), .groups = 'drop')
    
    # Create plot
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
  
  # ===================================================================
  # CAGE ENTROPY PLOTS (PHASE-BASED)
  # ===================================================================
  message("Generating phase-based cage entropy plots...")
  
  # Active phases
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
  
  # Inactive phases
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
  
    # ===================================================================
  # HALF-HOUR ENTROPY PLOTS (with consecutive numbering)
  # ===================================================================
  if (analyze_by_halfhour) {
    message("Generating half-hour entropy plots with consecutive numbering...")
    
    # Calculate cage change boundaries
    cage_boundaries <- consec_animal_entropy_halfhour %>%
      group_by(CageChange) %>%
      summarise(max_consec = max(ConsecHalfHour, na.rm = TRUE)) %>%
      pull(max_consec) %>%
      cumsum()
    
    # Remove the last boundary (end of data)
    if (length(cage_boundaries) > 0) {
      cage_boundaries <- cage_boundaries[-length(cage_boundaries)]
    }
    
    # Create data frame for annotations (if we have boundaries)
    if (length(cage_boundaries) > 0) {
      boundary_labels <- data.frame(
        x = cage_boundaries,
        label = paste("CC", 1:length(cage_boundaries), "→", 2:(length(cage_boundaries) + 1))
      )
    }
    
    # ===================================================================
    # PLOT 1: Animal Entropy Over Continuous Time
    # ===================================================================
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
    
    # Add vertical lines and labels if we have boundaries
    if (length(cage_boundaries) > 0) {
      animal_halfhour_plot_consec <- animal_halfhour_plot_consec +
        geom_vline(xintercept = cage_boundaries, linetype = "dashed", color = "black", alpha = 0.5) +
        geom_text(data = boundary_labels, 
                 aes(x = x, y = Inf, label = label),
                 inherit.aes = FALSE,
                 vjust = 1.5, size = 3, angle = 90, color = "black")
    }
    
    # ===================================================================
    # PLOT 2: Cage Entropy Over Continuous Time
    # ===================================================================
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
    
    # Add vertical lines and labels if we have boundaries
    if (length(cage_boundaries) > 0) {
      cage_halfhour_plot_consec <- cage_halfhour_plot_consec +
        geom_vline(xintercept = cage_boundaries, linetype = "dashed", color = "black", alpha = 0.5) +
        geom_text(data = boundary_labels, 
                 aes(x = x, y = Inf, label = label),
                 inherit.aes = FALSE,
                 vjust = 1.5, size = 3, angle = 90, color = "black")
    }
    
    # ===================================================================
    # PLOT 3: Animal Entropy by Cage Change (Original Style)
    # ===================================================================
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
    
    # ===================================================================
    # PLOT 4: Cage Entropy by Cage Change (Original Style)
    # ===================================================================
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
    
    # ===================================================================
    # SAVE HALF-HOUR PLOTS
    # ===================================================================
    if (save_plots) {
      # Save consecutive time plots
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
      
      # Save by cage change plots
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
      
      message("   ✓ Saved 4 half-hour entropy plots (2 consecutive + 2 by cage change)")
    }
    
    # Show plots if requested
    if (show_plots) {
      print(animal_halfhour_plot_consec)
      print(cage_halfhour_plot_consec)
      print(animal_halfhour_plot_bycage)
      print(cage_halfhour_plot_bycage)
    }
  }
}

message("=======================================================================")
message("ENTROPY ANALYSIS COMPLETE!")
message("=======================================================================")
