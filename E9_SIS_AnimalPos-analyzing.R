#' @title Analysis of Animal Positions and Proximity
#' @description This script analyzes the positions and proximity of animals in different systems (cages) over various phases (Active and Inactive). 
#' It processes preprocessed data files, generates heatmaps and plots for social proximity and positions, and saves the results.
#' The script iterates through batches and cage changes, performing the following steps:
#' 1. Load necessary libraries and set customizable variables.
#' 2. Define paths and source custom functions.
#' 3. Initialize result lists for heatmaps and plots.
#' 4. Iterate through each batch and cage change:
#'    a. Read the preprocessed data file.
#'    b. Define unique systems, days, phases, and animal IDs.
#'    c. Initialize result tibbles for proximity and movement.
#'    d. Iterate through each system:
#'       i. Filter data for the current system.
#'       ii. Initialize variables and lists for analysis.
#'       iii. Iterate through each phase (Active and Inactive):
#'            - Filter data for the current phase.
#'            - Initialize variables for analysis.
#'            - Perform calculations for proximity, positions, and movements.
#'            - Generate heatmaps for proximity and positions.
#'            - Update result tibbles with proximity and movement data.
#'       iv. Generate total proximity plots.
#'       v. Save heatmaps in general lists.
#'    e. Save plots and tables if specified.
#' 5. Optionally, display plots in R.
#'
#' @details
#' The script requires a specific file structure in the working directory:
#' - A folder called "preprocessed_data" containing the preprocessed data files.
#' - A folder called "plots" for saving the generated plots (if save_plots == TRUE).
#'
#' Customizable variables include:
#' - show_plots: Boolean to display plots in R.
#' - save_plots: Boolean to save the generated plots.
#' - save_tables: Boolean to save the generated tables.
#' - working_directory: Path to the current working directory.
#' - saving_directory: Path to the directory where plots and tables will be saved.
#' - batches: Vector of batch identifiers to process.
#' - cageChanges: Vector of cage change identifiers to process.
#'
#' @param show_plots Boolean to display plots in R.
#' @param save_plots Boolean to save the generated plots.
#' @param save_tables Boolean to save the generated tables.
#' @param excludeHomecage Boolean to exclude the homecage position from the analysis.
#' @param working_directory Path to the current working directory.
#' @param saving_directory Path to the directory where plots and tables will be saved.
#' @param batches Vector of batch identifiers to process.
#' @param cageChanges Vector of cage change identifiers to process.
#' @param processed_data Data frame containing the preprocessed data.
#' @param uniqueSystems Vector of unique system identifiers.
#' @param days Vector of unique days.
#' @param phases Vector of phases ("Active" and "Inactive").
#' @param animal_ids Vector of unique animal IDs.
#' @param systemData Data frame filtered for the current system.
#' @param systemPhaseData Data frame filtered for the current phase.
#' @param animal_list List containing the current state of animals.
#' @param count_proximity_list List storing proximity counts for each pair of mice.
#' @param count_position_list List storing position counts for each cage location.
#' @param total_proximity_list List storing total proximity counts for each mouse.
#' @param count_movement_list List storing movement counts for each mouse.
#' @param result_total_proximity Tibble to store total proximity results.
#' @param result_total_movement Tibble to store total movement results.
#' @param systemHeatmaps_proximity List to store heatmaps for the current system's social proximity.
#' @param systemHeatmaps_positions List to store heatmaps for the current system's positions.
#' @param allHeatmaps_proximity List to store all social proximity heatmaps.
#' @param allHeatmaps_positions List to store all position heatmaps.
#' @param all_plots_total_proximity List to store all total proximity plots.
#' @param system_animal_ids Tibble storing which mouse is in which system for each cage change.
#' @param filename String representing the current CSV filename.
#' @param csvFilePath String representing the path of the current CSV file.
#' @param active_phases_number Vector of unique active phase numbers.
#' @param inactive_phases_number Vector of unique inactive phase numbers.
#' @param phases_column Vector of phases in the correct chronological order.
#' @param system_complete Boolean indicating if the system data is complete.
#' @param timeTemp Temporary variable for time.
#' @param lineTemp Temporary variable for line number.
#' @param secTemp Temporary variable for time difference in seconds.
#' @param lastRow Integer representing the last row number in the data.
#' @param old_animal_list List containing the previous state of animals.
#' @param heatmap_proximity Heatmap for social proximity.
#' @param heatmap_positions Heatmap for positions.
#' @param plot_total_proximity Plot for total proximity.
#' @param p String representing the current phase identifier.
#' @param row Integer representing the row number in the result tibble.
#' @param mouse String representing the current mouse.
#' @param title String representing the title of the heatmap.
#' @param pattern String representing the pattern to identify the system substring.
#' @param system_substring String representing the system substring.
#' @param batch String representing the current batch identifier.
#' @param cageChange String representing the current cage change identifier.
#' @param system_id String representing the current system identifier.
#' @param phase String representing the current phase.
#' @param phase_number Integer representing the current phase number.

# Load necessary libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readr, dplyr, lubridate, tibble, purrr, ggplot2, reshape2, scales, stringr)

# Customizable variables
show_plots = FALSE
save_plots = FALSE
save_tables = TRUE
excludeHomecage = TRUE

# Paths
working_directory <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability"
saving_directory <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability"
plots_directory <- "/plots/noHomeCage"
tables_directory <- "/tables/noHomeCage"

# Create directories if they don't exist
if (!dir.exists(paste0(saving_directory, plots_directory))) {
  dir.create(paste0(saving_directory, plots_directory), recursive = TRUE)
}

if (!dir.exists(paste0(saving_directory, tables_directory))) {
  dir.create(paste0(saving_directory, tables_directory), recursive = TRUE)
}

# Source custom functions
source(paste0("C:/Users/topohl/Documents/GitHub/MMMSociability/E9_SIS_AnimalPos-functions.R"))

# Define batches and cage changes
batches <- c("B1", "B2", "B3", "B4", "B5", "B6")  # Add more batches if needed
cageChanges <- c("CC1", "CC2", "CC3", "CC4")  # Add more cage changes if needed

# Initialize result lists for heatmaps and plots
allHeatmaps_proximity <- list()
allHeatmaps_positions <- list()

# Initialize result tibbles for proximity and movement
all_plots_total_proximity <- list()

# Iterate through each batch and cage change
for (batch in batches) {
  for (cageChange in cageChanges) {
    print (paste(batch, cageChange))
    
    # Current csv filename
    filename <- paste0("E9_SIS_", batch, "_", cageChange, "_AnimalPos")
    # Path of csv file
    csvFilePath <-  paste0(working_directory,"/preprocessed_data/", filename, "_preprocessed.csv")
    
    # Read csv file in tibble
    processed_data <- as_tibble(read_delim(csvFilePath,delim = ",", show_col_types = FALSE))
    
    # Note: The data in "processed_data" has already been preprocessed in another script.
    
    # ---------------------------------------------------
    # Define unique systems, days, phases, and animal IDs
    # ---------------------------------------------------
    
    # Define unique systems
    uniqueSystems <- unique(processed_data$System)
    uniqueSystems <- str_sort(uniqueSystems)
    
    # Define days
    days <- unique(format(processed_data$DateTime, "%D"))
    
    # Define phases
    phases <- c("Active", "Inactive")
    
    # This section analyzes the number of active and inactive phases in the recorded data
    active_phases_number <- unique(processed_data$ConsecActive)
    inactive_phases_number <- unique(processed_data$ConsecInactive)
    
    # Remove the the first active and last inactive phase
    # This is done to exclude the first and last phases as they are incomplete
    # Remove the phase numbers 0, 1, and the maximum phase number from the inactive phases
    active_phases_number <- active_phases_number[! active_phases_number %in% 0]
    inactive_phases_number <- inactive_phases_number[! inactive_phases_number %in% c(0, 1, max(inactive_phases_number))]
    
    # Define Animal IDs
    animal_ids <- unique(processed_data$AnimalID)
    
    #' Initialize and Organize Result Lists
    #'
    #' @description
    #' This section handles the generation and organization of heatmaps across systems, 
    #' phases, and batches, ensuring structured visualization and storage.
    #'
    #' @details
    #' - Heatmaps are generated for each system, with segmentation based on phase and batch.
    #' - Each system-specific heatmap is stored in a dedicated list for modular organization.
    #' - All individual system lists are combined into a master list for holistic management.
    #' - The master list is used to produce a comprehensive grid of heatmaps spanning all systems.
    #' - The final heatmap grid is exported as a single image file for visualization purposes.
    #' - This process is iterated for every batch and cage change to ensure complete coverage.
    #'
    #' @note
    #' Ensure proper input data structure for accurate heatmap generation and alignment.
    #'
    #' @output
    #' A consolidated grid of heatmaps saved as an image and organized result lists for further analysis.
    systemHeatmaps_proximity <- list()
    systemHeatmaps_positions <- list()

    # Create a vector of phases in the correct chronological order as they appear in the data
    phases_column <- c()
    for (i in 1:max(length(inactive_phases_number), length(active_phases_number))) {
      if (i <= length(active_phases_number)) {
        phases_column <- c(phases_column, paste0("A", i))
      }
      if (i <= length(inactive_phases_number)) {
        phases_column <- c(phases_column, paste0("I", i + 1)) # we start with the second inactive phase
      }
    }
    
    # Tibble for proximity:
    # Create a tibble with a Phase column
    result_total_proximity <- tibble("Phase" = phases_column)

    # Add columns for each mouse present in the current batch
    for (mouse in animal_ids) {
      result_total_proximity[[mouse]] <- NA
    }
    
    # Tibble for movement:
    # This section creates a tibble with a Phase column
    result_total_movement <- tibble("Phase" = phases_column)
    
    # This section adds columns for each mouse present in the current batch.
    for (mouse in animal_ids) {
      result_total_movement[[mouse]] <- NA
    }

    # Add columns for each system to capture general movement within the system
    for (system in uniqueSystems) {
      result_total_movement[[system]] <- NA
    }
    
    # This section creates a tibble for storing animal positions.
    # This section creates a tibble with a Phase column and columns for each system to capture general movement within the system.
    result_total_positions <- tibble("Phase" = phases_column)
    for (system in uniqueSystems) {
      result_total_positions[[system]] <- NA
    }
    
    #tibble of information which mouse is in which system in this particular cagechange
    system_animal_ids <- tibble(
      "sys.1" = rep(NA, times = 4),
      "sys.2" = rep(NA, times = 4),
      "sys.3" = rep(NA, times = 4),
      "sys.4" = rep(NA, times = 4),
      "sys.5" = rep(NA, times = 4)
    )
    
    # ----------------------------
    # ANALYSIS
    # ----------------------------
    # Iterate through each system (representing different mouse cages).
    # There are 4 experimental systems and 1 control system.
    
    # Iterate through each system
    for (system_id in uniqueSystems) {
      
      print (system_id)
      # Filter the processed data to include only the current system
      systemData <- processed_data %>%
        filter(System==system_id) %>%
        as_tibble()
      
      # Define the animal IDs for the current system
      animal_ids <- unique(systemData$AnimalID)
      
      # Define if the system is complete (contains all 4 animals)
      system_complete <- length(animal_ids) >= 4
      
      # Fill the animal_ids vector with NA values if the system is incomplete
      while (length(animal_ids) < 4) {
        animal_ids <- append(animal_ids, NA)
      }
      
      # Add the animal IDs to the system_animal_ids tibble
      system_animal_ids[[system_id]] <- c(animal_ids)
      
      # Iterate through each phase (Active and Inactive) for the current system.
      for (phase in phases) {        
        # Determine the number of phases based on the current phase (Active or Inactive)
        ifelse (phase == "Active", number_of_phases <- active_phases_number, number_of_phases <- inactive_phases_number)
        
        # Iterate through each phase number for the current phase (Active or Inactive)
        for (phase_number in number_of_phases) {
          print (paste0("System: ", system_id, ", ", phase, " phase number: ", phase_number))
          
          # Filter the system data to include only the current phase
          systemPhaseData <- systemData %>%
            filter(ConsecActive == ifelse(phase == "Active", phase_number, 0)) %>%
            filter(ConsecInactive == ifelse(phase == "Inactive", phase_number, 0)) %>%
            as_tibble()

          #' Initialize Variables for Analysis
          #'
          #' @description
          #' Sets up the data structures and variables required for analyzing animal behavior 
          #' and proximity. This includes initializing lists to track individual animal states, 
          #' proximity counts, position counts, and movement counts.
          #'
          #' @details
          #' - **Animal State Initialization:**  
          #'   - Individual animal states are represented as lists containing the animal's name, 
          #'     start time, and start position.  
          #'   - These are consolidated into a master list (`animal_list`) along with temporary 
          #'     data for time difference and line number tracking.  
          #' - **Proximity Tracking:**  
          #'   - `count_proximity_list` records the cumulative seconds each pair of animals spent in proximity.  
          #'   - `total_proximity_list` tracks the total time each animal spent in proximity across all phases.  
          #' - **Position Tracking:**  
          #'   - `count_position_list` maintains a cumulative count of seconds spent at each position.  
          #'   - Positions can be dynamically adjusted by adding or removing entries for specific IDs.  
          #' - **Movement Tracking:**  
          #'   - `count_movement_list` logs the number of movements for each animal, including 
          #'     system-wide movements.
          #'
          #' @param animal_ids Vector. Unique identifiers for the animals being analyzed.
          #' @param system_id String. Identifier for the system under analysis.
          #'
          #' @return
          #' A set of initialized lists:
          #' - `animal_list`: Tracks the state of each animal and temporary analysis data.
          #' - `count_proximity_list`: Tracks seconds of proximity for each animal pair.
          #' - `count_position_list`: Tracks cumulative time at each position.
          #' - `total_proximity_list`: Tracks total proximity time for each animal.
          #' - `count_movement_list`: Tracks movement counts for each animal and system.
          #'
          #' @note
          #' Ensure that `animal_ids` contains at least four unique identifiers and that `system_id` 
          #' is appropriately defined before initializing variables.
          #'
          #' @examples
          #' # Initialize lists:
          #' animal_1 <- list(name = "", time = "", position = 0)
          #' animal_1 <- list(name = "", time = "", position = 0)
          #' animal_1 <- list(name = "", time = "", position = 0)
          #' animal_1 <- list(name = "", time = "", position = 0)
          #' 
          #' # Temporary data holder for time difference and line number in the analysis.
          #' data_temp <- list(secTemp = 0, lineTemp = 0)
          #' 
          #' animal_list <- list(
          #'   "animal_1" = animal_1,
          #'   "animal_2" = animal_2,
          #'   "animal_3" = animal_3,
          #'   "animal_4" = animal_4,
          #'   "data_temp" = data_temp
          #' )
          #' count_proximity_list <- list(
          #'   m1 = c(0, 0, 0, 0),
          #'   m2 = c(0, 0, 0, 0),
          #'   m3 = c(0, 0, 0, 0),
          #'   m4 = c(0, 0, 0, 0)
          #' )
          #' count_position_list <- list(
          #'   c(1, 0), c(2, 0), c(3, 0), c(4, 0), 
          #'   c(5, 0), c(6, 0), c(7, 0), c(8, 0)
          #' )
          #' total_proximity_list <- list(
          #'  c(animal_ids[1], 0), 
          #'  c(animal_ids[2], 0),
          #'  c(animal_ids[3], 0),
          #'  c(animal_ids[4], 0)
          #' )
          #' count_movement_list <- list(
          #'  c(animal_ids[1], 0), 
          #'  c(animal_ids[2], 0),
          #'  c(animal_ids[3], 0),
          #'  c(animal_ids[4], 0),
          #'  c(system_id, 0)
          #' )
          animal_1 <- list(name = "", time = "", position = 0)
          animal_2 <- list(name = "", time = "", position = 0)
          animal_3 <- list(name = "", time = "", position = 0)
          animal_4 <- list(name = "", time = "", position = 0)
          
          # Temporary data holder for time difference and line number in the analysis.
          data_temp <- list(secTemp = 0, lineTemp = 0
          )
          animal_list <- list(
            "animal_1" = animal_1,
            "animal_2" = animal_2,
            "animal_3" = animal_3,
            "animal_4" = animal_4,
            "data_temp" = data_temp
          )
          count_proximity_list <- list(
            m1 = c(0, 0, 0, 0),
            m2 = c(0, 0, 0, 0),
            m3 = c(0, 0, 0, 0),
            m4 = c(0, 0, 0, 0)
          )
          count_position_list <- list(
            c(1, 0), c(2, 0), c(3, 0), c(4, 0), 
            c(5, 0), c(6, 0), c(7, 0), c(8, 0)
          )
          total_proximity_list <- list(
            c(animal_ids[1], 0), 
            c(animal_ids[2], 0),
            c(animal_ids[3], 0),
            c(animal_ids[4], 0)
          )
          count_movement_list <- list(
            c(animal_ids[1], 0), 
            c(animal_ids[2], 0),
            c(animal_ids[3], 0),
            c(animal_ids[4], 0),
            c(system_id, 0)
          )

          #' Perform Calculations for Proximity and Position Analysis
          #'
          #' @description
          #' This section initializes variables and defines the logic required for proximity 
          #' and position analysis. It prepares the necessary data structures and sets up 
          #' parameters for iterating through the dataset.
          #'
          #' @details
          #' - **Initialization:**  
          #'   - The animal list is populated with the first recorded time and position 
          #'     for each animal using `find_first_pos_and_time()`.  
          #'   - A reference start time (`timeTemp`) is assigned based on the initial entry 
          #'     of the first animal in the list.  
          #' - **While Loop Setup:**  
          #'   - The starting line (`lineTemp`) is set to the first data row to be processed 
          #'     after the initial metadata lines.  
          #'   - The time difference between consecutive entries (`secTemp`) is initialized to 0.  
          #'   - The loop termination condition is defined as reaching the last row of the dataset.  
          #'
          #' @param animal_ids Vector. Unique identifiers for animals in the dataset.
          #' @param systemPhaseData Data frame. Contains the phase-specific data for the system, 
          #' including time and positional information for all animals.
          #' @param animal_list List. Stores information about each animal, including initial 
          #' time and position.
          #'
          #' @return
          #' - Updated `animal_list` with initial time and position data for each animal.  
          #' - Initialized variables (`timeTemp`, `lineTemp`, `secTemp`, and `lastRow`) 
          #'   for proximity and position calculations.
          #'
          #' @note
          #' Ensure that `systemPhaseData` is pre-processed and contains no missing values. 
          #' The `animal_ids` vector should correspond to unique identifiers present in `systemPhaseData`.
          #'
          #' @examples
          #' # Prepare data for analysis:
          #' animal_ids <- c("A1", "A2", "A3")
          #' systemPhaseData <- data.frame(Time = c("2025-01-01 12:00", "2025-01-01 12:01"), Position = c("P1", "P2"))
          #' animal_list <- list()
          #'
          #' # Initialize calculations:
          #' animal_list <- find_first_pos_and_time(animal_ids, systemPhaseData, animal_list)
          #' timeTemp <- animal_list[[1]][[2]]
          #' lineTemp <- 5
          #' secTemp <- 0
          #' lastRow <- nrow(systemPhaseData) + 1
          
          # Initialize calculations:
          message ("Calculations")
          animal_list <- find_first_pos_and_time(animal_ids, systemPhaseData, animal_list)
          timeTemp <- animal_list[[1]][[2]]
          lineTemp <- 5
          secTemp <- 0
          lastRow <- nrow(systemPhaseData) + 1

          #' @description This script iterates through system data to update the animal list and perform various analyses, including social proximity, cage location usage, total proximity per mouse, and movement counts. The results are stored in dedicated lists and used to generate heatmaps and plots.
          #'
          #' @details
          #' 1. Iterate through the system data to update the animal list and perform the analysis.
          #' 2. Create a copy of the previous state of the animal list to compare with the updated list.
          #' 3. Update the animal_list with new positions and times from the next row of systemPhaseData.
          #' 4. Update the time difference between the current and previous entries.
          #' 5. Use analysis functions to update the result lists:
          #'    - Compute social proximity between each pair of mice.
          #'    - Track the usage of cage locations.
          #'    - Aggregate the total social proximity for each mouse.
          #'    - Calculate the movement count for each mouse.
          #' 6. Update lineTemp and timeTemp for the next iteration.
          #' 7. Print results if needed.
          #' 8. Generate heatmaps for social proximity and positions.
          #' 9. Save heatmaps in specific system lists.
          #' 10. Integrate total proximity data into the result tibble.
          #' 11. Integrate movement data into the result tibble.
          #' 12. Generate total proximity plot for each system.
          #' 13. Save all heatmaps from one system in the general list.
          #' 14. Clear the current system's lists for the next iteration.
          #'
          #' @param animal_ids List of animal IDs.
          #' @param animal_list List containing the current state of animals.
          #' @param systemPhaseData Data frame containing system phase data.
          #' @param timeTemp Temporary variable for time.
          #' @param lineTemp Temporary variable for line number.
          #' @param lastRow The last row number in the data.
          #' @param system_complete Boolean indicating if the system data is complete.
          #' @param count_proximity_list List storing proximity counts for each pair of mice.
          #' @param count_position_list List storing position counts for each cage location.
          #' @param total_proximity_list List storing total proximity counts for each mouse.
          #' @param count_movement_list List storing movement counts for each mouse.
          #' @param batch Batch identifier.
          #' @param cageChange Cage change identifier.
          #' @param system_id System identifier.
          #' @param phase Current phase.
          #' @param phase_number Current phase number.
          #' @param result_total_proximity Tibble to store total proximity results.
          #' @param result_total_movement Tibble to store total movement results.
          #' @param all_plots_total_proximity List to store all total proximity plots.
          #' @param allHeatmaps_proximity List to store all social proximity heatmaps.
          #' @param allHeatmaps_positions List to store all position heatmaps.
          #' @param systemHeatmaps_proximity List to store heatmaps for the current system's social proximity.
          #' @param systemHeatmaps_positions List to store heatmaps for the current system's positions.
          
          # Iterate through the system data to update the animal list and perform the analysis.
          while (lineTemp != lastRow && lineTemp < lastRow) {
            
            # Create a copy of the previous state of the animal list 
            # to compare with the updated list (differences between two consecutive rows in the data)
            old_animal_list <- animal_list
            
            # Update the animal_list with new positions and times from the next row of systemPhaseData
            animal_list <- update_animal_list(animal_ids, animal_list, systemPhaseData, timeTemp, lineTemp)
            
            # Update the time difference between the current and previous entries
            secTemp <- animal_list[["data_temp"]][["secTemp"]]
            
            ##################################################################################################
            # Use analysis functions to update the result lists
            # These functions process the system data to compute social proximity, cage location usage, 
            # total proximity per mouse, and movement counts. Results are stored in dedicated lists.
            ##################################################################################################

            # Update the social proximity data between each pair of mice.
            # This function calculates the time mice spend in proximity to each other and updates the list.
            if (system_complete) {
              count_proximity_list <- compute_proximity(
                old_animal_list,  # The state of animals from the previous time step
                animal_list,      # The current state of animals
                count_proximity_list,  # List storing proximity counts for each pair of mice
                secTemp           # Time difference in seconds between the current and previous state
              )
            }

            # Update the usage of cage locations.
            # This function tracks how much time each mouse spends in different cage positions.
            count_position_list <- compute_position(
              old_animal_list,  # The state of animals from the previous time step
              animal_list,      # The current state of animals
              count_position_list,  # List storing position counts for each cage location
              secTemp           # Time difference in seconds between the current and previous state
            )

            # Update the total social proximity for each mouse.
            # This function aggregates the total time each mouse spends in proximity to others.
            if (system_complete) {
              total_proximity_list <- compute_total_proximity(
                old_animal_list,  # The state of animals from the previous time step
                animal_list,      # The current state of animals
                total_proximity_list,  # List storing total proximity counts for each mouse
                secTemp           # Time difference in seconds between the current and previous state
              )
            }

            # Update the movement count for each mouse.
            # This function calculates the number of coil crossings performed by each mouse.
            count_movement_list <- compute_movement(
              old_animal_list,  # The state of animals from the previous time step
              animal_list,      # The current state of animals
              count_movement_list,  # List storing movement counts for each mouse
              secTemp           # Time difference in seconds between the current and previous state
            )

            # Update lineTemp and timeTemp for the next iteration
            # These variables control the loop progression and are updated based on the current data row.
            lineTemp <- animal_list[["data_temp"]][["lineTemp"]]
            timeTemp <- animal_list[[1]][[2]]
          }
          
          # --------------------
          # PRINT RESULTS
          # --------------------

          if (FALSE) {
          message ("results: ")
          message ("count_proximity_list")
          print (count_proximity_list)
          message ("count_position_list")
          print (count_position_list)
          message ("total_proximity_list")
          print (total_proximity_list)
          }
          
          #print (paste("heatmap", system_id))
          message ("Generate heatmaps...")
          
          # Generate heatmaps
          if (system_complete) {
            heatmap_proximity <- generateHeatMapproximity(count_proximity_list, batch, cageChange, system_id, animal_ids, phase, phase_number)
          }

          # Positions
          heatmap_positions <- generateHeatMapPositions(count_position_list, batch, cageChange, system_id, phase, phase_number)
          #print (heatmap_positions)
          
          # Add the generated heatmaps to the system lists
          # Social Proximity
          if (system_complete) {
            systemHeatmaps_proximity <- c(systemHeatmaps_proximity, list(heatmap_proximity))
          }

          # Positions
          systemHeatmaps_positions <- c(systemHeatmaps_positions, list(heatmap_positions))
          
          if (system_complete) {
            message("Entering total proximity data into tibble")

            # Integrate total_proximity_list data into the result tibble
            for (i in 1:4) { # For each mouse in the system
              mouse <- total_proximity_list[[i]][[1]] # Current mouse
              phase_identifier <- paste0(substr(phase, 1, 1), phase_number) # Current phase identifier
              
              # Find the row corresponding to the current phase
              row <- which(result_total_proximity$Phase == phase_identifier)
              result_total_proximity[[mouse]][[row]] <- total_proximity_list[[i]][[2]]
            }
          }
            
          message ("Entering total movement data into tibble")

          # This section of the code integrates the count_movement_list information into a comprehensive result tibble.
          # The result tibble is structured to include data for every phase and system, which will be utilized for subsequent plotting.

          # Concatenate the first character of the phase variable with the phase_number variable to create the current phase identifier.
          p <- paste0(substr(phase, 1, 1), phase_number)

          # Find the row corresponding to the current phase
          row <- which(result_total_movement$Phase == p)

          # Update the result_total_movement tibble with movement data for each mouse
          for (i in 1:4) {
            mouse <- count_movement_list[[i]][[1]]
            if (!is.na(mouse)) {
              result_total_movement[[mouse]][[row]] <- count_movement_list[[i]][[2]]
            } else {
                message(paste("Animal value from", system_id, "in", batch, "in", cageChange, "is not available."))
            }
          }

          # Update the result_total_movement tibble with movement data for the system
          result_total_movement[[system_id]][[row]] <- count_movement_list[[5]][[2]]
        }
      }

      if (system_complete) {
        message ("Generating total proximity plot")  
        # Generate total proximity plot for each system
        # Parameters: system, mouse names, result tibble
        plot_total_proximity <- generateGraph(result_total_proximity, batch, cageChange, animal_ids, system_id)
        
        # Add the generated plot to the list of plots
        all_plots_total_proximity <- c(all_plots_total_proximity, list(plot_total_proximity) )
      }
      
      # Save all heatmaps from one system in the general list
      # Social Proximity
      if (system_complete) {
        allHeatmaps_proximity <- c(allHeatmaps_proximity, list(systemHeatmaps_proximity))
      }

      # Positions
      # allHeatmaps_positions <- c(allHeatmaps_positions, list(systemHeatmaps_positions))
      message ("Save heatmaps from one system in the general list")
      print (batch)
      allHeatmaps_positions[[batch]][[cageChange]][[system_id]] <- systemHeatmaps_positions
      
      # Clear the current system's lists
      systemHeatmaps_proximity <- list()
      systemHeatmaps_positions <- list()
    }
    
    #' Save Plots and Tables
    #'
    #' @description
    #' This section handles the saving of plots and tables to user-specified directories. 
    #' The process is performed for each batch and cage change, ensuring all relevant 
    #' visualizations and data summaries are stored systematically.
    #'
    #' @details
    #' - **Plots:**  
    #'   - Total proximity plots are generated and saved for each system.  
    #'   - Heatmaps for proximity and position analysis are saved, with filenames 
    #'     including system identifiers for clarity.  
    #'   - Layouts of heatmaps are managed using `gridExtra` for better visualization.  
    #' - **Tables:**  
    #'   - Social proximity, movement, and system animal ID data are exported as CSV files.  
    #'   - Files are named to include batch and cage change details for easy identification.
    #'
    #' @param save_plots Logical. If `TRUE`, plots are saved to the specified directory.
    #' @param save_tables Logical. If `TRUE`, tables are saved to the specified directory.
    #' @param saving_directory String. The root directory where plots and tables will be saved.
    #' @param plots_directory String. Subdirectory for saving plot files.
    #' @param tables_directory String. Subdirectory for saving table files.
    #' @param batch String. Identifier for the current batch being processed.
    #' @param cageChange String. Identifier for the current cage change being processed.
    #' @param all_plots_total_proximity List. Contains total proximity plots for each system.
    #' @param allHeatmaps_proximity List. Contains proximity heatmaps for each system.
    #' @param allHeatmaps_positions List. Contains position heatmaps for each system.
    #' @param result_total_proximity Data frame. Results of total proximity analysis.
    #' @param result_total_movement Data frame. Results of total movement analysis.
    #' @param system_animal_ids Data frame. Animal IDs corresponding to each system.
    #'
    #' @output
    #' - Saved `.svg` files for total proximity plots and heatmaps (proximity and position).  
    #' - Exported `.csv` files for proximity, movement, and animal ID data.
    #'
    #' @note
    #' Ensure `save_plots` and `save_tables` flags are appropriately set before running this section.
    #' Directories specified in `saving_directory`, `plots_directory`, and `tables_directory` 
    #' should exist or be created beforehand to avoid errors.
    #'
    #' @examples
    #' # Save plots and tables for a given batch and cage change:
    #' save_plots <- TRUE
    #' save_tables <- TRUE
    #' saving_directory <- "output"
    #' plots_directory <- "plots"
    #' tables_directory <- "tables"
    #' # Assume data and plots are pre-defined:
    #' save_results(save_plots, save_tables, saving_directory, plots_directory, tables_directory, batch, cageChange, 
    #'              all_plots_total_proximity, allHeatmaps_proximity, allHeatmaps_positions, 
    #'              result_total_proximity, result_total_movement, system_animal_ids)

    # Generate and save plots if the `save_plots` flag is set to TRUE.
    if (save_plots == TRUE) {
      message("Saving plots...")
      
      # Save total proximity plots.
      # Each plot corresponds to a system's total proximity data.
      for (i in seq_along(all_plots_total_proximity)) {
        print(i)  # Log the plot index being saved.
        ggsave(
          filename = paste0(
            saving_directory, plots_directory, "/total_proximity_", batch, "_", cageChange, "_sys.", i, ".svg"),
          plot = all_plots_total_proximity[[i]],
          width = 5,
          height = 2
        )
      }
      
      # Save heatmaps for proximity and position analysis.
      # Iterate through each system's heatmap in the `allHeatmaps_proximity` and `allHeatmaps_positions` lists.
      for (i in seq_along(allHeatmaps_proximity)) {
        print(i)  # Log the heatmap index being processed.
        
        # Extract the system identifier from the title of the proximity heatmap.
        title <- allHeatmaps_proximity[[i]][[1]]$labels$title
        pattern <- "sys.."  # Define the pattern to identify the system substring.
        system_substring <- ifelse(
          (match <- regexec(pattern, title))[[1]][1] > 0,
          regmatches(title, match)[[1]],
          "Pattern not found."
        )
        
        # Save heatmap for proximity analysis.
        ggsave(
          filename = paste0(
            saving_directory, plots_directory, "/allHeatmaps_proximity_", batch, "_", cageChange, "_", system_substring, ".svg"),
          plot = gridExtra::arrangeGrob(
            grobs = allHeatmaps_proximity[[i]],
            ncol = 2,
            layout_matrix = rbind(c(1, 5), c(2, 6), c(3, 7), c(4, 8), c(NA, 9))),
          width = 12,
          height = 8
        )
        
        # Extract the system identifier from the title of the position heatmap.
        title <- allHeatmaps_positions[[i]][[1]]$labels$title
        system_substring <- ifelse(
          (match <- regexec(pattern, title))[[1]][1] > 0,
          regmatches(title, match)[[1]],
          "Pattern not found."
        )
        
        # Save heatmap for position analysis.
        ggsave(
          filename = paste0(
            saving_directory, plots_directory, "/allHeatmaps_positions_", batch, "_", cageChange, "_", system_substring, ".svg"),
          plot = gridExtra::arrangeGrob(
            grobs = allHeatmaps_positions[[i]],
            ncol = 2,
            layout_matrix = rbind(c(1, 5), c(2, 6), c(3, 7), c(4, 8), c(NA, 9))),
          width = 12,
          height = 8
        )
      }
    }

    # Save tables
    if (save_tables == TRUE) {
      message ("save tables")
      
      # Saving Social Proximity tables
      write.csv(result_total_proximity, file = paste0(saving_directory, tables_directory,"/", batch, "_", cageChange, "_total_proximity.csv"), row.names = FALSE)
      # Saving Movement tables
      write.csv(result_total_movement, file = paste0(saving_directory, tables_directory,"/", batch, "_", cageChange, "_total_movement.csv"), row.names = FALSE)
      # Saving System Animal ID tables
      write.csv(system_animal_ids, file = paste0(saving_directory, tables_directory,"/", batch, "_", cageChange, "_animal_ids.csv"), row.names = FALSE)
    }
  }
}

# ----------------------------
# Display plots in R
# ----------------------------

if (show_plots == TRUE) {  
  # Create a grid of the total proximity plots
  gridExtra::grid.arrange(grobs = all_plots_total_proximity, ncol = 2)
  # Create a grid of the Heatmaps of each system for each analysis
  for (batch in seq_along(allHeatmaps_positions)) {
    for (cc in seq_along(allHeatmaps_positions[[batch]])) {
      for (sys in seq_along(allHeatmaps_positions[[batch]][[cc]])) {
        if (cc != "CC4") {
          gridExtra::grid.arrange(grobs = allHeatmaps_positions[[batch]][[cc]][[sys]], ncol = 2, layout_matrix = rbind(c(1, 8), c(2, 5), c(3, 6), c(4, 7)))
        } else {
          gridExtra::grid.arrange(grobs = allHeatmaps_positions[[batch]][[cc]][[sys]], ncol = 2, layout_matrix = rbind(c(1, 10), c(2, 6), c(3, 7), c(4, 8), c(5, 9)))
        }
      }  
    }
  }
}