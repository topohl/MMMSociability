#' @title Analysis of Animal Positions
#' @description This script processes preprocessed animal position data, generating heatmaps, proximity, and movement analyses.
#' 
#' @details 
#' This script assumes a specific folder structure:
#' - `preprocessed_data/` must contain preprocessed CSV files.
#' - `plots/` must exist for saving plots if `save_plots == TRUE`.
#' - `tables/` must exist for saving tables if `save_tables == TRUE`.
#' 
#' Customizable parameters:
#' - `show_plots`: Display plots in R.
#' - `save_plots`: Save generated plots.
#' - `save_tables`: Save generated tables.
#' - `working_directory`: Set the working directory.
#' - `batch` and `cage change`: Define the filenames for processing.
#' 
#' @note 
#' The data in `data_preprocessed` has already been preprocessed in another script (`E9_SIS_AnimalPos-preprocessing_parallell.R`).
#' 
#' @authors 
#' Tobias Pohl, Anja Magister
#' 
#' @date 
#' February 2025

# Load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readr, dplyr, lubridate, tibble, purrr, ggplot2, reshape2, scales, stringr)

# Customizable variables
show_plots <- FALSE  # Set to TRUE to display plots in R
save_plots <- FALSE  # Set to TRUE to save generated plots
save_tables <- FALSE # Set to TRUE to save generated tables
exclude_homecage <- FALSE

# Paths
working_directory <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability"
saving_directory <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability"
plots_directory <- "/plots/noHomeCage"
tables_directory <- "/tables/noHomeCage"

# Source required functions
source(paste0(working_directory, "/E9_SIS_AnimalPos-functions.R"))

# Define batch and cage change identifiers
batches <- c("B1", "B2", "B3", "B4", "B5", "B6")
cageChanges <- c("CC1", "CC2", "CC3", "CC4")

# Initialize lists to store heatmaps and plots
allHeatmaps_proximity <- list()
allHeatmaps_positions <- list()
all_plots_total_proximity <- list()

# Loop through each batch and cage change
for(batch in batches) {
  
  for (cageChange in cageChanges) {
    print(paste(batch, cageChange))
   
    # Construct the filename for the current CSV file
    filename <- paste0("E9_SIS_", batch, "_", cageChange, "_AnimalPos")
    # Define the path to the CSV file
    csvFilePath <- paste0(working_directory, "/preprocessed_data/", filename, "_preprocessed.csv")
      
    # Read the preprocessed data from the CSV file into a tibble
    data_preprocessed <- as_tibble(read_delim(csvFilePath, delim = ",", show_col_types = FALSE))
  
    # Note: The data in "data_preprocessed" has already been preprocessed in another script.
    # ---------------------------------------------------
    # Define the unique systems, days, phases, and animals
    # ---------------------------------------------------

    unique_systems <- unique(data_preprocessed$System)
    # Sort the systems in a natural order
    unique_systems <- str_sort(unique_systems)
    
    # Extract unique days from the DateTime column in the preprocessed data
    days <- unique(format(data_preprocessed$DateTime, "%D"))
    
    # Define the phases of activity
    phases <- c("Active", "Inactive")
    
    # Identify the number of active and inactive phases in the current dataset
    active_phases_number <- unique(data_preprocessed$ConsecActive)
    inactive_phases_number <- unique(data_preprocessed$ConsecInactive)
        
    # Remove the phase number 0 as it is not relevant for analysis
    active_phases_number <- active_phases_number[active_phases_number != 0]
    inactive_phases_number <- inactive_phases_number[!inactive_phases_number %in% c(0, 1, max(inactive_phases_number))]
    
    # Extract unique animal IDs present in the current batch of data
    unique_animals <- unique(data_preprocessed$AnimalID)
    
    # Initialize lists to store heatmaps for each system
    systemHeatmaps_proximity <- list()
    systemHeatmaps_positions <- list()
    
    # Create a vector of phases in the correct order as they appear in reality
    phases_column <- c()
    for (i in 1:max(length(inactive_phases_number), length(active_phases_number))) {
      # Add active phases to the vector
      if (i <= length(active_phases_number)) {
        phases_column <- c(phases_column, paste0("A", i))
      }
        # Add inactive phases to the vector, starting with the second inactive phase
      if (i <= length(inactive_phases_number)) {
        phases_column <- c(phases_column, paste0("I", i + 1))
      }
    }
    
    # Create a tibble to store total proximity data for each animal in each phase
    result_total_proximity <- tibble("Phase" = phases_column)

    # Add columns for each unique animal ID to store proximity data
    for (animal in unique_animals) {
      result_total_proximity[[animal]] <- NA
    }
          
    # Create a tibble to store total movement data for each animal in each phase
    result_total_movement <- tibble("Phase" = phases_column)
          
    # Add columns for each unique animal ID to store movement data
    for (animal in unique_animals) {
      result_total_movement[[animal]] <- NA
    }
          
    # Add columns for each unique system ID to track overall movement within the system
    for (system in unique_systems) {
      result_total_movement[[system]] <- NA
    }
          
    # Create a tibble to store positional data for each system in each phase
    result_total_positions <- tibble("Phase" = phases_column)
        
    # Add columns for each unique system ID to track positional data within the system
    for (system in unique_systems) {
      result_total_positions[[system]] <- NA
    }
          
    # Create a tibble to store information about animal IDs in each system for the current cage change
    system_animal_ids <- tibble(
      "sys.1" = rep(NA, times = 4), 
      "sys.2" = rep(NA, times = 4), 
      "sys.3" = rep(NA, times = 4), 
      "sys.4" = rep(NA, times = 4), 
      "sys.5" = rep(NA, times = 4)
    )
    
    #' Loop Through Each System and Process Phase Data
    #' 
    #' This function processes the data for each system, calculating proximity, movement, and position for each phase (Active/Inactive)
    #' and generates the corresponding heatmaps and graphs.
    #' 
    #' @param data_preprocessed A tibble containing preprocessed data for all systems, animals, and phases.
    #' @param active_phases_number A vector of phase numbers for the Active phases.
    #' @param inactive_phases_number A vector of phase numbers for the Inactive phases.
    #' @param systems A vector of system identifiers.
    #' @param batch The current batch identifier.
    #' @param cageChange The identifier for the cage change condition.
    #' @param phases A vector containing the phase names ("Active", "Inactive").
    #' @param system_animal_ids A tibble containing the list of animal IDs for each system.
    #' @param result_total_proximity A tibble to store the total proximity data for each animal in each phase.
    #' @param result_total_movement A tibble to store the total movement data for each animal in each phase.
    #' @param result_total_positions A tibble to store the total position data for each system in each phase.
    #' @param allHeatmaps_proximity A list to store all heatmaps for proximity.
    #' @param allHeatmaps_positions A list to store all heatmaps for positions.
    #' @param all_plots_total_proximity A list to store all total closeness plots.
    #' @param save_plots A logical value indicating whether to save the generated plots.
    #' @param save_tables A logical value indicating whether to save the generated tables.
    #' @param show_plots A logical value indicating whether to show the generated plots in R.
    #' @param elapsed_seconds A numeric value to store the current time in seconds.
    #' @param data_temp A list to store temporary data.
    #' 
    #' @return A list containing all heatmaps for proximity and positions for all systems in the current batch.
    # Loop through each system
    for(system_id in unique_systems) {
      
      print(system_id)
      # Filter data for the current system
      data_system <- data_preprocessed %>%
        filter(System == system_id) %>%
        as_tibble()
      
      # Define animal IDs for the current system
      animal_ids <- unique(data_system$AnimalID)
      # Check if the system is complete (contains 4 animals)
      system_complete <- ifelse(length(animal_ids) < 4, FALSE, TRUE)
      # Fill vector with NAs if the system is incomplete
      while(length(animal_ids) < 4) { animal_ids <- append(animal_ids, NA) }
      
      # Add animal IDs to the system_animal_ids tibble
      system_animal_ids[[system_id]] <- c(animal_ids)
      
      # Loop through each phase (Active/Inactive)
      for(phase in phases) {
        
        print(phase) 
        
        # Select phase numbers based on the current phase type (Active/Inactive)
        phase_numbers <- ifelse(phase == "Active", active_phases_number, inactive_phases_number)
        print(phase_numbers)
        
        # Iterate through each phase number
        for(phase_number in phase_numbers) {
          print(paste0("System: ", system_id, ", ", phase, " phase phase_number: ", phase_number))
          
          # Filter data for the current phase and system
          data_system_phase <- data_system %>%
            filter(ConsecActive == ifelse(phase == "Active", phase_number, 0)) %>%
            filter(ConsecInactive == ifelse(phase == "Inactive", phase_number, 0)) %>%
            as_tibble()
          
          # Initialize lists for storing animal positions and temporary data
          animal_list <- list(
            "animal_1" = list(name = "", time = "", position = 0),
            "animal_2" = list(name = "", time = "", position = 0),
            "animal_3" = list(name = "", time = "", position = 0),
            "animal_4" = list(name = "", time = "", position = 0),
            "data_temp" = list(elapsed_seconds = 0, current_row = 0)
          )
          
          # Initialize lists for storing proximity, position, and movement results
          count_proximity_list <- list(m1 = c(0, 0, 0, 0), m2 = c(0, 0, 0, 0), m3 = c(0, 0, 0, 0), m4 = c(0, 0, 0, 0))
          count_position_list <- list(c(1, 0), c(2, 0), c(3, 0), c(4, 0), c(5, 0), c(6, 0), c(7, 0), c(8, 0))
          total_proximity_list <- list(c(animal_ids[1], 0), c(animal_ids[2], 0), c(animal_ids[3], 0), c(animal_ids[4], 0))
          count_movement_list <- list(c(animal_ids[1], 0), c(animal_ids[2], 0), c(animal_ids[3], 0), c(animal_ids[4], 0), c(system_id, 0))
          
          # Perform calculations for the current phase
          message("Calculating proximity, movement, and positional data for the current phase")
          
          # Initialize animal positions based on the current phase data
          animal_list <- initialize_animal_positions(animal_ids, data_system_phase, animal_list)
          
          # Assign initial time and line number for the while loop
          initial_time <- animal_list[[1]][[2]]
          current_row <- 5
          elapsed_seconds <- 0
          total_rows <- nrow(data_system_phase) + 1
          
          # Iterate through the rows of the current phase data
          while(current_row != total_rows && current_row < total_rows) {
            
            # Make a copy of the previous animal list
            previous_animal_list <- animal_list
            
            # Update animal list with new positions from the next row
            animal_list <- update_animal_list(animal_ids, animal_list, data_system_phase, initial_time, current_row)
            
            # Update elapsed_seconds
            elapsed_seconds <- animal_list[["data_temp"]][["elapsed_seconds"]]
            
            # Update result lists using analysis functions
            if(system_complete) { 
              count_proximity_list <- update_proximity(previous_animal_list, animal_list, count_proximity_list, elapsed_seconds)
            }

            count_position_list <- update_position(previous_animal_list, animal_list, count_position_list, elapsed_seconds)
            
            if(system_complete) {
              total_proximity_list <- update_total_proximity(previous_animal_list, animal_list, total_proximity_list, elapsed_seconds)
            }

            count_movement_list <- update_movement(previous_animal_list, animal_list, count_movement_list, elapsed_seconds)
            
            # Update current_row and initial_time for the next iteration
            current_row <- animal_list[["data_temp"]][["current_row"]]
            initial_time <- animal_list[[1]][[2]]
          }
          
          # Generate heatmaps
          message("Generating heatmaps for proximity and positional data")
          if(system_complete) {
            heatmap_closeness <- generateHeatMapProximity(count_proximity_list, batch, cageChange, system_id, animal_ids, phase, phase_number) 
          }

          heatmap_positions <- generateHeatMapPositions(count_position_list, batch, cageChange, system_id, phase, phase_number)
          
          # Save heatmaps in specific system list
          if(system_complete) {
            systemHeatmaps_proximity <- c(systemHeatmaps_proximity, list(heatmap_closeness))
          }

          systemHeatmaps_positions <- c(systemHeatmaps_positions, list(heatmap_positions))
          
          # Record total proximity data in the result tibble
          if(system_complete) {
            message("Recording total proximity data in tibble")
            for(i in 1:4) {
              animal <- total_proximity_list[[i]][[1]]
              p <- paste0(substr(phase, 1, 1), phase_number)
              row <- which(result_total_proximity$Phase == p)
              result_total_proximity[[animal]][[row]] <- total_proximity_list[[i]][[2]]
            }
          }
          
          # Enter movement data in the result tibble
          # Record movement data in the result tibble
          message("Recording movement data in tibble")
          p <- paste0(substr(phase, 1, 1), phase_number)
          row <- which(result_total_movement$Phase == p)
          for(i in 1:4) {
            animal <- count_movement_list[[i]][[1]]
            if(!is.na(animal)) {
              result_total_movement[[animal]][[row]] <- count_movement_list[[i]][[2]]
            } else {
              cat("Animal value from ", system_id, " in ", batch, " in ", cageChange, " is not available.")
            }
          }
          result_total_movement[[system_id]][[row]] <- count_movement_list[[5]][[2]]
          
        }
      }
      
      # Generate total proximity graph
      if(system_complete) {
        message("Generating total proximity graph")
        plot_total_proximity <- generateGraph(result_total_proximity, batch, cageChange, animal_ids, system_id)
        all_plots_total_proximity <- c(all_plots_total_proximity, list(plot_total_proximity))
      }
      
      # Save all heatmaps from one system in the general list
      if(system_complete) { allHeatmaps_proximity <- c(allHeatmaps_proximity, list(systemHeatmaps_proximity)) }
      message("Saving heatmaps from one system in the general list")
      print(batch)
      allHeatmaps_positions[[batch]][[cageChange]][[system_id]] <- systemHeatmaps_positions
      
      # Clear the lists of the current system
      systemHeatmaps_proximity <- list()
      systemHeatmaps_positions <- list()
    }

    # ---------------------------------------------------
    # Save Plots and Tables
    # ---------------------------------------------------

    # Save plots if save_plots is TRUE
    if (save_plots == TRUE) {
      message("Saving plots")
      
      # Save total closeness plots
      for (i in seq_along(all_plots_total_proximity)) {
        print(i)
        ggsave(
          filename = paste0(saving_directory, plots_directory, "/total_proximity_", batch, "_", cageChange, "_sys.", i, ".svg"), 
          plot = all_plots_total_proximity[[i]], 
          width = 5, 
          height = 2
        )
      }
      
      # Save heatmaps for each system
      for (i in seq_along(allHeatmaps_proximity)) {
        print(i)
        
        # Save proximity heatmaps
        title <- allHeatmaps_proximity[[i]][[1]]$labels$title
        pattern <- "sys..."
        system_substring <- ifelse((match <- regexec(pattern, title))[[1]][1] > 0, regmatches(title, match)[[1]], "Pattern not found.")
        ggsave(
          filename = paste0(saving_directory, plots_directory, "/allHeatmaps_proximity_", batch, "_", cageChange, "_", system_substring, ".png"), 
          plot = gridExtra::arrangeGrob(grobs = allHeatmaps_proximity[[i]], ncol = 2, layout_matrix = rbind(c(1, 5), c(2, 6), c(3, 7), c(4, 8), c(NA, 9))), 
          width = 12, 
          height = 8
        )
        
        # Save positional heatmaps
        title <- allHeatmaps_positions[[i]][[1]]$labels$title
        system_substring <- ifelse((match <- regexec(pattern, title))[[1]][1] > 0, regmatches(title, match)[[1]], "Pattern not found.")
        ggsave(
          filename = paste0(saving_directory, plots_directory, "/allHeatmaps_positions_", batch, "_", cageChange, "_", system_substring, ".png"), 
          plot = gridExtra::arrangeGrob(grobs = allHeatmaps_positions[[i]], ncol = 2, layout_matrix = rbind(c(1, 5), c(2, 6), c(3, 7), c(4, 8), c(NA, 9))), 
          width = 12, 
          height = 8
        )
      }
    }
    
    # Save tables as CSV files if save_tables is TRUE
    if (save_tables == TRUE) {
      message("Saving tables")
      
      # Save Social Proximity tables
      write.csv(result_total_proximity, 
                file = paste0(saving_directory, tables_directory, "/", batch, "_", cageChange, "_total_proximity.csv"), 
                row.names = FALSE)
      
      # Save Movement tables
      write.csv(result_total_movement, 
                file = paste0(saving_directory, tables_directory, "/", batch, "_", cageChange, "_total_movement.csv"), 
                row.names = FALSE)
      
      # Save System Animal ID tables
      write.csv(system_animal_ids, 
                file = paste0(saving_directory, tables_directory, "/", batch, "_", cageChange, "_animal_ids.csv"), 
                row.names = FALSE)
    }
  }
}

# ---------------------------------------------------
# Display Plots in R
# ---------------------------------------------------

if(show_plots==TRUE) {
  message("show Plots")

  # Create a grid of the total closeness plots
  gridExtra::grid.arrange(grobs = all_plots_total_proximity, ncol = 2)
  
  # Create a grid of the Heatmaps of each system for each analysis
  for (batch in seq_along(allHeatmaps_positions)) {

    for(cc in seq_along(allHeatmaps_positions[[batch]])) {

      for(sys in seq_along(allHeatmaps_positions[[batch]][[cc]])) {
        print(paste(sys))

        if(cc != "CC4") {
          gridExtra::grid.arrange(grobs = allHeatmaps_positions[[batch]][[cc]][[sys]], ncol = 2, layout_matrix = rbind(c(1,8), c(2,5), c(3,6), c(4,7)))
        } else {
          gridExtra::grid.arrange(grobs = allHeatmaps_positions[[batch]][[cc]][[sys]], ncol = 2, layout_matrix = rbind(c(1,10), c(2,6), c(3,7), c(4,8), c(5,9)))
        }
      }
    }
  }
}
