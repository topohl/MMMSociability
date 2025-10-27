#' @title Analysis of Animal Positions
#' @description This script processes preprocessed animal position data,
#' generating heatmaps, proximity, and movement analyses.
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
#' - `analyze_by_halfhour`: Set to TRUE to also analyze by half-hour periods
#'
#' @note
#' The data in `data_preprocessed` has already been preprocessed in
#' another script (`E9_SIS_AnimalPos-preprocessing_parallell.R`).
#'
#' @authors
#' Tobias Pohl, Anja Magister
#'
#' @date
#' October 2025

# Load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readr, dplyr, lubridate, tibble, purrr, ggplot2, reshape2,
               scales, stringr)

# Customizable variables
show_plots <- FALSE  # Set to TRUE to display plots in R
save_plots <- TRUE  # Set to TRUE to save generated plots
save_tables <- TRUE # Set to TRUE to save generated tables
exclude_homecage <- TRUE  # Set to TRUE to exclude home cage positions
analyze_by_halfhour <- TRUE  # Set to TRUE to also analyze by half-hour periods

# Paths
#working_directory <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability"
#saving_directory <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability"

working_directory <- "D:/MMMSociability"
saving_directory <- "D:/MMMSociability"

# Define base directories
plots_base <- paste0("/plots/", ifelse(exclude_homecage, "noHomeCage", "withHomeCage"))
tables_base <- paste0("/tables/", ifelse(exclude_homecage, "noHomeCage", "withHomeCage"))

# Source required functions
#source(paste0(working_directory, "/E9_SIS_AnimalPos-functions.R"))
source("C:/Users/Tobias Pohl/Documents/GitHub/MMMSociability/Functions/E9_SIS_AnimalPos-functions.R")

# Define batch and cage change identifiers
batches <- c("B1", "B2", "B3", "B4", "B5", "B6")
cageChanges <- c("CC1", "CC2", "CC3", "CC4")

# Initialize lists to store heatmaps and plots
allHeatmaps_proximity <- list()
allHeatmaps_positions <- list()
all_plots_total_proximity <- list()

# Loop through each batch and cage change
for (batch in batches) {
  for (cageChange in cageChanges) {
    print(paste(batch, cageChange))

    # ===================================================================
    # CREATE FOLDER STRUCTURE FOR CURRENT BATCH AND CAGE CHANGE
    # ===================================================================
    # Define subdirectories for this specific batch and cage change
    current_plots_dir <- paste0(plots_base, "/", batch, "/", cageChange)
    current_tables_dir <- paste0(tables_base, "/", batch, "/", cageChange)
    
    # Create phase-based analysis directories
    plots_dir_phase <- paste0(current_plots_dir, "/by_phase")
    tables_dir_phase <- paste0(current_tables_dir, "/by_phase")
    
    # Create half-hour analysis directories
    plots_dir_halfhour <- paste0(current_plots_dir, "/by_halfhour")
    tables_dir_halfhour <- paste0(current_tables_dir, "/by_halfhour")
    
    # Create metadata directory
    tables_dir_metadata <- paste0(current_tables_dir, "/metadata")
    
    # Create all directories
    dir.create(paste0(saving_directory, plots_dir_phase), recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0(saving_directory, tables_dir_phase), recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0(saving_directory, tables_dir_metadata), recursive = TRUE, showWarnings = FALSE)
    
    if (analyze_by_halfhour) {
      dir.create(paste0(saving_directory, plots_dir_halfhour), recursive = TRUE, showWarnings = FALSE)
      dir.create(paste0(saving_directory, tables_dir_halfhour), recursive = TRUE, showWarnings = FALSE)
    }
    
    message("✓ Created directory structure for ", batch, " ", cageChange)

    # Construct the filename for the current CSV file
    filename <- paste0("E9_SIS_", batch, "_", cageChange, "_AnimalPos")
    # Define the path to the CSV file
    csvFilePath <- paste0(working_directory, "/preprocessed_data_test/", filename, "_preprocessed.csv")

    # Read the preprocessed data from the CSV file into a tibble
    data_preprocessed <- as_tibble(read_delim(csvFilePath, delim = ",",
                                              show_col_types = FALSE))

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
    # Preprocessing has already removed phase 0, first inactive, and last inactive
    active_phases_number <- unique(data_preprocessed$ConsecActive)
    inactive_phases_number <- unique(data_preprocessed$ConsecInactive)

    # remove active and inactive phases > 2 from CC4
    if (cageChange == "CC4") {
      active_phases_number <- active_phases_number[active_phases_number <= 2]
      inactive_phases_number <- inactive_phases_number[inactive_phases_number <= 2]
    } 

    # Remove any remaining 0 values (safety check)
    active_phases_number <- active_phases_number[active_phases_number != 0]
    inactive_phases_number <- inactive_phases_number[inactive_phases_number != 0]


    # Sort to ensure correct order
    active_phases_number <- sort(active_phases_number)
    inactive_phases_number <- sort(inactive_phases_number)

    print(paste("Active phases in data:", paste(active_phases_number, collapse = ", ")))
    print(paste("Inactive phases in data:", paste(inactive_phases_number, collapse = ", ")))

    # Get unique half-hour periods
    halfhour_periods <- sort(unique(data_preprocessed$HalfHoursElapsed))
    print(paste("Half-hour periods in data:", paste(halfhour_periods, collapse = ", ")))

    # Extract unique animal IDs present in the current batch of data
    unique_animals <- unique(data_preprocessed$AnimalID)

    # Initialize lists to store heatmaps for each system
    systemHeatmaps_proximity <- list()
    systemHeatmaps_positions <- list()

    # Create a vector of phases using ACTUAL phase numbers from the data
    phases_column <- c()
    max_phases <- max(length(active_phases_number), length(inactive_phases_number))

    for (i in 1:max_phases) {
      # Add active phase with its actual phase number
      if (i <= length(active_phases_number)) {
        phases_column <- c(phases_column, paste0("A", active_phases_number[i]))
      }
      # Add inactive phase with its actual phase number
      if (i <= length(inactive_phases_number)) {
        phases_column <- c(phases_column, paste0("I", inactive_phases_number[i]))
      }
    }
    print(paste("Phase column created:", paste(phases_column, collapse = ", ")))

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

    # ===================================================================
    # CREATE HALF-HOUR ANALYSIS TIBBLES
    # ===================================================================
    if (analyze_by_halfhour) {
      # Create column labels for half-hour periods (H0, H1, H2, etc.)
      halfhour_column <- paste0("H", halfhour_periods)
      
      # Create tibbles for half-hour analysis
      result_halfhour_proximity <- tibble("HalfHour" = halfhour_column)
      result_halfhour_movement <- tibble("HalfHour" = halfhour_column)
      
      # Add columns for each animal
      for (animal in unique_animals) {
        result_halfhour_proximity[[animal]] <- NA
        result_halfhour_movement[[animal]] <- NA
      }
      
      # Add columns for each system in movement
      for (system in unique_systems) {
        result_halfhour_movement[[system]] <- NA
      }
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
    for (system_id in unique_systems) {

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
      while (length(animal_ids) < 4) {
        animal_ids <- append(animal_ids, NA)
      }

      # Add animal IDs to the system_animal_ids tibble
      system_animal_ids[[system_id]] <- c(animal_ids)

      # ===================================================================
      # ANALYSIS BY PHASE (Active/Inactive)
      # ===================================================================
      # Loop through each phase (Active/Inactive)
      for (phase in phases) {
        print(phase)

        # Select phase numbers based on the current phase type (Active/Inactive)
        if (phase == "Active") {
          phase_numbers <- active_phases_number
        } else {
          phase_numbers <- inactive_phases_number
        }

        print(paste("Phase numbers to process:", paste(phase_numbers, collapse = ", ")))

        # Iterate through each phase number
        for (phase_number in phase_numbers) {
          print(paste0("System: ", system_id, ", ", phase, " phase phase_number: ", phase_number))

          # Filter data for the current phase and system
          data_system_phase <- data_system %>%
            filter(ConsecActive == ifelse(phase == "Active", phase_number, 0)) %>%
            filter(ConsecInactive == ifelse(phase == "Inactive", phase_number, 0)) %>%
            as_tibble()
          
          # Skip if no data exists for this phase
          if (nrow(data_system_phase) == 0) {
            message(paste0("No data found for System: ", system_id, ", ", phase, " phase number: ", phase_number, ". Skipping..."))
            next
          }

          # Initialize lists for storing animal positions and temporary data
          animal_list <- list(
            "animal_1" = list(name = "", time = "", position = 0),
            "animal_2" = list(name = "", time = "", position = 0),
            "animal_3" = list(name = "", time = "", position = 0),
            "animal_4" = list(name = "", time = "", position = 0),
            "data_temp" = list(elapsed_seconds = 0, current_row = 0)
          )

          # Initialize lists storing proximity, position, and movement results
          count_proximity_list <- list(m1 = c(0, 0, 0, 0),
                                       m2 = c(0, 0, 0, 0),
                                       m3 = c(0, 0, 0, 0),
                                       m4 = c(0, 0, 0, 0))

          count_position_list <- list(c(1, 0), c(2, 0), c(3, 0), c(4, 0),
                                      c(5, 0), c(6, 0), c(7, 0), c(8, 0))

          total_proximity_list <- list(c(animal_ids[1], 0), c(animal_ids[2], 0),
                                       c(animal_ids[3], 0), c(animal_ids[4], 0))

          count_movement_list <- list(c(animal_ids[1], 0), c(animal_ids[2], 0),
                                      c(animal_ids[3], 0), c(animal_ids[4], 0),
                                      c(system_id, 0))

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
          while (current_row != total_rows && current_row < total_rows) {

            # Make a copy of the previous animal list
            previous_animal_list <- animal_list

            # Update animal list with new positions from the next row
            animal_list <- update_animal_list(animal_ids,
                                              animal_list,
                                              data_system_phase,
                                              initial_time,
                                              current_row)

            # Update elapsed_seconds
            elapsed_seconds <- animal_list[["data_temp"]][["elapsed_seconds"]]

            # Update result lists using analysis functions
            if (system_complete) {
              count_proximity_list <- update_proximity(previous_animal_list,
                                                       animal_list,
                                                       count_proximity_list,
                                                       elapsed_seconds)
            }

            count_position_list <- update_position(previous_animal_list,
                                                   animal_list,
                                                   count_position_list,
                                                   elapsed_seconds)

            if (system_complete) {
              total_proximity_list <- update_total_proximity(previous_animal_list,
                                                             animal_list,
                                                             total_proximity_list,
                                                             elapsed_seconds)
            }

            count_movement_list <- update_movement(previous_animal_list,
                                                   animal_list,
                                                   count_movement_list,
                                                   elapsed_seconds)

            # Update current_row and initial_time for the next iteration
            current_row <- animal_list[["data_temp"]][["current_row"]]
            initial_time <- animal_list[[1]][[2]]
          }

          # Generate heatmaps
          message("Generating heatmaps for proximity and positional data")
          if(system_complete) {
            heatmap_closeness <- generateHeatMapProximity(count_proximity_list,
                                                          batch,
                                                          cageChange,
                                                          system_id,
                                                          animal_ids,
                                                          phase,
                                                          phase_number)
          }

          heatmap_positions <- generateHeatMapPositions(count_position_list,
                                                        batch,
                                                        cageChange,
                                                        system_id,
                                                        phase,
                                                        phase_number)

          # Save heatmaps in specific system list
          if (system_complete) {
            systemHeatmaps_proximity <- c(systemHeatmaps_proximity,
                                          list(heatmap_closeness))
          }

          systemHeatmaps_positions <- c(systemHeatmaps_positions,
                                        list(heatmap_positions))

          # Record total proximity data in the result tibble
          if (system_complete) {
            message("Recording total proximity data in tibble")
            for (i in 1:4) {
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

      # ===================================================================
      # ANALYSIS BY HALF-HOUR PERIODS
      # ===================================================================
      if (analyze_by_halfhour) {
        message("Starting half-hour period analysis")
        
        for (halfhour_period in halfhour_periods) {
          print(paste0("System: ", system_id, ", Half-hour period: ", halfhour_period))
          
          # Filter data for the current half-hour period
          data_system_halfhour <- data_system %>%
            filter(HalfHoursElapsed == halfhour_period) %>%
            as_tibble()
          
          # Skip if no data exists for this period
          if (nrow(data_system_halfhour) == 0) {
            message(paste0("No data found for System: ", system_id, ", Half-hour period: ", halfhour_period, ". Skipping..."))
            next
          }
          
          # Initialize lists for storing animal positions and temporary data
          animal_list <- list(
            "animal_1" = list(name = "", time = "", position = 0),
            "animal_2" = list(name = "", time = "", position = 0),
            "animal_3" = list(name = "", time = "", position = 0),
            "animal_4" = list(name = "", time = "", position = 0),
            "data_temp" = list(elapsed_seconds = 0, current_row = 0)
          )
          
          # Initialize lists storing proximity and movement results
          total_proximity_list <- list(c(animal_ids[1], 0), c(animal_ids[2], 0),
                                       c(animal_ids[3], 0), c(animal_ids[4], 0))
          
          count_movement_list <- list(c(animal_ids[1], 0), c(animal_ids[2], 0),
                                      c(animal_ids[3], 0), c(animal_ids[4], 0),
                                      c(system_id, 0))
          
          count_proximity_list <- list(m1 = c(0, 0, 0, 0),
                                       m2 = c(0, 0, 0, 0),
                                       m3 = c(0, 0, 0, 0),
                                       m4 = c(0, 0, 0, 0))
          
          # Perform calculations for the current half-hour period
          message("Calculating proximity and movement data for half-hour period")
          
          # Initialize animal positions
          animal_list <- initialize_animal_positions(animal_ids, data_system_halfhour, animal_list)
          
          # Assign initial time and line number for the while loop
          initial_time <- animal_list[[1]][[2]]
          current_row <- 5
          elapsed_seconds <- 0
          total_rows <- nrow(data_system_halfhour) + 1
          
          # Iterate through the rows of the current half-hour period data
          while (current_row != total_rows && current_row < total_rows) {
            
            # Make a copy of the previous animal list
            previous_animal_list <- animal_list
            
            # Update animal list with new positions from the next row
            animal_list <- update_animal_list(animal_ids,
                                              animal_list,
                                              data_system_halfhour,
                                              initial_time,
                                              current_row)
            
            # Update elapsed_seconds
            elapsed_seconds <- animal_list[["data_temp"]][["elapsed_seconds"]]
            
            # Update result lists using analysis functions
            if (system_complete) {
              count_proximity_list <- update_proximity(previous_animal_list,
                                                       animal_list,
                                                       count_proximity_list,
                                                       elapsed_seconds)
              
              total_proximity_list <- update_total_proximity(previous_animal_list,
                                                             animal_list,
                                                             total_proximity_list,
                                                             elapsed_seconds)
            }
            
            count_movement_list <- update_movement(previous_animal_list,
                                                   animal_list,
                                                   count_movement_list,
                                                   elapsed_seconds)
            
            # Update current_row and initial_time for the next iteration
            current_row <- animal_list[["data_temp"]][["current_row"]]
            initial_time <- animal_list[[1]][[2]]
          }
          
          # Record half-hour proximity data in the result tibble
          if (system_complete) {
            message("Recording half-hour proximity data in tibble")
            for (i in 1:4) {
              animal <- total_proximity_list[[i]][[1]]
              h <- paste0("H", halfhour_period)
              row <- which(result_halfhour_proximity$HalfHour == h)
              result_halfhour_proximity[[animal]][[row]] <- total_proximity_list[[i]][[2]]
            }
          }
          
          # Record half-hour movement data in the result tibble
          message("Recording half-hour movement data in tibble")
          h <- paste0("H", halfhour_period)
          row <- which(result_halfhour_movement$HalfHour == h)
          for(i in 1:4) {
            animal <- count_movement_list[[i]][[1]]
            if(!is.na(animal)) {
              result_halfhour_movement[[animal]][[row]] <- count_movement_list[[i]][[2]]
            }
          }
          result_halfhour_movement[[system_id]][[row]] <- count_movement_list[[5]][[2]]
        }
      }

      # Generate total proximity graph
      if (system_complete) {
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
    # Save Tables and Plots
    # ---------------------------------------------------
    
    if (save_tables == TRUE || save_plots == TRUE) {
      message("=======================================================================")
      message("SAVING RESULTS FOR ", batch, " ", cageChange)
      message("=======================================================================")
    }
    
    # ===================================================================
    # 1. SAVE METADATA
    # ===================================================================
    if (save_tables == TRUE) {
      message("1. Saving metadata tables...")
      
      # Save System Animal ID tables
      write.csv(system_animal_ids,
                file = paste0(saving_directory, tables_dir_metadata, "/", batch, "_", cageChange, "_animal_ids.csv"), 
                row.names = FALSE)
      message("   ✓ Saved: ", batch, "_", cageChange, "_animal_ids.csv")
      message("      Location: ", tables_dir_metadata)
    }
    
    # ===================================================================
    # 2. SAVE PHASE-BASED ANALYSIS (Active/Inactive)
    # ===================================================================
    if (save_tables == TRUE) {
      message("2. Saving phase-based analysis tables...")
      
      # Save Social Proximity tables (by phase)
      write.csv(result_total_proximity,
                file = paste0(saving_directory, tables_dir_phase, "/", batch, "_", cageChange, "_total_proximity.csv"), 
                row.names = FALSE)
      message("   ✓ Saved: ", batch, "_", cageChange, "_total_proximity.csv")
      
      # Save Movement tables (by phase)
      write.csv(result_total_movement,
                file = paste0(saving_directory, tables_dir_phase, "/", batch, "_", cageChange, "_total_movement.csv"), 
                row.names = FALSE)
      message("   ✓ Saved: ", batch, "_", cageChange, "_total_movement.csv")
      message("      Location: ", tables_dir_phase)
    }
    
    if (save_plots == TRUE) {
      message("3. Saving phase-based analysis plots...")
      
      # Save total proximity plots (one per complete system)
      if (length(all_plots_total_proximity) > 0) {
        
        # Apply theme and save each plot
        for (i in seq_along(all_plots_total_proximity)) {
          # Extract plot information for title
          original_plot <- all_plots_total_proximity[[i]]
          plot_title <- if(!is.null(original_plot$labels$title)) original_plot$labels$title else paste("Total Proximity -", batch, cageChange)
          plot_subtitle <- if(!is.null(original_plot$labels$subtitle)) original_plot$labels$subtitle else paste("System", i)
          
          # Apply ultra-modern minimalist theme
          proximity_lineplot <- generate_proximity_lineplot(original_plot, plot_title, plot_subtitle, batch, cageChange)
          
          ggsave(
        filename = paste0(saving_directory, plots_dir_phase, "/total_proximity_", batch, "_", cageChange, "_sys.", i, ".svg"), 
        plot = proximity_lineplot,
        width = 8,
        height = 8,
        dpi = 300
          )
        }
        message("   ✓ Saved ", length(all_plots_total_proximity), " total proximity plots with ultra-modern minimal design")
        
        # Save proximity heatmaps (one per complete system, multiple phases per system)
        if (length(allHeatmaps_proximity) > 0) {
        saved_count <- 0
        for (i in seq_along(allHeatmaps_proximity)) {
          
          # Check if heatmap list exists and has elements
          if (length(allHeatmaps_proximity[[i]]) == 0) next
          
          # Check if the first heatmap exists and has a title
          if (is.null(allHeatmaps_proximity[[i]][[1]]) || 
            is.null(allHeatmaps_proximity[[i]][[1]]$labels) ||
            is.null(allHeatmaps_proximity[[i]][[1]]$labels$title)) next
          
          # Extract system name from title
          title <- allHeatmaps_proximity[[i]][[1]]$labels$title
          pattern <- "sys\\.[0-9]"
          
          match_result <- regexec(pattern, title)
          if (match_result[[1]][1] > 0) {
          system_substring <- regmatches(title, match_result)[[1]]
          } else {
          system_substring <- paste0("system_", i)
          }
          
          # Separate Active and Inactive phase heatmaps by checking subtitle
          active_plots <- list()
          inactive_plots <- list()
          
          for (plot in allHeatmaps_proximity[[i]]) {
          if (!is.null(plot$labels$subtitle)) {
            subtitle <- as.character(plot$labels$subtitle)
            
            # Debug: print subtitle to see format
            # print(paste("Proximity subtitle:", subtitle))
            
            # Check if it starts with "Phase: Active" or "Phase: Inactive"
            if (grepl("^Phase:\\s*Active", subtitle, ignore.case = TRUE)) {
            active_plots <- c(active_plots, list(plot))
            } else if (grepl("^Phase:\\s*Inactive", subtitle, ignore.case = TRUE)) {
            inactive_plots <- c(inactive_plots, list(plot))
            } else {
            # Fallback: if subtitle contains the phase number pattern
            message("   Warning: Could not classify proximity plot with subtitle: ", subtitle)
            }
          }
          }
          
          n_active <- length(active_plots)
          n_inactive <- length(inactive_plots)
          n_rows <- max(n_active, n_inactive)
          
          message("   Arranging ", n_active, " active and ", n_inactive, " inactive proximity heatmaps for ", system_substring)
          
          # If no plots were classified, fall back to original list
          if (n_active == 0 && n_inactive == 0) {
          message("   Warning: No proximity plots classified by phase type. Using sequential layout.")
          combined_plots <- allHeatmaps_proximity[[i]]
          n_plots <- length(allHeatmaps_proximity[[i]])
          n_rows <- ceiling(n_plots / 2)
          layout_mat <- matrix(1:(n_rows * 2), ncol = 2, byrow = TRUE)
          if (n_plots %% 2 == 1) layout_mat[n_rows, 2] <- NA
          } else {
          # Create layout matrix: Active on left, Inactive on right
          layout_mat <- matrix(NA, nrow = n_rows, ncol = 2)
          
          # Fill left column with active plots (1, 2, 3, 4...)
          if (n_active > 0) {
            for (j in 1:n_active) {
            layout_mat[j, 1] <- j
            }
          }
          
          # Fill right column with inactive plots (n_active+1, n_active+2, ...)
          if (n_inactive > 0) {
            for (j in 1:n_inactive) {
            layout_mat[j, 2] <- n_active + j
            }
          }
          
          # Combine plots in order: all active first, then all inactive
          combined_plots <- c(active_plots, inactive_plots)
          }
          
          # Set plot dimensions
          plot_width <- 14
          plot_height <- 5 * n_rows
          
          # 
          # Save proximity heatmap grid as SVG
          tryCatch({
          ggsave(
            filename = paste0(saving_directory, plots_dir_phase, "/allHeatmaps_proximity_", batch, "_", cageChange, "_", system_substring, ".svg"),
            plot = gridExtra::arrangeGrob(grobs = combined_plots,
                          ncol = 2,
                          layout_matrix = layout_mat),
            width = plot_width,
            height = plot_height
          )
          saved_count <- saved_count + 1
          }, error = function(e) {
          message("   ✗ Error saving proximity heatmap for ", system_substring, ": ", e$message)
          })
        }
        message("   ✓ Saved ", saved_count, " proximity heatmap grids as SVG")
        }
        
        # Save positional heatmaps (one per system, multiple phases per system)
        if (length(allHeatmaps_positions) > 0 && 
          !is.null(allHeatmaps_positions[[batch]]) && 
          !is.null(allHeatmaps_positions[[batch]][[cageChange]])) {
        
        saved_count <- 0
        
        for (sys_name in names(allHeatmaps_positions[[batch]][[cageChange]])) {
          heatmap_list <- allHeatmaps_positions[[batch]][[cageChange]][[sys_name]]
          
          # Check if heatmap list exists and has elements
          if (is.null(heatmap_list) || length(heatmap_list) == 0) {
          message("   Skipping empty position heatmap for ", sys_name)
          next
          }
          
          # Check if the first heatmap exists and has a title
          if (is.null(heatmap_list[[1]]) || 
            is.null(heatmap_list[[1]]$labels) ||
            is.null(heatmap_list[[1]]$labels$title)) {
          message("   Skipping position heatmap ", sys_name, " - no valid title found")
          next
          }
          
          # Separate Active and Inactive phase heatmaps by checking subtitle
          active_plots <- list()
          inactive_plots <- list()
          
          for (plot in heatmap_list) {
          if (!is.null(plot$labels$subtitle)) {
            subtitle <- as.character(plot$labels$subtitle)
            
            # Debug: print subtitle to see format
            # print(paste("Subtitle:", subtitle))
            
            # Check if it starts with "Phase: Active" or "Phase: Inactive"
            if (grepl("^Phase:\\s*Active", subtitle, ignore.case = TRUE)) {
            active_plots <- c(active_plots, list(plot))
            } else if (grepl("^Phase:\\s*Inactive", subtitle, ignore.case = TRUE)) {
            inactive_plots <- c(inactive_plots, list(plot))
            } else {
            # Fallback: if subtitle contains the phase number pattern
            message("   Warning: Could not classify plot with subtitle: ", subtitle)
            }
          }
          }
          
          n_active <- length(active_plots)
          n_inactive <- length(inactive_plots)
          n_rows <- max(n_active, n_inactive)
          
          message("   Arranging ", n_active, " active and ", n_inactive, " inactive position heatmaps for ", sys_name)
          
          # If no plots were classified, fall back to original list
          if (n_active == 0 && n_inactive == 0) {
          message("   Warning: No plots classified by phase type. Using sequential layout.")
          combined_plots <- heatmap_list
          n_plots <- length(heatmap_list)
          n_rows <- ceiling(n_plots / 2)
          layout_mat <- matrix(1:(n_rows * 2), ncol = 2, byrow = TRUE)
          if (n_plots %% 2 == 1) layout_mat[n_rows, 2] <- NA
          } else {
          # Create layout matrix: Active on left, Inactive on right
          layout_mat <- matrix(NA, nrow = n_rows, ncol = 2)
          
          # Fill left column with active plots (1, 2, 3, 4...)
          if (n_active > 0) {
            for (i in 1:n_active) {
            layout_mat[i, 1] <- i
            }
          }
          
          # Fill right column with inactive plots (n_active+1, n_active+2, ...)
          if (n_inactive > 0) {
            for (i in 1:n_inactive) {
            layout_mat[i, 2] <- n_active + i
            }
          }
          
          # Combine plots in order: all active first, then all inactive
          combined_plots <- c(active_plots, inactive_plots)
          }
          
          # Set plot dimensions
          plot_width <- 14
          plot_height <- 5 * n_rows
          
          # Save positional heatmap grid as SVG
          tryCatch({
          ggsave(
            filename = paste0(saving_directory, plots_dir_phase, "/allHeatmaps_positions_", batch, "_", cageChange, "_", sys_name, ".svg"),
            plot = gridExtra::arrangeGrob(grobs = combined_plots,
                          ncol = 2,
                          layout_matrix = layout_mat),
            width = plot_width,
            height = plot_height
          )
          saved_count <- saved_count + 1
          }, error = function(e) {
          message("   ✗ Error saving position heatmap for ", sys_name, ": ", e$message)
          })
        }
        
        if (saved_count > 0) {
          message("   ✓ Saved ", saved_count, " position heatmap grids as SVG")
          message("      Location: ", plots_dir_phase)
        }
        }
      }
      
      # ===================================================================
      # 3. SAVE HALF-HOUR ANALYSIS
      # ===================================================================
      if (analyze_by_halfhour && save_tables == TRUE) {
        message("4. Saving half-hour analysis tables...")
        
        # Save Social Proximity tables (by half-hour)
        write.csv(result_halfhour_proximity,
            file = paste0(saving_directory, tables_dir_halfhour, "/", batch, "_", cageChange, "_halfhour_proximity.csv"), 
            row.names = FALSE)
        message("   ✓ Saved: ", batch, "_", cageChange, "_halfhour_proximity.csv")
        
        # Save Movement tables (by half-hour)
        write.csv(result_halfhour_movement,
            file = paste0(saving_directory, tables_dir_halfhour, "/", batch, "_", cageChange, "_halfhour_movement.csv"), 
            row.names = FALSE)
        message("   ✓ Saved: ", batch, "_", cageChange, "_halfhour_movement.csv")
        message("      Location: ", tables_dir_halfhour)
      }
      
      # ===================================================================
      # SUMMARY
      # ===================================================================
      if (save_tables == TRUE || save_plots == TRUE) {
        message("=======================================================================")
        message("FINISHED SAVING ", batch, " ", cageChange)
        message("=======================================================================")
        message("")
      }
      }
    }
  }

# ---------------------------------------------------
# Display Plots in R
# ---------------------------------------------------

if (show_plots == TRUE) {
  message("show Plots")

  # Create a grid of the total closeness plots
  gridExtra::grid.arrange(grobs = all_plots_total_proximity, ncol = 2)

  # Create a grid of the Heatmaps of each system for each analysis
  for (batch in seq_along(allHeatmaps_positions)) {
    for (cc in seq_along(allHeatmaps_positions[[batch]])) {
      for (sys in seq_along(allHeatmaps_positions[[batch]][[cc]])) {
        print(paste(sys))
        if (cc != "CC4") {
          gridExtra::grid.arrange(grobs = allHeatmaps_positions[[batch]][[cc]][[sys]],
                                  ncol = 2,
                                  layout_matrix = rbind(c(1,8),
                                                        c(2,5),
                                                        c(3,6),
                                                        c(4,7)))
        } else {
          gridExtra::grid.arrange(grobs = allHeatmaps_positions[[batch]][[cc]][[sys]],
                                  ncol = 2,
                                  layout_matrix = rbind(c(1,10),
                                                        c(2,6),
                                                        c(3,7),
                                                        c(4,8),
                                                        c(5,9)))
        }
      }
    }
  }
}
