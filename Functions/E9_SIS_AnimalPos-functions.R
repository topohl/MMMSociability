#' @file E9_SIS_AnimalPos-functions.R
#' 
#' @path /c:/Users/topohl/Documents/GitHub/MMMSociability/
#' 
#' @date 01/2025
#' 
#' @author Tobias Pohl, Anja Magister
#' 
#' @title Analysis of Animal Positions - Functions
#' 
#' @description
#' This script contains functions used for the analysis of animal positions
#' within the MMMSociability project. The functions are designed to process,
#' analyze, and visualize positional data of animals in various experimental
#' settings.
#' 
#' @details
#' The functions included in this script are essential for the preprocessing
#' and analysis of positional data. They facilitate tasks such as data
#' cleaning, transformation, and statistical analysis. The script is part of
#' the E9_SIS module within the MMMSociability project repository.
#' 
#' @usage
#' Source this script to access the functions for analyzing animal positions.
#' Ensure that all dependencies and required packages are installed and loaded
#' before using the functions.

#' @title Preprocess a single data file
#'
#' @description 
#' This function preprocesses a single CSV file for animal position data. It reads the file, performs data cleaning, 
#' extracts relevant columns, adds phase and position information, and saves the preprocessed file.
#' 
#' @details
#' The function reads the raw data from the specified file path, cleans the data by removing unnecessary columns,
#' converts the DateTime column to POSIXct format, splits the Animal column into AnimalID and System, adds PositionID
#' based on the xPos and yPos coordinates, and computes phase transitions. It also defines active and inactive phases
#' based on the time of day and calculates consecutive active and inactive phases. The preprocessed data is then saved
#' as a CSV file in the output directory.
#' 
#' @param batch Character. Name of the batch directory containing the file.
#' @param change Character. Change condition associated with the file.
#' @param exclAnimals Character vector. IDs of animals to exclude from the analysis.
#' @param outputDir Character. Path to the output directory for saving the preprocessed data.
#' @param find_id Function. A function to find the PositionID based on coordinates.
#'
#' @return None. The function saves the preprocessed data as a CSV file in the output directory.
#' @export
preprocess_file <- function(batch, change, exclAnimals) {
  
  # Construct the file path for the raw data
  filename <- paste0("E9_SIS_", batch, "_", change, "_AnimalPos")
  csvFilePath <- file.path(getwd(), "raw_data", batch, paste0(filename, ".csv"))
  
  # Check if the file exists; skip if not found
  if (!file.exists(csvFilePath)) {
    warning(paste("File", csvFilePath, "does not exist. Skipping to next file."))
    return(NULL)
  }
  
  # Read the CSV file into a tibble
  data <- as_tibble(read_delim(csvFilePath, delim = ";", show_col_types = FALSE))
  
  # Preprocessing steps
  data <- data %>%
    # Remove unnecessary columns
    select(-c(RFID, AM, zPos)) %>%
    # Convert DateTime to POSIXct format
    mutate(DateTime = as.POSIXct(DateTime, format = "%d.%m.%Y %H:%M:%S", tz = "UTC")) %>%
    # Split the Animal column into AnimalID and System
    separate(Animal, into = c("AnimalID", "System"), sep = "[-_]")
  
  # Define Position Mapping Table
  position_ids <- tibble(
    PositionID = 1:8, 
    xPos = c(0, 100, 200, 300, 0, 100, 200, 300), 
    yPos = c(0, 0, 0, 0, 116, 116, 116, 116)
  )
  
  # Add PositionID based on xPos and yPos
  data <- data %>%
    rowwise() %>%
    mutate(PositionID = find_id(xPos, yPos, position_id = position_ids)) %>%
    # Retain relevant columns and exclude specified animals
    select(DateTime, AnimalID, System, PositionID) %>%
    filter(!AnimalID %in% exclAnimals) %>%
    arrange(DateTime)
  
  # Add phase information and phase_transitions
  data <- data %>%
    # Add phase transition phase_transitions
    compute_phase_transitions() %>%
    # Define active and inactive phases
    mutate(Phase = ifelse(
      format(DateTime, "%H:%M", tz = "UTC") >= "18:30" | 
      format(DateTime, "%H:%M", tz = "UTC") < "06:30", 
      "Active", "Inactive"
    )) %>%
    # Initialize consecutive phase columns
    mutate(ConsecActive = 0, ConsecInactive = 0) %>%
    # Update consecutive phases and add day phase_transitions
    count_phases() %>%
    compute_day_transitions()
  
  # Define the output file path
  outputFileDir <- file.path(outputDir, paste0(filename, "_preprocessed.csv"))
  
  # Save the preprocessed data
  write_csv(data, outputFileDir)
  
  # Log message upon successful save
  message(paste("File saved at", outputFileDir))
}

#' Find Position ID Based on Coordinates
#'
#' This function searches for a corresponding position ID based on a given pair of x and y coordinates 
#' in a lookup table. The coordinates are adjusted to specific predefined ranges before the search.
#'
#' @param x_Pos A numeric value representing the x-coordinate of the cage.
#' @param y_Pos A numeric value representing the y-coordinate of the cage.
#' @param position_id A tibble containing the lookup table with columns `xPos`, `yPos`, and `PositionID`.
#'        The `xPos` and `yPos` columns represent the coordinates, and `PositionID` is the ID associated 
#'        with each coordinate pair.
#' 
#' @return A numeric value representing the `PositionID` corresponding to the given coordinates, 
#'         or `NA` if no match is found.
#' 
#' @details 
#' The function first adjusts the x and y coordinates to predefined ranges before searching the 
#' `position_id` for the matching position. If a match is found, the corresponding `PositionID` is returned. 
#' If no match is found, the function prints the input coordinates and returns `NA`.
#' 
#' The x-coordinate is adjusted as follows:
#' - Values less than 100 are set to 0
#' - Values between 100 and 199 are set to 100
#' - Values between 200 and 299 are set to 200
#' - Values greater than or equal to 300 are set to 300
#'
#' The y-coordinate is adjusted as follows:
#' - Values less than 116 are set to 0
#' - Values greater than or equal to 116 are set to 116
#'
find_id <- function(x_Pos, y_Pos, position_id) {
  
  # Adjust y-coordinate based on predefined ranges
  if(y_Pos < 116) {y_Pos <- 0} 
  else if(y_Pos >= 116) {y_Pos <- 116}
  
  # Adjust x-coordinate based on predefined ranges
  if(x_Pos < 100) {x_Pos <- 0} 
  else if(x_Pos < 200) {x_Pos <- 100} 
  else if(x_Pos < 300) {x_Pos <- 200} 
  else if(x_Pos >= 300) {x_Pos <- 300}
  
  # Search for position in the lookup table
  result <- position_id %>%
    filter(xPos == x_Pos, yPos == y_Pos) %>%
    select(PositionID)
  
  # If a match is found, return the PositionID, otherwise return NA
  if (nrow(result) > 0) {
    return(result$PositionID)
  } else {
    # Print the input coordinates and return NA if no match is found
    cat("xpos: ", x_Pos, " , ypos: ", y_Pos, "\n")
    print("NA!")
    return(NA) # Return NA if no match found
  }
}

#' Compute Phase Transitions in Experimental Data
#'
#' Adds additional rows to the dataset to represent RFID information at the phase transition points (06:30 AM and 06:30 PM) for each animal in each system.
#'
#' @param preprocessed_data A tibble containing the experimental data, which includes columns for `DateTime`, `System`, and `AnimalID`.
#'
#' @return A tibble with additional rows representing phase transitions. Each added row provides a starting value for each animal ID in every system during each phase.
#'
#' @details
#' - This function identifies the unique dates in the dataset and determines phase transition times for each date (06:30 AM and 06:30 PM).
#' - For each transition point, it ensures there is an entry for every animal in each system by duplicating the closest pre-transition row and updating the timestamp to the transition time.
#' - The resulting dataset is sorted by `DateTime` for consistency.
#' - The last 20 rows are removed to account for unnecessary phase transition rows at the end of the data.
#'
#' @param experiment_dates A Date vector containing the dates of the experiment.
#' @param unique_dates A Date vector containing the unique dates in the experiment.
#' @param unique_systems A character vector containing the unique system names in the experiment.
#' @param date A Date representing the current date in the iteration.
#' @param active_phase_end POSIXct. The end time of the active phase.
#' @param inactive_phase_start POSIXct. The start time of the inactive phase.
#' @param inactive_phase_end POSIXct. The end time of the inactive phase.
#' @param active_phase_start POSIXct. The start time of the active phase.
#' @param filtered_time A tibble containing the filtered data based on the phase transition time.
#' @param filtered_system A tibble containing the filtered data based on the system.
#' @param animal_ids A character vector containing the unique animal IDs in the current system.
#' @param filtered_row A tibble containing the filtered row for the specific animal.
#' @param system A character representing the system name for the current iteration.
#' @param animal_id A character representing the animal ID for the current iteration.
#' 
compute_phase_transitions <- function(preprocessed_data) {
  
  # Save dates of the experiment
  experiment_dates <- as.Date(preprocessed_data$DateTime)
  unique_dates <- unique(experiment_dates)
  
  # Save existing systems
  unique_systems <- unique(preprocessed_data$System)
  
  for(i in 1:length(unique_dates)) {
    date <- unique_dates[i]
  
    # Define phase transition points depending on the actual date
    # "GMT" means "UTC"
    active_phase_end <- as.POSIXct(paste(date, "06:29:59"), tz = "GMT")
    inactive_phase_start <- as.POSIXct(paste(date, "06:30:00"), tz = "GMT")
    inactive_phase_end <- as.POSIXct(paste(date, "18:29:59"), tz = "GMT")
    active_phase_start <- as.POSIXct(paste(date, "18:30:00"), tz = "GMT")
    
    for(phase_transition in c(active_phase_end, inactive_phase_start, inactive_phase_end, active_phase_start)) {
   
      # Filter rows which are earlier than the phase transition
      filtered_time <- preprocessed_data %>%
      as_tibble() %>%
      filter(DateTime < phase_transition)
    
      # We need a phase transition for every existing system and for every animal in each system
      for(system in unique_systems) {
    
        # Filter rows from the specific system 
        filtered_system <- filtered_time %>%
          as_tibble() %>%
          filter(System == system)
    
        # Check which animals exist in this system
        animal_ids <- unique(filtered_system$AnimalID)
    
        for(animal_id in animal_ids) {
          # Filter specific animal row from the system with the nearest DateTime to the phase transition but still earlier
          filtered_row <- filtered_system %>%
          as_tibble() %>%
          filter(AnimalID == animal_id) %>%
          slice_max(order_by = DateTime)
      
          # If time on this day is recorded
          if(length(filtered_row) != 0) {
            # Use copied row as a new one and change the date to the phase transition date
            filtered_row <- filtered_row %>%
            mutate(DateTime = as.POSIXct(phase_transition, tz = "UTC"))
      
            # Add new row to preprocessed_data
            preprocessed_data <- add_row(preprocessed_data, filtered_row)
          }
        }
      }
    }
  }
  
  # Sort tibble again by DateTime to bring new entries to the correct position
  preprocessed_data <- preprocessed_data %>%
  as_tibble() %>%
  arrange(., DateTime)
  
  # Remove the last 20 rows to eliminate redundant phase transition entries from the final day and late phase transitions
  preprocessed_data <- preprocessed_data %>% filter(row_number() <= n() - 20)
  
  return(preprocessed_data)
}

#' @title Compute Day Transitions in Experimental Data
#' 
#' @description
#' This function adds additional rows to the dataset to represent RFID information at the day transition points (midnight) for each animal in each system.
#' 
#' @param preprocessed_data A tibble containing the experimental data, which includes columns for `DateTime`, `System`, and `AnimalID`.
#' @param filtered_time A tibble containing the filtered data based on the phase transition time.
#' @param filtered_system A tibble containing the filtered data based on the system.
#' @param animal_ids A vector containing the unique animal IDs in the current system.
#' @param filtered_row A tibble containing the filtered row for the specific animal.
#' @param date Date. The date of the experiment.
#' @param phase_transition POSIXct. The phase transition time for the current date.
#' @param unique_dates Date vector. Unique dates in the experiment.
#' @param unique_systems Character vector. Unique system names in the experiment.
#' @param system Character. The system name for the current iteration.
#' @param animal_id Character. The animal ID for the current iteration.
#' 
compute_day_transitions <- function(preprocessed_data) {
  # Save dates of the experiment
  experiment_dates <- as.Date(preprocessed_data$DateTime)
  unique_dates <- unique(experiment_dates)
  
  # Save existing systems
  unique_systems <- unique(preprocessed_data$System)
  
  for(i in 1:length(unique_dates)) {
    date <- unique_dates[i]
    
    # Define phase transition point based on the actual date
    # "GMT" refers to "UTC"
    phase_transition <- as.POSIXct(paste(date, "00:00:00"), tz = "GMT")  #start of new day
    
    # Filter rows that are earlier than the phase transition
    filtered_time <- preprocessed_data %>%
      as_tibble() %>%
      filter(DateTime < phase_transition)
    
    for(system in unique_systems) {
      
      # Filter rows from the specific system 
      filtered_system <- filtered_time %>%
        as_tibble() %>%
        filter(System == system)
      
      # Identify which animals are present in the current system
      animal_ids <- unique(filtered_system$AnimalID)
      
      for(animal_id in animal_ids) {
        # Filter the specific animal row from the system with the nearest DateTime to the phase transition but still earlier
        filtered_row <- filtered_system %>%
          as_tibble() %>%
          filter(AnimalID == animal_id) %>%
          slice_max(order_by = DateTime)
        
        # If the time on this day is recorded
        if(length(filtered_row) != 0) {
          # Duplicate the row and update the DateTime to the phase transition time
          filtered_row <- filtered_row %>%
            mutate(DateTime = as.POSIXct(phase_transition, tz = "UTC"))
          
          # Append the new row to the preprocessed_data tibble
          preprocessed_data <- add_row(preprocessed_data, filtered_row)
        }
      }
    }
  }

  # Sort the tibble by DateTime to ensure new entries are correctly positioned
  preprocessed_data <- preprocessed_data %>%
    as_tibble() %>%
    arrange(., DateTime)
  
  return(preprocessed_data)
}

#' @title Count Phases in Experimental Data
#' 
#' @description
#' This function counts the number of consecutive active and inactive phases in the experimental data.
#' It updates the `ConsecActive` and `ConsecInactive` columns based on the phase transitions in the data.
#' 
#' @details
#' The function iterates over the data and counts the number of consecutive active and inactive phases.
#' It updates the `ConsecActive` and `ConsecInactive` columns based on the phase transitions in the data.
#' The function ensures that the `ConsecActive` and `ConsecInactive` columns are correctly updated for each row.
#' 
#' @param data A tibble containing the experimental data, which includes columns for `DateTime`, `Phase`, `ConsecActive`, and `ConsecInactive`.
#' @param current_row A tibble containing the current row of the data.
#' @param previous_row A tibble containing the previous row of the data.
#' @param active_phases Numeric. The count of active phases in the data.
#' @param inactive_phases Numeric. The count of inactive phases in the data.
#' @param i Integer. The current iteration index.
#' 
count_phases <- function(data) {
  
  #initialize Phase counter
  active_phases <- 0
  inactive_phases <- 0
  
  for (i in 1:nrow(data)) {
    
    current_row <- data[i,]
    
    if (i == 1) {
      
      if(current_row$Phase == "Active") {
        current_row$ConsecActive <- 1
        current_row$ConsecInactive <- 0
        active_phases <- active_phases + 1
      } else {
        current_row$ConsecActive <- 0
        current_row$ConsecInactive <- 1
        inactive_phases <- inactive_phases + 1
      }
    } else {#i>1
      
      #previous row from data
      previous_row <- data[(i - 1),]
      
      if (current_row$Phase != previous_row$Phase) {
        if (current_row$Phase == "Active") {
          active_phases <- active_phases + 1  # change means the counter has to set higher
          current_row$ConsecActive <- active_phases
          current_row$ConsecInactive <-  0
        } else {    #current_row$Phase == "Inactive"
          inactive_phases <- inactive_phases + 1
          current_row$ConsecActive <- 0
          current_row$ConsecInactive <- inactive_phases
        }
      } else {  #current_row$Phase == previous_row$Phase
        current_row$ConsecActive <- previous_row$ConsecActive      #values stay the same number as before
        current_row$ConsecInactive <- previous_row$ConsecInactive
      }
      
      #write new information back into data
      data[(i - 1),] <- previous_row
    }
    
    #write new information back into data
    data[i,] <- current_row  
  }
  
  return(data)
}

#' @title Initialize Animal Positions
#' 
#' @description
#' This function initializes the positions of animals in the system based on the first recorded positions in the data.
#' It extracts the initial time and position of each animal from the data and stores this information in a list.
#' If the initial times of the animals differ, the function sets all initial times to the latest shared time to ensure
#' that all animals start at the same time point for proximity tracking.
#' 
#' @details
#' The function iterates over the list of animal IDs in the system and searches for the first recorded entry for each animal.
#' It extracts the initial time and position of each animal and stores this information in a list. If the initial times of
#' the animals differ, the function sets all initial times to the latest shared time. This ensures that all animals start at
#' the same time point for proximity tracking.
#' 
#' @param system_animal_ids A list containing the IDs of animals in the current system.
#' @param data A tibble containing the positional data for the current system.
#' @param animal_list A list storing the animal IDs, initial time, and initial position for each animal.
#' @param times_vec A numeric vector to store the initial times of each animal for comparison.
#' @param i An integer representing the current iteration index.
#' @param system_animal_id Character. The ID of the current animal in the system.
#' @param first_entry A tibble containing the first recorded entry for the current animal.
#' @param start_time POSIXct. The initial time of the current animal.
#' @param start_position Numeric. The initial position of the current animal.
#' @param end_time POSIXct. The shared end time for all animals in the system.
#' 
initialize_animal_positions <- function(system_animal_ids, data, animal_list) {
  
  # Initialize an empty vector with four variables.
  times_vec <- rep(NA, times = 4)
  
  print(system_animal_ids)
  for (i in 1:length(system_animal_ids)) { #i=1-4
    
    # This section defines the current animal identification name (animal_id_name), 
    # the initial position (start_position), and the initial timestamp (start_time).
    system_animal_id <- system_animal_ids[[i]]
    if(is.na(system_animal_id)) {
      system_animal_id <- paste0("lost_", i)
      start_time <- 0        #time value should be adapted at the end of this function
      start_position <- (-1)  #add non existing position for non existing animal_id
    } else {
      # Search for the first entry in the entire dataset.
      first_entry <- data %>%
        filter(AnimalID == system_animal_id) %>%
        slice(1) #first row

      start_time <- first_entry$DateTime
      start_position <- first_entry$PositionID
    }
    
    # Record the name, position, and timestamp of each animal into the animal_list.
    animal_list[i][[1]] <- system_animal_id
    animal_list[[i]][[3]] <- start_position
    animal_list[[i]][[2]] <- start_time
    
    # Record the current time into a vector for subsequent comparisons.
    times_vec[i] <- start_time
  }
  
  # Check if all initial times are identical
  # If not, update all initial times to the latest shared time
  if(length(unique(times_vec)) != 1) {
    end_time <- max(times_vec)
    #print(end_time)
    for(i in 1:4) {
      animal_list[[i]][[2]] <- end_time
    }
  }
  
  return(animal_list)
}

#' @title Update Social Proximity Tracking in Animal Behavior Analysis
#' 
#' @description
#' This function tracks the proximity of animals over time by comparing their spatial positions. 
#' It updates a result list that records the number of seconds each pair of animals spent in social contact.
#' 
#' @param previous_animal_positions A list containing the recorded positions of animals from the previous time step 
#'        until the second before the current update. Each entry is assumed to be a sublist where the third element 
#'        represents the numerical position identifier.
#' @param current_animal_positions A list containing the updated positions of animals at the current time step. 
#'        The structure is identical to `previous_animal_positions`, where the third element corresponds to the 
#'        spatial position.
#' @param count_proximity_list A nested list tracking the number of seconds each pair of animals spent in proximity.
#'        Each entry represents a matrix-like structure where the row and column indices correspond to animal IDs, 
#'        and the values store the accumulated social contact time in seconds.
#' @param elapsed_seconds A numeric value representing the time interval (in seconds) that has passed since the last update.
#' 
#' @return
#' A nested list (`count_proximity_list`) updated with new proximity durations for each pair of animals.
#' 
#' @details
#' The function performs pairwise comparisons of animal positions to track proximity:
#' - It first compares `previous_animal_positions` to determine proximity duration over the past `(elapsed_seconds - 1)` seconds.
#' - It then compares `current_animal_positions` to update proximity for the most recent second.
#' - If two animals share the same position (excluding `-1` values), their contact time is incremented in `count_proximity_list`.
#' - The function ensures symmetry in the `count_proximity_list` matrix, where `count_proximity_list[[i]][[j]]` 
#'   and `count_proximity_list[[j]][[i]]` are updated simultaneously when different individuals are in contact.
#'
update_proximity <- function(previous_animal_positions, current_animal_positions, count_proximity_list, elapsed_seconds) {
  
  # Compare the positions of each pair of animals (third element in the list represents the position)
  # If the positions are the same, update the proximity count in count_proximity_list
  for (i in 1:4) {
    for (j in i:4) {

      # Comparison of previous_animal_positions for the last (elapsed_seconds - 1) seconds
      if (previous_animal_positions[[i]][[3]] == previous_animal_positions[[j]][[3]]) {
        count_proximity_list[[i]][[j]] <- count_proximity_list[[i]][[j]] + (elapsed_seconds - 1)
        # Enter new contact second in proximity list if the animals are at the same position and are two different individuals
        if (j != i) {
          count_proximity_list[[j]][[i]] <- count_proximity_list[[j]][[i]] + (elapsed_seconds - 1)
        }
      }

      # Comparison of current_animal_positions for the first new second
      if (current_animal_positions[[i]][[3]] == current_animal_positions[[j]][[3]]) {
        count_proximity_list[[i]][[j]] <- count_proximity_list[[i]][[j]] + 1
        # Enter new contact second in proximity list if the animals are at the same position and are two different individuals
        if (j != i) {
          count_proximity_list[[j]][[i]] <- count_proximity_list[[j]][[i]] + 1
        }
      }
    }
  }
  
  # return updated list of animals that are close to each other
  return(count_proximity_list)
}

#' @title Update Position Occupancy in Spatial Tracking
#' 
#' @description
#' This function updates the occupancy duration of spatial positions based on the movement of tracked animals. 
#' It ensures that each position is counted only once per second, even if occupied by multiple animals simultaneously.
#' 
#' @param previous_animal_positions A list containing the previous positions of tracked animals. 
#'        Each entry is expected to be a sublist where the third element corresponds to the numerical position identifier.
#' @param current_animal_positions A list containing the current positions of tracked animals.
#'        The structure is identical to `previous_animal_positions`, where the third element represents the position identifier.
#' @param count_position_list A list tracking the occupancy duration of each position.
#'        Each position index contains a sublist where the second element stores the cumulative occupancy time.
#' @param elapsed_seconds A numeric value representing the time interval (in seconds) that has passed since the last update.
#' 
#' @return
#' A list (`count_position_list`) with updated occupancy durations for each spatial position. 
#' The second element of each sublist is incremented based on the time an animal occupied the position.
#' 
#' @details
#' This function iterates over a fixed number of tracked animals (assumed to be four) and updates 
#' the duration for which each position was occupied. The function differentiates between:
#' - `previous_animal_positions`: The positions occupied before the update.
#' - `current_animal_positions`: The positions occupied at the time of the update.
#' - If a position is occupied in `previous_animal_positions` but not already recorded for this time period,
#'   its duration is increased by `(elapsed_seconds - 1)`.
#' - If a position is occupied in `current_animal_positions` and not previously recorded in this time step,
#'   its duration is incremented by 1 second.
#' - Positions marked as `-1` are ignored, as they do not correspond to valid spatial locations.
#' 
update_position <- function(previous_animal_positions, current_animal_positions, count_position_list, elapsed_seconds) {
  
  # Vectors to track unique positions occupied by animals within the given time period
  previous_positions <- c()
  current_positions <- c()
  
  # Iterate through the list of tracked animals (assumed to be four in total)
  for (i in 1:4) {
    previous_position <- as.numeric(previous_animal_positions[[i]][[3]])
    current_position <- as.numeric(current_animal_positions[[i]][[3]])

    # Ensure each occupied position is counted only once per second,
    # even if multiple animals are present at the same location.
    # This prevents redundant increments for the same position.
    
    # Process the previous position if it has not already been added and is valid (i.e., not -1)
    if (!previous_position %in% previous_positions & previous_position != (-1)) {
      
      # Update the duration the position was occupied, adjusting for elapsed time
      count_position_list[[previous_position]][[2]] <- count_position_list[[previous_position]][[2]] + (elapsed_seconds - 1)
      
      # Store the position in the tracking vector to prevent double counting
      previous_positions <- append(previous_positions, previous_position)
    }

    # Process the current position if it has not already been added and is valid (i.e., not -1)
    if (!current_position %in% current_positions & current_position != (-1)) {
      
      # Increment the occupancy duration for the current position
      count_position_list[[current_position]][[2]] <- count_position_list[[current_position]][[2]] + 1
      
      # Store the position in the tracking vector to prevent double counting
      current_positions <- append(current_positions, current_position)
    }
  }
  
  # Return the updated list tracking the duration each position was occupied
  return(count_position_list)
}

#' Update Total Proximity for Animals
#'
#' This function calculates the total amount of time spent in closer contact between animals. It compares the positions of
#' each animal at the previous and current timepoints and updates the total proximity list accordingly.
#'
#' @param previous_animal_positions A list of previous positions for each animal, including their IDs and positions.
#' @param current_animal_positions A list of current positions for each animal, including their IDs and positions.
#' @param total_proximity_list A list that tracks the total time spent in proximity for each animal.
#' @param elapsed_seconds The time difference (in seconds) between the previous and current timepoints.
#'
#' @return The updated `total_proximity_list` with the total proximity time for each animal.
#' @export
#'
#' @examples
#' # Assuming previous_animal_positions, current_animal_positions, total_proximity_list, and elapsed_seconds are defined:
#' updated_total_proximity <- update_total_proximity(previous_animal_positions, current_animal_positions, total_proximity_list, elapsed_seconds)
update_total_proximity <- function(previous_animal_positions, current_animal_positions, total_proximity_list, elapsed_seconds) {

  # Loop over each animal ID (assuming there are 4 animals)
  for (animal_id in 1:4) {

    # Define the list of cohabiting animals (excluding the current animal)
    cohabiting_animals <- c(1, 2, 3, 4)
    cohabiting_animals <- cohabiting_animals[-animal_id]

    # Get the current animal ID from the previous animal positions
    current_animal_id <- previous_animal_positions[[animal_id]][[1]]

    ## OLD: REMAINING SECONDS (check if animal was on the same position as another animal)
    if (previous_animal_positions[[animal_id]][[3]] == previous_animal_positions[[cohabiting_animals[1]]][[3]] ||
        previous_animal_positions[[animal_id]][[3]] == previous_animal_positions[[cohabiting_animals[2]]][[3]] ||
        previous_animal_positions[[animal_id]][[3]] == previous_animal_positions[[cohabiting_animals[3]]][[3]]) {

      # Add the remaining proximity seconds (since last entry) into the total proximity list for the current animal
      if (total_proximity_list[[animal_id]][[1]] == current_animal_id) {
        total_proximity_list[[animal_id]][[2]] <- as.numeric(total_proximity_list[[animal_id]][[2]]) + (elapsed_seconds - 1)
      } else {
        stop("Animal IDs are not the same in check_total_proximity")
      }
    }

    ## FIRST NEW SECOND (check if animal is on the same position as another animal in the current second)
    if (current_animal_positions[[animal_id]][[3]] == current_animal_positions[[cohabiting_animals[1]]][[3]] ||
        current_animal_positions[[animal_id]][[3]] == current_animal_positions[[cohabiting_animals[2]]][[3]] ||
        current_animal_positions[[animal_id]][[3]] == current_animal_positions[[cohabiting_animals[3]]][[3]]) {

      # Add the proximity second into the total proximity list for the current animal
      if (total_proximity_list[[animal_id]][[1]] == current_animal_id) {
        total_proximity_list[[animal_id]][[2]] <- as.numeric(total_proximity_list[[animal_id]][[2]]) + 1
      } else {
        stop("Animal ID mismatch in check_total_proximity")
      }
    }
  }

  # Return the updated total proximity list
  return(total_proximity_list)
}

#' Track Animal Movement Between Cage Positions
#'
#' This function tracks the movement of animals between cage positions over time, updating a movement counter list.
#' The function identifies whether any animal has moved between positions in a single time step and updates the movement count for each animal and the system as a whole.
#'
#' @param previous_animal_positions List containing the previous positions of animals in the cage.
#' @param current_animal_positions List containing the current positions of animals in the cage.
#' @param count_movement_list A list that stores movement counts for each animal and the system as a whole.
#' @param elapsed_seconds Integer representing the current time step (not directly used in this function, but retained for consistency with broader analysis).
#' @return Updated `count_movement_list` reflecting movements recorded in the current time step.
#' @details
#' - The function only tracks valid movements, ignoring invalid positions denoted by `-1`.
#' - Updates are made to both individual animal movement counts and the total system movement count.
#' - If no movement occurs in the current time step, a message is printed to indicate this.
update_movement <- function(previous_animal_positions, current_animal_positions, count_movement_list, elapsed_seconds) {
  
  # Variable to track if any movement occurred in this time step
  movement_count <- 0
  
  # Iterate through each animal's position data
  for (i in 1:4) {
    
    # Extract previous and current positions for the animal
    previous_position <- as.numeric(previous_animal_positions[[i]][[3]])
    current_position <- as.numeric(current_animal_positions[[i]][[3]])
    
    # Check if the animal has moved, excluding invalid positions (-1)
    if (previous_position != (-1) & current_position != (-1) & previous_position != current_position) {
      
      # Update the movement count for the specific animal
      count_movement_list[[i]][[2]] <- as.numeric(count_movement_list[[i]][[2]]) + 1
      
      # Update the total system movement count
      count_movement_list[[5]][[2]] <- as.numeric(count_movement_list[[5]][[2]]) + 1
      
      # Increment the movement counter
      movement_count <- movement_count + 1
    }
  }
  
  # Print a message if no movement occurred in this time step
  if (movement_count == 0) {
    message("track_movement: No movement detected in this time step.")
  }
  
  # Return the updated movement count list
  return(count_movement_list)
}

#' @title sec_shift
#' 
#' @description
#' This function shifts the time by one second forward, updating the current timepoint.
#' 
#' @param previous_timepoint The previous timepoint that will be shifted forward by one second.
#' @param current_timepoint The current timepoint after shifting one second forward.
#' 
#' @return The updated timepoint after shifting one second forward.
 
sec_shift <- function(previous_timepoint) {
  # Shift the time by one second forward
  current_timepoint <- previous_timepoint %>%
    as.numeric() %>%
    + 1 %>%
    as.character()

  return(current_timepoint)
}

#' Update Animal List with New Time and Position Information
#' 
#' This function updates the `animal_list` with the most recent time and position information for each animal
#' in the dataset. It checks if the current timepoint matches consecutive rows and updates the positions accordingly.
#' The function is similar to `initialize_animal_positions` but handles the dynamic updating of the list.
#' 
#' @param system_animal_ids A vector containing the IDs/names of the mice in the current system (not explicitly used in the function).
#' @param animal_list A list containing information about each animal, including its ID, current time, and position.
#' @param data A tibble or data frame containing the dataset with position and time information.
#' @param time The initial timepoint that is being compared to the new time data in the dataset.
#' @param line The current row index in the `data` that is being processed.
#' 
#' @return The updated `animal_list` with the new time and position information for the animals.
#' @export
#'
#' @examples
#' # Assuming system_animal_ids, animal_list, data, time, and line are defined:
#' updated_animal_list <- update_animal_list(system_animal_ids, animal_list, data, time, line)
update_animal_list <- function(system_animal_ids, animal_list, data, time, line) {
  
  # Extract the current timepoint as numeric value from the data
  current_timepoint <- as.numeric(data[line, "DateTime"])
  
  # Calculate and store the time difference (in seconds) between the current and previous timepoints
  animal_list[["data_temp"]][["elapsed_seconds"]] <- current_timepoint - as.numeric(time)
  
  # Update the time for each animal in the animal list
  for (i in 1:4) {
    animal_list[[i]][[2]] <- current_timepoint
  }
  
  # While the time is the same as the current timepoint, continue processing the data
  while (as.numeric(data[line, "DateTime"]) == current_timepoint) {
    
    # Update the position for the animal based on AnimalID match
    for (i in 1:4) {
      if (animal_list[[i]][[1]] == as.character(data[line, "AnimalID"])) {
        animal_list[[i]][[3]] <- as.numeric(data[line, "PositionID"])
      }
    }
    
    # If we have reached the last line, move to the next line and break the loop
    if (line == nrow(data)) {
      line <- line + 1
      break
    }
    
    # Otherwise, move to the next row in the data
    line <- line + 1
  }
  
  # Update the line number in the temporary data section of animal_list
  animal_list[["data_temp"]][["current_row"]] <- line
  
  # Return the updated animal_list
  return(animal_list)
}

#' Function: compute_rank
#'
#' Compute Ranks Based on Hourly Values for a Specific System, Cage Change, and Batch
#'
#' This function sorts a numeric vector and assigns ranks to a tibble based on the sorted values.
#' The largest value in the vector is assigned rank 1. The rank is then updated in the `rank_tibble`
#' where the values match specific criteria, such as the system, cage change, and batch.
#'
#' @param rank_tibble A tibble containing the data to which ranks will be assigned. This tibble
#'                    must include columns for `System`, `CageChange`, `Batch`, and the column 
#'                    specified in `hours_column_name`.
#' @param vector A numeric vector containing the values that will be ranked. These values correspond
#'               to the `hours_column_name` in the `rank_tibble`.
#' @param system A character string or factor representing the system that will be used to filter
#'               the rows in the tibble.
#' @param cageChange A character string or factor representing the cage change that will be used
#'                   to filter the rows in the tibble.
#' @param batch A character string or factor representing the batch that will be used to filter
#'              the rows in the tibble.
#' @param hours_column_name The name of the column containing the hourly values in the `rank_tibble`
#'                          that will be ranked.
#' @param rank_column_name The name of the column where the rank values will be stored in the tibble.
#' 
#' @return A tibble with the updated `rank_column_name` reflecting the assigned ranks.
#' 
#' @details The function sorts the `vector` in ascending order (smallest value gets rank 1) and assigns
#'          ranks to the corresponding rows in the `rank_tibble`. The ranks are assigned based on a match 
#'          between the values in the `hours_column_name` and the sorted `vector`. The assignment is restricted
#'          by the specified `system`, `cageChange`, and `batch`.
#' @export
compute_rank <- function(rank_tibble, vector, system, cageChange, batch, hours_column_name, rank_column_name) {
  
  # Sort the vector in ascending order so the smallest value gets rank 1
  vector <- sort(vector)
  
  # Assign ranks based on the sorted vector
  for (i in 1:length(vector)) {
    rank_tibble <- rank_tibble %>%
      mutate({{rank_column_name}} := ifelse(
        (Batch == batch) &
        (CageChange == cageChange) &
        (System == system) &
        (.data[[hours_column_name]] == vector[i]), i, .data[[rank_column_name]]
      )) 
  }
  
  # Return the tibble with the updated rank column
  return(rank_tibble)
}

#' Translate Rank Vector into Score Vector
#'
#' This function converts a vector of ranks into a corresponding vector of scores. 
#' The scoring system is as follows:
#' - Rank 1 = Score 4
#' - Rank 2 = Score 3
#' - Rank 3 = Score 2
#' - Rank 4 = Score 1
#' 
#' @param rank_vec A numeric vector containing ranks (values 1 to 4).
#' @return A numeric vector of scores corresponding to the input ranks.
#' @details If a rank outside the range 1 to 4 is encountered, it will be assigned `NA` in the output.
#' @examples
#' # Example usage:
#' rank_vec <- c(1, 2, 3, 4, 1, 3)
#' score_vec <- translate_rank_in_score_vec(rank_vec)
#' print(score_vec) # Output: c(4, 3, 2, 1, 4, 2)
#' @export
convert_rank_to_score <- function(rank_vec) {
  # Iterate over each rank in the rank_vec to convert it into a corresponding score
  for (i in 1:length(rank_vec)) {
    # Map rank to score using conditional statements
    score <- ifelse(rank_vec[i] == 1, 4, 
                    ifelse(rank_vec[i] == 2, 3, 
                    ifelse(rank_vec[i] == 3, 2, 
                    ifelse(rank_vec[i] == 4, 1, NA))))
    
    # Replace rank with the corresponding score
    rank_vec[i] <- score
  }
  
  score_vec <- rank_vec
  return(score_vec)
}  

#================================================
# FUNCTIONS FOR SHANNON ENTROPY CALCULATION
#================================================

#' Update Cage Position Probability
#'
#' Updates the probability distribution of animal positions in a cage based on previous and current position data. 
#' This function calculates the positional distribution of animals for each frame and updates the 
#' cumulative probability distribution for all positions over a specified time window.
#'
#' @param previous_animal_positions List containing the positions of animals from the previous frame. 
#'   Each sublist contains:
#'   \describe{
#'     \item{[[1]]}{Animal ID (character or numeric).}
#'     \item{[[2]]}{Additional metadata (not used in this function).}
#'     \item{[[3]]}{Animal's position in the previous frame (numeric).}
#'   }
#' @param current_animal_positions List containing the positions of animals from the current frame. 
#'   The structure is the same as `previous_animal_positions`.
#' @param cage_position_probability List to store cumulative probability and time spent at each cage position. 
#'   Each sublist contains:
#'   \describe{
#'     \item{[[1]]}{Position index (integer).}
#'     \item{[[2]]}{Cumulative time spent at this position (in seconds).}
#'     \item{[[3]]}{Cumulative probability for this position (numeric).}
#'   }
#' @param elapsed_seconds Integer. Time window (in seconds) to consider when updating the probability distribution.
#'
#' @return List. Updated `cage_position_probability` with cumulative probabilities and time spent at each position.
#'
#' @details 
#' The function calculates the proportion of animals at each position for both the previous and current frames. 
#' These proportions are used to update the cumulative probability and time spent at each position over the 
#' time window defined by `elapsed_seconds`.
update_cage_position_probability <- function(previous_animal_positions, current_animal_positions, cage_position_probability, elapsed_seconds) {
  
  # Create vectors to store position distributions for 8 positions (initialized to 0)
  previous_position_distribution <- integer(8) # Stores the counts of animals in each position for the previous frame
  current_position_distribution <- integer(8)  # Stores the counts of animals in each position for the current frame
  
  # Calculate the distribution of animals across positions
  for (i in 1:4) {
    # Get the positions from the previous and current frames
    previous_position <- as.numeric(previous_animal_positions[[i]][[3]])  # Extract the position from the previous frame
    current_position <- as.numeric(current_animal_positions[[i]][[3]])    # Extract the position from the current frame
    
    # Increment the respective position counters
    previous_position_distribution[previous_position] <- previous_position_distribution[previous_position] + 1
    current_position_distribution[current_position] <- current_position_distribution[current_position] + 1
  }
  
  # Normalize the position distributions (divide by 4 since there are 4 animals)
  previous_position_distribution <- previous_position_distribution / 4
  current_position_distribution <- current_position_distribution / 4
  
  # Update cumulative probabilities for each of the 8 positions
  for (i in 1:8) {
    # Add the probability for the last second of the previous frame
    cage_position_probability[[i]][[2]] <- cage_position_probability[[i]][[2]] + 1  # Increment the time spent at position i
    cage_position_probability[[i]][[3]] <- cage_position_probability[[i]][[3]] + previous_position_distribution[i]  # Add the proportion of animals
    
    # Add probabilities for the remaining seconds in the time window
    for (second in 1:(elapsed_seconds - 1)) {
      cage_position_probability[[i]][[2]] <- cage_position_probability[[i]][[2]] + 1  # Increment the time spent
      cage_position_probability[[i]][[3]] <- cage_position_probability[[i]][[3]] + previous_position_distribution[i]  # Add the proportion
    }
  }
  
  return(cage_position_probability)
}

#' Update Animal Position Probability
#'
#' Updates the probability distribution of animal positions based on previous and current position data for each animal.
#' The function calculates the number of observed seconds and the cumulative probability percentage for each animal's position.
#'
#' @param previous_animal_positions List containing the positions of animals from the previous frame. 
#'   Each sublist contains:
#'   \describe{
#'     \item{[[1]]}{Animal ID (numeric or character).}
#'     \item{[[2]]}{Additional information (not used in this function).}
#'     \item{[[3]]}{Animal's previous position (numeric).}
#'   }
#' @param current_animal_positions List containing the positions of animals from the current frame. 
#'   The structure is the same as `previous_animal_positions`.
#' @param animal_position_probability Data frame storing cumulative position probabilities and observed seconds for each animal and position. 
#'   It must include the following columns:
#'   \describe{
#'     \item{AnimalID}{Unique identifier for each animal.}
#'     \item{Position}{Position index (numeric).}
#'     \item{SumPercentage}{Cumulative percentage probability for the position (numeric).}
#'     \item{Seconds}{Cumulative seconds the animal has been observed (numeric).}
#'   }
#' @param elapsed_seconds Time window (in seconds) to update the probability distribution for the current frame.
#'
#' @return Updated `animal_position_probability` data frame with updated cumulative percentages and seconds for each animal's position.
#'
#' @details 
#' The function iterates over all animals, updates the observed seconds for their positions, and calculates the new cumulative percentage probability based on their previous and current positions.
#' If an animal is not tracked (indicated by `NA` in `animal_ids`), the function skips that animal.
update_animal_position_probability <- function(previous_animal_positions, current_animal_positions, animal_position_probability, elapsed_seconds) {
  
  # Loop through each animal
  for (i in 1:4) {
    
    # Skip if the animal ID is not tracked (incomplete system)
    if (is.na(animal_ids[i])) { 
      next 
    }
    
    # Define the previous and current positions for the animal
    previous_position <- as.numeric(previous_animal_positions[[i]][[3]])
    current_position <- as.numeric(current_animal_positions[[i]][[3]])
    
    # Find the row corresponding to the animal's previous position
    row_old <- which(
      animal_position_probability$AnimalID == animal_ids[i] & 
      animal_position_probability$Position == previous_position
    )
    
    # Find the row corresponding to the animal's current position
    row_new <- which(
      animal_position_probability$AnimalID == animal_ids[i] & 
      animal_position_probability$Position == current_position
    )
    
    # Update the cumulative percentage for the previous position
    animal_position_probability[["SumPercentage"]][[row_old]] <- 
      animal_position_probability[["SumPercentage"]][[row_old]] + 1
    
    # Update the cumulative percentage for the current position
    animal_position_probability[["SumPercentage"]][[row_new]] <- 
      animal_position_probability[["SumPercentage"]][[row_new]] + (elapsed_seconds - 1)
    
    # Update the observed seconds for the animal across all positions
    animal_position_probability <- animal_position_probability %>%
      mutate(Seconds = ifelse(AnimalID == animal_ids[i], Seconds + elapsed_seconds, Seconds))
  }
  
  return(animal_position_probability)
}

#' Calculate Shannon Entropy
#'
#' Calculates the Shannon entropy of a probability vector representing the distribution of animals across positions in a cage.
#' Shannon entropy is a measure of uncertainty or disorder in a probability distribution, calculated as the negative 
#' sum of the probabilities multiplied by their log base 2.
#'
#' @param prob_vec A numeric vector of length 8 representing the probabilities of animals being at each of 8 positions.
#'   The values in this vector must be between 0 and 1, and the sum of the values should equal 1.
#' 
#' @return Numeric. The Shannon entropy of the probability distribution.
#' 
#' @details 
#' The function iterates through the 8 positions and calculates the Shannon entropy using the formula:
#' 
#' H(X) = - \sum_{i=1}^{n} p(x_i) \log_2(p(x_i))
#'
#' Where \( p(x_i) \) is the probability of being in position \( i \). If a position has a probability of 0, 
#' its contribution to the entropy is ignored, as \( \log_2(0) \) is undefined.
calc_shannon_entropy <- function(prob_vec) {
  
  # Check if the probability vector has exactly 8 elements
  if (length(prob_vec) != 8) {
    print(prob_vec)
    stop("Error in Shannon calculation: prob_vec not the right size")
  }
  
  # Initialize variable to store Shannon entropy
  shannon_entropy <- 0
  
  # Loop through the probability vector and calculate entropy
  for (i in 1:8) {
    # Handle case where the probability is 0 (log2(0) is undefined)
    if (prob_vec[i] == 0) {
      shannon_entropy <- shannon_entropy + 0  # Add 0 for this position
    } else {
      # Add the term for non-zero probabilities
      shannon_entropy <- shannon_entropy + prob_vec[i] * log2(prob_vec[i])
    }
  }
  
  # Negate the sum to calculate entropy
  shannon_entropy <- -shannon_entropy
  
  return(shannon_entropy)
}

#' Generate Heatmap for Proximity Data
#'
#' This function generates a heatmap to visualize the proximity data of animals over a specified time period.
#' The heatmap represents the amount of time animals spent in close proximity to each other, converted from 
#' seconds to hours. The heatmap is generated using ggplot2, with a custom color gradient.
#'
#' @param count_proximity_list A list of lists, where each sublist contains proximity data (in seconds) between pairs of animals. 
#'   Each entry corresponds to the time (in seconds) that animals spent in close proximity.
#' @param batch A character string representing the batch or experimental group.
#' @param cageChange A numeric or character value indicating a change in the cage setup.
#' @param systemNum A numeric identifier for the experimental system.
#' @param animal_ids A vector of animal IDs corresponding to the animals in the experiment.
#' @param phase A character string indicating the phase of the experiment (e.g., "active" or "inactive").
#' @param phase_number An integer representing the phase number.
#'
#' @return A ggplot2 object representing the heatmap of animal proximity.
#'
#' @details 
#' The proximity data is first converted from seconds to hours. Then, the list of lists is flattened into a matrix, 
#' with animal IDs used as both row and column names. The matrix is melted into a long-format data frame suitable 
#' for ggplot2, where `Var1` and `Var2` represent the animal IDs, and the value represents the number of hours 
#' spent in proximity.
#'
#' The heatmap is created with a color gradient ranging from light blue to dark blue, with specified breaks and 
#' labels for clarity. The title of the plot includes batch, cage change, system number, phase, and phase number.
#'
#' @examples
#' # Example usage
#' count_proximity_list <- list(list(1000, 2000, 3000, 4000), list(1500, 2500, 3500, 4500))
#' animal_ids <- c("A1", "A2", "A3", "A4")
#' ggp <- generateHeatMapProximity(count_proximity_list, "Batch1", "Change1", 1, animal_ids, "active", 1)
#' print(ggp)
#' 
#' @export
generateHeatMapProximity <- function(count_proximity_list, batch, cageChange, systemNum, animal_ids, phase, phase_number) {
  
  # Convert proximity data from seconds to hours
  count_proximity_list_hours <- lapply(count_proximity_list, function(x) ifelse(x != 0, x / 3600, x))
  
  # Convert the list of lists into a matrix
  matrix_data <- do.call(rbind, count_proximity_list_hours)
  
  # Set animal IDs as row and column names
  dimnames(matrix_data) <- list(animal_ids, animal_ids)
  
  # Melt the data: Convert the matrix into a long-format data frame for ggplot
  data_melt <- melt(matrix_data, as.is = TRUE, value.name = "hours")
  
  # Create the heatmap plot using ggplot2
  ggp <- ggplot(data_melt, aes(Var1, Var2)) +
    geom_tile(aes(fill = hours)) +
    scale_fill_gradientn(colors = c("lightblue", "blue", "darkblue"), 
                         limits = c(0, 15),
                         breaks = c(2, 4, 6, 8, 10, 12, 14, ifelse(max(data_melt$hours) < 15, 15, warning("Higher scale required in heatmap for ", systemNum))),
                         labels = c("2", "4", "6", "8", "10", "12", "14", "15")) +
    labs(title = paste(batch, cageChange, systemNum, ": close contact, Phase: ", phase, phase_number), 
         x = "ID", y = "ID")
  
  return(ggp)
}

#' Generate Heatmap for Animal Positions
#'
#' This function generates a heatmap to visualize the positions where animals spent time during the experiment.
#' The heatmap shows the amount of time animals spent at each position, with the time represented in hours.
#' The heatmap is created using ggplot2, with a custom color gradient.
#'
#' @param count_position_list A list of lists containing position IDs and corresponding time spent at those positions 
#'   (in seconds). Each list entry corresponds to one animal's data.
#' @param batch A character string representing the batch or experimental group.
#' @param cageChange A numeric or character value indicating a change in the cage setup.
#' @param systemNum A numeric identifier for the experimental system.
#' @param phase A character string indicating the phase of the experiment (e.g., "active" or "inactive").
#' @param phase_number An integer representing the phase number.
#'
#' @return A ggplot2 object representing the heatmap of animal positions.
#'
#' @details 
#' The list of position data is first converted into a data frame. Position IDs are then mapped to coordinates using a 
#' predefined tibble. The time spent at each position is converted from seconds to hours. The data is then merged 
#' with position coordinates and used to create the heatmap. The x and y coordinates represent the positions on a 
#' grid, and the color of each tile corresponds to the amount of time spent at each position.
#'
#' The heatmap is generated with a custom color gradient from yellow to dark red, with specified breaks and labels.
#' The title of the plot includes batch, cage change, system number, phase, and phase number.
#'
#' @examples
#' # Example usage
#' count_position_list <- list(list(1, 100), list(2, 200), list(3, 300))
#' ggp <- generateHeatMapPositions(count_position_list, "Batch1", "Change1", 1, "active", 1)
#' print(ggp)
#' 
#' @export
generateHeatMapPositions <- function(count_position_list, batch, cageChange, systemNum, phase, phase_number) {
  
  # Convert list of position data into a data frame
  df_positions <- as.data.frame(do.call(rbind, count_position_list))
  
  # Create a tibble with position IDs and corresponding x and y coordinates
  Positions_tibble <- tibble(PositionID = c(1:8), 
                             xPos = c(0, 100, 200, 300, 0, 100, 200, 300), 
                             yPos = c(0, 0, 0, 0, 116, 116, 116, 116))
  
  # Merge position data with coordinates
  merged_df <- merge(df_positions, Positions_tibble, by.x = "V1", by.y = "PositionID", all = TRUE)
  
  # Convert time from seconds to hours
  hour_df <- merged_df %>%
    mutate(V2 = ifelse(V2 != 0, V2 / 3600, V2))
  
  # Rename columns for clarity
  hour_df <- hour_df %>%
    rename(ID = V1) %>%
    rename(hours = V2)
  
  # Create the heatmap using ggplot2
  heatmap <- ggplot(hour_df, aes(x = xPos, y = yPos, fill = hours)) +
    geom_tile() +                                                        # Create tiles based on hours
    scale_x_continuous(breaks = c(0, 100, 200, 300), labels = c("0", "100", "200", "300")) +
    scale_y_continuous(breaks = c(0, 116), labels = c("0", "116")) +
    scale_fill_gradientn(colors = c("yellow", "orange", "red", "darkred", "#290000"), # Custom color palette
                         limits = c(0, 12), 
                         breaks = c(0, 2, 4, 6, 8, 10, ifelse(max(hour_df$hours) < 12, 12, warning("Higher scale required in heatmap for ", systemNum))),
                         labels = c("0", "2", "4", "6", "8", "10", "12")) + 
    labs(title = paste(batch, cageChange, systemNum, ": used positions, Phase: ", phase, phase_number),
         x = "x-axis", 
         y = "y-axis")   # Add title and axis labels
  
  return(heatmap)  # Return the generated ggplot object
}

#' Generate Line Plot for Animal Data Across Phases
#'
#' This function generates a line plot to visualize the total close contact time of animals across different experimental phases.
#' The time spent in close contact is plotted on the y-axis, and the experimental phases are shown on the x-axis. Each animal's
#' data is represented by a separate line.
#'
#' @param data A data frame containing the experimental data with phase and animal-specific time data.
#' @param batch A character string representing the batch or experimental group.
#' @param cageChange A numeric or character value indicating a change in the cage setup.
#' @param animal_ids A vector containing animal IDs for the subjects involved in the experiment.
#' @param system A numeric identifier for the experimental system.
#'
#' @return A ggplot2 object representing the line plot of animal contact time across experimental phases.
#'
#' @details 
#' The function selects the relevant columns from the input data, melts it into a long format, and converts the time data from
#' seconds to hours. The line plot is generated with different lines for each animal, showing the total close contact time
#' across different phases. The color of each line represents a different animal.
#'
#' The x-axis represents the experimental phases, and the y-axis shows the total time spent in close contact, converted to hours.
#' The legend is labeled with the batch, cage change, and system number.
#'
#' @examples
#' # Example usage
#' data <- data.frame(Phase = c("I1", "A1", "I2"), Animal1 = c(3600, 7200, 5400), Animal2 = c(4500, 6000, 3000))
#' plot <- generateGraph(data, "Batch1", "Change1", c("Animal1", "Animal2"), 1)
#' print(plot)
#' 
#' @export
generateGraph <- function(data, batch, cageChange, animal_ids, system) {
  
  # Filter data to select only the relevant columns (Phase and animal IDs)
  subset_data <- select(data, Phase, animal_ids[1], animal_ids[2], animal_ids[3], animal_ids[4])

  # Melt the data into a long format
  long_data <- melt(subset_data, id = 'Phase')
  
  # Rename columns for clarity
  names(long_data) <- c('Phase', 'animal_id', 'time')
  
  # Convert time from seconds to hours
  hour_data <- long_data %>%
    mutate(time = as.integer(time)) %>%
    mutate(time = ifelse(time != 0, time / 3600, time))  # Convert non-zero times to hours
  
  # Create the line plot using ggplot2
  plot <- ggplot(data = hour_data, aes(x = Phase, y = time, color = animal_id)) +
    geom_line(aes(group = animal_id)) +         # Plot lines for each animal
    geom_point() +                             # Add points at each phase
    scale_y_continuous("Total Close Contact in Hours") +   # Y-axis label
    scale_x_discrete(limits = c("I1", "A1", "I2", "A2", "I3", "A3", "I4", "A4", "I5")) +   # Set X-axis limits
    scale_color_discrete(name = paste("Mice from", batch, cageChange, system))   # Legend label
  
  return(plot)  # Return the generated plot
}

#' Create a Plot Comparing Cage Change Data
#'
#' This function generates a plot comparing the total close contact time for different
#' cage change phases for a given animal. The data from four cage change conditions
#' (CC1, CC2, CC3, CC4) are merged and plotted, with time converted from seconds to hours.
#'
#' @param CC1 A tibble containing proximity data for the first cage change condition.
#' @param CC2 A tibble containing proximity data for the second cage change condition.
#' @param CC3 A tibble containing proximity data for the third cage change condition.
#' @param CC4 A tibble containing proximity data for the fourth cage change condition.
#' @param current_animal_id A character string representing the ID of the animal being analyzed.
#' @param batch A character string representing the batch information for the experiment.
#' @param stress_condition A character string representing the stress condition for the experiment.
#'
#' @return A `ggplot` object showing the comparison of total close contact time across cage change conditions.
#' @export
#'
#' @examples
#' plot <- create_joined_table(CC1, CC2, CC3, CC4, "Mouse1", "Batch1", "Control")
#' print(plot)
create_joined_table <- function(CC1, CC2, CC3, CC4, current_animal_id, batch, stress_condition) {
  
  # Filter each tibble to the current mouse ID and rename columns for clarity
  filtered_CC1 <- CC1 %>%
    select(Phase, current_animal_id) %>%
    rename(CC1 = current_animal_id)
  
  filtered_CC2 <- CC2 %>%
    select(Phase, current_animal_id) %>%
    rename(CC2 = current_animal_id)
  
  filtered_CC3 <- CC3 %>%
    select(Phase, current_animal_id) %>%
    rename(CC3 = current_animal_id)
  
  filtered_CC4 <- CC4 %>%
    select(Phase, current_animal_id) %>%
    rename(CC4 = current_animal_id)
  
  # Join the filtered tibbles (one for each cage change) into a single table by 'Phase'
  proximity_table_join <- Reduce(function(...) { merge(..., by = "Phase", all = TRUE) },
                                 list(filtered_CC1, filtered_CC2, filtered_CC3, filtered_CC4))
  
  # Sort the 'Phase' column numerically (ignoring any non-numeric characters)
  proximity_table_join <- proximity_table_join %>%
    arrange(as.integer(sub("[^0-9]", "", Phase)))
  
  # Melt the table to long format for easier plotting
  long_data <- melt(proximity_table_join, id = 'Phase')
  
  # Rename columns for easier understanding
  names(long_data) <- c('Phase', 'cageChange', 'time')
  
  # Convert time from seconds to hours
  hour_data <- long_data %>%
    mutate(time = as.integer(time)) %>%
    mutate(time = ifelse(time != 0, time / 3600, time))  # Avoid division by zero
  
  # Create the plot with ggplot2
  plot <- ggplot(data = hour_data, aes(x = Phase, y = time, color = cageChange)) +
    geom_line(aes(group = cageChange), na.rm = TRUE) +  # Lines connecting points for each cageChange
    geom_point(na.rm = TRUE) +  # Points for each data entry
    scale_color_manual(values = c("CC1" = "aquamarine2", "CC2" = "deepskyblue2", "CC3" = "deepskyblue4", "CC4" = "darkblue"),
                       name = paste("Animal ID:", current_animal_id, batch, stress_condition)) +
    scale_y_continuous("Total Close Contact in Hours") +  # Y-axis label
    scale_x_discrete(limits = c("I1", "A1", "I2", "A2", "I3", "A3", "I4", "A4", "I5")) +  # X-axis limits for phases
    facet_grid(~cageChange)  # Facet by cageChange (one plot per condition)
  
  return(plot)
}

#' Perform Normality Test and Statistical Test for Each Variable and Phase
#'
#' This function filters the data based on the specified phase and sex, checks if the specified variable is numeric, 
#' performs a normality test (Shapiro-Wilk) for each group, and then performs either a Wilcoxon rank-sum test (for two groups)
#' or an ANOVA/Kruskal-Wallis test (for more than two groups). The appropriate post-hoc test is performed as well.
#'
#' @param data A data frame containing the data to be analyzed.
#' @param value A string representing the value or label for the test.
#' @param variableName A character string specifying the name of the dependent variable in the data.
#' @param phase A character string specifying the phase for filtering the data.
#' @param sex A character string specifying the sex for filtering the data.
#'
#' @return A list containing:
#'   \item{testResults}{A list with the test details including p-values, significance levels, and normality results.}
#'   \item{plot}{A plot visualizing the results of the test.}
#'   \item{posthocResults}{A data frame with the post-hoc test results (if applicable).}
testAndPlotVariable <- function(data, value, variableName, phase, sex) {

  filteredData <- data %>%
    filter(if('Phase' %in% colnames(data))Phase == phase else TRUE) %>%   
    filter(Sex == sex)
  
  filteredData$Group <- as.factor(filteredData$Group)
  
  uniqueGroups <- unique(filteredData$Group)
  numGroups <- length(uniqueGroups)
  
  # Check if the variable is numeric (if not, return Null)
  if (is.numeric(filteredData[[variableName]])) {
    if (numGroups == 2) {
      
      group1 <- uniqueGroups[1]
      group2 <- uniqueGroups[2]
      
      dataGroup1 <- filteredData[[variableName]][filteredData$Group == group1]
      dataGroup2 <- filteredData[[variableName]][filteredData$Group == group2]
      
      # Check if there are at least 3 non-NA values in each group
      if (sum(!is.na(dataGroup1)) >= 3 && sum(!is.na(dataGroup2)) >= 3) {
        wilcoxRes <- performWilcoxonTest(dataGroup1, dataGroup2)
        
        # If the Wilcoxon test was successful, return the results and plot
        if (!is.null(wilcoxRes)) {
          testResults <- list(
            Testvalue = value,
            Variable = variableName,
            Phase = phase,
            Sex = sex,
            Test = "Wilcoxon rank-sum test",
            CON_Normality = NA,
            RES_Normality = NA,
            SUS_Normality = NA,
            P_Value = wilcoxRes$p.value,
            Significance_Level = sprintf("%.3f", wilcoxRes$p.value)
          )
          
          # Generate plot for Wilcoxon test
          p <- generatePlot(filteredData, value, variableName, phase, sex)
          
          return(list(testResults = testResults, plot = p, posthocResults = NULL))
        }
      }
    } else {  # numGroups /= 2
      ## NORMALIZATION
      # normalize group CON and RES with Shapiro-Wilk Normality Test
      conNorm <- shapiro.test(filteredData[[variableName]][filteredData$Group == "con"])
      resNorm <- shapiro.test(filteredData[[variableName]][filteredData$Group == "res"])
      
      # normalize group SUS, if it exists, else declare it to 1
      susGroupExists <- any(filteredData$Group == "SUS")
      if (susGroupExists) {
        susNorm <- shapiro.test(filteredData[[variableName]][filteredData$Group == "sus"])
      } else {
        susNorm <- list(p.value = 1)
      }
      #if one of the p values of the normalizations is empty, return NULL
      if (is.na(conNorm$p.value) || is.na(resNorm$p.value) || is.na(susNorm$p.value)) {
        return(NULL)
      }
      # register Normalization test results in testResults for return
      testResults <- list(
        Testvalue = value,
        Variable = variableName,
        Phase = phase,
        Sex = sex,
        CON_Normality = conNorm$p.value,
        RES_Normality = resNorm$p.value,
        SUS_Normality = susNorm$p.value
      )
      
      # initialize empty dataframes for ANOVA
      testResultsDf <- data.frame()
      posthocResultsDf <- data.frame()
      
      print("test6")
      # ANOVA or Kruskal-Wallis test
      if (numGroups > 2) {
        
        if (conNorm$p.value >= 0.05 && resNorm$p.value >= 0.05 && susNorm$p.value >= 0.05) {
          
          anovaTest <- aov(as.formula(paste(variableName, "~ Group")), data = filteredData)
          testResults$Test <- "ANOVA"
          testResults$P_Value <- summary(anovaTest)[[1]][["Pr(>F)"]][1]
          testResults$Significance_Level <- sprintf("%.3f", testResults$P_Value)
          #posthoc for ANOVA
          posthocResultsDf <- performPosthocAnova(filteredData, variableName)
        } else {
          kruskalTest <- kruskal.test(as.formula(paste(variableName, "~ Group")), data = filteredData)
          testResults$Test <- "Kruskal-Wallis"
          testResults$P_Value <- kruskalTest$p.value
          testResults$Significance_Level <- sprintf("%.3f", p.adjust(testResults$P_Value, method = "BH"))
          
          posthocResultsDf <- performPosthocKruskal(filteredData, variableName)
        }
        
        print("test7")
        if (!is.null(posthocResultsDf) && ncol(posthocResultsDf) > 0) {
          if (identical(names(testResultsDf), names(posthocResultsDf))) {
            testResultsDf <- bind_rows(testResultsDf, posthocResultsDf)
          } else {
            testResultsDf <- plyr::rbind.fill(testResultsDf, posthocResultsDf)
          }
        }
      }
      
      p <- generatePlot(filteredData, value, variableName, phase, sex)
      
      return(list(testResults = testResults, plot = p, posthocResults = testResultsDf))
    }
  } else {
    cat("Variable", variableName, "is not numeric and will be skipped.\n")
    return(NULL)
  }
}

# Function to generate plots for each variable and phase
generatePlot <- function(preprocessed_data, value, variableName, phase, sex) {
  filteredData <- preprocessed_data #%>%
  #  filter(if (include_phase) Phase == phase else TRUE) %>%  # Include/exclude "Phase" based on the variable
  #  filter(if (include_sex) Sex == sex else TRUE)  # Include/exclude "Sex" based on the variable
  
  p <- ggplot(filteredData, aes(Group, .data[[variableName]], color = Group)) +
    # Customize plot aesthetics and labels
    scale_x_discrete(name = NULL, expand = c(0.3, 0.1)) + 
    scale_y_continuous("avg rank", expand = c(0.1, 0.1)) +
    geom_jitter(aes(fill = Group), size = 4, alpha = 0.7, width = 0.2, shape = 16, na.rm = TRUE) +
    stat_summary(
      fun.min = function(z) {quantile(z, 0.25)},
      fun.max = function(z) {quantile(z, 0.75)},
      fun = median,
      color = "black",
      size = 0.8,
      shape = 16,
      na.rm = TRUE
    ) +
    labs(title = bquote(~bold(.(value))),
         subtitle = paste("(", phase, " , ", sex, ")", sep = ""),
         caption = "",
         x = NULL,
         y = "z score [a.u.]") +
    scale_color_manual(name = NULL, values = c("#1e3791", "#76A2E8", "#F79719")) +
    scale_fill_manual(name = NULL, values = c("#1e3791", "#76A2E8", "#F79719")) +
    cowplot::theme_minimal_hgrid(12, rel_small = 1) +
    theme(plot.title = element_text(hjust = 0.5, face = "plain"),
          plot.subtitle = element_text(hjust = 0.5, size = 10, face = "plain"),
          legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_text(),
          axis.ticks.x = element_blank())
  
  return(p)
}

#' Perform Post-Hoc Pairwise Tests for ANOVA
#'
#' This function performs ANOVA followed by pairwise post-hoc tests for the specified variable.
#' It adjusts p-values using the Bonferroni method and removes duplicate comparisons.
#'
#' @param data A data frame containing the data to be analyzed.
#' @param variableName A character string specifying the name of the dependent variable in the data.
#'
#' @return A data frame with pairwise comparison results including the adjusted p-values.
#' @export
#'
#' @examples
#' posthoc_results <- performPosthocAnova(my_data, "score")
performPosthocAnova <- function(data, variableName) {
  
  anovaTest <- aov(as.formula(paste(variableName, "~ Group")), data = data)
  
  pairwiseResults <- rstatix::pairwise_t_test(as_data_frame(data), formula = as.formula(paste(variableName, "~ Group")),
                                     p.adjust.method = "bonferroni")
  
  # Remove duplicate comparisons
  pairwiseResults <- pairwiseResults[!duplicated(pairwiseResults[, c("group1", "group2")]), ]
  pairwiseResults$GroupComparison <- paste(pairwiseResults$group1, "vs.", pairwiseResults$group2)
  
  return(pairwiseResults)
}

#' Perform Post-Hoc Pairwise Tests for Kruskal-Wallis Test
#'
#' This function performs Kruskal-Wallis test followed by pairwise post-hoc tests using Dunn's test.
#' P-values are adjusted using the Holm method.
#'
#' @param data A data frame containing the data to be analyzed.
#' @param variableName A character string specifying the name of the dependent variable in the data.
#'
#' @return A data frame with pairwise comparison results including the adjusted p-values.
#' @export
#'
#' @examples
#' posthoc_results <- performPosthocKruskal(my_data, "score")
performPosthocKruskal <- function(data, variableName) {
  kruskalTest <- kruskal.test(as.formula(paste(variableName, "~ Group")), data = data)
  pairwiseResults <- dunn_test(data, formula = as.formula(paste(variableName, "~ Group")),
                               p.adjust.method = "holm")
  pairwiseResults$GroupComparison <- paste(pairwiseResults$group1, "vs.", pairwiseResults$group2)
  
  return(pairwiseResults)
}

#' Perform Wilcoxon Rank-Sum Test
#'
#' This function performs a Wilcoxon rank-sum test (Mann-Whitney U test) to compare two groups.
#' It checks if each group contains at least 3 observations before performing the test.
#'
#' @param dataGroup1 A numeric vector containing the data for the first group.
#' @param dataGroup2 A numeric vector containing the data for the second group.
#'
#' @return A Wilcoxon test result object or NULL if there are fewer than 3 observations in either group.
#' @export
#'
#' @examples
#' wilcoxon_result <- performWilcoxonTest(group1_data, group2_data)
performWilcoxonTest <- function(dataGroup1, dataGroup2) {
  if (length(dataGroup1) >= 3 && length(dataGroup2) >= 3) {
    return(wilcox.test(dataGroup1, dataGroup2))
  } else {
    return(NULL)
  }
}
