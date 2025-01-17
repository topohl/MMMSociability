
## 11/2023
## Anja Magister
## ANALYSIS OF ANIMAL POSITIONS - FUNCTIONS ##
##

##############################################################################################################
# input:  x_Pos - x coordinate of cage
#         y_Pos - y coordinate of cage
#         lookup_tibble - table of defined position IDs for every possible coordinate couple
#
# output: Position Id
# effect: finds corresponding PositionID from two coordinates(x_Pos, y_pos) in lookup_tibble and returns ID 

find_id <- function(x_Pos, y_Pos, lookup_tibble){
  
  #give non-standard positions a standard position
  # y value
  if(y_Pos < 116) {y_Pos <- 0}
  else if(y_Pos >= 116) {y_Pos <- 116}
  
  # x value
  if(x_Pos<100) {x_Pos <- 0}
  else if(x_Pos<200) {x_Pos <- 100}
  else if(x_Pos<300) {x_Pos <- 200}
  else if(x_Pos>=300) {x_Pos <- 300}

  # search position in lookup table
  result <- lookup_tibble %>%
    filter(xPos == x_Pos, yPos == y_Pos) %>%
    select(PositionID)
  
  # return new ID of specific position
  if (nrow(result) > 0) {
    return(result$PositionID)
  } else {
    cat("xpos: ", x_Pos, " , ypos: ", y_Pos, "\n")
    print("NA!")
    return(NA) # Return NA if no match found
  }
}

################################################################################################################
# input:  processed_data - the tibble with the whole data that we want to change 
#         
#        
#
# output: processed_data
# effect: - adds extra rows of RFID information to the time when a Phase changes(6h30 a.m. &6h30 p.m.)
#         - we need a row for every mouse in every system bc the data will later filtered by the phases 
#         and we need a start value of every mouse for every phase
#         - the new rows will be added to the processed_data tibble

addPhaseMarkerRows <- function(processed_data){
  
  #save dates of experiment(the date of the days)
  onlyDates <- as.Date(processed_data$DateTime)
  uniqueDates <- unique(onlyDates)
  #print(uniqueDates)
  
  #save existing systems
  uniqueSystems <- unique(processed_data$System)
  #print(uniqueSystems)
  
  for(i in 1:length(uniqueDates)) {
    date <- uniqueDates[i]

    #define marker points depending on actual date
    # "GMT" means "UTC"
    early_marker_end_active <- as.POSIXct(paste(date,"06:29:59"), tz = "GMT")  #end active phase
    early_marker_start_inactive <- as.POSIXct(paste(date,"06:30:00"), tz = "GMT")  #start inactive phase
    early_marker_end_inactive <- as.POSIXct(paste(date,"18:29:59"), tz = "GMT")  #end inactive phase
    late_marker_start_active <- as.POSIXct(paste(date,"18:30:00"), tz = "GMT")   #start active phase

    #we have 2 early markers and 2 late markers
    for(marker in c(early_marker_end_active, early_marker_start_inactive, early_marker_end_inactive, late_marker_start_active)){
     
      # filter rows which are earlier than the marker
      filtered_time <- processed_data %>%
        as_tibble() %>%
        filter(DateTime < marker)
      
      #we need a marker for every existing system and for every mouse in every system
      for(system in uniqueSystems){

        # filter rows from specific system 
        filtered_system <- filtered_time %>%
          as_tibble() %>%
          filter(System == system)

        #check which mice exists in this system
        mice_names <- unique(filtered_system$AnimalID)

        for(mouse in mice_names){
          # filter specific mouse row from system with nearest DateTime to 18:30:00 or 06:30:00 but still earlier
          filtered_row <- filtered_system %>%
            as_tibble() %>%
            filter(AnimalID == mouse) %>%
            slice_max(order_by = DateTime)

          #if time on this day is recorded
          if(length(filtered_row) != 0){
            #use copied row as new one and change the date to marker date
            filtered_row <- filtered_row %>%
              mutate(DateTime = as.POSIXct(marker, tz = "UTC"))
            
            #add new row to processed_data
            processed_data <- add_row(processed_data,filtered_row)
          }
        }
      }
    }
  }
  
  #sort tibble again by dateTime to bring new entries to correct position
  processed_data <- processed_data %>%
    as_tibble() %>%
    arrange(., DateTime)

  #remove last 20 rows because they are 20 unnecessary marker rows of the last day and the late marker
  processed_data <- processed_data %>% filter(row_number() <= n()-20)

  return(processed_data)
}

################################################################################################################
# input:  processed_data - the tibble with the whole data that we want to change 
#         
#        
#
# output: processed_data 
# effect: - adds extra rows of RFID information to the time when a day changes(midnight)
addDayMarkerRows <- function(processed_data) {
  #save dates of experiment
  onlyDates <- as.Date(processed_data$DateTime)
  uniqueDates <- unique(onlyDates)
  
  #save existing systems
  uniqueSystems <- unique(processed_data$System)
  
  for(i in 1:length(uniqueDates)) {
    date <- uniqueDates[i]
    
    #define marker point depending on actual date
    # "GMT" means "UTC"
    marker <- as.POSIXct(paste(date, "00:00:00"), tz = "GMT")  #start of new day
    
    # filter rows which are earlier than the marker
    filtered_time <- processed_data %>%
      as_tibble() %>%
      filter(DateTime < marker)
    
    for(system in uniqueSystems){
      
      # filter rows from specific system 
      filtered_system <- filtered_time %>%
        as_tibble() %>%
        filter(System == system)

      #check which mice exists in this system
      mice_names <- unique(filtered_system$AnimalID)
      
      for(mouse in mice_names){
        # filter specific mouse row from system with nearest DateTime to 18:30:00 or 06:30:00 but still earlier
        filtered_row <- filtered_system %>%
          as_tibble() %>%
          filter(AnimalID == mouse) %>%
          slice_max(order_by = DateTime)

        #if time on this day is recorded
        if(length(filtered_row) != 0){
          #use copied row as new one and change the date to marker date
          filtered_row <- filtered_row %>%
            mutate(DateTime = as.POSIXct(marker, tz= "UTC"))
          
          #add new row to processed_data
          processed_data <- add_row(processed_data, filtered_row)
        }
      }
    }
  }
  #sort tibble again by dateTime to bring new entries to correct position
  processed_data <- processed_data %>%
    as_tibble() %>%
    arrange(., DateTime)
  
  return(processed_data)
}

#' Count Consecutive Phases
#' 
#' This function counts the consecutive phases of activity and inactivity in the data.
#' It adds two new columns to the input tibble: `ConsecActive` and `ConsecInactive`,
#' which represent the count of consecutive active and inactive phases, respectively.
#' 
#' @param data A tibble with a column named `Phase` containing the phases of activity ("Active" or "Inactive").
#' 
#' @return A tibble with the original data and two additional columns: `ConsecActive` and `ConsecInactive`.
#' @export 
#' 
#' @examples
#' data <- tibble::tibble(Phase = c("Active", "Active", "Inactive", "Active", "Inactive", "Inactive"))
#' count_consecPhases(data)
#' return (data)

consecPhases <- function(data){
  
  # initialize counters
  active_phases <- 0
  inactive_phases <- 0
  
  for(i in 1:nrow(data)){
    #select row from data
    current_row <- data[i,]

    #first row
    if(i == 1){
      
      if(current_row$Phase == "Active"){
        current_row$ConsecActive <- 1
        current_row$ConsecInactive <- 0
        active_phases <- active_phases + 1
      } else {    #current_row$Phase == "Inactive"
        current_row$ConsecActive <- 0
        current_row$ConsecInactive <- 1
        inactive_phases <- inactive_phases + 1
      }
    } else {
      
      #previous row from data
      previous_row <- data[(i - 1),]
      
      if(current_row$Phase != previous_row$Phase){
        if(current_row$Phase == "Active") {
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
  return (data)
}

##############################################################################################################
# input:  system_animal_ids -  the ids/names of the mice in the current system
#         data -  the tibble with the data of the current system
#         animal_list - the list that contains the information of the mouse ID and the current time and position
#         
#
# output: animal_list - updated
# effect: - the mice-list will be entered here for the first time, information gained from the data
#         - takes the first entry of every mouse in data and copies the time and position that we know the start position of each mouse

# find the FIRST TIME where mouse is tracked in the cage
# aka first value of mouse in processed_data_final
find_first_pos_and_time <- function(system_animal_ids, data, animal_list){
  
  # create empty vector with 4 variables
  times_vec <- rep(NA, times = 4)
  
  print(system_animal_ids)
  for (i in 1:length(system_animal_ids)){ #i=1-4
    
    #define current animal_id, first_position and first_time
    animal_id <- system_animal_ids[[i]]
    if(is.na(animal_id)){
      animal_id <- paste0("lost_", i)
      first_time <- 0        #time value should be adapted at the end of this function
      first_position <- (-1)  #add non existing position for non existing mouse
    } else {
      #search first entry in whole data
      first_entry <- data %>%
        filter(AnimalID == animal_id) %>%
        slice(1) #first row

      first_time <- first_entry$DateTime
      first_position <- first_entry$PositionID
    }
    
    #write name, position and time into animal_list
    animal_list[i][[1]] <- animal_id
    #print(animal_id)
    
    animal_list[[i]][[3]] <- first_position
    #print(first_position)
    
    animal_list[[i]][[2]] <- first_time
    #print(first_time)
    
    #enter time in times vec fot later to compare
    times_vec[i] <- first_time
  }
  
  #check, if all times are similar
  #if not, change them all to the first shared time
  if(length(unique(times_vec)) != 1) {
    latest_time <- max(times_vec)
    #print(latest_time)
    for(i in 1:4) {
      animal_list[[i]][[2]] <- latest_time
    }
  }
  
  #print("animal_list inside the first pos function: \n")
  #print(animal_list)
  return(animal_list)
}

#' Compute Proximity Between Mice
#'
#' This function calculates the social proximity between mice based on their positions. It updates the `count_proximity_list` 
#' to include the duration (in seconds) that each pair of mice remained in the same position.
#'
#' @param old_animal_list A list containing the previous positions of the mice from the last recorded time until the second 
#'   before the current time.
#' @param new_animal_list A list containing the updated positions of the mice for the current second.
#' @param count_proximity_list A list to store the results, tracking the number of seconds each pair of mice remained 
#'   in proximity.
#' @param secTemp An integer indicating the time difference (in seconds) between the `old_animal_list` and the current time.
#'
#' @return The updated `count_proximity_list` with the total seconds of social contact between mice.
#' @details
#' The function performs the following steps:
#' - Compares positions from `old_animal_list` for the time period defined by `secTemp - 1` and updates proximity counts.
#' - Compares positions from `new_animal_list` for the current second and updates proximity counts.
#' - Ensures proximity is only counted if two different individuals are in the same position.
#'
#' @examples
#' # Assuming old_animal_list, new_animal_list, and count_proximity_list are defined:
#' count_proximity_list <- compute_proximity(old_animal_list, new_animal_list, count_proximity_list, secTemp)
compute_proximity <- function(old_animal_list, new_animal_list, count_proximity_list, secTemp) {
  
  # Loop through all pairs of mice (4 assumed in both old and new lists)
  for (i in 1:4) {
    for (j in i:4) {
      # Compare positions in old_animal_list for (secTemp - 1) seconds
      if (old_animal_list[[i]][[3]] == old_animal_list[[j]][[3]]) {
        count_proximity_list[[i]][[j]] <- count_proximity_list[[i]][[j]] + (secTemp - 1)
        if (j != i) {
          count_proximity_list[[j]][[i]] <- count_proximity_list[[j]][[i]] + (secTemp - 1)
        }
      }
      # Compare positions in new_animal_list for the current second
      if (new_animal_list[[i]][[3]] == new_animal_list[[j]][[3]]) {
        count_proximity_list[[i]][[j]] <- count_proximity_list[[i]][[j]] + 1
        if (j != i) {
          count_proximity_list[[j]][[i]] <- count_proximity_list[[j]][[i]] + 1
        }
      }
    }
  }
  
  # Return the updated proximity list
  return(count_proximity_list)
}

#' Update Mouse Position Counts
#'
#' This function updates the count of seconds that mice have been standing on specific positions.
#'
#' @param old_animal_list A tibble containing the old positions of the mice from the last entered time until the second before the current time.
#' @param new_animal_list A tibble containing the new updated positions for the current second.
#' @param count_position_list A list to store the results, containing the amount of seconds a mouse (or multiple mice) were standing on a position.
#' @param secTemp The time difference between the new and old mice lists.
#'
#' @return A list containing the updated count of seconds that mice have been standing on specific positions.
#' @export
#'
#' @examples
#' old_animal_list <- tibble::tibble(id = 1:3, position = c("A", "B", "C"))
#' new_animal_list <- tibble::tibble(id = 1:3, position = c("A", "C", "B"))
#' count_position_list <- list(A = 5, B = 3, C = 2)
#' secTemp <- 1
#' update_position_counts(old_animal_list, new_animal_list, count_position_list, secTemp)
compute_position <- function(old_animal_list, new_animal_list, count_position_list, secTemp) {
  
  #vectors to check which positions were used in this time period
  old_positions <- c()
  new_positions <- c()
  
  #for every mouse in the mice list
  for (i in 1:4) {
    old_pos <- as.numeric(old_animal_list[[i]][[3]])
    new_pos <- as.numeric(new_animal_list[[i]][[3]])

    # Count used position only once per second, even if multiple mice are on it
    # Check if position is already added to the list from another mouse
    if(!old_pos %in% old_positions & old_pos != (-1)) {

      #add new seconds to count_position_list
      count_position_list[[old_pos]][[2]] <- count_position_list[[old_pos]][[2]] + (secTemp - 1)

      #add old position into vector of used positions
      old_positions <- append(old_positions, old_pos)
    }
    if(!new_pos %in% new_positions& new_pos != (-1)) {
      #print("add new pos")
      
      #add new seconds to count_position_list
      count_position_list[[new_pos]][[2]] <- count_position_list[[new_pos]][[2]] + 1
      
      #add new position into vector of used positions
      new_positions <- append(new_positions, new_pos)
    }
  }
  
  # return updated list of mice that are close to each other
  return(count_position_list)
}

#' Compute Total Social Proximity Time for Each Mouse
#'
#' This function calculates the total time (in seconds) each mouse spends in close social contact with other mice.
#' The results are stored in the `total_proximity_list`, which is updated with the proximity duration for each mouse.
#'
#' @param old_animal_list A list containing the previous positions of the mice from the last recorded time 
#'   until the second before the current time. Each sublist contains details about a mouse, including its ID and position.
#' @param new_animal_list A list containing the updated positions of the mice for the current second. 
#'   The format is consistent with `old_animal_list`.
#' @param total_proximity_list A list to store the total proximity time for each mouse. Each sublist contains the mouse ID 
#'   and the cumulative proximity duration in seconds.
#' @param secTemp An integer representing the time difference (in seconds) between the previous entry 
#'   and the current entry.
#' @param excludeHomecage Logical. If `TRUE`, positions with ID 5 (homecage) will be excluded from proximity calculations.
#'
#' @return The updated `total_proximity_list`, with the total proximity time in seconds for each mouse.
#' @details
#' The function performs the following operations:
#' - For each mouse, checks if it was in proximity with other mice during the time between the previous entry 
#'   (`old_animal_list`) and the current time.
#' - Updates the proximity count for each mouse if it shares a position with other mice.
#' - Ensures the IDs in `total_proximity_list` match those in `old_animal_list` for consistency.
#'
#' @examples
#' # Assuming old_animal_list, new_animal_list, and total_proximity_list are defined:
#' total_proximity_list <- compute_total_proximity(old_animal_list, new_animal_list, total_proximity_list, secTemp, TRUE)
compute_total_proximity <- function(old_animal_list, new_animal_list, total_proximity_list, secTemp) {

  # Iterate through each mouse in the list
  for (mouse in 1:4) {
    
    # Create a list of other mice indices, excluding the current mouse
    other_mice <- setdiff(1:4, mouse)
    
    # Retrieve the current mouse's ID from old_animal_list
    mouse_id <- old_animal_list[[mouse]][[1]]
    
    # Conditional for excluding homecage (PositionID 5)
    exclude_condition_old <- if (excludeHomecage) old_animal_list[[mouse]][[3]] != 5 else TRUE
    exclude_condition_new <- if (excludeHomecage) new_animal_list[[mouse]][[3]] != 5 else TRUE

    # Check proximity in old positions (for secTemp - 1 seconds)
    if (exclude_condition_old && 
        (old_animal_list[[mouse]][[3]] %in% sapply(other_mice, function(idx) old_animal_list[[idx]][[3]]))) {
      
      # Add proximity time to total_proximity_list for the current mouse
      if (total_proximity_list[[mouse]][[1]] == mouse_id) { # Verify matching IDs
        total_proximity_list[[mouse]][[2]] <- as.numeric(total_proximity_list[[mouse]][[2]]) + (secTemp - 1)
      } else {
        stop("Mouse IDs do not match in compute_total_proximity")
      }
    }
    
    # Check proximity in new positions (for the current second)
    if (exclude_condition_new && 
        (new_animal_list[[mouse]][[3]] %in% sapply(other_mice, function(idx) new_animal_list[[idx]][[3]]))) {
      
      # Add one second of proximity time to total_proximity_list for the current mouse
      if (total_proximity_list[[mouse]][[1]] == mouse_id) { # Verify matching IDs
        total_proximity_list[[mouse]][[2]] <- as.numeric(total_proximity_list[[mouse]][[2]]) + 1
      } else {
        stop("Mouse IDs do not match in compute_total_proximity")
      }
    }
  }
  
  # Return the updated total proximity list
  return(total_proximity_list)
}


##############################################################################################################
# input: old_animal_list,new_animal_list,count_movement_list, secTemp
#         
#        
#
# output: count_movement_list
# effect: 

# 
compute_movement <- function(old_animal_list,new_animal_list,count_movement_list, secTemp) {
  #variable check if anyone moved at all
  moves <- 0
  
  #for every mouse in the mice list
  for (i in 1:4) {
    old_pos <- as.numeric(old_animal_list[[i]][[3]])
    #print("old_pos")
    #print(old_pos)
    new_pos <- as.numeric(new_animal_list[[i]][[3]])
    #print("new_pos")
    #print(new_pos)
    
    #check if mouse changed place in this second, how long is not important
    #also checking for lost_mouse with position -1(not valuable)
    if(old_pos!=(-1) & new_pos != (-1) & old_pos!=new_pos) {
      #enter one movement to current mouse
      count_movement_list[[i]][[2]] <- as.numeric(count_movement_list[[i]][[2]]) + 1
      #enter one movement to whole system
      count_movement_list[[5]][[2]] <- as.numeric(count_movement_list[[5]][[2]]) + 1
      
      #add move to moves
      moves <- moves + 1
    }
    
    
  }
  
  if(moves == 0) {
    message("compute_movement: No movement detected in this interval.")
  }
  
  # return updated list of mice that are close to each other
  return(count_movement_list)
}
##############################################################################################################
# input:  old_time
#         
#        
#
# output: new_time
# effect: function to do a shift in time(one second forward)

# 
sec_shift <- function(old_time) {
  #put one second on top of old_time
  new_time <- old_time %>%
    as.numeric() %>%
    + 1 %>%
    as.character()
  #cat("old time: ", old_time, "\n")
  #cat("new time: ", new_time, "\n")
  return(new_time)
}

##############################################################################################################
# input:  system_animal_ids -  the ids/names of the mice in the current system
#         animal_list - the list containing the information of the mouse ID and the current time and position
#         data -  the tibble with the data of the current system
#         time - the time that we are looking at, can be written in several lines bc we have several mice
#         line - a counter for the actual line in our data that we are looking at
#
# output: 
# effect: - update animal_list(if its possible) and return it
#         - takes next entry in data and updates the new information in the mice list
#         - similarity to find_first_pos_and_time
update_animal_list <- function(system_animal_ids, animal_list, data, time, line) {
  
  #next_second <- sec_shift(time)
  new_time <- as.numeric(data[line,"DateTime"])

  # write sec difference between new and old time into secTemp
  animal_list[["data_temp"]][["secTemp"]] <- new_time-as.numeric(time)
  
  # write new time into every mouse information
  for(i in 1:4) {animal_list[[i]][[2]] <- new_time}
 
  # while line(and especially the next lines) is still same time
  while(as.numeric(data[line,"DateTime"]) == new_time) {
    # write new position into special mouse
    for(i in 1:4) {
      if(animal_list[[i]][[1]] == as.character(data[line,"AnimalID"])) {animal_list[[i]][[3]] <- as.numeric(data[line,"PositionID"])}
    }
    #if line is not the last line, check next line, else break while loop
    if(line == nrow(data)) {
      line <- line + 1
      break
    }
    #continue with the while condition
    line <- line + 1
  }
  
  # write new line into animal_list
  animal_list[["data_temp"]][["lineTemp"]] <- line
  
  return(animal_list)
}


########################################################################################
## functions for comparing-file

# input:  
#         
#        
#
# output: 
# effect:


compute_rank <- function(rank_tibble, vector, system, cageChange, batch, hours_column_name, rank_column_name) {
  
  #sort the total hour system vec by size, largest value gets rank 1
  vector <- sort(vector)
  #print(vector)
  #enter ranks in rank tibble
  for(i in 1:length(vector)) {
    #enter rank in corresponding line in rank tibble
    rank_tibble <- rank_tibble %>%
      mutate({{rank_column_name}} := ifelse((Batch == batch) & (CageChange == cageChange) & (System == system) & (.data[[hours_column_name]] == vector[i]), i,.data[[rank_column_name]]))
  }
  return(rank_tibble)
}


####################################

# input:  
#         
#        
#
# output: 
# effect:

translate_rank_in_score_vec <- function(rank_vec) {
  #change rank into score
  #rank = score
  # 1   =   4
  # 2   =   3
  # 3   =   2
  # 4   =   1
  
  for (i in 1:length(rank_vec)) {
    #change rank into score 
    score <- ifelse(rank_vec[i] == 1, 4, ifelse(rank_vec[i] == 2, 3, ifelse(rank_vec[i] == 3, 2, ifelse(rank_vec[i] == 4, 1, NA))))
    #overwrite rank with score
    rank_vec[i] <- score
  }
  
  score_vec <- rank_vec
  return(score_vec)
}  


############################################################################################
## FUNCTIONS FROM ANALYZING SHANNON ##
############################################################################################


# input:  
#         
#        
#
# output: 
# effect:
check_cage_prob <- function(old_animal_list,new_animal_list,cage_prob_list,secTemp) {
  #create a vector of 8 ints , all 0
  old_cage_positions <- c(integer(8))
  new_cage_positions <- c(integer(8))
  for(i in 1:4) {
    old_pos <- as.numeric(old_animal_list[[i]][[3]])
    new_pos <- as.numeric(new_animal_list[[i]][[3]])
    
    # add +1 to the position number the mouse sits on (pos nr is the index of the position vector)
    old_cage_positions[old_pos] <- old_cage_positions[old_pos] + 1
    new_cage_positions[new_pos] <- new_cage_positions[new_pos] + 1
    
  }
  #print(old_cage_positions)
  #print(new_cage_positions)
  #divide the whole position vectors throug 4(4 mice)
  old_cage_positions <- old_cage_positions / 4
  new_cage_positions <- new_cage_positions / 4
  
  #enter probability of positions into cage_prob_list:
  for(i in 1:8) {
    #one time for the last second of the old list
    cage_prob_list[[i]][[2]] <- cage_prob_list[[i]][[2]] + 1 #enter the amount of seconds from this position
    cage_prob_list[[i]][[3]] <- cage_prob_list[[i]][[3]] + old_cage_positions[i]
    #and (secTemp-1)-times for the new list 
    for(second in 1:(secTemp - 1)) {
      cage_prob_list[[i]][[2]] <- cage_prob_list[[i]][[2]] + 1 
      cage_prob_list[[i]][[3]] <- cage_prob_list[[i]][[3]] + old_cage_positions[i]  #add the probability (secTemp-1)-times to the prob-number, bc we look at (secTemp-1)different seconds
    }
  }
  return(cage_prob_list)
}


# input:  
#         
#        
#
# output: 
# effect:
check_mice_prob <- function(old_animal_list, new_animal_list, mice_prob_tibble, secTemp) {
  
  for(i in 1:4) {
    
    #if mouse is not tracked(incomplete system)
    if(is.na(animal_ids[i])) {next}

    #define old and new(current) position
    old_pos <- as.numeric(old_animal_list[[i]][[3]])
    new_pos <- as.numeric(new_animal_list[[i]][[3]])
    
    #row in which mice and position fits
    #animal_ids[i], old_pos
    row_old <- which(mice_prob_tibble$AnimalID == animal_ids[i] & mice_prob_tibble$Position == old_pos)
    row_new <- which(mice_prob_tibble$AnimalID == animal_ids[i] & mice_prob_tibble$Position == new_pos)
    
    #enter probability of positions into mice_prob_tibble: 
    mice_prob_tibble[["SumPercentage"]][[row_old]] <- mice_prob_tibble[["SumPercentage"]][[row_old]] + 1
    mice_prob_tibble[["SumPercentage"]][[row_new]] <- mice_prob_tibble[["SumPercentage"]][[row_new]] + (secTemp - 1)
    #enter amount of observed seconds into mice_prob_tibble: 
    mice_prob_tibble <- mice_prob_tibble %>%
      mutate(Seconds=ifelse(AnimalID == animal_ids[i], Seconds + secTemp, Seconds))
  }
  return(mice_prob_tibble)
}

# input:  
#         
#        
#
# output: 
# effect:
calc_shannon_entropy <- function(prob_vec) {
  if(length(prob_vec) != 8) {
    print(prob_vec)
    stop("error in shannon calculation, prob_vec not the rigth size")
  }
  #calculate shannon entropy with the known formula
  shannon_entropy <- numeric()
  for(i in 1:8) {
    if(prob_vec[i] == 0) {
      shannon_entropy <- sum(shannon_entropy, 0) #when no one was in position i, i will be 0 and log2(0) is not defined. thats why I declare this sum as 0
    } else {
      shannon_entropy <- sum(shannon_entropy, prob_vec[i] * log2(prob_vec[i]))}
  }
  shannon_entropy <- -shannon_entropy
  return(shannon_entropy)
}

############################################################################################
## PLOTS ##
######################################################################
# HEATMAPS



# input:  
#         
#        
#
# output: 
# effect:

generateHeatMapproximity <- function(count_proximity_list, batch, cageChange, systemNum, animal_ids, phase, phase_number) {
  # calculate second entrys to hour entrys
  count_proximity_list_hours <- lapply(count_proximity_list, function(x) ifelse(x != 0, x / 3600, x))

  # convert list of lists into a matrix
  matrix_data <- do.call(rbind, count_proximity_list_hours)
  
  # names of the mice Ids
  dimnames(matrix_data) <- list(animal_ids, animal_ids)
  
  # melt the data, means create values combinations out of the matrix
  data_melt <- melt(matrix_data, as.is = TRUE, value.name = "hours")
  
  #create the plot
  ggp <- ggplot(data_melt, aes(Var1, Var2)) +
    geom_tile(aes(fill = hours)) +
    scale_fill_gradientn(colors = c("lightblue","blue","darkblue"), 
                         limits = c(0, 15), 
                         breaks = c(2, 4, 6, 8, 10, 12, 14, ifelse(max(data_melt$hours)<15,15,warning("Higher scale required in heatmap for ", systemNum))),
                         labels = c("2", "4", "6", "8", "10", "12", "14", "15"))+
    labs(title = paste(batch, cageChange, systemNum,":close contact, Phase: ", phase, phase_number), x = "ID", y = "ID")
  
  return(ggp)            
}

# input:  
#         
#        
#
# output: 
# effect:
generateHeatMapPositions <- function(count_position_list, batch, cageChange, systemNum, phase, phase_number) {
  #list into dataframe
  df_positions <- as.data.frame(do.call(rbind, count_position_list))
  #create old posizions tibble(translates position ids into coordinates)
  Positions_tibble <- tibble(PositionID = c(1:8), xPos = c(0, 100, 200, 300, 0, 100, 200, 300), yPos = c(0, 0, 0, 0, 116, 116, 116, 116))
  #merge both dataframes
  merged_df <- merge(df_positions, Positions_tibble, by.x = "V1", by.y = "PositionID", all = TRUE)
  # calculate second entrys to hour entrys
  hour_df <-  merged_df %>%
    mutate(V2 = ifelse(V2 != 0, V2 / 3600, V2))
  
  #rename the columns
  hour_df <- hour_df %>%
    rename(ID = V1) %>%
    rename(hours = V2)

  # create heatmap with ggplot2
  heatmap <- ggplot(hour_df, aes(x = xPos, y = yPos, fill = hours)) +
    geom_tile() +
    scale_x_continuous(breaks = c(0, 100, 200, 300), labels = c("0", "100", "200", "300")) +
    scale_y_continuous(breaks = c(0, 116), labels = c("0", "116")) +
    scale_fill_gradientn(colors = c("yellow","orange", "red", "darkred", "#290000"), # colour palette
                         limits = c(0, 12), 
                         breaks = c(0, 2, 4, 6, 8, 10, ifelse(max(hour_df$hours) < 12, 12, warning("Higher scale required in heatmap for ", systemNum))),
                         labels = c("0", "2", "4", "6", "8", "10", "12")) + 
                         #limits are the borders of the scale,breaks are actual value breaks, labels are names for breakpoints
    labs(title = paste(batch, cageChange, systemNum, ": used positions, Phase: ", phase, phase_number),
         x = "x-axis",
         y = "y-axis")

  return(heatmap)            
}

#graph

#' Generate Graph of Animal Positions Over Phases
#'
#' This function generates a line graph showing the total close contact time 
#' of animals over different phases. The time is converted from seconds to hours.
#'
#' @param data A data frame containing the data to be plotted.
#' @param batch A character string indicating the batch of the experiment.
#' @param cageChange A character string indicating the cage change information.
#' @param animal_ids A character vector of length 4 containing the IDs of the animals.
#' @param system A character string indicating the system used in the experiment.
#'
#' @return A ggplot object representing the graph.
#'
#' @import dplyr
#' @import ggplot2
#' @import reshape2
#'
#' @examples
#' \dontrun{
#' data <- read.csv("path_to_data.csv")
#' batch <- "Batch1"
#' cageChange <- "CageChange1"
#' animal_ids <- c("Animal1", "Animal2", "Animal3", "Animal4")
#' system <- "System1"
#' plot <- generateGraph(data, batch, cageChange, animal_ids, system)
#' print(plot)
#' }

generateGraph <- function(data, batch, cageChange, animal_ids, system) {

  # Filter data to include only the necessary columns
  subset_data <- select(data, Phase, animal_ids[1], animal_ids[2], animal_ids[3], animal_ids[4])

  # Melt the data to long format
  long_data <- melt(subset_data, id = 'Phase')
  
  # Rename columns for clarity
  names(long_data) <- c('Phase', 'mouse', 'time')

  # Convert time from seconds to hours
  hour_data <- long_data %>%
    mutate(time = as.integer(time)) %>%
    mutate(time = ifelse(time != 0, time / 3600, time))

  # Generate the plot
  plot <- ggplot(data = hour_data, aes(x = Phase, y = time, color = mouse)) +
    geom_line(aes(group = mouse)) +
    geom_point() +
    scale_y_continuous("total close contact in h") +
    scale_x_discrete(limits = c("I1", "A1", "I2", "A2", "I3", "A3", "I4", "A4", "I5")) + 
    scale_color_discrete(name = paste("mice from", batch, cageChange, system))
  
  return(plot)
}

#######################################
# plotting function for comparing file
# input:  
#         
#        
#
# output: 
# effect:
create_joined_table <- function(CC1, CC2, CC3, CC4, mouse_id, batch, stress_condition) {
  
  # filter every tibble to the current mouse id
  filtered_CC1 <-  CC1 %>%
    select(Phase, mouse_id) %>%
    rename(CC1 = mouse_id)
  print(filtered_CC1)
  
  filtered_CC2 <- CC2 %>%
    select(Phase, mouse_id) %>%
    rename(CC2 = mouse_id)
  
  filtered_CC3 <- CC3 %>%
    select(Phase, mouse_id) %>%
    rename(CC3 = mouse_id)
  
  filtered_CC4 <- CC4 %>%
    select(Phase, mouse_id) %>%
    rename(CC4 = mouse_id)
  
  # join the four filtered tibble to one 
  proximity_table_join <- Reduce(function (...) { merge(..., by = "Phase", all = TRUE) },   # Full join of reduced tibbles 1,2,3 and 4
                                 list(filtered_CC1, filtered_CC2, filtered_CC3, filtered_CC4))
  #sort the phase column by int
  proximity_table_join <- proximity_table_join%>%
    arrange(as.integer(sub("[^0-9]", "", Phase)))
  
  print(proximity_table_join)
  message("plotting")
  #data melt
  long_data <- melt(proximity_table_join, id = 'Phase')
  
  #rename columns
  names(long_data) <- c('Phase', 'cageChange', 'time')
  
  #change seconds to hours
  hour_data <- long_data %>%
    mutate(time = as.integer(time)) %>%
    mutate(time = ifelse(time != 0, time / 3600, time))
  
  plot <- ggplot(data = hour_data, aes(x = Phase, y = time, color = cageChange))+
    geom_line(aes(group = cageChange),na.rm=TRUE) +
    geom_point(na.rm = TRUE) +
    scale_color_manual(values = c("CC1" = "aquamarine2", "CC2" = "deepskyblue2", "CC3" = "deepskyblue4", "CC4" = "darkblue"), name = paste(mouse, batch, stress_condition))+
    scale_y_continuous("total close contact in h")+
    scale_x_discrete(limits = c("I1", "A1", "I2", "A2", "I3", "A3", "I4", "A4", "I5"))+
    facet_grid(~cageChange)
  return(plot)
}

############################################################################################
### statistic functions from tobi: ###
########################################################################################

# Function to perform normality test and appropriate statistical test for each variable and phase
testAndPlotVariable <- function(data, value, variableName, phase, sex) {
  #filtering specific phase or sex if needed 
  filteredData <- data %>%
    filter(if('Phase' %in% colnames(data))Phase == phase else TRUE) %>%   
    filter(Sex == sex)
  
  # Convert Group to factor
  filteredData$Group <- as.factor(filteredData$Group)

  # Get unique groups and number of groups
  uniqueGroups <- unique(filteredData$Group)
  numGroups <- length(uniqueGroups)
  
  # Check if the variable is numeric (if not, return Null)
  if (is.numeric(filteredData[[variableName]])) {
    if (numGroups == 2) {
      
      group1 <- uniqueGroups[1]
      group2 <- uniqueGroups[2]
      
      dataGroup1 <- filteredData[[variableName]][filteredData$Group == group1]
      dataGroup2 <- filteredData[[variableName]][filteredData$Group == group2]
      
      # if there are more than 3 columns containing the name of group1 AND group 2 , perform Wilcoxon test
      # WILCOXON
      if (sum(!is.na(dataGroup1)) >= 3 && sum(!is.na(dataGroup2)) >= 3) {
        wilcoxRes <- performWilcoxonTest(dataGroup1, dataGroup2)
        
        #if results of Wilcoxontest are not empty add results to testResults
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
          
          #generate plot of wilcoxontest results
          p <- generatePlot(filteredData, value, variableName, phase, sex)
          
          return(list(testResults = testResults, plot = p, posthocResults = NULL))
        }
      }
    } else {  # numGroups != 2
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
      ## ANOVA
      if (numGroups > 2) {
        
        #every p value has to be not significant(>=0,05)
        if (conNorm$p.value >= 0.05 && resNorm$p.value >= 0.05 && susNorm$p.value >= 0.05) {
          
          anovaTest <- aov(as.formula(paste(variableName, "~ Group")), data = filteredData)
          testResults$Test <- "ANOVA"
          testResults$P_Value <- summary(anovaTest)[[1]][["Pr(>F)"]][1]
          testResults$Significance_Level <- sprintf("%.3f", testResults$P_Value)
          #posthoc for ANOVA
          posthocResultsDf <- performPosthocAnova(filteredData, variableName)
        } else {#Kruskal Wallis
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
      
      #generate plots in p for return
      p <- generatePlot(filteredData, value, variableName, phase, sex)
      
      return(list(testResults = testResults, plot = p, posthocResults = testResultsDf))
    }
  } else {
    cat("Variable", variableName, "is not numeric and will be skipped.\n")
    return(NULL)
  }
}

# Function to generate plots for each variable and phase
generatePlot <- function(processed_data, value, variableName, phase, sex) {
  filteredData <- processed_data #%>%
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
### MATHEMATIC TESTS IN FUNCTIONS: ###

## parametric
# Function to perform post hoc pairwise tests for ANOVA
performPosthocAnova <- function(data, variableName) {
  anovaTest <- aov(as.formula(paste(variableName, "~ Group")), data = data)
  pairwiseResults <- rstatix::pairwise_t_test(as_data_frame(data), formula = as.formula(paste(variableName, "~ Group")),
                                     p.adjust.method = "bonferroni")
  # Remove duplicate comparisons
  pairwiseResults <- pairwiseResults[!duplicated(pairwiseResults[, c("group1", "group2")]), ]
  pairwiseResults$GroupComparison <- paste(pairwiseResults$group1, "vs.", pairwiseResults$group2)
  return(pairwiseResults)
}

# Function to perform post hoc pairwise tests for Kruskal-Wallis
performPosthocKruskal <- function(data, variableName) {
  kruskalTest <- kruskal.test(as.formula(paste(variableName, "~ Group")), data = data)
  pairwiseResults <- dunn_test(data, formula = as.formula(paste(variableName, "~ Group")),
                               p.adjust.method = "holm")
  pairwiseResults$GroupComparison <- paste(pairwiseResults$group1, "vs.", pairwiseResults$group2)
  return(pairwiseResults)
}

# non parametric
# Function to perform Wilcoxon rank-sum test
performWilcoxonTest <- function(dataGroup1, dataGroup2) {
  if (length(dataGroup1) >= 3 && length(dataGroup2) >= 3) {
    return(wilcox.test(dataGroup1, dataGroup2))
  } else {
    return(NULL)
  }
}

############################################################################################
### saving functions: ###
########################################################################################
# input:  
#         
#        
#
# output: 
# effect:
