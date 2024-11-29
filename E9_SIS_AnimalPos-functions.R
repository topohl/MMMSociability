
## 11/2023
## Anja Magister
## ANALYSIS OF ANIMAL POSITIONS - FUNCTIONS ##
##

########### Functions for analyzing shannon entropy of animal positions ###########

analyze_batch <- function(batch, change, working_directory, uniqueSystems, phases) {
  # Define the sex based on the batch
  sex <- ifelse(batch %in% c("B3", "B4", "B6"), "female", "male")
  
  # Read the overall data for the batch and cage change
  filename <- paste0("E9_SIS_", batch, "_", change, "_AnimalPos")
  csvFilePath <- file.path(working_directory, "preprocessed_data", paste0(filename, "_preprocessed.csv"))
  
  if (!file.exists(csvFilePath)) {
    warning(paste("File not found:", csvFilePath))
    return(NULL)
  }
  
  overallData <- read_delim(csvFilePath, delim = ",", show_col_types = FALSE) %>% as_tibble()
  
  # Initialize results
  cagePosProb <- tibble()
  cagePosEntropy <- tibble()
  animalPosEntropy <- tibble()
  
  # Loop through each system
  for (system_id in uniqueSystems) {
    print(system_id)
    
    # Filter data for the current system
    systemData <- overallData %>%
      filter(System == system_id) %>%
      as_tibble()
    
    # Identify unique animals in the system
    animal_ids <- unique(systemData$AnimalID)
    system_complete <- length(animal_ids) >= 4
    
    # Pad with NA if the system is incomplete
    while (length(animal_ids) < 4) {
      animal_ids <- append(animal_ids, NA)
    }
    
    # Loop through each phase
    for (phase in phases) {
      print(phase)
      
      # Handle special case for `CC4`
      if (change == "CC4") {
        active_phases_number <- c(1, 2)
        inactive_phases_number <- c(2)
      } else {
        active_phases_number <- c(1, 2, 3)
        inactive_phases_number <- c(1, 2)
      }
      
      number_of_phases <- ifelse(phase == "Active", active_phases_number, inactive_phases_number)
      print(number_of_phases)
      
      # Loop through each phase number
      for (nr in number_of_phases) {
        print(paste0(batch, ", System: ", syste_id, ", ", change, ", ", phase, " phase nr: ", nr))
        
        # Filter data for the current phase
        system_data_phase <- system_data %>%
          filter(ConsecActive == ifelse(phase == "Active", nr, 0)) %>%
          filter(ConsecInactive == ifelse(phase == "Inactive", nr, 0)) %>%
          as_tibble()
        
        # Initialize data structures
        animal_list <- initialize_animal_list(animal_ids)
        cage_prob_list <- initialize_cage_prob_list()
        animal_prob_tibble <- initialize_animal_prob_tibble(animal_ids)
        
        # While loop through the system_data_phase rows
        animal_list <- process_rows(system_data_phase, animal_list, cage_prob_list, animal_prob_tibble, system_complete)
        
        # Calculate probabilities and entropy
        cagePosProb <- update_cage_probs(cage_prob_list, cagePosProb, batch, system_id, change, phase, nr)
        cagePosEntropy <- calculate_cage_entropy(cagePosProb, cagePosEntropy, batch, system_id, change, phase, nr, sex)
        animalPosEntropy <- calculate_animal_entropy(animal_prob_tibble, animalPosEntropy, batch, system_id, change, phase, nr, animal_ids, sex)
      }
    }
  }
  
  # Return all results
  list(
    cagePosProb = cagePosProb,
    cagePosEntropy = cagePosEntropy,
    animalPosEntropy = animalPosEntropy
  )
}

########### Support function to initialise animal list ###########
initialize_animal_list <- function(animal_ids) {
  list(
    animalOne = list(name = "", time = "", position = 0),
    animalTwo = list(name = "", time = "", position = 0),
    animalThree = list(name = "", time = "", position = 0),
    animalFour = list(name = "", time = "", position = 0),
    tempData = list(secTemp = 0, lineTemp = 0)
  )
}

########### Support function to initialise cage probability list ###########
initialize_cage_prob_list <- function() {
  lapply(1:8, function(i) c(i, 0, 0, 0))
}

########### Support function to initialise animal probability tibble ###########
initialize_animal_prob_tibble <- function(animal_ids) {
  tibble(
    AnimalID = rep(animal_ids, each = 8),
    Position = rep(1:8, length.out = 32),
    Seconds = 0,
    SumPercentage = 0,
    Prob = 0
  )
}

########### Support function to process rows of data ###########
process_rows <- function(system_data_phase, animal_list, cage_prob_list, animal_prob_tibble, system_complete) {
  timeTemp <- animal_list[[1]][["time"]]
  lineTemp <- 5
  theEnd <- nrow(system_data_phase) + 1
  
  while (lineTemp < theEnd) {
    old_animal_list <- animal_list
    animal_list <- update_animal_list(animal_ids, animal_list, system_data_phase, timeTemp, lineTemp)
    secTemp <- animal_list[["tempData"]][["secTemp"]]
    
    if (system_complete) {
      cage_prob_list <- check_cage_prob(old_animal_list, animal_list, cage_prob_list, secTemp)
    }
    animal_prob_tibble <- check_animal_prob(old_animal_list, animal_list, animal_prob_tibble, secTemp)
    lineTemp <- animal_list[["tempData"]][["lineTemp"]]
    timeTemp <- animal_list[[1]][["time"]]
  }
  
  animal_list
}

## Preprocessing
# Define function to preprocess a single file
preprocess_file <- function(batch, change, exclAnimals) {
  filename <- paste0("E9_SIS_", batch, "_", change, "_AnimalPos")
  csvFilePath <- file.path(getwd(), "raw_data", batch, paste0(filename, ".csv"))
  
  if (!file.exists(csvFilePath)) {
    warning(paste("File", csvFilePath, "does not exist. Skipping to next file."))
    return(NULL)
  }
  
  # Read CSV file
  data <- as_tibble(read_delim(csvFilePath, delim = ";", show_col_types = FALSE))
  
  # Preprocessing steps
  data <- data %>%
    select(-c(RFID, AM, zPos)) %>%
    mutate(DateTime = as.POSIXct(DateTime, format = "%d.%m.%Y %H:%M:%S", tz = "UTC")) %>%
    separate(Animal, into = c("AnimalID", "System"), sep = "[-_]") 
  
  # Define Position Mapping Table
  tblPosition <- tibble(
    PositionID = 1:8, 
    xPos = c(0, 100, 200, 300, 0, 100, 200, 300), 
    yPos = c(0, 0, 0, 0, 116, 116, 116, 116)
  )
  
  # Add PositionID to the data
  data <- data %>%
    rowwise() %>%
    mutate(PositionID = find_id(xPos, yPos, lookup_tibble = tblPosition)) %>%
    select(DateTime, AnimalID, System, PositionID) %>%
    filter(!AnimalID %in% exclAnimals) %>%
    arrange(DateTime)
  
  # Add phase information and markers
  data <- data %>%
    addPhaseTransitionMarkers() %>%
    mutate(Phase = ifelse(format(DateTime, "%H:%M", tz = "UTC") >= "18:30" | 
                          format(DateTime, "%H:%M", tz = "UTC") < "06:30", 
                          "Active", "Inactive")) %>%
    mutate(ConsecActive = 0, ConsecInactive = 0) %>%
    consecPhases() %>%
    addDayMarkerRows()
  
  # Save the preprocessed data
  outputFileDir <- file.path(outputDir, paste0(filename, "_preprocessed.csv"))
  write_csv(data, outputFileDir)
  
  message(paste("File saved at", outputFileDir))
}

##############################################################################################################
# input:  x_Pos - x coordinate of cage
#         y_Pos - y coordinate of cage
#         lookup_tibble - table of defined position IDs for every possible coordinate couple
#
# output: Position Id
# effect: finds corresponding PositionID from two coordinates(x_Pos, y_pos) in lookup_tibble and returns ID 

#' Find Position ID based on x and y coordinates
#'
#' This function standardizes animal positions by converting non-standard 
#' positions into a predefined standard format and then looks up the 
#' corresponding Position ID from a provided lookup table.
#'
#' @param x_Pos Numeric value representing the x-coordinate of the position.
#' @param y_Pos Numeric value representing the y-coordinate of the position.
#' @param lookup_tibble A tibble containing the lookup table with columns 
#'        `xPos`, `yPos`, and `PositionID`.
#'
#' @return The Position ID corresponding to the standardized x and y coordinates.
#'         If no match is found, the function returns NA.
#'
#' @details The function first standardizes the x and y coordinates:
#'          - y_Pos is set to 0 if it is less than 116, otherwise it is set to 116.
#'          - x_Pos is set to 0 if it is less than 100, 100 if it is between 100 and 199,
#'            200 if it is between 200 and 299, and 300 if it is 300 or greater.
#'          The function then searches for the standardized position in the lookup table.
#'          If a match is found, the corresponding Position ID is returned.
#'          If no match is found, the function prints the standardized coordinates and 
#'          returns NA.
#'
#' @examples
#' # Example usage:
#' lookup_tibble <- tibble::tibble(
#'   xPos = c(0, 100, 200, 300),
#'   yPos = c(0, 116),
#'   PositionID = 1:4
#' )
#' find_id(150, 120, lookup_tibble)
#'
#' @export
find_id <- function(x_Pos, y_Pos, lookup_tibble) {
  
  # This script contains functions to standardize animal positions.
  # It converts non-standard positions into a predefined standard format.
  # y value
  if(y_Pos<116){y_Pos <- 0}
  else if(y_Pos>=116){y_Pos <- 116}
  
  # X value
  if(x_Pos<100){x_Pos <- 0}
  else if(x_Pos<200){x_Pos <- 100}
  else if(x_Pos<300){x_Pos <- 200}
  else if(x_Pos>=300){x_Pos <- 300}
  
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
# Function: addPhaseTransitionMarkers
# 
# Description:
#   Inserts additional rows into the dataset to mark phase transitions at specific times (6:30 AM and 6:30 PM). 
#   These markers are essential for tracking the start and end of active/inactive phases for each animal in each system.
#
# Input:
#   data - A tibble containing the experiment data, including timestamped observations for various subjects across systems.
#
# Output:
#   A tibble with added rows corresponding to the phase transition times. The dataset is sorted by DateTime, 
#   and the additional rows provide a start value for every subject at the beginning of each phase.
#
# Effects:
#   - Inserts rows at 6:30 AM and 6:30 PM for each date in the dataset.
#   - Ensures there is a record for each subject in every system at these critical phase change times.
#
################################################################################################################

addPhaseTransitionMarkers <- function(data) {
  
  #save dates of experiment(the date of the days)
  onlyDates <- as.Date(data$DateTime)
  uniqueDates <- unique(onlyDates)
  #print(uniqueDates)
  
  #uniqueExperimentDates <- unique(as.Date(data$DateTime))
  
  #save existing systems
  uniqueSystems <- unique(data$System)
  #print(uniqueSystems)
  
  
  for(i in 1:length(uniqueDates)){
    date <- uniqueDates[i]
  #  
  #  #print(date)

  #for (experimentDate in uniqueExperimentDates) {  
    
    #define marker points depending on actual date
    # "GMT" means "UTC"
    # Define the specific phase transition markers for the current date
    phaseMarkers <- list(
      endOfActivePhaseAM = as.POSIXct(paste(date, "06:29:59"), tz = "GMT"),
      startOfInactivePhaseAM = as.POSIXct(paste(date, "06:30:00"), tz = "GMT"),
      endOfInactivePhasePM = as.POSIXct(paste(date, "18:29:59"), tz = "GMT"),
      startOfActivePhasePM = as.POSIXct(paste(date, "18:30:00"), tz = "GMT")
    )
    
    #we have 2 early markers and 2 late markers
    for(phaseMarker in phaseMarkers){
     
      # filter rows which are earlier than the marker
      filtered_time <- data %>%
        as_tibble()%>%
        filter(DateTime < phaseMarker)
      
      #we need a marker for every existing system and for every animal in every system
      for(system in uniqueSystems){
        
        #print(system)
        
        # filter rows from specific system 
        filtered_system <- filtered_time%>%
          as_tibble()%>%
          filter(System==system)
        
        #check which animal exists in this system
        animal_names <- unique(filtered_system$AnimalID)
        
        for(animal in animal_names){
          # filter specific animal row from system with nearest DateTime to 18:30:00 or 06:30:00 but still earlier
          filtered_row <- filtered_system%>%
            as_tibble()%>%
            filter(AnimalID==animal)%>%
            slice_max(order_by = DateTime)
          
          #print(filtered_row)
          
          #if time on this day is recorded
          if(length(filtered_row) != 0){
            #use copied row as new one and change the date to marker date
            filtered_row <- filtered_row%>%
              mutate(DateTime = as.POSIXct(phaseMarker, tz="UTC"))
            
            #add new row to data
            data <- add_row(data,filtered_row)
          }
        
        
        }
        
      }
      
    }
    
  }
  
  #sort tibble again by dateTime to bring new entries to correct position
  data <- data%>%
    as_tibble()%>%
    arrange(., DateTime)
  
  
  #remove last 20 rows because they are 20 unnecessary marker rows of the last day and the late marker
  data <- data %>% filter(row_number() <= n()-20)
  
  return(data)
}



################################################################################################################
# Function: addDayMarkerRows
#
# Description:
#   Adds extra rows of RFID information at the start of each new day (midnight) for each animal in each system.
#   This is essential for tracking the start of each day for every animal in every system.
#
# Input:
#   data - A tibble containing the experiment data, including timestamped observations for various subjects across systems.
#
# Output:
#   A tibble with added rows corresponding to the start of each new day. The dataset is sorted by DateTime, 
#   and the additional rows provide a start value for every subject at the beginning of each day.
#
# Effects:
#   - Inserts rows at midnight for each date in the dataset.
#   - Ensures there is a record for each subject in every system at the start of each day.
#
################################################################################################################

addDayMarkerRows <- function(data) {
  
  # Extract unique dates from the data
  uniqueDates <- unique(as.Date(data$DateTime))
  
  # Extract unique systems from the data
  uniqueSystems <- unique(data$System)
  
  for (date in uniqueDates) {
  
  # Define the marker point for the start of the new day (midnight)
  marker <- as.POSIXct(paste(date, "00:00:00"), tz = "GMT")
  
  # Filter rows which are earlier than the marker
  filtered_time <- data %>%
    filter(DateTime < marker)
  
  for (system in uniqueSystems) {
    
    # Filter rows for the specific system
    filtered_system <- filtered_time %>%
    filter(System == system)
    
    # Identify unique animals in the system
    animal_names <- unique(filtered_system$AnimalID)
    
    for (animal in animal_names) {
    
    # Filter the specific animal row with the nearest DateTime to midnight but still earlier
    filtered_row <- filtered_system %>%
      filter(AnimalID == animal) %>%
      slice_max(order_by = DateTime)
    
    # If time on this day is recorded
    if (nrow(filtered_row) != 0) {
      
      # Use the copied row as a new one and change the date to the marker date
      filtered_row <- filtered_row %>%
      mutate(DateTime = as.POSIXct(marker, tz = "UTC"))
      
      # Add the new row to the data
      data <- add_row(data, filtered_row)
    }
    }
  }
  }
  
  # Sort the tibble again by DateTime to bring new entries to the correct position
  data <- data %>%
  arrange(DateTime)
  
  return(data)
}

##############################################################################################################
# input:  data  - a tibble with data that we want to change 
#         
#        
#
# output: data
# effect: - edits the columns of the phases(ConsecActive, ConsecInactive)
#         - counts seperately how many active and inactive phases happened in that data, adds the counted int into the column
#         - when the active phase is happening , the ConsecInactive is 0 and vice versa
consecPhases <- function(data){
  
  #initialize Phase counter
  active_phases <- 0
  inactive_phases <- 0
  
  for(i in 1:nrow(data)){
    #select row from data
    current_row <- data[i,]
    
    #print(i)
    
    #first row
    if(i==1){
      
      if(current_row$Phase == "Active"){
        current_row$ConsecActive <- 1
        current_row$ConsecInactive <- 0
        active_phases <- active_phases+1
      }else{    #current_row$Phase == "Inactive"
        current_row$ConsecActive <- 0
        current_row$ConsecInactive <- 1
        inactive_phases <- inactive_phases+1
      }
    }
    else{#i>1
      
      #previous row from data
      previous_row <- data[(i-1),]
      
      if(current_row$Phase != previous_row$Phase){
        if(current_row$Phase == "Active"){
          active_phases <- active_phases+1  # change means the counter has to set higher
          current_row$ConsecActive <- active_phases
          current_row$ConsecInactive <-  0
        }else{    #current_row$Phase == "Inactive"
          inactive_phases <- inactive_phases+1
          current_row$ConsecActive <- 0
          current_row$ConsecInactive <- inactive_phases
        }
      }else{  #current_row$Phase == previous_row$Phase
        current_row$ConsecActive <- previous_row$ConsecActive      #values stay the same number as before
        current_row$ConsecInactive <- previous_row$ConsecInactive
      }
      
      #write new information back into data
      data[(i-1),] <- previous_row
    }
    
    #write new information back into data
    data[i,] <- current_row  
  }
  
  return (data)
}

##############################################################################################################
# input:  system_animal_names -  the ids/names of the animal in the current system
#         data -  the tibble with the data of the current system
#         animal_list - the list that contains the information of the animal ID and the current time and position
#         
#
# output: animal_list - updated
# effect: - the animal-list will be entered here for the first time, information gained from the data
#         - takes the first entry of every animal in data and copies the time and position that we know the start position of each animal

# find the FIRST TIME where animal is tracked in the cage
# aka first value of animal in data_final
find_first_pos_and_time <- function(system_animal_ids, data, animal_list){
  
  # create empty vector with 4 variables
  times_vec <- rep(NA, times=4)
  
  print(system_animal_ids)
  for (i in 1:length(system_animal_ids)){ #i=1-4
    
    #define current animal_name, first_position and first_time
    animal_id <- system_animal_ids[[i]]
    if(is.na(animal_id)){
      animal_id <- paste0("lost_",i)
      first_time <- 0        #time value should be adapted at the end of this function
      first_position <- (-1)  #add non existing position for non existing animal
    }
    else{
      #search first entry in whole data
      first_entry <- data%>%
        filter(AnimalID == animal_id)%>%
        slice(1) #first row
      
      #cat("first entry: \n")
      #print(first_entry)
      
      first_time <- first_entry$DateTime
      
      first_position <- first_entry$PositionID
      
      #cat("first position: ", first_position)
      
    }
    
    
    #write name, position and time into animal_list
    animal_list[i][[1]] <- animal_id
    #print(animal_name)
    
    animal_list[[i]][[3]] <- first_position
    #print(first_position)
    
    animal_list[[i]][[2]] <- first_time
    #print(first_time)
    
    #enter time in times vec fot later to compare
    times_vec[i] <- first_time
    
  }
  
  #check, if all times are similar
  #if not, change them all to the first shared time
  if(length(unique(times_vec)) != 1){
    latest_time <- max(times_vec)
    #print(latest_time)
    for(i in 1:4){
      animal_list[[i]][[2]] <- latest_time
    }
  }
  
  #print("animal_list inside the first pos function: \n")
  #print(animal_list)
  
  
  
  return(animal_list)
}


##############################################################################################################
# input:  old_animal_list - contains old positions of the animal from the last entered time in tibble until the second before the current time
#         new_animal_list - contains the new updated positions for the current second(time)
#         count_closeness_list - list for the results, contains the amount of seconds with animal in social contact
#         secTemp - time difference between new and old animal-list
#
# output: count_closeness_list
# effect: - reads the position information from the old positions and calculates the second difference to the actual time, adds seconds to results-list
#         - also reads the new positions and adds one second for every position into the results-list

# function to check for closeness
# compare every sublist(4)to each other
check_closeness <- function(old_animal_list, new_animal_list, count_closeness_list, secTemp){
  
  # compare the third value of every couple(which is the position of the animal)
  # if the position is the same, save in count_closeness_list list
  for (i in 1:4) {
    for (j in i:4) {
      #comparisation of old_animal_list for last (secTemp-1)-missing seconds
      if(old_animal_list[[i]][[3]]==old_animal_list[[j]][[3]]){
        count_closeness_list[[i]][[j]] <- count_closeness_list[[i]][[j]]+(secTemp-1)  #secTemp contains x-1 seconds of old positions and one sec of new pos

        #enter new contact second in closeness list if the animal are at the same position and are two different individuals
        if(j!=i){
          count_closeness_list[[j]][[i]] <- count_closeness_list[[j]][[i]]+(secTemp-1)
        }
      }

      #comparisation of new_animal_list for first new second
      if(new_animal_list[[i]][[3]]==new_animal_list[[j]][[3]]){
        count_closeness_list[[i]][[j]] <- count_closeness_list[[i]][[j]]+1

        #enter new contact second in closeness list if the animal are at the same position and are two different individuals
        if(j!=i){
          count_closeness_list[[j]][[i]] <- count_closeness_list[[j]][[i]]+1
        }
      }
    }
  }
  
  # return updated list of animal that are close to each other
  return(count_closeness_list)
}

##############################################################################################################
# input:  old_animal_list - contains old positions of the animal from the last entered time in tibble until the second before the current time
#         new_animal_list - contains the new updated positions for the current second(time)
#         count_position_list - list for the results, contains the amount of seconds a animal(or multiple animal)were standing on a position
#         secTemp - time difference between new and old animal-list
#         
#
# output: count_position_list
# effect:

# function to check for position


check_position <- function(old_animal_list, new_animal_list, count_position_list, secTemp){
  
  #vectors to check which positions were used in this time period
  old_positions <- c()
  new_positions <- c()
  
  #for every animal in the animal list
  for (i in 1:4) {
    old_pos <- as.numeric(old_animal_list[[i]][[3]])
    #print("old_pos")
    #print(old_pos)
    new_pos <- as.numeric(new_animal_list[[i]][[3]])
    #print("new_pos")
    #print(new_pos)
    
    #count used position only one time per second, even if multiple animal are on it
    #check if position is already added to list from another animal
    #(is actual position already in the vector of used positions)
    
    if(!old_pos %in% old_positions & old_pos!=(-1)){
      #print("add old pos")
      
      #add new seconds to count_position_list
      count_position_list[[old_pos]][[2]] <- count_position_list[[old_pos]][[2]]+(secTemp-1)
      
      #add old position into vector of used positions
      old_positions <- append(old_positions, old_pos)
    }
    if(!new_pos %in% new_positions& new_pos!=(-1)){
      #print("add new pos")
      
      #add new seconds to count_position_list
      count_position_list[[new_pos]][[2]] <- count_position_list[[new_pos]][[2]]+1
      
      #add new position into vector of used positions
      new_positions <- append(new_positions, new_pos)
    }
    
  }
  
  # return updated list of animal that are close to each other
  return(count_position_list)
}



##############################################################################################################
# input:  old_animal_list, new_animal_list, animal_ids
#         
#        
#
# output: total_closeness_list
# effect: calculate total amount of social(time in closer contact) seconds per animal 


check_total_closeness <- function(old_animal_list, new_animal_list, total_closeness_list, secTemp){

  for(animal in 1:4){
    
    other_animal <- c(1,2,3,4)
    #remove current animal index from other_animal
    other_animal <- other_animal[-animal]
    
    #get actual animal id from animal_list
    animal_id <- old_animal_list[[animal]][[1]]
    
    ##OLD; REMAINING SECONDS
    #animal_list[[x]][[3]]is looking for position
    #check if current animal was on the same position as one or multiple animal in the time between last entry and current entry(secTemp, old_animal_list)
    if(old_animal_list[[animal]][[3]]==old_animal_list[[other_animal[1]]][[3]] ||
       old_animal_list[[animal]][[3]]==old_animal_list[[other_animal[2]]][[3]] ||
       old_animal_list[[animal]][[3]]==old_animal_list[[other_animal[3]]][[3]] )
    {
      #add remaining closeness seconds(since the last entry, see secTemp) into total_closeness_list for current animal
      if(total_closeness_list[[animal]][[1]]==animal_id){#double check if id is the right one
        total_closeness_list[[animal]][[2]] <- as.numeric(total_closeness_list[[animal]][[2]])+(secTemp-1)
      }else{throw("animal ids are not the same in check_total_closeness")}
    }
    
    ##FIRST NEW SECOND
    #check if current animal was on the same position as one or multiple animal in the current second
    if(new_animal_list[[animal]][[3]]==new_animal_list[[other_animal[1]]][[3]] ||
       new_animal_list[[animal]][[3]]==new_animal_list[[other_animal[2]]][[3]] ||
       new_animal_list[[animal]][[3]]==new_animal_list[[other_animal[3]]][[3]] )
    {
      #add closeness second into total_closeness_list for current animal
      if(total_closeness_list[[animal]][[1]]==animal_id){#double check if id is the right one
        total_closeness_list[[animal]][[2]] <- as.numeric(total_closeness_list[[animal]][[2]])+1
      }else{throw("animal ids are not the same in check_total_closeness")}
    }
    
  }
  
  
  #browser()
  return(total_closeness_list)
}
##############################################################################################################
# input: old_animal_list,new_animal_list,count_movement_list, secTemp
#         
#        
#
# output: count_movement_list
# effect: 

# 
check_movement <- function(old_animal_list,new_animal_list,count_movement_list, secTemp){
  #variable check if anyone moved at all
  moves <- 0
  
  #for every animal in the animal list
  for (i in 1:4) {
    old_pos <- as.numeric(old_animal_list[[i]][[3]])
    #print("old_pos")
    #print(old_pos)
    new_pos <- as.numeric(new_animal_list[[i]][[3]])
    #print("new_pos")
    #print(new_pos)
    
    
    #check if animal changed place in this second, how long is not important
    #also checking for lost_animal with position -1(not valuable)
    if(old_pos!=(-1) & new_pos!=(-1) & old_pos!=new_pos){
      #enter one movement to current animal
      count_movement_list[[i]][[2]] <- as.numeric(count_movement_list[[i]][[2]])+1
      #enter one movement to whole system
      count_movement_list[[5]][[2]] <- as.numeric(count_movement_list[[5]][[2]])+1
      
      #add move to moves
      moves <- moves+1
    }
    
    
  }
  
  if(moves==0){print("check_movement: no one moved in this row :)")}
  
  # return updated list of animal that are close to each other
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
sec_shift <- function( old_time){
  #put one second on top of old_time
  new_time <- old_time%>%
    as.numeric()%>%
    +1%>%
    as.character()
  #cat("old time: ", old_time, "\n")
  #cat("new time: ", new_time, "\n")
  return(new_time)
}

##############################################################################################################
# input:  system_animal_names -  the ids/names of the animal in the current system
#         animal_list - the list containing the information of the animal ID and the current time and position
#         data -  the tibble with the data of the current system
#         time - the time that we are looking at, can be written in several lines bc we have several animal
#         line - a counter for the actual line in our data that we are looking at
#
# output: 
# effect: - update animal_list(if its possible) and return it
#         - takes next entry in data and updates the new information in the animal list
#         - similarity to find_first_pos_and_time
update_animal_list <- function(system_animal_ids, animal_list, data, time, line){
  
  #next_second <- sec_shift(time)
  new_time <- as.numeric(data[line,"DateTime"])

  # write sec difference between new and old time into secTemp
  animal_list[["tempData"]][["secTemp"]] <- new_time-as.numeric(time)
  
  # write new time into every animal information
  for(i in 1:4){animal_list[[i]][[2]] <- new_time}

  # while line(and especially the next lines) is still same time
  while(as.numeric(data[line,"DateTime"])==new_time){
    # write new position into special animal
    for(i in 1:4){
      if(animal_list[[i]][[1]]==as.character(data[line,"AnimalID"])){animal_list[[i]][[3]] <- as.numeric(data[line,"PositionID"])}
    }
    #if line is not the last line, check next line, else break while loop
    if(line==nrow(data)){
      line <- line+1
      break
    }
    #continue with the while condition
    line <- line+1
  }
  
  # write new line into animal_list
  animal_list[["tempData"]][["lineTemp"]] <- line
  
  return(animal_list)
}


########################################################################################
## functions for comparing-file

# input:  
#   rank_tibble - the tibble containing the ranks
#   vector - the vector to be ranked
#   system - the system identifier
#   change - the cage change identifier
#   batch - the batch identifier
#   hours_column_name - the name of the column containing the hours
#   rank_column_name - the name of the column containing the ranks
#
# output: 
#   rank_tibble - the updated tibble with ranks
#
# effect:
#   Computes the ranks for the given vector and updates the rank_tibble with the ranks for the corresponding system, change, and batch.

compute_rank <- function(rank_tibble, vector, system, change, batch, hours_column_name, rank_column_name) {
  
  # Sort the total hour system vector by size, with the largest value getting rank 1
  vector <- sort(vector)
  
  # Enter ranks in rank tibble
  for (i in 1:length(vector)) {
  # Enter rank in corresponding line in rank tibble
  rank_tibble <- rank_tibble %>%
    mutate({{rank_column_name}} := ifelse((Batch == batch) & (change == change) & (System == system) & (.data[[hours_column_name]] == vector[i]), i, .data[[rank_column_name]])) # New R version & instead of &&
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

translate_rank_in_score_vec <- function(rank_vec){
  #change rank into score
  #rank = score
  # 1   =   4
  # 2   =   3
  # 3   =   2
  # 4   =   1
  
  for (i in 1:length(rank_vec)){
    #change rank into score 
    score <- ifelse(rank_vec[i]==1, 4, ifelse(rank_vec[i]==2, 3, ifelse(rank_vec[i]==3, 2, ifelse(rank_vec[i]==4, 1, NA))))
    #overwrite rank with score
    rank_vec[i] <- score
  }
  
  score_vec <- rank_vec
  return(score_vec)
}  


############################################################################################
## FUNCTIONS FOR ANALYZING SHANNON ENTROPY METRICS ##
############################################################################################

# Function: check_cage_prob
# -------------------------
# Description:
# Updates a probability list (`cage_prob_list`) by calculating the positional 
# probabilities of animals in a cage across time intervals.
#
# Input:
# - old_animal_list: List containing the previous positions of animals.
# - new_animal_list: List containing the updated positions of animals.
# - cage_prob_list: List storing positional probabilities and counts for each cage position.
# - secTemp: Number of seconds for which new positions should be considered.
#
# Output:
# - Updated `cage_prob_list` with recalculated probabilities and counts for each position.
#
# Effect:
# The function modifies the `cage_prob_list` by incorporating data from both old and new
# animal position lists, weighted by the time intervals specified.

check_cage_prob <- function(old_animal_list, new_animal_list, cage_prob_list, secTemp) {
  # Initialize position counters for old and new positions (8 possible positions in the cage)
  old_cage_positions <- integer(8) # Vector of zeros, length 8
  new_cage_positions <- integer(8) # Vector of zeros, length 8

  # Iterate over the first 4 animals in the lists to update position counters
  for (i in 1:4) {
    # Retrieve position indices for old and new animal positions
    old_pos <- as.numeric(old_animal_list[[i]][[3]])
    new_pos <- as.numeric(new_animal_list[[i]][[3]])

    # Increment the position counters at the respective indices
    old_cage_positions[old_pos] <- old_cage_positions[old_pos] + 1
    new_cage_positions[new_pos] <- new_cage_positions[new_pos] + 1
  }

  # Normalize position counters by dividing by 4 (total number of animals)
  old_cage_positions <- old_cage_positions / 4
  new_cage_positions <- new_cage_positions / 4

  # Update the cage_prob_list with the computed positional probabilities
  for (i in 1:8) {
    # Increment the total count for the current position in the old list
    cage_prob_list[[i]][[2]] <- cage_prob_list[[i]][[2]] + 1
    # Add the probability of this position from the old positions
    cage_prob_list[[i]][[3]] <- cage_prob_list[[i]][[3]] + old_cage_positions[i]

    # Update probabilities for (secTemp - 1) additional seconds using the old positions
    for (second in 1:(secTemp - 1)) {
      cage_prob_list[[i]][[2]] <- cage_prob_list[[i]][[2]] + 1
      cage_prob_list[[i]][[3]] <- cage_prob_list[[i]][[3]] + old_cage_positions[i]
    }
  }

  # Return the updated probability list
  return(cage_prob_list)
}

# Function: check_animal_prob
# ---------------------------
# Description:
# Updates the `animal_prob_tibble` to track the time spent by animals at specific positions 
# and their cumulative probabilities based on positional changes over a given time interval.
#
# Input:
# - old_animal_list: List containing the previous positions of animals.
# - new_animal_list: List containing the updated positions of animals.
# - animal_prob_tibble: Tibble containing information about animal positions, including:
#   - AnimalID: Unique identifier for each animal.
#   - Position: The cage position for the animal.
#   - SumPercentage: Cumulative probability for each animal at each position.
#   - Seconds: Total time spent by each animal.
# - secTemp: Number of seconds to update for new positions.
#
# Output:
# - Updated `animal_prob_tibble` with recalculated probabilities and time data.
#
# Effect:
# Modifies `animal_prob_tibble` to reflect updated positional probabilities and 
# time spent for each animal, based on their old and new positions.

check_animal_prob <- function(old_animal_list, new_animal_list, animal_prob_tibble, secTemp) {
  for (i in 1:4) { # Loop through each of the 4 tracked animals
    
    # Skip updating if the current animal is not tracked (AnimalID is NA)
    if (is.na(animal_ids[i])) {
      next
    }

    # Retrieve old and new position indices for the current animal
    old_pos <- as.numeric(old_animal_list[[i]][[3]])
    new_pos <- as.numeric(new_animal_list[[i]][[3]])

    # Locate the row in `animal_prob_tibble` matching the animal's ID and position
    row_old <- which(animal_prob_tibble$AnimalID == animal_ids[i] & animal_prob_tibble$Position == old_pos)
    row_new <- which(animal_prob_tibble$AnimalID == animal_ids[i] & animal_prob_tibble$Position == new_pos)

    # Increment the cumulative probability for the old position by 1
    animal_prob_tibble[["SumPercentage"]][[row_old]] <- animal_prob_tibble[["SumPercentage"]][[row_old]] + 1

    # Increment the cumulative probability for the new position by (secTemp - 1)
    animal_prob_tibble[["SumPercentage"]][[row_new]] <- animal_prob_tibble[["SumPercentage"]][[row_new]] + (secTemp - 1)

    # Update the total seconds for the current animal across all positions
    animal_prob_tibble <- animal_prob_tibble %>%
      mutate(Seconds = ifelse(AnimalID == animal_ids[i], Seconds + secTemp, Seconds))
  }

  # Return the updated tibble with revised positional probabilities and time data
  return(animal_prob_tibble)
}

# Function: calc_shannon_entropy
# ------------------------------
# Description:
# Calculates the Shannon entropy of a given probability vector (`prob_vec`), 
# which quantifies the uncertainty or diversity in a system.
#
# Input:
# - prob_vec: A numeric vector of length 8, representing probabilities for 8 positions.
#
# Output:
# - A numeric value representing the calculated Shannon entropy.
#
# Effect:
# Throws an error if the input vector does not have exactly 8 elements, ensuring data integrity.
# Handles cases where probabilities are 0 to avoid undefined log2(0) operations.

calc_shannon_entropy <- function(prob_vec) {
  
  # Validate the input vector: ensure it has exactly 8 elements
  if (length(prob_vec) != 8) {
    print(prob_vec)
    stop("Error in Shannon entropy calculation: prob_vec is not the correct size.")
  }

  # Initialize the Shannon entropy value
  shannon_entropy <- numeric()

  # Compute Shannon entropy using the formula: -Σ(p * log2(p))
  for (i in 1:8) {
    if (prob_vec[i] == 0) {
      # If probability is 0, add 0 to the sum (log2(0) is undefined)
      shannon_entropy <- sum(shannon_entropy, 0)
    } else {
      # Add the term p * log2(p) to the sum
      shannon_entropy <- sum(shannon_entropy, prob_vec[i] * log2(prob_vec[i]))
    }
  }

  # Negate the result to obtain the Shannon entropy
  shannon_entropy <- -shannon_entropy

  # Return the calculated Shannon entropy
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
generateHeatMapCloseness <- function(count_closeness_list, batch, change, systemNum, animal_ids, phase, nr){
  # calculate second entrys to hour entrys
  count_closeness_list_hours <- lapply(count_closeness_list, function(x) ifelse(x!=0,x/3600,x))
  
  
  
  # convert list of lists into a matrix
  matrix_data <- do.call(rbind, count_closeness_list_hours)
  
  # names of the animal Ids
  dimnames(matrix_data) <- list(animal_ids, animal_ids)
  
  # melt the data, means create values combinations out of the matrix
  data_melt <- melt(matrix_data, as.is = TRUE, value.name = "hours")                                          # Reorder data
  #head(data_melt) 
  
  #print("System:")
  #print(systemNum)
  
  #create the plot
  ggp <- ggplot(data_melt, aes(Var1, Var2)) +                                 # Create heatmap with ggplot2
    geom_tile(aes(fill = hours))+
    scale_fill_gradientn(colors = c("lightblue","blue","darkblue"), 
                         limits = c(0, 15), 
                         breaks = c(2, 4, 6, 8, 10, 12, 14, ifelse(max(data_melt$hours)<15,15,warning("Higher scale required in heatmap for ", systemNum))),
                         labels = c("2", "4", "6", "8", "10", "12", "14", "15"))+
    labs(title = paste(batch, change, systemNum,":close contact, Phase: ", phase, nr), x = "ID", y = "ID")   # add labels and caption
  
  return(ggp)            
}


# input:  
#         
#        
#
# output: 
# effect:
generateHeatMapPositions <- function(count_position_list, batch, change, systemNum, phase, nr){
  #list into dataframe
  df_positions <- as.data.frame(do.call(rbind, count_position_list))
  #create old posizions tibble(translates position ids into coordinates)
  Positions_tibble <- tibble(PositionID = c(1:8), xPos = c(0,100,200,300,0,100,200,300), yPos = c(0,0,0,0,116,116,116,116))
  #merge both dataframes
  merged_df <- merge(df_positions, Positions_tibble, by.x = "V1", by.y = "PositionID", all = TRUE)
  # calculate second entrys to hour entrys
  hour_df <-  merged_df%>%
    mutate(V2= ifelse(V2!=0,V2/3600,V2))
  
  #rename the columns
  hour_df <- hour_df%>%
    rename(ID=V1)%>%
    rename(hours=V2)
  
  #print(hour_df)
  #print("max value: ")
  #print(max(hour_df$hours))
  
  
  # create heatmap with ggplot2
  heatmap <- ggplot(hour_df, aes(x = xPos, y = yPos, fill = hours)) +
    geom_tile() +
    scale_x_continuous(breaks = c(0, 100, 200, 300), labels = c("0", "100", "200", "300")) +
    scale_y_continuous(breaks = c(0, 116), labels = c("0", "116")) +
    scale_fill_gradientn(colors = c("yellow","orange", "red", "darkred", "#290000"), # colour palette
                         limits = c(0, 12), 
                         breaks = c(0, 2, 4, 6, 8, 10, ifelse(max(hour_df$hours)<12,12,warning("Higher scale required in heatmap for ", systemNum))),
                         labels = c("0", "2", "4", "6", "8", "10", "12")) + 
                         #limits are the borders of the scale,breaks are actual value breaks, labels are names for breakpoints
    labs(title = paste(batch, change, systemNum, ": used positions, Phase: ", phase, nr),
         x = "x-axis",
         y = "y-axis")

  return(heatmap)            
}



#graph

# input:  
#         
#        
#
# output: 
# effect:
generateGraph <- function(data, batch, change,animal_ids,system){
  
  #print(data)
  #print(colnames(data))
  
  #filter data to needed columns
  subset_data <- select(data, Phase, animal_ids[1], animal_ids[2], animal_ids[3], animal_ids[4])
  
  
  #print(subset_data)
  
  
  #data melt
  long_data <- melt(subset_data, id='Phase')
  
  #rename columns
  names(long_data) <- c('Phase', 'animal', 'time')
  
  #print(long_data)
  
  #print(typeof(long_data$time))
  #change seconds to hours
  hour_data <- long_data%>%
    mutate(time = as.integer(time))%>%
    mutate(time= ifelse(time!=0,time/3600,time))
  
  #print(typeof(hour_data$time))
  #print(hour_data)
  
  
  plot <- ggplot(data=hour_data, aes(x=Phase, y=time, color=animal))+
    geom_line(aes(group=animal))+
    geom_point()+
    scale_y_continuous("total close contact in h")+
    scale_x_discrete(limits = c("I1", "A1", "I2", "A2", "I3", "A3", "I4", "A4", "I5"))+ 
    scale_color_discrete(name = paste("animal from", batch, change, system))
    
 
    
    
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
create_joined_table <- function(CC1, CC2, CC3, CC4, animal_id, batch, stress_condition){
  
  # filter every tibble to the current animal id
  filtered_CC1 <-  CC1%>%
    select(Phase, animal_id)%>%
    rename(CC1=animal_id)
  print(filtered_CC1)
  
  filtered_CC2 <- CC2%>%
    select(Phase, animal_id)%>%
    rename(CC2=animal_id)
  
  filtered_CC3 <- CC3%>%
    select(Phase, animal_id)%>%
    rename(CC3=animal_id)
  
  filtered_CC4 <- CC4%>%
    select(Phase, animal_id)%>%
    rename(CC4=animal_id)
  
  # join the four filtered tibble to one 
  closeness_table_join <- Reduce(function (...) { merge(..., by = "Phase", all = TRUE) },   # Full join of reduced tibbles 1,2,3 and 4
                                 list(filtered_CC1, filtered_CC2, filtered_CC3, filtered_CC4))
  #sort the phase column by int
  closeness_table_join <- closeness_table_join%>%
    arrange(as.integer(sub("[^0-9]", "", Phase)))
  
  print(closeness_table_join)
  
  
  message("plotting")
  #data melt
  long_data <- melt(closeness_table_join, id='Phase')
  
  #rename columns
  names(long_data) <- c('Phase', 'change', 'time')
  
  #change seconds to hours
  hour_data <- long_data%>%
    mutate(time = as.integer(time))%>%
    mutate(time= ifelse(time!=0,time/3600,time))
  
  plot <- ggplot(data=hour_data, aes(x=Phase, y=time, color=change))+
    geom_line(aes(group=change),na.rm=TRUE)+
    geom_point(na.rm=TRUE)+
    #geom_path(linewidth = 10)+
    scale_color_manual(values = c("CC1" = "aquamarine2", "CC2" = "deepskyblue2", "CC3" = "deepskyblue4", "CC4" = "darkblue"),name = paste(animal, batch, stress_condition))+
    scale_y_continuous("total close contact in h")+
    scale_x_discrete(limits = c("I1", "A1", "I2", "A2", "I3", "A3", "I4", "A4", "I5"))+
    facet_grid(~change)
  return(plot)
}


############################################################################################
### statistic functions from tobi: ###
########################################################################################






# Function to perform normality test and appropriate statistical test for each variable and phase
testAndPlotVariable <- function(data, value, variableName, phase, sex) {
  #filtering specific phase or sex if needed 
  filteredData <- data %>%
    filter(if('Phase'%in% colnames(data))Phase == phase else TRUE) %>%   
    filter(Sex == sex)
  
  #factorize column group
  filteredData$Group <- as.factor(filteredData$Group)
  
  
  # save unique group names
  uniqueGroups <- unique(filteredData$Group)  #SUS,RES,CON...
  # number of different groups in data
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
generatePlot <- function(data, value, variableName, phase, sex) {
  filteredData <- data #%>%
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
