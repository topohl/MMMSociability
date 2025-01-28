#' @title Analysis of Shannon Entropy
#' @description This script performs an analysis of Shannon entropy for animal positions.
#' @details 
#' This script is part of the MMMSociability project and is used to analyze the Shannon entropy 
#' of animal positions. The analysis is based on data collected and processed in previous steps 
#' of the project.
#' 
#' @file /C:/Users/topohl/Documents/GitHub/MMMSociability/E9_SIS_AnimalPos-analyzing-shannon_AnjasVersion.r
#' @date 01/2025
#' @author Tobias Pohl, Anja Magister
#' 
#' @section Notes:
#' - Ensure that all required libraries and data files are available before running the script.
#' - The script is designed to be run in an R environment.
#' - The results of the analysis will be saved in the specified output directory.
#' 
#' @section References:
#' - Please refer to the project documentation for more details on the methodology and data collection.
#' 
#' @keywords Shannon entropy, animal positions, MMMSociability

# Load required packages using pacman
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(readr, dplyr, lubridate, tibble, purrr, ggplot2, reshape2, scales, stringr)

# Customizable variables
show_plots <- FALSE
save_plots <- FALSE
save_tables <- TRUE

# Paths
working_directory <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability"
saving_directory <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability"
plots_directory <- "/plots"
tables_directory <- "/tables"

# Load custom functions
source(paste0(working_directory, "/E9_SIS_AnimalPos-functions_AnjasVersion.R"))

# Load lists of susceptible and control animals
sus_animals <- readLines(paste0(working_directory, "/raw_data/sus_animals.csv"))
con_animals <- readLines(paste0(working_directory, "/raw_data/con_animals.csv"))

#' Initialize Result Tibbles for Analysis
#'
#' This script initializes three empty tibbles that will store the results of the analysis:
#' - `cagePosProb`: Stores the probability of each position in the cage for each phase and system.
#' - `cagePosEntropy`: Stores the Shannon entropy of the cage for each phase and system.
#' - `animalPosEntropy`: Stores the Shannon entropy of each animal for each phase, system, and cage change.
#' 
#' The tibbles are initialized with the appropriate column names and types, but are empty initially.
#'
#' @details
#' - **`cagePosProb`** contains columns:
#'   - `Batch`: Batch identifier (character).
#'   - `System`: System identifier (character).
#'   - `CageChange`: Indicator of cage change (character).
#'   - `Phase`: Phase of the experiment (character).
#'   - `Position`: Position in the cage (numeric).
#'   - `Probability`: Probability associated with each position (numeric).
#'
#' - **`cagePosEntropy`** contains columns:
#'   - `Batch`: Batch identifier (character).
#'   - `Sex`: Sex of the animals (character).
#'   - `System`: System identifier (character).
#'   - `CageChange`: Indicator of cage change (character).
#'   - `Phase`: Phase of the experiment (character).
#'   - `CageEntropy`: Shannon entropy of the cage (numeric).
#'
#' - **`animalPosEntropy`** contains columns:
#'   - `Batch`: Batch identifier (character).
#'   - `Sex`: Sex of the animals (character).
#'   - `System`: System identifier (character).
#'   - `CageChange`: Indicator of cage change (character).
#'   - `Phase`: Phase of the experiment (character).
#'   - `AnimalID`: Unique identifier for each animal (character).
#'   - `animalEntropy`: Shannon entropy of the animal's position data (numeric).
#'
#' @return Three empty tibbles: `cagePosProb`, `cagePosEntropy`, and `animalPosEntropy`.
#' @examples
#' # Example usage
#' print(cagePosProb)  # Check the initialized tibble
#' print(cagePosEntropy)  # Check the initialized tibble
#' print(animalPosEntropy)  # Check the initialized tibble

# Initialize the result tibble for cage position probabilities
cagePosProb <- tibble(
  Batch = character(),       # Batch identifier
  System = character(),      # System identifier
  CageChange = character(),  # Cage change indicator
  Phase = character(),       # Experimental phase
  Position = numeric(),      # Cage position
  Probability = numeric()    # Probability of the position
)

# Initialize the result tibble for cage entropy
cagePosEntropy <- tibble(
  Batch = character(),       # Batch identifier
  Sex = character(),         # Sex of animals
  System = character(),      # System identifier
  CageChange = character(),  # Cage change indicator
  Phase = character(),       # Experimental phase
  CageEntropy = numeric()    # Shannon entropy of the cage
)

# Initialize the result tibble for animal entropy
animalPosEntropy <- tibble(
  Batch = character(),         # Batch identifier
  Sex = character(),           # Sex of animals
  System = character(),        # System identifier
  CageChange = character(),    # Cage change indicator
  Phase = character(),         # Experimental phase
  AnimalID = character(),      # Unique animal ID
  animalEntropy = numeric()    # Shannon entropy for each animal
)

# Define batches and cage changes
batches <- c("B1", "B2", "B3", "B4", "B5", "B6")
cageChanges <- c("CC1", "CC2", "CC3", "CC4")

for(batch in batches) {
  
  # Define sex based on batch
  sex <- ifelse(batch %in% c("B3", "B4", "B6"), "female", "male")
  
  # Iterate over cage changes
  for(cageChange in cageChanges) {
    cat(batch, cageChange)
    
    #' @title Read Preprocessed Data
    #' @description This function constructs the filename for the current CSV file based on the batch and cageChange variables, defines the full path to the CSV file, and reads the preprocessed data from the CSV file into a tibble.
    #' @param batch A variable representing the batch identifier.
    #' @param cageChange A variable representing the cage change identifier.
    #' @param working_directory A string representing the working directory where the preprocessed data is stored.
    #' @return A tibble containing the preprocessed data read from the CSV file.
    # Construct the filename for the current CSV file
    filename <- paste0("E9_SIS_", batch, "_", cageChange, "_AnimalPos")
    # Define the full path to the CSV file
    csvFilePath <- paste0(working_directory,"/preprocessed_data/", filename, "_preprocessed.csv")
    
    # Read the preprocessed data from the CSV file into a tibble
    data_preprocessed <- read_delim(csvFilePath,delim = ",", show_col_types = FALSE) %>%
      as_tibble()
    
    #' @title Define Systems and Phases
    #' @description This section defines the unique systems and phases present in the preprocessed data.
    #' It also identifies the number of active and inactive phases for the current batch.
    #' @param data_preprocessed A tibble containing the preprocessed data for the current batch and cage change.
    #' @return A list containing unique systems, phases, active phase count, inactive phase count, and unique animals.
    
    # Define unique systems
    unique_systems <- unique(data_preprocessed$System)
    unique_systems <- str_sort(unique_systems)
    
    # Define phases
    phases <- c("Active", "Inactive")
    
    # Identify the number of active and inactive phases in the preprocessed data
    active_phase_count <- unique(data_preprocessed$ConsecActive)
    inactive_phase_count <- unique(data_preprocessed$ConsecInactive)
    
    # Remove zero from phase counts
    active_phase_count <- active_phase_count[!active_phase_count %in% 0]
    inactive_phase_count <- inactive_phase_count[!inactive_phase_count %in% c(0, 1, max(inactive_phase_count))] # Remove first and last phase as they are incomplete
    
    # Define unique animals in the data batch
    unique_animals <- unique(data_preprocessed$AnimalID)
   
    #================================================
    # ANALYSIS
    #================================================
    # for every system (different animal cages)
    # one system has 4 animals, one batch has 5 systems
    # there are 4 experimental systems and one control system
    
    ## FOR LOOP ## (goes through every of the 5 systems)
    for(system_id in unique_systems) {
            
      print(system_id)
      
      # Filter the preprocessed data for the current system
      system_data <- data_preprocessed %>%
        filter(System == system_id) %>%
        as_tibble()
      
      # Define animal IDs for the current system
      animal_ids <- unique(system_data$AnimalID)
      
      # Determine if the current system is complete (i.e., has 4 animals)
      system_complete <- ifelse(length(animal_ids) < 4, FALSE, TRUE)
      
      # Ensure the animal_ids vector has 4 values, filling with NA if the system is incomplete
      while (length(animal_ids) < 4) {
        animal_ids <- append(animal_ids, NA)
      }
      
      # Loop through the phases (Active and Inactive)
      for(phase in phases) {
        
        print(phase) 
        
        # Special handling for CC4 after A2 (grid in cage) - not usable for social analysis
        # Use only A1, I2, A2 for the computation of the rank from CC4
        if(cageChange == "CC4") {
          active_phase_count <- c(1, 2)
          inactive_phase_count <- c(2)
        }
        
        # Determine the number of phases based on the current phase type
        ifelse(phase == "Active", phase_numbers <- active_phase_count, phase_numbers <- inactive_phase_count)
        print(phase_numbers)
        
        # Loop through the phase numbers
        for(phase_number in phase_numbers) {
          # Print the current batch, system, cage change, phase, and phase number for debugging
          print(paste0(batch, ", System: ", system_id, ", ", cageChange, ", ", phase, " phase number: ", phase_number))
          
          # Filter the data for the current phase and phase number
          data_system_phase <- system_data %>%
          filter(ConsecActive == ifelse(phase == "Active", phase_number, 0)) %>%   
          filter(ConsecInactive == ifelse(phase == "Inactive", phase_number, 0)) %>% 
            as_tibble()
          
          #================================================
          # INITIALIZATIONS
          #================================================

          # Initialize animal lists with empty name, start time, and start position for each animal in one system (4 animals together)
          animal_1 <- list(name = "", time = "", position = 0)
          animal_2 <- list(name = "", time = "", position = 0)
          animal_3 <- list(name = "", time = "", position = 0)
          animal_4 <- list(name = "", time = "", position = 0)
          data_temp <- list(elapsed_seconds = 0, current_row = 0)
            
          # Combine individual animal lists into a list of lists
          animal_list <- list(
            "animal_1" = animal_1,
            "animal_2" = animal_2,
            "animal_3" = animal_3,
            "animal_4" = animal_4,
            "data_temp" = data_temp
          )
            
          #' @title Initialize Probability Lists for Cage and Animal Positions
          #' 
          #' @description
          #' This script initializes two probability lists: one for cage positions and one for animal positions.
          #' 
          #' @details
          #' - `cage_position_probability`: A list of vectors where each vector represents a cage position.
          #'   Each vector contains:
          #'   - Cage position (1 to 8)
          #'   - Total observed seconds (initially 0)
          #'   - Sum of probabilities that an animal was at this position at each second (initially 0)
          #'   - Space for the calculated probability (initially 0)
          #'   The last element of each vector will be divided by the second element to get the average probability.
          #' 
          #' - `animal_position_probability`: A tibble containing the following columns:
          #'   - `AnimalID`: Repeated animal IDs for each of the 8 positions
          #'   - `Position`: Cage positions (1 to 8)
          #'   - `Seconds`: Total observed seconds (initially 0)
          #'   - `SumPercentage`: Sum of probabilities that an animal was at this position at each second (initially 0)
          #'   - `Prob`: Calculated probability (initially 0)
          
          # Cage position probability
          cage_position_probability <- list(
            c(1, 0, 0, 0), 
            c(2, 0, 0, 0), 
            c(3, 0, 0, 0), 
            c(4, 0, 0, 0), 
            c(5, 0, 0, 0), 
            c(6, 0, 0, 0), 
            c(7, 0, 0, 0), 
            c(8, 0, 0, 0)
          )
            
          # Animal position probability
          animal_position_probability <- tibble(
            AnimalID = rep(animal_ids, each = 8),
            Position = rep(1:8, length.out = 32),
            Seconds = 0,
            SumPercentage = 0,
            Prob = 0
          )
          
          #================================================
          # ANALYSIS
          #================================================

          message("Analysing data...")
          
          # Update animal_list to first time and first position
          animal_list <- initialize_animal_positions(animal_ids, data_system_phase, animal_list)
          
          initial_time <- animal_list[[1]][[2]]
          current_row <- 5
          elapsed_seconds <- 0
          total_rows <- nrow(data_system_phase) + 1
          
          # WHILE LOOP: Iterate through the data rows of one data_system_phase
          while (current_row != total_rows && current_row < total_rows) {
          
            # Create a copy of the old version of the animals list for comparison to new list (to calculate differences between the gaps of two rows)
            previous_animal_positions <- animal_list
          
            # Update animals list with new positions from the next row
            animal_list <- update_animal_list(animal_ids, animal_list, data_system_phase, initial_time, current_row)
          
            # Update elapsed_seconds
            elapsed_seconds <- animal_list[["data_temp"]][["elapsed_seconds"]]
          
            ##################################################################################################
          
            # Use ANALYSIS FUNCTIONS to update result lists
          
            if (system_complete) {
              # Update the probability of positions in the cage
              cage_position_probability <- update_cage_position_probability(previous_animal_positions, animal_list, cage_position_probability, elapsed_seconds)
            }
            
            # Update the probability of positions in the cage for each animal
            animal_position_probability <- update_animal_position_probability(previous_animal_positions, animal_list, animal_position_probability, elapsed_seconds)
          
            ##################################################################################################
          
            # Update current_row
            current_row <- animal_list[["data_temp"]][["current_row"]]
            # Update initial_time
            initial_time <- animal_list[[1]][[2]]
          }
          
          #================================================
          # RESULTS
          #================================================
          
          # Calculate the probability for each cage position
          for(i in 1:8) {
            cage_position_probability[[i]][[4]] <- cage_position_probability[[i]][[3]] / cage_position_probability[[i]][[2]]
          }
          
          # Calculate the probability for each animal position
          animal_position_probability$Prob <- animal_position_probability$SumPercentage / animal_position_probability$Seconds
          
          # Print the results
          message("Results: ")
          message("Animal Position Probability")
          print(animal_position_probability)
          
          if(system_complete) {
          message("Entering data from this phase into the total result tibble")
            # Enter information into the result tibble for each phase and system
            for(i in 1:8) { # For each position in the cage
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
            
          #================================================
          # SHANNON ENTROPY
          #================================================
            
          if(system_complete) {
            # Cage
            # Extract vector with the probabilities of every position
            cage_prob_vec <- cagePosProb %>%
              filter(Batch == batch) %>%
              filter(System == system_id) %>%
              filter(CageChange == cageChange) %>%
              filter(Phase == paste0(substr(phase, 1, 1), phase_number)) %>%
              pull(Probability)
            
            # Calculate the Shannon entropy for the cage
            cage_shannon_entropy <- calc_shannon_entropy(cage_prob_vec)
            
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
          for(animal in animal_ids) { # Should add 4 new rows, if system incomplete then less
              
            # Skip if animal is not tracked (incomplete system)
            if(is.na(animal)) {next}
              
            # Extract probability vector for the current animal
            animal_prob_vec <- animal_position_probability %>%
              filter(AnimalID == animal) %>%
              pull(Prob)
              
            # Calculate Shannon entropy for the current animal
            animal_shannon_entropy <- calc_shannon_entropy(animal_prob_vec)
              
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
  }
}

############### save tables ########################################################################

#tables as csv data
if(save_tables == TRUE) {
  message("save analysis tables")
  
  #result cage entropy
  write.csv(cagePosEntropy, file = paste0(saving_directory, tables_directory,"/all_batches", "_all_cageChanges", "_cagePosEntropy.csv"), row.names = FALSE)
  #result animals entropy
  write.csv(animalPosEntropy, file = paste0(saving_directory, tables_directory,"/all_batches", "_all_cageChanges", "_animalPosEntropy.csv"), row.names = FALSE)
}

############### read saved tables into a tibble(if you skip the analysis part on top) ########################################################################

cagePosEntropy <- as_tibble(read_delim(paste0(saving_directory, tables_directory,"/all_batches", "_all_cageChanges", "_cagePosEntropy.csv"),delim = ",", show_col_types = FALSE))

animalPosEntropy <- as_tibble(read_delim(paste0(saving_directory, tables_directory,"/all_batches", "_all_cageChanges", "_animalPosEntropy.csv"),delim = ",", show_col_types = FALSE))

####################################### CONSEC COLUMNS ##########################################################################################
#PREPROCESS
#rename tibbles to make a copy
consec_cage_entropy <- cagePosEntropy
consec_animal_entropy <- animalPosEntropy

#definitions for dynamic processing
entropy_tibble_names <- c("consec_cage_entropy", "consec_animal_entropy")
values <- c('CageEntropy', 'animalEntropy')

for (i in seq_along(entropy_tibble_names)) {
  name <- entropy_tibble_names[i]
  ifelse(name=="consec_cage_entropy", entropy_list <- consec_cage_entropy, entropy_list <- consec_animal_entropy)
  
  #create consec colums:
  #split Phase column
  entropy_list[c('Phase', 'Consec')] <- str_split_fixed(entropy_list$Phase, '', 2)
  
  #change consec column to ConsecAct and Inact
  #and change content of Phase column to active and inactive
  entropy_list <- entropy_list%>%
    filter(Batch!= "B6") %>%
    mutate(ConsecActive = ifelse(Phase=="A", as.numeric(Consec), 0)) %>%
    mutate(ConsecInactive = ifelse(Phase=="I", as.numeric(Consec)-1, 0)) %>%
    mutate(Phase = ifelse(Phase=="A", "active", "inactive"))
  
  if(name=="consec_animal_entropy")entropy_list <- entropy_list%>%mutate(Group = ifelse(AnimalID %in% sus_animals, "sus", ifelse(AnimalID %in% con_animals, "con", "res")))
  
  for(change in c("CC1","CC2", "CC3", "CC4")) {
    
    if(change!="CC1") {
      entropy_list <- entropy_list %>%
        mutate(ConsecActive = ifelse(CageChange == change & Phase == "active", ConsecActive + max_consecAct, ConsecActive)) %>%
        mutate(ConsecInactive = ifelse(CageChange == change & Phase == "inactive", ConsecInactive + max_consecInact, ConsecInactive))
      
    }
    
    max_consecAct <- entropy_list %>%
      filter(CageChange==change) %>%
      pull(ConsecActive) %>%
      unique() %>%
      max()
    #print(max_consecAct)
    
    max_consecInact <- entropy_list %>%
      filter(CageChange==change) %>%
      pull(ConsecInactive) %>%
      unique() %>%
      max()
    #print(max_consecInact)
  }
  
  #bring columns into right order
  if(name == "consec_animal_entropy")entropy_list <- entropy_list[c('CageChange', 'Batch', 'System', 'AnimalID', 'Sex', 'Group', 'Phase', 'ConsecActive', 'ConsecInactive', 'animalEntropy')]
  if(name == "consec_cage_entropy")entropy_list <- entropy_list[c('CageChange', 'Batch', 'System', 'Sex', 'Phase', 'ConsecActive', 'ConsecInactive', 'CageEntropy')]
  
  #overwrite old tibbles with new processed tibbles
  ifelse(name == "consec_cage_entropy", consec_cage_entropy <- entropy_list, consec_animal_entropy <- entropy_list)
}

############### save consec tables ########################################################################

#tables as csv data
if(save_tables==TRUE) {
  message("save consec tables")
  
  #result cage entropy
  write.csv(consec_cage_entropy, file = paste0(saving_directory, tables_directory,"/all_batches", "_all_cageChanges", "_consec_cage_entropy.csv"), row.names = FALSE)
  #result animals entropy
  write.csv(consec_animal_entropy, file = paste0(saving_directory, tables_directory,"/all_batches", "_all_cageChanges", "_consec_animal_entropy.csv"), row.names = FALSE)
}

############### show plots in R with new table structure ########################################################################

## outdated code, maybe use for statistics file!

### animals ENTROPY ###
#PREPROCESS
#create consec colums:
#split Phase column
cagePosEntropy[c('Phase', 'Consec')] <- str_split_fixed(cagePosEntropy$Phase, '', 2)
animalPosEntropy[c('Phase', 'Consec')] <- str_split_fixed(animalPosEntropy$Phase, '', 2)

#change consec column to ConsecAct and Inact
#and change content of Phase column to active and inactive
cagePosEntropy <- cagePosEntropy%>%
  filter(Batch!= "B6") %>%
  mutate(ConsecActive = ifelse(Phase=="A", Consec, 0)) %>%
  mutate(ConsecInactive = ifelse(Phase=="I", Consec, 0)) %>%
  mutate(Phase = ifelse(Phase=="A", "active", "inactive"))

animalPosEntropy <- animalPosEntropy%>%
  filter(Batch!= "B6") %>%
  mutate(ConsecActive = ifelse(Phase == "A", Consec, 0)) %>%
  mutate(ConsecInactive = ifelse(Phase == "I", Consec, 0)) %>%
  mutate(Phase = ifelse(Phase=="A", "active", "inactive")) %>%
  mutate(Group = ifelse(AnimalID %in% sus_animals, "sus", ifelse(AnimalID %in% con_animals, "con", "res"))) #group column(add maybe already before in analysis part( for cage entropy))

#GROUP AND PLOT IN LOOP

columns_to_group <- list(c("Sex","AnimalID","Group"), c("Sex","AnimalID","CageChange","Group"), c("Sex","AnimalID","Phase","Group"), c("Sex","AnimalID","Phase","Group"))#phase group 2* bc of 2 phases
x_axis <- list("Group","CageChange", "Group", "Group")
#counter to distinguish the two phases
phasecount <- 1

###
# Initialisieren einer Liste für die Ergebnisse
results <- list()
plots <- list()

# For-Loop
for (i in seq_along(columns_to_group)) {
  
  #variables from dynamic lists
  group_vars <- columns_to_group[[i]]
  x_var <- x_axis[[i]]
  #working copy of data(needed)
  data <- animalPosEntropy
  
  #dynamic filtering/grouping:
  
  #when looking at phases we want to divide the tibble
  if("Phase" %in% group_vars) data <- filter(data, Phase == ifelse(phasecount==1,"active", "inactive"))
  #increase phasecount for the second phase
  if("Phase" %in% group_vars) phasecount <- phasecount+1
  
  #group tibble and calculate mean
  result <- data %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(Mean_Entropy = mean(animalEntropy, na.rm = TRUE), .groups = 'drop')
  
  
  print(paste("Grouping by:", paste(group_vars, collapse = ", ")))
  print(paste("X axis:", x_var))
  
  
  # dynamic scatterplot
  p <- ggplot(data = result, aes(x = !!sym(x_var), y = Mean_Entropy, color = Group, group = Group)) + 
    geom_jitter(aes(fill=Group), 
      size=4, 
      alpha=0.7, 
      shape=16, 
      position = position_dodge(width = ifelse("CageChange" %in% group_vars, 0.75, 0))
    ) +
    labs(title = paste("Animal-Entropy-Plot\n for Grouping by", paste(group_vars, collapse = ", ")),
         x = x_var,
         y = "Mean Entropy") +
    scale_color_manual(values = c("sus" = "tomato", "res" = "darkgreen", "con" = "deepskyblue4")) +
    facet_grid(ifelse("Phase" %in% group_vars, "Phase~Sex", "~Sex")) +
    stat_summary(
      fun.min=function(z) {quantile(z, 0.25)},
      fun.max = function(z) {quantile(z, 0.75)},
      fun=median,
      color="black",
      size=0.8,
      shape=16, 
      position = position_dodge(width = ifelse("CageChange" %in% group_vars, 0.75, 0))
    ) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          title = element_text(size = 20),
          legend.key.size = unit(3, "lines"),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          axis.text = element_text(size = 20), 
          axis.title = element_text(size = 22), 
          strip.text = element_text(size = 20)
    )
  
  #print(p)
  
  plots[[i]] <- p
}

# Print all generated plots
for (i in seq_along(plots)) {
  print(plots[[i]])
}

##plots outside of loop for consec, not perfect consec tibble as base
##consec plots active/inactive
#filter phases
animalPosEntropy_act <- animalPosEntropy %>%
  filter(Phase == "active")

animalPosEntropy_inact <- animalPosEntropy %>%
  filter(Phase == "inactive")

p <- ggplot(data = animalPosEntropy_inact, aes(x = Consec, y = animalEntropy, color = Group, group = Group)
  ) + 
  geom_jitter(aes(fill = Group),
    size = 4, 
    alpha=0.7, 
    shape=16, 
    position = position_dodge(width = 0.75)
  ) +
  labs(title = paste("Animal-Entropy-Plot\n Inactive Phases")) +
  scale_color_manual(values = c("sus" = "tomato", "res" = "darkgreen", "con" = "deepskyblue4")) +
  facet_grid(Sex~CageChange) +
  stat_summary(
    fun.min=function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun=median,
    color="black",
    size=0.8,
    shape=16, 
    position=position_dodge(width = 0.75)) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        title = element_text(size = 20),
        legend.key.size = unit(3, "lines"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 20), 
        axis.title = element_text(size = 22), 
        strip.text = element_text(size = 20))


############### show plots in R ########################################################################
if(show_plots==TRUE) {
  message("show Plots")
}

### CAGE ENTROPY ###
cage_ent_plot <- ggplot(data=cagePosEntropy, aes(x = Phase, y = CageEntropy, color = System)) +
  geom_point() +
  geom_line(aes(group = System)) +  # connects the dots between the phases with a line
  scale_x_discrete(limits = c("I1", "A1", "I2", "A2", "I3", "A3", "I4", "A4", "I5")) +
  facet_grid(~CageChange)

#divide into active/inactive
cagePosEntropy_act <- cagePosEntropy%>%
  filter(grepl("^A\\d+", Phase))

cagePosEntropy_inact<- cagePosEntropy%>%
  filter(grepl("^I\\d+", Phase))

cage_ent_plot_act <- ggplot(data=cagePosEntropy_act, aes(x = Phase, y = CageEntropy, color=System)) +
  geom_jitter(aes(fill=System), size=4, alpha=0.7, width=0.2, shape=16) +
  scale_y_continuous("Cage Entropy all batches") +
  scale_x_discrete("active Phases") +
  labs(title = paste("Cage Entropy active Phases")) +
  stat_summary(
    fun.min = function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun=median,
    color="black",
    size=0.8,
    shape=16) +
  #theme_bw() +
  facet_grid(Sex~CageChange) +
  theme(title = element_text(size = 20),
        legend.key.size = unit(3, "lines"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 20), 
        axis.title = element_text(size = 22), 
        strip.text = element_text(size = 20))

cage_ent_plot_inact <- ggplot(data=cagePosEntropy_inact, aes(x = Phase, y = CageEntropy, color=System)) +
  geom_jitter(aes(fill=System), size=4, alpha=0.7, width=0.2, shape=16) +
  scale_y_continuous("Cage Entropy all batches") +
  scale_x_discrete("inactive Phases") +
  labs(title = paste("Cage Entropy inactive Phases")) +
  stat_summary(
    fun.min=function(z) {quantile(z, 0.25)},
    fun.max = function(z) {quantile(z, 0.75)},
    fun = median,
    color = "black",
    size = 0.8,
    shape = 16
  ) +
  #theme_bw() +
  facet_grid(Sex~CageChange) +
  theme(title = element_text(size = 20),
        legend.key.size = unit(3, "lines"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 20), 
        axis.title = element_text(size = 22), 
        strip.text = element_text(size = 20))