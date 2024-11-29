## 05/2024
## Anja Magister
## ANALYSIS OF SHANNON ENTROPY ##
##
##

# Load required packages
required_packages <- c("readr", "dplyr", "lubridate", "tibble", "purrr", "ggplot2", "reshape2", "scales", "stringr")

# Install required packages if not already installed
for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  } 
  library(package, character.only = TRUE)
}

# customisable variables
show_plots = TRUE
save_plots = TRUE
save_tables= TRUE

# paths

#setwd("C:/Users/topohl/Documents/GitHub/DLCAnalyzer")
#source('R/DLCAnalyzer_Functions_final.R')

working_directory <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability"

saving_directory <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability"
plots_directory <- "/plots"
#plots_directory <- "/R-plots"
tables_directory <- "/tables"

#functions
source(paste0("C:/Users/topohl/Documents/GitHub/MMMSociability/E9_SIS_AnimalPos-functions.R"))

#sus and con animals
sus_animals <- readLines(paste0(working_directory,"/raw_data/sus_animals.csv"))
##csv path for con animals
con_animals <- readLines(paste0(working_directory,"/raw_data/con_animals.csv"))


################################################################################################################################
## INITIALIZE RESULT TIBBLES ##
# These tibbles will store the results of the analysis.

# cagePosProb: Stores the probability of each position in the cage for each phase and system.
# cagePosEntropy: Stores the Shannon entropy of the cage for each phase and system.
# animalPosEntropy: Stores the Shannon entropy of each animal for each phase, system, and cage change.

# The tibbles are initialized with empty columns to be filled with the analysis results.

# Initialize the result tibbles for cage probability, cage entropy, and animal entropy
cagePosProb <- tibble(Batch = character(),
                      System = character(),
                      CageChange = character(),
                      Phase = character(),
                      Position = numeric(),
                      Probability = numeric())

cagePosEntropy <- tibble(Batch = character(),
                         Sex = character(),
                         System = character(),
                         CageChange = character(),
                         Phase = character(),
                         CageEntropy = numeric())

animalPosEntropy <- tibble(Batch = character(),
                          Sex = character(),
                          System = character(),
                          CageChange = character(),
                          Phase = character(),
                          AnimalID = character(),
                          animalEntropy = numeric())

################################################################################################################################

# define batch and cage change
batches <- c("B1", "B2", "B3", "B4", "B5", "B6")
cageChanges <- c("CC1", "CC2", "CC3", "CC4")

for(batch in batches){
  
  # Define sex based on batch
  sex <- ifelse(batch %in% c("B3", "B4", "B6"), "female", "male")
  
  for(cageChange in cageChanges){
    cat(batch, cageChange)
    
    ################################################################################################################################
    # Current CSV filename
    filename <- paste0("E9_SIS_", batch, "_", cageChange, "_AnimalPos")
    # Path of CSV file 
    csvFilePath <-  paste0(working_directory,"/preprocessed_data/", filename, "_preprocessed.csv")
    
    # Read preprocessed data (CSV file) into tibble
    preprocessed_data <- as_tibble(read_delim(csvFilePath, delim = ",", show_col_types = FALSE))
    
    ################################################################################################################################
    ## DEFINITIONS ##
    
    # Define systems
    uniqueSystems <- unique(preprocessed_data$System)
    uniqueSystems <- str_sort(uniqueSystems)
    # Define phases
    phases <- c("Active", "Inactive")
    
    # Determine the number of active and inactive phases in this record data
    active_phases_number <- unique(preprocessed_data$ConsecActive)
    inactive_phases_number <- unique(preprocessed_data$ConsecInactive)
    # Remove the 0
    active_phases_number <- active_phases_number[! active_phases_number %in% 0]
    inactive_phases_number <- inactive_phases_number[! inactive_phases_number %in% c(0,1,max(inactive_phases_number))] # Remove first and last phase as well because they are incomplete
    
    # Every animal in the data batch
    unique_animal <- unique(preprocessed_data$AnimalID)
    
    ################################################################################################################################
    ## ANALYSIS ##
    ##########################################################
    # For every system (different animal cages)
    # One system has 4 animal, one batch has 5 systems
    # There are 4 experimental systems and one control system
    
    ## FOR LOOP ## (goes through every of the 5 systems)
    for(system_id in uniqueSystems){
      
      #print(system_id)
      # Filter preprocessed_data to the actual system
      systemData <- preprocessed_data %>%
        filter(System == system_id) %>%
        as_tibble()
      
      # Define animal names of the system
      animal_ids <- unique(systemData$AnimalID)
      # Define boolean value for completeness of current system
      system_complete <- ifelse(length(animal_ids) < 4, FALSE, TRUE)
      # Fill vector with NAs if incomplete system (we need 4 values in vector)
      while(length(animal_ids) < 4){animal_ids <- append(animal_ids, NA)} # If there are less than 4 animal in a system (lost chip ...)
      
      ## FOR LOOP ## (difference between 2 phases)
      for(phase in phases){
        
        #print(phase) 
        
        # Special treatment for every CC4 after A2 (grid in cage) -> not usable for social analysis
        # Use only A1, I2, A2 for the computation of the rank from CC4
        if(cageChange == "CC4"){
          active_phases_number <- c(1, 2)
          inactive_phases_number <- c(2)
        }
        
        # Depending on the phase take the number of phases
        ifelse(phase == "Active",  number_of_phases <- active_phases_number,  number_of_phases <- inactive_phases_number)
        #print(number_of_phases)
        
        ## FOR LOOP ## (difference between the existing number of phases)
        for(nr in number_of_phases){
          print(paste0(batch, ", System: ", system_id, ", ", cageChange, ", ", phase, " phase #",nr))
          
          # Filter system data to the actual phase
          systemPhaseData <- systemData %>%
            filter(ConsecActive == ifelse(phase == "Active", nr, 0)) %>%   # Depending on phase the special number of the phase has to be selected
            filter(ConsecInactive == ifelse(phase == "Inactive", nr, 0)) %>%
            as_tibble()
          
          ## INITIALIZATIONS ##
          # Initialize animal lists with empty name, start time and start position of every animal in one system (4 animal together)
          animalOne    <- list(name = "", time = "", position = 0)
          animalTwo    <- list(name = "", time = "", position = 0)
          animalThree  <- list(name = "", time = "", position = 0)
          animalFour   <- list(name = "", time = "", position = 0)
          tempData    <- list(secTemp = 0, lineTemp = 0)
          
          # Combine them to a list of lists
          animal_list <- list(
            "animalOne" = animalOne,
            "animalTwo" = animalTwo,
            "animalThree" = animalThree,
            "animalFour" = animalFour,
            "tempData" = tempData)
          
          # Initialize probability lists
          # Cage probability
          # c(cagePosition, number of general observed seconds, number of added probabilities that a animal was on this position at this second, space for later calculated probability)
          # Last number of each vector has to be divided through the second number -> sum of probabilities of every second/number of seconds = average
          cage_prob_list <- list(c(1, 0, 0, 0), 
                                 c(2, 0, 0, 0), 
                                 c(3, 0, 0, 0), 
                                 c(4, 0, 0, 0), 
                                 c(5, 0, 0, 0), 
                                 c(6, 0, 0, 0), 
                                 c(7, 0, 0, 0), 
                                 c(8, 0, 0, 0))
          
          # animal probability
          animal_prob_tibble <- tibble(AnimalID = rep(animal_ids, each = 8),
                                     Position = rep(1:8, length.out = 32),
                                     Seconds = 0,
                                     SumPercentage = 0,
                                     Prob = 0)
          
          ## CALCULATIONS ##

          # Display a message indicating the start of animal ID processing
          message("Processing animal IDs:")

          # Update `animal_list` to initialize with the first position and start time of each animal
          animal_list <- find_first_pos_and_time(animal_ids, systemPhaseData, animal_list)

          # Initialize `timeTemp` with the start time from the first animal in the list
          # This is used to track the current time during the while loop
          timeTemp <- animal_list[[1]][[2]]

          # Initialize `lineTemp` with the first line index to process after the four initial lines
          # This is used to keep track of the current row being processed in `systemPhaseData`
          lineTemp <- 5

          # Initialize `secTemp` as 0 to store the time difference between consecutive rows during processing
          secTemp <- 0

          # Define the ending condition for the while loop based on the number of rows in `systemPhaseData`
          theEnd <- nrow(systemPhaseData) + 1

          ## WHILE LOOP ## 
          # Iterates through the rows of `systemPhaseData` to process animal positions and update probability data
          while (lineTemp != theEnd && lineTemp < theEnd) {

          # Create a copy of the current `animal_list` for comparison with the updated list
          # This is necessary for calculating changes between consecutive rows
          old_animal_list <- animal_list

          # Update `animal_list` with new positions and times from the current row
          animal_list <- update_animal_list(animal_ids, animal_list, systemPhaseData, timeTemp, lineTemp)

          # Update `secTemp` with the time difference (in seconds) between consecutive rows
          secTemp <- animal_list[["tempData"]][["secTemp"]]

          ##################################################################################################

          ## Use ANALYSIS FUNCTIONS to update result lists ##

          if (system_complete) {
            # Update cage probability list to reflect the likelihood of positions within the cage
            cage_prob_list <- check_cage_prob(old_animal_list, animal_list, cage_prob_list, secTemp)
          }

          # Update the probability table for individual animal positions
          animal_prob_tibble <- check_animal_prob(old_animal_list, animal_list, animal_prob_tibble, secTemp)

          ##################################################################################################

          # Update `lineTemp` with the next row index to process
          lineTemp <- animal_list[["tempData"]][["lineTemp"]]

          # Update `timeTemp` with the current time for the first animal
          timeTemp <- animal_list[[1]][[2]]

          } ## END WHILE LOOP ## (Processes rows of one `systemPhaseData`)

          ## RESULTS ##
          
          # Calculate the probability out of the given numbers (space 4 in the vectors of the list)
          ## Cage
          for(i in 1:8){
            cage_prob_list[[i]][[4]] <- cage_prob_list[[i]][[3]] / cage_prob_list[[i]][[2]]
          }
          ## animal
          animal_prob_tibble$Prob <- animal_prob_tibble$SumPercentage / animal_prob_tibble$Seconds
          
          # Print list results
          #message("Results: ")
          #message("animal_prob_tibble")
          # print(animal_prob_tibble)
          
          if(system_complete){
            message("enter data from this phase in total result tibble")
            # Enter information in big result tibble for every phase and system (later used for plotting)
            for(i in 1:8){ # For every position in cage
              p <- paste0(substr(phase, 1, 1), nr) # The current phase
              cagePosProb <- cagePosProb %>%
                add_row(Batch = batch,
                        System = system_id,
                        CageChange = cageChange, 
                        Phase = p, 
                        Position = i, 
                        Probability = cage_prob_list[[i]][[4]])
            }
          }
          
          ###### Calculate Shannon entropy: ######
          
          if(system_complete){
            ## Cage
            # Extract vector with the probabilities of every position
            cage_prob_vec <- cagePosProb %>%
              filter(Batch == batch) %>%
              filter(System == system_id) %>%
              filter(CageChange == cageChange) %>%
              filter(Phase == paste0(substr(phase, 1, 1), nr)) %>%
              pull(Probability)
            # Enter vector with the rest of information into calc function and calculate the Shannon entropy
            cage_shannon_entropy <- calc_shannon_entropy(cage_prob_vec)
            # Enter new row with Shannon entropy into result tibble
            cagePosEntropy <- cagePosEntropy %>%
              add_row(Batch = batch,
                      Sex = sex,
                      System = system_id,
                      CageChange = cageChange,
                      Phase = paste0(substr(phase, 1, 1), nr),
                      CageEntropy = cage_shannon_entropy)
          }  
          
          ## animal
          # Also probability vector
          for(animal in animal_ids){ # Should add 4 new rows, if system incomplete then less
            
            # If animal is not tracked (incomplete system)
            if(is.na(animal)) {
              next
            }
            
            animal_prob_vec <- animal_prob_tibble %>%
              filter(AnimalID == animal) %>%
              pull(Prob)
            # Calculate entropy
            animal_shannon_entropy <- calc_shannon_entropy(animal_prob_vec)
            # Add row
            animalPosEntropy <- animalPosEntropy %>%
              add_row(Batch = batch,
                      Sex = sex,
                      System = system_id,
                      CageChange = cageChange,
                      Phase = paste0(substr(phase, 1, 1), nr),
                      AnimalID = animal,
                      animalEntropy = animal_shannon_entropy)
          }
        } ## END FOR LOOP ## (nr of phases)
      } ## END FOR LOOP ## (phases)
    } ## END FOR LOOP ## (systems)
  }
}

############### save tables ########################################################################
#tables as csv data
if(save_tables==TRUE){
  message("save analysis tables")
  
  #result cage entropy
  write.csv(cagePosEntropy, file = paste0(saving_directory, tables_directory,"/all_batches", "_all_cageChanges", "_cagePosEntropy.csv"), row.names = FALSE)
  #result animal entropy
  write.csv(animalPosEntropy, file = paste0(saving_directory, tables_directory,"/all_batches", "_all_cageChanges", "_animalPosEntropy.csv"), row.names = FALSE)
}

############### read saved tables into a tibble(if you skip the analysis part on top) ########################################################################
cagePosEntropy <- as_tibble(read_delim(paste0(saving_directory, tables_directory,"/all_batches", "_all_cageChanges", "_cagePosEntropy.csv"),delim = ",", show_col_types = FALSE))
animalPosEntropy <- as_tibble(read_delim(paste0(saving_directory, tables_directory,"/all_batches", "_all_cageChanges", "_animalPosEntropy.csv"),delim = ",", show_col_types = FALSE))

####################################### CONSEC COLUMNS ##########################################################################################
#PREPROCESS
#rename tibbles to make a copy
consecCagePosEntropy <- cagePosEntropy
consecanimalPosEntropy <- animalPosEntropy

#definitions for dynamic processing
tibbleNames <- c("consecCagePosEntropy ", "consecanimalPosEntropy ")
values <- c('CageEntropy', 'animalEntropy')

for (i in seq_along(tibbleNames)) {
  name <- tibbleNames[i]
  if (name == "consecCagePosEntropy ") {
    dataTibble <- consecCagePosEntropy
  } else {
    dataTibble <- consecanimalPosEntropy
  }
  
  #create consec colums:
  #split Phase column
  dataTibble[c('Phase', 'Consec')] <- str_split_fixed(dataTibble$Phase, '', 2)
  
  # Change the consec column to ConsecActive and ConsecInactive
  # and change the content of the Phase column to "active" and "inactive"
  dataTibble <- dataTibble %>%
    # filter(Batch != "B6") %>% # filter out batch B6
    mutate(ConsecActive = ifelse(Phase == "A", as.numeric(Consec), 0)) %>%
    mutate(ConsecInactive = ifelse(Phase == "I", as.numeric(Consec) - 1, 0)) %>%
    mutate(Phase = ifelse(Phase == "A", "active", "inactive"))
  
  if (name == "consecanimalPosEntropy ") {
    dataTibble <- dataTibble %>% mutate(Group = ifelse(AnimalID %in% sus_animals, "sus", ifelse(AnimalID %in% con_animals, "con", "res")))
  }
  
  for (cc in c("CC1", "CC2", "CC3", "CC4")) {
    
    if (cc != "CC1") {
      dataTibble <- dataTibble %>%
        mutate(ConsecActive = ifelse(CageChange == cc & Phase == "active", ConsecActive + max_consecAct, ConsecActive)) %>%
        mutate(ConsecInactive = ifelse(CageChange == cc & Phase == "inactive", ConsecInactive + max_consecInact, ConsecInactive))
    }
    
    max_consecAct <- dataTibble %>%
      filter(CageChange == cc) %>%
      pull(ConsecActive) %>%
      unique() %>%
      max()
    #print(max_consecAct)
    
    max_consecInact <- dataTibble %>%
      filter(CageChange == cc) %>%
      pull(ConsecInactive) %>%
      unique() %>%
      max()
    #print(max_consecInact)
  }
  
  #bring columns into right order
  if (name == "consecanimalPosEntropy ") {
    dataTibble <- dataTibble[c('CageChange', 'Batch', 'System', 'AnimalID', 'Sex', 'Group', 'Phase', 'ConsecActive', 'ConsecInactive', 'animalEntropy')]
  }
  if (name == "consecCagePosEntropy ") {
    dataTibble <- dataTibble[c('CageChange', 'Batch', 'System', 'Sex', 'Phase', 'ConsecActive', 'ConsecInactive', 'CageEntropy')]
  }
  
  #overwrite old tibbles with new processed tibbles
  if (name == "consecCagePosEntropy ") {
    consecCagePosEntropy <- dataTibble
  } else {
    consecanimalPosEntropy <- dataTibble
  }
}

############### save consec tables ########################################################################
#tables as csv data
if(save_tables==TRUE){
  message("save consec tables")
  
  #result cage entropy
  write.csv(consecCagePosEntropy , file = paste0(saving_directory, tables_directory,"/all_batches", "_all_cageChanges", "_consecCagePosEntropy .csv"), row.names = FALSE)
  #result animal entropy
  write.csv(consecanimalPosEntropy , file = paste0(saving_directory, tables_directory,"/all_batches", "_all_cageChanges", "_consecanimalPosEntropy .csv"), row.names = FALSE)
}

############### show plots in R with new table structure ########################################################################
## outdated code, maybe use for statistics file!

### animal ENTROPY ###
#PREPROCESS
#create consec colums:
#split Phase column
cagePosEntropy[c('Phase', 'Consec')] <- str_split_fixed(cagePosEntropy$Phase, '', 2)
animalPosEntropy[c('Phase', 'Consec')] <- str_split_fixed(animalPosEntropy$Phase, '', 2)

# Change consec column to ConsecAct and ConsecInact
# Change content of Phase column to "active" and "inactive"
cagePosEntropy <- cagePosEntropy%>%
  # filter(Batch!= "B6")%>% # filter out batch B6
  mutate(ConsecActive = ifelse(Phase=="A", Consec, 0))%>%
  mutate(ConsecInactive = ifelse(Phase=="I", Consec, 0))%>%
  mutate(Phase = ifelse(Phase=="A", "active", "inactive"))

animalPosEntropy <- animalPosEntropy%>%
  # filter(Batch!= "B6")%>% # filter out batch B6
  mutate(ConsecActive = ifelse(Phase=="A", Consec, 0))%>%
  mutate(ConsecInactive = ifelse(Phase=="I", Consec, 0))%>%
  mutate(Phase = ifelse(Phase=="A", "active", "inactive"))%>%
  mutate(Group = ifelse(AnimalID %in% sus_animals, "sus", ifelse(AnimalID %in% con_animals, "con", "res"))) #group column(add maybe already before in analysis part( for cage entropy))

#GROUP AND PLOT IN LOOP
columns_to_group <- list(c("Sex","AnimalID","Group"), c("Sex","AnimalID","CageChange","Group"), c("Sex","AnimalID","Phase","Group"), c("Sex","AnimalID","Phase","Group"))#phase group 2* bc of 2 phases
x_axis <- list("Group","CageChange", "Group", "Group")
#counter to distinguish the two phases
phasecount <- 1

###
# Initialize lists to store the results and plots
results <- list()
plots <- list()

# Loop through the columns to group by
for (i in seq_along(columns_to_group)) {
  
  # Get the group variables and the x-axis variable
  group_vars <- columns_to_group[[i]]
  x_var <- x_axis[[i]]
  # Filter the data based on the group variables
  data <- animalPosEntropy
  
  # When looking at phases, we want to divide the tibble
  
  # Filter the data based on the group variables
  if("Phase" %in% group_vars) data <- filter(data, Phase == ifelse(phasecount==1,"active", "inactive"))
  # Count the number of phases
  if("Phase" %in% group_vars) phasecount <- phasecount+1
  
  #group tibble and calculate mean
  result <- data %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(Mean_Entropy = mean(animalEntropy, na.rm = TRUE), .groups = 'drop')
  
  print(paste("Grouping by:", paste(group_vars, collapse = ", ")))
  print(paste("x_axis:", x_var))
  
  # dynamic scatterplot
  p <- ggplot(data = result, aes(x = !!sym(x_var), y = Mean_Entropy, color = Group, group = Group)) + 
    geom_jitter(aes(fill=Group), size=4, alpha=0.7, shape=16, position=position_dodge(width = ifelse("CageChange" %in% group_vars, 0.75, 0))) +
    labs(title = paste("animal-Entropy-Plot\n for Grouping by", paste(group_vars, collapse = ", ")),
         x = x_var,
         y = "Mean Entropy") +
    scale_color_manual(values = c("sus" = "tomato", "res" = "darkgreen", "con" = "deepskyblue4")) +
    facet_grid(ifelse("Phase" %in% group_vars, "Phase~Sex", "~Sex"))+
    stat_summary(
      fun.min=function(z){quantile(z, 0.25)},
      fun.max=function(z){quantile(z, 0.75)},
      fun=median,
      color="black",
      size=0.8,
      shape=16, 
      position=position_dodge(width = ifelse("CageChange" %in% group_vars, 0.75, 0)))+
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          title = element_text(size = 20),
          legend.key.size = unit(3, "lines"),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          axis.text = element_text(size = 20), 
          axis.title = element_text(size = 22), 
          strip.text = element_text(size = 20))
  
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
animalPosEntropy_act <- animalPosEntropy%>%
  filter(Phase == "active")

animalPosEntropy_inact <- animalPosEntropy%>%
  filter(Phase == "inactive")

p <- ggplot(data = animalPosEntropy_inact, aes(x = Consec, y = animalEntropy, color = Group, group = Group)) + 
  geom_jitter(aes(fill=Group), size=4, alpha=0.7, shape=16, position=position_dodge(width = 0.75)) +
  labs(title = paste("animal-Entropy-Plot\n Inactive Phases")) +
  scale_color_manual(values = c("sus" = "tomato", "res" = "darkgreen", "con" = "deepskyblue4")) +
  facet_grid(Sex~CageChange)+
  stat_summary(
    fun.min=function(z){quantile(z, 0.25)},
    fun.max=function(z){quantile(z, 0.75)},
    fun=median,
    color="black",
    size=0.8,
    shape=16, 
    position=position_dodge(width = 0.75))+
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
if(show_plots==TRUE){
  message("show Plots")
}

### CAGE ENTROPY ###
cage_ent_plot <- ggplot(data=cagePosEntropy, aes(x=Phase, y=CageEntropy, color=System))+
  geom_point()+
  geom_line(aes(group = System)) +  # connects the dots between the phases with a line
  scale_x_discrete(limits = c("I1", "A1", "I2", "A2", "I3", "A3", "I4", "A4", "I5"))+
  facet_grid(~CageChange)

#divide into active/inactive
cagePosEntropy_act <- cagePosEntropy%>%
  filter(grepl("^A\\d+", Phase))

cagePosEntropy_inact<- cagePosEntropy%>%
  filter(grepl("^I\\d+", Phase))

cage_ent_plot_act <- ggplot(data=cagePosEntropy_act, aes(x=Phase, y=CageEntropy, color=System))+
  geom_jitter(aes(fill=System), size=4, alpha=0.7, width=0.2, shape=16)+
  scale_y_continuous("Cage Entropy all batches")+
  scale_x_discrete("active Phases")+
  labs(title = paste("Cage Entropy active Phases"))+
  stat_summary(
    fun.min=function(z){quantile(z, 0.25)},
    fun.max=function(z){quantile(z, 0.75)},
    fun=median,
    color="black",
    size=0.8,
    shape=16)+
  #theme_bw()+
  facet_grid(Sex~CageChange)+
  theme(title = element_text(size = 20),
        legend.key.size = unit(3, "lines"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 20), 
        axis.title = element_text(size = 22), 
        strip.text = element_text(size = 20))

cage_ent_plot_inact <- ggplot(data=cagePosEntropy_inact, aes(x=Phase, y=CageEntropy, color=System))+
  geom_jitter(aes(fill=System), size=4, alpha=0.7, width=0.2, shape=16)+
  scale_y_continuous("Cage Entropy all batches")+
  scale_x_discrete("inactive Phases")+
  labs(title = paste("Cage Entropy inactive Phases"))+
  stat_summary(
    fun.min=function(z){quantile(z, 0.25)},
    fun.max=function(z){quantile(z, 0.75)},
    fun=median,
    color="black",
    size=0.8,
    shape=16)+
  #theme_bw()+
  facet_grid(Sex~CageChange)+
  theme(title = element_text(size = 20),
        legend.key.size = unit(3, "lines"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 20), 
        axis.title = element_text(size = 22), 
        strip.text = element_text(size = 20))