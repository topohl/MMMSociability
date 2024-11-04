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
source(paste0(working_directory,"/E9_SIS_AnimalPos-functions.R"))

#sus and con animals
sus_animals <- readLines(paste0(working_directory,"/raw_data/sus_animals.csv"))
##csv path for con animals
con_animals <- readLines(paste0(working_directory,"/raw_data/con_animals.csv"))


################################################################################################################################
## INITIALIZE RESULT TIBBLES ##
# These tibbles will store the results of the analysis.

# cagePosProb: Stores the probability of each position in the cage for each phase and system.
# cagePosEntropy: Stores the Shannon entropy of the cage for each phase and system.
# mousePosEntropy: Stores the Shannon entropy of each mouse for each phase, system, and cage change.

# The tibbles are initialized with empty columns to be filled with the analysis results.

# Initialize the result tibbles for cage probability, cage entropy, and mouse entropy
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

mousePosEntropy <- tibble(Batch = character(),
                Sex = character(),
                System = character(),
                CageChange = character(),
                Phase = character(),
                AnimalID = character(),
                MiceEntropy = numeric())

################################################################################################################################

# define batch and cage change
batches <- c("B1", "B2", "B3", "B4", "B5", "B6")
cageChanges <- c("CC1", "CC2", "CC3", "CC4")

for(batch in batches){
  
  #define sex
  sex <- ifelse(batch %in% c("B3", "B4", "B6"), "female", "male")
  
  for(cageChange in cageChanges){
  #for(cageChange in c("CC1")){
    cat(batch, cageChange)
    
    ################################################################################################################################
    #current csv filename
    filename <- paste0("E9_SIS_", batch, "_", cageChange, "_AnimalPos")
    #path of csv file 
    csvFilePath <-  paste0(working_directory,"/preprocessed_data/", filename, "_preprocessed.csv")
    
    # read preprocessed data(csv file) in tibble
    overallData <- as_tibble(read_delim(csvFilePath,delim = ",", show_col_types = FALSE))
    
    ################################################################################################################################
    ## DEFINITIONS ##
    
    #define systems
    uniqueSystems <- unique(overallData$System)
    uniqueSystems <- str_sort(uniqueSystems)
    #define phases
    phases <- c("Active", "Inactive")
    
    #look for number of active and inactive phases in this record data
    active_phases_number <- unique(overallData$ConsecActive)
    inactive_phases_number <- unique(overallData$ConsecInactive)
    #delete the 0
    active_phases_number <- active_phases_number[! active_phases_number %in% 0]
    inactive_phases_number <- inactive_phases_number[! inactive_phases_number %in% c(0,1,max(inactive_phases_number))]#delete first and last phase as well bc they are incomplete
    
    #every mouse in the data batch
    unique_mice <- unique(overallData$AnimalID)
    
    ################################################################################################################################
    ## ANALYSIS ##
    ##########################################################
    # for every system (different mouse cages)
    # one system has 4 mice, one batch has 5 systems
    # there are 4 experimental systems and one control system
    
    ## FOR LOOP ## (goes through every of the 5 systems)
    for(systemName in uniqueSystems){
      
      print(systemName)
      #sort overallData to the actual system
      systemData <- overallData%>%
        filter(System==systemName)%>%
        as_tibble()
      #print(system)
      
      #define mouse names of the system
      mouse_names <- unique(systemData$AnimalID)
      #define boolean value for completeness of current system
      system_complete <- ifelse(length(mouse_names)<4,FALSE,TRUE)
      #fill vector with NAs if incomplete system(we need 4 values in vector)
      while(length(mouse_names)<4){mouse_names <- append(mouse_names,NA)}#if there are less than 4 mice in a system(lost chip ...)
      
      ## FOR LOOP ## (difference between 2 phases)
      for(phase in phases){
        
        print(phase) 
        
        #special treatment for every CC4 after A2(grid in cage)-> not usable for social analysis
        #use only A1,I2,A2 for the computation of the rank from CC4
        if(cageChange=="CC4"){
          active_phases_number <- c(1,2)
          inactive_phases_number <- c(2)
        }
        
        #depending on the phase take the number of phases
        ifelse(phase == "Active",  number_of_phases <- active_phases_number,  number_of_phases <- inactive_phases_number)
        print(number_of_phases)
        
        ## FOR LOOP ## (difference between the existing number of phases)
        for(nr in number_of_phases){
          print(paste0( batch, ", System: ", systemName, ", ", cageChange, ", ", phase, " phase nr: ", nr))
          
          #sort system data to the actual phase
          systemPhaseData <- systemData%>%
            filter(ConsecActive == ifelse(phase == "Active",nr, 0))%>%   #depending on phase the special number of the phase has to be selected
            filter(ConsecInactive == ifelse(phase == "Inactive",nr, 0))%>%
            as_tibble()
          #print(systemPhaseData)
          
          
          ##INITIALIZATIONS##
          # initialize mice lists with empty name, start time and start position of every mouse in one system(4mice together)
          mouseOne    <- list(name="", time="", position=0)
          mouseTwo    <- list(name="", time="", position=0)
          mouseThree  <- list(name="", time="", position=0)
          mouseFour   <- list(name="", time="", position=0)
          tempData    <- list(secTemp=0, lineTemp=0)
          
          # combine them to a list of lists
          mice_list <- list(
            "mouseOne" = mouseOne,
            "mouseTwo" = mouseTwo,
            "mouseThree" = mouseThree,
            "mouseFour" = mouseFour,
            "tempData" = tempData)
          
          
          #initialize probability lists
          #cage probability
          #c(cagePosition, number of general observed seconds,number of added probabilitys that a mouse was on this position at this second, space for later calclated probability)
          #last number of each vector has to be divided through the second number-> sum of probabilitys of every second/number of seconds= average
          cage_prob_list <- list( c(1, 0, 0, 0), 
                                  c(2, 0, 0, 0), 
                                  c(3, 0, 0, 0), 
                                  c(4, 0, 0, 0), 
                                  c(5, 0, 0, 0), 
                                  c(6, 0, 0, 0), 
                                  c(7, 0, 0, 0), 
                                  c(8, 0, 0, 0))
          
          #mice probability
          mice_prob_tibble <- tibble(AnimalID=rep(mouse_names,each=8),
                                     Position=rep(1:8, length.out = 32),
                                     Seconds=0,
                                     SumPercentage=0,
                                     Prob=0)
          
          
          ##CALCULATIONS##
          message("calculates")
          
          
          #update mice_list to first time and first position
          mice_list <- find_first_pos_and_time(mouse_names, systemPhaseData, mice_list)
          
          #assign start time(choose one of the mices start time) for while loop
          #contains object of DateTime
          #for mice list update
          timeTemp <- mice_list[[1]][[2]]
          
          #assign first line number(after the four initial lines) for while loop
          #contains an int
          #needed for while loop ending 
          lineTemp <- 5
          
          #assign first seconds difference between two entrys in data for while loop
          #contains an int
          secTemp <- 0
          
          #define ending for while loop
          theEnd <- nrow(systemPhaseData)+1
          
          ## WHILE LOOP ## (goes through the data rows of one systemPhaseData)
          while(lineTemp!=theEnd && lineTemp<theEnd){
            
            #create a copy of the old version of the mice list for comparison to new list(to calculate differences between the gaps of two rows)
            old_mice_list <- mice_list
            
            #update mice list with new positions from the next row
            mice_list <- update_mice_list(mouse_names, mice_list, systemPhaseData, timeTemp, lineTemp)
            
            #update secTemp
            secTemp <-  mice_list[["tempData"]][["secTemp"]]
            
            ##################################################################################################
            
            ## use ANALYSIS FUNCTIONS to update result lists ##
            
            if(system_complete){
              ## probability of positions in cage 
              cage_prob_list <- check_cage_prob(old_mice_list,mice_list,cage_prob_list,secTemp)
            }
            
            ## probability of positions in cage of every mouse
            mice_prob_tibble <- check_mice_prob(old_mice_list,mice_list,mice_prob_tibble,secTemp)
            
            ##################################################################################################
            
            
            #update lineTemp
            lineTemp <- mice_list[["tempData"]][["lineTemp"]]
            #update timeTemp
            timeTemp <- mice_list[[1]][[2]]
          }## END WHILE LOOP ##(rows of one systemPhaseData)
          
          ## RESULTS ##
          
          #calculate the probability out of the given numbers (space 4 in the vectors of the list)
          ##cage
          for(i in 1:8){
            cage_prob_list[[i]][[4]] <- cage_prob_list[[i]][[3]]/cage_prob_list[[i]][[2]]
          }
          ##mice
          mice_prob_tibble$Prob <- mice_prob_tibble$SumPercentage/mice_prob_tibble$Seconds
          
          #print list results
          message("results: ")
          #message("cage_prob_list")
          #print(cage_prob_list)
          message("mice_prob_tibble")
          print(mice_prob_tibble)
          
          if(system_complete){
            message("enter data from this phase in total result tibble")
            #enter information in big result tibble for every phase and system(later used for plotting)
            for(i in 1:8){#for every position in cage
              p <- paste0(substr(phase, 1, 1),nr) #the current phase
              cagePosProb <- cagePosProb%>%
                add_row(Batch=batch,
                        System=systemName,
                        CageChange = cageChange, 
                        Phase = p, 
                        Position=i, 
                        Probability=cage_prob_list[[i]][[4]] )
              
            }
          }
            
            ###### calculate shannon entropy: ######
            
          if(system_complete){
            ##cage
            #extract vector with the probs of every position
            cage_prob_vec <- cagePosProb%>%
              filter(Batch==batch)%>%
              filter(System==systemName)%>%
              filter(CageChange==cageChange)%>%
              filter(Phase==paste0(substr(phase,1,1),nr))%>%
              pull(.,Probability)
            #enter vector with the rest of inform. into calc function and calcuate the shannon entropy
            cage_shannon_entropy <- calc_shannon_entropy(cage_prob_vec)
            #enter new row with shannon entropy into result tibble
            cagePosEntropy <- cagePosEntropy%>%
              add_row(Batch=batch,
                      Sex=sex,
                      System=systemName,
                      CageChange=cageChange,
                      Phase=paste0(substr(phase,1,1),nr),
                      CageEntropy=cage_shannon_entropy)
          }  
          
            
            ##mice
            #also prob vec
            for(mouse in mouse_names){#should add 4 new rows, if system incomplete then less
              
              #if mouse is not tracked(incomplete system)
              if(is.na(mouse)){next}
              
              mice_prob_vec <- mice_prob_tibble%>%
                filter(AnimalID==mouse)%>%
                pull(Prob)
              #calc entropy
              mice_shannon_entropy <- calc_shannon_entropy(mice_prob_vec)
              #add row
              mousePosEntropy<- mousePosEntropy%>%
                add_row(Batch=batch,
                        Sex=sex,
                        System=systemName,
                        CageChange=cageChange,
                        Phase=paste0(substr(phase,1,1),nr),
                        AnimalID=mouse,
                        MiceEntropy=mice_shannon_entropy)
            }
        }## END FOR LOOP ##(nr of phases)
      }## END FOR LOOP ##(phases)
    }
    ## END FOR LOOP ##(systems)
  }
}



############### save tables ########################################################################


#tables as csv data
if(save_tables==TRUE){
  message("save analysis tables")
  
  #result cage entropy
  write.csv(cagePosEntropy, file = paste0(saving_directory, tables_directory,"/all_batches", "_all_cageChanges", "_cagePosEntropy.csv"), row.names = FALSE)
  #result mice entropy
  write.csv(mousePosEntropy, file = paste0(saving_directory, tables_directory,"/all_batches", "_all_cageChanges", "_mousePosEntropy.csv"), row.names = FALSE)
}

############### read saved tables into a tibble(if you skip the analysis part on top) ########################################################################

cagePosEntropy <- as_tibble(read_delim(paste0(saving_directory, tables_directory,"/all_batches", "_all_cageChanges", "_cagePosEntropy.csv"),delim = ",", show_col_types = FALSE))

mousePosEntropy <- as_tibble(read_delim(paste0(saving_directory, tables_directory,"/all_batches", "_all_cageChanges", "_mousePosEntropy.csv"),delim = ",", show_col_types = FALSE))



####################################### CONSEC COLUMNS ##########################################################################################
#PREPROCESS
#rename tibbles to make a copy
consecCagePosEntropy <- cagePosEntropy
consecMousePosEntropy <- mousePosEntropy

#definitions for dynamic processing
tibbleNames <- c("consecCagePosEntropy ", "consecMousePosEntropy ")
values <- c('CageEntropy', 'MiceEntropy')

for (i in seq_along(tibbleNames)) {
  name <- tibbleNames[i]
  if (name == "consecCagePosEntropy ") {
    dataTibble <- consecCagePosEntropy
  } else {
    dataTibble <- consecMousePosEntropy
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
  
  if (name == "consecMousePosEntropy ") {
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
  if (name == "consecMousePosEntropy ") {
    dataTibble <- dataTibble[c('CageChange', 'Batch', 'System', 'AnimalID', 'Sex', 'Group', 'Phase', 'ConsecActive', 'ConsecInactive', 'MiceEntropy')]
  }
  if (name == "consecCagePosEntropy ") {
    dataTibble <- dataTibble[c('CageChange', 'Batch', 'System', 'Sex', 'Phase', 'ConsecActive', 'ConsecInactive', 'CageEntropy')]
  }
  
  #overwrite old tibbles with new processed tibbles
  if (name == "consecCagePosEntropy ") {
    consecCagePosEntropy <- dataTibble
  } else {
    consecMousePosEntropy <- dataTibble
  }
}

############### save consec tables ########################################################################

#tables as csv data
if(save_tables==TRUE){
  message("save consec tables")
  
  #result cage entropy
  write.csv(consecCagePosEntropy , file = paste0(saving_directory, tables_directory,"/all_batches", "_all_cageChanges", "_consecCagePosEntropy .csv"), row.names = FALSE)
  #result mice entropy
  write.csv(consecMousePosEntropy , file = paste0(saving_directory, tables_directory,"/all_batches", "_all_cageChanges", "_consecMousePosEntropy .csv"), row.names = FALSE)
}

############### show plots in R with new table structure ########################################################################

## outdated code, maybe use for statistics file!

### MICE ENTROPY ###
#PREPROCESS
#create consec colums:
#split Phase column
cagePosEntropy[c('Phase', 'Consec')] <- str_split_fixed(cagePosEntropy$Phase, '', 2)
mousePosEntropy[c('Phase', 'Consec')] <- str_split_fixed(mousePosEntropy$Phase, '', 2)

# Change consec column to ConsecAct and ConsecInact
# Change content of Phase column to "active" and "inactive"
cagePosEntropy <- cagePosEntropy%>%
  # filter(Batch!= "B6")%>% # filter out batch B6
  mutate(ConsecActive = ifelse(Phase=="A", Consec, 0))%>%
  mutate(ConsecInactive = ifelse(Phase=="I", Consec, 0))%>%
  mutate(Phase = ifelse(Phase=="A", "active", "inactive"))

mousePosEntropy <- mousePosEntropy%>%
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
  data <- mousePosEntropy
  
  # When looking at phases, we want to divide the tibble
  
  # Filter the data based on the group variables
  if("Phase" %in% group_vars) data <- filter(data, Phase == ifelse(phasecount==1,"active", "inactive"))
  # Count the number of phases
  if("Phase" %in% group_vars) phasecount <- phasecount+1
  
  #group tibble and calculate mean
  result <- data %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(Mean_Entropy = mean(MiceEntropy, na.rm = TRUE), .groups = 'drop')
  
  print(paste("Grouping by:", paste(group_vars, collapse = ", ")))
  print(paste("X axis:", x_var))
  
  # dynamic scatterplot
  p <- ggplot(data = result, aes(x = !!sym(x_var), y = Mean_Entropy, color = Group, group = Group)) + 
    geom_jitter(aes(fill=Group), size=4, alpha=0.7, shape=16, position=position_dodge(width = ifelse("CageChange" %in% group_vars, 0.75, 0))) +
    labs(title = paste("Mice-Entropy-Plot\n for Grouping by", paste(group_vars, collapse = ", ")),
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
mousePosEntropy_act <- mousePosEntropy%>%
  filter(Phase == "active")

mousePosEntropy_inact <- mousePosEntropy%>%
  filter(Phase == "inactive")

p <- ggplot(data = mousePosEntropy_inact, aes(x = Consec, y = MiceEntropy, color = Group, group = Group)) + 
  geom_jitter(aes(fill=Group), size=4, alpha=0.7, shape=16, position=position_dodge(width = 0.75)) +
  labs(title = paste("Mice-Entropy-Plot\n Inactive Phases")) +
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