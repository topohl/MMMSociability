## 11/2023
## Anja Magister
## ANALYSIS OF ANIMAL POSITIONS - ANALYZING ##
##
##
## NEEDED FILE STRUCTURE IN WORKING DIRECTORY
## - there must be a folder called "preprocessed_data" in which our preprocessed data is located from the previous code
## - there must be a folder called "plots", because the results are saved in this folder (if save_plots == TRUE)
##
## CUSTOMISABLE VARIABLES:
## - show_plots: do you want to see the plots in R?
## - save_plots_ do you want the generated plots to be saved?
## - working directory: choose your current directory
## - batch and cage change: on which file are you working->filename will be created with batch and cc information


# libraries
library(readr)        # load readr package for reading csv files
library(dplyr)
library(lubridate)    # for rounding time, time operations in general
library(tibble)       #important for tibble operations
library(purrr)
library(ggplot2)      #for plots
library(reshape2)     #for heatmap plot
library(scales)       # for heatmap rescale
library("stringr")    # for sorting vector of strings   

# customisable variables
show_plots = FALSE
save_plots = FALSE
save_tables= FALSE



# paths
working_directory <- "S:/Lab_Member/Anja/Git/MDC_Bachelor/E9_SIS_AnimalPos"
#working_directory <- "/home/anja/Dokumente/FU BERLIN/BA/Git/MDC_Bachelor/E9_SIS_AnimalPos"


#functions
source(paste0(working_directory,"/E9_SIS_AnimalPos-functions.R"))

# define batch and cage change
batches <- c("B1", "B2", "B3", "B4", "B5")#, "B6"
cageChanges <- c("CC1", "CC2", "CC3", "CC4")


#initialize result heatmap lists
#allHeatmaps <- list()
allHeatmaps_closeness <- list()
allHeatmaps_positions <- list()

# for total closeness plots
all_plots_total_closeness <- list()





for(batch in batches){
  
  for(cageChange in cageChanges){
    print(paste(batch, cageChange))
    
    #current csv filename
    filename <- paste0("E9_SIS_", batch, "_", cageChange, "_AnimalPos")
    #path of csv file 
    csvFilePath <-  paste0(working_directory,"/preprocessed_data/", filename, "_preprocessed.csv")
    
    # read preprocessed data(csv file) in tibble
    overallData <- as_tibble(read_delim(csvFilePath,delim = ",", show_col_types = FALSE))
    
    #convert dateTime from utc to utc+1 (Berlin)?
    
    ################################################################################################################################
    ## given data "overallData" has already been processed in other file ##
    ################################################################################################################################
    ## DEFINITIONS ##
    
    #define systems
    uniqueSystems <- unique(overallData$System)
    uniqueSystems <- str_sort(uniqueSystems)
    
    #define days
    #filter the unique day dates
    days <- unique(format(overallData$DateTime, "%D"))
    
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
    ## INITIALIZE RESULT LISTS ##

    
    #list of heatmaps every system
    systemHeatmaps_closeness <- list()
    systemHeatmaps_positions <- list()
    
  
    
    #create vector of phases
    #in the correct order of the phases how they appear in real
    phases_column <- c()
    for(i in 1:max(length(inactive_phases_number),length(active_phases_number))){
      if(i<=length(active_phases_number)){
        phases_column <- c(phases_column, paste0("A",i))
      }
      if(i<=length(inactive_phases_number)){
        phases_column <- c(phases_column, paste0("I",i+1))#we start with the second inactive phase
      }
    }
    
    # tibble for closeness:
    #create tibble with Phase column
    result_total_closeness <- tibble("Phase"=phases_column)
    #add columns for every existing mouse in this batch
    for(mouse in unique_mice){
      result_total_closeness[[mouse]] <- NA
    }
    
    # tibble for movement:
    #create tibble with Phase column
    result_total_movement <- tibble("Phase"=phases_column)
    #add columns for every existing mouse in this batch
    for(mouse in unique_mice){
      result_total_movement[[mouse]] <- NA
    }
    #add columns for every system(for general movement in system)
    for(system in uniqueSystems){
      result_total_movement[[system]] <- NA
    }
    
    # tibble for positions:
    #create tibble with Phase column
    result_total_positions <- tibble("Phase"=phases_column)
    #add columns for every system(for general movement in system)
    for(system in uniqueSystems){
      result_total_positions[[system]] <- NA
    }
    
    #tibble of information which mouse is in which system in this particular cagechange
    mouse_which_system_tibble <- tibble("sys.1"=rep(NA, times=4), "sys.2"=rep(NA, times=4), "sys.3"=rep(NA, times=4), "sys.4"=rep(NA, times=4), "sys.5"=rep(NA, times=4))
    
    
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
      
      #add the mice of this system information in the tibble
      mouse_which_system_tibble[[systemName]] <-  c(mouse_names)
      
      ## FOR LOOP ## (difference between 2 phases)
      for(phase in phases){
        
        print(phase) 
        
        #depending on the phase take the number of phases
        ifelse(phase == "Active",  number_of_phases <- active_phases_number,  number_of_phases <- inactive_phases_number)
        print(number_of_phases)
        
        ## FOR LOOP ## (difference between the existing number of phases)
        for(nr in number_of_phases){
          print(paste0("System: ", systemName, ", ", phase, " phase nr: ", nr))
          
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
          
          
          #initialize mice closeness result
          #m1 on the third int means number of seconds together from m1 and m3
          count_closeness_list <- list(   m1=c(0,0,0,0),
                                          m2=c(0,0,0,0),
                                          m3=c(0,0,0,0),
                                          m4=c(0,0,0,0))
          #initialize mice closeness result
          #c(positionID, number of seconds)
          count_position_list <- list( c(1, 0), 
                                       c(2, 0), 
                                       c(3, 0), 
                                       c(4, 0), 
                                       c(5, 0), 
                                       c(6, 0), 
                                       c(7, 0), 
                                       c(8, 0))
          
          #initialize total closeness result
          #c(name, total seconds of closeness)
          total_closeness_list <- list( c(mouse_names[1], 0), 
                                        c(mouse_names[2], 0),
                                        c(mouse_names[3], 0),
                                        c(mouse_names[4], 0))
          
          #initialize movement list
          #c(name, number of movements in phase)
          count_movement_list <- list( c(mouse_names[1], 0), 
                                       c(mouse_names[2], 0),
                                       c(mouse_names[3], 0),
                                       c(mouse_names[4], 0),
                                       c(systemName, 0))
          
          
          
          ##CALCULATIONS##
          message("calculates")
          
          ## Definitions for closeness and postition analysis ##
          
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
          
          #define end for while loop
          theEnd <- nrow(systemPhaseData)+1
          
          ## WHILE LOOP ## (goes through the rows of current systemPhaseData)
          while(lineTemp!=theEnd && lineTemp<theEnd){
            
            #make a copy of the old version of the mice list to compare with the new list(differences between two consecutive lines in data)
            old_mice_list <- mice_list
            
            #update mice list with new positions from the next row
            mice_list <- update_mice_list(mouse_names, mice_list, systemPhaseData, timeTemp, lineTemp)
            
            #update secTemp
            secTemp <-  mice_list[["tempData"]][["secTemp"]]
            
            ##################################################################################################
            
            ## use ANALYSIS FUNCTIONS to update result lists ##
            
            ## social closeness between each mouse
            if(system_complete){count_closeness_list <- check_closeness(old_mice_list,mice_list,count_closeness_list, secTemp)}
            ## usage of cage locations
            count_position_list <- check_position(old_mice_list,mice_list,count_position_list, secTemp)
            ## total social closeness of one individual mouse
            if(system_complete){total_closeness_list <- check_total_closeness(old_mice_list, mice_list, total_closeness_list, secTemp)}
            ## number of coil crossings for each individual mouse
            count_movement_list <- check_movement(old_mice_list,mice_list,count_movement_list, secTemp)
            
            ##################################################################################################
            
            
            #update lineTemp
            lineTemp <- mice_list[["tempData"]][["lineTemp"]]
            #update timeTemp
            timeTemp <- mice_list[[1]][[2]]
          }## END WHILE LOOP ##(all rows of current systemPhaseData)
          
          ## RESULTS ##
          
          
          
          if(FALSE){
          #print list results
          message("results: ")
            
          message("count_closeness_list")
          print(count_closeness_list)
          message("count_position_list")
          print(count_position_list)
          message("total_closeness_list")
          print(total_closeness_list)
          }
          
          
          #print(paste("heatmap", systemName))
          message("generate heatmaps")
          
          #generate heatmaps
          if(system_complete){heatmap_closeness <- generateHeatMapCloseness(count_closeness_list, batch, cageChange, systemName, mouse_names, phase, nr)}
          heatmap_positions <- generateHeatMapPositions(count_position_list, batch, cageChange, systemName, phase, nr)
          #print(heatmap_positions)
          
          #save heatmaps in specific system list
          #social prox
          if(system_complete){systemHeatmaps_closeness <- c(systemHeatmaps_closeness, list(heatmap_closeness))}
          #positions
          systemHeatmaps_positions <- c(systemHeatmaps_positions, list(heatmap_positions))
          
          if(system_complete){
            message("enter total_closeness data in tibble")
            #enter total_closeness_list-information in big result tibble for every phase and system(later used for plotting)
            for(i in 1:4){#for every mouse in the system
              
              mouse <- total_closeness_list[[i]][[1]] #the current mouse
              p <- paste0(substr(phase, 1, 1),nr) #the current phase
              
              #find the row with the current phase we are in
              row <- which(result_total_closeness$Phase == p)
              result_total_closeness[[mouse]][[row]] <- total_closeness_list[[i]][[2]]
            }
          }
          
          
          message("enter count_movement data in tibble")
          #enter count_movement_list-information in big result tibble for every phase and system(later used for plotting)
          p <- paste0(substr(phase, 1, 1),nr) #the current phase
          row <- which(result_total_movement$Phase == p)#find the row with the current phase we are in
          #->mice
          for(i in 1:4){#for every mouse in the system
            mouse <- count_movement_list[[i]][[1]] #the current mouse
            if(!is.na(mouse)){
              result_total_movement[[mouse]][[row]] <- count_movement_list[[i]][[2]]
            }else{cat("mouse value from ", systemName, " in ", batch, " in ", cageChange, "is not available.")}
            
          }
          #->system
          result_total_movement[[systemName]][[row]] <- count_movement_list[[5]][[2]]
          
          #####################
#          message("enter count_position data in tibble")
#          #enter count_position_list-information in big result tibble for every phase and system(later used for plotting)
#          p <- paste0(substr(phase, 1, 1),nr) #the current phase
#          row <- which(result_total_positions$Phase == p)#find the row with the current phase we are in
#          #->mice
#          for(i in 1:4){#for every mouse in the system
#            mouse <- count_movement_list[[i]][[1]] #the current mouse
#            if(!is.na(mouse)){
#              result_total_positions[[mouse]][[row]] <- count_movement_list[[i]][[2]]##
#            }else{cat("mouse value from ", systemName, " in ", batch, " in ", cageChange, "is not available.")}
#            
#          }
#          #->system
#          result_total_positions[[systemName]][[row]] <- count_movement_list[[5]][[2]]###
          #########################
          
        }## END FOR LOOP ##(nr of phases)
        
        
        
        
      }## END FOR LOOP ##(phases)
      
      
      if(system_complete){
        message("generate total closeness graph")  
        #generate total closeness plot for every system
        #system, mousenames, result tibble
        plot_total_closeness <- generateGraph(result_total_closeness, batch, cageChange, mouse_names, systemName)
        
        #add generated plot into list of plots
        all_plots_total_closeness <- c(all_plots_total_closeness, list(plot_total_closeness) )
      }
      
      #save all heatmaps from one system in general list
      #social prox
      if(system_complete){allHeatmaps_closeness <- c(allHeatmaps_closeness, list(systemHeatmaps_closeness))}
      #positions
      #allHeatmaps_positions <- c(allHeatmaps_positions, list(systemHeatmaps_positions))
      message("save heatmaps from one system in general list")
      print(batch)
      allHeatmaps_positions[[batch]][[cageChange]][[systemName]] <- systemHeatmaps_positions
      
      
      #clear the lists of the current system
      systemHeatmaps_closeness <- list()
      systemHeatmaps_positions <- list()
    }
    ## END FOR LOOP ##(systems)
    
    
    
    
    
    ############## saving results ##########################################################################
    
    
    saving_directory <- "S:/Lab_Member/Anja/Git/MDC_Bachelor/E9_SIS_AnimalPos"
    #saving_directory <- "/home/anja/Dokumente/FU BERLIN/BA/Git/MDC_Bachelor/E9_SIS_AnimalPos"
    plots_directory <- "/plots"
    #plots_directory <- "/R-plots"
    tables_directory <- "/tables"
    
    
    #plots
    if(save_plots==TRUE){
      message("save plots")
      
      
      ###total_closeness###
      for (i in seq_along(all_plots_total_closeness)) {
        print(i)
        ggsave(filename = paste0(saving_directory, plots_directory, "/total_closeness","_", batch, "_", cageChange, "_sys.",  i, ".png"), plot = all_plots_total_closeness[[i]], width = 5, height = 2)  
        
      }
      
      for (i in seq_along(allHeatmaps_closeness)) {               #for every system
        print(i)
        
        ###allHeatmaps_closeness###
        title <- allHeatmaps_closeness[[i]][[1]]$labels$title     #we need the information from one title of the system
        #find substring of current title which contains the system information
        pattern <- "sys.."
        system_substring <- (ifelse((match <- regexec(pattern, title))[[1]][1] > 0, regmatches(title, match)[[1]], "Pattern not found."))
        #save
        ggsave(filename = paste0(saving_directory, plots_directory, "/allHeatmaps_closeness_", batch, "_", cageChange, "_", system_substring, ".png"), plot = gridExtra::arrangeGrob(grobs = allHeatmaps_closeness[[i]], ncol = 2, layout_matrix = rbind(c(1,5), c(2,6), c(3,7), c(4,8), c(NA,9))), width = 12, height = 8)  
        
        
        ###allHeatmaps_positions####
        title <- allHeatmaps_positions[[i]][[1]]$labels$title 
        system_substring <- (ifelse((match <- regexec(pattern, title))[[1]][1] > 0, regmatches(title, match)[[1]], "Pattern not found."))
        ggsave(filename = paste0(saving_directory, plots_directory, "/allHeatmaps_positions_", batch, "_", cageChange, "_", system_substring, ".png"), plot = gridExtra::arrangeGrob(grobs = allHeatmaps_positions[[i]], ncol = 2, layout_matrix = rbind(c(1,5), c(2,6), c(3,7), c(4,8), c(NA,9))), width = 12, height = 8)  
      }
      
    }
    
    #tables as csv data
    if(save_tables==TRUE){
      message("save tables")
      
      #total closeness
      write.csv(result_total_closeness, file = paste0(saving_directory, tables_directory,"/", batch, "_", cageChange, "_total_closeness_table.csv"), row.names = FALSE)
      #total movement
      write.csv(result_total_movement, file = paste0(saving_directory, tables_directory,"/", batch, "_", cageChange, "_total_movement_table.csv"), row.names = FALSE)
      #which mouse in which system
      write.csv(mouse_which_system_tibble, file = paste0(saving_directory, tables_directory,"/", batch, "_", cageChange, "_mouse_which_system_tibble.csv"), row.names = FALSE)
    }
    
    
  }
}

############### show plots in R ########################################################################
if(show_plots==TRUE){
  message("show Plots")
  # Create a grid of the total closeness plots
  gridExtra::grid.arrange(grobs = all_plots_total_closeness, ncol = 2)
  
  # Create a grid of the Heatmaps of each system for each analysis
  for (batch in seq_along(allHeatmaps_positions)) {
    for(cc in seq_along(allHeatmaps_positions[[batch]])){
      for(sys in seq_along(allHeatmaps_positions[[batch]][[cc]])){
        print(paste(sys))
        
        if(cc != "CC4") {
          #gridExtra::grid.arrange(grobs = allHeatmaps_closeness[[batch]][[cc]][[sys]], ncol = 2, layout_matrix = rbind(c(1,8), c(2,5), c(3,6), c(4,7)))
          gridExtra::grid.arrange(grobs = allHeatmaps_positions[[batch]][[cc]][[sys]], ncol = 2, layout_matrix = rbind(c(1,8), c(2,5), c(3,6), c(4,7)))
        }else{
          #gridExtra::grid.arrange(grobs = allHeatmaps_closeness[[batch]][[cc]][[sys]], ncol = 2, layout_matrix = rbind(c(1,10), c(2,6), c(3,7), c(4,8), c(5,9)))
          gridExtra::grid.arrange(grobs = allHeatmaps_positions[[batch]][[cc]][[sys]], ncol = 2, layout_matrix = rbind(c(1,10), c(2,6), c(3,7), c(4,8), c(5,9)))
        }
      }
      
    }
  }
}



