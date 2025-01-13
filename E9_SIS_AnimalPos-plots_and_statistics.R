# This script performs data analysis and visualization for the E9 Social Stress experiment.
# It loads necessary libraries, sets up directories for saving plots and statistics, 
# and defines paths for working and saving directories.
# Load required libraries using pacman. If pacman is not installed, install it first.
# The libraries used include readr, tibble, dplyr, reshape2, stringr, ggplot2, forcats, plotrix, and lmerTest.
# Set flags to determine whether to save plots and statistics.
# Define paths for the working directory and saving directory.
# Uncomment and modify the paths as needed for different environments (e.g., local machine, server).
# Define subdirectories for saving plots and tables.

# Load packages using pacman.
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readr, tibble, dplyr, reshape2, stringr, ggplot2, forcats, plotrix, lmerTest, tidyr, gridExtra, ggpubr, cowplot, writexl)

save_plots = TRUE
save_statistics= TRUE

# Paths
working_directory <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability"
saving_directory <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability"
plots_directory <- "/plots"
tables_directory <- "/tables"
consec_plot_dir <- paste0(working_directory, "/plots/consec_plots")
consec_lme_result_dir <- paste0(working_directory, "/lme/consec-mixed_model_results")
cc_lme_result_dir <- paste0(working_directory, "/lme/cc-mixed_model_results")

dirs <- c(consec_plot_dir, consec_lme_result_dir, cc_lme_result_dir)
lapply(dirs, function(dir) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
})

# Load external functions script
source(paste0("C:/Users/topohl/Documents/GitHub/MMMSociability/E9_SIS_AnimalPos-functions.R"))

## general definitions and read data ###########################################################
# Define batches and cage changes
batches <- c("B1", "B2", "B3", "B4", "B5")
cageChanges <- c("CC1", "CC2", "CC3", "CC4") 

# Read con and sus animals
con_animals <- readLines(paste0(working_directory,"/raw_data/con_animals.csv"))
sus_animals <- readLines(paste0(working_directory,"/raw_data/sus_animals.csv"))

####################################################################################
#### CONSEC TIBBLES ####
# create from analysis data similar tibble with a consecutive Phase row
# nice for lme statistics and consecutive plots
####################################################################################
message("## create consec phase mean tibble ##")

###### SOCIAL  PROXIMITY AND COIL CROSSING ######
#initialize empty tibble
social_prox_consec <- tibble(CageChange = character(),
                             Batch = character(),
                             System = character(),
                             AnimalID = character(),
                             Sex = character(),
                             Group = character(),
                             Phase = character(),
                             ConsecActive = numeric(),
                             ConsecInactive = numeric(),
                             SocialProx = numeric())

coil_crossing_consec <- tibble(CageChange = character(),
                             Batch = character(),
                             System = character(),
                             AnimalID = character(),
                             Sex = character(),
                             Group = character(),
                             Phase = character(),
                             ConsecActive = numeric(),
                             ConsecInactive = numeric(),
                             CoilCrossing = numeric())

tibbles_to_process <- c("total_closeness_table", "total_movement_table")
values <- c('SocialProx', 'CoilCrossing')

for (i in seq_along(tibbles_to_process)) {
  # Define the current tibble name
  tibblename <- tibbles_to_process[i]
  
  for (batch in batches) {
    max_consecAct <- 0
    max_consecInact <- 0
    
    for (cageChange in cageChanges) {
      sex <- ifelse(batch %in% c("B3", "B4", "B6"), "female", "male")
      
      # Load CSV files for each cage change
      system_animal_ids <- tibble(read_delim(paste0(working_directory, "/tables/", batch, "_", cageChange, "_animal_ids.csv"), delim = ",", show_col_types = FALSE))
      if (tibblename == "total_closeness_table") {
        cc_tibble <- tibble(read_delim(paste0(working_directory, "/tables/", batch, "_", cageChange, "_total_proximity.csv"), delim = ",", show_col_types = FALSE))
      } else if (tibblename == "total_movement_table") {
        cc_tibble <- tibble(read_delim(paste0(working_directory, "/tables/", batch, "_", cageChange, "_total_movement.csv"), delim = ",", show_col_types = FALSE))
      } else {
        stop("Invalid tibble name")
      }
      
      # Special handling for CC4 after A2 (grid in cage) - not usable for social analysis
      if (cageChange == "CC4") {
        cc_tibble <- cc_tibble %>% filter(Phase %in% c("A1", "I2", "A2"))
      }
      
      # Preprocess the data
      cc_tibble <- melt(cc_tibble, id = 'Phase')
      names(cc_tibble) <- c('Phase', 'AnimalID', values[i])
      
      # Create consecutive columns
      cc_tibble[c('Phase', 'Consec')] <- str_split_fixed(cc_tibble$Phase, '', 2)
      cc_tibble <- cc_tibble %>%
        mutate(ConsecActive = ifelse(Phase == "A", as.numeric(Consec), 0)) %>%
        mutate(ConsecInactive = ifelse(Phase == "I", as.numeric(Consec), 0)) %>%
        mutate(Phase = ifelse(Phase == "A", "active", "inactive"))
      
      if (tibblename == "total_closeness_table") {
        cc_tibble <- cc_tibble %>% mutate(SocialProx = ifelse(SocialProx != 0, SocialProx / 3600, SocialProx))
      }
      
      # Adjust consecutive values for higher cage changes
      if (cageChange != "CC1") {
        batch_table <- if (tibblename == "total_closeness_table") {
          social_prox_consec %>% filter(Batch == batch)
        } else if (tibblename == "total_movement_table") {
          coil_crossing_consec %>% filter(Batch == batch)
        } else {
          stop("Invalid tibble name")
        }
        
        max_consecAct <- batch_table %>% pull(ConsecActive) %>% unique() %>% max()
        max_consecInact <- batch_table %>% pull(ConsecInactive) %>% unique() %>% max()
        
        cc_tibble <- cc_tibble %>%
          mutate(ConsecActive = ifelse(Phase == "active", max_consecAct + ConsecActive, 0)) %>%
          mutate(ConsecInactive = ifelse(Phase == "inactive", max_consecInact - 1 + ConsecInactive, 0))
      }
      
      # Add new columns to the tibble
      cc_tibble <- cc_tibble %>%
        mutate(Group = ifelse(AnimalID %in% sus_animals, "sus", ifelse(AnimalID %in% con_animals, "con", "res"))) %>%
        mutate(Sex = sex) %>%
        mutate(Batch = batch) %>%
        mutate(AnimalID = as.character(AnimalID)) %>%
        mutate(System = purrr::map_chr(AnimalID, ~ {
          index <- which(system_animal_ids == .x, arr.ind = TRUE)[2]
          if (length(index) > 0) { 
            names(system_animal_ids)[index]
            } else { 
              NA 
              }
        })) %>%
        mutate(CageChange = cageChange)
      
      # Reorder columns
      cc_tibble <- cc_tibble[c('CageChange', 'Batch', 'System', 'AnimalID', 'Sex', 'Group', 'Phase', 'ConsecActive', 'ConsecInactive', values[i])]
      
      # Combine with the general tibble for all cage changes and batches
      if (tibblename == "total_closeness_table") {
        social_prox_consec <- bind_rows(social_prox_consec, cc_tibble)
      } else if (tibblename == "total_movement_table") {
        coil_crossing_consec <- bind_rows(coil_crossing_consec, cc_tibble)
      }
    }
  }
  
  # Remove rows with NA values
  if (tibblename == "total_closeness_table") {
    social_prox_consec <- na.omit(social_prox_consec)
  } else if (tibblename == "total_movement_table") {
    coil_crossing_consec <- na.omit(coil_crossing_consec)
  }
}

###### CAGE AND MICE ENTROPY ######
# read saved tables into a tibble 
cage_entropy_consec <- as_tibble(read_delim(paste0(saving_directory, tables_directory,"/all_batches", "_all_cageChanges", "_consec_cage_entropy.csv"),delim = ",", show_col_types = FALSE))
mice_entropy_consec <- as_tibble(read_delim(paste0(saving_directory, tables_directory,"/all_batches", "_all_cageChanges", "_consec_mice_entropy.csv"),delim = ",", show_col_types = FALSE))

##all tibbles as a list together##
consec_tibbles <- list(social_prox_consec, coil_crossing_consec, mice_entropy_consec, cage_entropy_consec)

###############################################################################################################################################################
#### LME AND PLOTS OF CONSEC TIBBLES ####
###############################################################################################################################################################
message("## create consec phase mean plots and statistics ##")

#dynamic plot and lme for all consec tibbles 

#definitions for dynamic processing
values <- c('SocialProx', 'CoilCrossing', 'MiceEntropy', 'CageEntropy')
phases <- c("active", "inactive")
consec_phases <- c('ConsecActive', 'ConsecInactive')

# Initialisieren einer Liste für die Ergebnisse
#statistics <- list()
plots <- list()

for(i in seq_along(consec_tibbles)){
  value <- values[i]
  color_by <- ifelse(value=='CageEntropy', 'System', 'Group')
  ifelse(value == 'CageEntropy', color_palette <- c("forestgreen", "goldenrod", "steelblue", "firebrick", "darkorchid"), color_palette <- c("con" = "#1e3791", "res" = "#8aacdb", "sus" = "#f49620"))
  
  message(value)
  
  for(j in seq_along(phases)){
    consec_data <- consec_tibbles[[i]]
    #print(consec_data)
    phase <- phases[j]
    consec_phase <- consec_phases[j]
    
    message(phase)
    message(consec_phase)
    
    #divide into act or inact phases
    consec_data <- consec_data %>%
      filter(Phase == phase)
    
    ## PLOT  ##
    p <- ggplot(data = consec_data, aes(x = fct_inorder(as_factor(.data[[consec_phase]])), y = .data[[value]], color = !!sym(color_by), group = !!sym(color_by))) + 
      stat_summary(aes(group = !!sym(color_by)), 
                   fun = median, 
                   geom = "line", 
                   linewidth = 1) +
      stat_summary(aes(fill =!!sym(color_by)), 
                   fun.min = function(z){quantile(z, 0.25)},
                   fun.max = function(z){quantile(z, 0.75)},
                   fun = median,
                   geom = "ribbon",                              
                   alpha = 0.2, color = NA) +
      labs(title = paste(value), subtitle = paste(phase, "phases")) +
      scale_fill_manual(values = color_palette) +
      scale_color_manual(values = color_palette) +
      scale_y_continuous(paste(value)) +
      scale_x_discrete(name = consec_phase) + 
      facet_grid(Sex~.) +
      theme_minimal(base_size = 14) +
      theme(
            panel.grid.major = element_blank(),  # Remove major gridlines
            panel.grid.minor = element_blank(),  # Remove minor gridlines
            plot.title = element_text(hjust = 0.5, face = "bold", size = 18),  # Center title
            plot.subtitle = element_text(hjust = 0.5, size = 14, face = "italic"),  # Italic for subtitle
            axis.text.x = element_text(size = 12),  # Adjust text sizes for clarity
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 14, face = "bold"),
            axis.ticks.x = element_line(size = 0.5),
            legend.position = "top",  # Place legend inside the plot area
            panel.background = element_blank())  # Clean panel background
    
    #save plot in plot list
    if(j == 1) plots[[2 * i - 1]] <- p
    else plots[[2 * i]] <- p
    
    #save plot in file
    if(save_plots) ggsave(filename = paste0(working_directory, "/plots", "/consec_plots/", "consec-",value,"-",consec_phase, ".svg"), width = 5, height = 5)
    
    # Create dynamic formula for linear mixed effects model
    if(value=='CageEntropy') {
      formula <- as.formula(paste(value, "~ System * Sex + ", consec_phase, " + (1 | Batch)"))
    } else {
      formula <- as.formula(paste(value, "~ Group * Sex + ", consec_phase, " + (1 | AnimalID)"))
    }

    # Perform linear mixed effects model
    message(paste("consec lme for", value,", ", phase))
    print(formula)
    model <- lmerTest::lmer(formula, data = consec_data)
    print(summary(model))
    
    # save mixed model results to a file
    mixed_model_results_file <- paste0(working_directory, "/lme/consec-mixed_model_results/", paste0("lme-",value,"-",consec_phase, ".txt"))
    model_summary <- capture.output(summary(model))
    if(save_statistics) writeLines(model_summary, con = mixed_model_results_file)
  }
}

## PRINT PLOTS IN R ##
# Print all generated plots
#separate
for (i in seq_along(plots)) {
  print(plots[[i]])
}
#together
gridExtra::grid.arrange(grobs = plots, ncol = 2)

#####################################################################################################################################
## CageChange Tibble ##
#####################################################################################################################################
#plots of every analysis in perspective of every cagechange
message("## create cage change mean tibble ##")

#initialize rank tibble
social_prox_CC <- tibble(CageChange = character(),
                      System = character(),
                      AnimalID = character(),
                      Group = character(),
                      Batch = character(),
                      Sex = character(),
                      Sum_act = numeric(),
                      Sum_inact = numeric(),
                      Avg_act = numeric(),
                      Avg_inact = numeric(),
                      SystemRank_act = numeric(),
                      SystemRank_inact = numeric())

coil_crossing_CC <- tibble(CageChange=character(),
                         System=character(),
                         AnimalID=character(),
                         Group=character(),
                         Batch=character(),
                         Sex=character(),
                         Sum_act=numeric(),
                         Sum_inact=numeric(),
                         Avg_act=numeric(),
                         Avg_inact=numeric(),
                         SystemRank_act=numeric(),
                         SystemRank_inact=numeric())

mice_entropy_CC <- tibble(CageChange=character(),
                          System=character(),
                          AnimalID=character(),
                          Group=character(),
                          Batch=character(),
                          Sex=character(),
                          Sum_act=numeric(),
                          Sum_inact=numeric(),
                          Avg_act=numeric(),
                          Avg_inact=numeric(),
                          SystemRank_act=numeric(),
                          SystemRank_inact=numeric())

cage_entropy_CC <- tibble(CageChange=character(),
                          System=character(),
                          Batch=character(),
                          Sex=character(),
                          Sum_act=numeric(),
                          Sum_inact=numeric(),
                          Avg_act=numeric(),
                          Avg_inact=numeric(),
                          SystemRank_act=numeric(),
                          SystemRank_inact=numeric())

#definitions for dynamic processing
cc_tibbles <- list('SocialProx' = social_prox_CC, 'CoilCrossing' = coil_crossing_CC, 'MiceEntropy' = mice_entropy_CC, 'CageEntropy' = cage_entropy_CC)
values <- c('SocialProx', 'CoilCrossing', 'MiceEntropy', 'CageEntropy')

for(i in c(1,2,3,4)){
  consec_tibble <- consec_tibbles[[i]]
  #print(consec_tibble)
  value <- values[i]
  message(value)
  cc_tibble <- cc_tibbles[[value]]
  
  for(batch in batches){
    for(cageChange in cageChanges){
      message(paste(batch, cageChange))
     
      filtered_consec_tibble <- consec_tibble %>%
        filter(Batch==batch) %>%
        filter(CageChange==cageChange)
      
      #load csv
      system_animal_ids <-  tibble(read_delim(paste0(working_directory,"/tables/", batch, "_", cageChange, "_animal_ids.csv"), delim = ",", show_col_types = FALSE))
      
      for(system in colnames(system_animal_ids)){
        
        #message(system)
        
        animal_ids <- system_animal_ids %>%  # this will create a character vector
          select(all_of(system)) %>%
          pull()
        
        #if system is incomplete->skip system
        if(NA %in% animal_ids){
          cat("skip ", system, " in ", cageChange, "in", batch, "\n")
          #skip for loop
          next
        }
        
        #empty vector for determining rank:
        #(max value,active phase, of a mouse of a system/of a system)
        act_system_vec <- c()
        #(max value,inactive phase, of a mouse of a system/of a system)
        inact_system_vec <- c()
        
        #calculating/generating values of each column and add new row to tibble
        if(value == 'CageEntropy'){
          #message("columns for cageEntropy")
          sex <- ifelse(batch %in% c("B3", "B4", "B6"), "female", "male")
          
          #filter into act and inact phases of system and select all existing values
          act_tibble <- filtered_consec_tibble %>%
            filter(Phase == "active") %>%
            filter(System==system) %>%
            select(all_of(value))
          
          inact_tibble <- filtered_consec_tibble %>%
            filter(Phase == "inactive") %>%
            filter(System==system) %>%
            select(all_of(value))
          
          #sum of multiple phases of one CageChange
          sum_act <- act_tibble %>%
            sum()
          #cat("sum of active: ", sum_act, "\n")
          
          sum_inact <- inact_tibble %>%
            sum()
          #cat("sum of inactive: ", sum_inact, "\n")
          
          avg_act <- act_tibble %>%
            unlist() %>%
            mean()
          #cat("mean of inactive: ", avg_act, "\n")
          
          avg_inact <- inact_tibble %>%
            unlist() %>%
            mean()
          #cat("mean of inactive: ", avg_inact, "\n")
          
          #add row to tibble
          cc_tibble <- cc_tibble %>% 
            add_row(CageChange = cageChange,
                    System = system,
                    Batch = batch,
                    Sex = sex,
                    Sum_act = sum_act,
                    Sum_inact = sum_inact,
                    Avg_act = avg_act,
                    Avg_inact = avg_inact,
                    SystemRank_act = NA,
                    SystemRank_inact = NA)
          
          #FOR RANK COMPUTING:
          #add total hours to vector
          act_system_vec <- append(act_system_vec, sum_act)
          inact_system_vec <- append(inact_system_vec, sum_inact)
       
        } else {
          for(animal_id in animal_ids) {

            ## create variables for new row of tibble ##
            group <- ifelse(animal_id %in% sus_animals, "sus", ifelse(animal_id %in% con_animals, "con", "res"))
            #print(group)
            sex <- ifelse(batch %in% c("B3", "B4", "B6"), "female", "male")
            
            #filter into act and inact phases of the individual and select all existing values
            act_tibble <- filtered_consec_tibble %>%
              filter(Phase == "active") %>%
              filter(AnimalID == animal_id) %>%
              select(all_of(value))
            
            inact_tibble <- filtered_consec_tibble %>%
              filter(Phase == "inactive") %>%
              filter(AnimalID == animal_id) %>%
              select(all_of(value))
            
            #sum of multiple phases of one CageChange
            sum_act <- act_tibble %>%
              sum()
            #cat("sum of active: ", sum_act, "\n")
            
            sum_inact <- inact_tibble %>%
              sum()
            #cat("sum of inactive: ", sum_inact, "\n")
            
            avg_act <- act_tibble %>%
              unlist() %>%
              mean()
            #cat("mean of inactive: ", avg_act, "\n")
            
            avg_inact <- inact_tibble %>%
              unlist() %>%
              mean()
            #cat("mean of inactive: ", avg_inact, "\n")
            
            #add row to tibble
            cc_tibble <- cc_tibble %>% 
              add_row(CageChange = cageChange,
                      System = system,
                      AnimalID = animal_id,
                      Group = group,
                      Batch = batch,
                      Sex = sex,
                      Sum_act = sum_act,
                      Sum_inact = sum_inact,
                      Avg_act = avg_act,
                      Avg_inact = avg_inact,
                      SystemRank_act = NA,
                      SystemRank_inact = NA)
            
            #FOR RANK COMPUTING:
            #add total hours to vector
            act_system_vec <- append(act_system_vec, sum_act)
            inact_system_vec <- append(inact_system_vec, sum_inact)
          }
        }
        #RANK COMPUTING:
        cc_tibble <- compute_rank(cc_tibble, act_system_vec, system, cageChange, batch,'Sum_act', 'SystemRank_act')
        cc_tibble <- compute_rank(cc_tibble, inact_system_vec, system, cageChange, batch, 'Sum_inact', 'SystemRank_inact')
      }
    }
  }
  #overwrite empty tibble in list with filled one 
  cc_tibbles[[value]] <- cc_tibble
}

###############################################################################################################################################################
#### STATISTIC AND PLOTS OF CC TIBBLES ####
###############################################################################################################################################################
message("## create cage change mean plots and statistics ##")

#dynamic plot and statistic for all CC tibbles 

#definitions for dynamic processing
values <- c('SocialProx', 'CoilCrossing', 'MiceEntropy', 'CageEntropy')
phases <- c("active", "inactive")

# Initialisieren einer Liste für die Ergebnisse
#avg hors per cc
avg_plots <- list()
#avg rank in total
rank_plots <- list()

##statistics on rank plots(tobis code, plots again)
allTestResults <- list()
allPlots <- list()
allPosthocResults <- list()

for(i in seq_along(cc_tibbles)){
  value <- values[i]
  color_by <- ifelse(value == 'CageEntropy', 'System', 'Group')
  ifelse(value == 'CageEntropy', color_palette <- c("forestgreen", "goldenrod", "steelblue", "firebrick", "darkorchid"), color_palette <- c("con" = "#1e3791", "res" = "#8aacdb", "sus" = "#f49620"))
  
  message(value)
  
  for(j in seq_along(phases)){
    cc_data <- cc_tibbles[[value]]
    #print(cc_data)
    phase <- phases[j]
    message(phase)
    
    ## PLOTS  ##
    ## AVG PLOT
    avg_phase <- ifelse(phase == "active", 'Avg_act', 'Avg_inact')
    
    #style ribbon
    avg_cc_p <- ggplot(data = cc_data, aes(x = CageChange, y = !!sym(avg_phase), color = !!sym(color_by), group=!!sym(color_by))) + 
      stat_summary(aes(group = !!sym(color_by)), fun = median, geom = "line", linewidth = 1) +
      stat_summary(aes(fill = !!sym(color_by)),
                   fun.min = function(z){quantile(z, 0.25)},
                   fun.max = function(z){quantile(z, 0.75)},
                   fun = median,
                   geom = "ribbon",                              
                   alpha = 0.2, color = NA) +
      labs(title = paste(value, " Cage Change Plot"), subtitle = phase, "phases") +
      scale_fill_manual(values = color_palette) +
      scale_color_manual(values = color_palette) +
      scale_y_continuous(paste(value)) +
      scale_x_discrete(name = "Cage Change")+ 
      facet_grid(Sex~.) +
      theme_minimal(base_size = 14) +  # Increase base font size for readability
      theme(
            panel.grid.major = element_blank(),  # Remove major gridlines
            panel.grid.minor = element_blank(),  # Remove minor gridlines
            plot.title = element_text(hjust = 0.5, face = "bold", size = 18),  # Center title
            plot.subtitle = element_text(hjust = 0.5, size = 14, face = "italic"),  # Italic for subtitle
            axis.text.x = element_text(size = 12),  # Adjust text sizes for clarity
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 14, face = "bold"),
            axis.ticks.x = element_line(size = 0.5),
            legend.position = "top",  # Place legend inside the plot area
            panel.background = element_blank())  # Clean panel background
    
    #save plot in plot list
    if(j == 1) avg_plots[[2 * i - 1]] <- avg_cc_p
    else avg_plots[[2 * i]] <- avg_cc_p
    #save plot in file
    if(save_plots) ggsave(filename = paste0(working_directory, "/plots", "/cc_plots/", "avg-cc-",value,"-",phase,".svg"), plot = avg_cc_p, width = 5, height = 5)
    
    # LME-AVG 
    message(paste("cc lme for", value ,", ", avg_phase))
    
    #create dynamic formula
    if(value=='CageEntropy') {
      formula <- as.formula(paste(avg_phase, "~ System * Sex + CageChange + (1 | Batch)")) ####maybe * CageChange?
    } else {
      formula <- as.formula(paste(avg_phase, "~ Group * Sex + CageChange + (1 | AnimalID)"))
    }
    message(formula)
    
    model <- lmerTest::lmer(formula, data = cc_data)
    #print(summary(model))
    
    # save mixed model results to a file
    mixed_model_results_file <- paste0(working_directory, "/lme/cc-mixed_model_results/", paste0("avg-lme-",value,"-",phase,".txt"))
    model_summary <- capture.output(summary(model))
    if(save_statistics) writeLines(model_summary, con = mixed_model_results_file)
    if(value!='CageEntropy'){

      ## AVG RANK PLOT
      #define the current rank-value based o the current phase
      rank_phase <- ifelse(phase == "active", 'SystemRank_act', 'SystemRank_inact')
      
      #group data to summarize all Cage changes for each individual
      avg_data <- cc_data %>%
        group_by(Batch, Sex, Group, AnimalID) %>%
        summarise(avg_rank = mean(!!sym(rank_phase)))
      
      #exclude con group from visualisation if  socialProx (here not valuable)
      if(value=='SocialProx') avg_data <- avg_data %>% filter(Group!="con")
      
      avg_rank_p <- ggplot(data = avg_data, aes(x = Group, y = avg_rank, color = Group, shape = Batch)) +
      geom_jitter(size = 4, alpha = 0.7, width = 0.2, height = 0) +
        scale_shape_manual(values = c(3, 16, 17, 15, 5, 20)) + 
        scale_x_discrete("Animal ID") +
        scale_y_discrete("Rank", limits = c("1", "2", "3", "4")) +
        labs(title = paste("Median rank of each individual"), subtitle = paste(value, "-", phase)) +
        scale_color_manual(values = c("con" = "#1e3791", "res" = "#8aacdb", "sus" = "#f49620")) +
        stat_summary(
          fun.min = function(z) {quantile(z, 0.25)},
          fun.max = function(z) {quantile(z, 0.75)},
          fun = median,
          color = "black",
          size = 0.8,
          shape = 16) +
        facet_grid(Sex ~ .) +
        theme_minimal(base_size = 14) +  # Increase base font size for readability
        theme(
              panel.grid.major = element_blank(),  # Remove major gridlines
              panel.grid.minor = element_blank(),  # Remove minor gridlines
              plot.title = element_text(hjust = 0.5, face = "bold", size = 18),  # Center title
              plot.subtitle = element_text(hjust = 0.5, size = 14, face = "italic"),  # Italic for subtitle
              axis.text.x = element_text(size = 12),  # Adjust text sizes for clarity
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 14, face = "bold"),
              axis.ticks.x = element_line(size = 0.5),
              legend.position = "top",  # Place legend inside the plot area
              panel.background = element_blank())  # Clean panel background
      
      #save plot in plot list
      if(j == 1) rank_plots[[2 * i - 1]] <- avg_rank_p
      else rank_plots[[2 * i]] <- avg_rank_p
      #save plot in file
      if(save_plots) ggsave(filename = paste0(working_directory, "/plots", "/avg_rank_plots/", "avg_ranks-",value,"-", phase, ".svg"), plot = avg_rank_p, width = 5, height = 5)
      
      ##TOBIS STATISTICS CODE
      for(sex in c("female", "male")){
        
        message("statistics on rank avgs")
        print(sex)
        result <- testAndPlotVariable(avg_data, value, 'avg_rank', phase, sex)#avg_data or cc_data???
        # add result-list(containing the columns testResults, plot, posthocResults) to other fitting list
        if (!is.null(result)) {
          # posthocResults is always NULL for the Wilcoxon test
          if (is.null(result$posthocResults)) {
            allTestResults <- c(allTestResults, list(result$testResults))
          } else {
            allTestResults <- c(allTestResults, list(result$testResults))
            allPosthocResults <- c(allPosthocResults, list(result$posthocResults))
          }
          #add the plot
          allPlots <- c(allPlots, list(result$plot))
        }
      }
    }
  }
}

## PRINT PLOTS IN R ##
# Print all plots individually
for (i in seq_along(avg_plots)) {
  print(avg_plots[[i]])
}
# Arrange plots into grids
gridExtra::grid.arrange(grobs = avg_plots, ncol = 2)   # Average plots
gridExtra::grid.arrange(grobs = rank_plots, ncol = 2)  # Rank plots
gridExtra::grid.arrange(grobs = allPlots, ncol = 4)    # All plots

###################### saving the results #############################################################

# Convert the list of test results to a data frame
allTestResultsDf <- bind_rows(allTestResults)
# Save the test results data frame to a CSV file
# Create the statistics directory if it doesn't exist
statistics_directory <- paste0(working_directory, "/statistics")
if (!dir.exists(statistics_directory)) {
  dir.create(statistics_directory, recursive = TRUE)
}
write.csv(allTestResultsDf, file = paste0(working_directory,"/statistics/test_results.csv"), row.names = FALSE)
# Save the post hoc results to a CSV file
if (!is.null(allPosthocResults) && length(allPosthocResults) > 0) {
  allPosthocResultsDf <- bind_rows(allPosthocResults)
  write.csv(allPosthocResultsDf, file = paste0(working_directory,"/statistics/posthoc_results.csv"), row.names = FALSE)
}