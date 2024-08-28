# libraries
library(readr)        # load readr package for reading csv files
library(tibble)       #important for tibble operations
library(dplyr)        #changes in tibbles
library(reshape2)     #for melt()
library(stringr)
library(ggplot2)
library(forcats)
library(plotrix)      #for std.error
library(lmerTest)     # Perform mixed-effects models with pairwise comparisons

save_plots = FALSE
save_statistics= FALSE


# paths
working_directory <- "S:/Lab_Member/Anja/Git/MDC_Bachelor/E9_SIS_AnimalPos"
#working_directory <- "/home/anja/Dokumente/FU BERLIN/BA/Git/MDC_Bachelor/E9_SIS_AnimalPos"

saving_directory <- "S:/Lab_Member/Anja/Git/MDC_Bachelor/E9_SIS_AnimalPos"
#saving_directory <- "/home/anja/Dokumente/FU BERLIN/BA/Git/MDC_Bachelor/E9_SIS_AnimalPos"
plots_directory <- "/plots"
tables_directory <- "/tables"


#functions
source(paste0(working_directory,"/E9_SIS_AnimalPos-functions.R"))

## general definitions and read data ###########################################################
# define batch and cage change
batches <- c("B1", "B2", "B3", "B4", "B5")
cageChanges <- c("CC1", "CC2", "CC3", "CC4") 


#sus and con animals
sus_animals <- readLines(paste0(working_directory,"/raw_data/sus_animals.csv"))
##csv path for con animals
con_animals <- readLines(paste0(working_directory,"/raw_data/con_animals.csv"))


####################################################################################
#### CONSEC TIBBLES ####
# create from analysis data similar tibble with a consecutive Phase row
# nice for lme statistics and consecutive plots
####################################################################################
message("## create consec phase mean tibble ##")

###### SOCIAL  PROXIMITY AND COIL CROSSING ######
#initialize empty tibble
social_prox_consec <- tibble(CageChange=character(),
                             Batch=character(),
                             System=character(),
                             AnimalID=character(),
                             Sex=character(),
                             Group=character(),
                             Phase=character(),
                             ConsecActive=numeric(),
                             ConsecInactive=numeric(),
                             SocialProx=numeric())

coil_crossing_consec <- tibble(CageChange=character(),
                             Batch=character(),
                             System=character(),
                             AnimalID=character(),
                             Sex=character(),
                             Group=character(),
                             Phase=character(),
                             ConsecActive=numeric(),
                             ConsecInactive=numeric(),
                             CoilCrossing=numeric())


tibbles_to_process <- c("total_closeness_table", "total_movement_table")
values <- c('SocialProx', 'CoilCrossing')

for (i in seq_along(tibbles_to_process)) {
  #define current tibble
  tibblename <- tibbles_to_process[i]
  #print(tibblename)
  
  for(batch in batches){
    
    max_consecAct <- 0
    max_consecInact <-0
    
    for(cageChange in cageChanges){
      
      #define sex
      sex <- ifelse(batch %in% c("B3", "B4", "B6"), "female", "male")
      
      #load csv per cage change
      mouse_which_system_tibble <-  tibble(read_delim(paste0(working_directory,"/tables/", batch, "_", cageChange, "_mouse_which_system_tibble.csv"),delim = ",", show_col_types = FALSE))
      if(tibblename=="total_closeness_table")cc_tibble <- tibble(read_delim(paste0(working_directory,"/tables/", batch, "_", cageChange, "_total_closeness_table.csv"),delim = ",", show_col_types = FALSE))
      else if(tibblename=="total_movement_table")cc_tibble <- tibble(read_delim(paste0(working_directory,"/tables/", batch, "_", cageChange, "_total_movement_table.csv"),delim = ",", show_col_types = FALSE))
      else stop("no such tibble to process")
      
      
      #special treatment for every CC4 after A2(grid in cage)-> not usable for social analysis
      #use only A1,I2,A2 for the computation of the rank from CC4
      if(cageChange=="CC4"){
        cc_tibble <- cc_tibble%>%
          filter(Phase %in% c("A1","I2","A2"))
      }
      
      #PREPROCESS
      # melt tibble
      cc_tibble<- melt(cc_tibble, id='Phase')
      names(cc_tibble) <- c('Phase', 'AnimalID', values[i])
      
      
      #create consec colums:
      #split Phase column
      cc_tibble[c('Phase', 'Consec')] <- str_split_fixed(cc_tibble$Phase, '', 2)
      
      #change consec column to ConsecAct and Inact
      #and change content of Phase column to active and inactive
      cc_tibble <- cc_tibble%>%
        mutate(ConsecActive = ifelse(Phase=="A", as.numeric(Consec), 0))%>%
        mutate(ConsecInactive = ifelse(Phase=="I",as.numeric(Consec), 0))%>%
        mutate(Phase = ifelse(Phase=="A", "active", "inactive"))
      
      
      if(tibblename=="total_closeness_table"){
        #calculate the counted seconds to hours
        cc_tibble <- cc_tibble%>%
          mutate(SocialProx = ifelse(SocialProx!=0,SocialProx/3600,SocialProx))
      }
      
      
      #count the consec higher for higher cage changes
      if(cageChange!="CC1"){
        
        if(tibblename=="total_closeness_table")batch_table <- social_prox_consec%>%filter(Batch==batch)
        else if(tibblename=="total_movement_table")batch_table <- coil_crossing_consec%>%filter(Batch==batch)
        
        max_consecAct <- batch_table%>%
          pull(ConsecActive)%>%
          unique()%>%
          max()
        max_consecInact <- batch_table%>%
          pull(ConsecInactive)%>%
          unique()%>%
          max()
        
        cc_tibble <- cc_tibble%>%
          mutate(ConsecActive = ifelse(Phase=="active", max_consecAct+ConsecActive, 0))%>%
          mutate(ConsecInactive = ifelse(Phase=="inactive", max_consecInact-1+ConsecInactive, 0))#Inact phase numbers start always with 2...
      }  
      
      #MUTATE NEW COLUMNS  
      cc_tibble <- cc_tibble%>%
        mutate(Group = ifelse(AnimalID %in% sus_animals, "sus", ifelse(AnimalID %in% con_animals, "con", "res")))%>%
        mutate(Sex=sex)%>%
        mutate(Batch=batch)%>%
        mutate(AnimalID = as.character(AnimalID))%>%
        mutate(System = purrr::map_chr(AnimalID, ~ {  #searches system number in the tibble depending on the animalID in the row
          index <- which(mouse_which_system_tibble == .x, arr.ind = TRUE)[2]
          if (length(index) > 0) { names(mouse_which_system_tibble)[index]}
          else {NA}
        }))%>%
        mutate(CageChange=cageChange)
      
      
      #sort tibble into correct order of columns
      cc_tibble <- cc_tibble[c('CageChange', 'Batch', 'System', 'AnimalID', 'Sex', 'Group', 'Phase', 'ConsecActive', 'ConsecInactive', values[i])]
      
      ## combine with general tibble for all CCs and batches ##
      if(tibblename=="total_closeness_table")social_prox_consec <- bind_rows(social_prox_consec,cc_tibble)
      else if(tibblename=="total_movement_table")coil_crossing_consec <- bind_rows(coil_crossing_consec,cc_tibble)
      
    }
  }
  
  
  # delete NA rows in tibble (NA values)
  if(tibblename=="total_closeness_table")social_prox_consec <- na.omit(social_prox_consec)
  else if(tibblename=="total_movement_table")coil_crossing_consec <- na.omit(coil_crossing_consec)
  
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
  ifelse(value=='CageEntropy',manual_colors <- c("cornflowerblue", "yellow", "darkgreen", "orange", "purple"), manual_colors <- c("sus" = "tomato", "res" = "darkgreen", "con" = "deepskyblue4"))
  
  message(value)
  
  for(j in seq_along(phases)){
    consec_data <- consec_tibbles[[i]]
    #print(consec_data)
    phase <- phases[j]
    consec_phase <- consec_phases[j]
    
    message(phase)
    message(consec_phase)
    
    #divide into act or inact phases
    consec_data <- consec_data%>%
      filter(Phase == phase)
    
    
    ## PLOT  ##
    p <- ggplot(data = consec_data, aes(x = fct_inorder(as_factor(.data[[consec_phase]])), y = .data[[value]], color = !!sym(color_by), group = !!sym(color_by))) + 
      #geom_path(stat= "summary", fun.data = "mean_se", linewidth=1)+# path line shows mean of y values
      stat_summary(aes(group = !!sym(color_by)), 
                   fun = median, 
                   geom = "line", 
                   linewidth = 1) +  # Use stat_summary to plot median path
      stat_summary(aes(fill=!!sym(color_by)), 
                   fun.min=function(z){quantile(z, 0.25)},
                   fun.max=function(z){quantile(z, 0.75)},
                   fun=median,
                   geom="ribbon",                              
                   alpha=0.2, color=NA)+#alpha is intensity of color
      labs(title = paste(value, " Consecutive Plot\n", phase, " phases")) +
      scale_fill_manual(values = manual_colors) +
      scale_color_manual(values = manual_colors) +
      scale_y_continuous(paste("Mean of ", value))+
      scale_x_discrete(name= consec_phase)+ 
      facet_grid(Sex~.)+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            title = element_text(size = 20),
            legend.key.size = unit(3, "lines"),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 20),
            axis.text = element_text(size = 20), 
            axis.title = element_text(size = 22), 
            strip.text = element_text(size = 20))
    
    #save plot in plot list
    if(j==1) plots[[2*i-1]] <- p
    else plots[[2*i]] <- p
    
    #save plot in file
    if(save_plots) ggsave(filename = paste0(working_directory, "/plots", "/consec_plots/", "consec-",value,"-",consec_phase, ".png"), plot = p, width = 12, height = 8)
    
    
    ## LME ##
    #create dynamic formula
    if(value=='CageEntropy'){
      
      formula <- as.formula(paste(value, "~ System * Sex + ", consec_phase, " + (1 | Batch)"))
    }else{
      
      formula <- as.formula(paste(value, "~ Group * Sex + ", consec_phase, " + (1 | AnimalID)"))
    }
    
    #do the lme
    message(paste("consec lme for", value,", ", phase))
    print(formula)
    model <- lmerTest::lmer(formula, data = consec_data)
    print(summary(model))
   
    #save model in statistics list
    #statistics[[i]] <- summary(model)
    
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
social_prox_CC <- tibble(CageChange=character(),
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
cc_tibbles <- list('SocialProx'= social_prox_CC, 'CoilCrossing'= coil_crossing_CC, 'MiceEntropy' = mice_entropy_CC, 'CageEntropy' = cage_entropy_CC)
values <- c('SocialProx', 'CoilCrossing', 'MiceEntropy', 'CageEntropy')

#for(i in seq_along(consec_tibbles)){
for(i in c(1,2,3,4)){
  consec_tibble <- consec_tibbles[[i]]
  #print(consec_tibble)
  value <- values[i]
  message(value)
  cc_tibble <- cc_tibbles[[value]]
  
  
  for(batch in batches){
    for(cageChange in cageChanges){
      message(paste(batch, cageChange))
     
      filtered_consec_tibble <- consec_tibble%>%
        filter(Batch==batch)%>%
        filter(CageChange==cageChange)
      #message("filtered_consec_tibble")
      #print(filtered_consec_tibble)
      
      
      #load csv
      mouse_which_system_tibble <-  tibble(read_delim(paste0(working_directory,"/tables/", batch, "_", cageChange, "_mouse_which_system_tibble.csv"),delim = ",", show_col_types = FALSE))
      
      
      for(system in colnames(mouse_which_system_tibble)){
        
        #message(system)
        
        mice_of_system <- mouse_which_system_tibble%>%  # this will create a character vector
          select(all_of(system))%>%
          pull()
        
        #if system is incomplete->skip system
        if(NA %in% mice_of_system){
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
          act_tibble <- filtered_consec_tibble%>%
            filter(Phase == "active")%>%
            filter(System==system)%>%
            select(all_of(value))
          
          inact_tibble <- filtered_consec_tibble%>%
            filter(Phase == "inactive")%>%
            filter(System==system)%>%
            select(all_of(value))
          
          #sum of multiple phases of one CageChange
          sum_act <- act_tibble%>%
            sum()
          #cat("sum of active: ", sum_act, "\n")
          
          sum_inact <- inact_tibble%>%
            sum()
          #cat("sum of inactive: ", sum_inact, "\n")
          
          avg_act <- act_tibble%>%
            unlist()%>%
            mean()
          #cat("mean of inactive: ", avg_act, "\n")
          
          avg_inact <- inact_tibble%>%
            unlist()%>%
            mean()
          #cat("mean of inactive: ", avg_inact, "\n")
          
          
          #add row to tibble
          cc_tibble <- cc_tibble%>% 
            add_row(CageChange=cageChange,
                    System=system,
                    Batch=batch,
                    Sex=sex,
                    Sum_act=sum_act,
                    Sum_inact=sum_inact,
                    Avg_act=avg_act,
                    Avg_inact=avg_inact,
                    SystemRank_act=NA,
                    SystemRank_inact=NA)
          
          #FOR RANK COMPUTING:
          #add total hours to vector
          act_system_vec <- append(act_system_vec, sum_act)
          inact_system_vec <- append(inact_system_vec, sum_inact)
       
          
        }else{    # 'SocialProx', 'CoilCrossing', 'MiceEntropy'
          for(mouse in mice_of_system){
            
            #print(mouse)
            
            ## create variables for new row of tibble ##
            group <- ifelse(mouse %in% sus_animals, "sus", ifelse(mouse %in% con_animals, "con", "res"))
            #print(group)
            sex <- ifelse(batch %in% c("B3", "B4", "B6"), "female", "male")
            
            #filter into act and inact phases of the individual and select all existing values
            act_tibble <- filtered_consec_tibble%>%
              filter(Phase == "active")%>%
              filter(AnimalID==mouse)%>%
              select(all_of(value))
            
            inact_tibble <- filtered_consec_tibble%>%
              filter(Phase == "inactive")%>%
              filter(AnimalID==mouse)%>%
              select(all_of(value))
            
            #sum of multiple phases of one CageChange
            sum_act <- act_tibble%>%
              sum()
            #cat("sum of active: ", sum_act, "\n")
            
            sum_inact <- inact_tibble%>%
              sum()
            #cat("sum of inactive: ", sum_inact, "\n")
            
            avg_act <- act_tibble%>%
              unlist()%>%
              mean()
            #cat("mean of inactive: ", avg_act, "\n")
            
            avg_inact <- inact_tibble%>%
              unlist()%>%
              mean()
            #cat("mean of inactive: ", avg_inact, "\n")
            
            
            
            #add row to tibble
            cc_tibble <- cc_tibble%>% 
              add_row(CageChange=cageChange,
                      System=system,
                      AnimalID=mouse,
                      Group=group,
                      Batch=batch,
                      Sex=sex,
                      Sum_act=sum_act,
                      Sum_inact=sum_inact,
                      Avg_act=avg_act,
                      Avg_inact=avg_inact,
                      SystemRank_act=NA,
                      SystemRank_inact=NA)
            
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
  color_by <- ifelse(value=='CageEntropy', 'System', 'Group')
  ifelse(value=='CageEntropy',manual_colors <- c("cornflowerblue", "yellow", "darkgreen", "orange", "purple"), manual_colors <- c("sus" = "tomato", "res" = "darkgreen", "con" = "deepskyblue4"))
  
  message(value)
  
  for(j in seq_along(phases)){
    cc_data <- cc_tibbles[[value]]
    #print(cc_data)
    phase <- phases[j]
    message(phase)
    
    
    ## PLOTS  ##
    
    ## AVG PLOT
    avg_phase <- ifelse(phase=="active", 'Avg_act', 'Avg_inact')
    
    #style ribbon
    avg_cc_p <- ggplot(data = cc_data, aes(x = CageChange, y = !!sym(avg_phase), color = !!sym(color_by), group=!!sym(color_by))) + 
      stat_summary(aes(group = !!sym(color_by)), fun = median, geom = "line", linewidth = 1) +  # plot median line
      stat_summary(aes(fill=!!sym(color_by)),                                                   # plot ribbon (median+-quantils)
                   fun.min=function(z){quantile(z, 0.25)},
                   fun.max=function(z){quantile(z, 0.75)},
                   fun=median,
                   geom="ribbon",                              
                   alpha=0.2, color=NA)+
      labs(title = paste(value, " CageChange Plot\n", phase, " phases")) +
      scale_fill_manual(values = manual_colors) +
      scale_color_manual(values = manual_colors) +
      scale_y_continuous(paste("Mean of ", value))+
      scale_x_discrete(name= "CageChange")+ 
      facet_grid(Sex~.)+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            title = element_text(size = 20),
            legend.key.size = unit(3, "lines"),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 20),
            axis.text = element_text(size = 20), 
            axis.title = element_text(size = 22), 
            strip.text = element_text(size = 20))
    
    
    #save plot in plot list
    if(j==1) avg_plots[[2*i-1]] <- avg_cc_p
    else avg_plots[[2*i]] <- avg_cc_p
    #save plot in file
    if(save_plots) ggsave(filename = paste0(working_directory, "/plots", "/cc_plots/", "avg-cc-",value,"-", phase, ".png"), plot = avg_cc_p, width = 12, height = 8)
    
    # LME-AVG 
    message(paste("cc lme for", value ,", ", avg_phase))
    
    #create dynamic formula
    if(value=='CageEntropy'){
      
      formula <- as.formula(paste(avg_phase, "~ System * Sex + CageChange + (1 | Batch)")) ####maybe * CageChange?
    }else{
      
      formula <- as.formula(paste(avg_phase, "~ Group * Sex + CageChange + (1 | AnimalID)"))
    }
    message(formula)
    
    model <- lmerTest::lmer(formula, data = cc_data)
    #print(summary(model))
    
    # save mixed model results to a file
    mixed_model_results_file <- paste0(working_directory, "/lme/cc-mixed_model_results/", paste0("avg-lme-",value,"-",phase, ".txt"))
    model_summary <- capture.output(summary(model))
    if(save_statistics) writeLines(model_summary, con = mixed_model_results_file)
    

    
    if(value!='CageEntropy'){
      
      
      
      ## AVG RANK PLOT
      
      #define the current rank-value based o the current phase
      rank_phase <- ifelse(phase=="active", 'SystemRank_act', 'SystemRank_inact')
      
      #group data to summarize all Cage changes for each individual
      avg_data <- cc_data%>%
        group_by(Batch, Sex, Group, AnimalID)%>%
        summarise(avg_rank = mean(!!sym(rank_phase)))#
      
      #exclude con group from visualisation if  socialProx (here not valuable)
      if(value=='SocialProx') avg_data <- avg_data%>% filter(Group!="con")
      
      
      avg_rank_p <- ggplot(data = avg_data, aes(x = Group, y = avg_rank, color = Group, shape = Batch)) +
      geom_jitter(size = 4, alpha = 0.7, width = 0.2, height = 0) +
        scale_shape_manual(values = c(3, 16, 17, 15, 5, 20)) + 
        scale_x_discrete("Animal ID") +
        scale_y_discrete("Rank", limits = c("1", "2", "3", "4")) +
        labs(title = paste("Mean cage change rank of each individual"), subtitle = paste(value, "-", phase))+ #, subtitle = paste(value,"sum from each",phase, " phase in a cageChange:computed in ranks", "\n", "rank4 = the highest sum", "\n", "rank1 = the lowest sum")
        scale_color_manual(values = c("sus" = "tomato", "res" = "darkgreen", "con" = "deepskyblue4")) +
        stat_summary(
          fun.min = function(z) {quantile(z, 0.25)},
          fun.max = function(z) {quantile(z, 0.75)},
          fun = median,
          color = "black",
          size = 0.8,
          shape = 16
        ) +
        facet_grid(Sex ~ .)+
        theme(title = element_text(size = 20),
              legend.key.size = unit(3, "lines"),
              legend.title = element_text(size = 18), 
              legend.text = element_text(size = 18),
              axis.text = element_text(size = 20), 
              axis.title = element_text(size = 22), 
              strip.text = element_text(size = 20),
              plot.title = element_text(hjust = 0.5),#centralised title
              plot.subtitle = element_text(hjust = 0.5)) 
      
      #save plot in plot list
      if(j==1) rank_plots[[2*i-1]] <- avg_rank_p
      else rank_plots[[2*i]] <- avg_rank_p
      #save plot in file
      if(save_plots) ggsave(filename = paste0(working_directory, "/plots", "/avg_rank_plots/", "avg_ranks-",value,"-", phase, ".png"), plot = avg_rank_p, width = 8, height = 8)
      
      
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
        
        
        
      }##end tobis statistic code
      
    }
    
    
  }
}
    


## PRINT PLOTS IN R ##
# Print all generated avg_plots
#separate
for (i in seq_along(avg_plots)) {
  print(avg_plots[[i]])
}
#together
gridExtra::grid.arrange(grobs = avg_plots, ncol = 2)
gridExtra::grid.arrange(grobs = rank_plots, ncol = 2)



# Create a grid of plotsfrom tobis code
gridExtra::grid.arrange(grobs = allPlots, ncol = 4)


###################### saving the results #############################################################

# Convert the list of test results to a data frame
allTestResultsDf <- bind_rows(allTestResults)
# Save the test results data frame to a CSV file
write.csv(allTestResultsDf, file = paste0(working_directory,"/statistics/test_results.csv"), row.names = FALSE)
# Save the post hoc results to a CSV file
if (!is.null(allPosthocResults) && length(allPosthocResults) > 0) {
  allPosthocResultsDf <- bind_rows(allPosthocResults)
  write.csv(allPosthocResultsDf, file = paste0(working_directory,"/statistics/posthoc_results.csv"), row.names = FALSE)
}


