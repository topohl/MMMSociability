## 03/2024
## Anja Magister
## ANALYSIS OF ANIMAL POSITIONS - COMPARING ##
##
##
## NEEDED FILE STRUCTURE IN WORKING DIRECTORY
##
##
## CUSTOMISABLE VARIABLES:
##
##

# libraries
library(readr)        # load readr package for reading csv files
library(ggplot2)      #for plots
library(tibble)       #important for tibble operations
library(dplyr)        #changes in tibbles
library(reshape2)     #for heatmap plot

# paths
working_directory <- "S:/Lab_Member/Anja/Git/MDC_Bachelor/E9_SIS_AnimalPos"
#working_directory <- "/home/anja/Dokumente/FU BERLIN/BA/Git/MDC_Bachelor/E9_SIS_AnimalPos"

#functions
source(paste0(working_directory,"/E9_SIS_AnimalPos-functions.R"))

# define batch and cage change
batch <- "B3"
cageChanges <- c("CC1", "CC2", "CC3", "CC4")

#sus and con animals
sus_animals <- readLines(paste0(working_directory,"/raw_data/sus_animals.csv"))
##csv path for con animals
con_animals <- readLines(paste0(working_directory,"/raw_data/con_animals.csv"))


##############################################################################################################

#take one table as information for unique mice
CC1_total_closeness_table <- tibble(read_delim(paste0(working_directory,"/tables/", batch, "_", cageChanges[1], "_total_closeness_table.csv"),delim = ",", show_col_types = FALSE))
# get the animal Ids of every mouse in the batch from one of the CC lists(no matter which) 
unique_mice <- unique(colnames(CC1_total_closeness_table)) 
unique_mice <- unique_mice[unique_mice != "Phase"]          #delete the phase column

#define sus res and con groups
sus_in_batch <- intersect(unique_mice, sus_animals)
con_in_batch <- intersect(unique_mice, con_animals)
res_in_batch <- setdiff(setdiff(unique_mice,sus_in_batch), con_in_batch)

#define sex
sex <- ifelse(batch %in% c("B3", "B4", "B6"), "female", "male")

##########################################################################################################
#create rank tibble

#initialize rank tibble
rank_tibble <- tibble(CageChange=NA,
                      System=NA,
                      Animal_ID=NA,
                      Group=NA,
                      Sex=NA,
                      Social_h_act=NA,
                      Social_h_inact=NA,
                      Social_h_total=NA,
                      Avg_social_h_act=NA,
                      Avg_social_h_inact=NA,
                      SystemRank_act=NA,
                      SystemRank_inact=NA,
                      SystemRank_total=NA)

##fill rank_tibble of this batch##
for (cageChange in cageChanges) {
  #read csv
  mouse_which_system_tibble <-  tibble(read_delim(paste0(working_directory,"/tables/", batch, "_", cageChange, "_mouse_which_system_tibble.csv"),delim = ",", show_col_types = FALSE))
  total_closeness_table <- tibble(read_delim(paste0(working_directory,"/tables/", batch, "_", cageChange, "_total_closeness_table.csv"),delim = ",", show_col_types = FALSE))
  
  #special treatment for b3 1545...
  if (batch=="B3"&&cageChange=="CC3") {
    mouse_which_system_tibble[4,2] <- "1545"
    
    total_closeness_table <- total_closeness_table%>%
      rename('1545' = OR1545)
  }

  for (system in colnames(mouse_which_system_tibble)) {
    print(system)

    mice_of_system <- mouse_which_system_tibble%>%  # this will create a character vector
      select(all_of(system)) %>%
      pull()
    
    #if system is incomplete->skip system
    if (NA %in% mice_of_system) {
      cat("skip ", system, " in ", cageChange, "\n")
      #skip for loop
      next
    }
    
    #empty vector for determining rank:
    #(max social hours,active phase, of a mouse of a system)
    social_hours_act_system_vec <- c()
    #(max social hours,inactive phase, of a mouse of a system)
    social_hours_inact_system_vec <- c()
    #(max social hours of a mouse of a system)
    social_hours_total_system_vec <- c()
    
    for (mouse in mice_of_system) {
      
      ## create variables for new row of tibble ##
      group <- ifelse(mouse %in% res_in_batch, "res", ifelse(mouse %in% sus_in_batch, "sus", "con"))
      
      #sum of multiple phases of one CageChange
      social_h_act <- total_closeness_table%>%
        filter(Phase %in% c("A1", "A2", "A3", "A4")) %>%
        select(all_of(mouse)) %>%#select column with mouse id
        sum()
      #cat("sum of active: ", social_h_act, "\n")
      
      social_h_inact <- total_closeness_table%>%
        filter(Phase %in% c("I2", "I3", "I4")) %>%
        select(all_of(mouse)) %>%
        sum()
      #cat("sum of inactive: ", social_h_inact, "\n")
      
      social_h_total <- sum(social_h_act, social_h_inact)
      
      #average of a phase of one CageChange
      avg_social_h_act <- total_closeness_table%>%
        filter(Phase %in% c("A1", "A2", "A3", "A4")) %>%
        select(all_of(mouse)) %>%
        unlist() %>%#change tibble into vector
        mean()
      #cat("mean of active: ", avg_social_h_act, "\n")
      
      avg_social_h_inact <- total_closeness_table%>%
        filter(Phase %in% c("I2", "I3", "I4")) %>%
        select(all_of(mouse)) %>%
        unlist() %>%
        mean()
      #cat("mean of inactive: ", avg_social_h_inact, "\n")
      
      #add row to tibble
      rank_tibble <- rank_tibble %>% 
        add_row(CageChange = cageChange,
          System = system,
          Animal_ID = mouse,
          Group = group,
          Sex = sex,
          Social_h_act = social_h_act,
          Social_h_inact = social_h_inact,
          Social_h_total = social_h_total,
          Avg_social_h_act = avg_social_h_act,
          Avg_social_h_inact = avg_social_h_inact,
          SystemRank_act = NA,
          SystemRank_inact = NA,
          SystemRank_total = NA)
      
      #FOR RANK COMPUTING:
      #add total hours to vector
      social_hours_act_system_vec <- append(social_hours_act_system_vec, social_h_act)
      social_hours_inact_system_vec <- append(social_hours_inact_system_vec, social_h_inact)
      social_hours_total_system_vec <- append(social_hours_total_system_vec, social_h_total)
    }
    
    #RANK COMPUTING:
    rank_tibble <- compute_rank(rank_tibble, social_hours_act_system_vec, system, 'Social_h_act', 'SystemRank_act')
    rank_tibble <- compute_rank(rank_tibble, social_hours_inact_system_vec, system, 'Social_h_inact', 'SystemRank_inact')
    rank_tibble <- compute_rank(rank_tibble, social_hours_total_system_vec, system, 'Social_h_total', 'SystemRank_total')
  }
}
# delete first row in tibble (NA values)
rank_tibble <- na.omit(rank_tibble)

### plots for rank tibble ####################### 

#change counted seconds to hours
rank_tibble <- rank_tibble%>%
  mutate(Social_h_total = as.integer(Social_h_total)) %>%  ##sum
  mutate(Social_h_total = ifelse(Social_h_total!=0,Social_h_total/3600,Social_h_total)) %>%
  mutate(Social_h_act = as.integer(Social_h_act)) %>%
  mutate(Social_h_act = ifelse(Social_h_act!=0,Social_h_act/3600,Social_h_act)) %>%
  mutate(Social_h_inact = as.integer(Social_h_inact)) %>%
  mutate(Social_h_inact = ifelse(Social_h_inact!=0,Social_h_inact/3600,Social_h_inact)) %>%
  mutate(Avg_social_h_act = as.integer(Avg_social_h_act)) %>%  ##average
  mutate(Avg_social_h_act = ifelse(Avg_social_h_act!=0,Avg_social_h_act/3600,Avg_social_h_act)) %>%
  mutate(Avg_social_h_inact = as.integer(Avg_social_h_inact)) %>%
  mutate(Avg_social_h_inact = ifelse(Avg_social_h_inact!=0,Avg_social_h_inact/3600,Avg_social_h_inact))

#######
#generate plots
total_hours_plot <-  ggplot(data=rank_tibble, aes(x=Animal_ID, y=Social_h_total, color=Group))+
  geom_point(na.rm=TRUE)+
  scale_color_manual(values= c("sus"= "red", "res" = "green", "con" = "deepskyblue4"), name = paste("total hours, ", batch))+
  scale_y_continuous("total close contact in h")+
  scale_x_discrete(labels = NULL, breaks = NULL)+
  facet_grid(~CageChange)

rank_plot <- ggplot(data=rank_tibble, aes(x=factor(SystemRank_total), y=factor(Animal_ID, levels = unique(Animal_ID[order(Group)])), color=Group))+
  geom_point(na.rm=TRUE)+
  scale_color_manual(values= c("sus"= "red", "res" = "green", "con" = "deepskyblue4"), name = paste("rank overview, ", batch))+
  scale_y_discrete("Animal ID")+
  scale_x_discrete("rank", limits = c("1","2","3","4"), breaks = c("1","2","3","4")) +
  facet_grid(~CageChange) +
  facet_grid(~CageChange)+
  theme(strip.placement = "bottom")

## total_hours_plot with divided phases

##first attempt##
total_hours_phase2_plot <-  ggplot(data=rank_tibble)+
  geom_point(aes(x=Animal_ID, y=Social_h_act, color="active"), na.rm=TRUE)+
  geom_point(aes(x=Animal_ID, y=Social_h_inact, color="inactive"), na.rm=TRUE)+
  #scale_color_manual(values= c("sus"= "red", "res" = "green", "con" = "deepskyblue4"), name = paste("total hours, ", batch))+
  scale_y_continuous("total close contact in h")+
  scale_x_discrete(labels = NULL, breaks = NULL)+
  facet_grid(~CageChange)

##second attempt##
#select needed columns
phase_rank_tibble <- select(rank_tibble, c(CageChange,Animal_ID, Group, Social_h_act, Social_h_inact))
melt_rank_tibble <- melt(phase_rank_tibble, id=c('CageChange', 'Animal_ID', 'Group'))
names(melt_rank_tibble) <- c('CageChange', 'Animal_ID', 'Group', 'Phase', 'Phase_total_hours')

total_hours_phase_plot <-  ggplot(data = melt_rank_tibble, aes(x = Phase, y = Phase_total_hours, color=Group))+
  geom_jitter(na.rm=TRUE)+
  scale_color_manual(values= c("sus"= "red", "res" = "green", "con" = "deepskyblue4"), name = paste("total hours, ", batch)) +
  scale_y_continuous("total close contact in h") +
  scale_x_discrete() +
  facet_grid(~CageChange)

## third attempt ##
## boxplot 
total_hours_phase_boxplot <-  ggplot(data=melt_rank_tibble, aes(x = Phase, y = Phase_total_hours, fill = Group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("sus"= "tomato", "res" = "lightgreen", "con" = "deepskyblue4"), name = paste("total social hours, ", batch)) +
  scale_y_continuous("total close contact in h") +
  theme_bw() +
  facet_grid(~CageChange)#for general plot just leave out this line

########################################################################################################################
## social score ##
social_score_tibble <- tibble(
  Animal_ID = unique_mice,
  Group = NA,
  Social_score_act = NA,
  Social_score_inact = NA,
  Social_score_total = NA
)

for (mouse in unique_mice) {
  #define group of animal
  group <- ifelse(mouse %in% res_in_batch, "res", ifelse(mouse %in% sus_in_batch, "sus", "con"))
  
  #pull ranks of the mouse from every CC
  rank_act_vec <- rank_tibble %>%
    filter(Animal_ID == mouse) %>%
    select(SystemRank_act) %>%
    pull()
    #sum() 
  
  rank_inact_vec <- rank_tibble %>%
    filter(Animal_ID == mouse) %>%
    select(SystemRank_inact) %>%
    pull()
    #sum() 
  
  rank_total_vec <- rank_tibble %>%
    filter(Animal_ID == mouse) %>%
    select(SystemRank_total) %>%
    pull() 
    #sum() 
  
  #compute scores (number of ranks in all 4 CCs together)
  score_act <- translate_rank_in_score(rank_act_vec)
  score_inact <- translate_rank_in_score(rank_inact_vec) 
  score_total <- translate_rank_in_score(rank_total_vec) 
  
  social_score_tibble <- social_score_tibble %>%
    mutate(Group = ifelse(Animal_ID == mouse, group, Group)) %>%
    mutate(Social_score_act = ifelse(Animal_ID == mouse, score_act, Social_score_act)) %>%
    mutate(Social_score_inact = ifelse(Animal_ID == mouse, score_inact, Social_score_inact)) %>%
    mutate(Social_score_total = ifelse(Animal_ID == mouse, score_total, Social_score_total))
}

###################################################
## plots for social score

#define labels with sorted "Group" and sorted ids in vector
#label_x_axis_animalIDs <- social_score_tibble%>%
#  arrange(., Group) %>%
#  pull(Animal_ID)


social_score_plot <- ggplot(data=social_score_tibble, aes(x=Animal_ID, y=Social_score_total, color=Group))+
  geom_point()+
  scale_color_manual(values= c("sus"= "red", "res" = "green", "con" = "deepskyblue4"), name = paste("social score, ", batch))+
  scale_y_continuous("social score", breaks = seq(1, 20, by = 1), labels = seq(1, 20, by = 1))+
  scale_x_discrete("Animal ID")


###
#same plot with mean added as a line
# calculate grouped means

###total###
#exclude con group from plots(not interesting)
social_score_tibble <- social_score_tibble%>%
  filter(Group != "con")

group_means <- social_score_tibble %>%
  group_by(Group) %>%
  summarise(mean_social_score = mean(Social_score_total))

social_score_plot_with_mean <- ggplot()+
  geom_point(data=social_score_tibble, aes(x=Group, y=Social_score_total, color=Group)) +
  geom_line(data = group_means, aes(x = Group, y = mean_social_score, group = 1), linetype = "dashed") +
  #geom_segment(data = group_means, aes(x = 0, xend = Group, y = mean_social_score, yend = mean_social_score), linetype = "dashed") +
  scale_color_manual(values= c("sus"= "red", "res" = "green", "con" = "deepskyblue4"), name = paste("total social score, ", batch)) +
  scale_y_continuous("social score", breaks = seq(1, 20, by = 1), labels = seq(1, 20, by = 1)) +
  scale_x_discrete("Animal ID")

###active###
group_means_act <- social_score_tibble %>%
  group_by(Group) %>%
  summarise(mean_social_score = mean(Social_score_act))

act_social_score_plot_with_mean <- ggplot()+
  geom_point(data=social_score_tibble, aes(x=Group, y=Social_score_act, color=Group))+
  geom_line(data = group_means_act, aes(x = Group, y = mean_social_score, group = 1), linetype = "dashed") +
  #geom_segment(data = group_means, aes(x = 0, xend = Group, y = mean_social_score, yend = mean_social_score), linetype = "dashed") +
  scale_color_manual(values= c("sus"= "red", "res" = "green", "con" = "deepskyblue4"), name = paste("active social score, ", batch))+
  scale_y_continuous("social score", breaks = seq(1, 20, by = 1), labels = seq(1, 20, by = 1))+
  scale_x_discrete("Animal ID")

###inactive###
group_means_inact <- social_score_tibble %>%
  group_by(Group) %>%
  summarise(mean_social_score = mean(Social_score_inact))

inact_social_score_plot_with_mean <- ggplot()+
  geom_point(data=social_score_tibble, aes(x=Group, y=Social_score_inact, color=Group))+
  geom_line(data = group_means_inact, aes(x = Group, y = mean_social_score, group = 1), linetype = "dashed") +
  #geom_segment(data = group_means, aes(x = 0, xend = Group, y = mean_social_score, yend = mean_social_score), linetype = "dashed") +
  scale_color_manual(values= c("sus"= "red", "res" = "green", "con" = "deepskyblue4"), name = paste("inact social score, ", batch))+
  scale_y_continuous("social score", breaks = seq(1, 20, by = 1), labels = seq(1, 20, by = 1))+
  scale_x_discrete("Animal ID")

### barplot
total_social_score_barplot <- ggplot(data=social_score_tibble, aes(x = factor(Animal_ID, levels = unique(Animal_ID[order(Group)])), y = Social_score_total, fill = Group))+
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_discrete("Animal ID")

act_social_score_barplot <- ggplot(data=social_score_tibble, aes(x = factor(Animal_ID, levels = unique(Animal_ID[order(Group)])), y = Social_score_act, fill = Group))+
  geom_bar(stat = "identity")+
  scale_x_discrete("Animal ID")

inact_social_score_barplot <- ggplot(data = social_score_tibble, aes(x = factor(Animal_ID, levels = unique(Animal_ID[order(Group)])), y = Social_score_inact, fill = Group))+
  geom_bar(stat = "identity")+
  scale_x_discrete("Animal ID")

##boxplot
par(mar = c(0, 0, 0, 0))                  # Remove space around plot
par(bg = "#353436")  

boxplot(data=social_score_tibble,
        Social_score_total ~ Group,         
        xlab = "Group", ylab = "Social_score_total",
        col = "yellow",
        border = "yellow",
        pch = 16)
points(social_score_tibble$Group, social_score_tibble$Social_score_total,                       # Sophisticated overlay of jittered X variable 
       col = "#1b98e0",
       pch = 16,
       cex = 0.4)


###########################################################################################################################
#do a copy of the rank tibble and create another column for score of every mouse per cagechange
#just turns the rank-column into the social score
#
rank_tibble_with_SystemRankScore <- rank_tibble%>%
  mutate(SystemRankScore_total = ifelse(SystemRank_total == 1,4,ifelse(SystemRank_total == 2, 3, ifelse(SystemRank_total == 3, 2, ifelse(SystemRank_total == 4, 1, NA)))))

#rank plot-mean--WITH FACET FOR CCs
###total###
group_means <- rank_tibble_with_SystemRankScore %>%
  group_by(Group, CageChange) %>%
  #group_by(CageChange) %>%
  summarise(mean_social_score = mean(SystemRankScore_total))

social_score_plot_with_mean <- ggplot()+
  geom_point(data = rank_tibble_with_SystemRankScore, aes(x = Group, y= SystemRankScore_total, color = Group), position = position_jitter(w = 0.3, h = 0)) +
  geom_line(data = group_means, aes(x = Group, y = mean_social_score, group = 1), linetype = "dashed") +
  #geom_segment(data = group_means, aes(x = 0, xend = Group, y = mean_social_score, yend = mean_social_score), linetype = "dashed") +
  scale_color_manual(values= c("sus"= "red", "res" = "green", "con" = "deepskyblue4"), name = paste("total social score, ", batch)) +
  scale_y_continuous("social score", breaks = seq(1, 20, by = 1), labels = seq(1, 20, by = 1))+
  scale_x_discrete("Animal ID") +
  facet_grid(~CageChange)

###########################################################################################################################
#total_movement_tibble
#entropy

#generate tibble with everything together
all_CC_total_movement_tibble <- tibble()
for (cageChange in cageChanges) {
  print(cageChange)
  #read csv
  mouse_which_system_tibble <-  tibble(read_delim(paste0(working_directory, "/tables/", batch, "_", cageChange, "_mouse_which_system_tibble.csv"), delim = ",", show_col_types = FALSE))
  total_movement_tibble <- tibble(read_delim(paste0(working_directory, "/tables/", batch, "_", cageChange, "_total_movement_table.csv"), delim = ",", show_col_types = FALSE))
  
  ##alterate tibble
  total_movement_tibble <- total_movement_tibble%>%
    mutate(CageChange=cageChange) %>% #add column CageChange
    melt(id = c('CageChange', 'Phase'))
  
  #rename
  names(total_movement_tibble) <- c('CageChange', 'Phase', 'AnimalID', 'movement')
  
  if (length(all_CC_total_movement_tibble) == 0) {
    all_CC_total_movement_tibble <- total_movement_tibble
  } else {
    all_CC_total_movement_tibble <- bind_rows(all_CC_total_movement_tibble, total_movement_tibble)
  } 
}


#########PLOTS FOR MOVEMENT#########

###plot for every individual
#filter tibble without systems
all_CC_total_movement_tibble_individuals <- all_CC_total_movement_tibble%>%
  filter(!grepl("^sys\\.[1-9]$", AnimalID))

movement_plot_individual <- ggplot(data=all_CC_total_movement_tibble_individuals, aes(x=Phase, y=movement, color=AnimalID, group=AnimalID))+
  geom_point()+
  geom_line()+
  scale_y_continuous("number of movements")+
  scale_x_discrete(limits = c("I1", "A1", "I2", "A2", "I3", "A3", "I4", "A4", "I5"))+
  facet_grid(~CageChange)
  

###plot for every system
#filter tibble with systems
all_CC_total_movement_tibble_systems <- all_CC_total_movement_tibble%>%
  filter(grepl("^sys\\.[1-9]$", AnimalID))
#rename
names(all_CC_total_movement_tibble_systems) <- c('CageChange', 'Phase', 'System', 'movement')

movement_plot_systems <- ggplot(data=all_CC_total_movement_tibble_systems, aes(x=Phase, y=movement, color=System, group=System))+
  geom_point()+
  geom_line()+
  scale_y_continuous("number of movements")+
  scale_x_discrete(limits = c("I1", "A1", "I2", "A2", "I3", "A3", "I4", "A4", "I5"))+
  facet_grid(~CageChange)

###plot for every group with mean (divided in active/inactive phases)
#mutate group
all_CC_total_movement_tibble_groups <- all_CC_total_movement_tibble_individuals %>%
  mutate(Group = ifelse(AnimalID %in% res_in_batch, "res", ifelse(AnimalID %in% sus_in_batch, "sus", "con")))


#divide in active / inactive
all_CC_total_movement_tibble_act_phases <- all_CC_total_movement_tibble_groups%>%
  filter(grepl("^A\\d+", Phase))

all_CC_total_movement_tibble_inact_phases <- all_CC_total_movement_tibble_groups%>%
  filter(grepl("^I\\d+", Phase))
  

##active plot
group_means <- all_CC_total_movement_tibble_act_phases %>%
  group_by(Group, CageChange, Phase) %>%
  #group_by(CageChange) %>%
  summarise(mean_movement = mean(movement)) 

group_means_res <- group_means%>%
  filter(Group == "res")
group_means_sus <- group_means%>%
  filter(Group == "sus")
group_means_con <- group_means%>%
  filter(Group == "con")
  
movement_plot_groups <- ggplot()+
  geom_point(data=all_CC_total_movement_tibble_act_phases, aes(x=Phase, y=movement, color=Group, group=Group))+
  geom_line(data = group_means_res, aes(x = Phase, y = mean_movement, group = 1), color="lightgreen") +
  geom_line(data = group_means_sus, aes(x = Phase, y = mean_movement, group = 1), color="tomato") +
  geom_line(data = group_means_con, aes(x = Phase, y = mean_movement, group = 1), color="deepskyblue4") +
  scale_color_manual(values= c("sus"= "tomato", "res" = "lightgreen", "con" = "deepskyblue4"), name = paste("act. phase, number of movements ", batch))+
  scale_y_continuous()+
  scale_x_discrete(limits = c("I1", "A1", "I2", "A2", "I3", "A3", "I4", "A4", "I5"))+
  facet_grid(~CageChange)

##inactive plot
group_means <- all_CC_total_movement_tibble_inact_phases %>%
  group_by(Group, CageChange, Phase) %>%
  #group_by(CageChange) %>%
  summarise(mean_movement = mean(movement)) 

group_means_res <- group_means%>%
  filter(Group == "res")
group_means_sus <- group_means%>%
  filter(Group == "sus")
group_means_con <- group_means%>%
  filter(Group == "con")

movement_plot_groups <- ggplot()+
  geom_point(data=all_CC_total_movement_tibble_inact_phases, aes(x=Phase, y=movement, color=Group, group=Group))+
  geom_line(data = group_means_res, aes(x = Phase, y = mean_movement, group = 1), color="lightgreen") +
  geom_line(data = group_means_sus, aes(x = Phase, y = mean_movement, group = 1), color="tomato") +
  geom_line(data = group_means_con, aes(x = Phase, y = mean_movement, group = 1), color="deepskyblue4") +
  scale_color_manual(values= c("sus"= "tomato", "res" = "lightgreen", "con" = "deepskyblue4"), name = paste("inact phase, number of movements ", batch))+
  scale_y_continuous()+
  scale_x_discrete(limits = c("I1", "A1", "I2", "A2", "I3", "A3", "I4", "A4", "I5"))+
  facet_grid(~CageChange)
