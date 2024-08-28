## 05/2024
## Anja Magister
## ANALYSIS OF ANIMAL POSITIONS - COMPARING2 ##
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

batches <- c("B1", "B2", "B3", "B4", "B5")#B6
cageChanges <- c("CC1", "CC2", "CC3", "CC4")

#sus and con animals
sus_animals <- readLines(paste0(working_directory,"/raw_data/sus_animals.csv"))
##csv path for con animals
con_animals <- readLines(paste0(working_directory,"/raw_data/con_animals.csv"))


##########################################################################################################
#create rank tibble



#initialize rank tibble
rank_tibble <- tibble(CageChange=NA,
                      System=NA,
                      Animal_ID=NA,
                      Group=NA,
                      Batch=NA,
                      Sex=NA,
                      Social_h_act=NA,
                      Social_h_inact=NA,
                      Social_h_total=NA,
                      Avg_social_h_act=NA,
                      Avg_social_h_inact=NA,
                      SystemRank_act=NA,
                      SystemRank_inact=NA,
                      SystemRank_total=NA)

for(batch in batches){
  for(cageChange in cageChanges){
    #read csv
    mouse_which_system_tibble <-  tibble(read_delim(paste0(working_directory,"/tables/", batch, "_", cageChange, "_mouse_which_system_tibble.csv"),delim = ",", show_col_types = FALSE))
    total_closeness_table <- tibble(read_delim(paste0(working_directory,"/tables/", batch, "_", cageChange, "_total_closeness_table.csv"),delim = ",", show_col_types = FALSE))
    
    
    #special treatment for b3 1545...
    #if(batch=="B3"&&cageChange=="CC3"){
    #  mouse_which_system_tibble[4,2] <- "1545"
      
   #   total_closeness_table <- total_closeness_table%>%
   #     rename('1545' = OR1545)
   # }
    
    #special treatment for every CC4 after A2(grid in cage)-> not usable for social analysis
    #use only A1,I2,A2 for the computation of the rank from CC4
    if(cageChange=="CC4"){
      
      total_closeness_table <- total_closeness_table%>%
        filter(Phase %in% c("A1","I2","A2"))
    }
    
    for(system in colnames(mouse_which_system_tibble)){
      
      print(system)
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
      #(max social hours,active phase, of a mouse of a system)
      social_hours_act_system_vec <- c()
      #(max social hours,inactive phase, of a mouse of a system)
      social_hours_inact_system_vec <- c()
      #(max social hours of a mouse of a system)
      social_hours_total_system_vec <- c()
      
      for(mouse in mice_of_system){
        
        #print(mouse)
        ## create variables for new row of tibble ##
        group <- ifelse(mouse %in% sus_animals, "sus", ifelse(mouse %in% con_animals, "con", "res"))
        #print(group)
        sex <- ifelse(batch %in% c("B3", "B4", "B6"), "female", "male")
        
        #sum of multiple phases of one CageChange
        social_h_act <- total_closeness_table%>%
          filter(Phase %in% c("A1", "A2", "A3", "A4"))%>%
          select(all_of(mouse))%>%#select column with mouse id
          sum()
        cat("sum of active: ", social_h_act, "\n")
        
        social_h_inact <- total_closeness_table%>%
          filter(Phase %in% c("I2", "I3", "I4"))%>%
          select(all_of(mouse))%>%
          sum()
        #cat("sum of inactive: ", social_h_inact, "\n")
        
        social_h_total <- sum(social_h_act, social_h_inact)
        #average of a phase of one CageChange
        avg_social_h_act <- total_closeness_table%>%
          filter(Phase %in% c("A1", "A2", "A3", "A4"))%>%
          select(all_of(mouse))%>%
          unlist()%>%#change tibble into vector
          mean()
        #cat("mean of active: ", avg_social_h_act, "\n")
        
        avg_social_h_inact <- total_closeness_table%>%
          filter(Phase %in% c("I2", "I3", "I4"))%>%
          select(all_of(mouse))%>%
          unlist()%>%
          mean()
        #cat("mean of inactive: ", avg_social_h_inact, "\n")
        
        
        #add row to tibble
        rank_tibble <- rank_tibble%>% 
          add_row(CageChange=cageChange,
                  System=system,
                  Animal_ID=mouse,
                  Group=group,
                  Batch=batch,
                  Sex=sex,
                  Social_h_act=social_h_act,
                  Social_h_inact=social_h_inact,
                  Social_h_total=social_h_total,
                  Avg_social_h_act=avg_social_h_act,
                  Avg_social_h_inact=avg_social_h_inact,
                  SystemRank_act=NA,
                  SystemRank_inact=NA,
                  SystemRank_total=NA)
        
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
  
}

# delete first row in tibble (NA values)
rank_tibble <- na.omit(rank_tibble)

#change counted seconds to hours??
rank_tibble <- rank_tibble%>%
  mutate(Social_h_total = as.integer(Social_h_total))%>%  ##sum
  mutate(Social_h_total = ifelse(Social_h_total!=0,Social_h_total/3600,Social_h_total))%>%
  mutate(Social_h_act = as.integer(Social_h_act))%>%
  mutate(Social_h_act = ifelse(Social_h_act!=0,Social_h_act/3600,Social_h_act))%>%
  mutate(Social_h_inact = as.integer(Social_h_inact))%>%
  mutate(Social_h_inact = ifelse(Social_h_inact!=0,Social_h_inact/3600,Social_h_inact))%>%
  mutate(Avg_social_h_act = as.integer(Avg_social_h_act))%>%  ##average
  mutate(Avg_social_h_act = ifelse(Avg_social_h_act!=0,Avg_social_h_act/3600,Avg_social_h_act))%>%
  mutate(Avg_social_h_inact = as.integer(Avg_social_h_inact))%>%
  mutate(Avg_social_h_inact = ifelse(Avg_social_h_inact!=0,Avg_social_h_inact/3600,Avg_social_h_inact))



################## PLOTS ########################################

################# SOCIAL PROX IN HOURS #################
## ACTIVE VS INACTIVE PHASES ##
#select needed columns
phase_rank_tibble <- select(rank_tibble, c(CageChange,Animal_ID, Group, Sex, Social_h_act, Social_h_inact))
names(phase_rank_tibble) <- c('CageChange', 'Animal_ID', 'Group', 'Sex', 'active', 'inactive')

melt_rank_tibble <- melt(phase_rank_tibble, id=c('CageChange', 'Animal_ID', 'Group', 'Sex'))
names(melt_rank_tibble) <- c('CageChange', 'Animal_ID', 'Group', 'Sex', 'Phase', 'Hours')

## scatterplot 
total_hours_phase_scatterplot <-  ggplot(data=melt_rank_tibble, aes(x=Group, y=Hours, color=Phase))+
  geom_jitter(aes(fill=Phase), size=4, alpha=0.7, width=0.2, shape=16)+
  #scale_color_manual(values= c("sus"= "tomato", "res" = "lightgreen", "con" = "deepskyblue4"), name = paste("total social hours, ", batch))+
  scale_y_continuous("social proximity in h")+
  stat_summary(
    fun.min=function(z){quantile(z, 0.25)},
    fun.max=function(z){quantile(z, 0.75)},
    fun=median,
    color="black",
    size=0.8,
    shape=16,
    width=1)+
  theme_bw()+
  facet_grid(Phase~.)#for general plot just leave out this line

## boxplot 
total_hours_phase_boxplot <-  ggplot(data=melt_rank_tibble, aes(x=Phase, y=Hours, fill=Group))+
  geom_boxplot()+
  scale_fill_manual(values= c("sus"= "tomato", "res" = "lightgreen", "con" = "deepskyblue4"), name = "group")+
  scale_y_continuous("social proximity in h")+
  labs(title = paste("Social measurement in sum hours Batch 1-5"))+
  theme_bw()+
  facet_grid(Sex~CageChange)+
  theme(title = element_text(size = 20),
        legend.key.size = unit(5, "lines"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 20), 
        axis.title = element_text(size = 22), 
        strip.text = element_text(size = 20))

################# SOCIAL PROX IN AVG HOURS #################
## ACTIVE VS INACTIVE PHASES ##
#select needed columns
phase_rank_tibble <- select(rank_tibble, c(CageChange,Animal_ID, Group, Sex, Avg_social_h_act, Avg_social_h_inact))
names(phase_rank_tibble) <- c('CageChange', 'Animal_ID', 'Group', 'Sex', 'active', 'inactive')

melt_rank_tibble <- melt(phase_rank_tibble, id=c('CageChange', 'Animal_ID', 'Group', 'Sex'))
names(melt_rank_tibble) <- c('CageChange', 'Animal_ID', 'Group', 'Sex', 'Phase', 'Hours')

## boxplot 
total_hours_phase_boxplot <-  ggplot(data=melt_rank_tibble, aes(x=Phase, y=Hours, fill=Group))+
  geom_boxplot()+
  scale_fill_manual(values= c("sus"= "tomato", "res" = "lightgreen", "con" = "deepskyblue4"), name = "group")+
  scale_y_continuous("social proximity in avg h")+
  labs(title = paste("Social measurement in avg hours Batch 1-5"))+
  theme_bw()+
  facet_grid(Sex~CageChange)+
  theme(title = element_text(size = 20),
        legend.key.size = unit(5, "lines"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 20), 
        axis.title = element_text(size = 22), 
        strip.text = element_text(size = 20))

################# SOCIAL PROX IN RANKS #################
## ACTIVE VS INACTIVE PHASES ##
#select needed columns
phase_rank_tibble <- select(rank_tibble, c(CageChange,Animal_ID, Group, Sex, SystemRank_act, SystemRank_inact))
names(phase_rank_tibble) <- c('CageChange', 'Animal_ID', 'Group', 'Sex', 'active', 'inactive')
melt_rank_tibble <- melt(phase_rank_tibble, id=c('CageChange', 'Animal_ID', 'Group', 'Sex'))
names(melt_rank_tibble) <- c('CageChange', 'Animal_ID', 'Group', 'Sex', 'Phase', 'Rank')

#delete the rows with con group for this plot. not necessary
plot_rank_tibble <- melt_rank_tibble%>%
  filter(Group != "con")

##sex difference!!
rank_phase_scatterplot <-  ggplot(data=plot_rank_tibble, aes(x=Group, y=Rank, color=Phase))+
  geom_jitter(aes(fill=Phase), size=4, alpha=0.7, width=0.2, shape=16)+
  scale_y_continuous("Rank(1 = most social)")+
  labs(title = paste("Social measurement in Ranks"))+
  stat_summary(
    fun.min=function(z){quantile(z, 0.25)},
    fun.max=function(z){quantile(z, 0.75)},
    fun=median,
    color="black",
    size=0.8,
    shape=16)+
  theme_bw()+
  theme(legend.key.size = unit(5, "lines"), # Legendenliniengröße einstellen
        legend.title = element_text(size = 20), # Legendenüberschriftsgröße einstellen
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 20), # Achsenbeschriftungsgröße einstellen
        axis.title = element_text(size = 22), # Achsentitelgröße einstellen
        strip.text = element_text(size = 20))+ 
  facet_grid(Sex~CageChange)

#general between phases all together
rank_phase_scatterplot <-  ggplot(data=plot_rank_tibble, aes(x=Group, y=Phase_ranks, color=Phase))+
  geom_jitter(aes(fill=Phase), size=4, alpha=0.7, width=0.2, shape=16)+
  scale_y_continuous("social measurement in Rank(1 = most social)")+
  stat_summary(
    fun.min=function(z){quantile(z, 0.25)},
    fun.max=function(z){quantile(z, 0.75)},
    fun=median,
    color="black",
    size=0.8,
    shape=16)+
  theme_bw()

#grid between phases(not that important)
rank_phase_scatterplot <-  ggplot(data=plot_rank_tibble, aes(x=Group, y=Phase_ranks, color=Phase))+
  geom_jitter(aes(fill=Phase), size=4, alpha=0.7, width=0.2, shape=16)+
  scale_y_continuous("social measurement in Rank(1 = most social)")+
  stat_summary(
    fun.min=function(z){quantile(z, 0.25)},
    fun.max=function(z){quantile(z, 0.75)},
    fun=median,
    color="black",
    size=0.8,
    shape=16)+
  theme_bw()+
  facet_grid(Phase~.)


#number of sex entrys:
female <- rank_tibble%>%
  filter(Sex=="female")
#=144
male <- rank_tibble%>%
  filter(Sex=="male")
#=104

########################################################################################################################
## social score ##
#tibble for all batches together on basis of the rank tibble

unique_mice <- unique(rank_tibble$Animal_ID)

social_score_tibble <- tibble(Animal_ID=unique_mice,
                              Group=NA,
                              Batch=NA,
                              Sex=NA,
                              Social_score_act=NA,
                              Social_score_inact=NA,
                              Social_score_total=NA,
                              Std_dev_act=NA,
                              Std_dev_inact=NA,
                              Std_dev_total=NA) 


for(mouse in unique_mice){
  
  #define group of animal
  group <- ifelse(mouse %in% sus_animals, "sus", ifelse(mouse %in% con_animals, "con", "res"))
  
  #define batch of animal
  batch <- rank_tibble%>%
    filter(Animal_ID==mouse)%>%
    head(1)%>%
    select(Batch)%>%
    pull()
  #define sex of animal
  sex <- rank_tibble%>%
    filter(Animal_ID==mouse)%>%
    head(1)%>%
    select(Sex)%>%
    pull()
  
  #pull ranks of the mouse from every CC
  rank_act_vec <- rank_tibble%>%
    filter(Animal_ID==mouse)%>%
    select(SystemRank_act)%>%
    pull()
  #sum() 
  
  rank_inact_vec <- rank_tibble%>%
    filter(Animal_ID==mouse)%>%
    select(SystemRank_inact)%>%
    pull()
  #sum() 
  
  rank_total_vec <- rank_tibble%>%
    filter(Animal_ID==mouse)%>%
    select(SystemRank_total)%>%
    pull() 
  #sum() 
  
  cat("\n", mouse)
  print(group)
  cat("active ranks: ", rank_act_vec, ", ")
  cat("inactive ranks: ", rank_inact_vec, ", ")
  cat("total ranks: ", rank_total_vec, "\n")
  
  #test if mouse data from all 4 CCs was tracked(if the system was always complete)
  if(length(rank_act_vec) !=4) {#is not important which vector is getting tested, all should have the same length
    
    cat(batch, ", ", mouse, " the animal has not 4 tracked cage changes", "\n")
    #1: skip and go to the next iteration of for loop?
    next
    #2: add similar values until vector length is 4?
    
    #3: change calc function that you alway get the average as a score-> then we need a variable for num of CCs
  }
  #compute score_vecs (now rank 1 means score 4 and vice versa)
  score_act_vec <- translate_rank_in_score_vec(rank_act_vec)
  score_inact_vec <- translate_rank_in_score_vec(rank_inact_vec) 
  score_total_vec <- translate_rank_in_score_vec(rank_total_vec) 
  
  
  
  social_score_tibble <- social_score_tibble%>%
    mutate(Group = ifelse(Animal_ID == mouse, group, Group))%>%
    mutate(Batch = ifelse(Animal_ID == mouse, batch, Batch))%>%
    mutate(Sex = ifelse(Animal_ID == mouse, sex, Sex))%>%
    mutate(Social_score_act = ifelse(Animal_ID == mouse, sum(score_act_vec), Social_score_act))%>%
    mutate(Social_score_inact = ifelse(Animal_ID == mouse, sum(score_inact_vec), Social_score_inact))%>%
    mutate(Social_score_total = ifelse(Animal_ID == mouse, sum(score_total_vec), Social_score_total))%>%
    mutate(Std_dev_act = ifelse(Animal_ID == mouse, sd(score_act_vec), Std_dev_act))%>% #standard deviation with sample
    mutate(Std_dev_inact = ifelse(Animal_ID == mouse, sd(score_inact_vec), Std_dev_inact))%>% 
    mutate(Std_dev_total = ifelse(Animal_ID == mouse, sd(score_total_vec), Std_dev_total))
  
  
}

# delete rows with NA in tibble (incomplete systems)
social_score_tibble <- na.omit(social_score_tibble)


################## ANALYSIS OF THE DERIVATION ##################

#animals with low std. derivation(-> consistent rank):
low_deriv_tibble <- social_score_tibble%>%
  #filter(Social_score_total>12)%>%  #very consistent social 
  #filter(Social_score_total<8)%>%  #very consistent insocial
  filter(Std_dev_total<=0.75)%>%
  select(-c(Social_score_act, Social_score_inact, Std_dev_act, Std_dev_inact))

#animals with high std. derivation(-> inconsistent rank):
high_deriv_tibble <- social_score_tibble%>%
  filter(Std_dev_total>=1.4)%>%
  select(-c(Social_score_act, Social_score_inact, Std_dev_act, Std_dev_inact))


#active and inactive
low_deriv_tibble <- social_score_tibble%>%
  filter(Std_dev_act<=0.75)

high_deriv_tibble <- social_score_tibble%>%
  filter(Std_dev_inact>=1.4)

################################
#kendalls tau

#plots  graphic correlation
plot(social_score_tibble$Social_score_total, social_score_tibble$Std_dev_total)
#calculates corr with kendalls tau
cor(social_score_tibble$Social_score_total, social_score_tibble$Std_dev_total,method="kendall")
#0.1486019 -> means positive correlation

#one try for order of cagechange vs Ranks of one individual:
cor(c(1,2,3,4), c(4,4,1,3),method="kendall")
#-0.5477226 

###################################
#lme
library(lmerTest)
library(emmeans)

model <- lmerTest::lmer(Std_dev_total ~ Group * Sex + (1 | Group), data = social_score_tibble)
emmeans_obj <- emmeans(model, ~ Group, data = social_score_tibble, 
                       cov.reduce = FALSE, adjust = "sidak", pbkrtest.limit = 10000)
pairwise_results <- pairs(emmeans_obj, by = c("Group"), adjust = "bonferroni")

cat("Mixed-effects model results :\n")
print(summary(model))

cat("emmeans results :\n")
print(emmeans_obj)

# Perform pairwise comparisons with Tukey adjustment
cat("Pairwise day comparisons :\n")
print(pairwise_results)

################## PLOTS ##########

## SOCIAL PROX IN SOCIAL SCORE ##


#barplot single animals and their score, colored in groups -> shitty
total_social_score_barplot <- ggplot(data=social_score_tibble, aes(x = factor(Animal_ID, levels = unique(Animal_ID[order(Group)])), y = Social_score_total, fill=Group))+
  geom_bar(stat = "identity", position = "dodge")+
  scale_x_discrete("Animal ID")+
  facet_grid(Sex~.)


#delete the rows with con group for this plot. not necessary
filtered_social_score_tibble <- social_score_tibble%>%
  filter(Group != "con")#%>%
  ##filter(Batch != "B5")%>%
  #filter(Batch != "B6")
#scatterplot groups and their score, colored in groups
total_social_score_scatterplot <- ggplot(data = filtered_social_score_tibble, aes(x = Group, y = Social_score_total, color = Group, shape = Batch)) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2) +
  scale_shape_manual(values = c(3, 16, 17, 15, 5, 20)) + # definition of the forms(shapes) (z.B. 17 =triangle , 16 = circle)
  scale_x_discrete("Animal ID") +
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
  theme(legend.key.size = unit(5, "lines"), # Legendenliniengröße einstellen
        legend.title = element_text(size = 20), # Legendenüberschriftsgröße einstellen
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 20), # Achsenbeschriftungsgröße einstellen
        axis.title = element_text(size = 22), # Achsentitelgröße einstellen
        strip.text = element_text(size = 20)) 


#heatmap for sum-score in comparison to standard deviation of the single scores from the sum

#select needed columns and bring them in fitting order
selected_score_tibble <- select(social_score_tibble, c(Animal_ID, Group, Sex, Social_score_total, Std_dev_total))
melt_score_tibble <- melt(selected_score_tibble, id='Animal_ID')
#names(melt_score_tibble) <- c('Animal_ID', 'Value', 'amount')

total_social_score_heatmap <- ggplot(melt_score_tibble, aes(x = variable, y = Animal_ID, fill=value)) +
  geom_tile() +
  #scale_x_continuous(breaks = c(0, 100, 200, 300), labels = c("0", "100", "200", "300")) +
  #scale_y_continuous(breaks = c(0, 116), labels = c("0", "116")) +
  #scale_fill_gradientn(colors = c("yellow","orange", "red", "darkred", "#290000")) + 
  #limits are the borders of the scale,breaks are actual value breaks, labels are names for breakpoints
  labs(title = paste("Social Score in Comparison to Standard Deviation"))
#size of scale!!!


total_social_score_heatmap <- ggplot(melt_score_tibble, aes(x = variable, y = Animal_ID, fill=value)) +
  geom_tile() +
  facet_grid(~ variable, scales = "free", space = "free") +
  scale_fill_manual(values = c("male" = "blue", "female" = "pink", "sus" = "red", "res" = "green", "con" = "purple"), na.value = NA) +
  labs(title = paste("Social Score in Comparison to Standard Deviation"))

#################################################################################################only cont values
filter_melt <- select(social_score_tibble, c(Animal_ID, Social_score_total, Std_dev_total))
melt_score_tibble <- melt(filter_melt, id='Animal_ID')

total_social_score_heatmap <- ggplot(melt_score_tibble, aes(x = variable, y = Animal_ID, fill=value)) +
  geom_tile() +
  facet_grid(~ variable, scales = "free", space = "free") +
  scale_fill_gradientn(colors = c("yellow","orange", "red", "darkred", "#290000")) +
  labs(title = paste("Social Score in Comparison to Standard Deviation"))
####################################################################################################pheatmap
library(pheatmap)


# arrange columns
conditions_sorted <- selected_score_tibble%>% arrange( Sex, Group, Social_score_total)

#select std_dev
data_subset <- select(conditions_sorted, c(Animal_ID,  Std_dev_total))%>%
  column_to_rownames(var="Animal_ID") %>% 
  as.matrix()

# create legend for heatmap
df_legend <- as.data.frame(conditions_sorted)[, c("Animal_ID",
                                                  "Social_score_total",
                                                  "Group",
                                                  "Sex"),  
                                              drop = FALSE] %>%
  column_to_rownames(var="Animal_ID") %>% 
  `colnames<-`(c( "Social_score_total", "Group", "Sex"))


row_lables <-  base::paste(selected_score_tibble$Animal_ID)
#row_lables <-  base::paste(row.names(df_legend))

pheatmap::pheatmap(data_subset, 
                   cluster_rows = F, #originally T, maybe change again
                   cluster_cols = F,
                   clustering_method = T,
                   show_rownames = T,
                   legend = T,
                   fontsize_col = 8,
                   fontsize_row = 8,
                   annotation_row = df_legend, #data frame that specifies the annotations shown on left side of the heatmap. Each row defines the features for a specific row. The rows in the data and in the annotation are matched using corresponding row names. Note that color schemes takes into account if variable is continuous or discrete.
                   labels_row = row_lables
                   )

#annotation_colors	-list for specifying annotation_row and annotation_col track colors manually. It is possible to define the colors for only some of the features. Check examples for details.


####################################################################################################
#std dev plot

#scatterplot like the on for the social scor only for std_dev
# with con group

#scatterplot groups and their score, colored in groups
total_std_dev_scatterplot <- ggplot(data = social_score_tibble, aes(x = Group, y = Std_dev_total, color = Group, shape = Batch)) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2) +
  scale_shape_manual(values = c(3, 16, 17, 15, 5, 20)) + # definition of the forms(shapes) (z.B. 17 =triangle , 16 = circle)
  scale_x_discrete("Animal ID") +
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
  theme(legend.key.size = unit(5, "lines"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 20), 
        axis.title = element_text(size = 22), 
        strip.text = element_text(size = 20)) 


########################################################################################################################
## social rank avg##
#tibble for all batches together on basis of the rank tibble

unique_mice <- unique(rank_tibble$Animal_ID)

avg_social_rank_tibble <- tibble(Animal_ID=unique_mice,
                              Group=NA,
                              Batch=NA,
                              Sex=NA,
                              Avg_social_rank_act=NA,
                              Avg_social_rank_inact=NA,
                              Avg_social_rank_total=NA,
                              Std_dev_act=NA,
                              Std_dev_inact=NA,
                              Std_dev_total=NA,
                              CV_total=NA) 


for(mouse in unique_mice){
  
  #define group of animal
  group <- ifelse(mouse %in% sus_animals, "sus", ifelse(mouse %in% con_animals, "con", "res"))
  
  #define batch of animal
  batch <- rank_tibble%>%
    filter(Animal_ID==mouse)%>%
    head(1)%>%
    select(Batch)%>%
    pull()
  #define sex of animal
  sex <- rank_tibble%>%
    filter(Animal_ID==mouse)%>%
    head(1)%>%
    select(Sex)%>%
    pull()
  
  #pull ranks of the mouse from every CC
  rank_act_vec <- rank_tibble%>%
    filter(Animal_ID==mouse)%>%
    select(SystemRank_act)%>%
    pull()
  #sum() 
  
  rank_inact_vec <- rank_tibble%>%
    filter(Animal_ID==mouse)%>%
    select(SystemRank_inact)%>%
    pull()
  #sum() 
  
  rank_total_vec <- rank_tibble%>%
    filter(Animal_ID==mouse)%>%
    select(SystemRank_total)%>%
    pull() 
  #sum() 
  
 
  if(length(rank_act_vec) == 1 | length(rank_inact_vec) == 1 | length(rank_total_vec) == 1) {cat(batch, mouse, "\n")}
  
  avg_social_rank_tibble <- avg_social_rank_tibble%>%
    mutate(Group = ifelse(Animal_ID == mouse, group, Group))%>%
    mutate(Batch = ifelse(Animal_ID == mouse, batch, Batch))%>%
    mutate(Sex = ifelse(Animal_ID == mouse, sex, Sex))%>%
    mutate(Avg_social_rank_act = ifelse(Animal_ID == mouse, sum(rank_act_vec)/length(rank_act_vec), Avg_social_rank_act))%>%
    mutate(Avg_social_rank_inact = ifelse(Animal_ID == mouse, sum(rank_inact_vec)/length(rank_inact_vec), Avg_social_rank_inact))%>%
    mutate(Avg_social_rank_total = ifelse(Animal_ID == mouse, sum(rank_total_vec)/length(rank_total_vec), Avg_social_rank_total))%>%
    mutate(Std_dev_act = ifelse(Animal_ID == mouse, sd(rank_act_vec), Std_dev_act))%>% #standard deviation with sample
    mutate(Std_dev_inact = ifelse(Animal_ID == mouse, sd(rank_inact_vec), Std_dev_inact))%>% 
    mutate(Std_dev_total = ifelse(Animal_ID == mouse, sd(rank_total_vec), Std_dev_total))%>%
    mutate(CV_total = ifelse(Animal_ID == mouse, sd(rank_total_vec)/mean(rank_total_vec), CV_total))
  
  
}

# delete rows with NA in tibble (mice with only one rank, bc of incomplete systems)
avg_social_rank_tibble <- na.omit(avg_social_rank_tibble)
########################################################################################################################
### PLOTS AVG TIBBLE ###

#delete the rows with con group for this plot. not necessary
filtered_avg_social_rank_tibble <- avg_social_rank_tibble%>%
  filter(Group != "con")

#scatterplot groups and their score, colored in groups
total_avg_social_rank_tibble_scatterplot <- ggplot(data = filtered_avg_social_rank_tibble, aes(x = Group, y = Avg_social_rank_total, color = Group, shape = Batch)) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2) +
  scale_shape_manual(values = c(3, 16, 17, 15, 5, 20)) + # definition of the forms(shapes) (z.B. 17 =triangle , 16 = circle)
  scale_x_discrete("Animal ID") +
  labs(title = paste("Social measurement in avg of the total ranks from each cageChange"), subtitle = paste("rank1 = the most social", "\n", "rank4 = the least social"))+
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
        legend.key.size = unit(5, "lines"),
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 20), 
        axis.title = element_text(size = 22), 
        strip.text = element_text(size = 20)) 

#active vs inactive
act_avg_social_rank_tibble_scatterplot <- ggplot(data = filtered_avg_social_rank_tibble, aes(x = Group, y = Avg_social_rank_act, color = Group, shape = Batch)) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2) +
  scale_shape_manual(values = c(3, 16, 17, 15, 5, 20)) + 
  scale_x_discrete("Animal ID") +
  labs(title = paste("Social measurement in avg of the active ranks from each cageChange"), subtitle = paste("rank1 = the most social", "\n", "rank4 = the least social"))+
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
        legend.key.size = unit(5, "lines"),
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 20), 
        axis.title = element_text(size = 22), 
        strip.text = element_text(size = 20)) 

inact_avg_social_rank_tibble_scatterplot <- ggplot(data = filtered_avg_social_rank_tibble, aes(x = Group, y = Avg_social_rank_inact, color = Group, shape = Batch)) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2) +
  scale_shape_manual(values = c(3, 16, 17, 15, 5, 20)) + 
  scale_x_discrete("Animal ID") +
  labs(title = paste("Social measurement in avg of the inactive ranks from each cageChange"), subtitle = paste("rank1 = the most social", "\n", "rank4 = the least social"))+
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
        legend.key.size = unit(5, "lines"),
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 20), 
        axis.title = element_text(size = 22), 
        strip.text = element_text(size = 20)) 

###std dev plot
total_std_dev_scatterplot <- ggplot(data = avg_social_rank_tibble, aes(x = Group, y = Std_dev_total, color = Group, shape = Batch)) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2) +
  scale_shape_manual(values = c(3, 16, 17, 15, 5, 20)) + 
  scale_x_discrete("Animal ID") +
  labs(title = paste("Standard deviation\n of total ranks from each cageChange"), subtitle = paste("how stable is the rank\n of each individual?"))+
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
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        title = element_text(size = 16),
        legend.key.size = unit(3, "lines"),
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16),
        axis.text = element_text(size = 20), 
        axis.title = element_text(size = 22), 
        strip.text = element_text(size = 20))

### CV plot
total_CV_scatterplot <- ggplot(data = avg_social_rank_tibble, aes(x = Group, y = CV_total, color = Group, shape = Batch)) +
  geom_jitter(size = 4, alpha = 0.7, width = 0.2) +
  scale_shape_manual(values = c(3, 16, 17, 15, 5, 20)) + 
  scale_x_discrete("Animal ID") +
  labs(title = paste("Coefficient of variation\n of ranks from each cageChange"), subtitle = paste("how stable is the rank\n of each individual?"))+
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
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        title = element_text(size = 16),
        legend.key.size = unit(3, "lines"),
        legend.title = element_text(size = 16), 
        legend.text = element_text(size = 16),
        axis.text = element_text(size = 20), 
        axis.title = element_text(size = 22), 
        strip.text = element_text(size = 20))

