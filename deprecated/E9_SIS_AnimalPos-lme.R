# LME file
# Workin with one batch at a time, all the CageChanges

# package management via pacman
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  readr, tibble, dplyr, reshape2, forcats, lme4, nlme,
  lmerTest, emmeans,  ggplot2
)


# paths
#working_directory <- "S:/Lab_Member/Anja/Git/MDC_Bachelor/E9_SIS_AnimalPos"
#working_directory <- "/home/anja/Dokumente/FU BERLIN/BA/Git/MDC_Bachelor/E9_SIS_AnimalPos"
working_directory <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability"

# ensure output dir and compact file path setup
out_dir <- file.path(working_directory, "lme", "social_prox")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- setNames(
  file.path(out_dir, sprintf("%s_%s.txt",
                             c("mixed_model_results", "pairwise_results", "emmeans_results"),
                             PhaseValue)),
  c("mixed_model_results_file", "pairwise_results_file", "emmeans_results_file")
)

invisible(lapply(files, function(p) if (!file.exists(p)) file.create(p)))
list2env(as.list(files), envir = environment())

## general definitions and read data ###########################################################
# define batch and cage change
batches <- c("B1", "B2", "B3", "B4", "B5")
cageChanges <- c("CC1", "CC2", "CC3", "CC4")


#sus and con animals
sus_animals <- readLines(paste0(working_directory, "/raw_data/sus_animals.csv"))
##csv path for con animals
con_animals <- readLines(paste0(working_directory, "/raw_data/con_animals.csv"))


#define sus res and con groups
#sus_in_batch <- intersect(unique_mice, sus_animals)
#con_in_batch <- intersect(unique_mice, con_animals)
#res_in_batch <- setdiff(setdiff(unique_mice,sus_in_batch), con_in_batch)


####################################################################################
########## SOCIAL  PROXIMITY ##################
#create general tibble for the batch

#initialize empty tibbles
batch_prox_active <- tibble(AnimalID = character(),
                            Sex = character(),
                            System = character(),
                            Batch = character(),
                            Group = character(),
                            Phase = numeric(),
                            SocialProx = numeric())

batch_prox_inactive <- tibble(AnimalID = character(),
                              Sex = character(),
                              System = character(),
                              Batch = character(),
                              Group = character(),
                              Phase = numeric(),
                              SocialProx = numeric())

batch_prox <- tibble(AnimalID = character(),
                     Sex = character(),
                     System = character(),
                     Batch = character(),
                     Group = character(),
                     Phase = numeric(),
                     SocialProx = numeric())

for (batch in batches) {
  for (cageChange in cageChanges) {
    sex <- ifelse(batch %in% c("B3", "B4", "B6"), "female", "male")

    #read csv
    #mouse_which_system_tibble <-  tibble(read_delim(paste0(working_directory,"/tables/", batch, "_", cageChange, "_mouse_which_system_tibble.csv"), delim = ",", show_col_types = FALSE))
    mouse_which_system_tibble <-  tibble(read_delim(paste0(working_directory,"/tables/", batch, "_", cageChange, "_animal_ids.csv"),
                                                    delim = ",", show_col_types = FALSE))
    #total_closeness_table <- tibble(read_delim(paste0(working_directory,"/tables/", batch, "_", cageChange, "_total_closeness_table.csv"), delim = ",", show_col_types = FALSE))
    total_closeness_table <- tibble(read_delim(paste0(working_directory,"/tables/", batch, "_", cageChange, "_total_proximity.csv"),
                                               delim = ",", show_col_types = FALSE))

    #use only A1,I2,A2 for the computation of the rank from CC4
    if (cageChange == "CC4") {

      total_closeness_table <- total_closeness_table %>%
        filter(Phase %in% c("A1", "I2", "A2"))
    }

    # divide into act and inactive
    active_data <- total_closeness_table %>%
      filter(grepl("^A\\d+", Phase))
    inactive_data <- total_closeness_table %>%
      filter(grepl("^I\\d+", Phase))

    # melt tibbles
    active_data <- melt(active_data, id = "Phase")
    names(active_data) <- c("Phase", "AnimalID", "SocialProx")

    inactive_data <- melt(inactive_data, id = "Phase")
    names(inactive_data) <- c("Phase", "AnimalID", "SocialProx")

    #original table, without division
    total_closeness_table$Phase <- 1:nrow(total_closeness_table)#rename both phases into numbers
    total_closeness_table <- melt(total_closeness_table, id = "Phase")
    names(total_closeness_table) <- c("Phase", "AnimalID", "SocialProx")

    #add extra columns
    active_data <- active_data %>%
      mutate(Group = ifelse(AnimalID %in% sus_animals, "sus",
                            ifelse(AnimalID %in% con_animals,
                                   "con", "res"))) %>%
      mutate(Sex = sex) %>%
      mutate(Batch = batch) %>%
      mutate(AnimalID = as.character(AnimalID)) %>%
      mutate(System = purrr::map_chr(AnimalID, ~ {
        index <- which(mouse_which_system_tibble == .x, arr.ind = TRUE)[2]
        if (length(index) > 0) {
          names(mouse_which_system_tibble)[index]
        } else {NA}
      })) %>%
      mutate(Phase = as.numeric(substr(Phase, 2, 3))) #take second char from phase and make it numeric

    inactive_data <- inactive_data %>%
      mutate(Group = ifelse(AnimalID %in% sus_animals, "sus",
                            ifelse(AnimalID %in% con_animals,
                                   "con", "res"))) %>%
      mutate(Sex = sex) %>%
      mutate(Batch = batch) %>%
      mutate(AnimalID = as.character(AnimalID)) %>%
      mutate(System = purrr::map_chr(AnimalID, ~ {
        index <- which(mouse_which_system_tibble == .x, arr.ind = TRUE)[2]
        if (length(index) > 0) {
          names(mouse_which_system_tibble)[index]
        } else {NA}
      })) %>%
      mutate(Phase = as.numeric(substr(Phase, 2, 3))) %>% #take second char from phase and make it numeric
      mutate(Phase = Phase - 1) #renaming phases starting with number one(bc inact phases always started with 2)

    total_closeness_table <- total_closeness_table %>%
      mutate(Group = ifelse(AnimalID %in% sus_animals, "sus",
                            ifelse(AnimalID %in% con_animals,
                                   "con", "res"))) %>%
      mutate(Sex = sex) %>%
      mutate(Batch = batch) %>%
      mutate(AnimalID = as.character(AnimalID)) %>%
      mutate(System = purrr::map_chr(AnimalID, ~ {
        index <- which(mouse_which_system_tibble == .x, arr.ind = TRUE)[2]
        if (length(index) > 0) { names(mouse_which_system_tibble)[index]
        } else {NA}
      }))

    #sort tibble into correct orde of columns
    active_data <- active_data[c("AnimalID", "Sex", "System", "Group",
                                 "Batch", "Phase", "SocialProx")]
    inactive_data <- inactive_data[c("AnimalID", "Sex", "System", "Group",
                                     "Batch", "Phase", "SocialProx")]
    total_closeness_table <- total_closeness_table[c("AnimalID", "Sex", "System", "Group",
                                                     "Batch", "Phase", "SocialProx")]

    ## combine with general tibble ##
    #first check if data already in tibble (esp. the Phases)
    #phase number eventually needs to be changed for consecutive phase number
    unique_phase_number_act <- batch_prox_active %>%
      filter(Batch == batch) %>%
      pull(Phase) %>%
      unique()
    if (sex == "female") {print(unique_phase_number_act)}
    if (length(unique_phase_number_act) == 0) {#general tibble is empty, thus first cagChange and nothing has to be alterated
      #print("first cageChange")
      batch_prox_active <- bind_rows(batch_prox_active, active_data)

    } else {#not the first cageChange, phase number has to be alterated in current tibble for consecutivity
      last_consec <- max(unique_phase_number_act)
      #print(last_consec)
      active_data <- active_data %>%
        mutate(Phase = Phase + last_consec)
      #after that, add rows to general tibble
      batch_prox_active <- bind_rows(batch_prox_active, active_data)
    }

    #same again for inactive
    unique_phase_number_inact <- batch_prox_inactive %>%
      filter(Batch == batch) %>%
      pull(Phase) %>%
      unique()

    #print(unique_phase_number_inact)
    if (length(unique_phase_number_inact) == 0) {
     # print("first cageChange")
      batch_prox_inactive <- bind_rows(batch_prox_inactive, inactive_data)
    } else {
      last_consec <- max(unique_phase_number_inact)
      #print(last_consec)
      inactive_data <- inactive_data %>%
        mutate(Phase = Phase + last_consec)
      #print(inactive_data)
      #after that, add rows to general tibble
      batch_prox_inactive <- bind_rows(batch_prox_inactive, inactive_data)
    }

    #for original table
    unique_phase_number <- batch_prox %>%
      filter(Batch == batch) %>%
      pull(Phase) %>%
      unique()
    #print(unique_phase_number)
    if (length(unique_phase_number) == 0) {
      #print("first cageChange")
      #print(total_closeness_table)
      batch_prox <- bind_rows(batch_prox, total_closeness_table)
    } else {
      last_consec <- max(unique_phase_number)
      #print(last_consec)
      total_closeness_table <- total_closeness_table %>%
        mutate(Phase = Phase + last_consec)
      #print(total_closeness_table)
      #after that, add rows to general tibble
      batch_prox <- bind_rows(batch_prox, total_closeness_table)
    }
  }
}

###################### LME ANALYSIS #################################
### for one seperate batch -> change for loop with batches to only one batch ###
#https://www.youtube.com/watch?v=u1ePV1ntMNs&t=916s

library(lme4)
#with two fixed values
#Group*Phase
model6 <- lmer(SocialProx ~ Group * Phase +
                 (1 | AnimalID) +
                 (1 | System),
               data = batch_prox_active)

#same without group value
model6.null <- lmer(SocialProx ~ Phase +
                      (1 | AnimalID) +
                      (1 | System),
                    data = batch_prox_active)

#do an anova to test if group variable is significant
anova(model6.null, model6)
#value is not significant...0.4985

#same without Phase value
model6.null= lmer(SocialProx ~ Group +
                    (1 | AnimalID) +
                    (1 | System),
                  data = batch_prox_active)

#do an anova to test if group variable is significant
anova(model6.null, model6)
# value 0.1591 ...

#https://www.youtube.com/watch?v=VhMWPkTbXoY
library(nlme)
#first, a linear model with only fixed effects
model1 <- lm(SocialProx ~ Group * Phase, data = batch_prox_active)
summary(model1)
coef(model1)

#now also with random effects and random variation from systems
model2 <- lme(SocialProx ~ Group * Phase, data = batch_prox_active, random = ~ 1 | System / Group)
summary(model2)

#coefficients
coef(model2)
#plot of random effects
plot(ranef(model2))
#this plot should show basically nothing
# should show symmetrical spread of random effects around 0 (expected value is always 0)
# no particular pattern

#plot residuals
plot(model2)
#should ideally scattered around 0, symmetrical, with no obvious patterns

#other random effect??
model3 <- lme(SocialProx~Group*Phase, data=batch_prox_active, random=~1|AnimalID/Group)
summary(model3)

model4 <- lme(SocialProx~Group*Phase, data=batch_prox_inactive, random=~1|System/Group)
summary(model4)
plot(ranef(model4))
plot(model4)

###lme with both phases together
model5 <- lme(SocialProx ~ Group * Phase, data = batch_prox, random = ~ 1 | AnimalID / Group)
summary(model5)

#######lme on data with all batches together(when all batches included)#######
#library(nlme)
batch_prox_na_free <- na.omit(batch_prox)
batch_prox_na_free_active <- na.omit(batch_prox_active)
batch_prox_na_free_inactive <- na.omit(batch_prox_inactive)

model7 <- lme(SocialProx ~ Group * Phase * Sex, data = batch_prox_na_free, random = ~ 1 | AnimalID / Group)
summary(model7)

model8 <- lme(SocialProx ~ Group * Phase * Sex, data = batch_prox_na_free_active, random = ~ 1 | AnimalID / Group)
summary(model8)


model9 <- lme(SocialProx ~ Group * Phase * Sex, data = batch_prox_na_free_inactive, random = ~ 1 | AnimalID / Group)
summary(model9)
coef(model9)

#lmer, library(lme4)
model10 = lmer(SocialProx ~ Group * Phase + Sex + (1 | AnimalID), data = batch_prox_na_free_active)
model10

model11 = lmer(SocialProx ~ Group * Phase + Sex + (1 | AnimalID), data = batch_prox_na_free_inactive)
model11

########################################################
## SCRIPT FROM TOBI ##
##################### MIXED MODELS #####################
# Perform mixed-effects models with pairwise comparisons
library(lmerTest)
library(emmeans)

#model_results <- list()
#emmeans_results <- list()

for (PhaseValue in c("Active", "Inactive")) {
  #data_filtered_agg_Phase_subset <- data_filtered_agg %>%
  #  filter(Phase == PhaseValue)

  #data_filtered_agg_Phase_subset <- data_filtered_agg_Phase_subset %>%
  #  mutate(Group = ifelse(Group == "SIS" & SUS == FALSE, "RES", 
  #                        ifelse(Group == "SIS" & SUS == TRUE, "SUS", Group)))

  if (PhaseValue == "Active") {
    #model <- lmerTest::lmer(SocialProx ~ Group * Phase + (1 | AnimalID), data = batch_prox_na_free_active)
    #model <- lmerTest::lmer(SocialProx ~ Group * Phase + Sex + (1 | AnimalID), data = batch_prox_na_free_active) #sex added as factor
    model <- lmerTest::lmer(SocialProx ~ Group * Phase * Sex + (1 | AnimalID), data = batch_prox_na_free_active) # * means INTERACTION between all three factors
    emmeans_obj <- emmeans(model, ~ Group * Phase, data = batch_prox_na_free_active, 
                           cov.reduce = FALSE, adjust = "sidak", pbkrtest.limit = 10000)
    pairwise_results <- pairs(emmeans_obj, by = c("Phase"), adjust = "bonferroni")
  } else if (PhaseValue == "Inactive") {
    #model <- lmerTest::lmer(SocialProx ~ Group * Phase + (1 | AnimalID), data = batch_prox_na_free_inactive)
    #model <- lmerTest::lmer(SocialProx ~ Group * Phase + Sex + (1 | AnimalID), data = batch_prox_na_free_inactive)
    model <- lmerTest::lmer(SocialProx ~ Group * Phase * Sex + (1 | AnimalID), data = batch_prox_na_free_inactive)
    emmeans_obj <- emmeans(model, ~ Group * Phase, data = batch_prox_na_free_inactive, 
                           cov.reduce = FALSE, adjust = "sidak", pbkrtest.limit = 10000)
    pairwise_results <- pairs(emmeans_obj, by = c("Phase"), adjust = "bonferroni")
  }

  #model_results[[PhaseValue]] <- model

  # save mixed model results to a file
  
  model_summary <- capture.output(summary(model))
  writeLines(model_summary, con = mixed_model_results_file)

  cat("Mixed-effects model results for", PhaseValue, ":\n")
  print(summary(model))

  cat("emmeans results for", PhaseValue, ":\n")
  print(emmeans_obj)

  # Perform pairwise comparisons with Tukey adjustment
  cat("Pairwise day comparisons for", PhaseValue, ":\n")
  print(pairwise_results)

  # Save pairwise results to a file
  
  pairwise_summary <- capture.output(print(pairwise_results))
  writeLines(pairwise_summary, con = pairwise_results_file)

  # Save emmeans results to a file
  
  emmeans_summary <- capture.output(summary(emmeans_obj))
  writeLines(emmeans_summary, con = emmeans_results_file)
}

####################################################################################
### Consec Plot :Social hours per day ###

#change counted seconds to hours and remove na rows
plotData_batch_prox <- batch_prox %>%
  mutate(SocialProx = ifelse(SocialProx != 0, SocialProx / 3600, SocialProx)) %>%
  na.omit()

plotData_batch_prox_active <- batch_prox_active %>%
  mutate(SocialProx = ifelse(SocialProx != 0, SocialProx / 3600, SocialProx)) %>%
  na.omit() %>%
  group_by(Sex, Group, Phase) %>%
  summarise(Mean_prox = mean(SocialProx))

plotData_batch_prox_inactive <- batch_prox_inactive %>%
  mutate(SocialProx = ifelse(SocialProx != 0, SocialProx / 3600, SocialProx)) %>%
  na.omit() %>%
  group_by(Sex, Group, Phase) %>%
  summarise(Mean_prox = mean(SocialProx))

#rename phase columns
names(plotData_batch_prox)[names(plotData_batch_prox) == "Phase"] <- "ConsecPhase"
names(plotData_batch_prox_active)[names(plotData_batch_prox_active) == "Phase"] <- "ConsecActive"
names(plotData_batch_prox_inactive)[names(plotData_batch_prox_inactive) == "Phase"] <- "ConsecInactive"

#lineplot active phases
total_active_hours_per_day_lineplot <-  ggplot(data = plotData_batch_prox_active,
                                               aes(x = fct_inorder(as_factor(ConsecActive)),
                                               y = Mean_prox, group = Group, color = Group)) +
  geom_line() +
  scale_color_manual(values = c("sus" = "#E63946",
                                "res" = "grey60",
                                "con" = "#457B9D"),
                     name = "group") +
  scale_y_continuous("social proximity in h") +
  scale_x_discrete("ConsecActive") +
  labs(title = paste("Social measurement in hours per active phase(4days)")) +
  theme_bw() +
  facet_grid(Sex ~ .) +
  theme(title = element_text(size = 20),
        legend.key.size = unit(5, "lines"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        strip.text = element_text(size = 20))

#lineplot inactive phases
total_inactive_hours_per_day_lineplot <-  ggplot(data = plotData_batch_prox_inactive,
                                                 aes(x = fct_inorder(as_factor(ConsecInactive)),
                                                     y = Mean_prox, group = Group, color = Group)) +
  geom_line() +
  scale_color_manual(values = c("sus" = "#E63946",
                                "res" = "grey60",
                                "con" = "#457B9D"),
                     name = "group") +
  scale_y_continuous("social proximity in h") +
  scale_x_discrete("ConsecInactive") +
  labs(title = paste("Social measurement in hours per inactive phase(4days)")) +
  theme_bw() +
  facet_grid(Sex ~ .) +
  theme(title = element_text(size = 20),
        legend.key.size = unit(5, "lines"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        strip.text = element_text(size = 20))

####################################################################################
####### COIL CROSSING ##############
#create general tibble for the batch

#initialize empty tibble
batch_cross <- tibble(AnimalID = character(),
                      System = character(),
                      Group = character(),
                      Phase = numeric(),
                      Coilcross = numeric())

for (batch in batches) {
  for (cageChange in cageChanges) {
    # read csv
    mouse_which_system_tibble <- tibble(
      readr::read_delim(
        paste0(
          working_directory, "/tables/",
          batch, "_", cageChange, "_animal_ids.csv"
        ),
        delim = ",",
        show_col_types = FALSE
      )
    )
    total_movement_table <- tibble(
      readr::read_delim(
        paste0(
          working_directory, "/tables/",
          batch, "_", cageChange, "_total_movement.csv"
        ),
        delim = ",",
        show_col_types = FALSE
      )
    )

    # cut the last five columns
    total_movement_table <- dplyr::select(
      total_movement_table,
      -dplyr::starts_with("sys.")
    )

    # melt tibble
    total_movement_table$Phase <- seq_len(nrow(total_movement_table)) # rename phases
    total_movement_table <- reshape2::melt(total_movement_table, id = "Phase")
    names(total_movement_table) <- c("Phase", "AnimalID", "Coilcross")

    # add extra columns
    total_movement_table <- total_movement_table %>%
      dplyr::mutate(
        Group = ifelse(
          AnimalID %in% sus_animals, "sus",
          ifelse(AnimalID %in% con_animals, "con", "res")
        ),
        AnimalID = as.character(AnimalID),
        System = purrr::map_chr(
          AnimalID,
          ~ {
            index <- which(
              mouse_which_system_tibble == .x,
              arr.ind = TRUE
            )[2]
            if (length(index) > 0) {
              names(mouse_which_system_tibble)[index]
            } else {
              NA
            }
          }
        )
      )

    # sort tibble into correct order of columns
    total_movement_table <- total_movement_table[
      c("AnimalID", "System", "Group", "Phase", "Coilcross")
    ]

    # combine with general tibble
    unique_phase_number <- unique(batch_cross$Phase)
    if (length(unique_phase_number) == 0) {
      batch_cross <- dplyr::bind_rows(batch_cross, total_movement_table)
    } else {
      last_consec <- max(unique_phase_number)
      total_movement_table <- total_movement_table %>%
        dplyr::mutate(Phase = Phase + last_consec)
      batch_cross <- dplyr::bind_rows(batch_cross, total_movement_table)
    }
  }
}

### lme coil crossing
model6 <- nlme::lme(
  Coilcross ~ Group * Phase,
  data = batch_cross,
  random = ~ 1 | AnimalID / Group
)
summary(model6)
