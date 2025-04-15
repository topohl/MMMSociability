#' @title E9 Social Stress Experiment Data Analysis and Visualization
#' 
#' @description This script is designed to perform data analysis and visualization 
#' for the E9 Social Stress experiment. It includes functionality for loading 
#' necessary libraries, setting up directories for saving plots and statistics, 
#' and defining paths for working and saving directories.
#' 
#' @details 
#' - The script uses the `pacman` package to load required libraries. If `pacman` 
#'   is not installed, it will need to be installed first.
#' - Libraries used in this script include:
#'   - `readr`
#'   - `tibble`
#'   - `dplyr`
#'   - `reshape2`
#'   - `stringr`
#'   - `ggplot2`
#'   - `forcats`
#'   - `plotrix`
#'   - `lmerTest`
#' - Flags are set to determine whether to save plots and statistics.
#' - Paths for the working directory and saving directory are defined, with 
#'   options to uncomment and modify paths for different environments 
#'   (e.g., local machine, server).
#' - Subdirectories for saving plots and tables are also defined.
#' 
#' @note Ensure that the required libraries are installed and that the paths 
#' are correctly set up for your environment before running the script.
# improve comments professionally like data scientist

# Load required packages using pacman. Install pacman if not already installed.
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  readr, tibble, dplyr, reshape2, stringr, ggplot2, forcats, plotrix, lmerTest, 
  tidyr, gridExtra, ggpubr, cowplot, writexl, rstatix, lme4, purrr, stringi, scales
)

# Flags to control saving of plots and statistics
save_plots <- TRUE
save_statistics <- TRUE

# Define paths for working and saving directories
working_directory <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability"
saving_directory <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability"

# Define subdirectories for plots and tables
plots_directory <- "/plots"
tables_directory <- "/tables"

# Define specific directories for consecutive plots and mixed model results
consec_plot_dir <- paste0(working_directory, "/plots/consec_plots")
consec_lme_result_dir <- paste0(working_directory, "/lme/consec-mixed_model_results")
cc_lme_result_dir <- paste0(working_directory, "/lme/cc-mixed_model_results")

# Create directories if they do not already exist
dirs <- c(consec_plot_dir, consec_lme_result_dir, cc_lme_result_dir)
lapply(dirs, function(dir) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
})

# Load external functions from a separate script
source(paste0("C:/Users/topohl/Documents/GitHub/MMMSociability/Functions/E9_SIS_AnimalPos-functions.R"))

# General definitions and data loading ###########################################################

# Define experimental batches and cage changes
batches <- c("B1", "B2", "B3", "B4", "B5", "B6")
cageChanges <- c("CC1", "CC2", "CC3", "CC4") 

# Read lists of control (con) and susceptible (sus) animals from CSV files
con_animals <- readLines(paste0(working_directory, "/raw_data/con_animals.csv"))
sus_animals <- readLines(paste0(working_directory, "/raw_data/sus_animals.csv"))

# ------------------------------------------
# INITIALIZATION OF CONSECUTIVE PHASE TIBBLES
# ------------------------------------------
# This section initializes tibbles to store processed data for consecutive phases.
# These tibbles are structured to support:
# - Linear Mixed-Effects (LME) modeling
# - Visualization through plots
# The data includes metadata (e.g., CageChange, Batch, System, AnimalID, etc.)
# and metrics such as active/inactive phase durations and behavioral measures.

message("## Initializing tibbles for consecutive phase data ##")

# ------------------------------------------
# SOCIAL PROXIMITY DATA
# ------------------------------------------
# Initialize an empty tibble to store consecutive phase data for social proximity.
# This tibble will hold:
# - Metadata: CageChange, Batch, System, AnimalID, Sex, Group, Phase
# - Metrics: Consecutive active/inactive phase durations and social proximity values
social_prox_consec <- tibble(
  CageChange = character(),   # Identifier for the cage change
  Batch = character(),        # Experimental batch identifier
  System = character(),       # System identifier
  AnimalID = character(),     # Animal ID
  Sex = character(),          # Sex of the animal
  Group = character(),        # Group classification (e.g., con, sus, res)
  Phase = character(),        # Experimental phase (active/inactive)
  ConsecActive = numeric(),   # Duration of consecutive active phases
  ConsecInactive = numeric(), # Duration of consecutive inactive phases
  SocialProx = numeric()      # Social proximity value (in hours)
)

# ------------------------------------------
# COIL CROSSING DATA
# ------------------------------------------
# Initialize an empty tibble to store consecutive phase data for coil crossing.
# This tibble will hold:
# - Metadata: CageChange, Batch, System, AnimalID, Sex, Group, Phase
# - Metrics: Consecutive active/inactive phase durations and coil crossing values
coil_crossing_consec <- tibble(
  CageChange = character(),   # Identifier for the cage change
  Batch = character(),        # Experimental batch identifier
  System = character(),       # System identifier
  AnimalID = character(),     # Animal ID
  Sex = character(),          # Sex of the animal
  Group = character(),        # Group classification (e.g., con, sus, res)
  Phase = character(),        # Experimental phase (active/inactive)
  ConsecActive = numeric(),   # Duration of consecutive active phases
  ConsecInactive = numeric(), # Duration of consecutive inactive phases
  CoilCrossing = numeric()    # Coil crossing value (in counts)
)

# ------------------------------------------
# DATA CATEGORIES AND COLUMN MAPPING
# ------------------------------------------
# Define the data categories and their corresponding column names for processing.
# These categories will be used to dynamically process and populate the tibbles.
data_categories <- c("total_proximity", "total_movement")  # Data categories
values <- c("SocialProx", "CoilCrossing")                 # Corresponding column names

# Iterate over each data category (e.g., total proximity and total movement)
for (i in seq_along(data_categories)) {
  # Define the current data category being processed
  data_category <- data_categories[i]

  for (batch in batches) {
    max_consecAct <- 0
    max_consecInact <- 0

    for (cageChange in cageChanges) {
      sex <- ifelse(batch %in% c("B3", "B4", "B6"), "female", "male")

      # Load CSV files for each cage change
      system_animal_ids <- tibble(read_delim(paste0(working_directory, "/tables/", batch, "_", cageChange, "_animal_ids.csv"), delim = ",", show_col_types = FALSE))
      if (data_category == "total_proximity") {
        cage_change_data <- tibble(read_delim(paste0(working_directory, "/tables/", batch, "_", cageChange, "_total_proximity.csv"), delim = ",", show_col_types = FALSE))
      } else if (data_category == "total_movement") {
        cage_change_data <- tibble(read_delim(paste0(working_directory, "/tables/", batch, "_", cageChange, "_total_movement.csv"), delim = ",", show_col_types = FALSE))
      } else {
        stop("Invalid tibble name")
      }

      # Special handling for CC4 after A2 (grid in cage) - not usable for social analysis
      if (cageChange == "CC4") {
        cage_change_data <- cage_change_data %>% filter(Phase %in% c("A1", "I2", "A2"))
      }

      # Preprocess the data
      cage_change_data <- melt(cage_change_data, id = "Phase")
      names(cage_change_data) <- c("Phase", "AnimalID", values[i])

      # Create consecutive columns
      cage_change_data[c("Phase", "Consec")] <- str_split_fixed(cage_change_data$Phase, "", 2)
      cage_change_data <- cage_change_data %>%
        mutate(ConsecActive = ifelse(Phase == "A", as.numeric(Consec), 0)) %>%
        mutate(ConsecInactive = ifelse(Phase == "I", as.numeric(Consec), 0)) %>%
        mutate(Phase = ifelse(Phase == "A", "active", "inactive"))

      if (data_category == "total_proximity") {
        cage_change_data <- cage_change_data %>% mutate(SocialProx = ifelse(SocialProx != 0, SocialProx / 3600, SocialProx))
      }

      # Adjust consecutive values for higher cage changes
      if (cageChange != "CC1") {
        batch_table <- if (data_category == "total_proximity") {
          social_prox_consec %>% filter(Batch == batch)
        } else if (data_category == "total_movement") {
          coil_crossing_consec %>% filter(Batch == batch)
        } else {
          stop("Invalid tibble name")
        }

        max_consecAct <- batch_table %>% pull(ConsecActive) %>% unique() %>% max()
        max_consecInact <- batch_table %>% pull(ConsecInactive) %>% unique() %>% max()

        cage_change_data <- cage_change_data %>%
          mutate(ConsecActive = ifelse(Phase == "active", max_consecAct + ConsecActive, 0)) %>%
          mutate(ConsecInactive = ifelse(Phase == "inactive", max_consecInact - 1 + ConsecInactive, 0))
      }

      # Add new columns to the tibble
      cage_change_data <- cage_change_data %>%
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
      cage_change_data <- cage_change_data[c("CageChange", "Batch", "System", "AnimalID", "Sex", "Group", "Phase", "ConsecActive", "ConsecInactive", values[i])]

      # Combine with the general tibble for all cage changes and batches
      if (data_category == "total_proximity") {
        social_prox_consec <- bind_rows(social_prox_consec, cage_change_data)
      } else if (data_category == "total_movement") {
        coil_crossing_consec <- bind_rows(coil_crossing_consec, cage_change_data)
      }
    }
  }

  # Remove rows with NA values
  if (data_category == "total_proximity") {
    social_prox_consec <- na.omit(social_prox_consec)
  } else if (data_category == "total_movement") {
    coil_crossing_consec <- na.omit(coil_crossing_consec)
  }
}

# ---------------------------------------------
# CAGE AND MICE ENTROPY DATA LOADING
# ---------------------------------------------
# Load pre-saved tables for cage and mice entropy into tibbles.
# These tibbles contain consecutive phase data for all batches and cage changes.

# Load cage entropy data
cage_entropy_consec <- as_tibble(
  read_delim(
    paste0(saving_directory, tables_directory, "/all_batches", "_all_cageChanges", "_consec_cage_entropy.csv"),
    delim = ",",
    show_col_types = FALSE
  )
)

# Load mice (animal) entropy data
animal_entropy_consec <- as_tibble(
  read_delim(
    paste0(saving_directory, tables_directory, "/all_batches", "_all_cageChanges", "_consec_mice_entropy.csv"),
    delim = ",",
    show_col_types = FALSE
  )
)

# Combine all consecutive phase tibbles into a list for streamlined processing
# The list includes:
# - Social proximity data
# - Coil crossing data
# - Mice entropy data
# - Cage entropy data
consec_tibbles <- list(
  social_prox_consec,   # Social proximity data
  coil_crossing_consec, # Coil crossing data
  animal_entropy_consec, # Mice entropy data
  cage_entropy_consec   # Cage entropy data
)

# ------------------------------------------
# LINEAR MIXED-EFFECTS MODELS (LME) AND PLOTS FOR CONSECUTIVE TIBBLES
# ------------------------------------------
message("## Generating consecutive phase mean plots and performing statistical analysis ##")

# Define data categories for dynamic processing
values <- c(
  "SocialProx",      # Social proximity data
  "CoilCrossing",    # Coil crossing data
  "MiceEntropy",     # Individual animal entropy data
  "CageEntropy"      # Cage-level entropy data
)

# Define experimental phases and their corresponding consecutive phase columns
phases <- c("active", "inactive")  # Experimental phases (active and inactive)
consec_phases <- c("ConsecActive", "ConsecInactive")  # Consecutive phase columns

# Initialize a list to store generated plots for visualization
plots <- list()

# Iterate over all consecutive tibbles
for (i in seq_along(consec_tibbles)) {
  value <- values[i]  # Current data category
  color_by <- ifelse(value == "CageEntropy", "System", "Group")
  ifelse(value == "CageEntropy",
    color_palette <- c("forestgreen",
                       "goldenrod",
                       "steelblue",
                       "firebrick",
                       "darkorchid"),
    color_palette <- c("con" = "#1e3791",
                       "res" = "#8aacdb",
                       "sus" = "#f49620")
  )

  message(value)  # Log the current data category being processed

  # Iterate over active and inactive phases
  for (j in seq_along(phases)) {
    consec_data <- consec_tibbles[[i]]  # Extract the current tibble
    phase <- phases[j]  # Define the current phase (active or inactive)
    consec_phase <- consec_phases[j]  # Define the corresponding consecutive phase column

    message(phase)  # Log the current phase
    message(consec_phase)  # Log the current consecutive phase column

    # Filter the data for the current phase
    consec_data <- consec_data %>%
      filter(Phase == phase)

    ## Generate Plot ##
    # Create a ggplot object for visualizing the data
    p <- ggplot(data = consec_data,
                aes(x = fct_inorder(as_factor(.data[[consec_phase]])),
                    y = .data[[value]],
                    color = !!sym(color_by),
                    group = !!sym(color_by))) +
      stat_summary(aes(group = !!sym(color_by)),
                   fun = median,
                   geom = "line",
                   linewidth = 1) +  # Add a line representing the median
      stat_summary(aes(fill = !!sym(color_by)),
                   fun.min = function(z) {quantile(z, 0.25)},
                   fun.max = function(z) {quantile(z, 0.75)},
                   fun = median,
                   geom = "ribbon",
                   alpha = 0.2, color = NA) +
      labs(title = paste(value),
           subtitle = paste(phase, "phases")) +
      scale_fill_manual(values = color_palette) +
      scale_color_manual(values = color_palette) +
      scale_y_continuous(paste(value)) +
      scale_x_discrete(name = consec_phase) +
      facet_grid(Sex ~ .) +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 14, face = "italic"),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.ticks.x = element_line(size = 0.5),  # Style x-axis ticks
        legend.position = "top",  # Position the legend at the top
        panel.background = element_blank()  # Remove background color
      )

    # Save the generated plot to the list
    if (j == 1) plots[[2 * i - 1]] <- p
    else plots[[2 * i]] <- p

    # Save the plot to a file if the save_plots flag is enabled
    if (save_plots) ggsave(filename = paste0(working_directory, "/plots", "/consec_plots/", "consec-", value, "-", consec_phase, ".svg"), width = 5, height = 5)

    ## Linear Mixed Effects Model ##
    # Dynamically create the formula for the linear mixed-effects model
    if (value == "CageEntropy") {
      formula <- as.formula(paste(value, "~ System * Sex + ", consec_phase, " + (1 | Batch)"))
    } else {
      formula <- as.formula(paste(value, "~ Group * Sex + ", consec_phase, " + (1 | AnimalID)"))
    }

    # Perform the linear mixed-effects model analysis
    message(paste("consec lme for", value, ", ", phase))
    print(formula)
    model <- lmerTest::lmer(formula, data = consec_data)
    print(summary(model))

    # Save the LME results to a file if the save_statistics flag is enabled
    mixed_model_results_file <- paste0(working_directory, "/lme/consec-mixed_model_results/", paste0("lme-", value, "-", consec_phase, ".txt"))
    model_summary <- capture.output(summary(model))
    if (save_statistics) writeLines(model_summary, con = mixed_model_results_file)
  }
}

## Print Generated Plots ##
# Print all plots individually for visualization
for (i in seq_along(plots)) {
  print(plots[[i]])
}

# Arrange and display plots in grids
gridExtra::grid.arrange(grobs = plots, ncol = 2)

# ---------------------------------------------
# Cage Change Data Initialization
# ---------------------------------------------
# Create tibbles to store summarized data for each cage change.
# These tibbles include phase-specific statistics, system rankings,
# and relevant metadata for different analysis types.

# Social Proximity data tibble
social_prox_CC <- tibble(
  CageChange = character(),   # Cage change identifier
  System = character(),       # System identifier
  AnimalID = character(),     # Animal ID
  Group = character(),        # Group classification (e.g., con, sus, res)
  Batch = character(),        # Experimental batch
  Sex = character(),          # Sex of the animal
  Sum_act = numeric(),        # Total active phase value
  Sum_inact = numeric(),      # Total inactive phase value
  Avg_act = numeric(),        # Average active phase value
  Avg_inact = numeric(),      # Average inactive phase value
  SystemRank_act = numeric(), # Rank of the system during active phase
  SystemRank_inact = numeric()# Rank of the system during inactive phase
)

# Coil Crossing data tibble
coil_crossing_CC <- tibble(
  CageChange = character(),   # Cage change identifier
  System = character(),       # System identifier
  AnimalID = character(),     # Animal ID
  Group = character(),        # Group classification (e.g., con, sus, res)
  Batch = character(),        # Experimental batch
  Sex = character(),          # Sex of the animal
  Sum_act = numeric(),        # Total active phase value
  Sum_inact = numeric(),      # Total inactive phase value
  Avg_act = numeric(),        # Average active phase value
  Avg_inact = numeric(),      # Average inactive phase value
  SystemRank_act = numeric(), # Rank of the system during active phase
  SystemRank_inact = numeric()# Rank of the system during inactive phase
)

# Mice Entropy data tibble
mice_entropy_CC <- tibble(
  CageChange = character(),   # Cage change identifier
  System = character(),       # System identifier
  AnimalID = character(),     # Animal ID
  Group = character(),        # Group classification (e.g., con, sus, res)
  Batch = character(),        # Experimental batch
  Sex = character(),          # Sex of the animal
  Sum_act = numeric(),        # Total active phase value
  Sum_inact = numeric(),      # Total inactive phase value
  Avg_act = numeric(),        # Average active phase value
  Avg_inact = numeric(),      # Average inactive phase value
  SystemRank_act = numeric(), # Rank of the system during active phase
  SystemRank_inact = numeric()# Rank of the system during inactive phase
)

# Cage Entropy data tibble
cage_entropy_CC <- tibble(
  CageChange = character(),   # Cage change identifier
  System = character(),       # System identifier
  Batch = character(),        # Experimental batch
  Sex = character(),          # Sex of the system
  Sum_act = numeric(),        # Total active phase value
  Sum_inact = numeric(),      # Total inactive phase value
  Avg_act = numeric(),        # Average active phase value
  Avg_inact = numeric(),      # Average inactive phase value
  SystemRank_act = numeric(), # Rank of the system during active phase
  SystemRank_inact = numeric()# Rank of the system during inactive phase
)

# Create a list to dynamically process cage change data for each analysis type
# The list maps analysis types (e.g., SocialProx, CoilCrossing) to their respective tibbles.
cage_change_datas <- list(
  "SocialProx" = social_prox_CC,      # Social Proximity data
  "CoilCrossing" = coil_crossing_CC,  # Coil Crossing data
  "MiceEntropy" = mice_entropy_CC,    # Mice Entropy data
  "CageEntropy" = cage_entropy_CC     # Cage Entropy data
)

# Define the analysis types for dynamic processing
values <- c("SocialProx", "CoilCrossing", "MiceEntropy", "CageEntropy")

# Iterate over all consec tibbles and their corresponding values
for (i in c(1, 2, 3, 4)) {
  consec_tibble <- consec_tibbles[[i]]
  value <- values[i]
  message(value)
  cage_change_data <- cage_change_datas[[value]]

  # Loop through each batch and cage change
  for (batch in batches) {
    for (cageChange in cageChanges) {
      message(paste(batch, cageChange))

      # Filter the consec tibble for the current batch and cage change
      filtered_consec_tibble <- consec_tibble %>%
        filter(Batch == batch) %>%
        filter(CageChange == cageChange)

      # Load the system-animal mapping CSV for the current batch and cage change
      system_animal_ids <- tibble(read_delim(paste0(working_directory, "/tables/", batch, "_", cageChange, "_animal_ids.csv"), delim = ",", show_col_types = FALSE))

      # Iterate over each system in the system-animal mapping
      for (system in colnames(system_animal_ids)) {

        # Extract animal IDs for the current system
        animal_ids <- system_animal_ids %>%
          select(all_of(system)) %>%
          pull()

        # Skip the system if it contains incomplete data (NA values)
        if (NA %in% animal_ids) {
          cat("Skipping system", system, "in", cageChange, "of batch", batch, "\n")
          next
        }

        # Initialize vectors for rank computation
        act_system_vec <- c()  # Active phase values for the system
        inact_system_vec <- c()  # Inactive phase values for the system

        # Process CageEntropy separately from other values
        if (value == "CageEntropy") {
          sex <- ifelse(batch %in% c("B3", "B4", "B6"), "female", "male")

          # Filter active and inactive phases for the current system
          act_tibble <- filtered_consec_tibble %>%
            filter(Phase == "active") %>%
            filter(System == system) %>%
            select(all_of(value))

          inact_tibble <- filtered_consec_tibble %>%
            filter(Phase == "inactive") %>%
            filter(System == system) %>%
            select(all_of(value))

          # Compute sum and average for active and inactive phases
          sum_act <- sum(act_tibble)
          sum_inact <- sum(inact_tibble)
          avg_act <- mean(unlist(act_tibble))
          avg_inact <- mean(unlist(inact_tibble))

          # Add a new row to the cage change data tibble
          cage_change_data <- cage_change_data %>% 
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

          # Append total hours to vectors for rank computation
          act_system_vec <- append(act_system_vec, sum_act)
          inact_system_vec <- append(inact_system_vec, sum_inact)

        } else {
          # Process other values (e.g., SocialProx, CoilCrossing, etc.)
          for (animal_id in animal_ids) {

            # Determine group and sex for the current animal
            group <- ifelse(animal_id %in% sus_animals, "sus", ifelse(animal_id %in% con_animals, "con", "res"))
            sex <- ifelse(batch %in% c("B3", "B4", "B6"), "female", "male")

            # Filter active and inactive phases for the current animal
            act_tibble <- filtered_consec_tibble %>%
              filter(Phase == "active") %>%
              filter(AnimalID == animal_id) %>%
              select(all_of(value))

            inact_tibble <- filtered_consec_tibble %>%
              filter(Phase == "inactive") %>%
              filter(AnimalID == animal_id) %>%
              select(all_of(value))

            # Compute sum and average for active and inactive phases
            sum_act <- sum(act_tibble)
            sum_inact <- sum(inact_tibble)
            avg_act <- mean(unlist(act_tibble))
            avg_inact <- mean(unlist(inact_tibble))
            
            # Add a new row to the cage change data tibble
            cage_change_data <- cage_change_data %>%
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

            # Append total hours to vectors for rank computation
            act_system_vec <- append(act_system_vec, sum_act)
            inact_system_vec <- append(inact_system_vec, sum_inact)
          }
        }
        # Compute ranks for active and inactive phases
        cage_change_data <- compute_rank(cage_change_data, act_system_vec, system, cageChange, batch, "Sum_act", "SystemRank_act")
        cage_change_data <- compute_rank(cage_change_data, inact_system_vec, system, cageChange, batch, "Sum_inact", "SystemRank_inact")
      }
    }
  }
  # Update the cage change data list with the processed tibble
  cage_change_datas[[value]] <- cage_change_data
}

 # color_by <- ifelse(value == "CageEntropy", "System", "Group")
 # ifelse(value == "CageEntropy", color_palette <- c("forestgreen", "goldenrod", "steelblue", "firebrick", "darkorchid"), color_palette <- c("con" = "#1e3791", "res" = "#8aacdb", "sus" = "#f49620"))

#----------------------------------------------
# STATISTICS AND PLOTS FOR CAGE CHANGE (CC) TIBBLES
#----------------------------------------------

# Dynamic generation of plots and statistical analysis for all CC tibbles

# Definitions for dynamic processing
values <- c("SocialProx", "CoilCrossing", "MiceEntropy", "CageEntropy")  # Data categories
phases <- c("active", "inactive")  # Experimental phases

# Initialize lists to store results
avg_plots <- list()  # Average plots for each phase
rank_plots <- list()  # Rank plots for each phase
allTestResults <- list()  # Statistical test results
allPlots <- list()  # Combined plots
allPosthocResults <- list()  # Post-hoc test results

# Iterate over each data category (e.g., SocialProx, CoilCrossing, etc.)
for (i in seq_along(cage_change_datas)) {
  value <- values[i]  # Current data category
  color_by <- ifelse(value == "CageEntropy", "System", "Group")
  ifelse(value == "CageEntropy",
         color_palette <- c("forestgreen",
                            "goldenrod",
                            "steelblue",
                            "firebrick",
                            "darkorchid"),
         color_palette <- c("con" = "#1e3791",
                            "res" = "#8aacdb",
                            "sus" = "#f49620"))

  message(value)  # Log the current data category being processed

  # Iterate over active and inactive phases
  for (j in seq_along(phases)) {
    cc_data <- cage_change_datas[[value]]  # Extract the current CC tibble
    phase <- phases[j]  # Define the current phase (active or inactive)
    message(phase)  # Log the current phase

    ## PLOTS ##
    ## Average Plot ##
    avg_phase <- ifelse(phase == "active", "Avg_act", "Avg_inact")

    # Generate average plot with median and interquartile range
    avg_cc_p <- ggplot(data = cc_data, aes(x = CageChange,
                                           y = !!sym(avg_phase),
                                           color = !!sym(color_by),
                                           group = !!sym(color_by))) +
      stat_summary(aes(group = !!sym(color_by)),
                   fun = median, 
                   geom = "line",
                   linewidth = 1) +  # Median line
      stat_summary(aes(fill = !!sym(color_by)),
                   fun.min = function(z) {quantile(z, 0.25)},  # Lower quartile
                   fun.max = function(z) {quantile(z, 0.75)},  # Upper quartile
                   fun = median,
                   geom = "ribbon",
                   alpha = 0.2, color = NA) +  # Shaded ribbon for IQR
      labs(title = paste(value, "Cage Change Plot"),
           subtitle = paste(phase, "phases")) +
      scale_fill_manual(values = color_palette) +
      scale_color_manual(values = color_palette) +
      scale_y_continuous(paste(value)) +
      scale_x_discrete(name = "Cage Change") +
      facet_grid(Sex ~ .) +
      theme_minimal(base_size = 14) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5,
                                      face = "bold",
                                      size = 18),
            plot.subtitle = element_text(hjust = 0.5,
                                         size = 14,
                                         face = "italic"),
            axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 14,
                                        face = "bold"),
            axis.ticks.x = element_line(size = 0.5),
            legend.position = "top",
            panel.background = element_blank())

    # Save plot to the list
    if (j == 1) avg_plots[[2 * i - 1]] <- avg_cc_p
    else avg_plots[[2 * i]] <- avg_cc_p

    # Save plot to file if save_plots flag is enabled
    if (save_plots) ggsave(filename = paste0(working_directory, "/plots", "/cc_plots/", "avg-cc-", value, "-", phase, ".svg"), plot = avg_cc_p, width = 5, height = 5)

    ## Linear Mixed-Effects Model (LME) for Average Values ##
    message(paste("CC LME for", value, ", ", avg_phase))

    # Dynamically create formula for LME
    if (value == "CageEntropy") {
      formula <- as.formula(paste(avg_phase, "~ System * Sex + CageChange + (1 | Batch)"))
    } else {
      formula <- as.formula(paste(avg_phase, "~ Group * Sex + CageChange + (1 | AnimalID)"))
    }
    message(formula)

    # Fit the LME model
    model <- lmerTest::lmer(formula, data = cc_data)

    # Save LME results to a file
    mixed_model_results_file <- paste0(working_directory, "/lme/cc-mixed_model_results/", paste0("avg-lme-", value, "-", phase, ".txt"))
    model_summary <- capture.output(summary(model))
    if (save_statistics) writeLines(model_summary, con = mixed_model_results_file)

    # Generate rank plots and perform statistical tests for non-CageEntropy data
    if (value != "CageEntropy") {
      ## Rank Plot ##
      rank_phase <- ifelse(phase == "active", "SystemRank_act", "SystemRank_inact")  # Define rank column

      # Summarize data by grouping across all cage changes for each individual
      avg_data <- cc_data %>%
        group_by(Batch, Sex, Group, AnimalID) %>%
        summarise(avg_rank = mean(!!sym(rank_phase)))

      # Exclude "con" group for SocialProx data (not meaningful for visualization)
      if (value == "SocialProx") avg_data <- avg_data %>% filter(Group != "con")

      # Generate rank plot
      avg_rank_p <- ggplot(data = avg_data, aes(x = Group, y = avg_rank, color = Group, shape = Batch)) +
        geom_jitter(size = 4, alpha = 0.7, width = 0.2, height = 0) +
        scale_shape_manual(values = c(3, 16, 17, 15, 5, 20)) +
        scale_x_discrete("Animal ID") +
        scale_y_discrete("Rank", limits = c("1", "2", "3", "4")) +
        labs(title = paste("Median rank of each individual"),
             subtitle = paste(value, "-", phase)) +
        scale_color_manual(values = c("con" = "#1e3791",
                                      "res" = "#8aacdb",
                                      "sus" = "#f49620")) +
        stat_summary(
          fun.min = function(z) {quantile(z, 0.25)},  # Lower quartile
          fun.max = function(z) {quantile(z, 0.75)},  # Upper quartile
          fun = median,
          color = "black",
          size = 0.8,
          shape = 16) +
        facet_grid(Sex ~ .) +
        theme_minimal(base_size = 14) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
              plot.subtitle = element_text(hjust = 0.5, size = 14, face = "italic"),
              axis.text.x = element_text(size = 12),
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 14, face = "bold"),
              axis.ticks.x = element_line(size = 0.5),
              legend.position = "top",
              panel.background = element_blank())

      # Save rank plot to the list
      if (j == 1) rank_plots[[2 * i - 1]] <- avg_rank_p
      else rank_plots[[2 * i]] <- avg_rank_p

      # Save rank plot to file if save_plots flag is enabled
      if (save_plots) ggsave(filename = paste0(working_directory, "/plots", "/avg_rank_plots/", "avg_ranks-", value, "-", phase, ".svg"),
                             plot = avg_rank_p,
                             width = 5, height = 5)

      ## Statistical Analysis on Rank Averages ##
      for (sex in c("female", "male")) {
        message("Statistics on rank averages")
        print(sex)

        # Perform statistical tests and generate plots
        result <- testAndPlotVariable(avg_data, value, "avg_rank", phase, sex)

        # Append results to respective lists
        if (!is.null(result)) {
          if (is.null(result$posthocResults)) {
            allTestResults <- c(allTestResults, list(result$testResults))
          } else {
            allTestResults <- c(allTestResults, list(result$testResults))
            allPosthocResults <- c(allPosthocResults, list(result$posthocResults))
          }
          allPlots <- c(allPlots, list(result$plot))
        }
      }
    }
  }
}

# ----------------------------------------------
# SAVE PLOTS AND STATISTICS
# ----------------------------------------------
# Save all plots to a file if the save_plots flag is enabled

# Display all average plots individually for inspection
for (i in seq_along(avg_plots)) {
  print(avg_plots[[i]])
}

# Arrange and display average plots in a grid layout
gridExtra::grid.arrange(grobs = avg_plots, ncol = 2)   # Grid of average plots

# Arrange and display rank plots in a grid layout
gridExtra::grid.arrange(grobs = rank_plots, ncol = 2)  # Grid of rank plots

# Arrange and display all combined plots in a grid layout
gridExtra::grid.arrange(grobs = allPlots, ncol = 4)    # Grid of all plots

# ----------------------------------------------
# SAVE STATISTICAL RESULTS
# ----------------------------------------------
# Combine all statistical test results into a single data frame
# Ensure the statistics directory exists, creating it if necessary

# Combine all statistical test results into a single data frame
allTestResultsDf <- bind_rows(allTestResults)

# Save the combined test results to a CSV file
# Ensure the statistics directory exists, creating it if necessary
statistics_directory <- paste0(working_directory, "/statistics")
if (!dir.exists(statistics_directory)) {
  dir.create(statistics_directory, recursive = TRUE)
}
write.csv(allTestResultsDf, file = paste0(working_directory, "/statistics/test_results.csv"), row.names = FALSE)

# Save post hoc test results to a CSV file, if available
if (!is.null(allPosthocResults) && length(allPosthocResults) > 0) {
  allPosthocResultsDf <- bind_rows(allPosthocResults)
  write.csv(allPosthocResultsDf, file = paste0(working_directory, "/statistics/posthoc_results.csv"), row.names = FALSE)
}