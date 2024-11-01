## 12/2023
## Anja Magister
## ANALYSIS OF ANIMAL POSITIONS - PREPROCESSING ##

## NEEDED FILE STRUCTURE IN WORKING DIRECTORY
## - the data to be preprocessed must be in a folder called "raw_data"
## - there must be a folder called "preprocessed_data" so that this code can save its result there

## Load required packages
required_packages <- c("readr", "stringr", "dplyr", "lubridate", "tibble", "writexl", "foreach", "doParallel", "parallel", "tidyverse")

# Install required packages if not already installed
for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  } 
  library(package, character.only = TRUE)
}

# Set paths
setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability")
source('C:/Users/topohl/Documents/GitHub/MMMSociability/E9_SIS_AnimalPos-functions.R')

#working_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability"
#source_dir <- "C:/Users/topohl/Documents/GitHub/MMMSociability/"
output_dir <- file.path(getwd(), "preprocessed_data")

# Define batch and cage change
batches <- c("B1", "B2", "B3", "B4", "B5", "B6")
cageChanges <- c("CC1", "CC2", "CC3", "CC4")

# Read excluded animals
excl_animals <- readLines(file.path(getwd(), "raw_data", "excluded_animals.csv"))

for (batch in batches) {
  for (cageChange in cageChanges) {
    print(paste(batch, cageChange))
    
    # Current csv filename
    filename <- paste0("E9_SIS_", batch, "_", cageChange, "_AnimalPos")
    
    # Path of csv file 
    csvFilePath <-  paste0(getwd(), "/raw_data", "/", batch, "/", filename, ".csv")
    
    # Load functions
    # source(paste0(source_dir, "E9_SIS_AnimalPos-functions.R"))
    
    # Read csv file into tibble
    data <- as_tibble(read_delim(csvFilePath, delim = ";", show_col_types = FALSE))
    
    ####################################################################################################################################
    ###### PREPROCESSING OF DATA: ######

    # The columns 'RFID', 'AM', and 'zPos' are not required for the analysis and are therefore removed from the dataset.
    data <- select(data, -c(RFID, AM, zPos))

    # Convert 'DateTime' Column to POSIXct Format
    # The 'DateTime' column is converted to a POSIXct datetime object to facilitate time-based operations.
    # The format of the datetime is specified as "%d.%m.%Y %H:%M:%S".
    # The timezone is set to 'UTC' to standardize the data. Note that the actual experiment may be in UTC+1, but using UTC ensures
    # consistency, particularly when exporting the data to CSV.
    data$DateTime <- as.POSIXct(data$DateTime, format = "%d.%m.%Y %H:%M:%S", tz = "UTC")

    # Separate 'Animal' Column into 'AnimalID' and 'System'
    # The 'Animal' column, which contains a combination of the animal's ID and the system identifier, is split into two separate
    # columns: 'AnimalID' and 'System'. The split is based on either an underscore ('_') or hyphen ('-') delimiter.
    data[c('AnimalID', 'System')] <- str_split_fixed(data$Animal, '[_-]', 2)

    ###### Convert xPos and yPos into a Single Column Named "PositionID" ######
    message("Merging xPos and yPos into PositionID...")

    # Create a Position Mapping Table
    # A reference table ('tblPosition') is created that maps each unique combination of xPos and yPos coordinates to a unique PositionID.
    # This step ensures that each spatial position in the system is uniquely identified by an integer ID, simplifying downstream analyses.
    positions <- select(data, c(xPos, yPos))
    unique_positions <- unique(positions)
    tblPosition <- tibble(PositionID = c(1:8), 
                         xPos = c(0, 100, 200, 300, 0, 100, 200, 300), 
                         yPos = c(0, 0, 0, 0, 116, 116, 116, 116))

    # Assign PositionID to Each Row in the Dataset
    # The dataset is augmented with a new column 'PositionID', which replaces the original xPos and yPos coordinates.
    # The 'find_id' function maps each (xPos, yPos) pair to its corresponding PositionID from 'tblPosition'.
    data_ids <- data %>% rowwise() %>% mutate(PositionID = find_id(xPos, yPos, tblPosition))

    # Reorder Columns for Consistency
    # The columns are reordered to follow a logical sequence: 'DateTime', 'AnimalID', 'System', 'PositionID'.
    data <- data_ids[c('DateTime', 'AnimalID', 'System', 'PositionID')]

    # Exclude Data from Certain Animals
    # Any rows corresponding to animals listed in 'excl_animals' are removed from the dataset.
    # This step ensures that only data from animals of interest are retained for further analysis.
    data <- data %>% filter(!AnimalID %in% excl_animals)

    # Sort Data by DateTime
    # The dataset is sorted by the 'DateTime' column to ensure that all rows are in chronological order.
    data <- data %>% arrange(., DateTime)
    
    ################################################################################################################################
    #### Add new column "Phase" ####
    
    message("adding row for phase transition markers...")
    
    # Add rows with specific time when the active/inactive phases for the mice change (6h30 and 18h30)
    # Important for next step (active/inactive phase division)
    # Without the marker rows, we would have slightly incorrect assignments of the time of the phases
    data <- addPhaseTransitionMarkers(data)
    
    #### Add column "phase" and define values for every row depending on the time in the row ####
    data <- data %>% rowwise() %>% mutate(Phase = ifelse(format(DateTime, "%H:%M", tz = "UTC") >= "18:30" | format(DateTime, "%H:%M", tz = "UTC") < "06:30", "Active", "Inactive"))
    
    ### Add new columns "consecActive" and "consecInactive" (consecutive), tells us how often the phases change and at which phase we are
    message("adding consecutive phase counters...")
    
    # Initialize empty columns
    data <- data %>% rowwise() %>% mutate(ConsecActive = 0) %>% mutate(ConsecInactive = 0)
    
    # Edit the columns so that the counting of the phases works
    # Each start of a new phase, it counts a new phase
    # Important for comparison between the phases (day per day)
    data <- consecPhases(data)
    
    ################################################################################################################################
    #### Add necessary rows for every new day ####
    
    # A start time of every mouse in every cage is needed to divide data per day
    data <- addDayMarkerRows(data)
    
    ################################################################################################################################
    # SAVING #
    
    ## Saving data in csv file
    output_dir <- file.path(getwd(), "preprocessed_data")
    
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    write_csv(data, file.path(output_dir, paste0(filename, "_preprocessed.csv")))
    
    message(paste("saving file at", file.path(output_dir, paste0(filename, "_preprocessed.csv"))))
    }
}
