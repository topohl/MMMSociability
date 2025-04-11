# ANALYSIS OF ANIMAL POSITIONS - PREPROCESSING
#
# The script performs the following steps:
# 1. Loads the required R packages.
# 2. Installs any missing packages.
# 3. Sets the working directory and sources additional functions from an external R script.
# 4. Defines constants for the output directory, batches, and cage changes.
# 5. Reads a list of excluded animals from a CSV file.
# 6. Ensures the output directory exists, creating it if necessary.
# 7. Sets up parallel processing using available CPU cores minus one.
# 8. Processes all files in parallel for each batch and cage change, excluding specified animals.
# 9. Stops the parallel processing cluster after completion.
#
# Note:
# - Ensure that the "raw_data" folder contains the necessary data files.
# - The "excluded_animals.csv" file should be located in the "raw_data" folder.
# - The external functions script should be located at the specified path.
# - Adjust the working directory path as needed.
# - Adjust the path to the excluded animals file as needed.

requiredPackages <- c("readr", "stringr", "dplyr", "lubridate", "tibble", "writexl", "foreach", "doParallel", "parallel", "tidyverse")

# Load required packages using pacman
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", dependencies = TRUE)
}
pacman::p_load(readr, stringr, dplyr, lubridate, tibble, writexl, foreach, doParallel, parallel, tidyverse, tcltk)

# Set working directory and source functions
setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability")
source('C:/Users/topohl/Documents/GitHub/MMMSociability/E9_SIS_AnimalPos-functions.R')

# Check if preprocessed_data folder exists, if not create it
if (!dir.exists("preprocessed_data")) {
  dir.create("preprocessed_data")
}

# Define constants for output directory, batches, and cage changes
outputDir <- file.path(getwd(), "preprocessed_data")
batches <- c("B1", "B2", "B3", "B4", "B5", "B6")
cageChanges <- c("CC1", "CC2", "CC3", "CC4")

# Read excluded animals from CSV file
excludedAnimalsDir <- file.path(getwd(), "raw_data", "excluded_animals.csv")
if (!file.exists(excludedAnimalsDir)) {
  stop("The excluded animals file does not exist at the expected path.")
}
exclAnimals <- readLines(excludedAnimalsDir)

# Ensure the output directory exists, create if it doesn't
if (!dir.exists(outputDir)) {
  dir.create(outputDir, recursive = TRUE)
}

# Setup parallel processing using available cores minus one
numCores <- detectCores() - 1
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Process all files in parallel for each batch and cage change
foreach(batch = batches, .packages = requiredPackages) %:%
  foreach(cageChange = cageChanges) %dopar% {
    preprocess_file(batch, cageChange, exclAnimals)
  }

# Stop the parallel processing cluster
stopCluster(cl)
