## ANALYSIS OF ANIMAL POSITIONS - PREPROCESSING ##

## NEEDED FILE STRUCTURE IN WORKING DIRECTORY
## - the data to be preprocessed must be in a folder called "raw_data"
## - there must be a folder called "preprocessed_data" so that this code can save its result there

# Load required packages
requiredPackages <- c("readr", "stringr", "dplyr", "lubridate", "tibble", "writexl", "foreach", "doParallel", "parallel", "tidyverse")

# Install and load required packages
lapply(requiredPackages, function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, dependencies = TRUE)
  }
  library(package, character.only = TRUE)
})

# Set working directory and source functions
setwd("S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability")
source('C:/Users/topohl/Documents/GitHub/MMMSociability/E9_SIS_AnimalPos-functions.R')

# Define constants
outputDir <- file.path(getwd(), "preprocessed_data")
batches <- c("B1", "B2", "B3", "B4", "B5", "B6")
cageChanges <- c("CC1", "CC2", "CC3", "CC4")

# Read excluded animals
excludedAnimalsDir <- file.path(getwd(), "raw_data", "excluded_animals.csv")
if (!file.exists(excludedAnimalsDir)) {
  stop("The excluded animals file does not exist at the expected path.")
}
exclAnimals <- readLines(excludedAnimalsDir)

# Ensure output directory exists
if (!dir.exists(outputDir)) {
  dir.create(outputDir, recursive = TRUE)
}

# Setup parallel processing
numCores <- detectCores() - 1
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Process all files in parallel
foreach(batch = batches, .packages = requiredPackages) %:%
  foreach(cageChange = cageChanges) %dopar% {
    preprocess_file(batch, cageChange, exclAnimals)
  }

# Stop the cluster
stopCluster(cl)
