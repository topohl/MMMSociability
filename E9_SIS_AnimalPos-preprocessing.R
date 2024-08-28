## 12/2023
## Anja Magister
## ANALYSIS OF ANIMAL POSITIONS -PREPROCESSING ##
##
## NEEDED FILE STRUCTURE IN WORKING DIRECTORY
## - the data to be preprocessed must be in a folder called "raw_data"
## - there must be a folder called "preprocessed_data" so that this code can save its result there
##
## CUSTOMISABLE VARIABLES:
## - working directory: choose your current directory
## - batch and cage change: on which file are you working->filename will be created with batch and cc information

# libraries
library(readr)        # load readr package for reading csv files
library(stringr)      #for splitting strings in tibble
library(dplyr)
library(lubridate)    # for rounding time, time operations in general
library(tibble)       #important for tibble operations
library(writexl)      #save as csv file

# paths
working_directory <- "S:/Lab_Member/Anja/Git/MDC_Bachelor/E9_SIS_AnimalPos"
#working_directory <- "/home/anja/Dokumente/FU BERLIN/BA/Git/MDC_Bachelor/E9_SIS_AnimalPos"


# define batch and cage change
batches <- c("B1", "B2", "B3", "B4", "B5", "B6")
cageChanges <- c("CC1", "CC2", "CC3", "CC4")


#excluded animals
excl_animals <- readLines(paste0(working_directory,"/raw_data/excluded_animals.csv"))


for(batch in batches){
  
  for(cageChange in cageChanges){
    print(paste(batch, cageChange))
    
    #current csv filename
    filename <- paste0("E9_SIS_", batch, "_", cageChange, "_AnimalPos")
    
    #path of csv file 
    csvFilePath <-  paste0(working_directory,"/raw_data", "/", batch, "/", filename, ".csv")
    
    #functions
    source(paste0(working_directory,"/E9_SIS_AnimalPos-functions.R"))
    
    # read csv file in tibble
    overallData <- as_tibble(read_delim(csvFilePath,delim = ";", show_col_types = FALSE))
    
    ####################################################################################################################################
    ###### PREPROCESSING OF overallData: ######
    
    # delete unnecessary columns
    overallData <- select(overallData, -c(RFID, AM, zPos))
    # convert the DateTime column to a datetime format(also rounds the DateTime)
    # timezone is very important! here in utc, but we are in utc+1. UTC is important for write_csv...
    overallData$DateTime <- as.POSIXct(overallData$DateTime, format = "%d.%m.%Y %H:%M:%S", tz="UTC")
    
    # separate Animal into his ID an his system
    overallData[c('AnimalID', 'System')] <- str_split_fixed(overallData$Animal, '[_-]', 2)
    
    
    
    
    ###### convert xPos and yPos into one column named "PositionID" ######
    message("translate positions")
    #create Positions_tibble that contains every possible combination of our coordinates together with an ID
    positions <- select(overallData, c(xPos,yPos))
    unique_positions <- unique(positions)
    Positions_tibble <- tibble(PositionID = c(1:8), xPos = c(0,100,200,300,0,100,200,300), yPos = c(0,0,0,0,116,116,116,116))
    
    # Adding column PositionID to overallData instead of two colums with x and y coordinates
    overallData_ids <- overallData %>% rowwise() %>%
      mutate(PositionID = find_id(xPos, yPos, Positions_tibble))
    
    
    
    ##### sort columns #####
    overallData <- overallData_ids[c('DateTime', 'AnimalID', 'System', 'PositionID')]
    
    ##### delete rows of excluded animals #####
    overallData <- overallData%>%
      filter(!AnimalID %in% excl_animals)
    
    ##### sort by Date Time #####
    overallData <- overallData%>%
      arrange(., DateTime)
    
    ################################################################################################################################
    #### add new column "Phase" ####
    
    message("add PhaseMarker rows")
    # add rows with specific time when the active/inactive phases for the mice change(6h30 and 18h30)
    # important for next step (active/inactive phase division)
    # without the marker rows we would have slightly incorrect assignments of the time of the phases
    overallData <- addPhaseMarkerRows(overallData)
    
    #### add column "phase" and define values for every row depending on the time in the row ####
    overallData <- overallData %>%
      rowwise() %>%
      mutate(Phase = ifelse(format(DateTime, "%H:%M", tz="UTC") >= "18:30" | format(DateTime, "%H:%M", tz="UTC") < "06:30","Active","Inactive"))
    
    ### add new columns "consecActive" and "consecInactive"(consecutive), tells  us how often the phases change and at which phase we are
    message("add Consec Columns")
    # initialize empty columns
    overallData <- overallData %>%
      rowwise() %>%
      mutate(ConsecActive = 0)%>%
      mutate(ConsecInactive = 0)
    
    # edit the colums so that the counting of the phases works
    # each start of a new phase it counts a new phase
    # important for comparison between the phases(day per day)
    overallData <- consecPhases(overallData)
    
    
    ################################################################################################################################
    #### add necessary rows for every new day ####
    # a start time of every mouse in every cage is needed to divide data per day
    overallData <- addDayMarkerRows(overallData)
    
    ################################################################################################################################
    # SAVING #
    
    ## saving data in csv file
    write_csv(overallData, paste0(working_directory,"/preprocessed_data/", filename, "_preprocessed.csv"))
  }
}
