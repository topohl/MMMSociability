###############################################
# Data Loading and Formatting Script
###############################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

# -------------------------------------------------
# Paths
# -------------------------------------------------
#base_dir <- "D:/MMMSociability/tables/noHomeCage_test"
#raw_data_dir <- "D:/MMMSociability/raw_data"
#output_dir <- "D:/MMMSociability/processed_data/data_lme_format"

base_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability/tables/noHomeCage"
raw_data_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability/raw_data"
output_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability/processed_data/data_lme_format"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -------------------------------------------------
# Load group assignments
# -------------------------------------------------
con_animals <- read_csv(file.path(raw_data_dir, "con_animals.csv"), col_names = FALSE)
sus_animals <- read_csv(file.path(raw_data_dir, "sus_animals.csv"), col_names = FALSE)

con_list <- con_animals[[1]]
sus_list <- sus_animals[[1]]

# Function to assign group
assign_group <- function(animal_id) {
  if (animal_id %in% con_list) return("CON")
  if (animal_id %in% sus_list) return("SUS")
  return("RES")
}

# -------------------------------------------------
# Metadata mappings
# -------------------------------------------------
batch_sex_map <- c(
  "B1" = "m", "B2" = "m", "B5" = "m",
  "B3" = "f", "B4" = "f", "B6" = "f"
)

# Function to determine phase from HalfHour
assign_phase <- function(halfhour_name) {
  h_num <- as.integer(sub("H", "", halfhour_name))
  block <- floor(h_num / 24)
  ifelse(block %% 2 == 0, "Active", "Inactive")
}

# Calculate HalfHourElapsed
calc_halfhour_elapsed <- function(halfhour_name) {
  as.integer(sub("H", "", halfhour_name))
}

# Calculate TH (time in hours, half-hour increments)
calc_th <- function(halfhour_elapsed) {
  halfhour_elapsed * 0.5
}

# -------------------------------------------------
# Main processing function
# -------------------------------------------------
process_batch_change <- function(batch, change) {
  cat(sprintf("Processing %s - %s...\n", batch, change))

  # Construct file paths with by_halfhour subdirectory
  data_dir <- file.path(base_dir, batch, change, "by_halfhour")
  movement_file <- file.path(data_dir, paste0(batch, "_", change, "_halfhour_movement.csv"))
  proximity_file <- file.path(data_dir, paste0(batch, "_", change, "_halfhour_proximity.csv"))

  # Check if files exist
  if (!file.exists(movement_file)) {
    cat(sprintf("  Warning: Movement file not found: %s\n", movement_file))
    return(NULL)
  }
  if (!file.exists(proximity_file)) {
    cat(sprintf("  Warning: Proximity file not found: %s\n", proximity_file))
    return(NULL)
  }

  # Read data
  movement_df <- read_csv(movement_file, show_col_types = FALSE)
  proximity_df <- read_csv(proximity_file, show_col_types = FALSE)

  # Get sex for this batch
  sex <- batch_sex_map[batch]

  # Process movement data - remove system columns
  movement_df <- movement_df %>%
    select(-matches("^sys\\."))

  # Convert to long format
  movement_long <- movement_df %>%
    pivot_longer(
      cols = -HalfHour,
      names_to = "AnimalNum",
      values_to = "Movement"
    ) %>%
    mutate(
      HalfHourElapsed = calc_halfhour_elapsed(HalfHour),
      Phase = assign_phase(HalfHour)
    )

  # Process proximity data
  proximity_long <- proximity_df %>%
    pivot_longer(
      cols = -HalfHour,
      names_to = "AnimalNum",
      values_to = "Proximity"
    ) %>%
    mutate(
      HalfHourElapsed = calc_halfhour_elapsed(HalfHour)
    )

  # Merge movement and proximity
  combined_df <- movement_long %>%
    left_join(
      proximity_long %>% select(HalfHour, AnimalNum, Proximity),
      by = c("HalfHour", "AnimalNum")
    )

  # Add metadata and calculate derived variables
  combined_df <- combined_df %>%
    mutate(
      Batch = batch,
      Change = change,
      Sex = sex,
      Group = sapply(AnimalNum, assign_group),
      TH = calc_th(HalfHourElapsed)
    ) %>%
    select(
      Batch, Change, Sex, Phase, Group, AnimalNum,
      HalfHour, HalfHourElapsed, TH,
      Movement, Proximity
    ) %>%
    arrange(AnimalNum, HalfHourElapsed)

  cat(sprintf("  Processed %d observations\n", nrow(combined_df)))
  return(combined_df)
}

# -------------------------------------------------
# Process all batches and changes
# -------------------------------------------------
batches <- c("B1", "B2", "B3", "B4", "B5", "B6")
changes <- c("CC1", "CC2", "CC3", "CC4")

all_data <- list()

for (batch in batches) {
  for (change in changes) {
    result <- process_batch_change(batch, change)
    if (!is.null(result)) {
      all_data[[paste(batch, change, sep = "_")]] <- result
    }
  }
}

# Combine all data
if (length(all_data) > 0) {
  data_filtered_agg <- bind_rows(all_data)

  # Convert to factors with appropriate levels
  data_filtered_agg <- data_filtered_agg %>%
    mutate(
      Batch = factor(Batch),
      Change = factor(Change),
      Sex = factor(Sex, levels = c("f", "m")),
      Phase = factor(Phase, levels = c("Active", "Inactive")),
      Group = factor(Group, levels = c("CON", "RES", "SUS")),
      AnimalNum = factor(AnimalNum),
      HalfHourElapsed = as.integer(HalfHourElapsed)
    )

  # remove all inactive and active phases > 2 from CC4
  data_filtered_agg <- data_filtered_agg %>%
    filter(!(Change == "CC4" & 
               ((Phase == "Active" & HalfHourElapsed >= 96) | 
                  (Phase == "Inactive" & HalfHourElapsed >= 72))))

  # Summary statistics
  cat("\n=== Data Summary ===\n")
  cat(sprintf("Total observations: %d\n", nrow(data_filtered_agg)))
  cat(sprintf("Unique animals: %d\n", length(unique(data_filtered_agg$AnimalNum))))
  cat(sprintf("Batches: %s\n", paste(levels(data_filtered_agg$Batch), collapse = ", ")))
  cat(sprintf("Changes: %s\n", paste(levels(data_filtered_agg$Change), collapse = ", ")))
  cat("\nGroup distribution:\n")
  print(table(data_filtered_agg$Group))
  cat("\nSex distribution:\n")
  print(table(data_filtered_agg$Sex))
  cat("\nPhase distribution:\n")
  print(table(data_filtered_agg$Phase))

  # Check for missing values
  cat("\n=== Missing Value Summary ===\n")
  missing_summary <- data_filtered_agg %>%
    summarise(
      Movement_NA = sum(is.na(Movement)),
      Proximity_NA = sum(is.na(Proximity))
    )
  print(missing_summary)

  # Save processed data
  output_file <- file.path(output_dir, "data_filtered_agg.csv")
  write_csv(data_filtered_agg, output_file)
  cat(sprintf("\nData saved to: %s\n", output_file))

  # Also save as RDS for faster loading in R
  output_rds <- file.path(output_dir, "data_filtered_agg.rds")
  saveRDS(data_filtered_agg, output_rds)
  cat(sprintf("Data also saved as RDS: %s\n", output_rds))

} else {
  cat("No data files found. Please check your directory structure.\n")
}

cat("\n=== Processing Complete ===\n")
