if (!require("pacman")) install.packages("pacman")
pacman::p_load(readr, dplyr, lubridate, tibble, tidyr, purrr, ggplot2, reshape2, scales, stringr)

show_plots <- FALSE
save_plots <- FALSE
save_tables <- FALSE
exclude_homecage <- FALSE

working_directory <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability"
saving_directory <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability"
plots_directory <- "/plots/noHomeCage"
tables_directory <- "/tables/noHomeCage"

cat("Checking rep function identity:\n")
print(rep)
cat("Checking class(rep): ", class(rep), "\n")
cat("Checking identical rep and base::rep: ", identical(rep, base::rep), "\n")
cat("Testing static rep outside any block:\n")
test1 <- rep(NA, times=4)
cat("test1:"); print(test1)
test2 <- rep(NA, 4)
cat("test2:"); print(test2)
cat("\n--- End pre-checks ---\n")

source(paste0(working_directory, "/E9_SIS_AnimalPos-functions.R"))

batches <- c("B1", "B2", "B3", "B4", "B5", "B6")
cageChanges <- c("CC1", "CC2", "CC3", "CC4")

for (batch in batches) {
  for (cageChange in cageChanges) {
    cat("\n=== Batch:", batch, " CageChange:", cageChange, "===\n")
    filename <- paste0("E9_SIS_", batch, "_", cageChange, "_AnimalPos")
    csvFilePath <- paste0(working_directory, "/preprocessed_data/", filename, "_preprocessed.csv")
    cat("Reading:", csvFilePath, "\n")
    data_preprocessed <- as_tibble(read_delim(csvFilePath, delim = ",", show_col_types = FALSE))
    cat("Read rows:", nrow(data_preprocessed), "\n")

    data_preprocessed$DateTime <- as.POSIXct(data_preprocessed$DateTime, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

    # Floor DateTime to 30-min bins
    data_preprocessed$Bin30min <- lubridate::floor_date(data_preprocessed$DateTime, unit = "30 minutes")

    unique_systems <- str_sort(unique(data_preprocessed$System))
    cat("Unique systems:", unique_systems, "\n")
    phases <- c("Active", "Inactive")
    active_phases_number <- unique(data_preprocessed$ConsecActive)
    inactive_phases_number <- unique(data_preprocessed$ConsecInactive)
    active_phases_number <- active_phases_number[active_phases_number != 0]
    inactive_phases_number <- inactive_phases_number[!inactive_phases_number %in% c(0, 1, max(inactive_phases_number))]
    cat("Used phase numbers: active:", active_phases_number, " inactive:", inactive_phases_number, "\n")
    unique_animals <- unique(data_preprocessed$AnimalID)
    cat("Unique animals:", unique_animals, "\n")

    cat("Forward-filling positions for last known location carry-forward...\n")
    data_preprocessed <- data_preprocessed %>%
      arrange(AnimalID, DateTime) %>%
      group_by(AnimalID) %>%
      tidyr::fill(PositionID, .direction = "down") %>%
      tidyr::fill(System, .direction = "down") %>%
      tidyr::fill(Phase, .direction = "down") %>%
      tidyr::fill(ConsecActive, .direction = "down") %>%
      tidyr::fill(ConsecInactive, .direction = "down") %>%
      ungroup()

    # Create bin x animal expanded grid for consistent coverage
    all_bins <- sort(unique(data_preprocessed$Bin30min))
    expand_bin_animal <- expand.grid(Bin30min = all_bins, AnimalID = unique_animals, stringsAsFactors = FALSE)

    # Merge forward-filled data into the expanded grid
    filled_bins <- left_join(expand_bin_animal, data_preprocessed, by = c("Bin30min", "AnimalID"))
    cat("filled_bins rows:", nrow(filled_bins), " Preview:\n")
    print(head(filled_bins, 5))

    # Initialize results tables
    phases_column <- c()
    for (i in 1:max(length(inactive_phases_number), length(active_phases_number))) {
      if (i <= length(active_phases_number)) phases_column <- c(phases_column, paste0("A", i))
      if (i <= length(inactive_phases_number)) phases_column <- c(phases_column, paste0("I", i + 1))
    }

    result_total_proximity <- tibble("Phase" = phases_column)
    for (animal in unique_animals) { result_total_proximity[[animal]] <- NA }
    result_total_movement <- tibble("Phase" = phases_column)
    for (animal in unique_animals) { result_total_movement[[animal]] <- NA }
    for (system in unique_systems) { result_total_movement[[system]] <- NA }
    result_total_positions <- tibble("Phase" = phases_column)
    for (system in unique_systems) { result_total_positions[[system]] <- NA }

    # Initialize system animal IDs
    system_animal_ids <- tibble(
      "sys.1" = rep(NA, times=4),
      "sys.2" = rep(NA, times=4),
      "sys.3" = rep(NA, times=4),
      "sys.4" = rep(NA, times=4),
      "sys.5" = rep(NA, times=4)
    )

    # Create a tibble with all animals and bins
all_bins <- sort(unique(data_preprocessed$Bin30min))
expand_bin_animal <- expand.grid(Bin30min = all_bins, AnimalID = unique_animals, stringsAsFactors = FALSE)

# Left join the data with this full animal-bin grid to have rows for all animal-bin combinations
# Forward fill positions to fill missing data (i.e. position when animal did not move in that bin)
filled_bins <- expand_bin_animal %>%
  left_join(data_preprocessed, by = c("Bin30min", "AnimalID")) %>%
  arrange(AnimalID, Bin30min, DateTime) %>%
  group_by(AnimalID) %>%
  tidyr::fill(PositionID, .direction = "down") %>%
  tidyr::fill(System, .direction = "down") %>%
  tidyr::fill(Phase, .direction = "down") %>%
  tidyr::fill(ConsecActive, .direction = "down") %>%
  tidyr::fill(ConsecInactive, .direction = "down") %>%
  ungroup()


    for (system_id in unique_systems) {
      cat("\n--- SYSTEM:", system_id, "---\n")
      data_system <- filled_bins %>% filter(System == system_id) %>% as_tibble()
      animal_ids <- unique(data_system$AnimalID)
      system_complete <- ifelse(length(animal_ids) < 4, FALSE, TRUE)
      while (length(animal_ids) < 4) { animal_ids <- append(animal_ids, NA) }
      system_animal_ids[[system_id]] <- animal_ids

      systemHeatmaps_proximity <- list()
      systemHeatmaps_positions <- list()

      for (phase in phases) {
        phase_numbers <- ifelse(phase == "Active", active_phases_number, inactive_phases_number)
        for (phase_number in phase_numbers) {
          data_system_phase <- data_system %>%
            filter(ConsecActive == ifelse(phase == "Active", phase_number, 0)) %>%
            filter(ConsecInactive == ifelse(phase == "Inactive", phase_number, 0)) %>%
            as_tibble()

          animal_list <- list(
            "animal_1" = list(name = "", time = "", position = 0),
            "animal_2" = list(name = "", time = "", position = 0),
            "animal_3" = list(name = "", time = "", position = 0),
            "animal_4" = list(name = "", time = "", position = 0),
            "data_temp" = list(elapsed_seconds = 0, current_row = 0)
          )
          count_proximity_list <- list(m1 = c(0,0,0,0), m2 = c(0,0,0,0), m3 = c(0,0,0,0), m4 = c(0,0,0,0))
          count_position_list <- list(c(1,0),c(2,0),c(3,0),c(4,0),c(5,0),c(6,0),c(7,0),c(8,0))
          total_proximity_list <- list(c(animal_ids[1],0),c(animal_ids[2],0),c(animal_ids[3],0),c(animal_ids[4],0))
          count_movement_list <- list(c(animal_ids[1],0),c(animal_ids[2],0),c(animal_ids[3],0),c(animal_ids[4],0),c(system_id,0))

          animal_list <- initialize_animal_positions(animal_ids, data_system_phase, animal_list)

          initial_time <- animal_list[[1]][[2]]
          current_row <- 5
          total_rows <- nrow(data_system_phase) + 1

          while(current_row < total_rows){
            previous_animal_list <- animal_list
            animal_list <- update_animal_list(animal_ids, animal_list, data_system_phase, initial_time, current_row)
            elapsed_seconds <- animal_list[["data_temp"]][["elapsed_seconds"]]
            if(system_complete){
              count_proximity_list <- update_proximity(previous_animal_list, animal_list, count_proximity_list, elapsed_seconds)
              total_proximity_list <- update_total_proximity(previous_animal_list, animal_list, total_proximity_list, elapsed_seconds)
            }
            count_position_list <- update_position(previous_animal_list, animal_list, count_position_list, elapsed_seconds)
            count_movement_list <- update_movement(previous_animal_list, animal_list, count_movement_list, elapsed_seconds)
            current_row <- animal_list[["data_temp"]][["current_row"]]
            initial_time <- animal_list[[1]][[2]]
          }

          if(system_complete){
            heatmap_closeness <- generateHeatMapProximity(count_proximity_list, batch, cageChange, system_id, animal_ids, phase, phase_number)
            systemHeatmaps_proximity <- c(systemHeatmaps_proximity, list(heatmap_closeness))
          }
          heatmap_positions <- generateHeatMapPositions(count_position_list, batch, cageChange, system_id, phase, phase_number)
          systemHeatmaps_positions <- c(systemHeatmaps_positions, list(heatmap_positions))

          if(system_complete){
            p <- paste0(substr(phase,1,1), phase_number)
            row <- which(result_total_proximity$Phase == p)
            for(i in 1:4){
              animal <- total_proximity_list[[i]][[1]]
              result_total_proximity[[animal]][[row]] <- total_proximity_list[[i]][[2]]
            }
            for(i in 1:4){
              animal <- count_movement_list[[i]][[1]]
              if(!is.na(animal)){
                result_total_movement[[animal]][[row]] <- count_movement_list[[i]][[2]]
              }
            }
            result_total_movement[[system_id]][[row]] <- count_movement_list[[5]][[2]]
          }

          unique_bins <- sort(unique(data_system_phase$Bin30min))
          for(bin_time in unique_bins){
            bin_data <- data_system_phase %>% filter(Bin30min == bin_time)
            if(nrow(bin_data) == 0){
              bin_data <- filled_bins %>% filter(Bin30min == bin_time, System == system_id)
            }
            if(nrow(bin_data) == 0) next

            animal_list_bin <- list(
              "animal_1"=list(name="", time="", position=0),
              "animal_2"=list(name="", time="", position=0),
              "animal_3"=list(name="", time="", position=0),
              "animal_4"=list(name="", time="", position=0),
              "data_temp"=list(elapsed_seconds=0, current_row=0)
            )
            count_prox_bin <- list(m1=c(0,0,0,0), m2=c(0,0,0,0), m3=c(0,0,0,0), m4=c(0,0,0,0))
            total_prox_bin <- list(c(animal_ids[1],0), c(animal_ids[2],0), c(animal_ids[3],0), c(animal_ids[4],0))
            count_move_bin <- list(c(animal_ids[1],0), c(animal_ids[2],0), c(animal_ids[3],0), c(animal_ids[4],0), c(system_id,0))

            animal_list_bin <- initialize_animal_positions(animal_ids, bin_data, animal_list_bin)

            initial_time_bin <- animal_list_bin[[1]][[2]]
            current_row_bin <- 5
            total_rows_bin <- nrow(bin_data) + 1

            while(current_row_bin < total_rows_bin){
              prev_animal_list_bin <- animal_list_bin
              animal_list_bin <- update_animal_list(animal_ids, animal_list_bin, bin_data, initial_time_bin, current_row_bin)
              elapsed_seconds_bin <- animal_list_bin[["data_temp"]][["elapsed_seconds"]]
              if(system_complete){
                count_prox_bin <- update_proximity(prev_animal_list_bin, animal_list_bin, count_prox_bin, elapsed_seconds_bin)
                total_prox_bin <- update_total_proximity(prev_animal_list_bin, animal_list_bin, total_prox_bin, elapsed_seconds_bin)
              }
              count_move_bin <- update_movement(prev_animal_list_bin, animal_list_bin, count_move_bin, elapsed_seconds_bin)
              current_row_bin <- animal_list_bin[["data_temp"]][["current_row"]]
              initial_time_bin <- animal_list_bin[[1]][[2]]
            }

            for(i in 1:4){
              animal <- animal_ids[i]
              if(!is.na(animal)){
                bin_results <- add_row(bin_results,
                                       Batch=batch,
                                       CageChange=cageChange,
                                       System=system_id,
                                       Phase=phase,
                                       PhaseNumber=phase_number,
                                       Bin30min=as.POSIXct(bin_time),
                                       AnimalID=animal,
                                       TotalProximity=as.numeric(total_prox_bin[[i]][[2]]),
                                       TotalMovement=as.numeric(count_move_bin[[i]][[2]])
                )
              }
            }
          }
        }
      }
      if(system_complete){
        plot_total_proximity <- generateGraph(result_total_proximity, batch, cageChange, animal_ids, system_id)
        all_plots_total_proximity <- c(all_plots_total_proximity, list(plot_total_proximity))
      }
      allHeatmaps_positions[[batch]][[cageChange]][[system_id]] <- systemHeatmaps_positions
      systemHeatmaps_proximity <- list()
      systemHeatmaps_positions <- list()
    }

    if (save_plots == TRUE) {
      message("Saving plots")
      for (i in seq_along(all_plots_total_proximity)) {
        print(i)
        ggsave(
          filename = paste0(saving_directory, plots_directory, "/total_proximity_", batch, "_", cageChange, "_sys.", i, ".svg"),
          plot = all_plots_total_proximity[[i]],
          width = 5,
          height = 2
        )
      }
      for (i in seq_along(allHeatmaps_proximity)) {
        print(i)
        title <- allHeatmaps_proximity[[i]][[1]]$labels$title
        pattern <- "sys..."
        system_substring <- ifelse((match <- regexec(pattern, title))[[1]][1] > 0,
          regmatches(title, match)[[1]], "Pattern not found.")
        ggsave(
          filename = paste0(saving_directory, plots_directory, "/allHeatmaps_proximity_", batch, "_", cageChange, "_", system_substring, ".png"),
          plot = gridExtra::arrangeGrob(
            grobs = allHeatmaps_proximity[[i]], ncol = 2,
            layout_matrix = rbind(c(1, 5), c(2, 6), c(3, 7), c(4, 8), c(NA, 9))
          ),
          width = 12, height = 8
        )
        title <- allHeatmaps_positions[[i]][[1]]$labels$title
        system_substring <- ifelse((match <- regexec(pattern, title))[[1]][1] > 0,
          regmatches(title, match)[[1]], "Pattern not found.")
        ggsave(
          filename = paste0(saving_directory, plots_directory, "/allHeatmaps_positions_", batch, "_", cageChange, "_", system_substring, ".png"),
          plot = gridExtra::arrangeGrob(
            grobs = allHeatmaps_positions[[i]], ncol = 2,
            layout_matrix = rbind(c(1, 5), c(2, 6), c(3, 7), c(4, 8), c(NA, 9))
          ),
          width = 12, height = 8
        )
      }
    }

    if (save_tables == TRUE) {
      message("Saving tables")
      write.csv(result_total_proximity, file = paste0(saving_directory, tables_directory, "/", batch, "_", cageChange, "_total_proximity.csv"), row.names = FALSE)
      write.csv(result_total_movement, file = paste0(saving_directory, tables_directory, "/", batch, "_", cageChange, "_total_movement.csv"), row.names = FALSE)
      write.csv(system_animal_ids, file = paste0(saving_directory, tables_directory, "/", batch, "_", cageChange, "_animal_ids.csv"), row.names = FALSE)
      write.csv(bin_results, file = paste0(saving_directory, tables_directory, "/", batch, "_", cageChange, "_binned_proximity_movement.csv"), row.names = FALSE)
    }
  }
}

if (show_plots == TRUE) {
  message("show Plots")
  gridExtra::grid.arrange(grobs = all_plots_total_proximity, ncol = 2)
  for (batch in seq_along(allHeatmaps_positions)) {
    for (cc in seq_along(allHeatmaps_positions[[batch]])) {
      for (sys in seq_along(allHeatmaps_positions[[batch]][[cc]])) {
        print(paste(sys))
        if (cc != "CC4") {
          gridExtra::grid.arrange(
            grobs = allHeatmaps_positions[[batch]][[cc]][[sys]], ncol = 2,
            layout_matrix = rbind(c(1, 8), c(2, 5), c(3, 6), c(4, 7))
          )
        } else {
          gridExtra::grid.arrange(
            grobs = allHeatmaps_positions[[batch]][[cc]][[sys]], ncol = 2,
            layout_matrix = rbind(c(1, 10), c(2, 6), c(3, 7), c(4, 8), c(5, 9))
          )
        }
      }
    }
  }
}
