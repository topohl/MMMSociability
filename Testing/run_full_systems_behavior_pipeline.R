scripts_to_run <- c(
  "Testing/check_behavioral_dynamics_structure.R",
  "Analysis/06_burstiness_temporal_instability.R",
  "Analysis/07_behavioral_state_space.R",
  "Analysis/08_early_prediction_models.R",
  "Analysis/09_dynamic_social_networks.R",
  "Analysis/10_hmm_behavioral_states.R",
  "Analysis/11_gamm_trajectory_features.R",
  "Analysis/12_behavior_proteomics_integration.R"
)

run_status <- list()

for (scr in scripts_to_run) {

  cat("\n================================================\n")
  cat("Running:", scr, "\n")
  cat("================================================\n")

  start_time <- Sys.time()

  result <- tryCatch(
    {
      source(scr, local = new.env())
      list(success = TRUE, error = NA_character_)
    },
    error = function(e) {
      list(success = FALSE, error = conditionMessage(e))
    }
  )

  elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)

  run_status[[scr]] <- data.frame(
    script = scr,
    success = result$success,
    runtime_seconds = elapsed,
    error = result$error,
    stringsAsFactors = FALSE
  )

  if (isTRUE(result$success)) {
    cat("SUCCESS (", elapsed, " sec)\n", sep = "")
  } else {
    cat("FAILED (", elapsed, " sec)\n", sep = "")
    cat(result$error, "\n")
  }
}

status_tbl <- do.call(rbind, run_status)

out_dir <- "analysis_ready/06_behavioral_dynamics"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(
  status_tbl,
  file.path(out_dir, "full_systems_behavior_pipeline_status.csv")
)

print(status_tbl)
