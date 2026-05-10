# ================================================================
# Behavioral Dynamics Pipeline Runner
# MMMSociability
# ================================================================
# Runs:
#   1) structure checks
#   2) burstiness / instability analysis
#   3) behavioral state-space analysis
#   4) early prediction modeling
#
# Recommended workflow:
#   source("Testing/run_behavioral_dynamics_pipeline.R")
# ================================================================

cat("\n================================================\n")
cat("MMMSociability behavioral dynamics pipeline\n")
cat("================================================\n\n")

scripts_to_run <- c(
  "Testing/check_behavioral_dynamics_structure.R",
  "Analysis/06_burstiness_temporal_instability.R",
  "Analysis/07_behavioral_state_space.R",
  "Analysis/08_early_prediction_models.R"
)

run_status <- list()

for (scr in scripts_to_run) {

  cat("\n------------------------------------------------\n")
  cat("Running:", scr, "\n")
  cat("------------------------------------------------\n")

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
    cat("Error:\n")
    cat(result$error, "\n")
  }
}

status_tbl <- do.call(rbind, run_status)

out_dir <- "analysis_ready/06_behavioral_dynamics"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

readr::write_csv(
  status_tbl,
  file.path(out_dir, "behavioral_dynamics_pipeline_status.csv")
)

cat("\n================================================\n")
cat("Pipeline finished\n")
cat("================================================\n\n")

print(status_tbl)
