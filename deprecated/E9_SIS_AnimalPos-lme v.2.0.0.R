# ---- 1. Package management ----
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  readr, tibble, dplyr, tidyr, purrr, lme4, lmerTest, emmeans, tcltk
)

# ---- 2. Paths and directory setup ----
working_directory <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability"
results_dir <- file.path(working_directory, "lme")
dirs <- list(
  logs        = file.path(results_dir, "logs"),
  models      = file.path(results_dir, "models"),
  tables_sum  = file.path(results_dir, "tables", "summary"),
  tables_emm  = file.path(results_dir, "tables", "emmeans"),
  tables_contr= file.path(results_dir, "tables", "contrasts"),
  plots_line  = file.path(results_dir, "plots", "line"),
  plots_pub   = file.path(results_dir, "plots", "publication")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

batches <- c("B1", "B2", "B3", "B4", "B5")
cageChanges <- c("CC1", "CC2", "CC3", "CC4")
sus_animals <- readLines(file.path(working_directory, "raw_data", "sus_animals.csv"))
con_animals <- readLines(file.path(working_directory, "raw_data", "con_animals.csv"))

get_sex <- function(batch) {
  ifelse(batch %in% c("B3", "B4", "B6"), "female", "male")
}

animal_to_system <- function(animal_ids, system_tbl) {
  map_chr(animal_ids, ~ {
    idx <- which(system_tbl == .x, arr.ind = TRUE)[2]
    if (length(idx) > 0) names(system_tbl)[idx] else NA
  })
}

# Helper: melt wide proximity data to long
melt_matrix_csv <- function(tbl, value_col) {
  tbl %>%
    tidyr::pivot_longer(-Phase, names_to = "AnimalID", values_to = value_col)
}

# Helper: Add metadata columns to long tibble
add_group <- function(df, sex, batch, change, sus, con, sys_tbl) {
  df %>%
    dplyr::mutate(
      Group = dplyr::case_when(
        AnimalID %in% sus ~ "sus",
        AnimalID %in% con ~ "con",
        TRUE              ~ "res"
      ),
      Sex = sex,
      Batch = batch,
      CageChange = change,
      System = animal_to_system(AnimalID, sys_tbl)
    )
}

# Read and merge all proximity data files into one tibble
merge_all_data_long <- function() {
  message("Start merging all proximity data for all batches and cage changes")
  prox_list <- list()
  for (batch in batches) {
    sex <- get_sex(batch)
    for (change in cageChanges) {
      animal_ids_path <- file.path(working_directory, "tables", paste0(batch, "_", change, "_animal_ids.csv"))
      if(!file.exists(animal_ids_path)) next
      system_tbl <- readr::read_csv(animal_ids_path, show_col_types = FALSE)
      prox_path <- file.path(working_directory, "tables", paste0(batch, "_", change, "_total_proximity.csv"))
      if (!file.exists(prox_path)) { next }

      prox_tbl <- readr::read_csv(prox_path, show_col_types = FALSE)
      if(change=="CC4") prox_tbl <- dplyr::filter(prox_tbl, Phase %in% c("A1", "I2", "A2"))
      prox_long <- melt_matrix_csv(prox_tbl, "SocialProx")
      prox_long <- add_group(prox_long, sex, batch, change, sus_animals, con_animals, system_tbl)

      prox_list[[paste(batch, change)]] <- prox_long
      message(sprintf("Merged data for batch %s, change %s: %d rows", batch, change, nrow(prox_long)))
    }
  }
  full_prox <- dplyr::bind_rows(prox_list)
  message(sprintf("Full merged proximity data rows: %d", nrow(full_prox)))
  return(full_prox)
}

proximity_data <- merge_all_data_long()

# GUI for user selections
get_user_input <- function() {
  tt <- tktoplevel()
  tkwm.title(tt, "Analysis Settings")

  includeChange_var <- tclVar("0")
  includeSex_var    <- tclVar("0")
  includePhase_var  <- tclVar("0")

  frame <- tkframe(tt, relief = "groove", borderwidth = 2, padx = 20, pady = 15)
  tkpack(frame, fill = "both", expand = TRUE)

  tklabel(frame, text = "Select factors to include in analysis", font = c("Helvetica", 14, "bold")) %>% tkpack(pady = c(0, 15))

  add_checkbox <- function(label_text, var) {
    frame_inner <- tkframe(frame)
    tkpack(frame_inner, fill = "x", pady = 5)
    tklabel(frame_inner, text = label_text, font = c("Helvetica", 12)) %>% tkpack(side = "left", anchor = "w")
    tkcheckbutton(frame_inner, variable = var) %>% tkpack(side = "right")
  }

  add_checkbox("Include Cage Change (Batch)", includeChange_var)
  add_checkbox("Include Sex", includeSex_var)
  add_checkbox("Include Phase", includePhase_var)

  # OK button
  ok_button <- tkbutton(frame, text = "OK", font = c("Helvetica", 12, "bold"),
                       bg = "#007AFF", fg = "white", relief = "flat",
                       command = function() tkdestroy(tt))
  tkpack(ok_button, pady = 20, fill = "x")

  tkwait.window(tt)

  list(
    includeChange = as.logical(as.integer(tclvalue(includeChange_var))),
    includeSex = as.logical(as.integer(tclvalue(includeSex_var))),
    includePhase = as.logical(as.integer(tclvalue(includePhase_var)))
  )
}
ui <- get_user_input()
includeChange <- ui$includeChange
includeSex <- ui$includeSex
includePhase <- ui$includePhase
message(sprintf("User selected: includeChange=%s, includeSex=%s, includePhase=%s", includeChange, includeSex, includePhase))

# Helper to build model formula rhs dynamically
build_fixed_effects <- function(includeSex, includePhase) {
  rhs <- "Group"
  if (includePhase) rhs <- paste0(rhs, " * Phase")
  if (includeSex) rhs <- paste0(rhs, " * Sex")
  return(rhs)
}

# Adaptive model fitting with fallback on simpler random effects if needed
fit_lmer_adaptive <- function(response_var, fixed_rhs, df, time_var = "Phase", animal = "AnimalID", allow_slope = TRUE) {
  ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  formula_str <- function(random_part) paste0(response_var, " ~ ", fixed_rhs, " + ", random_part)
  if (!allow_slope) {
    f <- as.formula(formula_str(paste0("(1 | ", animal, ")")))
    return(lmerTest::lmer(f, data = df, control = ctrl))
  }
  f1 <- as.formula(formula_str(paste0("(1 + ", time_var, " | ", animal, ")")))
  fit1 <- try(lmerTest::lmer(f1, data = df, control = ctrl), silent = TRUE)
  if (!inherits(fit1, "try-error") && !lme4::isSingular(fit1, tol = 1e-5)) return(fit1)
  f2 <- as.formula(formula_str(paste0("(1 | ", animal, ") + (0 + ", time_var, " | ", animal, ")")))
  fit2 <- try(lmerTest::lmer(f2, data = df, control = ctrl), silent = TRUE)
  if (!inherits(fit2, "try-error") && !lme4::isSingular(fit2, tol = 1e-5)) return(fit2)
  f3 <- as.formula(formula_str(paste0("(1 | ", animal, ")")))
  return(lmerTest::lmer(f3, data = df, control = ctrl))
}

# Compute phase average means and contrasts via emmeans
compute_phase_avg_for_model <- function(model, Change, Sex, Phase, adjust_method = "holm") {
  emm_grid <- emmeans::emmeans(model, ~ Group)
  emm_df <- as.data.frame(emm_grid)
  avg_df <- emm_df %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise( emmean_avg = mean(emmean, na.rm=TRUE), SE_avg = sqrt(mean(SE^2, na.rm=TRUE)) ) %>%
    dplyr::mutate(Change = Change, Sex = Sex, Phase = Phase)
  # Contrast vectors
  con_mat <- list(
    "RES-CON" = c(-1,1,0),
    "SUS-CON" = c(-1,0,1),
    "RES-SUS" = c(0,1,-1)
  )
  contr_list <- lapply(names(con_mat), function(nm) {
    L <- matrix(con_mat[[nm]], nrow=1)
    est <- as.numeric(L %*% avg_df$emmean_avg)
    se  <- sqrt(as.numeric(L %*% diag(avg_df$SE_avg^2) %*% t(L)))
    data.frame(contrast=nm, estimate=est, SE=se, Change=Change, Sex=Sex, Phase=Phase)
  })
  contr_df <- dplyr::bind_rows(contr_list)
  avg_df <- avg_df %>%
    dplyr::mutate(lwr = emmean_avg - 1.96 * SE_avg, upr = emmean_avg + 1.96 * SE_avg)
  contr_df <- contr_df %>%
    dplyr::mutate(lwr = estimate - 1.96 * SE, upr = estimate + 1.96 * SE)
  return(list(means = avg_df, contrasts = contr_df))
}

# ----- Main analysis loop -----
change_vals <- if(includeChange) unique(proximity_data$Batch) else "all"
sex_vals <- if(includeSex) unique(proximity_data$Sex) else "all"
phase_vals <- if(includePhase) unique(proximity_data$Phase) else "all"

phase_avg_means_all <- list()
phase_avg_contrasts_all <- list()

for(chg in change_vals) {
  for(sx in sex_vals) {
    for(ph in phase_vals) {
      subset_df <- proximity_data %>%
        dplyr::filter(
          (chg == "all" | Batch == chg),
          (sx == "all"  | Sex == sx),
          (ph == "all"  | Phase == ph)
        )
      message(sprintf("Processing subset: Change=%s, Sex=%s, Phase=%s, Rows=%d", chg, sx, ph, nrow(subset_df)))
      if(nrow(subset_df) < 10) { message("Skipping: insufficient rows"); next }
      if(length(unique(subset_df$Group)) < 2) { message("Skipping: insufficient groups"); next }
      if(includeSex && length(unique(subset_df$Sex)) < 2) { message("Skipping: insufficient sex levels"); next }
      if(includePhase && length(unique(subset_df$Phase)) < 2) { message("Skipping: insufficient phases"); next }
      fe_rhs <- build_fixed_effects(includeSex, includePhase)
      mod <- fit_lmer_adaptive("SocialProx", fe_rhs, subset_df, time_var = ifelse(includePhase,"Phase","1"), animal = "AnimalID")
      stats <- compute_phase_avg_for_model(mod, chg, sx, ph)
      phase_avg_means_all[[length(phase_avg_means_all)+1]] <- stats$means
      phase_avg_contrasts_all[[length(phase_avg_contrasts_all)+1]] <- stats$contrasts
    }
  }
}

# Bind all results together
phase_avg_means_tbl <- dplyr::bind_rows(phase_avg_means_all)
phase_avg_contrasts_tbl <- dplyr::bind_rows(phase_avg_contrasts_all)

# Save outputs
readr::write_csv(phase_avg_means_tbl, file.path(dirs$tables_emm, "socialprox_phaseavg_means.csv"))
readr::write_csv(phase_avg_contrasts_tbl, file.path(dirs$tables_contr, "socialprox_phaseavg_contrasts.csv"))

message("Analysis complete! Output saved in:")
message(dirs$tables_emm)
message(dirs$tables_contr)
