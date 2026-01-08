###############################################
# Mixed models with random time slopes + fixed-grid phase-edge windows
# Modified to analyze both Movement and Proximity metrics
###############################################

suppressPackageStartupMessages({
  packages <- c(
    "ggplot2","dplyr","tidyr","gridExtra","lme4",
    "lmerTest","cowplot","lsmeans","emmeans","Matrix",
    "tcltk","openxlsx","readr","stringr","patchwork","scales",
    "broom.mixed","pbkrtest","rlang","tibble"
  )
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) install.packages(package, dependencies = TRUE)
    library(package, character.only = TRUE)
  }
})
options(contrasts = c("contr.sum","contr.poly"))

# -------------------------------------------------
# Options and thresholds
# -------------------------------------------------
use_group_time_interaction <- TRUE
min_animals <- 2
min_groups  <- 2
min_obs     <- 8
min_time_lv <- 3

# Suppress specific emmeans notes for large datasets
options(
  emmeans = list(
    pbkrtest.limit = 10000,  # Increase limit for Kenward-Roger df adjustments
    lmerTest.limit = 10000,  # Increase limit for Satterthwaite df adjustments
    msg.interaction = FALSE  # Suppress interaction involvement warnings
  )
)

# -------------------------------------------------
# Base Paths
# -------------------------------------------------
#base_results_dir <- "D:/MMMSociability/statistics/lme/"
base_results_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability/statistics/lme/"
if (!dir.exists(base_results_dir)) dir.create(base_results_dir, recursive = TRUE)
dir_create_safe <- function(path) { if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE) }

# -------------------------------------------------
# Helpers
# -------------------------------------------------
scope_all <- function(includeChange, includeSex, includePhase) paste0(
  if (includeChange) "subsetChanges" else "allChanges", "_",
  if (includeSex)    "subsetSexes"  else "allSexes",  "_",
  if (includePhase)  "subsetPhases" else "allPhases"
)

quiet_singular <- function(expr) {
  withCallingHandlers(expr,
    warning = function(w) {
      if (grepl("boundary \\(singular\\) fit", conditionMessage(w))) invokeRestart("muffleWarning")
    }
  )
}

has_cols <- function(x, cols) all(cols %in% names(x))

# -------------------------------------------------
# LOAD DATA
# -------------------------------------------------
cat("Loading data...\n")
#data_file <- "D:/MMMSociability/processed_data/data_filtered_agg.csv"
data_file <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability/processed_data/data_lme_format/data_filtered_agg.csv"

if (file.exists(data_file)) {
  data_filtered_agg <- read_csv(data_file, show_col_types = FALSE)
  cat(sprintf("Data loaded: %d rows, %d columns\n", nrow(data_filtered_agg), ncol(data_filtered_agg)))
} else {
  # Alternative: load from RDS (faster)
  data_file_rds <- "D:/MMMSociability/processed_data/data_filtered_agg.rds"
  if (file.exists(data_file_rds)) {
    data_filtered_agg <- readRDS(data_file_rds)
    cat(sprintf("Data loaded from RDS: %d rows, %d columns\n", nrow(data_filtered_agg), ncol(data_filtered_agg)))
  } else {
    stop("Data file not found! Please run the data formatting script first.")
  }
}

# -------------------------------------------------
# Modified fit_lmer_adaptive - response_var is now required
# -------------------------------------------------
fit_lmer_adaptive <- function(fixed_rhs, df, response_var, time_var = "TH", animal = "AnimalNum", allow_slope = TRUE) {
  ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  if (!allow_slope) {
    f <- as.formula(paste0(response_var, " ~ ", fixed_rhs, " + (1 | ", animal, ")"))
    return(lmerTest::lmer(f, data = df, control = ctrl))
  }
  f1 <- as.formula(paste0(response_var, " ~ ", fixed_rhs, " + (1 + ", time_var, " | ", animal, ")"))
  fit1 <- try(lmerTest::lmer(f1, data = df, control = ctrl), silent = TRUE)
  if (!inherits(fit1, "try-error") && !lme4::isSingular(fit1, tol = 1e-5)) return(fit1)
  f2 <- as.formula(paste0(response_var, " ~ ", fixed_rhs, " + (1 | ", animal, ") + (0 + ", time_var, " | ", animal, ")"))
  fit2 <- try(lmerTest::lmer(f2, data = df, control = ctrl), silent = TRUE)
  if (!inherits(fit2, "try-error") && !lme4::isSingular(fit2, tol = 1e-5)) return(fit2)
  f3 <- as.formula(paste0(response_var, " ~ ", fixed_rhs, " + (1 | ", animal, ")"))
  lmerTest::lmer(f3, data = df, control = ctrl)
}

compute_phase_avg_for_model <- function(model, Change, Sex, Phase, adjust_method = "holm") {
  th_step <- 0.5
  th_grid <- seq(0, 12, by = th_step)
  
  emm_grid <- emmeans::emmeans(model, ~ Group, at = list(TH = th_grid))
  emm_grid_df <- as.data.frame(emm_grid)
  
  avg_df <- emm_grid_df %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(
      emmean_avg = mean(emmean, na.rm = TRUE),
      SE_avg     = sqrt(mean(SE^2, na.rm = TRUE) / length(th_grid)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(Change = as.character(Change),
                  Sex    = as.character(Sex),
                  Phase  = as.character(Phase))
  
  avg_df$Group <- factor(avg_df$Group, levels = c("CON","RES","SUS"))
  get_contrast <- function(L) {
    est <- as.numeric(L %*% avg_df$emmean_avg)
    V   <- diag(avg_df$SE_avg^2, nrow = length(avg_df$SE_avg))
    se  <- sqrt(as.numeric(L %*% V %*% t(L)))
    c(estimate = est, SE = se)
  }
  L_list <- list(
    "RES - CON" = c(-1, 1, 0),
    "SUS - CON" = c(-1, 0, 1),
    "RES - SUS" = c( 0, 1, -1)
  )
  contr_rows <- lapply(names(L_list), function(nm) {
    L <- matrix(L_list[[nm]], nrow = 1)
    est <- get_contrast(L)
    data.frame(
      contrast = nm, estimate = est["estimate"], SE = est["SE"],
      Change = as.character(Change), Sex = as.character(Sex), Phase = as.character(Phase),
      stringsAsFactors = FALSE
    )
  })
  contr_df <- dplyr::bind_rows(contr_rows)
  
  df_rep <- tryCatch({
    sm <- summary(model)
    as.numeric(sm$coefficients[,"df"][1])
  }, error = function(e) NA_real_)
  if (!is.finite(df_rep)) df_rep <- 200
  
  tcrit <- qt(0.975, df = df_rep)
  avg_df <- avg_df %>%
    dplyr::mutate(
      lwr = emmean_avg - tcrit * SE_avg,
      upr = emmean_avg + tcrit * SE_avg
    )
  contr_df <- contr_df %>%
    dplyr::mutate(
      t.ratio = estimate / SE,
      df = df_rep,
      p.value = 2 * pt(-abs(t.ratio), df = df_rep),
      p.adjust = p.adjust(p.value, method = adjust_method),
      lwr = estimate - tcrit * SE,
      upr = estimate + tcrit * SE
    )
  
  list(means = avg_df, contrasts = contr_df)
}

# -------------------------------------------------
# GUI
# -------------------------------------------------
get_user_input <- function() {
  tt <- tktoplevel(); tkwm.title(tt, "Analysis Settings"); tkconfigure(tt, bg = "#F7F9FC", padx = 10, pady = 20)
  includeChange_var <- tclVar(0); includeSex_var <- tclVar(0); includePhase_var <- tclVar(0)
  frame <- tkframe(tt, bg = "#F7F9FC"); tkgrid(frame, padx = 10, pady = 10)
  tkgrid(tklabel(frame, text = "Set Inclusion", bg = "#F7F9FC", font = c("Helvetica", 14)))
  tkgrid(tklabel(frame, text = "Parameters",   bg = "#F7F9FC", font = c("Helvetica", 14)))
  tkgrid(tklabel(frame, text="Change", bg="#F7F9FC", font=c("Helvetica",12)),
         tkcheckbutton(frame, variable = includeChange_var, bg="#F7F9FC"), padx=5, pady=5)
  tkgrid(tklabel(frame, text="Sex", bg="#F7F9FC", font=c("Helvetica",12)),
         tkcheckbutton(frame, variable = includeSex_var, bg="#F7F9FC"), padx=5, pady=5)
  tkgrid(tklabel(frame, text="Phase", bg="#F7F9FC", font=c("Helvetica",12)),
         tkcheckbutton(frame, variable = includePhase_var, bg="#F7F9FC"), padx=5, pady=5)
  onOK <- function() tkdestroy(tt)
  ok_button <- tkbutton(tt, text="OK", command=onOK, bg="#007AFF", fg="white", relief="flat", font=c("Helvetica",12))
  tkgrid(ok_button, padx=10, pady=10); tkwait.window(tt)
  list(
    includeChange = as.logical(as.numeric(tclvalue(includeChange_var))),
    includeSex    = as.logical(as.numeric(tclvalue(includeSex_var))),
    includePhase  = as.logical(as.numeric(tclvalue(includePhase_var)))
  )
}
ui <- get_user_input()
includeChange <- ui$includeChange
includeSex    <- ui$includeSex
includePhase  <- ui$includePhase

# -------------------------------------------------
# Data prep - check for Movement and Proximity
# -------------------------------------------------
stopifnot(all(c("Change","Sex","Phase","Group","HalfHourElapsed","Movement","Proximity","AnimalNum") %in% names(data_filtered_agg)))
data_filtered_agg <- data_filtered_agg %>%
  mutate(
    Change = factor(Change),
    Sex    = factor(Sex),
    Phase  = factor(Phase, levels = c("Active","Inactive")),
    Group  = factor(Group, levels = c("CON","RES","SUS")),
    HalfHourElapsed = as.integer(HalfHourElapsed)
  )

# -------------------------------------------------
# Main analysis wrapper function
# -------------------------------------------------
run_analysis_for_metric <- function(metric_name, data, includeChange, includeSex, includePhase) {
  
  cat(sprintf("\n\n========================================\n"))
  cat(sprintf("ANALYZING METRIC: %s\n", metric_name))
  cat(sprintf("========================================\n\n"))
  
  # Create metric-specific directories
  results_dir <- file.path(base_results_dir, metric_name)
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
  
  dirs <- list(
    logs        = file.path(results_dir, "logs"),
    models      = file.path(results_dir, "models"),
    tables_sum  = file.path(results_dir, "tables", "summary"),
    tables_emm  = file.path(results_dir, "tables", "emmeans"),
    plots_line  = file.path(results_dir, "plots", "line"),
    plots_heat  = file.path(results_dir, "plots", "heatmap"),
    plots_fixed = file.path(results_dir, "plots", "fixed_effects"),
    plots_pub   = file.path(results_dir, "plots", "publication")
  )
  invisible(lapply(dirs, dir_create_safe))
  dirs$tables_contr <- file.path(results_dir, "tables", "contrasts")
  dirs$tables_inv   <- file.path(results_dir, "tables", "inventory")
  dirs$tables_man   <- file.path(results_dir, "tables", "manifest")
  dirs$tables_phase <- file.path(results_dir, "tables", "phase_average")
  invisible(lapply(dirs[c("tables_contr","tables_inv","tables_man","tables_phase")], dir_create_safe))
  
  dirs$edgeAggAll <- file.path(results_dir, "plots", "edge_emm")
  dir_create_safe(dirs$edgeAggAll)
  
  # Directory helper functions
  plot_dir_line  <- function(Change, Sex, includeSex) { 
    sub <- file.path(dirs$plots_line, paste0("CC-", Change)); 
    dir_create_safe(sub); 
    if (isTRUE(includeSex)) { 
      sub <- file.path(sub, paste0("Sex-", Sex)); 
      dir_create_safe(sub) 
    }
    sub 
  }
  plot_dir_heat  <- function(Change, Sex, includeSex) { 
    sub <- file.path(dirs$plots_heat, paste0("CC-", Change)); 
    dir_create_safe(sub); 
    if (isTRUE(includeSex)) { 
      sub <- file.path(sub, paste0("Sex-", Sex)); 
      dir_create_safe(sub) 
    }
    sub 
  }
  plot_dir_fixed <- function(Change, Sex) { 
    sub <- file.path(dirs$plots_fixed, paste0("CC-", Change)); 
    dir_create_safe(sub); 
    sub2 <- file.path(sub, paste0("Sex-", as.character(Sex))); 
    dir_create_safe(sub2); 
    sub2 
  }
  
  cat(sprintf("Results will be saved to: %s\n", results_dir))
  cat(sprintf("Data summary for %s:\n", metric_name))
  cat(sprintf("  - Mean: %.2f\n", mean(data[[metric_name]], na.rm = TRUE)))
  cat(sprintf("  - SD: %.2f\n", sd(data[[metric_name]], na.rm = TRUE)))
  cat(sprintf("  - Range: [%.2f, %.2f]\n", min(data[[metric_name]], na.rm = TRUE), max(data[[metric_name]], na.rm = TRUE)))
  cat(sprintf("  - Missing values: %d\n\n", sum(is.na(data[[metric_name]]))))
  
  # -------------------------------------------------
  # Line and heatmap plots
  # -------------------------------------------------
  data_by_change <- split(data, data$Change)
  
  for (Change in sort(unique(data$Change))) {
    df_change <- data_by_change[[as.character(Change)]]
    for (SexValue in c("m","f")) {
      df_sex <- if (includeSex) dplyr::filter(df_change, Sex == SexValue) else df_change
      pal <- c("#457B9D", "#C6C3BB", "#E63946")
      p <- ggplot(df_sex, aes(HalfHourElapsed, .data[[metric_name]], color = Group, group = Group)) +
        geom_path(stat = "summary", fun = mean, linewidth = 0.8) +
        stat_summary(aes(fill = Group), fun = mean,
                     fun.min = function(x) mean(x) - sd(x),
                     fun.max = function(x) mean(x) + sd(x),
                     geom = "ribbon", alpha = 0.3, colour = NA) +
        scale_color_manual(values = pal) + scale_fill_manual(values = pal) +
        scale_x_continuous(labels = function(x) x/2, breaks = scales::pretty_breaks(6)) +
        labs(title = paste("Circadian", metric_name), x = "Time Elapsed [h]", y = paste(metric_name, "[a.u.]")) +
        theme_minimal(base_size = 14) +
        theme(panel.grid.minor = element_blank(), legend.position = "top", plot.title = element_text(face="bold", hjust=0.5))
      out_dir <- plot_dir_line(as.character(Change), SexValue, includeSex)
      ggsave(file.path(out_dir, paste0("line_", Change, "_", if (includeSex) SexValue else "allSexes", ".svg")), p, width = 7, height = 4)
    }
  }
  
  for (Change in sort(unique(data$Change))) {
    df_change <- data_by_change[[as.character(Change)]]
    for (SexValue in c("m","f")) {
      df_sex <- if (includeSex) dplyr::filter(df_change, Sex == SexValue) else df_change
      p <- ggplot(df_sex, aes(HalfHourElapsed, reorder(Group, -as.numeric(Group)), fill = .data[[metric_name]])) +
        geom_tile() +
        scale_fill_gradient(low = "white", high = "darkblue", limits = c(0, max(data[[metric_name]]))) +
        scale_x_continuous(labels = function(x) x/2, breaks = scales::pretty_breaks(6)) +
        facet_grid(Sex ~ ., scales = "free") +
        cowplot::theme_minimal_hgrid(12, rel_small = 1) +
        labs(title = bquote(~bold(.(paste(metric_name)))),
             subtitle = paste("Cage Change:", Change), x = "Time elapsed [h]", y = "Group") +
        theme(legend.position = "top", plot.title = element_text(hjust = 0.5))
      out_dir <- plot_dir_heat(as.character(Change), SexValue, includeSex)
      ggsave(file.path(out_dir, paste0("heatmap_", Change, "_", if (includeSex) SexValue else "allSexes", ".svg")), p, width = 5, height = 5)
    }
  }
  
  # -------------------------------------------------
  # Modeling (main slices)
  # -------------------------------------------------
  changes <- if (includeChange) unique(data$Change) else "allChanges"
  sexes   <- if (includeSex)    unique(data$Sex)    else "allSexes"
  phases  <- if (includePhase)  unique(data$Phase)  else "allPhases"
  
  summary_results <- tibble::tibble()
  emmeans_results <- tibble::tibble()
  log_file <- file.path(dirs$logs, paste0("fit_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
  
  total_iterations <- length(changes) * length(sexes) * length(phases)
  pb <- tkProgressBar(title = paste("Model Fitting Progress:", metric_name), min = 0, max = total_iterations, width = 300)
  i <- 0
  
  # Collectors
  phase_avg_means_all <- list()
  phase_avg_contr_all <- list()
  
  for (Change in changes) {
    for (Sex in sexes) {
      for (Phase in phases) {
        i <- i + 1
        setTkProgressBar(pb, i, label = paste("Progress:", round(i/total_iterations*100,2), "%"))
        dsub <- data %>%
          dplyr::filter(
            (!includeChange | Change == !!Change),
            (!includeSex    | Sex    == !!Sex),
            (!includePhase  | Phase  == !!Phase)
          )
        if (nrow(dsub) == 0) next
        
        n_anim  <- dplyr::n_distinct(dsub$AnimalNum)
        n_group <- dplyr::n_distinct(dsub$Group)
        if (n_anim < min_animals || n_group < min_groups || nrow(dsub) < min_obs) {
          cat(paste(Sys.time(), Change, Sex, Phase, "SKIP main: few_animals/groups/obs\n"), file = log_file, append = TRUE)
          next
        }
        
        t0 <- min(dsub$HalfHourElapsed, na.rm = TRUE)
        dsub <- dsub %>% dplyr::mutate(TH = (HalfHourElapsed - t0)/2)
        use_int <- use_group_time_interaction && (dplyr::n_distinct(dsub$TH) >= min_time_lv)
        fixed_rhs <- if (use_int) "Group * TH" else "Group + TH"
        per_animal <- dsub %>% dplyr::count(AnimalNum, name = "n_i")
        allow_slope <- all(per_animal$n_i >= 3) && (dplyr::n_distinct(dsub$TH) >= min_time_lv)
        
        model <- quiet_singular(fit_lmer_adaptive(fixed_rhs, dsub, response_var = metric_name, time_var = "TH", allow_slope = allow_slope))
        ms <- summary(model)
        cat(paste(Sys.time(), Change, Sex, Phase, "MAIN form:", deparse(formula(model)),
                  "singular:", lme4::isSingular(model), "\n"), file = log_file, append = TRUE)

        # Save model object
        model_filename <- paste0(
        "model_", 
        metric_name, "_",
        "Change-", Change, "_",
        "Sex-", Sex, "_",
        "Phase-", Phase, 
        ".rds"
        )
        saveRDS(model, file.path(dirs$models, model_filename))

        # Save model summary as text
        summary_filename <- paste0(
        "summary_", 
          metric_name, "_",
          "Change-", Change, "_",
          "Sex-", Sex, "_",
          "Phase-", Phase, 
          ".txt"
        )
        sink(file.path(dirs$models, summary_filename))
        print(summary(model))
        cat("\n\nFormula:\n")
        print(formula(model))
        cat("\n\nRandom effects:\n")
        print(VarCorr(model))
        sink()
        
        fx <- tibble::tibble(
          Fixed_effect = rownames(ms$coefficients),
          Estimate     = unname(ms$coefficients[, "Estimate"]),
          Std.Error    = unname(ms$coefficients[, "Std. Error"]),
          df           = unname(ms$coefficients[, "df"]),
          t.value      = unname(ms$coefficients[, "t value"]),
          p.value      = unname(ms$coefficients[, grep("Pr\\(>|t\\)", colnames(ms$coefficients))]),
          p.round      = round(p.value, 3),
          p.sign       = dplyr::case_when(
            p.value < 0.001 ~ "***",
            p.value < 0.01  ~ "**",
            p.value < 0.05  ~ "*",
            p.value < 0.1   ~ "T",
            TRUE ~ "ns"
          ),
          Change = as.character(Change), Phase = as.character(Phase), Sex = as.character(Sex)
        )
        summary_results <- dplyr::bind_rows(summary_results, fx)
        
        if (!any(grepl("^Group", names(fixef(model))))) next
        
        emm_lvl <- emmeans::emmeans(model, pairwise ~ Group, at = list(TH = 0))
        emm_lvl_df <- as.data.frame(summary(emm_lvl$contrasts, infer = TRUE, adjust = "holm")) %>%
          dplyr::mutate(contrast_type = "level@start",
                        Change = as.character(Change), Phase = as.character(Phase), Sex = as.character(Sex))
        emm_lvl_df <- emm_lvl_df[!is.na(emm_lvl_df$p.value), ]
        emm_lvl_df$p.adjust <- p.adjust(emm_lvl_df$p.value, "holm")
        if (nrow(emm_lvl_df) > 0) emmeans_results <- dplyr::bind_rows(emmeans_results, emm_lvl_df)
        
        if (use_int && any(grepl("Group:TH", names(fixef(model))))) {
          emm_tr <- emmeans::emtrends(model, pairwise ~ Group, var = "TH")
          emm_tr_df <- as.data.frame(summary(emm_tr$contrasts, infer = TRUE, adjust = "holm")) %>%
            dplyr::mutate(contrast_type = "slope",
                          Change = as.character(Change), Phase = as.character(Phase), Sex = as.character(Sex))
          emm_tr_df <- emm_tr_df[!is.na(emm_tr_df$p.value), ]
          emm_tr_df$p.adjust <- p.adjust(emm_tr_df$p.value, "holm")
          if (nrow(emm_tr_df) > 0) emmeans_results <- dplyr::bind_rows(emmeans_results, emm_tr_df)
        }
        if (exists("model") && inherits(model, "lmerMod")) {
          stat_k <- compute_phase_avg_for_model(model, Change, Sex, Phase, adjust_method = "holm")
          if (!is.null(stat_k$means)     && nrow(stat_k$means))     phase_avg_means_all[[length(phase_avg_means_all)+1]]   <- stat_k$means
          if (!is.null(stat_k$contrasts) && nrow(stat_k$contrasts)) phase_avg_contr_all[[length(phase_avg_contr_all)+1]] <- stat_k$contrasts
        }
      }
    }
  }
  close(pb)
  
  run_scope <- scope_all(includeChange, includeSex, includePhase)
  summary_xlsx <- file.path(dirs$tables_sum, paste0("lme_summary_results_", run_scope, ".xlsx"))
  emmeans_xlsx <- file.path(dirs$tables_emm, paste0("emmeans_results_", run_scope, ".xlsx"))
  openxlsx::write.xlsx(summary_results, summary_xlsx, rowNames = FALSE)
  openxlsx::write.xlsx(emmeans_results, emmeans_xlsx, rowNames = FALSE)
  
  phase_avg_means_tbl <- if (length(phase_avg_means_all)) dplyr::bind_rows(phase_avg_means_all) else tibble::tibble()
  phase_avg_contr_tbl <- if (length(phase_avg_contr_all)) dplyr::bind_rows(phase_avg_contr_all) else tibble::tibble()
  
  if (nrow(phase_avg_means_tbl) == 0 && nrow(phase_avg_contr_tbl) == 0) {
    message("No phase-average results were generated for the current selection.")
  } else {
    if (nrow(phase_avg_means_tbl)) {
      phase_avg_means_tbl <- phase_avg_means_tbl %>%
        dplyr::mutate(
          Change = factor(Change, levels = unique(Change)),
          Phase  = factor(Phase,  levels = c("Active","Inactive")),
          Sex    = factor(Sex,    levels = unique(Sex)),
          Group  = factor(Group,  levels = c("CON","RES","SUS"))
        ) %>%
        dplyr::distinct(Change, Sex, Phase, Group, .keep_all = TRUE)
    }
    if (nrow(phase_avg_contr_tbl)) {
      phase_avg_contr_tbl <- phase_avg_contr_tbl %>%
        dplyr::mutate(
          Change = factor(Change, levels = unique(Change)),
          Phase  = factor(Phase,  levels = c("Active","Inactive")),
          Sex    = factor(Sex,    levels = unique(Sex))
        )
    }
    
    # Save tables
    openxlsx::write.xlsx(phase_avg_means_tbl, file.path(dirs$tables_emm, paste0("phaseAvg_emmeans_", run_scope, ".xlsx")), rowNames = FALSE)
    openxlsx::write.xlsx(phase_avg_contr_tbl, file.path(dirs$tables_emm, paste0("phaseAvg_contrasts_", run_scope, ".xlsx")), rowNames = FALSE)
    
    # Canonicalize contrast function
    canonicalize_contrast <- function(df, keep_pairs = c("RES-CON","SUS-CON","RES-SUS")) {
      if (!nrow(df)) return(df)
      parts <- strsplit(df$contrast, "\\s*-\\s*")
      lhs <- vapply(parts, `[`, character(1), 1)
      rhs <- vapply(parts, `[`, character(1), 2)
      lhs <- trimws(lhs); rhs <- trimws(rhs)
      
      flip_vs_con <- lhs == "CON" & rhs %in% c("RES","SUS")
      if (any(flip_vs_con)) {
        df$estimate[flip_vs_con] <- -df$estimate[flip_vs_con]
        if (all(c("lwr","upr") %in% names(df))) {
          lwr_new <- -df$upr[flip_vs_con]
          upr_new <- -df$lwr[flip_vs_con]
          df$lwr[flip_vs_con] <- pmin(lwr_new, upr_new, na.rm = TRUE)
          df$upr[flip_vs_con] <- pmax(lwr_new, upr_new, na.rm = TRUE)
        }
        df$contrast[flip_vs_con] <- paste(rhs[flip_vs_con], "- CON")
      }
      
      df$pair <- dplyr::case_when(
        grepl("^RES\\s*-\\s*CON$", df$contrast) ~ "RES-CON",
        grepl("^SUS\\s*-\\s*CON$", df$contrast) ~ "SUS-CON",
        grepl("^RES\\s*-\\s*SUS$", df$contrast) ~ "RES-SUS",
        TRUE ~ NA_character_
      )
      
      df <- df %>% dplyr::filter(pair %in% keep_pairs)
      
      if (!all(c("lwr","upr") %in% names(df)) && all(c("SE","df") %in% names(df))) {
        df$tcrit <- qt(0.975, df = df$df)
        df$lwr   <- df$estimate - df$tcrit * df$SE
        df$upr   <- df$estimate + df$tcrit * df$SE
      }
      
      df %>% dplyr::distinct(Change, Sex, Phase, pair, .keep_all = TRUE)
    }
    
    phase_avg_contr_tbl <- canonicalize_contrast(phase_avg_contr_tbl, keep_pairs = c("RES-CON","SUS-CON","RES-SUS"))
    
    # Forest plot
    if (nrow(phase_avg_contr_tbl) > 0) {
      forest_phaseavg <- ggplot(
        phase_avg_contr_tbl,
        aes(x = Change, y = estimate, ymin = lwr, ymax = upr, color = pair, shape = pair)
      ) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
        geom_pointrange(position = position_dodge(width = 0.55), linewidth = 0.6) +
        geom_text(
          aes(label = ifelse(p.adjust < 0.001, "***",
            ifelse(p.adjust < 0.01, "**",
              ifelse(p.adjust < 0.05, "*", "")
              )
            )
          ),
          position = position_dodge(width = 0.55),
          vjust = -1.0,
          size = 3.6,
          color = "#333333"
        ) +
        scale_color_manual(values = c("RES-CON" = "grey60", "SUS-CON" = "#E63946", "RES-SUS" = "#457B9D")) +
        scale_shape_manual(values = c("RES-CON" = 16, "SUS-CON" = 16, "RES-SUS" = 16)) +
        labs(title = paste("Phase-averaged contrasts:", metric_name),
             x = "Cage Change", y = "Estimate (difference)") +
        facet_grid(Phase ~ Sex, scales = "free_x") +
        theme_minimal(base_size = 12) +
        theme(panel.grid.minor = element_blank(), legend.position = "top",
              plot.title = element_text(face = "bold", hjust = 0.5)) +
        coord_flip()
      
      ggsave(file.path(dirs$plots_pub, paste0("panelA_forest_phaseAvg_with_RES-SUS_", run_scope, ".svg")),
             forest_phaseavg, width = 5, height = 9)
    }
  }
  
  # -------------------------------------------------
# Save to proper directories + generate manifest & inventory
# -------------------------------------------------

# 1. Move phase-average to correct folder
if (nrow(phase_avg_means_tbl) > 0) {
  openxlsx::write.xlsx(phase_avg_means_tbl, 
                       file.path(dirs$tables_phase, paste0("phaseAvg_emmeans_", run_scope, ".xlsx")), 
                       rowNames = FALSE)
}

if (nrow(phase_avg_contr_tbl) > 0) {
  openxlsx::write.xlsx(phase_avg_contr_tbl, 
                       file.path(dirs$tables_phase, paste0("phaseAvg_contrasts_", run_scope, ".xlsx")), 
                       rowNames = FALSE)
}

# 2. Save contrasts separately
if (nrow(emmeans_results) > 0) {
  contrasts_only <- emmeans_results %>%
    dplyr::filter(!is.na(contrast))
  
  if (nrow(contrasts_only) > 0) {
    openxlsx::write.xlsx(contrasts_only, 
                         file.path(dirs$tables_contr, paste0("contrasts_", run_scope, ".xlsx")), 
                         rowNames = FALSE)
  }
}

# 3. Generate inventory (data summary per slice)
inventory_list <- list()

for (Change in changes) {
  for (Sex in sexes) {
    for (Phase in phases) {
      dsub <- data %>%
        dplyr::filter(
          (!includeChange | Change == !!Change),
          (!includeSex    | Sex    == !!Sex),
          (!includePhase  | Phase  == !!Phase)
        )
      
      if (nrow(dsub) > 0) {
        inventory_list[[length(inventory_list) + 1]] <- data.frame(
          Change = as.character(Change),
          Sex = as.character(Sex),
          Phase = as.character(Phase),
          n_observations = nrow(dsub),
          n_animals = dplyr::n_distinct(dsub$AnimalNum),
          n_groups = dplyr::n_distinct(dsub$Group),
          n_timepoints = dplyr::n_distinct(dsub$HalfHourElapsed),
          mean_Movement = mean(dsub$Movement, na.rm = TRUE),
          sd_Movement = sd(dsub$Movement, na.rm = TRUE),
          mean_Proximity = mean(dsub$Proximity, na.rm = TRUE),
          sd_Proximity = sd(dsub$Proximity, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

inventory_df <- dplyr::bind_rows(inventory_list)

if (nrow(inventory_df) > 0) {
  openxlsx::write.xlsx(inventory_df, 
                       file.path(dirs$tables_inv, paste0("data_inventory_", run_scope, ".xlsx")), 
                       rowNames = FALSE)
}

# 4. Generate manifest (analysis metadata)
manifest <- data.frame(
  metric = metric_name,
  analysis_date = Sys.time(),
  run_scope = run_scope,
  includeChange = includeChange,
  includeSex = includeSex,
  includePhase = includePhase,
  use_group_time_interaction = use_group_time_interaction,
  min_animals = min_animals,
  min_groups = min_groups,
  min_obs = min_obs,
  min_time_lv = min_time_lv,
  n_models_fitted = nrow(summary_results) / length(unique(summary_results$Fixed_effect)),
  n_contrasts = nrow(emmeans_results),
  n_phase_avg_contrasts = nrow(phase_avg_contr_tbl),
  total_observations = nrow(data),
  total_animals = dplyr::n_distinct(data$AnimalNum),
  data_file = data_file,
  results_dir = results_dir,
  stringsAsFactors = FALSE
)

openxlsx::write.xlsx(manifest, 
                     file.path(dirs$tables_man, paste0("analysis_manifest_", metric_name, "_", run_scope, ".xlsx")), 
                     rowNames = FALSE)

cat(sprintf("\n✓ Additional tables saved:\n"))
cat(sprintf("  - Inventory: %s\n", dirs$tables_inv))
cat(sprintf("  - Manifest: %s\n", dirs$tables_man))
cat(sprintf("  - Phase average: %s\n", dirs$tables_phase))
cat(sprintf("  - Contrasts: %s\n", dirs$tables_contr))


  # -------------------------------------------------
  # Additional publication plots (volcano, panels A/B/C)
  # -------------------------------------------------
  
  # Load phase-average result tables
  phase_avg_contr_tbl <- openxlsx::read.xlsx(file.path(dirs$tables_emm, paste0("phaseAvg_contrasts_", run_scope, ".xlsx")))
  phase_avg_means_tbl <- openxlsx::read.xlsx(file.path(dirs$tables_emm, paste0("phaseAvg_emmeans_", run_scope, ".xlsx")))
  
  # Convert to factors
  phase_avg_contr_tbl <- phase_avg_contr_tbl %>%
    mutate(
      contrast = factor(contrast, levels = c("RES - CON", "SUS - CON", "RES - SUS")),
      Change = factor(Change),
      Sex = factor(Sex),
      Phase = factor(Phase, levels = c("Active", "Inactive"))
    )
  
  phase_avg_means_tbl <- phase_avg_means_tbl %>%
    mutate(
      Group = factor(Group, levels = c("CON","RES","SUS")),
      Phase = factor(Phase, levels = c("Active", "Inactive")),
      Change = factor(Change),
      Sex = factor(Sex)
    )
  
  # Volcano plot
plot_volcano_phaseavg <- function(df, n_labels = 10, pval_thresh = 0.05, effect_thresh = 1) {
    df <- df %>%
        mutate(
            sig = p.adjust < pval_thresh,
            large_effect = abs(estimate) > effect_thresh,
            label_flag = sig & large_effect,
            contrast = factor(contrast, levels = c("RES - CON", "SUS - CON", "RES - SUS")),
            Phase = factor(Phase, levels = c("Active", "Inactive")),
            Sex = factor(Sex)
        )
    
    # Calculate x-axis range based on data
    est_range <- range(df$estimate, na.rm = TRUE)
    est_min <- floor(est_range[1])
    est_max <- ceiling(est_range[2])
    x_breaks <- scales::pretty_breaks(n = 6)(c(est_min, est_max))
    
    ggplot(
        df,
        aes(
            x = estimate,
            y = -log10(p.adjust),
            color = contrast,
            shape = Phase
        )
    ) +
        geom_vline(
            xintercept = c(-effect_thresh, effect_thresh),
            linetype = "dotted",
            color = "gray50"
        ) +
        geom_hline(
            yintercept = -log10(0.05),
            linetype = "dashed",
            color = "gray70"
        ) +
        geom_point(size = 5, alpha = 0.8) +
        scale_color_manual(values = c(
            "RES - CON" = "#457B9D",
            "SUS - CON" = "#E63946",
            "RES - SUS" = "#C6C3BB"
        )) +
        scale_shape_manual(values = c("Active" = 16, "Inactive" = 1), name = "Phase") +
        facet_grid(. ~ Sex) +
        theme_minimal(base_size = 14) +
        theme(
            panel.grid.major.x = element_line(linewidth = 0.2),
            panel.grid.major.y = element_line(linewidth = 0.2),
            legend.position = "top"
        ) +
        scale_x_continuous(breaks = x_breaks, minor_breaks = NULL) +
        labs(
            x = "Effect size (estimate)",
            y = expression(-log[10](adjusted~p~value)),
            color = "Contrast",
            title = paste("Volcano Plot:", metric_name)
        )
}

  volcano_plot_phaseavg <- plot_volcano_phaseavg(phase_avg_contr_tbl)
  ggsave(file.path(dirs$plots_pub, paste0("volcano_phaseAvg_", run_scope, ".svg")), 
         volcano_plot_phaseavg, width = 8, height = 6)
  
  # Forest plot (additional)
  plot_forest_phaseavg <- function(df, contrast_levels = levels(df$contrast)) {
    sub <- df %>% filter(contrast %in% contrast_levels)
    ggplot(sub, aes(y = reorder(contrast, estimate), x = estimate, xmin = lwr, xmax = upr, color = contrast)) +
      geom_point(size = 3) +
      geom_errorbarh(height = 0.2) +
      facet_grid(Phase ~ Change) +
      theme_minimal(base_size = 14) +
      labs(x = "Estimate (phase avg, with 95% CI)", y = "Contrast", 
           title = paste("Forest Plot:", metric_name))
  }
  
  forest_plot_phaseavg <- plot_forest_phaseavg(phase_avg_contr_tbl)
  ggsave(file.path(dirs$plots_pub, paste0("forest_phaseAvg_", run_scope, ".svg")), 
         forest_plot_phaseavg, width = 9, height = 7)
  
  # Dot estimation plot
  plot_dot_estimation_phaseavg <- function(df) {
    ggplot(df, aes(x = Group, y = emmean_avg, color = Sex)) +
      geom_point(position = position_jitter(width = 0.2), size = 3) +
      geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
      facet_grid(Phase ~ Change) +
      theme_minimal(base_size = 14) +
      labs(x = "Group", y = "Estimated Phase Average", 
           title = paste("Dot/Estimation Plot:", metric_name))
  }
  
  dot_estimation_plot_phaseavg <- plot_dot_estimation_phaseavg(phase_avg_means_tbl)
  ggsave(file.path(dirs$plots_pub, paste0("dotEstimation_phaseAvg_", run_scope, ".svg")), 
         dot_estimation_plot_phaseavg, width = 8, height = 6)
  
  # -------------------------------------------------
  # Publication panels A/B/C
  # -------------------------------------------------
  theme_pub_small <- theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "top", 
          plot.title = element_text(face = "bold", hjust = 0.5))
  
  emm_df <- openxlsx::read.xlsx(emmeans_xlsx)
  parse_contrast <- function(x) {
    x <- gsub("\\|.*$", "", x); x <- gsub("\\s+", "", x); x <- tolower(x)
    parts <- strsplit(x, "-", fixed = TRUE)[[1]]
    if (length(parts) != 2) return(c(NA_character_, NA_character_))
    parts
  }
  pc <- t(vapply(emm_df$contrast, parse_contrast, FUN.VALUE = character(2)))
  colnames(pc) <- c("lhs","rhs")
  emm_df <- cbind(emm_df, pc) |>
    as.data.frame() |>
    dplyr::mutate(
      pair = dplyr::case_when(
        (lhs == "sus" & rhs == "con") | (lhs == "con" & rhs == "sus") ~ "SUS-CON",
        (lhs == "res" & rhs == "con") | (lhs == "con" & rhs == "res") ~ "RES-CON",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::filter(!is.na(pair), !is.na(Sex), !is.na(Change), !is.na(Phase)) |>
    dplyr::mutate(
      tcrit = qt(0.975, df = df),
      ymin = estimate - tcrit * SE,
      ymax = estimate + tcrit * SE,
      Sex = droplevels(factor(Sex)),
      ChangeCollapsed = if (!includeChange) "allChanges" else as.character(Change),
      Phase = factor(Phase, levels = c("Active","Inactive"))
    )
  
  emm_df <- emm_df %>%
    dplyr::mutate(
      XFacet = if (includeChange)
        factor(as.character(Change), levels = sort(unique(as.character(Change))))
      else
        factor("allChanges", levels = "allChanges")
    )
  
  # Panel A: Forest plot
  file_panelA <- file.path(dirs$plots_pub, paste0("panelA_forest_phase_", run_scope, ".svg"))
  if (nrow(emm_df) > 0 && length(unique(emm_df$Sex)) > 0) {
    forest_plot <- ggplot(
      data = emm_df,
      aes(x = XFacet, y = estimate, ymin = ymin, ymax = ymax, color = pair, shape = pair)
    ) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_pointrange(position = position_dodge(width = 0.5), linewidth = 0.6) +
      geom_text(aes(label = ifelse(pair == "SUS-CON" & p.adjust < 0.05, "★", "")),
                position = position_dodge(width = 0.5), vjust = -1.0, size = 3.8, color = "#333333") +
      scale_color_manual(values = c("RES-CON" = "grey70", "SUS-CON" = "#E63946")) +
      scale_shape_manual(values = c("RES-CON" = 16, "SUS-CON" = 17)) +
      labs(title = paste("EMM contrasts vs CON by Phase:", metric_name), 
           x = if (includeChange) "Cage change" else "allChanges", 
           y = "Estimate (difference vs CON)") +
      facet_grid(Phase ~ Sex, scales = "free_x") +
      theme_pub_small
    ggsave(file_panelA, forest_plot, width = 8.5, height = 4.8)
  }
  
  # Panel B: Ribbons
  cc_starts <- data %>% 
    dplyr::group_by(Change) %>% 
    dplyr::summarise(t0 = min(HalfHourElapsed, na.rm = TRUE), .groups = "drop")
  
  df_plot_B <- data %>% 
    dplyr::filter(
      (!includeChange | Change %in% unique(Change)), 
      (!includeSex | Sex %in% unique(Sex)), 
      (!includePhase | Phase %in% unique(Phase))
    )
  
  tc_df <- df_plot_B %>% 
    dplyr::left_join(cc_starts, by = "Change") %>% 
    dplyr::mutate(rel_time_h = (HalfHourElapsed - t0) / 2) %>%
    dplyr::group_by(Change, Sex, Group, rel_time_h) %>% 
    dplyr::summarise(
      mean_act = mean(.data[[metric_name]], na.rm = TRUE), 
      sd_act = sd(.data[[metric_name]], na.rm = TRUE), 
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Group = factor(Group, levels = c("CON","RES","SUS")), 
      Sex = droplevels(factor(Sex)), 
      Change = droplevels(factor(Change))
    )
  
  file_panelB <- file.path(dirs$plots_pub, paste0("panelB_ribbons_rel_phaseCollapsed_", run_scope, ".svg"))
  if (nrow(tc_df) > 0) {
    max_rel <- max(tc_df$rel_time_h, na.rm = TRUE)
    ribbons_plot <- ggplot(tc_df, aes(x = rel_time_h, y = mean_act, color = Group, fill = Group)) +
      geom_ribbon(aes(ymin = mean_act - sd_act, ymax = mean_act + sd_act), alpha = 0.22, colour = NA) +
      geom_line(size = 0.8) +
      scale_color_manual(values = c("CON"="#C6C3BB", "RES"="#457B9D", "SUS"="#E63946")) +
      scale_fill_manual(values = c("CON"="#C6C3BB", "RES"="#457B9D", "SUS"="#E63946")) +
      scale_x_continuous(breaks = scales::breaks_pretty(), limits = c(0, max_rel)) +
      labs(title = paste(metric_name, "over CC-relative time (mean ± SD)"), 
           x = "Time since CC start [h]", 
           y = paste(metric_name, "[a.u.]")) +
      facet_grid(Sex ~ Change, scales = "free_y") +
      theme_pub_small
    ggsave(file_panelB, ribbons_plot, width = 10.5, height = 4.6)
  }
  
  # Panel C: Heatmap
  df_plot_C <- data %>% 
    dplyr::filter(
      Group %in% c("CON","SUS"), 
      (!includeChange | Change %in% unique(Change)), 
      (!includeSex | Sex %in% unique(Sex)), 
      (!includePhase | Phase %in% unique(Phase))
    )
  
  phase_starts <- df_plot_C %>% 
    dplyr::group_by(Change, Sex, Phase) %>% 
    dplyr::summarise(t0_phase = min(HalfHourElapsed, na.rm = TRUE), .groups = "drop")
  
  diff_df <- df_plot_C %>% 
    dplyr::left_join(phase_starts, by = c("Change","Sex","Phase")) %>% 
    dplyr::mutate(rel_phase_h = (HalfHourElapsed - t0_phase) / 2) %>%
    dplyr::group_by(Change, Sex, Phase, Group, rel_phase_h) %>% 
    dplyr::summarise(mean_act = mean(.data[[metric_name]], na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Group, values_from = mean_act) %>% 
    dplyr::mutate(
      diff_SC = SUS - CON, 
      Sex = droplevels(factor(Sex)), 
      Change = droplevels(factor(Change)), 
      Phase = factor(Phase, levels = c("Active","Inactive"))
    ) %>%
    dplyr::filter(!is.na(diff_SC))
  
  file_panelC <- file.path(dirs$plots_pub, paste0("panelC_heat_phaseAligned_", run_scope, ".svg"))
  if (nrow(diff_df) > 0) {
    thresh <- stats::quantile(abs(diff_df$diff_SC), probs = 0.95, na.rm = TRUE)
    diff_df <- diff_df %>% dplyr::mutate(sig_bin = abs(diff_SC) >= as.numeric(thresh))
    heat_plot <- ggplot(diff_df, aes(x = rel_phase_h, y = factor(Change, levels = sort(unique(as.character(Change)))), fill = diff_SC)) +
      geom_tile() +
      scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c", midpoint = 0) +
      geom_point(data = subset(diff_df, sig_bin),
                 aes(x = rel_phase_h, y = factor(Change, levels = sort(unique(as.character(Change))))),
                 shape = 4, size = 0.8, color = "black", inherit.aes = FALSE) +
      labs(title = paste("SUS−CON difference:", metric_name), 
           x = "Time since phase start [h]", 
           y = "Cage change", 
           fill = paste("Δ", metric_name)) +
      facet_grid(Phase ~ Sex, scales = "free_y") +
      theme_pub_small
    ggsave(file_panelC, heat_plot, width = 10.5, height = 4.8)
  }
  
  # Combined plots
  if (exists("forest_plot")) {
    if (exists("ribbons_plot") && exists("heat_plot")) {
      combined_emm_abc <- (forest_plot / ribbons_plot / heat_plot) +
        patchwork::plot_layout(heights = c(4.8, 4.6, 4.8), guides = "collect") &
        theme(legend.position = "top")
      ggsave(file.path(dirs$plots_pub, paste0("Combined_EMM_ABC_", run_scope, ".svg")), 
             combined_emm_abc, width = 10.5, height = 14.2)
    } else if (exists("ribbons_plot")) {
      combined_emm_ab <- (forest_plot / ribbons_plot) +
        patchwork::plot_layout(heights = c(4.8, 4.6), guides = "collect") &
        theme(legend.position = "top")
      ggsave(file.path(dirs$plots_pub, paste0("Combined_EMM_AB_", run_scope, ".svg")), 
             combined_emm_ab, width = 10.5, height = 9.4)
    }
  }

  cat(sprintf("\n=== Analysis complete for %s ===\n", metric_name))
  return(list(dirs = dirs, metric = metric_name))
}

# =========================
# Cross-metric correlation analysis
# =========================

cat("\n=== CROSS-METRIC CORRELATION ANALYSIS ===\n")

# Create correlation directory
correlation_dir <- file.path(base_results_dir, "Correlation_Analysis")
dir_create_safe(correlation_dir)
dirs_corr <- list(
  plots = file.path(correlation_dir, "plots"),
  tables = file.path(correlation_dir, "tables")
)
invisible(lapply(dirs_corr, dir_create_safe))

# Prepare data with both metrics
corr_data <- data_filtered_agg %>%
  dplyr::select(Change, Sex, Phase, Group, AnimalNum, HalfHourElapsed, TH, Movement, Proximity) %>%
  dplyr::filter(!is.na(Movement), !is.na(Proximity)) %>%
  dplyr::mutate(
    Change = factor(Change),
    Sex = factor(Sex),
    Phase = factor(Phase, levels = c("Active", "Inactive")),
    Group = factor(Group, levels = c("CON", "RES", "SUS"))
  )

cat(sprintf("Correlation data: %d observations with both metrics\n", nrow(corr_data)))

# Overall correlation
overall_cor <- cor.test(corr_data$Movement, corr_data$Proximity, method = "pearson")
overall_spearman <- cor.test(corr_data$Movement, corr_data$Proximity, method = "spearman")

cat(sprintf("Overall Pearson correlation: r = %.3f, p = %.3e\n", 
            overall_cor$estimate, overall_cor$p.value))
cat(sprintf("Overall Spearman correlation: rho = %.3f, p = %.3e\n", 
            overall_spearman$estimate, overall_spearman$p.value))

# Correlation by Group, Phase, Sex, Change
corr_by_group <- corr_data %>%
  dplyr::group_by(Change, Sex, Phase, Group) %>%
  dplyr::summarise(
    n = n(),
    pearson_r = cor(Movement, Proximity, method = "pearson", use = "complete.obs"),
    spearman_rho = cor(Movement, Proximity, method = "spearman", use = "complete.obs"),
    .groups = "drop"
  )

# Calculate p-values separately (cleaner and faster)
corr_by_group$pearson_p <- NA_real_
corr_by_group$spearman_p <- NA_real_

for (i in 1:nrow(corr_by_group)) {
  d <- corr_data %>% 
    dplyr::filter(
      Change == corr_by_group$Change[i], 
      Sex == corr_by_group$Sex[i], 
      Phase == corr_by_group$Phase[i], 
      Group == corr_by_group$Group[i]
    )
  
  if (nrow(d) >= 3) {
    # Pearson test
    pearson_test <- tryCatch(
      cor.test(d$Movement, d$Proximity, method = "pearson"),
      error = function(e) NULL
    )
    if (!is.null(pearson_test)) {
      corr_by_group$pearson_p[i] <- pearson_test$p.value
    }
    
    # Spearman test (with exact = FALSE to avoid ties warning)
    spearman_test <- tryCatch(
      suppressWarnings(
        cor.test(d$Movement, d$Proximity, method = "spearman", exact = FALSE)
      ),
      error = function(e) NULL
    )
    if (!is.null(spearman_test)) {
      corr_by_group$spearman_p[i] <- spearman_test$p.value
    }
  }
}


# Save correlation statistics
openxlsx::write.xlsx(
  corr_by_group, 
  file.path(dirs_corr$tables, "correlation_by_group_phase_sex_change.xlsx"), 
  rowNames = FALSE
)

# Overall scatter plot with regression line
p_overall <- ggplot(corr_data, aes(x = Movement, y = Proximity)) +
  geom_point(aes(color = Group), alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", color = "black", linewidth = 1.2) +
  scale_color_manual(values = c("CON" = "#457B9D", "RES" = "#C6C3BB", "SUS" = "#E63946")) +
  labs(
    title = sprintf("Movement vs Proximity\nr = %.3f, p < %.3e", 
                    overall_cor$estimate, overall_cor$p.value),
    x = "Movement [a.u.]",
    y = "Proximity [a.u.]"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave(
  file.path(dirs_corr$plots, "correlation_overall.svg"),
  p_overall, width = 7, height = 6
)

# Scatter plot faceted by Group
p_by_group <- ggplot(corr_data, aes(x = Movement, y = Proximity, color = Group)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("CON" = "#457B9D", "RES" = "#C6C3BB", "SUS" = "#E63946")) +
  facet_wrap(~ Group, ncol = 3) +
  labs(
    title = "Movement vs Proximity by Group",
    x = "Movement [a.u.]",
    y = "Proximity [a.u.]"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave(
  file.path(dirs_corr$plots, "correlation_by_group.svg"),
  p_by_group, width = 10, height = 4
)

# Scatter plot faceted by Phase and Sex
p_phase_sex <- ggplot(corr_data, aes(x = Movement, y = Proximity, color = Group)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
  scale_color_manual(values = c("CON" = "#457B9D", "RES" = "#C6C3BB", "SUS" = "#E63946")) +
  facet_grid(Phase ~ Sex, scales = "free") +
  labs(
    title = "Movement vs Proximity by Phase and Sex",
    x = "Movement [a.u.]",
    y = "Proximity [a.u.]"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave(
  file.path(dirs_corr$plots, "correlation_by_phase_sex.svg"),
  p_phase_sex, width = 8, height = 7
)

# Correlation heatmap by Group/Phase/Sex
corr_heatmap_data <- corr_by_group %>%
  dplyr::mutate(
    label = sprintf("r=%.2f\np=%.3f", pearson_r, pearson_p),
    sig = ifelse(pearson_p < 0.05, "*", "")
  )

p_heatmap <- ggplot(
  corr_heatmap_data,
  aes(x = interaction(Phase, Sex), y = Group, fill = pearson_r)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f%s", pearson_r, sig)), size = 3.5, color = "black") +
  scale_fill_gradient2(
    low = "#2c7bb6", 
    mid = "white", 
    high = "#d7191c", 
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson r"
  ) +
  facet_wrap(~ Change, ncol = 2) +
  labs(
    title = "Correlation Heatmap: Movement vs Proximity",
    subtitle = "* indicates p < 0.05",
    x = "Phase × Sex",
    y = "Group"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave(
  file.path(dirs_corr$plots, "correlation_heatmap.svg"),
  p_heatmap, width = 10, height = 8
)

# Time-resolved correlation (sliding window)
if ("TH" %in% names(corr_data)) {
  cat("Computing time-resolved correlations...\n")
  
  # Calculate correlation for each time point by Group
  time_corr <- corr_data %>%
    dplyr::group_by(Change, Sex, Phase, Group, HalfHourElapsed) %>%
    dplyr::summarise(
      n = n(),
      pearson_r = {
        if (n() >= 3) {
          # Check for zero variance
          var_movement <- var(Movement, na.rm = TRUE)
          var_proximity <- var(Proximity, na.rm = TRUE)
          
          if (!is.na(var_movement) && !is.na(var_proximity) && 
              var_movement > 0 && var_proximity > 0) {
            cor(Movement, Proximity, method = "pearson", use = "complete.obs")
          } else {
            NA_real_
          }
        } else {
          NA_real_
        }
      },
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(pearson_r))
  
  # Plot time-resolved correlation
  p_time_corr <- ggplot(time_corr, aes(x = HalfHourElapsed/2, y = pearson_r, color = Group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(linewidth = 0.8, alpha = 0.8) +
    geom_point(size = 1.5, alpha = 0.6) +
    scale_color_manual(values = c("CON" = "#457B9D", "RES" = "#C6C3BB", "SUS" = "#E63946")) +
    facet_grid(Phase ~ Sex, scales = "free_x") +
    labs(
      title = "Time-resolved Correlation: Movement vs Proximity",
      x = "Time Elapsed [h]",
      y = "Pearson r"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  ggsave(
    file.path(dirs_corr$plots, "correlation_time_resolved.svg"),
    p_time_corr, width = 10, height = 7
  )
  
  # Save time-resolved correlation data
  openxlsx::write.xlsx(
    time_corr, 
    file.path(dirs_corr$tables, "correlation_time_resolved.xlsx"), 
    rowNames = FALSE
  )
  
  cat(sprintf("  - Time-resolved correlations: %d valid time points (out of %d total)\n", 
              nrow(time_corr), 
              nrow(corr_data %>% dplyr::group_by(Change, Sex, Phase, Group, HalfHourElapsed) %>% dplyr::summarise(n = n(), .groups = "drop"))))
}


# Correlation by animal (individual level) WITH STATISTICS
# Correlation by animal (individual level) WITH COMPREHENSIVE STATISTICS
animal_corr <- corr_data %>%
  dplyr::group_by(AnimalNum, Group, Sex, Change) %>%
  dplyr::summarise(
    n = n(),
    pearson_r = if(n() >= 3) cor(Movement, Proximity, method = "pearson", use = "complete.obs") else NA_real_,
    spearman_rho = if(n() >= 3) cor(Movement, Proximity, method = "spearman", use = "complete.obs") else NA_real_,
    .groups = "drop"
  ) %>%
  dplyr::filter(!is.na(pearson_r))

# Save animal-level correlations
openxlsx::write.xlsx(
  animal_corr, 
  file.path(dirs_corr$tables, "correlation_by_animal.xlsx"), 
  rowNames = FALSE
)

# Statistical tests: Compare correlations between groups
stat_results_list <- list()
all_posthoc_list <- list()

for (sex_val in unique(animal_corr$Sex)) {
  animal_corr_sex <- animal_corr %>% dplyr::filter(Sex == sex_val)
  
  if (nrow(animal_corr_sex) < 3) next
  
  # Test normality
  normality_test <- shapiro.test(animal_corr_sex$pearson_r)
  is_normal <- normality_test$p.value > 0.05
  
  # Descriptive stats per group
  group_stats <- animal_corr_sex %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(
      n = n(),
      mean = mean(pearson_r, na.rm = TRUE),
      sd = sd(pearson_r, na.rm = TRUE),
      .groups = "drop"
    )
  
  cat(sprintf("\n=== Sex: %s ===\n", sex_val))
  cat(sprintf("Normality test: p = %.4f (%s)\n", normality_test$p.value, ifelse(is_normal, "Normal", "Non-normal")))
  print(group_stats)
  
  # Choose appropriate test
  if (is_normal && length(unique(animal_corr_sex$Group)) > 2) {
    # ANOVA
    anova_result <- aov(pearson_r ~ Group, data = animal_corr_sex)
    anova_summary <- summary(anova_result)
    overall_p <- anova_summary[[1]]$`Pr(>F)`[1]
    F_stat <- anova_summary[[1]]$`F value`[1]
    df1 <- anova_summary[[1]]$Df[1]
    df2 <- anova_summary[[1]]$Df[2]
    test_type <- "ANOVA"
    
    cat(sprintf("ANOVA: F(%d,%d) = %.3f, p = %.4f\n", df1, df2, F_stat, overall_p))
    
    # ALWAYS run Tukey HSD (even if not significant, for completeness)
    tukey_result <- TukeyHSD(anova_result)
    posthoc_df <- as.data.frame(tukey_result$Group)
    posthoc_df$comparison <- rownames(posthoc_df)
    posthoc_df <- posthoc_df %>%
      dplyr::mutate(
        Sex = sex_val,
        test = "Tukey HSD",
        group1 = sapply(strsplit(comparison, "-"), function(x) trimws(x[1])),
        group2 = sapply(strsplit(comparison, "-"), function(x) trimws(x[2])),
        p_value = `p adj`,
        sig_label = dplyr::case_when(
          p_value < 0.001 ~ "***",
          p_value < 0.01  ~ "**",
          p_value < 0.05  ~ "*",
          TRUE ~ "ns"
        ),
        overall_sig = overall_p < 0.05
      ) %>%
      dplyr::select(Sex, test, comparison, group1, group2, diff, lwr, upr, p_value, sig_label, overall_sig)
    
    cat("\nPost-hoc (Tukey HSD):\n")
    print(posthoc_df %>% dplyr::select(comparison, diff, p_value, sig_label))
    
  } else {
    # Kruskal-Wallis (non-parametric)
    kw_result <- kruskal.test(pearson_r ~ Group, data = animal_corr_sex)
    overall_p <- kw_result$p.value
    chi_stat <- kw_result$statistic
    df_kw <- kw_result$parameter
    test_type <- "Kruskal-Wallis"
    
    cat(sprintf("Kruskal-Wallis: χ²(%d) = %.3f, p = %.4f\n", df_kw, chi_stat, overall_p))
    
    # ALWAYS run pairwise Wilcoxon (even if not significant, for completeness)
    pairwise_result <- pairwise.wilcox.test(
      animal_corr_sex$pearson_r, 
      animal_corr_sex$Group, 
      p.adjust.method = "holm",
      exact = FALSE
    )
    
    # Convert matrix to dataframe
    p_matrix <- pairwise_result$p.value
    group_names <- c(rownames(p_matrix), colnames(p_matrix)[ncol(p_matrix)])
    
    posthoc_rows <- list()
    for (i in 1:nrow(p_matrix)) {
      for (j in 1:ncol(p_matrix)) {
        if (!is.na(p_matrix[i, j])) {
          group1 <- rownames(p_matrix)[i]
          group2 <- colnames(p_matrix)[j]
          p_val <- p_matrix[i, j]
          
          posthoc_rows[[length(posthoc_rows) + 1]] <- data.frame(
            Sex = sex_val,
            test = "Wilcoxon",
            comparison = paste(group2, group1, sep = " - "),  # Match Tukey format
            group1 = group2,
            group2 = group1,
            p_value = p_val,
            sig_label = dplyr::case_when(
              p_val < 0.001 ~ "***",
              p_val < 0.01  ~ "**",
              p_val < 0.05  ~ "*",
              TRUE ~ "ns"
            ),
            overall_sig = overall_p < 0.05,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    posthoc_df <- dplyr::bind_rows(posthoc_rows)
    
    cat("\nPost-hoc (Pairwise Wilcoxon with Holm correction):\n")
    print(posthoc_df %>% dplyr::select(comparison, p_value, sig_label))
  }
  
  # Store overall test result
  overall_result <- data.frame(
    Sex = sex_val,
    test_type = test_type,
    overall_p = overall_p,
    overall_sig = overall_p < 0.05,
    n_animals = nrow(animal_corr_sex),
    normality_p = normality_test$p.value,
    stringsAsFactors = FALSE
  )
  
  stat_results_list[[sex_val]] <- overall_result
  all_posthoc_list[[sex_val]] <- posthoc_df
}

# Combine results
overall_stats <- dplyr::bind_rows(stat_results_list)
posthoc_stats <- dplyr::bind_rows(all_posthoc_list)

# Save statistical results
openxlsx::write.xlsx(
  overall_stats, 
  file.path(dirs_corr$tables, "animal_correlation_overall_stats.xlsx"), 
  rowNames = FALSE
)

openxlsx::write.xlsx(
  posthoc_stats, 
  file.path(dirs_corr$tables, "animal_correlation_posthoc_stats.xlsx"), 
  rowNames = FALSE
)

# Create summary statistics per group
animal_corr_summary <- animal_corr %>%
  dplyr::group_by(Group, Sex) %>%
  dplyr::summarise(
    n = n(),
    mean_r = mean(pearson_r, na.rm = TRUE),
    sd_r = sd(pearson_r, na.rm = TRUE),
    se_r = sd_r / sqrt(n),
    median_r = median(pearson_r, na.rm = TRUE),
    Q1_r = quantile(pearson_r, 0.25, na.rm = TRUE),
    Q3_r = quantile(pearson_r, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

openxlsx::write.xlsx(
  animal_corr_summary, 
  file.path(dirs_corr$tables, "animal_correlation_summary.xlsx"), 
  rowNames = FALSE
)

# Enhanced violin plot with post-hoc comparisons (CORRECTED for facets)
p_animal_corr <- ggplot(animal_corr, aes(x = Group, y = pearson_r, fill = Group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 2) +
  scale_fill_manual(values = c("CON" = "#457B9D", "RES" = "#C6C3BB", "SUS" = "#E63946")) +
  facet_wrap(~ Sex) +
  labs(
    title = "Individual Animal Correlations: Movement vs Proximity",
    x = "Group",
    y = "Pearson r (per animal)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# Add post-hoc brackets - CORRECTED: Create data frames with Sex for proper faceting
for (sex_val in unique(posthoc_stats$Sex)) {
  posthoc_sex <- posthoc_stats %>% dplyr::filter(Sex == sex_val)
  
  if (nrow(posthoc_sex) > 0) {
    # Get y range for this sex
    animal_corr_sex <- animal_corr %>% dplyr::filter(Sex == sex_val)
    max_y <- max(animal_corr_sex$pearson_r, na.rm = TRUE)
    min_y <- min(animal_corr_sex$pearson_r, na.rm = TRUE)
    y_range <- max_y - min_y
    
    for (i in 1:nrow(posthoc_sex)) {
      row <- posthoc_sex[i, ]
      
      # Map group names to x positions
      group_levels <- c("CON", "RES", "SUS")
      x1 <- which(group_levels == row$group1)
      x2 <- which(group_levels == row$group2)
      
      if (length(x1) > 0 && length(x2) > 0) {
        # Calculate y position
        y_pos <- max_y + 0.08 * y_range + (i - 1) * 0.12 * y_range
        
        # Color based on significance
        bracket_color <- ifelse(row$sig_label != "ns", "black", "gray60")
        bracket_alpha <- ifelse(row$sig_label != "ns", 1, 0.4)
        label_color <- ifelse(row$sig_label != "ns", "black", "gray50")
        
        # Create data frames with Sex column for faceting
        segment_horizontal <- data.frame(
          Sex = sex_val,
          x = x1, xend = x2,
          y = y_pos, yend = y_pos
        )
        
        segment_left <- data.frame(
          Sex = sex_val,
          x = x1, xend = x1,
          y = y_pos - 0.02, yend = y_pos
        )
        
        segment_right <- data.frame(
          Sex = sex_val,
          x = x2, xend = x2,
          y = y_pos - 0.02, yend = y_pos
        )
        
        label_df <- data.frame(
          Sex = sex_val,
          x = (x1 + x2) / 2,
          y = y_pos + 0.03,
          label = row$sig_label
        )
        
        # Add segments using geom_segment with data argument
        p_animal_corr <- p_animal_corr +
          geom_segment(data = segment_horizontal,
                       aes(x = x, xend = xend, y = y, yend = yend),
                       color = bracket_color, linewidth = 0.5, alpha = bracket_alpha,
                       inherit.aes = FALSE) +
          geom_segment(data = segment_left,
                       aes(x = x, xend = xend, y = y, yend = yend),
                       color = bracket_color, linewidth = 0.5, alpha = bracket_alpha,
                       inherit.aes = FALSE) +
          geom_segment(data = segment_right,
                       aes(x = x, xend = xend, y = y, yend = yend),
                       color = bracket_color, linewidth = 0.5, alpha = bracket_alpha,
                       inherit.aes = FALSE) +
          geom_text(data = label_df,
                    aes(x = x, y = y, label = label),
                    size = 4.5, fontface = "bold", color = label_color,
                    alpha = bracket_alpha, inherit.aes = FALSE)
      }
    }
  }
}

# Add caption with statistics
caption_text <- ""
for (sex_val in unique(overall_stats$Sex)) {
  sex_stat <- overall_stats %>% dplyr::filter(Sex == sex_val)
  posthoc_sex <- posthoc_stats %>% dplyr::filter(Sex == sex_val)
  
  if (nrow(sex_stat) > 0) {
    p_val <- sex_stat$overall_p[1]
    test_name <- sex_stat$test_type[1]
    
    caption_text <- paste0(
      caption_text,
      sex_val, ": ", test_name, " p = ", 
      ifelse(p_val < 0.001, "< 0.001", sprintf("%.3f", p_val))
    )
    
    if (sex_stat$overall_sig[1]) {
      caption_text <- paste0(caption_text, " (significant)")
      
      sig_comparisons <- posthoc_sex %>% dplyr::filter(sig_label != "ns")
      if (nrow(sig_comparisons) > 0) {
        caption_text <- paste0(caption_text, "; Post-hoc: ")
        for (j in 1:nrow(sig_comparisons)) {
          comp_row <- sig_comparisons[j, ]
          caption_text <- paste0(
            caption_text,
            comp_row$comparison, " p=",
            ifelse(comp_row$p_value < 0.001, "<0.001", sprintf("%.3f", comp_row$p_value)),
            comp_row$sig_label
          )
          if (j < nrow(sig_comparisons)) caption_text <- paste0(caption_text, ", ")
        }
      }
    } else {
      caption_text <- paste0(caption_text, " (ns)")
    }
    caption_text <- paste0(caption_text, "\n")
  }
}

p_animal_corr <- p_animal_corr +
  labs(caption = caption_text) +
  theme(
    plot.caption = element_text(hjust = 0, size = 9, face = "italic"),
    plot.margin = margin(10, 10, 30, 10)
  )

ggsave(
  file.path(dirs_corr$plots, "correlation_by_animal.svg"),
  p_animal_corr, width = 10, height = 7
)

# Bar plot version (CORRECTED for facets)
p_animal_corr_bars <- ggplot(animal_corr_summary, aes(x = Group, y = mean_r, fill = Group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_bar(stat = "identity", alpha = 0.7, width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_r - se_r, ymax = mean_r + se_r),
    width = 0.2,
    linewidth = 0.8
  ) +
  scale_fill_manual(values = c("CON" = "#457B9D", "RES" = "#C6C3BB", "SUS" = "#E63946")) +
  facet_wrap(~ Sex) +
  labs(
    title = "Mean Correlation by Group (± SE)",
    x = "Group",
    y = "Mean Pearson r",
    caption = caption_text
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.caption = element_text(hjust = 0, size = 9, face = "italic"),
    plot.margin = margin(10, 10, 30, 10)
  )

# Add brackets to bar plot
for (sex_val in unique(posthoc_stats$Sex)) {
  posthoc_sex <- posthoc_stats %>% dplyr::filter(Sex == sex_val)
  
  if (nrow(posthoc_sex) > 0) {
    # Get y range for bar plot
    summary_sex <- animal_corr_summary %>% dplyr::filter(Sex == sex_val)
    max_y_bar <- max(summary_sex$mean_r + summary_sex$se_r, na.rm = TRUE)
    y_range_bar <- max_y_bar - min(summary_sex$mean_r - summary_sex$se_r, na.rm = TRUE)
    
    for (i in 1:nrow(posthoc_sex)) {
      row <- posthoc_sex[i, ]
      
      group_levels <- c("CON", "RES", "SUS")
      x1 <- which(group_levels == row$group1)
      x2 <- which(group_levels == row$group2)
      
      if (length(x1) > 0 && length(x2) > 0) {
        y_pos_bar <- max_y_bar + 0.05 * y_range_bar + (i - 1) * 0.12 * y_range_bar
        
        bracket_color <- ifelse(row$sig_label != "ns", "black", "gray60")
        bracket_alpha <- ifelse(row$sig_label != "ns", 1, 0.4)
        label_color <- ifelse(row$sig_label != "ns", "black", "gray50")
        
        # Create data frames with Sex column
        segment_df <- data.frame(
          Sex = sex_val,
          x = x1, xend = x2,
          y = y_pos_bar, yend = y_pos_bar
        )
        
        label_df <- data.frame(
          Sex = sex_val,
          x = (x1 + x2) / 2,
          y = y_pos_bar + 0.02,
          label = row$sig_label
        )
        
        p_animal_corr_bars <- p_animal_corr_bars +
          geom_segment(data = segment_df,
                       aes(x = x, xend = xend, y = y, yend = yend),
                       color = bracket_color, linewidth = 0.5, alpha = bracket_alpha,
                       inherit.aes = FALSE) +
          geom_text(data = label_df,
                    aes(x = x, y = y, label = label),
                    size = 4.5, fontface = "bold", color = label_color,
                    alpha = bracket_alpha, inherit.aes = FALSE)
      }
    }
  }
}

ggsave(
  file.path(dirs_corr$plots, "correlation_by_animal_barplot.svg"),
  p_animal_corr_bars, width = 10, height = 7
)

cat("\n✓ Animal-level correlation statistics with facet-specific annotations saved\n")

cat(sprintf("  - Overall stats saved to: %s\n", file.path(dirs_corr$tables, "animal_correlation_overall_stats.xlsx")))
cat(sprintf("  - Post-hoc stats saved to: %s\n", file.path(dirs_corr$tables, "animal_correlation_posthoc_stats.xlsx")))
cat(sprintf("  - Summary stats saved to: %s\n", file.path(dirs_corr$tables, "animal_correlation_summary.xlsx")))



cat("\n=== CORRELATION ANALYSIS COMPLETE ===\n")
cat(sprintf("Results saved to: %s\n", correlation_dir))


# -------------------------------------------------
# Run analyses for both metrics
# -------------------------------------------------
cat("\n=== RUNNING ANALYSES FOR BOTH METRICS ===\n")

# Analysis for Movement
movement_results <- run_analysis_for_metric("Movement", data_filtered_agg, includeChange, includeSex, includePhase)

# Analysis for Proximity  
proximity_results <- run_analysis_for_metric("Proximity", data_filtered_agg, includeChange, includeSex, includePhase)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Both Movement and Proximity have been analyzed.\n")
cat(sprintf("Movement results directory: %s\n", movement_results$dirs$models))
cat(sprintf("Proximity results directory: %s\n", proximity_results$dirs$models))












# =========================
# Edge analysis for current metric
# This should be added at the END of run_analysis_for_metric() function
# OR run separately after both metrics are processed
# =========================

run_edge_analysis_for_metric <- function(metric_name, data, dirs, includeChange, includeSex, includePhase) {
  
  cat(sprintf("\n=== Running edge analysis for %s ===\n", metric_name))
  
  # First, create edge_animal_means_collapsed using current metric
  edges_all_segments <- data %>%
    dplyr::mutate(
      seg_id = dplyr::case_when(
        Phase == "Active"   ~ as.integer(ConsecActive),
        Phase == "Inactive" ~ as.integer(ConsecInactive),
        TRUE ~ NA_integer_
      )
    ) %>% 
    dplyr::filter(!is.na(seg_id), seg_id >= 1L)
  
  seg_bounds <- edges_all_segments %>%
    dplyr::group_by(Change, Sex, Phase, AnimalNum, seg_id) %>%
    dplyr::summarise(
      seg_h0 = min(HalfHourElapsed, na.rm = TRUE),
      seg_h1 = max(HalfHourElapsed, na.rm = TRUE),
      .groups = "drop"
    )
  
  edges_labeled_all <- edges_all_segments %>%
    dplyr::left_join(seg_bounds, by = c("Change","Sex","Phase","AnimalNum","seg_id")) %>%
    dplyr::mutate(
      edge_bin = dplyr::case_when(
        HalfHourElapsed >= seg_h0 & HalfHourElapsed <= seg_h0 + 3L ~ "start4",
        HalfHourElapsed >= seg_h1 - 3L & HalfHourElapsed <= seg_h1   ~ "end4",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(edge_bin)) %>%
    dplyr::select(Change, Sex, Phase, Group, AnimalNum, seg_id, HalfHourElapsed, 
                  ActivityEdge = all_of(metric_name), edge_bin)  # KEY CHANGE: rename metric to ActivityEdge
  
  edge_animal_means_all <- edges_labeled_all %>%
    dplyr::group_by(Change, Sex, Phase, Group, AnimalNum, seg_id, edge_bin) %>%
    dplyr::summarise(ActivityEdge = mean(ActivityEdge, na.rm = TRUE), .groups = "drop")
  
  edge_animal_means_collapsed <- edge_animal_means_all %>%
    dplyr::group_by(Change, Sex, Phase, Group, AnimalNum, edge_bin) %>%
    dplyr::summarise(ActivityEdge = mean(ActivityEdge, na.rm = TRUE), .groups = "drop")
  
  # EMMeans options
  emm_options(
    pbkrtest.limit = 100000,
    lmerTest.limit = 100000,
    lmer.df = "kenward-roger"
  )
  
  # Canonical levels
  lvl_change <- sort(unique(as.character(edge_animal_means_collapsed$Change)))
  lvl_phase  <- c("Active", "Inactive")
  lvl_group  <- c("CON", "RES", "SUS")
  lvl_sex    <- sort(unique(as.character(edge_animal_means_collapsed$Sex)))
  lvl_edge   <- c("start4", "end4")
  
  # Prepare data
  edge_dat <- edge_animal_means_collapsed %>%
    dplyr::transmute(
      ActivityEdge = as.numeric(ActivityEdge),
      Change_f = factor(as.character(Change), levels = lvl_change),
      Change   = as.character(Change),
      Sex_f    = factor(as.character(Sex), levels = lvl_sex),
      Sex      = as.character(Sex),
      Phase_f  = factor(as.character(Phase), levels = lvl_phase),
      Phase    = as.character(Phase),
      Group_f  = factor(as.character(Group), levels = lvl_group),
      Group    = as.character(Group),
      AnimalNum = as.factor(AnimalNum),
      edge_bin  = factor(as.character(edge_bin), levels = lvl_edge)
    )
  
  edge_start <- edge_dat %>% dplyr::filter(edge_bin == "start4")
  edge_end   <- edge_dat %>% dplyr::filter(edge_bin == "end4")
  
  # Mixed-model fitting function
  fit_emm_subset <- function(df) {
    if (dplyr::n_distinct(df$Group_f) < 2 ||
        dplyr::n_distinct(df$AnimalNum) < 2) {
      return(NULL)
    }
    ok_change <- dplyr::n_distinct(df$Change_f) >= 2
    ok_phase  <- dplyr::n_distinct(df$Phase_f)  >= 2
    form_try <- if (ok_change && ok_phase) {
      ActivityEdge ~ Group_f + Phase_f + Change_f +
        Group_f:Phase_f + Group_f:Change_f + (1 | AnimalNum)
    } else {
      ActivityEdge ~ Group_f + Phase_f + Change_f + (1 | AnimalNum)
    }
    m <- quiet_singular(
      lmerTest::lmer(
        form_try,
        data = df,
        control = lme4::lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 2e5)
        )
      )
    )
    emmeans::emmeans(m, specs = ~ Group_f | Phase_f + Change_f)
  }
  
  # Converter function
  emm_to_df <- function(emm_obj, sex_label, edge_lab) {
    smry <- try(summary(emm_obj, infer = TRUE), silent = TRUE)
    if (inherits(smry, "try-error")) {
      smry <- emm_obj
    }
    
    df_emm <- NULL
    if (is.data.frame(smry)) {
      df_emm <- smry
    } else if (is.list(smry)) {
      df_emm <- try(as.data.frame(smry, stringsAsFactors = FALSE), silent = TRUE)
      if (inherits(df_emm, "try-error") || !is.data.frame(df_emm)) {
        df_emm <- tibble::tibble(emmean = NA_real_)
      }
    } else {
      df_emm <- tibble::tibble(emmean = NA_real_)
    }
    
    if (!("emmean" %in% names(df_emm)) && ("lsmean" %in% names(df_emm))) {
      df_emm$emmean <- df_emm$lsmean
    }
    if (!("lower.CL" %in% names(df_emm)) && ("asymp.LCL" %in% names(df_emm))) {
      df_emm$lower.CL <- df_emm$asymp.LCL
    }
    if (!("upper.CL" %in% names(df_emm)) && ("asymp.UCL" %in% names(df_emm))) {
      df_emm$upper.CL <- df_emm$asymp.UCL
    }
    
    if ("Phase_f" %in% names(df_emm) && !("Phase" %in% names(df_emm))) {
      df_emm$Phase <- as.character(df_emm$Phase_f)
    }
    if ("Change_f" %in% names(df_emm) && !("Change" %in% names(df_emm))) {
      df_emm$Change <- as.character(df_emm$Change_f)
    }
    if ("Group_f" %in% names(df_emm) && !("Group" %in% names(df_emm))) {
      df_emm$Group <- as.character(df_emm$Group_f)
    }
    
    if (!("Phase" %in% names(df_emm)))  df_emm$Phase  <- NA_character_
    if (!("Change" %in% names(df_emm))) df_emm$Change <- NA_character_
    if (!("Group" %in% names(df_emm)))  df_emm$Group  <- NA_character_
    
    df_emm$Sex <- sex_label
    
    if (!("emmean" %in% names(df_emm)))   df_emm$emmean   <- NA_real_
    if (!("SE" %in% names(df_emm)))       df_emm$SE       <- NA_real_
    if (!("df" %in% names(df_emm)))       df_emm$df       <- NA_real_
    if (!("lower.CL" %in% names(df_emm))) df_emm$lower.CL <- NA_real_
    if (!("upper.CL" %in% names(df_emm))) df_emm$upper.CL <- NA_real_
    
    em_vec  <- suppressWarnings(as.numeric(df_emm[["emmean"]]))
    se_vec  <- suppressWarnings(as.numeric(df_emm[["SE"]]))
    df_vec  <- suppressWarnings(as.numeric(df_emm[["df"]]))
    lwr_vec <- suppressWarnings(as.numeric(df_emm[["lower.CL"]]))
    upr_vec <- suppressWarnings(as.numeric(df_emm[["upper.CL"]]))
    
    out <- tibble::tibble(
      Change   = factor(as.character(df_emm[["Change"]]), levels = lvl_change),
      Phase    = factor(as.character(df_emm[["Phase"]]),  levels = lvl_phase),
      Group    = factor(as.character(df_emm[["Group"]]),  levels = lvl_group),
      Sex      = factor(as.character(df_emm[["Sex"]]),    levels = lvl_sex),
      emmean   = em_vec,
      SE       = se_vec,
      df       = df_vec,
      lwr      = lwr_vec,
      upr      = upr_vec,
      edge_bin = factor(edge_lab, levels = lvl_edge)
    )
    out
  }
  
  # Build EMMeans per sex
  emm_list_start <- list()
  emm_list_end   <- list()
  for (sx in levels(edge_dat$Sex_f)) {
    cat("\n---- Processing Sex:", sx, "----\n")
    ds <- edge_start %>% dplyr::filter(Sex == !!sx)
    de <- edge_end   %>% dplyr::filter(Sex == !!sx)
    
    es <- fit_emm_subset(ds)
    ee <- fit_emm_subset(de)
    
    if (!is.null(es)) {
      emm_list_start[[sx]] <- emm_to_df(es, sx, "start4")
    }
    if (!is.null(ee)) {
      emm_list_end[[sx]] <- emm_to_df(ee, sx, "end4")
    }
  }
  
  df_emm_start <- if (length(emm_list_start)) {
    dplyr::bind_rows(emm_list_start)
  } else {
    tibble::tibble(
      Change   = factor(character(), levels = lvl_change),
      Phase    = factor(character(), levels = lvl_phase),
      Group    = factor(character(), levels = lvl_group),
      Sex      = factor(character(), levels = lvl_sex),
      emmean   = double(),
      SE       = double(),
      df       = double(),
      lwr      = double(),
      upr      = double(),
      edge_bin = factor(character(), levels = lvl_edge)
    )
  }
  
  df_emm_end <- if (length(emm_list_end)) {
    dplyr::bind_rows(emm_list_end)
  } else {
    tibble::tibble(
      Change   = factor(character(), levels = lvl_change),
      Phase    = factor(character(), levels = lvl_phase),
      Group    = factor(character(), levels = lvl_group),
      Sex      = factor(character(), levels = lvl_sex),
      emmean   = double(),
      SE       = double(),
      df       = double(),
      lwr      = double(),
      upr      = double(),
      edge_bin = factor(character(), levels = lvl_edge)
    )
  }
  
  # Normalize function
  normalize_emm_df <- function(x) {
    dplyr::mutate(
      x,
      Change   = factor(as.character(.data$Change), levels = lvl_change),
      Phase    = factor(as.character(.data$Phase),  levels = lvl_phase),
      Group    = factor(as.character(.data$Group),  levels = lvl_group),
      Sex      = factor(as.character(.data$Sex),    levels = lvl_sex),
      edge_bin = factor(as.character(.data$edge_bin), levels = lvl_edge),
      emmean   = as.numeric(.data$emmean),
      SE       = as.numeric(.data$SE),
      df       = as.numeric(.data$df),
      lwr      = as.numeric(.data$lwr),
      upr      = as.numeric(.data$upr)
    )
  }
  
  df_emm_start <- normalize_emm_df(df_emm_start)
  df_emm_end   <- normalize_emm_df(df_emm_end)
  df_emm_all   <- normalize_emm_df(dplyr::bind_rows(df_emm_start, df_emm_end))
  
  run_scope <- scope_all(includeChange, includeSex, includePhase)
  
  # Save table
  openxlsx::write.xlsx(
    df_emm_all,
    file = file.path(
      dirs$tables_emm,
      paste0("edge_EMMeans_linearPred_", run_scope, ".xlsx")
    ),
    rowNames = FALSE
  )
  
  # Plotting
  facet_formula_all <- Phase + Sex ~ Change
  pal_grp <- c("CON"="#457B9D","RES"="#C6C3BB","SUS"="#E63946")
  
  # START plot
  df_emm_start$Group <- factor(df_emm_start$Group, levels = c("CON","RES","SUS"))
  plot_start_all <- ggplot(df_emm_start, aes(x = .data$Group, y = .data$emmean)) +
    geom_pointrange(
      aes(ymin = .data$lwr, ymax = .data$upr, color = .data$Group),
      linewidth = 0.9, fatten = 1.8
    ) +
    geom_point(
      aes(fill = .data$Group),
      shape = 21, color = "white", size = 3.4, stroke = 1.0
    ) +
    scale_color_manual(values = pal_grp, guide = "none") +
    scale_fill_manual(values = pal_grp, name = "Group") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    facet_grid(facet_formula_all, scales = "free_y") +
    labs(
      title = paste("Edge START (first 2h) —", metric_name),
      x = "Group", y = "Linear Prediction (EMM)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "top",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  ggsave(
    file.path(
      dirs$edgeAggAll,
      paste0("EMM_edge_START_all_", run_scope, ".svg")
    ),
    plot_start_all,
    width = 8,
    height = 5.2
  )
  
  # END plot
  df_emm_end$Group <- factor(df_emm_end$Group, levels = c("CON","RES","SUS"))
  plot_end_all <- ggplot(df_emm_end, aes(x = .data$Group, y = .data$emmean)) +
    geom_pointrange(
      aes(ymin = .data$lwr, ymax = .data$upr, color = .data$Group),
      linewidth = 0.9, fatten = 1.8
    ) +
    geom_point(
      aes(fill = .data$Group),
      shape = 21, color = "white", size = 3.4, stroke = 1.0
    ) +
    scale_color_manual(values = pal_grp, guide = "none") +
    scale_fill_manual(values = pal_grp, name = "Group") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    facet_grid(facet_formula_all, scales = "free_y") +
    labs(
      title = paste("Edge END (last 2h) —", metric_name),
      x = "Group", y = "Linear Prediction (EMM)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "top",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  ggsave(
    file.path(
      dirs$edgeAggAll,
      paste0("EMM_edge_END_all_", run_scope, ".svg")
    ),
    plot_end_all,
    width = 8,
    height = 5.2
  )
  
  # Trajectory plots
  make_traj <- function(df_emm, title) {
    df_emm <- normalize_emm_df(df_emm)
    df_emm$Change <- factor(as.character(df_emm$Change), levels = lvl_change)
    ggplot(
      df_emm,
      aes(x = .data$Change, y = .data$emmean, color = .data$Group, group = .data$Group)
    ) +
      geom_errorbar(
        aes(ymin = .data$lwr, ymax = .data$upr),
        width = 0.15,
        alpha = 0.9,
        position = position_dodge2(width = 0.25, padding = 0.15)
      ) +
      geom_point(
        size = 2.6,
        position = position_dodge2(width = 0.25, padding = 0.15)
      ) +
      geom_line(
        linewidth = 0.8,
        position = position_dodge2(width = 0.25, padding = 0.15)
      ) +
      scale_color_manual(values = pal_grp) +
      facet_grid(Phase ~ Sex, scales = "free_y") +
      labs(title = title, x = "Cage Change", y = "Linear Prediction (EMM)") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5))
  }
  
  plot_start_traj <- make_traj(
    df_emm_start,
    paste("START (first 2h) —", metric_name)
  )
  ggsave(
    file.path(
      dirs$edgeAggAll,
      paste0("EMM_edge_START_trajectory_", run_scope, ".svg")
    ),
    plot_start_traj,
    width = 8,
    height = 5.5
  )
  
  plot_end_traj <- make_traj(
    df_emm_end,
    paste("END (last 2h) —", metric_name)
  )
  ggsave(
    file.path(
      dirs$edgeAggAll,
      paste0("EMM_edge_END_trajectory_", run_scope, ".svg")
    ),
    plot_end_traj,
    width = 8,
    height = 5.5
  )
  
  # Compute edge statistics
  compute_edge_stats <- function(df_edge, edge_lab) {
    out_means <- list()
    out_contr <- list()
    
    for (sx in levels(df_edge$Sex_f)) {
      d_sx <- df_edge %>% dplyr::filter(Sex == !!sx)
      if (nrow(d_sx) < 1) next
      
      emm_obj <- fit_emm_subset(d_sx)
      if (is.null(emm_obj)) next
      
      emm_tbl <- try(summary(emm_obj, infer = TRUE), silent = TRUE)
      if (!inherits(emm_tbl, "try-error")) {
        emm_df <- as.data.frame(emm_tbl)
        names(emm_df) <- gsub("^lsmean$", "emmean", names(emm_df))
        names(emm_df) <- sub("^asymp\\.LCL$", "lower.CL", names(emm_df))
        names(emm_df) <- sub("^asymp\\.UCL$", "upper.CL", names(emm_df))
        
        emm_df$Sex      <- sx
        emm_df$edge_bin <- edge_lab
        if ("Group_f" %in% names(emm_df)) emm_df$Group <- as.character(emm_df$Group_f)
        if ("Phase_f" %in% names(emm_df)) emm_df$Phase <- as.character(emm_df$Phase_f)
        if ("Change_f" %in% names(emm_df)) emm_df$Change <- as.character(emm_df$Change_f)
        
        out_means[[paste(edge_lab, sx, sep="_")]] <- emm_df %>%
          dplyr::select(Change, Phase, Sex, Group, emmean, SE, df, lower.CL, upper.CL, edge_bin)
      }
      
      contr_pairs <- try(
        summary(pairs(emm_obj, by = c("Phase_f","Change_f")), infer = TRUE, adjust = "holm"),
        silent = TRUE
      )
      if (!inherits(contr_pairs, "try-error")) {
        contr_df <- as.data.frame(contr_pairs)
        if ("Phase_f" %in% names(contr_df)) contr_df$Phase  <- as.character(contr_df$Phase_f)
        if ("Change_f" %in% names(contr_df)) contr_df$Change <- as.character(contr_df$Change_f)
        
        contr_df$Sex      <- sx
        contr_df$edge_bin <- edge_lab
        if (!("p.value" %in% names(contr_df))) {
          if (all(c("t.ratio","df") %in% names(contr_df))) {
            contr_df$p.value <- 2 * pt(-abs(contr_df$t.ratio), df = contr_df$df)
          } else {
            contr_df$p.value <- NA_real_
          }
        }
        contr_df$p.holm <- p.adjust(contr_df$p.value, "holm")
        contr_df$p.BH   <- p.adjust(contr_df$p.value, "BH")
        
        out_contr[[paste(edge_lab, sx, sep="_")]] <- contr_df %>%
          dplyr::select(Change, Phase, Sex, edge_bin, contrast, estimate, SE, df, t.ratio, p.value, p.holm, p.BH, lower.CL, upper.CL)
      }
    }
    
    list(
      means = dplyr::bind_rows(out_means),
      contr = dplyr::bind_rows(out_contr)
    )
  }
  
  stats_start <- compute_edge_stats(edge_start, "start4")
  stats_end   <- compute_edge_stats(edge_end,   "end4")
  
  edge_emm_means    <- dplyr::bind_rows(stats_start$means, stats_end$means) %>%
    dplyr::mutate(
      Change = factor(Change, levels = lvl_change),
      Phase  = factor(Phase,  levels = lvl_phase),
      Sex    = factor(Sex,    levels = lvl_sex),
      Group  = factor(Group,  levels = lvl_group),
      edge_bin = factor(edge_bin, levels = lvl_edge)
    )
  
  edge_emm_contrasts <- dplyr::bind_rows(stats_start$contr, stats_end$contr) %>%
    dplyr::mutate(
      Change = factor(Change, levels = lvl_change),
      Phase  = factor(Phase,  levels = lvl_phase),
      Sex    = factor(Sex,    levels = lvl_sex),
      edge_bin = factor(edge_bin, levels = lvl_edge)
    )
  
  edge_stats_inventory <- edge_emm_contrasts %>%
    dplyr::mutate(sig_holm = p.holm < 0.05, sig_BH = p.BH < 0.05)
  
  # Save edge statistics
  edge_means_xlsx <- file.path(dirs$tables_emm, paste0("edge_emm_means_", run_scope, ".xlsx"))
  edge_contr_xlsx <- file.path(dirs$tables_emm, paste0("edge_emm_contrasts_", run_scope, ".xlsx"))
  edge_contr_csv  <- file.path(dirs$tables_emm, paste0("edge_stats_inventory_", run_scope, ".csv"))
  
  openxlsx::write.xlsx(edge_emm_means, edge_means_xlsx, rowNames = FALSE)
  openxlsx::write.xlsx(edge_emm_contrasts, edge_contr_xlsx, rowNames = FALSE)
  readr::write_csv(edge_stats_inventory, edge_contr_csv)
  
  cat(sprintf("\n=== Edge analysis complete for %s ===\n", metric_name))
}

# -------------------------------------------------
# Call edge analysis for both metrics
# -------------------------------------------------

# Run edge analysis for Movement
run_edge_analysis_for_metric("Movement", data_filtered_agg, movement_results$dirs, includeChange, includeSex, includePhase)

# Run edge analysis for Proximity
run_edge_analysis_for_metric("Proximity", data_filtered_agg, proximity_results$dirs, includeChange, includeSex, includePhase)
