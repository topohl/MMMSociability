###############################################
# Mixed models with random time slopes + fixed-grid phase-edge windows
# Modified to analyze both Movement and Proximity metrics
###############################################

suppressPackageStartupMessages({
  packages <- c(
    "ggplot2","dplyr","tidyr","gridExtra","lme4",
    "lmerTest","cowplot","emmeans","Matrix","nlme",
    "tcltk","openxlsx","readr","stringr","patchwork","scales",
    "broom.mixed","pbkrtest","rlang","tibble","MuMIn",
    "mgcv","pROC", "reformulas",
    "parallel","doParallel","foreach"
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

# Parallel model fitting: set to 1 to disable, or e.g. parallel::detectCores(logical=FALSE) - 1
n_cores <- 1L

# ─────────────────────────────────────────────────────────────────────
# Biologically prioritized contrasts (subset-level multiple-testing)
# ─────────────────────────────────────────────────────────────────────
# BIOLOGICAL RATIONALE:
#   SUS-CON: Tests susceptible vs control (primary hypothesis about stress response)
#   SUS-RES: Tests susceptible vs resilient (secondary: examines individual differences)
#   RES-CON: Not primary (exploratory; tests if resilience differs from control)
#
# PURPOSE: Reduce Type I/II error trade-offs by focusing p-value correction on
#          the 2 a priori contrasts, while maintaining exploratory significance
#          columns (family-level and global BH/Holm) for transparency.
#
# RECOMMENDATION: Report primary p-values in main text; family and global in tables
#                 to show sensitivity and full discovery landscape.
# ─────────────────────────────────────────────────────────────────────
primary_contrast_set <- c("SUS-CON", "SUS-RES")

# Contrast configuration: centralized definitions to avoid duplication
# Change these in one place to update all contrast references throughout the script
contrast_config <- list(
  display_names = c("RES - CON", "SUS - CON", "SUS - RES"),
  vectors = list(
    c(-1, 1, 0),   # RES - CON
    c(-1, 0, 1),   # SUS - CON
    c(0, -1, 1)    # SUS - RES
  ),
  colors = c(
    "RES - CON" = "#3d3b6e",
    "SUS - CON" = "#e63947",
    "SUS - RES" = "#C6C3BB"
  ),
  description = list(
    "RES - CON" = "Resilient vs Control",
    "SUS - CON" = "Susceptible vs Control (primary)",
    "SUS - RES" = "Susceptible vs Resilient (primary)"
  )
)

# Centralized color palettes — update here to change colors everywhere
group_colors <- c("CON" = "#3d3b6e", "RES" = "#C6C3BB", "SUS" = "#e63947")
pair_colors  <- c("RES-CON" = "#3d3b6e", "SUS-CON" = "#e63947", "SUS-RES" = "#C6C3BB")

# Suppress specific emmeans notes for large datasets
options(
  emmeans = list(
    pbkrtest.limit = 20000,  # Increase limit for Kenward-Roger df adjustments
    lmerTest.limit = 20000,  # Increase limit for Satterthwaite df adjustments
    msg.interaction = FALSE  # Suppress interaction involvement warnings
  )
)

# -------------------------------------------------
# Base Paths
# -------------------------------------------------
#base_results_dir <- "D:/MMMSociability/statistics/lme/"
base_results_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/statistics/gamm_new/"
if (!dir.exists(base_results_dir)) dir.create(base_results_dir, recursive = TRUE)
dir_create_safe <- function(path) { if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE) }

# -------------------------------------------------
# Helpers
# -------------------------------------------------
scope_all <- function(includeChange, includeSex, includePhase) paste0(
  if (includeChange) "byCC"    else "allCC",    "_",
  if (includeSex)    "bySex"   else "allSex",   "_",
  if (includePhase)  "byPhase" else "allPhase"
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
# LOAD DATA & REFERENCES
# -------------------------------------------------
cat("Loading data and SUS reference list...\n")
data_file <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability/processed_data/data_lme_format/data_filtered_agg.csv"
data_file_activity <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Raw Data/Behavior/RFID/BatchAnalysis/processed_data/aggregated/data_filtered_agg_new.csv"
sus_file <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/sus_animals.csv"

if (all(file.exists(data_file, data_file_activity, sus_file))) {
  
  # 1. LOAD DATASETS
  df_primary <- read_csv(data_file, show_col_types = FALSE)
  df_activity <- read_csv(data_file_activity, show_col_types = FALSE)
  sus_ids_raw <- read_csv(sus_file, col_names = "AnimalNum", show_col_types = FALSE)

  # 2. STANDARDIZE ANIMAL IDs
  # Removes leading zeros to ensure "0003" matches "3"
  clean_id <- function(x) {
    x_chr <- as.character(x)
    x_num <- suppressWarnings(as.numeric(x_chr))
    out <- ifelse(is.na(x_num), x_chr, as.character(x_num))
    as.character(out)
  }
  
  df_primary  <- df_primary  %>% mutate(AnimalNum = clean_id(AnimalNum))
  df_activity <- df_activity %>% mutate(AnimalNum = clean_id(AnimalNum))
  sus_list    <- clean_id(sus_ids_raw$AnimalNum)

  # 3. REASSIGN GROUP (Preserving CON)
  # Uses the sus_animals.csv as the ground truth for SUS/RES
  df_primary <- df_primary %>%
    mutate(Group = case_when(
      Group == "CON"            ~ "CON",  
      AnimalNum %in% sus_list   ~ "SUS",  
      TRUE                      ~ "RES"   
    ))

  # 4. DEDUPLICATE PRIMARY DATA
  # Fixes the row duplication caused by the 0003 vs 3 issue
  df_primary <- df_primary %>%
    distinct(Batch, Change, Sex, Phase, Group, AnimalNum, HalfHourElapsed, .keep_all = TRUE)

  # 5. PERFORM JOIN
  # We use 'distinct' on the activity data to make sure we don't have 
  # multiple rows for the same animal/time that would cause duplication.
  join_keys <- c("Batch", "AnimalNum", "Phase", "HalfHourElapsed")

  data_filtered_agg <- df_primary %>%
    left_join(
      df_activity %>% 
        # CRITICAL: Select only unique combinations of our keys
        group_by(Batch, AnimalNum, Phase, HalfHourElapsed) %>%
        summarise(ActivityIndex = mean(ActivityIndex, na.rm = TRUE), .groups = 'drop'), 
      by = join_keys
    )

  # 6. FINAL LOGGING
  cat(sprintf("Success! Merged dataset: %d rows\n", nrow(data_filtered_agg)))
  cat("Missing ActivityIndex values:", sum(is.na(data_filtered_agg$ActivityIndex)), "\n")
  cat("Final Group distribution:\n")
  print(table(data_filtered_agg$Group))
  
} else {
  stop("Check paths: One or more files (data, activity, or sus_list) are missing.")
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

compute_phase_avg_for_model <- function(model, Change, Sex, Phase) {
  
  # 1. Estimate marginal means and contrasts
  # Uses emmeans with correct Satterthwaite degrees of freedom
  emm <- emmeans::emmeans(model, specs = ~ Group)
  
  # Contrasts in predefined order from contrast_config
  # Create contrast method list dynamically from centralized configuration
  contrast_method <- setNames(
    contrast_config$vectors,
    contrast_config$display_names
  )
  contr_res <- emmeans::contrast(emm, method = contrast_method, adjust = "none")
  
  # 2. Total variance estimation for Cohen's d calculation
  # In mixed models: Total SD = sqrt(random variance + residual variance)
  var_comp <- as.data.frame(lme4::VarCorr(model))
  total_variance <- sum(var_comp$vcov) + sigma(model)^2
  total_sd <- sqrt(total_variance)
  
  # 3. Helper function for robust DataFrame formatting
  # Handles variable column naming from emmeans output
  safe_format <- function(emm_object) {
    # infer=TRUE forces calculation of confidence intervals (LCL/UCL)
    df <- as.data.frame(summary(emm_object, infer = TRUE))
    
    # Dynamically find column names (emmeans naming varies)
    low_col  <- grep("lower|LCL", names(df), value = TRUE)
    upp_col  <- grep("upper|UCL", names(df), value = TRUE)
    stat_col <- grep("t.ratio|z.ratio", names(df), value = TRUE)
    
    # Add metadata columns
    df <- df %>%
      dplyr::mutate(
        Change = as.character(Change), 
        Sex = as.character(Sex), 
        Phase = as.character(Phase)
      )
    
    # Standardize column names
    if(length(low_col) > 0) df <- df %>% dplyr::rename(lwr = !!low_col[1])
    if(length(upp_col) > 0) df <- df %>% dplyr::rename(upr = !!upp_col[1])
    if(length(stat_col) > 0) df <- df %>% dplyr::rename(t.ratio = !!stat_col[1])
    
    return(df)
  }

  # 4. Assemble and return results
  # Mean estimates table
  avg_df <- safe_format(emm) %>%
    dplyr::rename(emmean_avg = emmean, SE_avg = SE)

  # Contrast table with p-values and Cohen's d effect sizes
  contr_df <- safe_format(contr_res) %>%
    dplyr::mutate(
      d = estimate / total_sd,
      SE_d = SE / total_sd,
      lwr_d = lwr / total_sd,
      upr_d = upr / total_sd
    )

  return(list(means = avg_df, contrasts = contr_df))
}

# Helper: enrich contrast table with BH-primary and Holm-sensitivity columns
add_dual_adjustment_columns <- function(df,
                                        family_cols = c("Change", "Sex", "Phase", "window"),
                                        primary_method = "BH",
                                        strict_method = "holm",
                                        primary_contrast_set = c("SUS-CON", "SUS-RES")) {
  if (nrow(df) == 0 || !"p.value" %in% names(df)) return(df)

  family_cols <- family_cols[family_cols %in% names(df)]
  if (length(family_cols) == 0) {
    df$contrast_family <- "all_rows"
  } else {
    df$contrast_family <- do.call(paste, c(df[family_cols], sep = "|"))
  }

  # Normalize contrast labels ("SUS - CON" -> "SUS-CON") for biological subset tagging
  if ("contrast" %in% names(df)) {
    df$contrast_norm <- gsub("\\s+", "", as.character(df$contrast))
    df$contrast_norm <- gsub("--", "-", df$contrast_norm)
    df$contrast_norm <- toupper(df$contrast_norm)
    
    # Check if each normalized contrast matches primary set (either direction)
    # e.g., "SUS-CON" matches "SUS-CON", and "CON-SUS" also matches (reversed)
    df$is_primary_contrast <- FALSE
    primary_upper <- toupper(primary_contrast_set)
    # Vectorized direct match
    df$is_primary_contrast <- df$contrast_norm %in% primary_upper
    # Vectorized reversed match (e.g. "CON-SUS" matches "SUS-CON")
    parts_split <- strsplit(df$contrast_norm, "-")
    two_part    <- lengths(parts_split) == 2L
    if (any(two_part)) {
      reversed <- vapply(parts_split[two_part], function(p) paste0(p[2], "-", p[1]), character(1))
      df$is_primary_contrast[two_part] <- df$is_primary_contrast[two_part] | (reversed %in% primary_upper)
    }
  } else {
    df$contrast_norm <- NA_character_
    df$is_primary_contrast <- FALSE
  }

  df %>%
    dplyr::mutate(
      p.value_raw = p.value
    ) %>%
    dplyr::group_by(contrast_family) %>%
    dplyr::mutate(
      n_tests_family = dplyr::n(),
      p.adjust_bh_family = p.adjust(p.value_raw, method = primary_method),
      p.adjust_holm_family = p.adjust(p.value_raw, method = strict_method)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(contrast_family) %>%
    dplyr::group_modify(~{
      d <- .x
      idx <- d$is_primary_contrast %in% TRUE
      d$p.adjust_bh_primary_family <- NA_real_
      d$p.adjust_holm_primary_family <- NA_real_
      d$n_tests_primary_family <- sum(idx, na.rm = TRUE)

      if (any(idx)) {
        d$p.adjust_bh_primary_family[idx] <- p.adjust(d$p.value_raw[idx], method = primary_method)
        d$p.adjust_holm_primary_family[idx] <- p.adjust(d$p.value_raw[idx], method = strict_method)
      }
      d
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      n_tests_global = dplyr::n(),
      p.adjust_bh_global = p.adjust(p.value_raw, method = primary_method),
      p.adjust_holm_global = p.adjust(p.value_raw, method = strict_method),
      q.value_family = p.adjust_bh_family,
      q.value_global = p.adjust_bh_global,
      q.value = q.value_family,
      p.adjust = q.value,
      p.adjust_strict = p.adjust_holm_family,
      p.adjust_primary = p.adjust_bh_primary_family,
      p.adjust_primary_strict = p.adjust_holm_primary_family,
      p.sign = dplyr::case_when(
        p.adjust < 0.001 ~ "***",
        p.adjust < 0.01  ~ "**",
        p.adjust < 0.05  ~ "*",
        p.adjust < 0.1   ~ "†",
        TRUE ~ "ns"
      ),
      p.sign_strict = dplyr::case_when(
        p.adjust_strict < 0.001 ~ "***",
        p.adjust_strict < 0.01  ~ "**",
        p.adjust_strict < 0.05  ~ "*",
        p.adjust_strict < 0.1   ~ "†",
        TRUE ~ "ns"
      ),
      p.sign_primary = dplyr::case_when(
        is.na(p.adjust_primary) ~ NA_character_,
        p.adjust_primary < 0.001 ~ "***",
        p.adjust_primary < 0.01  ~ "**",
        p.adjust_primary < 0.05  ~ "*",
        p.adjust_primary < 0.1   ~ "†",
        TRUE ~ "ns"
      ),
      reporting_notes = sprintf(
        "Primary=%s(n=%d SUS-focus); Strict=%s(n=%d family); Global=%s(n=%d all); †=p<0.1",
        toupper(primary_method),
        n_tests_primary_family,
        toupper(strict_method),
        n_tests_family,
        toupper(primary_method),
        n_tests_global
      )
    )
}

# Backward-compatible wrapper for phase-average tables
finalize_phase_contrast_pvals <- function(df, method = "BH", primary_contrast_set = c("SUS-CON", "SUS-RES")) {
  add_dual_adjustment_columns(
    df,
    family_cols = c("Change", "Sex", "Phase", "window"),
    primary_method = method,
    strict_method = "holm",
    primary_contrast_set = primary_contrast_set
  )
}

# Helper: generate summary statistics about primary vs exploratory contrast counts
summarize_contrast_stratification <- function(df) {
  if (nrow(df) == 0 || !"is_primary_contrast" %in% names(df)) {
    return(tibble::tibble(
      category = "Summary",
      n_primary = 0,
      n_exploratory = 0,
      n_total = 0,
      pct_primary = 0
    ))
  }
  
  primary_count <- sum(df$is_primary_contrast == TRUE, na.rm = TRUE)
  exploratory_count <- sum(df$is_primary_contrast == FALSE, na.rm = TRUE)
  total_count <- nrow(df)
  pct_primary <- if (total_count > 0) (primary_count / total_count) * 100 else 0
  
  tibble::tibble(
    category = "Overall Contrast Stratification",
    n_primary = primary_count,
    n_exploratory = exploratory_count,
    n_total = total_count,
    pct_primary = round(pct_primary, 1)
  )
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
# Data preparation and transformations
# -------------------------------------------------
# Applied transformations (documented to justify model assumptions):
#   Movement: log1p() transformation
#     REASON: Handles zero-inflation and right skewness typical of motor activity
#             log(x+1) avoids undefined log(0) and stabilizes variance for LMM
#   Proximity: No transformation applied
#     REASON: Approximately normally distributed; within-animal variation stable
#   ActivityIndex: No transformation applied
#     REASON: Calculated from normalized proxies; bounded scale [0,1] approx.
#
stopifnot(all(c("Change","Sex","Phase","Group","HalfHourElapsed","Movement","Proximity","AnimalNum", "ActivityIndex") %in% names(data_filtered_agg)))
data_filtered_agg <- data_filtered_agg %>%
  mutate(
    Change = factor(Change),
    Sex    = factor(Sex),
    Phase  = factor(Phase, levels = c("Active","Inactive")),
    Group  = factor(Group, levels = c("CON","RES","SUS")),
    HalfHourElapsed = as.integer(HalfHourElapsed),
    Movement = log1p(Movement)   # log(x + 1) variance-stabilizing transformation
  )

# -------------------------------------------------
# Publication-grade Excel writer
# -------------------------------------------------
write_pub_xlsx <- function(df, filepath, sheet_name = "Results",
                           p_cols        = c("p.value", "p.adjust", "p.adj_FDR"),
                           freeze_cols   = 0L,
                           note          = NULL) {
  # ── Guard ────────────────────────────────────────────────────────────────────
  if (nrow(df) == 0) return(invisible(NULL))

  # ── 1. Insert significance-star columns next to every detected p-column ─────
  sig_stars <- function(p) {
    dplyr::case_when(
      is.na(p)  ~ "",
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      p < 0.10  ~ "\u2020",   # † for trend
      TRUE      ~ "ns"
    )
  }

  df_out      <- df
  star_cols   <- character(0)
  active_pcols <- intersect(p_cols, names(df_out))

  for (pc in rev(active_pcols)) {   # rev so forward-insertions preserve order
    if (is.numeric(df_out[[pc]])) {
      sc <- paste0(pc, "_sig")
      idx <- which(names(df_out) == pc)
      df_out <- tibble::add_column(df_out, !!sc := sig_stars(df_out[[pc]]), .after = idx)
      star_cols <- c(star_cols, sc)
    }
  }

  # ── 2. Build workbook ────────────────────────────────────────────────────────
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheet_name, gridLines = TRUE)

  # ── 3. Style palette ─────────────────────────────────────────────────────────
  FN <- "Arial"

  # Header: modern blue, bold, centred, medium bottom border
  hdr <- openxlsx::createStyle(
    fontName = FN, fontSize = 10, textDecoration = "bold",
    fgFill   = "#1F4E79", fontColour = "#FFFFFF", halign = "CENTER", valign = "CENTER",
    border   = "Bottom", borderColour = "#173A5E", borderStyle = "medium",
    wrapText = FALSE
  )

  # Base text / number styles
  s_txt  <- openxlsx::createStyle(fontName = FN, fontSize = 9, halign = "LEFT",  valign = "CENTER")
  s_ctr  <- openxlsx::createStyle(fontName = FN, fontSize = 9, halign = "CENTER",valign = "CENTER")  # stars
  s_2dp  <- openxlsx::createStyle(fontName = FN, fontSize = 9, halign = "RIGHT", valign = "CENTER", numFmt = "0.00")
  s_3dp  <- openxlsx::createStyle(fontName = FN, fontSize = 9, halign = "RIGHT", valign = "CENTER", numFmt = "0.000")
  s_4dp  <- openxlsx::createStyle(fontName = FN, fontSize = 9, halign = "RIGHT", valign = "CENTER", numFmt = "0.0000")
  s_int  <- openxlsx::createStyle(fontName = FN, fontSize = 9, halign = "RIGHT", valign = "CENTER", numFmt = "0")
  s_sci  <- openxlsx::createStyle(fontName = FN, fontSize = 9, halign = "RIGHT", valign = "CENTER", numFmt = "0.00E+00")

  # Banded row fills (very subtle)
  s_odd  <- openxlsx::createStyle(fgFill = "#FFFFFF")
  s_even <- openxlsx::createStyle(fgFill = "#F7F7F7")
  # Whole-row highlighting tiers
  s_row_sig   <- openxlsx::createStyle(fgFill = "#EAF6EC")  # significant row (green tint)
  s_row_trend <- openxlsx::createStyle(fgFill = "#FFFDEB")  # trend row (light yellow)


  n_rows  <- nrow(df_out)
  n_cols  <- ncol(df_out)
  col_nms <- names(df_out)

  # ── 4. Write data ─────────────────────────────────────────────────────────────
  openxlsx::writeData(wb, sheet_name, df_out, headerStyle = hdr, rowNames = FALSE)

  # Row heights: header slightly taller
  openxlsx::setRowHeights(wb, sheet_name, rows = 1, heights = 20)

  # ── 5. Banded rows (vectorized: 2 calls instead of n_rows calls) ─────────────
  all_data_rows <- seq_len(n_rows) + 1L          # rows 2 .. n_rows+1
  odd_rows  <- all_data_rows[seq(1L, n_rows, by = 2L)]
  even_rows <- if (n_rows >= 2L) all_data_rows[seq(2L, n_rows, by = 2L)] else integer(0)
  if (length(odd_rows))  openxlsx::addStyle(wb, sheet_name, style = s_odd,  rows = odd_rows,  cols = 1:n_cols, gridExpand = TRUE, stack = FALSE)
  if (length(even_rows)) openxlsx::addStyle(wb, sheet_name, style = s_even, rows = even_rows, cols = 1:n_cols, gridExpand = TRUE, stack = FALSE)

  # Row-level trend/significance shading (vectorized: batch by tier instead of row-by-row)
  if (length(active_pcols) > 0) {
    p_mat <- suppressWarnings(
      matrix(apply(df[, active_pcols, drop = FALSE], 2L, as.numeric),
             nrow = n_rows)
    )
    min_p_per_row <- apply(p_mat, 1L, function(x) {
      x <- x[is.finite(x)]; if (!length(x)) NA_real_ else min(x)
    })
    sig_rows   <- which(!is.na(min_p_per_row) & min_p_per_row <  0.05) + 1L
    trend_rows <- which(!is.na(min_p_per_row) & min_p_per_row >= 0.05 & min_p_per_row < 0.10) + 1L
    if (length(sig_rows))   openxlsx::addStyle(wb, sheet_name, style = s_row_sig,   rows = sig_rows,   cols = 1:n_cols, gridExpand = TRUE, stack = TRUE)
    if (length(trend_rows)) openxlsx::addStyle(wb, sheet_name, style = s_row_trend, rows = trend_rows, cols = 1:n_cols, gridExpand = TRUE, stack = TRUE)
  }

  # ── 6. Smart per-column formatting ───────────────────────────────────────────
  p_pat   <- "^p\\.|p\\.value|p\\.adj|p_boot|p_BH|p_unadj|p_adj|pvalue"
  est_pat <- "estimate|Estimate|^beta$|^coef$|^b$|^Effect$"
  se_pat  <- "^SE$|std\\.error|std_error|^se\\."
  t_pat   <- "^t\\.ratio|^z\\.ratio|^F\\.ratio|statistic|^t$|^F$|^z$|^chi|^W$"
  d_pat   <- "^d$|cohen|^r$|eta2|omega2|R2|r2_|^rho|ICC|^AUC$|AUC_norm"
  ci_pat  <- "^ci_|CI_|^lwr$|^upr$|lower|upper|^conf|^asymp"
  int_pat <- "^n$|^N$|n_obs|n_models|count|n_"

  for (j in seq_len(n_cols)) {
    cn  <- col_nms[j]
    col <- df_out[[j]]

    if (cn %in% star_cols) {
      openxlsx::addStyle(wb, sheet_name, style = s_ctr,
                         rows = 2:(n_rows + 1), cols = j, gridExpand = FALSE, stack = TRUE)
      next
    }

    if (!is.numeric(col)) {
      openxlsx::addStyle(wb, sheet_name, style = s_txt,
                         rows = 2:(n_rows + 1), cols = j, gridExpand = FALSE, stack = TRUE)
      next
    }

    # Detect format by column name
    col_fmt <- dplyr::case_when(
      grepl(p_pat,   cn, ignore.case = TRUE)  ~ "p",
      grepl(int_pat, cn, ignore.case = TRUE)  ~ "int",
      grepl(t_pat,   cn, ignore.case = TRUE)  ~ "2dp",
      grepl(d_pat,   cn, ignore.case = TRUE)  ~ "3dp",
      grepl(ci_pat,  cn, ignore.case = TRUE)  ~ "3dp",
      grepl(est_pat, cn, ignore.case = TRUE)  ~ "3dp",
      grepl(se_pat,  cn, ignore.case = TRUE)  ~ "3dp",
      grepl("AIC|BIC|logLik|df_|^df$|^DF$", cn) ~ "2dp",
      grepl("edf|Ref_df|^k$",       cn)       ~ "2dp",
      grepl("auc|AUC",  cn, ignore.case = TRUE) ~ "3dp",
      TRUE                                     ~ "3dp"
    )

    # Apply style based on detected format (p cells get special cell-level handling below)
    s_apply <- switch(col_fmt,
      "p"   = s_4dp,
      "2dp" = s_2dp,
      "3dp" = s_3dp,
      "int" = s_int,
      s_3dp
    )
    openxlsx::addStyle(wb, sheet_name, style = s_apply,
                       rows = 2:(n_rows + 1), cols = j, gridExpand = FALSE, stack = TRUE)
  }

  # ── 7. Thin right border after the last label (non-numeric) column ───────────
  label_cols <- which(!sapply(df_out, is.numeric) & !(col_nms %in% star_cols))
  if (length(label_cols) > 0) {
    last_lbl <- max(label_cols)
    if (last_lbl < n_cols) {
      b_sty <- openxlsx::createStyle(border = "Right", borderColour = "#AAAAAA", borderStyle = "medium")
      openxlsx::addStyle(wb, sheet_name, style = b_sty,
                         rows = 1:(n_rows + 1), cols = last_lbl, gridExpand = FALSE, stack = TRUE)
    }
  }

  # ── 8. Significance legend footnote ──────────────────────────────────────────
  sig_note <- paste0(
    "*** p\u202f<\u202f0.001   ** p\u202f<\u202f0.01   * p\u202f<\u202f0.05",
    "   \u2020 p\u202f<\u202f0.10 (trend)   ns = not significant",
    "   | row tint: green = any p<0.05, yellow = any p<0.10",
    if (!is.null(note)) paste0("     \u2502 ", note) else ""
  )
  note_row <- n_rows + 3L
  openxlsx::writeData(wb, sheet_name, data.frame(x = sig_note),
                      startRow = note_row, startCol = 1, colNames = FALSE)
  openxlsx::addStyle(wb, sheet_name,
    style = openxlsx::createStyle(fontName = FN, fontSize = 8,
                                   fontColour = "#666666", textDecoration = "italic"),
    rows = note_row, cols = 1, stack = FALSE)
  if (n_cols >= 2)
    openxlsx::mergeCells(wb, sheet_name, cols = 1:min(n_cols, 10), rows = note_row)

  # ── 10. Freeze pane + auto column widths ─────────────────────────────────────
  openxlsx::freezePane(wb, sheet_name, firstRow = TRUE,
                       firstCol = as.integer(freeze_cols) > 0L)
  openxlsx::setColWidths(wb, sheet_name, cols = 1:n_cols, widths = "auto")

  openxlsx::saveWorkbook(wb, filepath, overwrite = TRUE)
  invisible(wb)
}

# Helper: extract model fit indices
get_model_fit <- function(model) {
  r2 <- tryCatch(MuMIn::r.squaredGLMM(model)[1, ], error = function(e) c(R2m = NA_real_, R2c = NA_real_))
  data.frame(
    AIC           = tryCatch(AIC(model),  error = function(e) NA_real_),
    BIC           = tryCatch(BIC(model),  error = function(e) NA_real_),
    logLik        = tryCatch(as.numeric(logLik(model)), error = function(e) NA_real_),
    R2_marginal   = unname(r2["R2m"]),
    R2_conditional = unname(r2["R2c"])
  )
}

# Helper: get Wald 95% CIs for fixed effects
get_fixed_ci <- function(model) {
  tryCatch(
    as.data.frame(confint(model, method = "Wald", parm = "beta_")),
    error = function(e) NULL
  )
}

# Helper: publication-grade fixed-effects table with robust alignment
build_fixed_effects_table <- function(model, ms, fit_idx, Change, Phase, Sex, converged, window_label) {
  coef_df <- as.data.frame(ms$coefficients)
  coef_df$Fixed_effect <- rownames(coef_df)

  # Decode Group sum-coding so terms like Group1/Group2 are interpretable
  mf <- tryCatch(model.frame(model), error = function(e) NULL)
  group_levels <- tryCatch(levels(mf$Group), error = function(e) character(0))
  group_contr <- tryCatch(stats::contrasts(mf$Group), error = function(e) NULL)

  # Human-readable algebra for 3-level contr.sum coding (common case: CON/RES/SUS)
  group_levels_txt <- if (length(group_levels)) paste(group_levels, collapse = ", ") else NA_character_
  coding_scheme <- if (!is.null(group_contr)) "contr.sum" else NA_character_
  implied_group_means <- NA_character_
  implied_pairwise <- NA_character_
  if (!is.null(group_contr) && length(group_levels) == 3 && ncol(group_contr) == 2) {
    g1 <- group_levels[1]; g2 <- group_levels[2]; g3 <- group_levels[3]
    implied_group_means <- paste0(
      g1, " = Intercept + Group1; ",
      g2, " = Intercept + Group2; ",
      g3, " = Intercept - Group1 - Group2"
    )
    implied_pairwise <- paste0(
      g2, "-", g1, " = Group2 - Group1; ",
      g3, "-", g1, " = -2*Group1 - Group2; ",
      g3, "-", g2, " = -Group1 - 2*Group2"
    )
  }

  decode_group_term <- function(term) {
    out <- list(term_label = term, coding_weights = NA_character_, interpretation = NA_character_)

    if (term == "(Intercept)") {
      out$term_label <- "Intercept (grand mean; sum-coded)"
      out$interpretation <- "Grand mean across Group levels (with contr.sum coding)"
      return(out)
    }

    if (term == "TH_scaled") {
      out$term_label <- "TH_scaled (time effect at grand mean)"
      out$interpretation <- "Time slope at the sum-coded grand mean across groups"
      return(out)
    }

    # Main group terms: Group1, Group2, ...
    if (!is.null(group_contr) && grepl("^Group\\d+$", term)) {
      idx <- suppressWarnings(as.integer(sub("^Group", "", term)))
      if (!is.na(idx) && idx >= 1 && idx <= ncol(group_contr)) {
        w <- group_contr[, idx]
        w_txt <- paste(sprintf("%s:%+g", rownames(group_contr), as.numeric(w)), collapse = "; ")
        c_name <- if (!is.null(colnames(group_contr)) && nzchar(colnames(group_contr)[idx])) colnames(group_contr)[idx] else paste0("C", idx)
        out$term_label <- sprintf("Group contrast %s (sum-coded across %s)", c_name, group_levels_txt)
        out$coding_weights <- w_txt
        out$interpretation <- "Linear contrast under contr.sum; not a simple pairwise difference"
        return(out)
      }
    }

    # Interaction terms: Group1:TH_scaled, Group2:TH_scaled, ...
    if (!is.null(group_contr) && grepl("^Group\\d+:TH_scaled$", term)) {
      idx <- suppressWarnings(as.integer(sub("^Group", "", sub(":TH_scaled$", "", term))))
      if (!is.na(idx) && idx >= 1 && idx <= ncol(group_contr)) {
        w <- group_contr[, idx]
        w_txt <- paste(sprintf("%s:%+g", rownames(group_contr), as.numeric(w)), collapse = "; ")
        c_name <- if (!is.null(colnames(group_contr)) && nzchar(colnames(group_contr)[idx])) colnames(group_contr)[idx] else paste0("C", idx)
        out$term_label <- sprintf("Group contrast %s × TH_scaled (levels: %s)", c_name, group_levels_txt)
        out$coding_weights <- w_txt
        out$interpretation <- "Difference-in-slope contrast under contr.sum coding"
        return(out)
      }
    }

    out
  }

  decoded <- lapply(coef_df$Fixed_effect, decode_group_term)
  coef_df$term_label <- vapply(decoded, function(x) x$term_label, character(1))
  coef_df$coding_weights <- vapply(decoded, function(x) x$coding_weights, character(1))
  coef_df$interpretation <- vapply(decoded, function(x) x$interpretation, character(1))

  # Robust p-value extraction (column name differs across model classes)
  p_col <- grep("Pr\\(>|t\\)|Pr\\(>|z\\)", colnames(coef_df), value = TRUE)
  if (length(p_col) == 0) {
    coef_df$p.value <- NA_real_
  } else {
    coef_df$p.value <- as.numeric(coef_df[[p_col[1]]])
  }

  # Robust df extraction
  coef_df$df_Satterthwaite <- if ("df" %in% colnames(coef_df)) as.numeric(coef_df$df) else NA_real_

  # Robust t/z statistic extraction
  stat_col <- if ("t value" %in% colnames(coef_df)) "t value" else if ("z value" %in% colnames(coef_df)) "z value" else NA_character_
  coef_df$stat_value <- if (!is.na(stat_col)) as.numeric(coef_df[[stat_col]]) else NA_real_

  # Align CIs by term name (never by row position)
  ci_fx <- get_fixed_ci(model)
  if (!is.null(ci_fx)) {
    ci_fx$Fixed_effect <- rownames(ci_fx)
    colnames(ci_fx)[1:2] <- c("CI_lower_95", "CI_upper_95")
    coef_df <- dplyr::left_join(coef_df, ci_fx[, c("Fixed_effect", "CI_lower_95", "CI_upper_95")], by = "Fixed_effect")
  } else {
    coef_df$CI_lower_95 <- NA_real_
    coef_df$CI_upper_95 <- NA_real_
  }

  coef_df <- coef_df %>%
    dplyr::mutate(
      Fixed_effect_display = dplyr::if_else(!is.na(term_label) & nzchar(term_label), term_label, Fixed_effect),
      Estimate = as.numeric(.data$Estimate),
      Std.Error = as.numeric(.data$`Std. Error`),
      p.adj_FDR = p.adjust(.data$p.value, method = "fdr"),
      p.sign = dplyr::case_when(
        p.adj_FDR < 0.001 ~ "***",
        p.adj_FDR < 0.01  ~ "**",
        p.adj_FDR < 0.05  ~ "*",
        p.adj_FDR < 0.1   ~ "†",
        TRUE              ~ "ns"
      ),
      p.sign_raw = dplyr::case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        p.value < 0.1   ~ "†",
        TRUE            ~ "ns"
      ),
      term_class = dplyr::case_when(
        Fixed_effect == "(Intercept)" ~ "intercept",
        grepl(":", Fixed_effect) ~ "interaction",
        grepl("TH_scaled", Fixed_effect) ~ "time",
        grepl("^Group", Fixed_effect) ~ "group",
        TRUE ~ "other"
      ),
      estimate_ci95 = ifelse(
        is.na(CI_lower_95) | is.na(CI_upper_95),
        sprintf("%.3f", Estimate),
        sprintf("%.3f [%.3f, %.3f]", Estimate, CI_lower_95, CI_upper_95)
      ),
      estimate_direction = dplyr::case_when(
        Estimate > 0 ~ "positive",
        Estimate < 0 ~ "negative",
        TRUE ~ "neutral"
      ),
      coding_scheme   = coding_scheme,
      group_levels    = group_levels_txt,
      implied_means   = implied_group_means,
      implied_pairwise_contrasts = implied_pairwise,
      AIC            = fit_idx$AIC,
      BIC            = fit_idx$BIC,
      logLik         = fit_idx$logLik,
      R2_marginal    = fit_idx$R2_marginal,
      R2_conditional = fit_idx$R2_conditional,
      model_formula  = paste(deparse(formula(model)), collapse = " "),
      singular       = lme4::isSingular(model),
      Change         = as.character(Change),
      Phase          = as.character(Phase),
      Sex            = as.character(Sex),
      converged      = converged,
      window         = window_label
    )

  tibble::as_tibble(coef_df %>%
    dplyr::select(
      Fixed_effect_display, Fixed_effect, term_label, term_class, coding_weights, interpretation,
      Estimate, CI_lower_95, CI_upper_95, estimate_ci95, estimate_direction,
      Std.Error, df_Satterthwaite,
      t.value = stat_value,
      p.value, p.adj_FDR, p.sign, p.sign_raw,
      coding_scheme, group_levels, implied_means, implied_pairwise_contrasts,
      AIC, BIC, logLik, R2_marginal, R2_conditional,
      model_formula, singular,
      Change, Phase, Sex, converged, window
    ))
}

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
    logs          = file.path(results_dir, "logs"),
    models        = file.path(results_dir, "models"),
    artifacts     = file.path(results_dir, "artifacts"),
    tables_sum    = file.path(results_dir, "tables", "summary"),
    tables_emm    = file.path(results_dir, "tables", "emmeans"),
    tables_contr  = file.path(results_dir, "tables", "contrasts"),
    tables_inv    = file.path(results_dir, "tables", "inventory"),
    tables_man    = file.path(results_dir, "tables", "manifest"),
    tables_phase  = file.path(results_dir, "tables", "phase_average"),
    tables_day1   = file.path(results_dir, "tables", "day1_active"),
    plots_line    = file.path(results_dir, "plots", "line"),
    plots_heat    = file.path(results_dir, "plots", "heatmap"),
    plots_fixed   = file.path(results_dir, "plots", "fixed_effects"),
    plots_pub     = file.path(results_dir, "plots", "publication"),
    edgeAggAll    = file.path(results_dir, "plots", "edge_emm"),
    plots_day1    = file.path(results_dir, "plots", "day1_active")
  )
  invisible(lapply(dirs, dir_create_safe))
  
  # Directory helper functions
  plot_subdir <- function(base_dir, Change, Sex, includeSex) {
    sub <- file.path(base_dir, paste0("CC-", Change))
    dir_create_safe(sub)
    if (isTRUE(includeSex)) {
      sub <- file.path(sub, paste0("Sex-", as.character(Sex)))
      dir_create_safe(sub)
    }
    sub
  }
  plot_dir_line  <- function(Change, Sex, includeSex) plot_subdir(dirs$plots_line,  Change, Sex, includeSex)
  plot_dir_heat  <- function(Change, Sex, includeSex) plot_subdir(dirs$plots_heat,  Change, Sex, includeSex)
  plot_dir_fixed <- function(Change, Sex)              plot_subdir(dirs$plots_fixed, Change, Sex, includeSex = TRUE)
  
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
  pal <- unname(group_colors)
  p <- ggplot(df_sex, aes(HalfHourElapsed, .data[[metric_name]], color = Group, group = Group)) +
    geom_path(stat = "summary", fun = mean, linewidth = 0.8, na.rm = TRUE) +
  stat_summary(aes(fill = Group), fun = mean,
       fun.min = function(x) mean(x) - sd(x),
       fun.max = function(x) mean(x) + sd(x),
      geom = "ribbon", alpha = 0.3, colour = NA, na.rm = TRUE) +
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
  # Modeling (main slices + optional first/last 2h per phase)
  # -------------------------------------------------
  changes <- if (includeChange) unique(data$Change) else "allChanges"
  sexes   <- if (includeSex)    unique(data$Sex)    else "allSexes"
  phases  <- if (includePhase)  unique(data$Phase)  else "allPhases"

  log_file <- file.path(dirs$logs, paste0("fit_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

  fit_model_slices <- function(window_label = "all") {
  summary_results_i <- tibble::tibble()
  emmeans_results_i <- tibble::tibble()
  phase_avg_means_i <- list()
  phase_avg_contr_i <- list()

  # Flat combination grid — enables both sequential (progress bar) and parallel dispatch
  combo_grid <- expand.grid(Change = changes, Sex = sexes, Phase = phases,
                            stringsAsFactors = FALSE)
  total_iterations <- nrow(combo_grid)

  # Register parallel backend if n_cores > 1
  if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    cat(sprintf("  [parallel] Using %d cores for %s (%s)\n", n_cores, metric_name, window_label))
    `%dispatch%` <- foreach::`%dopar%`
  } else {
    pb <- tkProgressBar(
      title = paste("Model Fitting Progress:", metric_name, "-", window_label),
      min = 0, max = total_iterations, width = 300
    )
    `%dispatch%` <- foreach::`%do%`
  }

  iter_results <- foreach::foreach(
    row_i = seq_len(total_iterations),
    .packages = c("lme4","lmerTest","emmeans","dplyr","tibble","MuMIn"),
    .export   = c("fit_lmer_adaptive","quiet_singular","compute_phase_avg_for_model",
                  "build_fixed_effects_table","get_model_fit","get_fixed_ci","contrast_config",
                  "use_group_time_interaction","min_animals","min_groups","min_obs","min_time_lv"),
    .errorhandling = "pass"
  ) %dispatch% {

    if (n_cores <= 1L) {
      setTkProgressBar(pb, row_i, label = paste("Progress:", round(row_i / total_iterations * 100, 2), "%"))
    }

    Change <- combo_grid$Change[row_i]
    Sex    <- combo_grid$Sex[row_i]
    Phase  <- combo_grid$Phase[row_i]

    log_msgs <- character(0)

      # Base filter
      dsub <- data %>%
      dplyr::filter(
      (!includeChange | Change == !!Change),
      (!includeSex    | Sex    == !!Sex),
      (!includePhase  | Phase  == !!Phase)
      )
      if (nrow(dsub) == 0) return(NULL)

      # Optional window per phase
      if (window_label %in% c("first2h", "last2h")) {
      t0_phase <- min(dsub$HalfHourElapsed, na.rm = TRUE)
      t1_phase <- max(dsub$HalfHourElapsed, na.rm = TRUE)

      if (window_label == "first2h") {
        dsub <- dsub %>%
        dplyr::mutate(rel_phase_h = (HalfHourElapsed - t0_phase) / 2) %>%
        dplyr::filter(!is.na(rel_phase_h), rel_phase_h >= 0, rel_phase_h <= 2)
      }

      if (window_label == "last2h") {
        dsub <- dsub %>%
        dplyr::mutate(rel_phase_h_from_end = (t1_phase - HalfHourElapsed) / 2) %>%
        dplyr::filter(!is.na(rel_phase_h_from_end), rel_phase_h_from_end >= 0, rel_phase_h_from_end <= 2)
      }

      if (nrow(dsub) == 0) {
        log_msgs <- c(log_msgs, paste(Sys.time(), Change, Sex, Phase, "SKIP", window_label, ": no rows in requested 2h window"))
        return(list(log_msgs = log_msgs))
      }
      }

      # Check minimum requirements
      n_anim  <- dplyr::n_distinct(dsub$AnimalNum)
      n_group <- dplyr::n_distinct(dsub$Group)
      if (n_anim < min_animals || n_group < min_groups || nrow(dsub) < min_obs) {
      log_msgs <- c(log_msgs, paste(Sys.time(), Change, Sex, Phase, "SKIP", window_label, ": few_animals/groups/obs"))
      return(list(log_msgs = log_msgs))
      }

      # Time variable
      t0 <- min(dsub$HalfHourElapsed, na.rm = TRUE)
      dsub <- dsub %>%
      dplyr::mutate(
        TH = (HalfHourElapsed - t0) / 2,
        TH_scaled = as.vector(scale(TH))
      ) %>%
      dplyr::filter(!is.na(TH_scaled), !is.na(.data[[metric_name]]))

      if (nrow(dsub) == 0) {
      log_msgs <- c(log_msgs, paste(Sys.time(), Change, Sex, Phase, "SKIP", window_label, ": no valid data after scaling"))
      return(list(log_msgs = log_msgs))
      }

      # Model structure
      use_int <- use_group_time_interaction && (dplyr::n_distinct(dsub$TH_scaled) >= min_time_lv)
      fixed_rhs <- if (use_int) "Group * TH_scaled" else "Group + TH_scaled"

      per_animal <- dsub %>% dplyr::count(AnimalNum, name = "n_i")
      allow_slope <- all(per_animal$n_i >= 3) && (dplyr::n_distinct(dsub$TH_scaled) >= min_time_lv)

      model <- quiet_singular(
      fit_lmer_adaptive(
      fixed_rhs, dsub,
      response_var = metric_name,
      time_var = "TH_scaled",
      allow_slope = allow_slope)
      )

      # Convergence
      converged <- TRUE
      if (inherits(model, "lmerMod")) {
      conv_code <- model@optinfo$conv$opt
      if (!is.null(conv_code) && conv_code != 0) {
        converged <- FALSE
        log_msgs <- c(log_msgs, paste(Sys.time(), Change, Sex, Phase, window_label,
          "WARNING: Model convergence issue - code:", conv_code))
      }
      }

      ms <- summary(model)
      log_msgs <- c(log_msgs, paste(Sys.time(), Change, Sex, Phase, window_label,
        "MAIN form:", deparse(formula(model)),
        "singular:", lme4::isSingular(model),
        "converged:", converged))

      # Save model
      model_filename <- paste0(
      "model_", metric_name, "_",
      "Change-", Change, "_",
      "Sex-", Sex, "_",
      "Phase-", Phase, "_",
      "Window-", window_label,
      ".rds"
      )
      saveRDS(model, file.path(dirs$models, model_filename))

      # Fixed effects — publication-grade (robust helper)
      fit_idx  <- get_model_fit(model)
      fx <- build_fixed_effects_table(
        model = model, ms = ms, fit_idx = fit_idx,
        Change = Change, Phase = Phase, Sex = Sex,
        converged = converged, window_label = window_label
      )

      iter_emm    <- tibble::tibble()
      iter_means  <- list()
      iter_contr  <- list()

      if (!converged) {
      log_msgs <- c(log_msgs, paste(Sys.time(), Change, Sex, Phase, window_label, "SKIP emmeans: convergence issue"))
      return(list(fx = fx, log_msgs = log_msgs))
      }

      # emmeans only if Group in model
      if (!any(grepl("^Group", names(fixef(model))))) {
        return(list(fx = fx, log_msgs = log_msgs))
      }

      emm_lvl <- suppressMessages(
      emmeans::emmeans(model, pairwise ~ Group, at = list(TH_scaled = 0))
      )
      emm_lvl_df <- as.data.frame(summary(emm_lvl$contrasts, infer = TRUE, adjust = "none")) %>%
      dplyr::mutate(
        contrast_type = "level@intercept",
        Change        = as.character(Change),
        Phase         = as.character(Phase),
        Sex           = as.character(Sex),
        window        = window_label
      )
      emm_lvl_df <- emm_lvl_df[!is.na(emm_lvl_df$p.value), ]
      emm_lvl_df$cohens_d <- tryCatch({
        emm_lvl_raw <- emmeans::emmeans(model, ~ Group, at = list(TH_scaled = 0))
        eff_d <- emmeans::eff_size(emm_lvl_raw, sigma = sigma(model), edf = df.residual(model))
        eff_df <- as.data.frame(eff_d)
        # match by contrast string
        idx <- match(emm_lvl_df$contrast, eff_df$contrast)
        eff_df$effect.size[idx]
      }, error = function(e) NA_real_)
      if (nrow(emm_lvl_df) > 0) iter_emm <- dplyr::bind_rows(iter_emm, emm_lvl_df)

      if (use_int && any(grepl("Group:TH_scaled", names(fixef(model))))) {
      emm_tr <- suppressMessages(
        emmeans::emtrends(model, pairwise ~ Group, var = "TH_scaled")
      )
      emm_tr_df <- as.data.frame(summary(emm_tr$contrasts, infer = TRUE, adjust = "none")) %>%
        dplyr::mutate(
          contrast_type = "slope",
          Change        = as.character(Change),
          Phase         = as.character(Phase),
          Sex           = as.character(Sex),
          window        = window_label
        )
      emm_tr_df <- emm_tr_df[!is.na(emm_tr_df$p.value), ]
      emm_tr_df$cohens_d <- NA_real_
      if (nrow(emm_tr_df) > 0) iter_emm <- dplyr::bind_rows(iter_emm, emm_tr_df)
      }

      # Phase-average stats
      if (inherits(model, "lmerMod") && converged) {
      stat_k <- compute_phase_avg_for_model(model, Change, Sex, Phase)
      if (!is.null(stat_k$means) && nrow(stat_k$means)) {
        stat_k$means$window <- window_label
        iter_means[[1]] <- stat_k$means
      }
      if (!is.null(stat_k$contrasts) && nrow(stat_k$contrasts)) {
        stat_k$contrasts$window <- window_label
        iter_contr[[1]] <- stat_k$contrasts
      }
      }

    list(fx = fx, emm = iter_emm, means = iter_means, contr = iter_contr, log_msgs = log_msgs)
  } # end foreach

  if (n_cores <= 1L) close(pb)

  # Flush collected log messages to file (safe: sequential write after all iterations)
  for (res in iter_results) {
    if (!is.null(res) && !is.null(res$log_msgs) && length(res$log_msgs) > 0) {
      cat(paste(res$log_msgs, collapse = "\n"), "\n", file = log_file, append = TRUE)
    }
  }

  # Combine results from all iterations
  for (res in iter_results) {
    if (is.null(res)) next
    if (!is.null(res$fx)    && nrow(res$fx)    > 0) summary_results_i <- dplyr::bind_rows(summary_results_i, res$fx)
    if (!is.null(res$emm)   && nrow(res$emm)   > 0) emmeans_results_i  <- dplyr::bind_rows(emmeans_results_i,  res$emm)
    if (length(res$means) > 0) phase_avg_means_i <- c(phase_avg_means_i, res$means)
    if (length(res$contr) > 0) phase_avg_contr_i <- c(phase_avg_contr_i, res$contr)
  }

  if (nrow(emmeans_results_i) > 0) {
    emmeans_results_i <- add_dual_adjustment_columns(
      emmeans_results_i,
      family_cols = c("Change", "Sex", "Phase", "window", "contrast_type"),
      primary_method = "BH",
      strict_method = "holm",
      primary_contrast_set = primary_contrast_set
    )
  }

  list(
    summary = summary_results_i,
    emmeans = emmeans_results_i,
    phase_means = if (length(phase_avg_means_i)) dplyr::bind_rows(phase_avg_means_i) else tibble::tibble(),
    phase_contr = if (length(phase_avg_contr_i)) dplyr::bind_rows(phase_avg_contr_i) else tibble::tibble()
  )
  }

  # -------------------------------------------------
  # Day 1 Active Phase Analysis (Change == "1", Phase == "Active")
  # -------------------------------------------------
  fit_model_day1_active <- function() {
  cat("\n--- Fitting model for Day 1 / First Active Phase ---\n")

  # Identify the first cage change level
  change_levels <- sort(unique(as.character(data$Change)))
  first_change <- change_levels[1]

  # Find the first Active phase block within the first Change
  # (there may be multiple Active phases across days)
  first_change_active <- data %>%
    dplyr::filter(
      as.character(Change) == first_change,
      as.character(Phase) == "Active",
      !is.na(.data[[metric_name]]),
      !is.na(AnimalNum),
      !is.na(Group)
    ) %>%
    dplyr::arrange(HalfHourElapsed)

  # Detect contiguous Active phase blocks by looking for gaps in HalfHourElapsed
  # A new block starts when the difference to the previous row > 1 (i.e., not consecutive)
  if (nrow(first_change_active) > 0) {
    first_change_active <- first_change_active %>%
      dplyr::arrange(HalfHourElapsed) %>%
      dplyr::mutate(
        time_diff = c(0, diff(HalfHourElapsed)),
        block_id  = cumsum(time_diff > 1)
      )
    first_block_id <- min(first_change_active$block_id)
    dsub_day1 <- first_change_active %>%
      dplyr::filter(block_id == first_block_id) %>%
      dplyr::select(-time_diff, -block_id)
  } else {
    dsub_day1 <- first_change_active
  }

  if (nrow(dsub_day1) == 0) {
    cat("SKIP day1_active: no data found for first Change + Active phase.\n")
    return(NULL)
  }

  sexes_day1 <- if (includeSex) unique(as.character(dsub_day1$Sex)) else "allSexes"

  day1_summary   <- tibble::tibble()
  day1_emmeans   <- tibble::tibble()
  day1_means_lst <- list()
  day1_contr_lst <- list()

  pb_day1 <- tkProgressBar(
    title = paste("Day1 Active Phase:", metric_name),
    min = 0, max = length(sexes_day1), width = 300
  )

  for (idx_s in seq_along(sexes_day1)) {
    Sex <- sexes_day1[idx_s]
    setTkProgressBar(pb_day1, idx_s, label = paste("Sex:", Sex))

    ds <- if (includeSex) {
    dsub_day1 %>% dplyr::filter(Sex == !!Sex)
    } else {
    dsub_day1
    }

    if (nrow(ds) == 0 ||
      dplyr::n_distinct(ds$AnimalNum) < min_animals ||
      dplyr::n_distinct(ds$Group) < min_groups ||
      nrow(ds) < min_obs) {
    cat(paste("SKIP day1_active Sex:", Sex, "- insufficient data\n"),
      file = log_file, append = TRUE)
    next
    }

    t0 <- min(ds$HalfHourElapsed, na.rm = TRUE)
    ds <- ds %>%
    dplyr::mutate(
      TH = (HalfHourElapsed - t0) / 2,
      TH_scaled = as.vector(scale(TH))
    ) %>%
    dplyr::filter(!is.na(TH_scaled), !is.na(.data[[metric_name]]))

    if (nrow(ds) == 0) next

    use_int <- use_group_time_interaction && (dplyr::n_distinct(ds$TH_scaled) >= min_time_lv)
    fixed_rhs <- if (use_int) "Group * TH_scaled" else "Group + TH_scaled"

    per_animal <- ds %>% dplyr::count(AnimalNum, name = "n_i")
    allow_slope <- all(per_animal$n_i >= 3) && (dplyr::n_distinct(ds$TH_scaled) >= min_time_lv)

    model_day1 <- quiet_singular(
    fit_lmer_adaptive(
      fixed_rhs, ds,
      response_var = metric_name,
      time_var = "TH_scaled",
      allow_slope = allow_slope
    )
    )

    converged <- TRUE
    if (inherits(model_day1, "lmerMod")) {
    conv_code <- model_day1@optinfo$conv$opt
    if (!is.null(conv_code) && conv_code != 0) {
      converged <- FALSE
      cat(paste(Sys.time(), "day1_active Sex:", Sex, "WARNING convergence code:", conv_code, "\n"),
        file = log_file, append = TRUE)
    }
    }

    ms <- summary(model_day1)
    cat(paste(Sys.time(), "day1_active Sex:", Sex,
        "form:", deparse(formula(model_day1)),
        "singular:", lme4::isSingular(model_day1),
        "converged:", converged, "\n"),
      file = log_file, append = TRUE)

    saveRDS(model_day1,
        file.path(dirs$models,
            paste0("model_", metric_name, "_day1_active_Sex-", Sex, ".rds")))

    fit_idx_day1 <- get_model_fit(model_day1)
    fx <- build_fixed_effects_table(
      model = model_day1, ms = ms, fit_idx = fit_idx_day1,
      Change = first_change, Phase = "Active", Sex = as.character(Sex),
      converged = converged, window_label = "day1_active"
    )
    day1_summary <- dplyr::bind_rows(day1_summary, fx)

    if (!converged) next
    if (!any(grepl("^Group", names(fixef(model_day1))))) next

    emm_lvl <- suppressMessages(
    emmeans::emmeans(model_day1, pairwise ~ Group, at = list(TH_scaled = 0))
    )
    emm_lvl_df <- as.data.frame(summary(emm_lvl$contrasts, infer = TRUE, adjust = "none")) %>%
    dplyr::mutate(
      contrast_type = "level@intercept",
      Change        = first_change,
      Phase         = "Active",
      Sex           = as.character(Sex),
      window        = "day1_active"
    )
    emm_lvl_df <- emm_lvl_df[!is.na(emm_lvl_df$p.value), ]
    emm_lvl_df$cohens_d <- tryCatch({
      emm_lvl_raw_d1 <- emmeans::emmeans(model_day1, ~ Group, at = list(TH_scaled = 0))
      eff_d1 <- emmeans::eff_size(emm_lvl_raw_d1, sigma = sigma(model_day1), edf = df.residual(model_day1))
      eff_df_d1 <- as.data.frame(eff_d1)
      idx <- match(emm_lvl_df$contrast, eff_df_d1$contrast)
      eff_df_d1$effect.size[idx]
    }, error = function(e) NA_real_)
    if (nrow(emm_lvl_df) > 0) day1_emmeans <- dplyr::bind_rows(day1_emmeans, emm_lvl_df)

    if (use_int && any(grepl("Group:TH_scaled", names(fixef(model_day1))))) {
    emm_tr <- suppressMessages(
      emmeans::emtrends(model_day1, pairwise ~ Group, var = "TH_scaled")
    )
    emm_tr_df <- as.data.frame(summary(emm_tr$contrasts, infer = TRUE, adjust = "none")) %>%
      dplyr::mutate(
        contrast_type = "slope",
        Change        = first_change,
        Phase         = "Active",
        Sex           = as.character(Sex),
        window        = "day1_active"
      )
    emm_tr_df <- emm_tr_df[!is.na(emm_tr_df$p.value), ]
    emm_tr_df$cohens_d <- NA_real_
    if (nrow(emm_tr_df) > 0) day1_emmeans <- dplyr::bind_rows(day1_emmeans, emm_tr_df)
    }

    stat_k <- compute_phase_avg_for_model(model_day1, first_change, Sex, "Active")
    if (!is.null(stat_k$means) && nrow(stat_k$means)) {
    stat_k$means$window <- "day1_active"
    day1_means_lst[[length(day1_means_lst) + 1]] <- stat_k$means
    }
    if (!is.null(stat_k$contrasts) && nrow(stat_k$contrasts)) {
    stat_k$contrasts$window <- "day1_active"
    day1_contr_lst[[length(day1_contr_lst) + 1]] <- stat_k$contrasts
    }
  }

  close(pb_day1)

  if (nrow(day1_emmeans) > 0) {
    day1_emmeans <- add_dual_adjustment_columns(
      day1_emmeans,
      family_cols = c("Change", "Sex", "Phase", "window", "contrast_type"),
      primary_method = "BH",
      strict_method = "holm",
      primary_contrast_set = primary_contrast_set
    )
  }

  list(
    summary     = day1_summary,
    emmeans     = day1_emmeans,
    phase_means = if (length(day1_means_lst)) dplyr::bind_rows(day1_means_lst) else tibble::tibble(),
    phase_contr = if (length(day1_contr_lst)) dplyr::bind_rows(day1_contr_lst) else tibble::tibble(),
    first_change = first_change,
    dsub         = dsub_day1
  )
  }

  # -------------------------------------------------
  # Day 1 Active Phase Plots
  # -------------------------------------------------
  plot_day1_active <- function(day1_res) {
  if (is.null(day1_res)) return(invisible(NULL))

  dsub   <- day1_res$dsub
  first_change <- day1_res$first_change
  run_scope_d1 <- paste0("Change", first_change, "_Active")

  # --- Raw time-course line plot ---
  t0_d1 <- min(dsub$HalfHourElapsed, na.rm = TRUE)
  dsub_plot <- dsub %>%
    dplyr::mutate(
    rel_h = (HalfHourElapsed - t0_d1) / 2,
    Group = factor(Group, levels = c("CON","RES","SUS"))
    )

  p_line_d1 <- ggplot(dsub_plot,
            aes(x = rel_h, y = .data[[metric_name]],
              color = Group, fill = Group, group = Group)) +
      geom_path(stat = "summary", fun = mean, linewidth = 0.9, na.rm = TRUE) +
    stat_summary(fun = mean,
           fun.min = function(x) mean(x) - sd(x),
           fun.max = function(x) mean(x) + sd(x),
        geom = "ribbon", alpha = 0.25, colour = NA, na.rm = TRUE) +
    scale_color_manual(values = group_colors) +
    scale_fill_manual(values  = group_colors) +
    labs(
    title    = paste0(metric_name, ": First Active Phase (CC", first_change, ")"),
    subtitle = "Mean ± SD",
    x        = "Time since phase start [h]",
    y        = paste(metric_name, "[a.u.]")
    ) +
    theme_classic(base_size = 11) +
    theme(
    legend.position = "top",
    plot.title      = element_text(face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(hjust = 0.5, size = 9)
    )

  if (includeSex && "Sex" %in% names(dsub_plot) && dplyr::n_distinct(dsub_plot$Sex) > 1) {
    p_line_d1 <- p_line_d1 + facet_wrap(~ Sex)
  }

  ggsave(file.path(dirs$plots_day1,
           paste0("timecourse_day1_", run_scope_d1, ".svg")),
       p_line_d1, width = 7, height = 4)

  # --- Phase-average forest plot (contrasts) ---
  if (!is.null(day1_res$phase_contr) && nrow(day1_res$phase_contr) > 0) {

    canonicalize_contrast_d1 <- function(df) {
    if (!nrow(df)) return(df)
    parts <- strsplit(df$contrast, "\\s*-\\s*")
    lhs <- vapply(parts, `[`, character(1), 1)
    rhs <- vapply(parts, `[`, character(1), 2)
    lhs <- trimws(lhs); rhs <- trimws(rhs)

    flip <- lhs == "CON" & rhs %in% c("RES","SUS")
    if (any(flip)) {
      df$estimate[flip] <- -df$estimate[flip]
      if (all(c("lwr","upr") %in% names(df))) {
      lwr_new <- -df$upr[flip]; upr_new <- -df$lwr[flip]
      df$lwr[flip] <- pmin(lwr_new, upr_new, na.rm = TRUE)
      df$upr[flip] <- pmax(lwr_new, upr_new, na.rm = TRUE)
      }
      df$contrast[flip] <- paste(rhs[flip], "- CON")
    }
    df$pair <- dplyr::case_when(
      grepl("^RES\\s*-\\s*CON$", df$contrast) ~ "RES-CON",
      grepl("^SUS\\s*-\\s*CON$", df$contrast) ~ "SUS-CON",
      grepl("^SUS\\s*-\\s*RES$", df$contrast) ~ "SUS-RES",
      TRUE ~ NA_character_
    )
    df %>%
      dplyr::filter(pair %in% c("RES-CON","SUS-CON","SUS-RES")) %>%
      dplyr::mutate(pair = factor(pair, levels = c("RES-CON","SUS-CON","SUS-RES")))
    }

    contr_d1 <- canonicalize_contrast_d1(day1_res$phase_contr)

    if (nrow(contr_d1) > 0 && all(c("estimate","lwr","upr","p.adjust","d") %in% names(contr_d1))) {

    p_forest_d1 <- ggplot(
      contr_d1,
      aes(x = estimate, y = pair, xmin = lwr, xmax = upr, color = pair)
    ) +
      geom_vline(xintercept = 0, linetype = "solid", color = "grey80", linewidth = 0.3) +
      geom_pointrange(size = 0.5, linewidth = 1.1) +
      geom_text(
      aes(x = upr,
        label = paste0(
          dplyr::case_when(
          p.adjust < 0.001 ~ "***",
          p.adjust < 0.01  ~ "**",
          p.adjust < 0.05  ~ "*",
          TRUE ~ ""
          ),
          " d=", sprintf("%.2f", d)
        )),
      hjust = -0.1, size = 2.8, color = "black", show.legend = FALSE
      ) +
      scale_color_manual(values = pair_colors) +
      labs(
      title    = paste0(metric_name, ": Phase-average contrasts — First Active Phase (CC", first_change, ")"),
      x        = paste0("\u0394 ", metric_name, " (a.u.)"),
      y        = NULL,
      color    = "Contrast"
      ) +
      theme_classic(base_size = 10) +
      theme(
      strip.background  = element_blank(),
      strip.text        = element_text(face = "bold"),
      axis.line.y       = element_blank(),
      axis.ticks.y      = element_blank(),
      axis.text.y       = element_blank(),
      legend.position   = "bottom",
      legend.title      = element_blank(),
      plot.title        = element_text(face = "bold", hjust = 0.5, size = 9),
      plot.margin       = margin(5, 50, 5, 5)
      ) +
      scale_x_continuous(expand = expansion(mult = c(0.1, 0.55)))

    if ("Sex" %in% names(contr_d1) && dplyr::n_distinct(contr_d1$Sex) > 1) {
      p_forest_d1 <- p_forest_d1 + facet_wrap(~ Sex)
    }

    ggsave(file.path(dirs$plots_day1,
             paste0("forest_day1_", run_scope_d1, ".svg")),
         p_forest_d1, width = 6, height = 3.5)
    }
  }

  # --- EMM dot plot (means per group) ---
  if (!is.null(day1_res$phase_means) && nrow(day1_res$phase_means) > 0) {
    means_d1 <- day1_res$phase_means %>%
    dplyr::mutate(Group = factor(Group, levels = c("CON","RES","SUS")))

    lwr_col <- grep("^lwr$|^lower|^LCL", names(means_d1), value = TRUE)
    upr_col <- grep("^upr$|^upper|^UCL", names(means_d1), value = TRUE)

    if (length(lwr_col) && length(upr_col)) {
    means_d1 <- means_d1 %>%
      dplyr::rename(lwr_plot = !!lwr_col[1], upr_plot = !!upr_col[1])
    } else if ("SE_avg" %in% names(means_d1)) {
    means_d1 <- means_d1 %>%
      dplyr::mutate(lwr_plot = emmean_avg - 1.96 * SE_avg,
            upr_plot = emmean_avg + 1.96 * SE_avg)
    } else {
    means_d1$lwr_plot <- NA_real_
    means_d1$upr_plot <- NA_real_
    }

    p_emm_d1 <- ggplot(means_d1,
              aes(x = Group, y = emmean_avg,
                ymin = lwr_plot, ymax = upr_plot,
                color = Group)) +
    geom_pointrange(size = 0.7, linewidth = 0.9) +
    scale_color_manual(values = group_colors) +
    labs(
      title = paste0(metric_name, ": Estimated Phase Average — First Active Phase (CC", first_change, ")"),
      x = "Group",
      y = paste0("EMM ", metric_name, " (± 95% CI)")
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 9)
    )

    if ("Sex" %in% names(means_d1) && dplyr::n_distinct(means_d1$Sex) > 1) {
    p_emm_d1 <- p_emm_d1 + facet_wrap(~ Sex)
    }

    ggsave(file.path(dirs$plots_day1,
             paste0("dot_emm_day1_", run_scope_d1, ".svg")),
       p_emm_d1, width = 5, height = 4)
  }

  # --- Animal-level jitter + boxplot ---
  p_jitter_d1 <- ggplot(
    dsub_plot %>%
    dplyr::group_by(AnimalNum, Group, Sex) %>%
    dplyr::summarise(mean_val = mean(.data[[metric_name]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Group = factor(Group, levels = c("CON","RES","SUS"))),
    aes(x = Group, y = mean_val, color = Group)
  ) +
    geom_boxplot(aes(fill = Group), alpha = 0.25, outlier.shape = NA, width = 0.45) +
    geom_jitter(width = 0.12, size = 2.2, alpha = 0.7) +
    scale_color_manual(values = group_colors) +
    scale_fill_manual(values  = group_colors) +
    labs(
    title = paste0(metric_name, ": Animal means — First Active Phase (CC", first_change, ")"),
    x = "Group",
    y = paste0("Mean ", metric_name, " [a.u.]")
    ) +
    theme_classic(base_size = 11) +
    theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 9)
    )

  if (includeSex && "Sex" %in% names(dsub_plot) && dplyr::n_distinct(dsub_plot$Sex) > 1) {
    p_jitter_d1 <- p_jitter_d1 + facet_wrap(~ Sex)
  }

  ggsave(file.path(dirs$plots_day1,
           paste0("jitter_animals_day1_", run_scope_d1, ".svg")),
       p_jitter_d1, width = 5, height = 4.5)

  cat(sprintf("  ✓ Day 1 active phase plots saved to: %s\n", dirs$plots_day1))
  }

  # -------------------------------------------------
  # Run Day 1 Active Phase Analysis
  # -------------------------------------------------
  cat("\n--- Running Day 1 First Active Phase Analysis ---\n")
  day1_res <- fit_model_day1_active()

  if (!is.null(day1_res) && nrow(day1_res$phase_contr) > 0) {
    day1_res$phase_contr <- finalize_phase_contrast_pvals(day1_res$phase_contr, method = "BH")
  }

  if (!is.null(day1_res)) {
  # Save tables
  if (nrow(day1_res$summary) > 0) {
    write_pub_xlsx(day1_res$summary,
               file.path(dirs$tables_day1, paste0("lme_fixedEffects_day1.xlsx")),
               sheet_name = "Fixed_Effects", p_cols = c("p.value", "p.adj_FDR"))
  }
  if (nrow(day1_res$emmeans) > 0) {
    write_pub_xlsx(day1_res$emmeans,
               file.path(dirs$tables_day1, paste0("emm_contrasts_day1.xlsx")),
               sheet_name = "EMM_Contrasts", p_cols = c("p.adjust_primary", "p.adjust", "p.adjust_strict", "p.value"))
  }
  if (nrow(day1_res$phase_means) > 0) {
    write_pub_xlsx(day1_res$phase_means,
               file.path(dirs$tables_day1, paste0("phaseAvg_emm_day1.xlsx")),
               sheet_name = "Phase_Avg_EMM", p_cols = c("p.adjust"))
  }
  if (nrow(day1_res$phase_contr) > 0) {
    write_pub_xlsx(day1_res$phase_contr,
               file.path(dirs$tables_day1, paste0("phaseAvg_contrasts_day1.xlsx")),
               sheet_name = "Phase_Avg_Contrasts", p_cols = c("p.adjust_primary", "p.adjust", "p.adjust_strict", "p.value_raw"))
  }
  cat(sprintf("  ✓ Day 1 active phase tables saved to: %s\n", dirs$tables_day1))

  # Generate plots
  plot_day1_active(day1_res)
  } else {
  cat("  - Day 1 active phase analysis skipped (no data).\n")
  }

  # Main (full phase duration)
  main_res <- fit_model_slices(window_label = "all")
  summary_results <- main_res$summary
  emmeans_results <- main_res$emmeans
  phase_avg_means_tbl <- main_res$phase_means
  phase_avg_contr_tbl <- main_res$phase_contr

  if (nrow(phase_avg_contr_tbl) > 0) {
    phase_avg_contr_tbl <- finalize_phase_contrast_pvals(phase_avg_contr_tbl, method = "BH")
  }

  # Optional: first/last 2h per phase
  first2h_res <- NULL
  last2h_res  <- NULL
  if (isTRUE(includePhase)) {
  first2h_res <- fit_model_slices(window_label = "first2h")
  last2h_res  <- fit_model_slices(window_label = "last2h")
  if (!is.null(first2h_res) && nrow(first2h_res$phase_contr) > 0) {
    first2h_res$phase_contr <- finalize_phase_contrast_pvals(first2h_res$phase_contr, method = "BH")
  }
  if (!is.null(last2h_res) && nrow(last2h_res$phase_contr) > 0) {
    last2h_res$phase_contr <- finalize_phase_contrast_pvals(last2h_res$phase_contr, method = "BH")
  }
  }

  run_scope <- scope_all(includeChange, includeSex, includePhase)
  summary_xlsx <- file.path(dirs$tables_sum, paste0("lme_fixedEffects_", run_scope, ".xlsx"))
  emmeans_xlsx <- file.path(dirs$tables_emm, paste0("emm_contrasts_", run_scope, ".xlsx"))
  
  write_pub_xlsx(summary_results, summary_xlsx, sheet_name = "Fixed_Effects",
                 p_cols = c("p.value", "p.adj_FDR"))

  # Add metadata sheet to fixed effects xlsx
  add_lme_metadata_sheet <- function(xlsx_path, summary_tbl, metric_name, run_scope) {
    if (!file.exists(xlsx_path) || nrow(summary_tbl) == 0) return(invisible(NULL))

    # Prepare metadata information
    model_types <- if ("model_engine" %in% names(summary_tbl)) unique(summary_tbl$model_engine) else "lmerTest::lmer"
    change_vals <- if ("Change" %in% names(summary_tbl)) unique(as.character(summary_tbl$Change)) else run_scope
    sex_vals    <- if ("Sex"    %in% names(summary_tbl)) unique(as.character(summary_tbl$Sex))    else "all"
    phase_vals  <- if ("Phase"  %in% names(summary_tbl)) unique(as.character(summary_tbl$Phase))  else "all"

    metadata_info <- tibble::tibble(
      Information_Type = c(
        "Analysis Type",
        "Metric",
        "Scope",
        "Model Engine(s)",
        "Formula",
        "Reference Group",
        "Data - Changes Included",
        "Data - Sexes Included",
        "Data - Phases Included",
        "Coefficient Interpretation - Intercept",
        "Coefficient Interpretation - Group Effects",
        "Coefficient Interpretation - Interaction Terms",
        "Statistical Testing",
        "Multiple Comparisons",
        "Publication Notes"
      ),
      Value = c(
        "Mixed-effects linear regression (lmerTest)",
        metric_name,
        run_scope,
        paste0(sort(unique(model_types)), collapse = "; "),
        "response ~ Group (+ optional: TH_scaled, Group:TH_scaled, + (1|AnimalNum))",
        "CON (Control group)",
        paste0(sort(change_vals), collapse = ", "),
        paste0(sort(sex_vals), collapse = ", "),
        paste0(sort(phase_vals), collapse = ", "),
        paste0("Baseline predicted value for CON group at TH_scaled=0. ",
               "Represents intercept of the mixed model."),
        paste0("Effect of each group (RES, SUS) relative to CON (reference). ",
               "Positive value = higher response than CON; negative = lower response."),
        paste0("If present (Group:TH_scaled terms): How the time trajectory differs by group. ",
               "E.g., GroupRES:TH_scaled = difference in slope between RES and CON."),
        paste0("All p-values from lmerTest (Satterthwaite DF approximation). ",
               "p.adj_FDR: Benjamini-Hochberg false discovery rate correction within family."),
        "For pairwise group comparisons, see emmeans_contrasts table.",
        paste0("Individual fixed effects show main effects and baseline differences vs. CON. ",
               "Positive estimate = higher response than reference group.")
      )
    )

    # Load workbook and add metadata sheet
    wb <- openxlsx::loadWorkbook(xlsx_path)

    if (!"Metadata & Interpretation" %in% openxlsx::sheets(wb)) {
      openxlsx::addWorksheet(wb, "Metadata & Interpretation")
    }

    openxlsx::writeData(wb, "Metadata & Interpretation", metadata_info, startRow = 1, startCol = 1)

    header_style <- openxlsx::createStyle(
      fgFill = "#E7E6E6", textDecoration = "bold", border = "TopBottomLeftRight"
    )
    openxlsx::addStyle(wb, "Metadata & Interpretation", header_style, rows = 1, cols = 1:2)
    openxlsx::setColWidths(wb, "Metadata & Interpretation", cols = 1:2, widths = c(35, 65))

    openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
    invisible(NULL)
  }

  # Apply metadata to main fixed effects file
  add_lme_metadata_sheet(summary_xlsx, summary_results, metric_name, run_scope)
  
  # EMM contrasts: p_cols order determines highlighting priority
  # - p.adjust_primary: focused on SUS-CON, SUS-RES (highlighted green if p<0.05)
  # - p.adjust: family-level BH (fallback for RES-CON)
  # - p.adjust_strict: family-level Holm (conservative sensitivity check)
  # - p.value: raw p-values (for transparency)
  write_pub_xlsx(emmeans_results, emmeans_xlsx, sheet_name = "EMM_Contrasts",
                 p_cols = c("p.adjust_primary", "p.adjust", "p.adjust_strict", "p.value"))

  # Append stratification summary as a second sheet
  stratif_summary <- summarize_contrast_stratification(emmeans_results)
  openxlsx::addWorksheet(wb <- openxlsx::loadWorkbook(emmeans_xlsx), "Summary_Stats")
  openxlsx::writeData(wb, "Summary_Stats", stratif_summary)
  openxlsx::saveWorkbook(wb, emmeans_xlsx, overwrite = TRUE)

  # Save optional first2h tables
  if (!is.null(first2h_res)) {
  first2h_xlsx <- file.path(dirs$tables_sum, paste0("lme_fixedEffects_", run_scope, "_first2h.xlsx"))
  write_pub_xlsx(
    first2h_res$summary,
    first2h_xlsx,
    sheet_name = "Fixed_Effects", p_cols = c("p.value", "p.adj_FDR")
  )
  add_lme_metadata_sheet(first2h_xlsx, first2h_res$summary, metric_name, paste0(run_scope, "_first2h"))
  
  write_pub_xlsx(
    first2h_res$emmeans,
    file.path(dirs$tables_emm, paste0("emm_contrasts_", run_scope, "_first2h.xlsx")),
    sheet_name = "EMM_Contrasts", p_cols = c("p.adjust_primary", "p.adjust", "p.adjust_strict", "p.value")
  )
  }

  # Save optional last2h tables
  if (!is.null(last2h_res)) {
  last2h_xlsx <- file.path(dirs$tables_sum, paste0("lme_fixedEffects_", run_scope, "_last2h.xlsx"))
  write_pub_xlsx(
    last2h_res$summary,
    last2h_xlsx,
    sheet_name = "Fixed_Effects", p_cols = c("p.value", "p.adj_FDR")
  )
  add_lme_metadata_sheet(last2h_xlsx, last2h_res$summary, metric_name, paste0(run_scope, "_last2h"))
  
  write_pub_xlsx(
    last2h_res$emmeans,
    file.path(dirs$tables_emm, paste0("emm_contrasts_", run_scope, "_last2h.xlsx")),
    sheet_name = "EMM_Contrasts", p_cols = c("p.adjust_primary", "p.adjust", "p.adjust_strict", "p.value")
  )
  }

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
  dplyr::distinct(Change, Sex, Phase, Group, window, .keep_all = TRUE)
  }
  if (nrow(phase_avg_contr_tbl)) {
  phase_avg_contr_tbl <- phase_avg_contr_tbl %>%
  dplyr::mutate(
    Change = factor(Change, levels = unique(Change)),
    Phase  = factor(Phase,  levels = c("Active","Inactive")),
    Sex    = factor(Sex,    levels = unique(Sex))
  )
  }

  # Save default (all window) for downstream compatibility
  write_pub_xlsx(
    dplyr::filter(phase_avg_means_tbl, window == "all"),
    file.path(dirs$tables_emm, paste0("phaseAvg_emm_", run_scope, ".xlsx")),
    sheet_name = "Phase_Avg_EMM", p_cols = c("p.adjust")
  )
  write_pub_xlsx(
    dplyr::filter(phase_avg_contr_tbl, window == "all"),
    file.path(dirs$tables_emm, paste0("phaseAvg_contrasts_", run_scope, ".xlsx")),
    sheet_name = "Phase_Avg_Contrasts", p_cols = c("p.adjust", "p.adjust_primary", "p.adjust_strict", "p.value_raw")
  )

  # Save optional first2h phase-average tables
  if (!is.null(first2h_res)) {
  if (nrow(first2h_res$phase_means) > 0) {
  write_pub_xlsx(
    first2h_res$phase_means,
    file.path(dirs$tables_emm, paste0("phaseAvg_emm_", run_scope, "_first2h.xlsx")),
    sheet_name = "Phase_Avg_EMM", p_cols = c("p.adjust")
  )
  }
  if (nrow(first2h_res$phase_contr) > 0) {
  write_pub_xlsx(
    first2h_res$phase_contr,
    file.path(dirs$tables_emm, paste0("phaseAvg_contrasts_", run_scope, "_first2h.xlsx")),
    sheet_name = "Phase_Avg_Contrasts", p_cols = c("p.adjust", "p.adjust_primary", "p.adjust_strict", "p.value_raw")
  )
  }
  }

  # Save optional last2h phase-average tables
  if (!is.null(last2h_res)) {
  if (nrow(last2h_res$phase_means) > 0) {
  write_pub_xlsx(
    last2h_res$phase_means,
    file.path(dirs$tables_emm, paste0("phaseAvg_emm_", run_scope, "_last2h.xlsx")),
    sheet_name = "Phase_Avg_EMM", p_cols = c("p.adjust")
  )
  }
  if (nrow(last2h_res$phase_contr) > 0) {
  write_pub_xlsx(
    last2h_res$phase_contr,
    file.path(dirs$tables_emm, paste0("phaseAvg_contrasts_", run_scope, "_last2h.xlsx")),
    sheet_name = "Phase_Avg_Contrasts", p_cols = c("p.adjust", "p.adjust_primary", "p.adjust_strict", "p.value_raw")
  )
  }
  }

  canonicalize_contrast <- function(df, keep_pairs = c("RES-CON","SUS-CON","SUS-RES")) {
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
  grepl("^SUS\\s*-\\s*RES$", df$contrast) ~ "SUS-RES",
  TRUE ~ NA_character_
  )

  df <- df %>% dplyr::filter(pair %in% keep_pairs)

  if (!all(c("lwr","upr") %in% names(df)) && all(c("SE","df") %in% names(df))) {
  df$tcrit <- qt(0.975, df = df$df)
  df$lwr   <- df$estimate - df$tcrit * df$SE
  df$upr   <- df$estimate + df$tcrit * df$SE
  }

  df %>% dplyr::distinct(Change, Sex, Phase, pair, window, .keep_all = TRUE)
  }

  phase_avg_contr_tbl <- canonicalize_contrast(phase_avg_contr_tbl, keep_pairs = c("RES-CON","SUS-CON","SUS-RES"))

  plot_phase_forest <- function(df_in, suffix = "", window_label_text = "all phase duration") {
  if (nrow(df_in) == 0) return(invisible(NULL))

  plot_data <- df_in %>%
  dplyr::mutate(
    Change = dplyr::recode(as.character(Change), "CC1" = "1", "CC2" = "2", "CC3" = "3", "CC4" = "4"),
    Change = factor(Change, levels = sort(unique(Change))),
    Sex = factor(Sex, levels = c("m", "f"), labels = c("Male", "Female")),
    Phase = factor(Phase, levels = c("Active", "Inactive")),
    pair = factor(pair, levels = c("RES-CON", "SUS-CON", "SUS-RES"))
  ) %>%
  dplyr::filter(!is.na(estimate), !is.na(lwr), !is.na(upr), !is.na(p.adjust), !is.na(d))

  if (nrow(plot_data) == 0) return(invisible(NULL))

  forest_ultra <- ggplot(
  plot_data,
  aes(x = estimate, y = pair, xmin = lwr, xmax = upr, color = pair)
  ) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey85", linewidth = 0.3) +
  geom_pointrange(size = 0.4, linewidth = 1.0) +
  geom_text(
    aes(
    x = upr,
    label = paste0(
    ifelse(p.adjust < 0.001, "***", ifelse(p.adjust < 0.01, "**", ifelse(p.adjust < 0.05, "*", ""))),
    " d=", sprintf("%.2f", d)
    )
    ),
    hjust = -0.1,
    size = 2.2,
    color = "black",
    show.legend = FALSE
  ) +
  facet_grid(Sex ~ Phase + Change, scales = "free_x",
       labeller = labeller(Change = function(x) rep("", length(x)))) +
  scale_color_manual(values = pair_colors) + 
  labs(
    title = paste0(metric_name, ": phase-average contrasts (", window_label_text, ")"),
    subtitle = "Stats shown are emmeans phase-average contrasts (adjusted p), not fixed-effects coefficient tests",
    x = paste0("\u0394 ", metric_name, " (a.u.)"),
    y = NULL
  ) +
  theme_classic(base_size = 8) +
  theme(
    panel.spacing = unit(0.8, "lines"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 8),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.3, "cm"),
    plot.margin = margin(5, 45, 5, 5),
    panel.grid.major.x = element_line(color = "grey98"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 9),
    plot.subtitle = element_text(hjust = 0.5, size = 7)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.6)))

  ggsave(
  file.path(dirs$plots_pub, paste0("forest_phaseAvg_", run_scope, suffix, ".svg")),
  forest_ultra,
  width = if (isTRUE(includeChange)) 3 * 2.5 else 2.5,
  height = 3,
  device = "svg"
  )
  }

  # Plot default + optional first2h + optional last2h
  plot_phase_forest(
    dplyr::filter(phase_avg_contr_tbl, window == "all"),
    suffix = "",
    window_label_text = "all phase duration"
  )
  if (!is.null(first2h_res) && nrow(first2h_res$phase_contr) > 0) {
  first2h_plot_tbl <- canonicalize_contrast(first2h_res$phase_contr, keep_pairs = c("RES-CON","SUS-CON","SUS-RES"))
  plot_phase_forest(
    first2h_plot_tbl,
    suffix = "_first2h",
    window_label_text = "first 2h per phase"
  )
  }
  if (!is.null(last2h_res) && nrow(last2h_res$phase_contr) > 0) {
  last2h_plot_tbl <- canonicalize_contrast(last2h_res$phase_contr, keep_pairs = c("RES-CON","SUS-CON","SUS-RES"))
  plot_phase_forest(
    last2h_plot_tbl,
    suffix = "_last2h",
    window_label_text = "last 2h per phase"
  )
  }

  # -------------------------------------------------
  # Save to proper directories + generate manifest & inventory
  # -------------------------------------------------
  if (nrow(phase_avg_means_tbl) > 0) {
  write_pub_xlsx(
    dplyr::filter(phase_avg_means_tbl, window == "all"),
    file.path(dirs$tables_phase, paste0("phaseAvg_emm_", run_scope, ".xlsx")),
    sheet_name = "Phase_Avg_EMM", p_cols = c("p.adjust")
  )
  }

  if (nrow(phase_avg_contr_tbl) > 0) {
  write_pub_xlsx(
    dplyr::filter(phase_avg_contr_tbl, window == "all"),
    file.path(dirs$tables_phase, paste0("phaseAvg_contrasts_", run_scope, ".xlsx")),
    sheet_name = "Phase_Avg_Contrasts", p_cols = c("p.adjust", "p.adjust_primary", "p.adjust_strict", "p.value_raw")
  )
  }

  if (!is.null(first2h_res)) {
  if (nrow(first2h_res$phase_means) > 0) {
  write_pub_xlsx(
    first2h_res$phase_means,
    file.path(dirs$tables_phase, paste0("phaseAvg_emm_", run_scope, "_first2h.xlsx")),
    sheet_name = "Phase_Avg_EMM", p_cols = c("p.adjust")
  )
  }
  if (nrow(first2h_res$phase_contr) > 0) {
  write_pub_xlsx(
    first2h_res$phase_contr,
    file.path(dirs$tables_phase, paste0("phaseAvg_contrasts_", run_scope, "_first2h.xlsx")),
    sheet_name = "Phase_Avg_Contrasts", p_cols = c("p.adjust", "p.adjust_primary", "p.adjust_strict", "p.value_raw")
  )
  }
  }

  if (!is.null(last2h_res)) {
  if (nrow(last2h_res$phase_means) > 0) {
  write_pub_xlsx(
    last2h_res$phase_means,
    file.path(dirs$tables_phase, paste0("phaseAvg_emm_", run_scope, "_last2h.xlsx")),
    sheet_name = "Phase_Avg_EMM", p_cols = c("p.adjust")
  )
  }
  if (nrow(last2h_res$phase_contr) > 0) {
  write_pub_xlsx(
    last2h_res$phase_contr,
    file.path(dirs$tables_phase, paste0("phaseAvg_contrasts_", run_scope, "_last2h.xlsx")),
    sheet_name = "Phase_Avg_Contrasts", p_cols = c("p.adjust", "p.adjust_primary", "p.adjust_strict", "p.value_raw")
  )
  }
  }

  if (nrow(emmeans_results) > 0) {
  contrasts_only <- emmeans_results %>% dplyr::filter(!is.na(contrast))
  if (nrow(contrasts_only) > 0) {
  write_pub_xlsx(
    contrasts_only,
    file.path(dirs$tables_contr, paste0("contrasts_", run_scope, ".xlsx")),
    sheet_name = "EMM_Contrasts", p_cols = c("p.adjust", "p.adjust_primary", "p.adjust_strict", "p.value")
  )
  }
  }

  if (!is.null(first2h_res) && nrow(first2h_res$emmeans) > 0) {
  contrasts_only_2h <- first2h_res$emmeans %>% dplyr::filter(!is.na(contrast))
  if (nrow(contrasts_only_2h) > 0) {
  write_pub_xlsx(
    contrasts_only_2h,
    file.path(dirs$tables_contr, paste0("contrasts_", run_scope, "_first2h.xlsx")),
    sheet_name = "EMM_Contrasts", p_cols = c("p.adjust", "p.adjust_primary", "p.adjust_strict", "p.value")
  )
  }
  }

  if (!is.null(last2h_res) && nrow(last2h_res$emmeans) > 0) {
  contrasts_only_last2h <- last2h_res$emmeans %>% dplyr::filter(!is.na(contrast))
  if (nrow(contrasts_only_last2h) > 0) {
  write_pub_xlsx(
    contrasts_only_last2h,
    file.path(dirs$tables_contr, paste0("contrasts_", run_scope, "_last2h.xlsx")),
    sheet_name = "EMM_Contrasts", p_cols = c("p.adjust", "p.adjust_primary", "p.adjust_strict", "p.value")
  )
  }
  }

  # Inventory
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
    mean_ActivityIndex = mean(dsub$ActivityIndex, na.rm = TRUE),
    sd_ActivityIndex = sd(dsub$ActivityIndex, na.rm = TRUE),
    stringsAsFactors = FALSE
    )
    }
  }
  }
  }

  inventory_df <- dplyr::bind_rows(inventory_list)
  if (nrow(inventory_df) > 0) {
  write_pub_xlsx(inventory_df,
  file.path(dirs$tables_inv, paste0("inventory_", run_scope, ".xlsx")),
  sheet_name = "Results", p_cols = character(0))
  }

  # Manifest
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
  n_models_fitted = ifelse(nrow(summary_results) > 0,
         nrow(summary_results) / length(unique(summary_results$Fixed_effect)), 0),
  n_models_fitted_first2h = ifelse(!is.null(first2h_res) && nrow(first2h_res$summary) > 0,
           nrow(first2h_res$summary) / length(unique(first2h_res$summary$Fixed_effect)), 0),
  n_models_fitted_last2h = ifelse(!is.null(last2h_res) && nrow(last2h_res$summary) > 0,
          nrow(last2h_res$summary) / length(unique(last2h_res$summary$Fixed_effect)), 0),
  n_models_fitted_day1_active = ifelse(!is.null(day1_res) && nrow(day1_res$summary) > 0,
          nrow(day1_res$summary) / max(1, length(unique(day1_res$summary$Fixed_effect))), 0),
  n_contrasts = nrow(emmeans_results),
  n_contrasts_first2h = ifelse(!is.null(first2h_res), nrow(first2h_res$emmeans), 0),
  n_contrasts_last2h = ifelse(!is.null(last2h_res), nrow(last2h_res$emmeans), 0),
  n_contrasts_day1_active = ifelse(!is.null(day1_res), nrow(day1_res$emmeans), 0),
  n_phase_avg_contrasts = nrow(dplyr::filter(phase_avg_contr_tbl, window == "all")),
  n_phase_avg_contrasts_first2h = ifelse(!is.null(first2h_res), nrow(first2h_res$phase_contr), 0),
  n_phase_avg_contrasts_last2h = ifelse(!is.null(last2h_res), nrow(last2h_res$phase_contr), 0),
  n_phase_avg_contrasts_day1_active = ifelse(!is.null(day1_res), nrow(day1_res$phase_contr), 0),
  total_observations = nrow(data),
  total_animals = dplyr::n_distinct(data$AnimalNum),
  data_file = data_file,
  results_dir = results_dir,
  stringsAsFactors = FALSE
  )

  write_pub_xlsx(manifest,
  file.path(dirs$tables_man, paste0("manifest_", run_scope, ".xlsx")),
  sheet_name = "Results", p_cols = character(0))

  # Full artifact snapshot for reproducibility
  snapshot_file_main <- file.path(dirs$artifacts, paste0("snapshot_main_", run_scope, ".rds"))
  saveRDS(
    list(
      metric = metric_name,
      run_scope = run_scope,
      includeChange = includeChange,
      includeSex = includeSex,
      includePhase = includePhase,
      dirs = dirs,
      summary_results = summary_results,
      emmeans_results = emmeans_results,
      phase_avg_means = phase_avg_means_tbl,
      phase_avg_contrasts = phase_avg_contr_tbl,
      day1_results = day1_res,
      first2h_results = first2h_res,
      last2h_results = last2h_res,
      manifest = manifest,
      inventory = inventory_df
    ),
    snapshot_file_main,
    compress = "gzip"   # faster than xz with only slightly larger files
  )

  cat(sprintf("\n✓ Additional tables saved:\n"))
  cat(sprintf("  - Inventory: %s\n", dirs$tables_inv))
  cat(sprintf("  - Manifest: %s\n", dirs$tables_man))
  cat(sprintf("  - Phase average: %s\n", dirs$tables_phase))
  cat(sprintf("  - Contrasts: %s\n", dirs$tables_contr))
  cat(sprintf("  - Day 1 active: %s\n", dirs$tables_day1))
  if (!is.null(first2h_res)) cat(sprintf("  - First 2h outputs saved with suffix: _first2h\n"))
  if (!is.null(last2h_res))  cat(sprintf("  - Last 2h outputs saved with suffix: _last2h\n"))
  }

  # -------------------------------------------------
  # Additional plots (volcano, panels A/B/C, change-development line plot)
  # -------------------------------------------------
  
  # Load phase-average result tables
  phase_contr_file <- file.path(dirs$tables_emm, paste0("phaseAvg_contrasts_", run_scope, ".xlsx"))
  phase_means_file <- file.path(dirs$tables_emm, paste0("phaseAvg_emm_", run_scope, ".xlsx"))

  phase_avg_contr_tbl <- if (file.exists(phase_contr_file)) {
  openxlsx::read.xlsx(phase_contr_file)
  } else {
  tibble::tibble()
  }

  phase_avg_means_tbl <- if (file.exists(phase_means_file)) {
  openxlsx::read.xlsx(phase_means_file)
  } else {
  tibble::tibble()
  }
  
  # Convert to factors only if data is present
  if (nrow(phase_avg_contr_tbl) > 0 && has_cols(phase_avg_contr_tbl, c("contrast","Change","Sex","Phase"))) {
  phase_avg_contr_tbl <- phase_avg_contr_tbl %>%
  mutate(
    contrast = factor(contrast, levels = contrast_config$display_names),
    Change = factor(Change),
    Sex = factor(Sex),
    Phase = factor(Phase, levels = c("Active", "Inactive"))
  )
  }

  if (nrow(phase_avg_means_tbl) > 0 && has_cols(phase_avg_means_tbl, c("Group","Phase","Change","Sex"))) {
  phase_avg_means_tbl <- phase_avg_means_tbl %>%
  mutate(
    Group = factor(Group, levels = c("CON","RES","SUS")),
    Phase = factor(Phase, levels = c("Active", "Inactive")),
    Change = factor(Change),
    Sex = factor(Sex)
  )
  }

  # ---------- EMMEANS development over cage changes ----------
  make_change_levels <- function(x) {
  lev <- unique(as.character(x))
  num <- suppressWarnings(as.numeric(stringr::str_extract(lev, "\\d+")))
  if (all(!is.na(num))) lev[order(num)] else sort(lev)
  }

  plot_emmeans_change_development <- function(df_means, suffix = "", subtitle_txt = "All phase duration") {
  if (!isTRUE(includeChange)) return(invisible(NULL))
  if (nrow(df_means) == 0) return(invisible(NULL))
  if (!all(c("Change","Group","emmean_avg") %in% names(df_means))) return(invisible(NULL))
  if (dplyr::n_distinct(df_means$Change) < 2) return(invisible(NULL))

  d <- df_means %>%
  dplyr::filter(!is.na(Change), !is.na(Group), !is.na(emmean_avg)) %>%
  dplyr::mutate(
  Change_chr = as.character(Change),
  Change = factor(Change_chr, levels = make_change_levels(Change_chr)),
  Change_num = as.numeric(Change),
  Group = factor(Group, levels = c("CON","RES","SUS"))
  )

  if (!("lwr" %in% names(d) && "upr" %in% names(d))) {
  if ("SE_avg" %in% names(d)) {
  d <- d %>%
    dplyr::mutate(
    lwr = emmean_avg - 1.96 * SE_avg,
    upr = emmean_avg + 1.96 * SE_avg
    )
  } else {
  d$lwr <- NA_real_
  d$upr <- NA_real_
  }
  }

  p <- ggplot(d, aes(x = Change_num, y = emmean_avg, color = Group, group = Group)) +
  geom_line(linewidth = 0.7, alpha = 0.95) +
  geom_point(size = 1.8, stroke = 0.2) +
  scale_color_manual(values = group_colors) +
  scale_x_continuous(
    breaks = seq_along(levels(d$Change)),
    labels = levels(d$Change),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    title = paste0(metric_name, ": estimated development over cage changes"),
    subtitle = subtitle_txt,
    x = "Cage change",
    y = paste0("Estimated ", metric_name, " (EMM)")
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 8, base_family = "Arial") +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size = 7),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7, color = "black"),
    axis.line = element_line(linewidth = 0.4, color = "black"),
    axis.ticks = element_line(linewidth = 0.4, color = "black"),
    axis.ticks.length = unit(1.5, "mm"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
    plot.subtitle = element_text(hjust = 0.5, size = 7),
    plot.margin = margin(4, 6, 4, 4)
  )

  if (all(!is.na(d$lwr)) && all(!is.na(d$upr))) {
  p <- p +
  geom_ribbon(
    aes(ymin = lwr, ymax = upr, fill = Group),
    alpha = 0.14, color = NA, inherit.aes = TRUE
  ) +
  scale_fill_manual(values = group_colors, guide = "none")
  }

  # Facet only if present
  if ("Phase" %in% names(d) && "Sex" %in% names(d) &&
  dplyr::n_distinct(d$Phase) > 1 && dplyr::n_distinct(d$Sex) > 1) {
  p <- p + facet_grid(Phase ~ Sex, scales = "free_y")
  } else if ("Phase" %in% names(d) && dplyr::n_distinct(d$Phase) > 1) {
  p <- p + facet_wrap(~ Phase, scales = "free_y")
  } else if ("Sex" %in% names(d) && dplyr::n_distinct(d$Sex) > 1) {
  p <- p + facet_wrap(~ Sex, scales = "free_y")
  }

  ggsave(
  file.path(dirs$plots_pub, paste0("line_emm_byCC_", run_scope, suffix, ".svg")),
  p, width = 6, height = 4
  )
  ggsave(
  file.path(dirs$plots_fixed, paste0("line_emm_byCC_", run_scope, suffix, ".svg")),
  p, width = 6, height = 4
  )
  }

  if (isTRUE(includeChange)) {
  plot_emmeans_change_development(phase_avg_means_tbl, suffix = "", subtitle_txt = "All phase duration")

  f_first2h <- file.path(dirs$tables_emm, paste0("phaseAvg_emm_", run_scope, "_first2h.xlsx"))
  if (file.exists(f_first2h)) {
  tmp_first <- openxlsx::read.xlsx(f_first2h)
  if (nrow(tmp_first) > 0) {
  plot_emmeans_change_development(tmp_first, suffix = "_first2h", subtitle_txt = "First 2h per phase")
  }
  }

  f_last2h <- file.path(dirs$tables_emm, paste0("phaseAvg_emm_", run_scope, "_last2h.xlsx"))
  if (file.exists(f_last2h)) {
  tmp_last <- openxlsx::read.xlsx(f_last2h)
  if (nrow(tmp_last) > 0) {
  plot_emmeans_change_development(tmp_last, suffix = "_last2h", subtitle_txt = "Last 2h per phase")
  }
  }
  cat("✓ Added EMM line plot(s) across cage changes.\n")
  } else {
  cat("Skipping EMM change-development plot (includeChange = FALSE).\n")
  }

  m_first_day_active <- file.path(dirs$tables_day1, "emm_contrasts_day1.xlsx")
  if (file.exists(m_first_day_active)) { 
  tmp_day1 <- openxlsx::read.xlsx(m_first_day_active)
  if (nrow(tmp_day1) > 0) {
  plot_emmeans_change_development(tmp_day1, suffix = "_day1_active", subtitle_txt = "Day 1 active phase")
  cat("✓ Added EMM line plot for Day 1 active phase.\n")
  } else {
  cat("No data in Day 1 active EMM results; skipping line plot.\n")
  }
  } else {  
  cat("No Day 1 active EMM results file found; skipping line plot.\n")
  }

  f_first_day_active <- file.path(dirs$tables_day1, "phaseAvg_emm_day1.xlsx")
  if (file.exists(f_first_day_active)) {
  tmp_day1_means <- openxlsx::read.xlsx(f_first_day_active)
  if (nrow(tmp_day1_means) > 0) { 
  plot_emmeans_change_development(tmp_day1_means, suffix = "_day1_active_means", subtitle_txt = "Day 1 active phase (means)")
  cat("✓ Added EMM line plot for Day 1 active phase means.\n")
  } else {
  cat("No data in Day 1 active phase means; skipping line plot.\n")
  }
  } else {
  cat("No Day 1 active phase means file found; skipping line plot.\n")
  }
  

  # ---------- END NEW BLOCK ----------

  # Volcano plot
  plot_volcano_phaseavg <- function(df, n_labels = 10, pval_thresh = 0.05, effect_thresh = 1) {
  req_cols <- c("p.adjust","estimate","contrast","Phase","Sex")
  if (nrow(df) == 0 || !has_cols(df, req_cols)) return(NULL)

  df <- df %>%
  mutate(
  sig = p.adjust < pval_thresh,
  large_effect = abs(estimate) > effect_thresh,
  label_flag = sig & large_effect,
  contrast = factor(contrast, levels = contrast_config$display_names),
  Phase = factor(Phase, levels = c("Active", "Inactive")),
  Sex = factor(Sex)
  )

  if (nrow(df) == 0 || dplyr::n_distinct(df$Sex) == 0) return(NULL)
  
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
  scale_color_manual(values = contrast_config$colors) +
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
  if (!is.null(volcano_plot_phaseavg)) {
  ggsave(file.path(dirs$plots_pub, paste0("volcano_phaseAvg_", run_scope, ".svg")), 
     volcano_plot_phaseavg, width = 8, height = 6)
  } else {
  cat("Skipping volcano plot: no valid phase-average contrast data.\n")
  }
  
  plot_forest_phaseavg <- function(df, contrast_levels = levels(df$contrast)) {
  req_cols <- c("contrast","estimate","lwr","upr","Phase","Change")
  if (nrow(df) == 0 || !has_cols(df, req_cols)) return(NULL)

  sub <- df %>% filter(contrast %in% contrast_levels)
  if (nrow(sub) == 0) return(NULL)

  ggplot(sub, aes(y = reorder(contrast, estimate), x = estimate, xmin = lwr, xmax = upr, color = contrast)) +
  geom_point(size = 3) +
  geom_errorbar(width = 0.2, orientation = "y") +
  facet_grid(Phase ~ Change) +
  theme_minimal(base_size = 14) +
  labs(x = "Estimate (phase avg, with 95% CI)", y = "Contrast", 
     title = paste("Forest Plot:", metric_name))
  }
  
  forest_plot_phaseavg <- plot_forest_phaseavg(phase_avg_contr_tbl)
  if (!is.null(forest_plot_phaseavg)) {
  ggsave(file.path(dirs$plots_pub, paste0("forest_phaseAvg_contrasts_", run_scope, ".svg")), 
     forest_plot_phaseavg, width = 9, height = 7)
  } else {
  cat("Skipping forest plot: no valid phase-average contrast data.\n")
  }
  
  plot_dot_estimation_phaseavg <- function(df) {
  req_cols <- c("Group","emmean_avg","Sex","lwr","upr","Phase","Change")
  if (nrow(df) == 0 || !has_cols(df, req_cols)) return(NULL)

  ggplot(df, aes(x = Group, y = emmean_avg, color = Sex)) +
  geom_point(position = position_jitter(width = 0.2), size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  facet_grid(Phase ~ Change) +
  theme_minimal(base_size = 14) +
  labs(x = "Group", y = "Estimated Phase Average", 
     title = paste("Dot/Estimation Plot:", metric_name))
  }

  dot_estimation_plot_phaseavg <- plot_dot_estimation_phaseavg(phase_avg_means_tbl)
  if (!is.null(dot_estimation_plot_phaseavg)) {
  ggsave(file.path(dirs$plots_pub, paste0("dot_phaseAvg_emm_", run_scope, ".svg")), 
     dot_estimation_plot_phaseavg, width = 8, height = 6)
  } else {
  cat("Skipping dot-estimation plot: no valid phase-average means data.\n")
  }
  
  # -------------------------------------------------
  # Panels A/B/C
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
  
  # Panel A
  file_panelA <- file.path(dirs$plots_pub, paste0("forest_phase_emm_", run_scope, ".svg"))
  if (nrow(emm_df) > 0 && length(unique(emm_df$Sex)) > 0) {
  forest_plot <- ggplot(
  data = emm_df,
  aes(x = XFacet, y = estimate, ymin = ymin, ymax = ymax, color = pair, shape = pair)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(position = position_dodge(width = 0.5), linewidth = 0.6) +
  geom_text(aes(label = ifelse(pair == "SUS-CON" & p.adjust < 0.05, "★", "")),
    position = position_dodge(width = 0.5), vjust = -1.0, size = 3.8, color = "#333333") +
  scale_color_manual(values = c("RES-CON" = "grey70", "SUS-CON" = "#e63947")) +
  scale_shape_manual(values = c("RES-CON" = 16, "SUS-CON" = 17)) +
  labs(title = paste("EMM contrasts vs CON by Phase:", metric_name), 
     x = if (includeChange) "Cage change" else "allChanges", 
     y = "Estimate (difference vs CON)") +
  facet_grid(Phase ~ Sex, scales = "free_x") +
  theme_pub_small
  ggsave(file_panelA, forest_plot, width = 8.5, height = 4.8)
  }
  
  # Panel B
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
  
  file_panelB <- file.path(dirs$plots_pub, paste0("ribbons_timecourse_", run_scope, ".svg"))
  if (nrow(tc_df) > 0) {
  max_rel <- max(tc_df$rel_time_h, na.rm = TRUE)
  ribbons_plot <- ggplot(tc_df, aes(x = rel_time_h, y = mean_act, color = Group, fill = Group)) +
  geom_ribbon(aes(ymin = mean_act - sd_act, ymax = mean_act + sd_act), alpha = 0.22, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_x_continuous(breaks = scales::breaks_pretty(), limits = c(0, max_rel)) +
  labs(title = paste(metric_name, "over CC-relative time (mean ± SD)"), 
     x = "Time since CC start [h]", 
     y = paste(metric_name, "[a.u.]")) +
  facet_grid(Sex ~ Change, scales = "free_y") +
  theme_pub_small
  ggsave(file_panelB, ribbons_plot, width = 10.5, height = 4.6)
  }
  
  # Panel C
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
  
  file_panelC <- file.path(dirs$plots_pub, paste0("heat_susCon_diff_", run_scope, ".svg"))
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
  ggsave(file.path(dirs$plots_pub, paste0("combined_ABC_", run_scope, ".svg")), 
     combined_emm_abc, width = 10.5, height = 14.2)
  } else if (exists("ribbons_plot")) {
  combined_emm_ab <- (forest_plot / ribbons_plot) +
  patchwork::plot_layout(heights = c(4.8, 4.6), guides = "collect") &
  theme(legend.position = "top")
  ggsave(file.path(dirs$plots_pub, paste0("combined_AB_", run_scope, ".svg")), 
     combined_emm_ab, width = 10.5, height = 9.4)
  }
  }

  # -------------------------------------------------
  # Phase-transition edge windows (last/first 2h)
  # -------------------------------------------------
  edge_plot_dir <- file.path(dirs$plots_pub, "phase_transition_edges")
  edge_tab_dir  <- file.path(dirs$tables_emm, "phase_transition_edges")
  dir_create_safe(edge_plot_dir)
  dir_create_safe(edge_tab_dir)

  edge_window_bins <- 4L  # 4 x 30min = 2h

  # Load base data from file and join ActivityIndex from the aggregated data
  phase_edge_base <- read_csv(data_file, show_col_types = FALSE) %>%
  dplyr::left_join(
    data %>% dplyr::select(AnimalNum, Change, HalfHourElapsed, ActivityIndex),
    by = c("AnimalNum", "Change", "HalfHourElapsed")
  ) %>%
  dplyr::filter(
    !is.na(.data[[metric_name]]),
    !is.na(Phase),
    !is.na(HalfHourElapsed),
    !is.na(AnimalNum),
    !is.na(Group),
    !is.na(Sex),
    !is.na(Change)
  ) %>%
  dplyr::mutate(Phase_chr = as.character(Phase)) %>%
  dplyr::arrange(AnimalNum, Change, HalfHourElapsed) %>%
  dplyr::group_by(AnimalNum, Change) %>%
  dplyr::mutate(row_id = dplyr::row_number()) %>%
  dplyr::ungroup()

  transitions_tbl <- phase_edge_base %>%
  dplyr::group_by(AnimalNum, Change) %>%
  dplyr::mutate(prev_phase = dplyr::lag(Phase_chr)) %>%
  dplyr::filter(!is.na(prev_phase), Phase_chr != prev_phase) %>%
  dplyr::transmute(
    AnimalNum,
    Change,
    trans_row = row_id,
    from_phase = prev_phase,
    to_phase = Phase_chr,
    transition = paste0(from_phase, " -> ", to_phase),
    transition_uid = paste(AnimalNum, Change, trans_row, sep = "_")
  ) %>%
  dplyr::ungroup()

  if (nrow(transitions_tbl) > 0) {
  phase_edge_long <- phase_edge_base %>%
    dplyr::select(AnimalNum, Change, Sex, Group, HalfHourElapsed, row_id, value = all_of(metric_name)) %>%
    dplyr::inner_join(transitions_tbl, by = c("AnimalNum", "Change"), relationship = "many-to-many") %>%
    dplyr::mutate(rel_bin = row_id - trans_row) %>%
    dplyr::filter(rel_bin >= -edge_window_bins, rel_bin <= (edge_window_bins - 1)) %>%
    dplyr::mutate(
    side = dplyr::if_else(rel_bin < 0, "before", "after"),
    side_label = factor(side, levels = c("before", "after"),
              labels = c("Before (last 2h)", "After (first 2h)"))
    )

  phase_edge_event <- phase_edge_long %>%
    dplyr::group_by(AnimalNum, Group, Sex, Change, transition, transition_uid, side, side_label) %>%
    dplyr::summarise(
    edge_mean = mean(value, na.rm = TRUE),
    n_bins = dplyr::n(),
    .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(edge_mean)) %>%
    dplyr::filter(n_bins == edge_window_bins)

  phase_edge_animal <- phase_edge_event %>%
    dplyr::group_by(AnimalNum, Group, Sex, Change, transition, side, side_label) %>%
    dplyr::summarise(
    edge_mean = mean(edge_mean, na.rm = TRUE),
    n_events = dplyr::n(),
    .groups = "drop"
    )

  phase_edge_summary <- phase_edge_animal %>%
    dplyr::group_by(Change, Sex, transition, Group, side, side_label) %>%
    dplyr::summarise(
    n_animals = dplyr::n_distinct(AnimalNum),
    mean_edge = mean(edge_mean, na.rm = TRUE),
    sd_edge = sd(edge_mean, na.rm = TRUE),
    se_edge = sd_edge / sqrt(n_animals),
    .groups = "drop"
    )

  write_pub_xlsx(phase_edge_event,
    file.path(edge_tab_dir, paste0("phaseTransition_eventMeans_", run_scope, ".xlsx")),
    sheet_name = "Results", p_cols = character(0))
  write_pub_xlsx(phase_edge_animal,
    file.path(edge_tab_dir, paste0("phaseTransition_animalMeans_", run_scope, ".xlsx")),
    sheet_name = "Results", p_cols = character(0))
  write_pub_xlsx(phase_edge_summary,
    file.path(edge_tab_dir, paste0("phaseTransition_summary_", run_scope, ".xlsx")),
    sheet_name = "Results", p_cols = character(0))

  p_phase_edge <- ggplot(
    phase_edge_summary,
    aes(x = side_label, y = mean_edge, color = Group, group = Group)
  ) +
    geom_line(position = position_dodge(width = 0.2), linewidth = 0.7) +
    geom_point(position = position_dodge(width = 0.2), size = 2.3) +
    geom_errorbar(
    aes(ymin = mean_edge - se_edge, ymax = mean_edge + se_edge),
    width = 0.12,
    position = position_dodge(width = 0.2),
    linewidth = 0.5
    ) +
    facet_grid(transition ~ Sex + Change, scales = "free_y") +
    scale_color_manual(values = c("CON"="#3e3c6f", "RES"="#c7c3bb", "SUS"="#e63a47")) +
    labs(
    title = paste("Phase-transition edge averages:", metric_name),
    x = NULL,
    y = paste(metric_name, "(mean ± SE)"),
    color = "Group"
    ) +
    theme_minimal(base_size = 12) +
    theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 15, hjust = 1),
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5)
    )

  ggsave(
    file.path(edge_plot_dir, paste0("phaseEdge_beforeAfter_", run_scope, ".svg")),
    p_phase_edge, width = 11, height = 6.5
  )

  # -------------------------------------------------
  # Single-timepoint line plot within the ±2h window
  # -------------------------------------------------
  phase_edge_tp_animal <- phase_edge_long %>%
    dplyr::group_by(AnimalNum, Group, Sex, Change, transition, rel_bin) %>%
    dplyr::summarise(
    tp_mean = mean(value, na.rm = TRUE),
    n_obs = dplyr::n(),
    .groups = "drop"
    )

  phase_edge_tp_summary <- phase_edge_tp_animal %>%
    dplyr::group_by(Change, Sex, transition, Group, rel_bin) %>%
    dplyr::summarise(
    n_animals = dplyr::n_distinct(AnimalNum),
    mean_tp = mean(tp_mean, na.rm = TRUE),
    sd_tp = sd(tp_mean, na.rm = TRUE),
    se_tp = sd_tp / sqrt(n_animals),
    .groups = "drop"
    )

  write_pub_xlsx(phase_edge_tp_animal,
    file.path(edge_tab_dir, paste0("phaseEdge_tp_animals_", run_scope, ".xlsx")),
    sheet_name = "Results", p_cols = character(0))
  write_pub_xlsx(phase_edge_tp_summary,
    file.path(edge_tab_dir, paste0("phaseEdge_tp_summary_", run_scope, ".xlsx")),
    sheet_name = "Results", p_cols = character(0))

  p_phase_edge_tp <- ggplot(
    phase_edge_tp_summary,
    aes(x = rel_bin, y = mean_tp, color = Group, group = Group)
  ) +
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "grey45") +
    geom_line(linewidth = 0.75) +
    geom_point(size = 1.8) +
    geom_errorbar(
    aes(ymin = mean_tp - se_tp, ymax = mean_tp + se_tp),
    width = 0.12,
    linewidth = 0.45
    ) +
    facet_grid(transition ~ Sex + Change, scales = "free_y") +
    scale_color_manual(values = c("CON"="#3e3c6f", "RES"="#c7c3bb", "SUS"="#e63a47")) +
    scale_x_continuous(
    breaks = seq(-edge_window_bins, edge_window_bins - 1, by = 1),
    labels = function(x) sprintf("%+.1f", x / 2)
    ) +
    labs(
    title = paste("Phase-transition timepoint profile (±2h):", metric_name),
    x = "Time relative to transition [h] (0 = first bin after transition)",
    y = paste(metric_name, "(mean ± SE)"),
    color = "Group"
    ) +
    theme_minimal(base_size = 12) +
    theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5)
    )

  ggsave(
    file.path(edge_plot_dir, paste0("phaseEdge_tp_line_", run_scope, ".svg")),
    p_phase_edge_tp, width = 11.5, height = 6.8
  )

  phase_edge_delta <- phase_edge_animal %>%
    dplyr::select(AnimalNum, Group, Sex, Change, transition, side, edge_mean) %>%
    tidyr::pivot_wider(names_from = side, values_from = edge_mean) %>%
    dplyr::filter(!is.na(before), !is.na(after)) %>%
    dplyr::mutate(delta_after_minus_before = after - before)

  phase_edge_delta_sum <- phase_edge_delta %>%
    dplyr::group_by(Change, Sex, transition, Group) %>%
    dplyr::summarise(
    n_animals = dplyr::n_distinct(AnimalNum),
    mean_delta = mean(delta_after_minus_before, na.rm = TRUE),
    sd_delta = sd(delta_after_minus_before, na.rm = TRUE),
    se_delta = sd_delta / sqrt(n_animals),
    .groups = "drop"
    )

  write_pub_xlsx(phase_edge_delta,
    file.path(edge_tab_dir, paste0("phaseTransition_deltaByAnimal_", run_scope, ".xlsx")),
    sheet_name = "Results", p_cols = character(0))
  write_pub_xlsx(phase_edge_delta_sum,
    file.path(edge_tab_dir, paste0("phaseEdge_delta_summary_", run_scope, ".xlsx")),
    sheet_name = "Results", p_cols = character(0))

  p_phase_edge_delta <- ggplot(
    phase_edge_delta_sum,
    aes(x = Group, y = mean_delta, ymin = mean_delta - se_delta, ymax = mean_delta + se_delta, color = Group)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_pointrange(linewidth = 0.6) +
    facet_grid(transition ~ Sex + Change, scales = "free_y") +
    scale_color_manual(values = c("CON"="#3e3c6f", "RES"="#c7c3bb", "SUS"="#e63a47")) +
    labs(
    title = paste("Δ across phase transition (After - Before):", metric_name),
    x = "Group",
    y = "Delta (mean ± SE)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
    )

  ggsave(
    file.path(edge_plot_dir, paste0("phaseEdge_delta_", run_scope, ".svg")),
    p_phase_edge_delta, width = 10.5, height = 6.5
  )

  # -------------------------------------------------
  # Stability plot — median across cage changes,
  # individual animal trajectories + group mean line (RAW values)
  # -------------------------------------------------
  stability_plot_dir <- file.path(dirs$plots_pub, "stability")
  dir_create_safe(stability_plot_dir)

  stability_animal <- phase_edge_animal %>%
    dplyr::group_by(AnimalNum, Group, Sex, transition, side, side_label) %>%
    dplyr::summarise(
    median_edge = median(edge_mean, na.rm = TRUE),
    .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(median_edge))

  write_pub_xlsx(stability_animal,
    file.path(edge_tab_dir, paste0("stability_median_", run_scope, ".xlsx")),
    sheet_name = "Results", p_cols = character(0))

  for (sex_val in unique(stability_animal$Sex)) {
    df_stab <- stability_animal %>%
    dplyr::filter(Sex == sex_val) %>%
    dplyr::mutate(Group = factor(Group, levels = c("CON", "RES", "SUS")))

    p_stability <- ggplot(
    df_stab,
    aes(x = side_label, y = median_edge, group = AnimalNum, color = Group)
    ) +
    geom_line(alpha = 0.3, linewidth = 0.3) +
    geom_point(alpha = 0.4, size = 1.2) +
    stat_summary(aes(group = Group), fun = mean, geom = "line", linewidth = 1.1, position = position_dodge(width = 0.15), na.rm = TRUE) +
    stat_summary(aes(group = Group), fun = mean, geom = "point", size = 2.5, position = position_dodge(width = 0.15), na.rm = TRUE) +
    stat_summary(aes(group = Group), fun.data = mean_se, geom = "errorbar", width = 0.08, linewidth = 0.6, position = position_dodge(width = 0.15), na.rm = TRUE) +
    scale_color_manual(values = group_colors) +
    facet_wrap(~ transition, scales = "free_y") +
    labs(
    title = paste0("Phenotypic Stability (", sex_val, "): ", metric_name),
    subtitle = "Median across cage changes (raw values)",
    x = NULL,
    y = paste(metric_name, "raw median"),
    color = "Group"
    ) +
    theme_minimal(base_size = 13) +
    theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 15, hjust = 1),
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
    )

    ggsave(
    file.path(stability_plot_dir, paste0("stability_median_", sex_val, "_", run_scope, ".svg")),
    p_stability, width = 8, height = 5
    )
  }

  df_stab_all <- stability_animal %>%
    dplyr::mutate(Group = factor(Group, levels = c("CON", "RES", "SUS")))

  p_stability_all <- ggplot(
    df_stab_all,
    aes(x = side_label, y = median_edge, group = AnimalNum, color = Group)
  ) +
    geom_line(alpha = 0.25, linewidth = 0.3) +
    geom_point(alpha = 0.35, size = 1.0) +
    stat_summary(aes(group = Group), fun = mean, geom = "line", linewidth = 1.1, position = position_dodge(width = 0.15), na.rm = TRUE) +
    stat_summary(aes(group = Group), fun = mean, geom = "point", size = 2.5, position = position_dodge(width = 0.15), na.rm = TRUE) +
    stat_summary(aes(group = Group), fun.data = mean_se, geom = "errorbar", width = 0.08, linewidth = 0.6, position = position_dodge(width = 0.15), na.rm = TRUE) +
    scale_color_manual(values = group_colors) +
    facet_grid(transition ~ Sex, scales = "free_y") +
    labs(
    title = paste0("Phenotypic Stability: ", metric_name),
    subtitle = "Median across cage changes (raw values)",
    x = NULL,
    y = paste(metric_name, "raw median"),
    color = "Group"
    ) +
    theme_minimal(base_size = 13) +
    theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 15, hjust = 1),
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
    )

  ggsave(
    file.path(stability_plot_dir, paste0("stability_median_allSex_", run_scope, ".svg")),
    p_stability_all, width = 9, height = 6
  )

  # -------------------------------------------------
  # Stability single-timepoint profile (±2h bins)
  # RAW timepoints (no median collapse across changes/events)
  # -------------------------------------------------
  stability_tp_raw <- phase_edge_long %>%
    dplyr::select(AnimalNum, Group, Sex, Change, transition, transition_uid, rel_bin, value) %>%
    dplyr::filter(!is.na(value))

  write_pub_xlsx(stability_tp_raw,
    file.path(edge_tab_dir, paste0("stability_tp_", run_scope, ".xlsx")),
    sheet_name = "Results", p_cols = character(0))

  for (sex_val in unique(stability_tp_raw$Sex)) {
    df_tp <- stability_tp_raw %>%
    dplyr::filter(Sex == sex_val) %>%
    dplyr::mutate(Group = factor(Group, levels = c("CON", "RES", "SUS")))

    p_stability_tp <- ggplot(
    df_tp,
    aes(x = rel_bin, y = value, group = interaction(AnimalNum, transition_uid), color = Group)
    ) +
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "grey45") +
    geom_line(alpha = 0.10, linewidth = 0.20) +
    stat_summary(aes(group = Group), fun = mean, geom = "line", linewidth = 1.0, na.rm = TRUE) +
    stat_summary(aes(group = Group), fun = mean, geom = "point", size = 2.1, na.rm = TRUE) +
    stat_summary(aes(group = Group), fun.data = mean_se, geom = "errorbar", width = 0.10, linewidth = 0.55, na.rm = TRUE) +
    scale_color_manual(values = group_colors) +
    scale_x_continuous(
    breaks = seq(-edge_window_bins, edge_window_bins - 1, by = 1),
    labels = function(x) sprintf("%+.1f", x / 2)
    ) +
    facet_wrap(~ transition, scales = "free_y") +
    labs(
    title = paste0("Phenotypic Stability timepoints (", sex_val, "): ", metric_name),
    subtitle = "Raw timepoints (no median collapse)",
    x = "Time relative to transition [h]",
    y = paste(metric_name, "raw value"),
    color = "Group"
    ) +
    theme_minimal(base_size = 7) +
    theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 7),
    plot.subtitle = element_text(hjust = 0.5, size = 6),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    strip.text = element_text(size = 7)
    )

    ggsave(
    file.path(stability_plot_dir, paste0("stability_tp_", sex_val, "_", run_scope, ".svg")),
    p_stability_tp, width = 3, height = 3 , units = "in", dpi = 300
    )
  }

  p_stability_tp_all <- ggplot(
    stability_tp_raw %>% dplyr::mutate(Group = factor(Group, levels = c("CON", "RES", "SUS"))),
    aes(x = rel_bin, y = value, group = interaction(AnimalNum, transition_uid), color = Group)
  ) +
    geom_vline(xintercept = -0.5, linetype = "dashed", color = "grey45") +
    geom_line(alpha = 0.08, linewidth = 0.20) +
    stat_summary(aes(group = Group), fun = mean, geom = "line", linewidth = 1.0, na.rm = TRUE) +
    stat_summary(aes(group = Group), fun = mean, geom = "point", size = 2.1, na.rm = TRUE) +
    stat_summary(aes(group = Group), fun.data = mean_se, geom = "errorbar", width = 0.10, linewidth = 0.55, na.rm = TRUE) +
    scale_color_manual(values = group_colors) +
    scale_x_continuous(
    breaks = seq(-edge_window_bins, edge_window_bins - 1, by = 1),
    labels = function(x) sprintf("%+.1f", x / 2)
    ) +
    facet_grid(transition ~ Sex, scales = "free_y") +
    labs(
    title = paste0("Phenotypic Stability timepoints: ", metric_name),
    subtitle = "Raw timepoints (no median collapse)",
    x = "Time relative to transition [h]",
    y = paste(metric_name, "raw value"),
    color = "Group"
    ) +
    theme_minimal(base_size = 13) +
    theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
    )

  ggsave(
    file.path(stability_plot_dir, paste0("stability_tp_allSex_", run_scope, ".svg")),
    p_stability_tp_all, width = 4, height = 4
  )

  } else {
  cat(sprintf("No phase transitions found for %s; phase-transition edge plots skipped.\n", metric_name))
  }

  # -------------------------------------------------
  # For ActivityIndex: create stability plots from raw phase data
  # -------------------------------------------------
  if (metric_name == "ActivityIndex" && nrow(transitions_tbl) == 0) {
  cat(sprintf("Generating ActivityIndex stability plots from phase data (no transitions found).\n"))
  
  stability_plot_dir <- file.path(dirs$plots_pub, "stability")
  dir_create_safe(stability_plot_dir)
  
  # Create artificial transitions for stability analysis (Active->Inactive and Inactive->Active)
  phase_edge_animal <- phase_edge_base %>%
    dplyr::group_by(AnimalNum, Group, Sex) %>%
    dplyr::mutate(
      prev_phase = dplyr::lag(Phase_chr),
      is_transition = !is.na(prev_phase) & Phase_chr != prev_phase,
      transition = dplyr::if_else(is_transition, 
                                  paste0(prev_phase, " -> ", Phase_chr), 
                                  NA_character_)
    ) %>%
    dplyr::filter(!is.na(transition)) %>%
    dplyr::mutate(
      side = dplyr::if_else(is_transition, "after", "after"),  # Mark as after for activity profile
      side_label = factor(side, levels = c("before", "after"),
                         labels = c("Before", "After"))
    ) %>%
    dplyr::group_by(AnimalNum, Group, Sex, transition, side, side_label) %>%
    dplyr::summarise(
      edge_mean = mean(.data[[metric_name]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(AnimalNum, Group, Sex, transition, side, side_label) %>%
    dplyr::summarise(
      edge_mean = mean(edge_mean, na.rm = TRUE),
      .groups = "drop"
    )

  # Generate stability plots for ActivityIndex

  for (sex_val in unique(phase_edge_animal$Sex)) {
    df_stab <- phase_edge_animal %>%
    dplyr::filter(Sex == sex_val) %>%
    dplyr::mutate(Group = factor(Group, levels = c("CON", "RES", "SUS")))

    p_stability <- ggplot(
    df_stab,
    aes(x = side_label, y = edge_mean, group = AnimalNum, color = Group)
    ) +
    geom_line(alpha = 0.3, linewidth = 0.3) +
    geom_point(alpha = 0.4, size = 1.2) +
    stat_summary(aes(group = Group), fun = mean, geom = "line", linewidth = 1.1, position = position_dodge(width = 0.15)) +
    stat_summary(aes(group = Group), fun = mean, geom = "point", size = 2.5, position = position_dodge(width = 0.15)) +
    stat_summary(aes(group = Group), fun.data = mean_se, geom = "errorbar", width = 0.08, linewidth = 0.6, position = position_dodge(width = 0.15)) +
    scale_color_manual(values = group_colors) +
    facet_wrap(~ transition, scales = "free_y") +
    labs(
    title = paste0("Phenotypic Stability (", sex_val, "): ", metric_name),
    subtitle = "Median across phase changes (raw values)",
    x = NULL,
    y = paste(metric_name, "raw value"),
    color = "Group"
    ) +
    theme_minimal(base_size = 13) +
    theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 15, hjust = 1),
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
    )

    ggsave(
    file.path(stability_plot_dir, paste0("stability_phaseChange_raw_", sex_val, "_", run_scope, ".svg")),
    p_stability, width = 8, height = 5
    )
  }

  df_stab_all <- phase_edge_animal %>%
    dplyr::mutate(Group = factor(Group, levels = c("CON", "RES", "SUS")))

  p_stability_all <- ggplot(
    df_stab_all,
    aes(x = side_label, y = edge_mean, group = AnimalNum, color = Group)
  ) +
    geom_line(alpha = 0.25, linewidth = 0.3) +
    geom_point(alpha = 0.35, size = 1.0) +
    stat_summary(aes(group = Group), fun = mean, geom = "line", linewidth = 1.1, position = position_dodge(width = 0.15)) +
    stat_summary(aes(group = Group), fun = mean, geom = "point", size = 2.5, position = position_dodge(width = 0.15)) +
    stat_summary(aes(group = Group), fun.data = mean_se, geom = "errorbar", width = 0.08, linewidth = 0.6, position = position_dodge(width = 0.15)) +
    scale_color_manual(values = group_colors) +
    facet_grid(transition ~ Sex, scales = "free_y") +
    labs(
    title = paste0("Phenotypic Stability: ", metric_name),
    subtitle = "Median across phase changes (raw values)",
    x = NULL,
    y = paste(metric_name, "raw value"),
    color = "Group"
    ) +
    theme_minimal(base_size = 13) +
    theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 15, hjust = 1),
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
    )

  ggsave(
    file.path(stability_plot_dir, paste0("stability_phaseChange_raw_allSexes_", run_scope, ".svg")),
    p_stability_all, width = 9, height = 6
  )
  
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

# Prepare data with all three metrics
corr_data <- data_filtered_agg %>%
  dplyr::select(Change, Sex, Phase, Group, AnimalNum, HalfHourElapsed, TH, Movement, Proximity, ActivityIndex) %>%
  dplyr::filter(!is.na(Movement), !is.na(Proximity), !is.na(ActivityIndex)) %>%
  dplyr::mutate(
    Change = factor(Change),
    Sex = factor(Sex),
    Phase = factor(Phase, levels = c("Active", "Inactive")),
    Group = factor(Group, levels = c("CON", "RES", "SUS"))
  )

cat(sprintf("Correlation data: %d observations with all three metrics\n", nrow(corr_data)))

# Overall correlations - all metric pairs
overall_cor_mp <- cor.test(corr_data$Movement, corr_data$Proximity, method = "pearson")
overall_spearman_mp <- cor.test(corr_data$Movement, corr_data$Proximity, method = "spearman", exact = FALSE)
overall_cor_ma <- cor.test(corr_data$Movement, corr_data$ActivityIndex, method = "pearson")
overall_spearman_ma <- cor.test(corr_data$Movement, corr_data$ActivityIndex, method = "spearman", exact = FALSE)
overall_cor_pa <- cor.test(corr_data$Proximity, corr_data$ActivityIndex, method = "pearson")
overall_spearman_pa <- cor.test(corr_data$Proximity, corr_data$ActivityIndex, method = "spearman", exact = FALSE)

cat(sprintf("Movement vs Proximity - Pearson: r = %.3f, p = %.3e\n", 
            overall_cor_mp$estimate, overall_cor_mp$p.value))
cat(sprintf("Movement vs Proximity - Spearman: rho = %.3f, p = %.3e\n", 
            overall_spearman_mp$estimate, overall_spearman_mp$p.value))
cat(sprintf("Movement vs ActivityIndex - Pearson: r = %.3f, p = %.3e\n", 
            overall_cor_ma$estimate, overall_cor_ma$p.value))
cat(sprintf("Movement vs ActivityIndex - Spearman: rho = %.3f, p = %.3e\n", 
            overall_spearman_ma$estimate, overall_spearman_ma$p.value))
cat(sprintf("Proximity vs ActivityIndex - Pearson: r = %.3f, p = %.3e\n", 
            overall_cor_pa$estimate, overall_cor_pa$p.value))
cat(sprintf("Proximity vs ActivityIndex - Spearman: rho = %.3f, p = %.3e\n", 
            overall_spearman_pa$estimate, overall_spearman_pa$p.value))

# Correlation by Group, Phase, Sex, Change - for all metric pairs
# MP = Movement vs Proximity; MA = Movement vs ActivityIndex; PA = Proximity vs ActivityIndex
corr_by_group_mp <- corr_data %>%
  dplyr::group_by(Change, Sex, Phase, Group) %>%
  dplyr::summarise(
    n = n(),
    metric_pair = "Movement_vs_Proximity",
    pearson_r = cor(Movement, Proximity, method = "pearson", use = "complete.obs"),
    spearman_rho = cor(Movement, Proximity, method = "spearman", use = "complete.obs"),
    .groups = "drop"
  )

corr_by_group_ma <- corr_data %>%
  dplyr::group_by(Change, Sex, Phase, Group) %>%
  dplyr::summarise(
    n = n(),
    metric_pair = "Movement_vs_ActivityIndex",
    pearson_r = cor(Movement, ActivityIndex, method = "pearson", use = "complete.obs"),
    spearman_rho = cor(Movement, ActivityIndex, method = "spearman", use = "complete.obs"),
    .groups = "drop"
  )

corr_by_group_pa <- corr_data %>%
  dplyr::group_by(Change, Sex, Phase, Group) %>%
  dplyr::summarise(
    n = n(),
    metric_pair = "Proximity_vs_ActivityIndex",
    pearson_r = cor(Proximity, ActivityIndex, method = "pearson", use = "complete.obs"),
    spearman_rho = cor(Proximity, ActivityIndex, method = "spearman", use = "complete.obs"),
    .groups = "drop"
  )

# Combine all correlation data
corr_by_group <- dplyr::bind_rows(corr_by_group_mp, corr_by_group_ma, corr_by_group_pa)

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
    metric_pair <- corr_by_group$metric_pair[i]
    
    # Determine which metrics to correlate
    if (metric_pair == "Movement_vs_Proximity") {
      var1 <- d$Movement
      var2 <- d$Proximity
    } else if (metric_pair == "Movement_vs_ActivityIndex") {
      var1 <- d$Movement
      var2 <- d$ActivityIndex
    } else {
      var1 <- d$Proximity
      var2 <- d$ActivityIndex
    }
    
    # Pearson test
    pearson_test <- tryCatch(
      cor.test(var1, var2, method = "pearson"),
      error = function(e) NULL
    )
    if (!is.null(pearson_test)) {
      corr_by_group$pearson_p[i] <- pearson_test$p.value
    }
    
    # Spearman test (with exact = FALSE to avoid ties warning)
    spearman_test <- tryCatch(
      suppressWarnings(
        cor.test(var1, var2, method = "spearman", exact = FALSE)
      ),
      error = function(e) NULL
    )
    if (!is.null(spearman_test)) {
      corr_by_group$spearman_p[i] <- spearman_test$p.value
    }
  }
}

# Save correlation statistics
write_pub_xlsx(corr_by_group,
  file.path(dirs_corr$tables, "correlation_by_group_phase_sex_change.xlsx"),
  sheet_name = "Results", p_cols = character(0))

# Overall scatter plots and by-group / phase-sex plots for each metric pair
corr_pairs <- list(
  list(x = "Movement",  y = "Proximity",    xlab = "Movement [a.u.]",    ylab = "Proximity [a.u.]",    tag = "MP", cor_obj = overall_cor_mp),
  list(x = "Movement",  y = "ActivityIndex", xlab = "Movement [a.u.]",   ylab = "Activity Index [a.u.]", tag = "MA", cor_obj = overall_cor_ma),
  list(x = "Proximity", y = "ActivityIndex", xlab = "Proximity [a.u.]",  ylab = "Activity Index [a.u.]", tag = "PA", cor_obj = overall_cor_pa)
)

for (cp in corr_pairs) {
  x_var <- cp$x; y_var <- cp$y
  xlab  <- cp$xlab; ylab <- cp$ylab
  tag   <- cp$tag;  co   <- cp$cor_obj
  pair_label <- paste(x_var, "vs", y_var)

  # Overall scatter
  p_overall <- ggplot(corr_data, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(aes(color = Group), alpha = 0.3, size = 0.8) +
    geom_smooth(method = "lm", color = "black", linewidth = 1.2) +
    scale_color_manual(values = group_colors) +
    labs(
      title = sprintf("%s\nr = %.3f, p < %.3e", pair_label, co$estimate, co$p.value),
      x = xlab, y = ylab
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5))

  ggsave(
    file.path(dirs_corr$plots, sprintf("correlation_overall_%s.svg", tag)),
    p_overall, width = 7, height = 6
  )

  # By Group
  p_by_group <- ggplot(corr_data, aes(x = .data[[x_var]], y = .data[[y_var]], color = Group)) +
    geom_point(alpha = 0.4, size = 1) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
    scale_color_manual(values = group_colors) +
    facet_wrap(~ Group, ncol = 3) +
    labs(title = sprintf("%s by Group", pair_label), x = xlab, y = ylab) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5))

  ggsave(
    file.path(dirs_corr$plots, sprintf("correlation_by_group_%s.svg", tag)),
    p_by_group, width = 10, height = 4
  )

  # By Phase × Sex
  p_phase_sex <- ggplot(corr_data, aes(x = .data[[x_var]], y = .data[[y_var]], color = Group)) +
    geom_point(alpha = 0.3, size = 0.8) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
    scale_color_manual(values = group_colors) +
    facet_grid(Phase ~ Sex, scales = "free") +
    labs(title = sprintf("%s by Phase and Sex", pair_label), x = xlab, y = ylab) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5))

  ggsave(
    file.path(dirs_corr$plots, sprintf("correlation_by_phase_sex_%s.svg", tag)),
    p_phase_sex, width = 8, height = 7
  )
}


# Correlation heatmap by Group/Phase/Sex — faceted by metric pair and Change
corr_heatmap_data <- corr_by_group %>%
  dplyr::mutate(
    label = sprintf("r=%.2f\np=%.3f", pearson_r, pearson_p),
    sig = ifelse(pearson_p < 0.05, "*", ""),
    metric_pair = factor(metric_pair, levels = c(
      "Movement_vs_Proximity", "Movement_vs_ActivityIndex", "Proximity_vs_ActivityIndex"
    ))
  )

p_heatmap <- ggplot(
  corr_heatmap_data,
  aes(x = interaction(Phase, Sex), y = Group, fill = pearson_r)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f%s", pearson_r, sig)), size = 3, color = "black") +
  scale_fill_gradient2(
    low = "#2c7bb6", 
    mid = "white", 
    high = "#d7191c", 
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson r"
  ) +
  facet_grid(metric_pair ~ Change) +
  labs(
    title = "Correlation Heatmap: All Metric Pairs",
    subtitle = "* indicates p < 0.05",
    x = "Phase × Sex",
    y = "Group"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(size = 9)
  )

ggsave(
  file.path(dirs_corr$plots, "correlation_heatmap.svg"),
  p_heatmap, width = 12, height = 10
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
    scale_color_manual(values = group_colors) +
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
  write_pub_xlsx(time_corr,
    file.path(dirs_corr$tables, "correlation_time_resolved.xlsx"),
    sheet_name = "Results", p_cols = character(0))
  
  cat(sprintf("  - Time-resolved correlations: %d valid time points (out of %d total)\n", 
              nrow(time_corr), 
              nrow(corr_data %>% dplyr::group_by(Change, Sex, Phase, Group, HalfHourElapsed) %>% dplyr::summarise(n = n(), .groups = "drop"))))
}

# Compute animal-level correlations
animal_corr_raw <- corr_data %>%
  dplyr::group_by(AnimalNum, Group, Sex, Change) %>%
  dplyr::summarise(
    pearson_r = if(n() >= 10) cor(Movement, Proximity, use = "complete.obs") else NA_real_,
    .groups = "drop"
  )

# Average correlations across changes for each animal
animal_corr <- animal_corr_raw %>%
  dplyr::group_by(AnimalNum, Group, Sex) %>%
  dplyr::summarise(
    pearson_r = mean(pearson_r, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::filter(!is.na(pearson_r))

# Save animal-level correlations
write_pub_xlsx(animal_corr,
  file.path(dirs_corr$tables, "correlation_by_animal.xlsx"),
  sheet_name = "Results", p_cols = character(0))

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
            comparison = paste(group2, group1, sep = " - "),
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
write_pub_xlsx(overall_stats,
  file.path(dirs_corr$tables, "animal_correlation_overall_stats.xlsx"),
  sheet_name = "Results", p_cols = character(0))

write_pub_xlsx(posthoc_stats,
  file.path(dirs_corr$tables, "animal_correlation_posthoc_stats.xlsx"),
  sheet_name = "Results", p_cols = character(0))

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

write_pub_xlsx(animal_corr_summary,
  file.path(dirs_corr$tables, "animal_correlation_summary.xlsx"),
  sheet_name = "Results", p_cols = character(0))

# Enhanced violin plot with post-hoc comparisons
p_animal_corr <- ggplot(animal_corr, aes(x = Group, y = pearson_r, fill = Group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 2) +
  scale_fill_manual(values = group_colors) +
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

# Add post-hoc brackets
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
  scale_fill_manual(values = group_colors) +
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

cat("\n ✓ Animal-level correlation statistics with facet-specific annotations saved\n")
cat(sprintf("  - Overall stats saved to: %s\n", file.path(dirs_corr$tables, "animal_correlation_overall_stats.xlsx")))
cat(sprintf("  - Post-hoc stats saved to: %s\n", file.path(dirs_corr$tables, "animal_correlation_posthoc_stats.xlsx")))
cat(sprintf("  - Summary stats saved to: %s\n", file.path(dirs_corr$tables, "animal_correlation_summary.xlsx")))
cat("\n === CORRELATION ANALYSIS COMPLETE ===\n")
cat(sprintf("Results saved to: %s\n", correlation_dir))

# -------------------------------------------------
# AR1 sensitivity analysis (separate output branch)
# -------------------------------------------------
fit_lme_ar1_adaptive <- function(fixed_rhs, df, response_var, time_var = "TH_scaled", animal = "AnimalNum", allow_slope = TRUE) {
  f_fixed <- as.formula(paste0(response_var, " ~ ", fixed_rhs))
  rand_slope <- as.formula(paste0("~ ", time_var, " | ", animal))
  rand_int <- as.formula(paste0("~ 1 | ", animal))
  cor_form <- as.formula(paste0("~ ", time_var, " | ", animal))

  ctrl <- nlme::lmeControl(
    opt = "optim",
    msMaxIter = 200,
    msMaxEval = 400,
    niterEM = 40,
    returnObject = TRUE
  )

  fit_try <- function(random_form, use_ar1 = TRUE) {
    tryCatch(
      nlme::lme(
        fixed = f_fixed,
        random = random_form,
        correlation = if (use_ar1) nlme::corAR1(form = cor_form) else NULL,
        data = df,
        control = ctrl,
        na.action = na.omit
      ),
      error = function(e) NULL
    )
  }

  # Fallback sequence: slope+AR1 -> intercept+AR1 -> slope(no AR1) -> intercept(no AR1)
  if (allow_slope) {
    m1 <- fit_try(rand_slope, use_ar1 = TRUE)
    if (!is.null(m1)) return(m1)
  }

  m2 <- fit_try(rand_int, use_ar1 = TRUE)
  if (!is.null(m2)) return(m2)

  if (allow_slope) {
    m3 <- fit_try(rand_slope, use_ar1 = FALSE)
    if (!is.null(m3)) return(m3)
  }

  m4 <- fit_try(rand_int, use_ar1 = FALSE)
  m4
}

compute_phase_avg_for_model_ar1 <- function(model, Change, Sex, Phase) {
  emm <- emmeans::emmeans(model, specs = ~ Group)
  contrast_method <- setNames(contrast_config$vectors, contrast_config$display_names)
  contr_res <- emmeans::contrast(emm, method = contrast_method, adjust = "none")

  safe_format <- function(emm_object) {
    df <- as.data.frame(summary(emm_object, infer = TRUE))
    low_col  <- grep("lower|LCL", names(df), value = TRUE)
    upp_col  <- grep("upper|UCL", names(df), value = TRUE)
    stat_col <- grep("t.ratio|z.ratio", names(df), value = TRUE)

    df <- df %>%
      dplyr::mutate(
        Change = as.character(Change),
        Sex = as.character(Sex),
        Phase = as.character(Phase)
      )

    if (length(low_col) > 0) df <- df %>% dplyr::rename(lwr = !!low_col[1])
    if (length(upp_col) > 0) df <- df %>% dplyr::rename(upr = !!upp_col[1])
    if (length(stat_col) > 0) df <- df %>% dplyr::rename(t.ratio = !!stat_col[1])

    df
  }

  avg_df <- safe_format(emm) %>% dplyr::rename(emmean_avg = emmean, SE_avg = SE)

  # Preferred effect size: emmeans::eff_size matched by contrast label
  # (fallback to estimate/sigma if eff_size fails)
  sigma_val <- tryCatch(as.numeric(stats::sigma(model)), error = function(e) NA_real_)
  edf_val <- tryCatch(as.numeric(stats::df.residual(model)), error = function(e) NA_real_)
  if (is.na(sigma_val) || sigma_val <= 0) sigma_val <- NA_real_

  eff_df <- tryCatch({
    eff_obj <- emmeans::eff_size(
      emm,
      sigma = sigma_val,
      edf = edf_val,
      method = contrast_method
    )
    as.data.frame(eff_obj)
  }, error = function(e) NULL)

  contr_df <- safe_format(contr_res)

  if (!is.null(eff_df) && nrow(eff_df) > 0 && "contrast" %in% names(eff_df)) {
    idx <- match(as.character(contr_df$contrast), as.character(eff_df$contrast))
    contr_df$d <- if ("effect.size" %in% names(eff_df)) as.numeric(eff_df$effect.size[idx]) else NA_real_
    contr_df$SE_d <- if ("SE" %in% names(eff_df)) as.numeric(eff_df$SE[idx]) else NA_real_
    if ("lower.CL" %in% names(eff_df)) {
      contr_df$lwr_d <- as.numeric(eff_df$lower.CL[idx])
    } else if ("asymp.LCL" %in% names(eff_df)) {
      contr_df$lwr_d <- as.numeric(eff_df$asymp.LCL[idx])
    } else {
      contr_df$lwr_d <- NA_real_
    }
    if ("upper.CL" %in% names(eff_df)) {
      contr_df$upr_d <- as.numeric(eff_df$upper.CL[idx])
    } else if ("asymp.UCL" %in% names(eff_df)) {
      contr_df$upr_d <- as.numeric(eff_df$asymp.UCL[idx])
    } else {
      contr_df$upr_d <- NA_real_
    }
  } else {
    # Fallback: standardized contrasts using residual sigma
    contr_df <- contr_df %>%
      dplyr::mutate(
        d = ifelse(is.na(sigma_val), NA_real_, estimate / sigma_val),
        SE_d = ifelse(is.na(sigma_val), NA_real_, SE / sigma_val),
        lwr_d = ifelse(is.na(sigma_val), NA_real_, lwr / sigma_val),
        upr_d = ifelse(is.na(sigma_val), NA_real_, upr / sigma_val)
      )
  }

  list(means = avg_df, contrasts = contr_df)
}

run_analysis_for_metric_ar1 <- function(metric_name, data, includeChange, includeSex, includePhase) {
  cat(sprintf("\n\n========================================\n"))
  cat(sprintf("AR1 SENSITIVITY ANALYSIS: %s\n", metric_name))
  cat(sprintf("========================================\n\n"))

  results_dir <- file.path(base_results_dir, "AR1_sensitivity", metric_name)
  dir_create_safe(results_dir)

  dirs <- list(
    logs        = file.path(results_dir, "logs"),
    models      = file.path(results_dir, "models"),
    artifacts   = file.path(results_dir, "artifacts"),
    tables_sum  = file.path(results_dir, "tables", "summary"),
    tables_emm  = file.path(results_dir, "tables", "emmeans"),
    tables_phase = file.path(results_dir, "tables", "phase_average"),
    plots_pub   = file.path(results_dir, "plots", "publication")
  )
  invisible(lapply(dirs, dir_create_safe))

  changes <- if (includeChange) unique(data$Change) else "allChanges"
  sexes   <- if (includeSex)    unique(data$Sex)    else "allSexes"
  phases  <- if (includePhase)  unique(data$Phase)  else "allPhases"

  log_file <- file.path(dirs$logs, paste0("fit_log_ar1_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

  fit_model_slices_ar1 <- function(window_label = "all") {
    summary_results_i <- tibble::tibble()
    emmeans_results_i <- tibble::tibble()
    phase_avg_means_i <- list()
    phase_avg_contr_i <- list()

    for (Change in changes) {
      for (Sex in sexes) {
        for (Phase in phases) {
          dsub <- data %>%
            dplyr::filter(
              (!includeChange | Change == !!Change),
              (!includeSex    | Sex    == !!Sex),
              (!includePhase  | Phase  == !!Phase)
            )
          if (nrow(dsub) == 0) next

          if (window_label %in% c("first2h", "last2h")) {
            t0_phase <- min(dsub$HalfHourElapsed, na.rm = TRUE)
            t1_phase <- max(dsub$HalfHourElapsed, na.rm = TRUE)

            if (window_label == "first2h") {
              dsub <- dsub %>%
                dplyr::mutate(rel_phase_h = (HalfHourElapsed - t0_phase) / 2) %>%
                dplyr::filter(!is.na(rel_phase_h), rel_phase_h >= 0, rel_phase_h <= 2)
            }

            if (window_label == "last2h") {
              dsub <- dsub %>%
                dplyr::mutate(rel_phase_h_from_end = (t1_phase - HalfHourElapsed) / 2) %>%
                dplyr::filter(!is.na(rel_phase_h_from_end), rel_phase_h_from_end >= 0, rel_phase_h_from_end <= 2)
            }
          }

          n_anim  <- dplyr::n_distinct(dsub$AnimalNum)
          n_group <- dplyr::n_distinct(dsub$Group)
          if (n_anim < min_animals || n_group < min_groups || nrow(dsub) < min_obs) next

          t0 <- min(dsub$HalfHourElapsed, na.rm = TRUE)
          dsub <- dsub %>%
            dplyr::mutate(
              TH = (HalfHourElapsed - t0) / 2,
              TH_scaled = as.vector(scale(TH)),
              AnimalNum = factor(AnimalNum)
            ) %>%
            dplyr::filter(!is.na(TH_scaled), !is.na(.data[[metric_name]]))

          if (nrow(dsub) == 0) next

          use_int <- use_group_time_interaction && (dplyr::n_distinct(dsub$TH_scaled) >= min_time_lv)
          fixed_rhs <- if (use_int) "Group * TH_scaled" else "Group + TH_scaled"

          per_animal <- dsub %>% dplyr::count(AnimalNum, name = "n_i")
          allow_slope <- all(per_animal$n_i >= 3) && (dplyr::n_distinct(dsub$TH_scaled) >= min_time_lv)

          model <- fit_lme_ar1_adaptive(
            fixed_rhs = fixed_rhs,
            df = dsub,
            response_var = metric_name,
            time_var = "TH_scaled",
            animal = "AnimalNum",
            allow_slope = allow_slope
          )

          if (is.null(model)) {
            cat(paste(Sys.time(), Change, Sex, Phase, window_label, "AR1 model failed\n"), file = log_file, append = TRUE)
            next
          }

          model_filename <- paste0(
            "model_ar1_", metric_name, "_",
            "Change-", Change, "_",
            "Sex-", Sex, "_",
            "Phase-", Phase, "_",
            "Window-", window_label,
            ".rds"
          )
          saveRDS(model, file.path(dirs$models, model_filename))

          ms <- summary(model)
          tt <- as.data.frame(ms$tTable)
          tt$Fixed_effect <- rownames(tt)
          p_col <- grep("p-value|Pr\\(>", names(tt), value = TRUE)
          if (length(p_col) == 0) tt$p.value <- NA_real_ else tt$p.value <- as.numeric(tt[[p_col[1]]])
          phi_ar1 <- tryCatch(as.numeric(nlme::coef(model$modelStruct$corStruct, unconstrained = FALSE)), error = function(e) NA_real_)

          fx_tbl <- tt %>%
            dplyr::mutate(
              Estimate = as.numeric(Value),
              Std.Error = as.numeric(`Std.Error`),
              df_Satterthwaite = as.numeric(DF),
              t.value = as.numeric(`t-value`),
              p.adj_FDR = p.adjust(p.value, method = "fdr"),
              AIC = tryCatch(AIC(model), error = function(e) NA_real_),
              BIC = tryCatch(BIC(model), error = function(e) NA_real_),
              logLik = tryCatch(as.numeric(logLik(model)), error = function(e) NA_real_),
              ar1_phi = phi_ar1,
              model_formula = paste(deparse(formula(model)), collapse = " "),
              model_engine = "nlme::lme (AR1 sensitivity)",
              Change = as.character(Change),
              Phase = as.character(Phase),
              Sex = as.character(Sex),
              window = window_label
            ) %>%
            dplyr::select(
              Fixed_effect, Estimate, Std.Error, df_Satterthwaite, t.value,
              p.value, p.adj_FDR, AIC, BIC, logLik, ar1_phi,
              model_formula, model_engine,
              Change, Phase, Sex, window
            )
          summary_results_i <- dplyr::bind_rows(summary_results_i, fx_tbl)

          if (!any(grepl("^Group", names(nlme::fixef(model))))) next

          emm_lvl <- suppressMessages(emmeans::emmeans(model, pairwise ~ Group, at = list(TH_scaled = 0)))
          emm_lvl_df <- as.data.frame(summary(emm_lvl$contrasts, infer = TRUE, adjust = "none")) %>%
            dplyr::mutate(
              contrast_type = "level@intercept",
              Change = as.character(Change),
              Phase = as.character(Phase),
              Sex = as.character(Sex),
              window = window_label,
              model_engine = "nlme::lme (AR1 sensitivity)"
            )
          emm_lvl_df <- emm_lvl_df[!is.na(emm_lvl_df$p.value), ]
          emm_lvl_df$cohens_d <- tryCatch({
            emm_lvl_raw <- emmeans::emmeans(model, ~ Group, at = list(TH_scaled = 0))
            eff_d <- emmeans::eff_size(emm_lvl_raw, sigma = sigma(model), edf = stats::df.residual(model))
            eff_df <- as.data.frame(eff_d)
            idx <- match(emm_lvl_df$contrast, eff_df$contrast)
            eff_df$effect.size[idx]
          }, error = function(e) NA_real_)
          if (nrow(emm_lvl_df) > 0) emmeans_results_i <- dplyr::bind_rows(emmeans_results_i, emm_lvl_df)

          if (use_int && any(grepl("Group:TH_scaled", names(nlme::fixef(model))))) {
            emm_tr <- suppressMessages(emmeans::emtrends(model, pairwise ~ Group, var = "TH_scaled"))
            emm_tr_df <- as.data.frame(summary(emm_tr$contrasts, infer = TRUE, adjust = "none")) %>%
              dplyr::mutate(
                contrast_type = "slope",
                Change = as.character(Change),
                Phase = as.character(Phase),
                Sex = as.character(Sex),
                window = window_label,
                model_engine = "nlme::lme (AR1 sensitivity)",
                cohens_d = NA_real_
              )
            emm_tr_df <- emm_tr_df[!is.na(emm_tr_df$p.value), ]
            if (nrow(emm_tr_df) > 0) emmeans_results_i <- dplyr::bind_rows(emmeans_results_i, emm_tr_df)
          }

          stat_k <- compute_phase_avg_for_model_ar1(model, Change, Sex, Phase)
          if (!is.null(stat_k$means) && nrow(stat_k$means)) {
            stat_k$means$window <- window_label
            stat_k$means$model_engine <- "nlme::lme (AR1 sensitivity)"
            phase_avg_means_i[[length(phase_avg_means_i) + 1]] <- stat_k$means
          }
          if (!is.null(stat_k$contrasts) && nrow(stat_k$contrasts)) {
            stat_k$contrasts$window <- window_label
            stat_k$contrasts$model_engine <- "nlme::lme (AR1 sensitivity)"
            phase_avg_contr_i[[length(phase_avg_contr_i) + 1]] <- stat_k$contrasts
          }
        }
      }
    }

    if (nrow(emmeans_results_i) > 0) {
      emmeans_results_i <- add_dual_adjustment_columns(
        emmeans_results_i,
        family_cols = c("Change", "Sex", "Phase", "window", "contrast_type"),
        primary_method = "BH",
        strict_method = "holm",
        primary_contrast_set = primary_contrast_set
      )
    }

    list(
      summary = summary_results_i,
      emmeans = emmeans_results_i,
      phase_means = if (length(phase_avg_means_i)) dplyr::bind_rows(phase_avg_means_i) else tibble::tibble(),
      phase_contr = if (length(phase_avg_contr_i)) dplyr::bind_rows(phase_avg_contr_i) else tibble::tibble()
    )
  }

  run_scope <- scope_all(includeChange, includeSex, includePhase)

  save_ar1_publication_plots <- function(res_obj, suffix = "") {
    # Volcano plot for phase-average contrasts
    if (!is.null(res_obj$phase_contr) && nrow(res_obj$phase_contr) > 0 &&
        all(c("estimate", "p.adjust", "contrast") %in% names(res_obj$phase_contr))) {
      vdf <- res_obj$phase_contr %>%
        dplyr::filter(!is.na(estimate), !is.na(p.adjust), p.adjust > 0) %>%
        dplyr::mutate(
          contrast = factor(contrast, levels = contrast_config$display_names),
          sig = p.adjust < 0.05,
          large_effect = abs(estimate) > 1,
          label_flag = sig & large_effect,
          Phase = factor(Phase, levels = c("Active", "Inactive")),
          Sex = factor(Sex)
        )

      if (nrow(vdf) > 0) {
        est_range <- range(vdf$estimate, na.rm = TRUE)
        est_min <- floor(est_range[1])
        est_max <- ceiling(est_range[2])
        x_breaks <- scales::pretty_breaks(n = 6)(c(est_min, est_max))

        p_volcano_ar1 <- ggplot(vdf, aes(x = estimate, y = -log10(p.adjust), color = contrast, shape = Phase)) +
          geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "gray50") +
          geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray70") +
          geom_point(size = 5, alpha = 0.8) +
          scale_color_manual(values = contrast_config$colors) +
          scale_shape_manual(values = c("Active" = 16, "Inactive" = 1), name = "Phase") +
          facet_grid(. ~ Sex) +
          scale_x_continuous(breaks = x_breaks, minor_breaks = NULL) +
          labs(
            title = paste("Volcano Plot:", metric_name),
            x = "Effect size (estimate)",
            y = expression(-log[10](adjusted~p~value)),
            color = "Contrast"
          ) +
          theme_minimal(base_size = 14) +
          theme(
            panel.grid.major.x = element_line(linewidth = 0.2),
            panel.grid.major.y = element_line(linewidth = 0.2),
            legend.position = "top"
          )

        ggsave(
          file.path(dirs$plots_pub, paste0("volcano_phaseAvg_ar1_", run_scope, suffix, ".svg")),
          p_volcano_ar1,
          width = 8,
          height = 6
        )
      }
    }

    # Forest plot for phase-average contrasts
    if (!is.null(res_obj$phase_contr) && nrow(res_obj$phase_contr) > 0 &&
        all(c("estimate", "lwr", "upr", "contrast") %in% names(res_obj$phase_contr))) {
      fdf <- res_obj$phase_contr %>%
        dplyr::filter(!is.na(estimate), !is.na(lwr), !is.na(upr)) %>%
        dplyr::mutate(
          contrast = factor(contrast, levels = contrast_config$display_names),
          Phase = factor(Phase, levels = c("Active", "Inactive")),
          Change = factor(Change)
        )

      if (nrow(fdf) > 0) {
        p_forest_ar1 <- ggplot(fdf, aes(y = contrast, x = estimate, xmin = lwr, xmax = upr, color = contrast)) +
          geom_point(size = 3) +
          geom_errorbar(width = 0.2, orientation = "y") +
          scale_color_manual(values = contrast_config$colors) +
          facet_grid(Phase ~ Change) +
          labs(
            title = paste("Forest Plot:", metric_name),
            x = "Estimate (phase avg, with 95% CI)",
            y = "Contrast",
            color = "Contrast"
          ) +
          theme_minimal(base_size = 14)

        ggsave(
          file.path(dirs$plots_pub, paste0("forest_phaseAvg_contrasts_ar1_", run_scope, suffix, ".svg")),
          p_forest_ar1,
          width = 9,
          height = 7
        )

        # Additional split plots: one file per Sex × Phase (2x2 design)
        fdf_split <- res_obj$phase_contr %>%
          dplyr::mutate(
            pair = dplyr::case_when(
              grepl("^RES\\s*-\\s*CON$", as.character(contrast)) ~ "RES-CON",
              grepl("^SUS\\s*-\\s*CON$", as.character(contrast)) ~ "SUS-CON",
              grepl("^SUS\\s*-\\s*RES$", as.character(contrast)) ~ "SUS-RES",
              TRUE ~ NA_character_
            )
          ) %>%
          dplyr::filter(pair %in% c("RES-CON", "SUS-CON", "SUS-RES")) %>%
          dplyr::mutate(
            pair = factor(pair, levels = c("RES-CON", "SUS-CON", "SUS-RES")),
            Change = dplyr::recode(as.character(Change), "CC1" = "1", "CC2" = "2", "CC3" = "3", "CC4" = "4"),
            Change = factor(Change, levels = sort(unique(Change))),
            Sex = factor(Sex),
            Phase = factor(Phase, levels = c("Active", "Inactive"))
          ) %>%
          dplyr::filter(!is.na(estimate), !is.na(lwr), !is.na(upr), !is.na(p.adjust))

        sex_levels <- unique(as.character(fdf_split$Sex))
        phase_levels <- unique(as.character(fdf_split$Phase))

        for (sx in sex_levels) {
          for (ph in phase_levels) {
            sub_sp <- fdf_split %>% dplyr::filter(as.character(Sex) == sx, as.character(Phase) == ph)
            if (nrow(sub_sp) == 0) next

            p_forest_sp <- ggplot(
              sub_sp,
              aes(x = estimate, y = pair, xmin = lwr, xmax = upr, color = pair)
            ) +
              geom_vline(xintercept = 0, linetype = "solid", color = "grey85", linewidth = 0.3) +
              geom_pointrange(size = 0.4, linewidth = 1.0) +
              geom_text(
                aes(
                  x = upr,
                  label = paste0(
                    ifelse(p.adjust < 0.001, "***", ifelse(p.adjust < 0.01, "**", ifelse(p.adjust < 0.05, "*", ""))),
                    if ("d" %in% names(sub_sp)) ifelse(is.finite(d), paste0(" d=", sprintf("%.2f", d)), "") else ""
                  )
                ),
                hjust = -0.1,
                size = 2.2,
                color = "black",
                show.legend = FALSE
              ) +
              facet_grid(. ~ Change, scales = "free_x", labeller = labeller(Change = function(x) rep("", length(x)))) +
              scale_color_manual(values = pair_colors) +
              labs(
                title = paste0(metric_name, ": phase-average contrasts (", ph, ", ", sx, ")"),
                subtitle = "Stats shown are emmeans phase-average contrasts (adjusted p)",
                x = paste0("\u0394 ", metric_name, " (a.u.)"),
                y = NULL
              ) +
              theme_classic(base_size = 8) +
              theme(
                panel.spacing = unit(0.8, "lines"),
                strip.background = element_blank(),
                strip.text = element_text(face = "bold", size = 8),
                axis.line.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.text.y = element_blank(),
                legend.position = "bottom",
                legend.title = element_blank(),
                legend.text = element_text(size = 7),
                legend.key.size = unit(0.3, "cm"),
                plot.margin = margin(5, 45, 5, 5),
                panel.grid.major.x = element_line(color = "grey98"),
                plot.title = element_text(face = "bold", hjust = 0.5, size = 9),
                plot.subtitle = element_text(hjust = 0.5, size = 7)
              ) +
              scale_x_continuous(expand = expansion(mult = c(0.1, 0.6)))

            sx_slug <- gsub("[^A-Za-z0-9]+", "-", tolower(sx))
            ph_slug <- gsub("[^A-Za-z0-9]+", "-", tolower(ph))
            ggsave(
              file.path(dirs$plots_pub, paste0("forest_phaseAvg_ar1_", run_scope, suffix, "_sex-", sx_slug, "_phase-", ph_slug, ".svg")),
              p_forest_sp,
              width = if (isTRUE(includeChange)) 3 * 2.5 else 2.5,
              height = 3,
              device = "svg"
            )
          }
        }
      }
    }

    # EMM means plot (phase-average means)
    if (!is.null(res_obj$phase_means) && nrow(res_obj$phase_means) > 0 &&
        all(c("Change", "Group", "emmean_avg", "SE_avg") %in% names(res_obj$phase_means))) {
      mdf <- res_obj$phase_means %>%
        dplyr::filter(!is.na(emmean_avg), !is.na(SE_avg)) %>%
        dplyr::mutate(
          Group = factor(Group, levels = c("CON", "RES", "SUS")),
          Change_chr = as.character(Change),
          Change = factor(Change_chr, levels = {
            lev <- unique(as.character(Change_chr))
            num <- suppressWarnings(as.numeric(stringr::str_extract(lev, "\\d+")))
            if (all(!is.na(num))) lev[order(num)] else sort(lev)
          }),
          Change_num = as.numeric(Change),
          Phase = factor(Phase, levels = c("Active", "Inactive")),
          Sex = factor(Sex)
        )

      if (nrow(mdf) > 0) {
        if (!("lwr" %in% names(mdf) && "upr" %in% names(mdf))) {
          mdf <- mdf %>% dplyr::mutate(lwr = emmean_avg - 1.96 * SE_avg, upr = emmean_avg + 1.96 * SE_avg)
        }

        p_emm_ar1 <- ggplot(mdf, aes(x = Change_num, y = emmean_avg, color = Group, group = Group)) +
          geom_line(linewidth = 0.7, alpha = 0.95) +
          geom_point(size = 1.8, stroke = 0.2) +
          geom_ribbon(aes(ymin = lwr, ymax = upr, fill = Group), alpha = 0.14, color = NA, inherit.aes = TRUE) +
          scale_color_manual(values = group_colors) +
          scale_fill_manual(values = group_colors, guide = "none") +
          scale_x_continuous(
            breaks = seq_along(levels(mdf$Change)),
            labels = levels(mdf$Change),
            expand = expansion(mult = c(0.02, 0.02))
          ) +
          coord_cartesian(clip = "off") +
          labs(
            title = paste0(metric_name, ": estimated development over cage changes"),
            subtitle = ifelse(suffix == "", "All phase duration", gsub("^_", "", suffix)),
            x = "Cage change",
            y = paste0("Estimated ", metric_name, " (EMM)")
          ) +
          theme_classic(base_size = 8, base_family = "Arial") +
          theme(
            legend.position = "top",
            legend.title = element_blank(),
            legend.key.size = unit(0.35, "cm"),
            legend.text = element_text(size = 7),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 7, color = "black"),
            axis.line = element_line(linewidth = 0.4, color = "black"),
            axis.ticks = element_line(linewidth = 0.4, color = "black"),
            axis.ticks.length = unit(1.5, "mm"),
            plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
            plot.subtitle = element_text(hjust = 0.5, size = 7),
            plot.margin = margin(4, 6, 4, 4)
          )

        if (dplyr::n_distinct(mdf$Phase) > 1 && dplyr::n_distinct(mdf$Sex) > 1) {
          p_emm_ar1 <- p_emm_ar1 + facet_grid(Phase ~ Sex, scales = "free_y")
        } else if (dplyr::n_distinct(mdf$Phase) > 1) {
          p_emm_ar1 <- p_emm_ar1 + facet_wrap(~ Phase, scales = "free_y")
        } else if (dplyr::n_distinct(mdf$Sex) > 1) {
          p_emm_ar1 <- p_emm_ar1 + facet_wrap(~ Sex, scales = "free_y")
        }

        ggsave(
          file.path(dirs$plots_pub, paste0("emm_phaseAvg_means_ar1_", run_scope, suffix, ".svg")),
          p_emm_ar1,
          width = 6,
          height = 4
        )
      }
    }
  }

  save_ar1_grouped_window_phase_forest <- function(main_res, first2h_res = NULL, last2h_res = NULL) {
    tbl_all <- list()
    if (!is.null(main_res$phase_contr) && nrow(main_res$phase_contr) > 0) {
      tbl_all[[length(tbl_all) + 1]] <- main_res$phase_contr %>% dplyr::mutate(window_group = "all")
    }
    if (!is.null(first2h_res) && !is.null(first2h_res$phase_contr) && nrow(first2h_res$phase_contr) > 0) {
      tbl_all[[length(tbl_all) + 1]] <- first2h_res$phase_contr %>% dplyr::mutate(window_group = "first2h")
    }
    if (!is.null(last2h_res) && !is.null(last2h_res$phase_contr) && nrow(last2h_res$phase_contr) > 0) {
      tbl_all[[length(tbl_all) + 1]] <- last2h_res$phase_contr %>% dplyr::mutate(window_group = "last2h")
    }

    if (length(tbl_all) == 0) return(invisible(NULL))

    gdf <- dplyr::bind_rows(tbl_all) %>%
      dplyr::mutate(
        pair = dplyr::case_when(
          grepl("^RES\\s*-\\s*CON$", as.character(contrast)) ~ "RES-CON",
          grepl("^SUS\\s*-\\s*CON$", as.character(contrast)) ~ "SUS-CON",
          grepl("^SUS\\s*-\\s*RES$", as.character(contrast)) ~ "SUS-RES",
          TRUE ~ NA_character_
        )
      ) %>%
      dplyr::filter(pair %in% c("RES-CON", "SUS-CON", "SUS-RES")) %>%
      dplyr::mutate(
        pair = factor(pair, levels = c("RES-CON", "SUS-CON", "SUS-RES")),
        window_group = factor(window_group, levels = c("all", "first2h", "last2h"),
                              labels = c("All", "First 2h", "Last 2h")),
        Phase = factor(Phase, levels = c("Active", "Inactive")),
        Change = dplyr::recode(as.character(Change), "CC1" = "1", "CC2" = "2", "CC3" = "3", "CC4" = "4"),
        Change = factor(Change, levels = sort(unique(Change))),
        Sex = factor(Sex)
      ) %>%
      dplyr::filter(!is.na(estimate), !is.na(lwr), !is.na(upr), !is.na(p.adjust))

    if (nrow(gdf) == 0) return(invisible(NULL))

    for (sx in unique(as.character(gdf$Sex))) {
      sub_sx <- gdf %>% dplyr::filter(as.character(Sex) == sx)
      if (nrow(sub_sx) == 0) next

      p_grouped <- ggplot(
        sub_sx,
        aes(x = estimate, y = pair, xmin = lwr, xmax = upr, color = pair)
      ) +
        geom_vline(xintercept = 0, linetype = "solid", color = "grey85", linewidth = 0.3) +
        geom_pointrange(size = 0.4, linewidth = 1.0) +
        geom_text(
          aes(
            x = upr,
            label = paste0(
              ifelse(p.adjust < 0.001, "***", ifelse(p.adjust < 0.01, "**", ifelse(p.adjust < 0.05, "*", ""))),
              if ("d" %in% names(sub_sx)) ifelse(is.finite(d), paste0(" d=", sprintf("%.2f", d)), "") else ""
            )
          ),
          hjust = -0.1,
          size = 2.1,
          color = "black",
          show.legend = FALSE
        ) +
        facet_grid(Phase ~ window_group + Change, scales = "free_x",
                   labeller = labeller(Change = function(x) rep("", length(x)))) +
        scale_color_manual(values = pair_colors) +
        labs(
          title = paste0(metric_name, ": grouped phase-average contrasts (", sx, ")"),
          subtitle = "Columns: window × cage change | Rows: phase",
          x = paste0("\u0394 ", metric_name, " (a.u.)"),
          y = NULL
        ) +
        theme_classic(base_size = 8) +
        theme(
          panel.spacing = unit(0.75, "lines"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 8),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 7),
          legend.key.size = unit(0.3, "cm"),
          plot.margin = margin(5, 55, 5, 5),
          panel.grid.major.x = element_line(color = "grey98"),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 9),
          plot.subtitle = element_text(hjust = 0.5, size = 7)
        ) +
        scale_x_continuous(expand = expansion(mult = c(0.1, 0.6)))

      sx_slug <- gsub("[^A-Za-z0-9]+", "-", tolower(sx))
      ggsave(
        file.path(dirs$plots_pub, paste0("forest_phaseAvg_ar1_groupedWindows_", run_scope, "_sex-", sx_slug, ".svg")),
        p_grouped,
        width = if (isTRUE(includeChange)) 2 * 2 * 2.3 else 2 * 2.3,
        height = 4,
        device = "svg"
      )
    }

    invisible(NULL)
  }

  main_res <- fit_model_slices_ar1("all")
  if (nrow(main_res$phase_contr) > 0) {
    main_res$phase_contr <- finalize_phase_contrast_pvals(main_res$phase_contr, method = "BH")
  }

  if (nrow(main_res$summary) > 0) {
    ar1_summary_file <- file.path(dirs$tables_sum, paste0("summary_ar1_", run_scope, ".xlsx"))
    write_pub_xlsx(main_res$summary, ar1_summary_file, "summary_ar1")
    # Add metadata sheet for AR1 analysis
    wb_ar1 <- openxlsx::loadWorkbook(ar1_summary_file)
    if (!"Metadata & Interpretation" %in% openxlsx::sheets(wb_ar1)) {
      ar1_metadata <- tibble::tibble(
        Information_Type = c(
          "Analysis Type", "Model Engine", "Formula",
          "Reference Group", "Correlation Structure",
          "Random Effects", "Coefficient Interpretation - Intercept",
          "Coefficient Interpretation - Group Effects",
          "Coefficient Interpretation - Interaction Terms",
          "Statistical Testing", "Multiple Comparisons",
          "Purpose", "Publication Notes"
        ),
        Value = c(
          "AR1 mixed-effects linear regression (nlme::lme)",
          "nlme::lme with adaptive random effects fallback sequence",
          "response ~ Group (+ optional: TH_scaled, Group:TH_scaled, + random intercept/slope)",
          "CON (Control group)",
          "Autoregressive AR1 correlation structure for repeated measures within animals",
          "Adaptive: random slope + AR1 → intercept + AR1 → slope (no AR1) → intercept (no AR1)",
          paste0("Baseline predicted value for CON group at TH_scaled=0."),
          paste0("Effect of each group (RES, SUS) relative to CON (reference). ",
                 "Similar interpretation to lmer but accounting for AR1 correlation."),
          paste0("If present (Group:TH_scaled terms): How the time trajectory differs by group"),
          "All p-values from nlme summary. df and t-values from tTable.",
          "FDR correction applied where indicated.",
          "SENSITIVITY CHECK: AR1 models account for temporal autocorrelation in within-animal measurements",
          paste0("Compare AR1 results to lmer results (see lme_fixedEffects table) to assess sensitivity to correlation structure. ",
                 "If AR1 and lmer conclusions agree, findings are robust.")
        )
      )
      openxlsx::addWorksheet(wb_ar1, "Metadata & Interpretation")
      openxlsx::writeData(wb_ar1, "Metadata & Interpretation", ar1_metadata)
      header_style <- openxlsx::createStyle(
        fgFill = "#E7E6E6", textDecoration = "bold", border = "TopBottomLeftRight"
      )
      openxlsx::addStyle(wb_ar1, "Metadata & Interpretation", header_style, rows = 1, cols = 1:2)
      openxlsx::setColWidths(wb_ar1, "Metadata & Interpretation", cols = 1:2, widths = c(35, 65))
      openxlsx::saveWorkbook(wb_ar1, ar1_summary_file, overwrite = TRUE)
    }
  }
  if (nrow(main_res$emmeans) > 0) {
    write_pub_xlsx(main_res$emmeans, file.path(dirs$tables_emm, paste0("emmeans_ar1_", run_scope, ".xlsx")), "emmeans_ar1")
    strat_tbl <- summarize_contrast_stratification(main_res$emmeans)
    wb_emm <- openxlsx::loadWorkbook(file.path(dirs$tables_emm, paste0("emmeans_ar1_", run_scope, ".xlsx")))
    openxlsx::addWorksheet(wb_emm, "contrast_stratification")
    openxlsx::writeData(wb_emm, "contrast_stratification", strat_tbl)
    openxlsx::saveWorkbook(wb_emm, file.path(dirs$tables_emm, paste0("emmeans_ar1_", run_scope, ".xlsx")), overwrite = TRUE)
  }
  if (nrow(main_res$phase_means) > 0) {
    write_pub_xlsx(main_res$phase_means, file.path(dirs$tables_phase, paste0("phaseAvg_means_ar1_", run_scope, ".xlsx")), "phaseAvg_means_ar1")
  }
  if (nrow(main_res$phase_contr) > 0) {
    write_pub_xlsx(main_res$phase_contr, file.path(dirs$tables_phase, paste0("phaseAvg_contrasts_ar1_", run_scope, ".xlsx")), "phaseAvg_contrasts_ar1")
  }

  save_ar1_publication_plots(main_res, suffix = "")

  if (includePhase) {
    first2h_res <- fit_model_slices_ar1("first2h")
    last2h_res  <- fit_model_slices_ar1("last2h")

    if (nrow(first2h_res$phase_contr) > 0) first2h_res$phase_contr <- finalize_phase_contrast_pvals(first2h_res$phase_contr, method = "BH")
    if (nrow(last2h_res$phase_contr)  > 0) last2h_res$phase_contr  <- finalize_phase_contrast_pvals(last2h_res$phase_contr, method = "BH")

    if (nrow(first2h_res$summary) > 0) {
      ar1_first2h_file <- file.path(dirs$tables_sum, paste0("summary_ar1_", run_scope, "_first2h.xlsx"))
      write_pub_xlsx(first2h_res$summary, ar1_first2h_file, "summary_ar1_first2h")
      # Add metadata sheet
      wb_ar1_f2h <- openxlsx::loadWorkbook(ar1_first2h_file)
      if (!"Metadata & Interpretation" %in% openxlsx::sheets(wb_ar1_f2h)) {
        ar1_f2h_metadata <- tibble::tibble(
          Information_Type = c("Time Window", "Analysis Type", "Correlation Structure", "Purpose"),
          Value = c(
            "First 2 hours of each phase",
            "AR1 mixed-effects linear regression (nlme::lme) - temporal subset",
            "Autoregressive AR1 correlation structure",
            "Examine early-phase behavioral dynamics. Compare with 'last2h' and full-phase analyses."
          )
        )
        openxlsx::addWorksheet(wb_ar1_f2h, "Metadata & Interpretation")
        openxlsx::writeData(wb_ar1_f2h, "Metadata & Interpretation", ar1_f2h_metadata)
        header_style <- openxlsx::createStyle(fgFill = "#E7E6E6", textDecoration = "bold", border = "TopBottomLeftRight")
        openxlsx::addStyle(wb_ar1_f2h, "Metadata & Interpretation", header_style, rows = 1, cols = 1:2)
        openxlsx::setColWidths(wb_ar1_f2h, "Metadata & Interpretation", cols = 1:2, widths = c(30, 60))
        openxlsx::saveWorkbook(wb_ar1_f2h, ar1_first2h_file, overwrite = TRUE)
      }
    }
    if (nrow(last2h_res$summary) > 0) {
      ar1_last2h_file <- file.path(dirs$tables_sum, paste0("summary_ar1_", run_scope, "_last2h.xlsx"))
      write_pub_xlsx(last2h_res$summary, ar1_last2h_file, "summary_ar1_last2h")
      # Add metadata sheet
      wb_ar1_l2h <- openxlsx::loadWorkbook(ar1_last2h_file)
      if (!"Metadata & Interpretation" %in% openxlsx::sheets(wb_ar1_l2h)) {
        ar1_l2h_metadata <- tibble::tibble(
          Information_Type = c("Time Window", "Analysis Type", "Correlation Structure", "Purpose"),
          Value = c(
            "Last 2 hours of each phase",
            "AR1 mixed-effects linear regression (nlme::lme) - temporal subset",
            "Autoregressive AR1 correlation structure",
            "Examine late-phase behavioral dynamics. Compare with 'first2h' and full-phase analyses."
          )
        )
        openxlsx::addWorksheet(wb_ar1_l2h, "Metadata & Interpretation")
        openxlsx::writeData(wb_ar1_l2h, "Metadata & Interpretation", ar1_l2h_metadata)
        header_style <- openxlsx::createStyle(fgFill = "#E7E6E6", textDecoration = "bold", border = "TopBottomLeftRight")
        openxlsx::addStyle(wb_ar1_l2h, "Metadata & Interpretation", header_style, rows = 1, cols = 1:2)
        openxlsx::setColWidths(wb_ar1_l2h, "Metadata & Interpretation", cols = 1:2, widths = c(30, 60))
        openxlsx::saveWorkbook(wb_ar1_l2h, ar1_last2h_file, overwrite = TRUE)
      }
    }

    if (nrow(first2h_res$emmeans) > 0) {
      write_pub_xlsx(first2h_res$emmeans, file.path(dirs$tables_emm, paste0("emmeans_ar1_", run_scope, "_first2h.xlsx")), "emmeans_ar1_first2h")
    }
    if (nrow(last2h_res$emmeans) > 0) {
      write_pub_xlsx(last2h_res$emmeans, file.path(dirs$tables_emm, paste0("emmeans_ar1_", run_scope, "_last2h.xlsx")), "emmeans_ar1_last2h")
    }

    if (nrow(first2h_res$phase_contr) > 0) {
      write_pub_xlsx(first2h_res$phase_contr, file.path(dirs$tables_phase, paste0("phaseAvg_contrasts_ar1_", run_scope, "_first2h.xlsx")), "phaseAvg_contrasts_ar1_first2h")
    }
    if (nrow(last2h_res$phase_contr) > 0) {
      write_pub_xlsx(last2h_res$phase_contr, file.path(dirs$tables_phase, paste0("phaseAvg_contrasts_ar1_", run_scope, "_last2h.xlsx")), "phaseAvg_contrasts_ar1_last2h")
    }

    save_ar1_publication_plots(first2h_res, suffix = "_first2h")
    save_ar1_publication_plots(last2h_res, suffix = "_last2h")
    save_ar1_grouped_window_phase_forest(main_res, first2h_res, last2h_res)
  }

  snapshot_file_ar1 <- file.path(dirs$artifacts, paste0("snapshot_ar1_", run_scope, ".rds"))
  saveRDS(
    list(
      metric = metric_name,
      run_scope = run_scope,
      includeChange = includeChange,
      includeSex = includeSex,
      includePhase = includePhase,
      dirs = dirs,
      main_results = main_res,
      first2h_results = if (exists("first2h_res")) first2h_res else NULL,
      last2h_results = if (exists("last2h_res")) last2h_res else NULL
    ),
    snapshot_file_ar1,
    compress = "xz"
  )

  cat(sprintf("AR1 sensitivity results saved to: %s\n", results_dir))
  invisible(list(dirs = dirs, run_scope = run_scope, snapshot = snapshot_file_ar1))
}

# -------------------------------------------------
# Run analyses
# -------------------------------------------------
cat("\n === RUNNING ANALYSES FOR BEHAVIORAL METRICS ===\n")

# Run analysis for each metric using loop (reduces code duplication)
# Note: ActivityIndex always analyzed by Phase as it is derived from Change × Phase
metrics_to_analyze <- c("Movement", "Proximity", "ActivityIndex")
analysis_results <- list()

for (metric in metrics_to_analyze) {
  cat(sprintf("\n  Processing: %s\n", metric))
  analysis_results[[metric]] <- run_analysis_for_metric(
    metric,
    data_filtered_agg,
    includeChange,
    includeSex,
    includePhase = includePhase  # Consistent parameter for all metrics
  )
}

# Optional AR1 sensitivity branch (saved separately from primary lmer outputs)
RUN_AR1_SENSITIVITY <- TRUE

# Additional analyses (all 10) — runs after main + AR1 pipelines
RUN_ANALYSES <- TRUE
ar1_results <- list()
if (RUN_AR1_SENSITIVITY) {
  cat("\n === RUNNING AR1 SENSITIVITY ANALYSES (SEPARATE OUTPUTS) ===\n")
  for (metric in metrics_to_analyze) {
    cat(sprintf("\n  AR1 sensitivity: %s\n", metric))
    ar1_results[[metric]] <- run_analysis_for_metric_ar1(
      metric,
      data_filtered_agg,
      includeChange,
      includeSex,
      includePhase = includePhase
    )
  }
}

# Global pipeline snapshot (main + AR1 branches)
pipeline_artifacts_dir <- file.path(base_results_dir, "artifacts")
dir_create_safe(pipeline_artifacts_dir)
pipeline_snapshot_file <- file.path(
  pipeline_artifacts_dir,
  paste0("pipeline_snapshot_", scope_all(includeChange, includeSex, includePhase), ".rds")
)
saveRDS(
  list(
    metrics_to_analyze = metrics_to_analyze,
    includeChange = includeChange,
    includeSex = includeSex,
    includePhase = includePhase,
    analysis_results = analysis_results,
    ar1_results = ar1_results
  ),
  pipeline_snapshot_file,
  compress = "xz"
)

# Extract results for backward compatibility
movement_results <- analysis_results$Movement
proximity_results <- analysis_results$Proximity
activity_index_results <- analysis_results$ActivityIndex

cat("\n=== ANALYSIS COMPLETE ===\n")
cat(sprintf("Analyzed %d metrics: %s\n", length(metrics_to_analyze), paste(metrics_to_analyze, collapse = ", ")))
for (metric in metrics_to_analyze) {
  cat(sprintf("  ✓ %s results: %s\n", metric, analysis_results[[metric]]$dirs$models))
}
if (RUN_AR1_SENSITIVITY) {
  for (metric in metrics_to_analyze) {
    cat(sprintf("  ✓ %s AR1 sensitivity results: %s\n", metric, ar1_results[[metric]]$dirs$models))
  }
}

# -------------------------------------------------
# Integrated robustness + cross-metric agreement
# -------------------------------------------------
cat("\n=== INTEGRATED ROBUSTNESS / AGREEMENT ANALYSIS ===\n")

integration_dir <- file.path(base_results_dir, "Integrated_robustness")
dir_create_safe(integration_dir)
dirs_int <- list(
  tables = file.path(integration_dir, "tables"),
  plots  = file.path(integration_dir, "plots"),
  artifacts = file.path(integration_dir, "artifacts")
)
invisible(lapply(dirs_int, dir_create_safe))

run_scope_current <- scope_all(includeChange, includeSex, includePhase)

normalize_pair <- function(x) {
  x <- gsub("\\s+", "", as.character(x))
  x <- toupper(x)
  dplyr::case_when(
    x %in% c("RES-CON", "CON-RES") ~ "RES-CON",
    x %in% c("SUS-CON", "CON-SUS") ~ "SUS-CON",
    x %in% c("SUS-RES", "RES-SUS") ~ "SUS-RES",
    TRUE ~ NA_character_
  )
}

read_phase_contrasts <- function(path, window_label, metric_name, model_label) {
  if (!file.exists(path)) return(tibble::tibble())
  df <- tryCatch(openxlsx::read.xlsx(path), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(tibble::tibble())

  required_cols <- c("contrast", "estimate", "p.adjust")
  if (!all(required_cols %in% names(df))) return(tibble::tibble())

  df %>%
    dplyr::mutate(
      metric = metric_name,
      model = model_label,
      window = window_label,
      pair = normalize_pair(contrast),
      Change = as.character(Change),
      Sex = as.character(Sex),
      Phase = as.character(Phase),
      d = if ("d" %in% names(df)) as.numeric(d) else NA_real_,
      p.adjust = as.numeric(p.adjust),
      estimate = as.numeric(estimate),
      sig = p.adjust < 0.05
    ) %>%
    dplyr::filter(!is.na(pair)) %>%
    dplyr::select(metric, model, window, Change, Sex, Phase, pair, estimate, p.adjust, d, sig)
}

collect_metric_contrasts <- function(metric_name) {
  out <- list()

  # LMER branch
  lmer_phase_dir <- analysis_results[[metric_name]]$dirs$tables_phase
  out[[length(out) + 1]] <- read_phase_contrasts(
    file.path(lmer_phase_dir, paste0("phaseAvg_contrasts_", run_scope_current, ".xlsx")),
    "all", metric_name, "lmer"
  )
  out[[length(out) + 1]] <- read_phase_contrasts(
    file.path(lmer_phase_dir, paste0("phaseAvg_contrasts_", run_scope_current, "_first2h.xlsx")),
    "first2h", metric_name, "lmer"
  )
  out[[length(out) + 1]] <- read_phase_contrasts(
    file.path(lmer_phase_dir, paste0("phaseAvg_contrasts_", run_scope_current, "_last2h.xlsx")),
    "last2h", metric_name, "lmer"
  )

  # AR1 branch (if available)
  if (RUN_AR1_SENSITIVITY && !is.null(ar1_results[[metric_name]])) {
    ar1_phase_dir <- ar1_results[[metric_name]]$dirs$tables_phase
    out[[length(out) + 1]] <- read_phase_contrasts(
      file.path(ar1_phase_dir, paste0("phaseAvg_contrasts_ar1_", run_scope_current, ".xlsx")),
      "all", metric_name, "ar1"
    )
    out[[length(out) + 1]] <- read_phase_contrasts(
      file.path(ar1_phase_dir, paste0("phaseAvg_contrasts_ar1_", run_scope_current, "_first2h.xlsx")),
      "first2h", metric_name, "ar1"
    )
    out[[length(out) + 1]] <- read_phase_contrasts(
      file.path(ar1_phase_dir, paste0("phaseAvg_contrasts_ar1_", run_scope_current, "_last2h.xlsx")),
      "last2h", metric_name, "ar1"
    )
  }

  dplyr::bind_rows(out)
}

integrated_tbl <- dplyr::bind_rows(lapply(metrics_to_analyze, collect_metric_contrasts))

if (nrow(integrated_tbl) > 0) {
  write_pub_xlsx(integrated_tbl,
    file.path(dirs_int$tables, paste0("integrated_phase_contrasts_", run_scope_current, ".xlsx")),
    sheet_name = "Results", p_cols = character(0))

  # ------------------------------
  # A) Movement vs ActivityIndex agreement
  # ------------------------------
  mov_act_tbl <- integrated_tbl %>%
    dplyr::filter(metric %in% c("Movement", "ActivityIndex")) %>%
    tidyr::pivot_wider(
      names_from = metric,
      values_from = c(estimate, p.adjust, d, sig)
    ) %>%
    dplyr::filter(!is.na(estimate_Movement), !is.na(estimate_ActivityIndex))

  if (nrow(mov_act_tbl) > 0) {
    mov_act_summary <- mov_act_tbl %>%
      dplyr::group_by(model, window) %>%
      dplyr::summarise(
        n = dplyr::n(),
        r_estimate = ifelse(dplyr::n() >= 3,
                            suppressWarnings(cor(estimate_Movement, estimate_ActivityIndex, use = "complete.obs")),
                            NA_real_),
        r_d = ifelse(sum(!is.na(d_Movement) & !is.na(d_ActivityIndex)) >= 3,
                     suppressWarnings(cor(d_Movement, d_ActivityIndex, use = "complete.obs")),
                     NA_real_),
        sign_match_rate = mean(sign(estimate_Movement) == sign(estimate_ActivityIndex), na.rm = TRUE),
        sig_match_rate = mean(sig_Movement == sig_ActivityIndex, na.rm = TRUE),
        .groups = "drop"
      )

    write_pub_xlsx(mov_act_summary,
      file.path(dirs_int$tables, paste0("movement_activity_similarity_summary_", run_scope_current, ".xlsx")),
      sheet_name = "Results", p_cols = character(0))

    p_mov_act <- ggplot(
      mov_act_tbl,
      aes(x = estimate_Movement, y = estimate_ActivityIndex, color = pair, shape = Phase)
    ) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey75") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey75") +
      geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black") +
      geom_point(size = 2.5, alpha = 0.85) +
      scale_color_manual(values = pair_colors) +
      facet_grid(model ~ window) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = "top",
        plot.title = element_text(face = "bold", hjust = 0.5)
      ) +
      labs(
        title = "Movement vs ActivityIndex agreement (phase-contrast estimates)",
        x = "Movement estimate",
        y = "ActivityIndex estimate",
        color = "Contrast",
        shape = "Phase"
      )

    ggsave(
      file.path(dirs_int$plots, paste0("movement_activity_similarity_scatter_", run_scope_current, ".svg")),
      p_mov_act,
      width = 8,
      height = 6
    )
  }

  # ------------------------------
  # B) LMER vs AR1 robustness
  # ------------------------------
  if (RUN_AR1_SENSITIVITY) {
    robust_tbl <- integrated_tbl %>%
      tidyr::pivot_wider(
        names_from = model,
        values_from = c(estimate, p.adjust, d, sig)
      ) %>%
      dplyr::filter(!is.na(estimate_lmer), !is.na(estimate_ar1)) %>%
      dplyr::mutate(
        delta_estimate = estimate_ar1 - estimate_lmer,
        abs_delta = abs(delta_estimate),
        sign_match = sign(estimate_ar1) == sign(estimate_lmer),
        sig_match = sig_ar1 == sig_lmer
      )

    if (nrow(robust_tbl) > 0) {
      robust_summary <- robust_tbl %>%
        dplyr::group_by(metric, window, pair) %>%
        dplyr::summarise(
          n = dplyr::n(),
          mean_abs_delta = mean(abs_delta, na.rm = TRUE),
          sign_match_rate = mean(sign_match, na.rm = TRUE),
          sig_match_rate = mean(sig_match, na.rm = TRUE),
          .groups = "drop"
        )

      write_pub_xlsx(robust_tbl,
        file.path(dirs_int$tables, paste0("lmer_vs_ar1_robustness_detail_", run_scope_current, ".xlsx")),
        sheet_name = "Results", p_cols = character(0))
      write_pub_xlsx(robust_summary,
        file.path(dirs_int$tables, paste0("lmer_vs_ar1_robustness_summary_", run_scope_current, ".xlsx")),
        sheet_name = "Results", p_cols = character(0))

      p_robust_scatter <- ggplot(
        robust_tbl,
        aes(x = estimate_lmer, y = estimate_ar1, color = pair)
      ) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey75") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey75") +
        geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "black") +
        geom_point(size = 2.4, alpha = 0.85) +
        scale_color_manual(values = pair_colors) +
        facet_grid(metric ~ window) +
        theme_minimal(base_size = 12) +
        theme(
          panel.grid.minor = element_blank(),
          legend.position = "top",
          plot.title = element_text(face = "bold", hjust = 0.5)
        ) +
        labs(
          title = "Robustness across model assumptions (LMER vs AR1)",
          x = "LMER estimate",
          y = "AR1 estimate",
          color = "Contrast"
        )

      ggsave(
        file.path(dirs_int$plots, paste0("lmer_vs_ar1_robustness_scatter_", run_scope_current, ".svg")),
        p_robust_scatter,
        width = 9,
        height = 6
      )

      p_robust_heat <- ggplot(
        robust_summary,
        aes(x = window, y = pair, fill = sign_match_rate)
      ) +
        geom_tile(color = "white", linewidth = 0.4) +
        geom_text(aes(label = sprintf("%.2f", sign_match_rate)), size = 3.2, fontface = "bold") +
        scale_fill_gradient(low = "#f3f3f3", high = "#2a9d8f", limits = c(0, 1)) +
        facet_wrap(~ metric) +
        theme_minimal(base_size = 12) +
        theme(
          panel.grid = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "right",
          plot.title = element_text(face = "bold", hjust = 0.5)
        ) +
        labs(
          title = "Sign-direction stability (LMER vs AR1)",
          x = "Window",
          y = "Contrast",
          fill = "Sign match"
        )

      ggsave(
        file.path(dirs_int$plots, paste0("lmer_vs_ar1_signMatch_heatmap_", run_scope_current, ".svg")),
        p_robust_heat,
        width = 8,
        height = 4.8
      )
    }
  }

  cat(sprintf("  ✓ Integrated robustness outputs saved to: %s\n", integration_dir))
}

integrated_snapshot_file <- file.path(dirs_int$artifacts, paste0("integrated_snapshot_", run_scope_current, ".rds"))
saveRDS(
  list(
    run_scope = run_scope_current,
    integrated_tbl = integrated_tbl,
    movement_activity_tbl = if (exists("mov_act_tbl")) mov_act_tbl else tibble::tibble(),
    movement_activity_summary = if (exists("mov_act_summary")) mov_act_summary else tibble::tibble(),
    robust_tbl = if (exists("robust_tbl")) robust_tbl else tibble::tibble(),
    robust_summary = if (exists("robust_summary")) robust_summary else tibble::tibble()
  ),
  integrated_snapshot_file,
  compress = "xz"
)


# =================================================================
# Additional analyses (all 10)
# =================================================================
if (RUN_ANALYSES) {
  cat("\n=== ADDITIONAL ANALYSES ===\n")

  analyses_dir <- file.path(base_results_dir, "analyses")
  dir_create_safe(analyses_dir)
  nat_dirs <- list(
    # Spline analyses (lmer and AR1)
    splines_lmer_tables = file.path(analyses_dir, "splines", "lmer", "tables"),
    splines_lmer_plots  = file.path(analyses_dir, "splines", "lmer", "plots"),
    splines_ar1_tables  = file.path(analyses_dir, "splines", "ar1",  "tables"),
    splines_ar1_plots   = file.path(analyses_dir, "splines", "ar1",  "plots"),
    # GAMM analyses
    gamm_tables = file.path(analyses_dir, "gamm", "tables"),
    gamm_plots  = file.path(analyses_dir, "gamm", "plots"),
    # General outputs
    tables = file.path(analyses_dir, "tables"),
    plots  = file.path(analyses_dir, "plots"),
    artifacts = file.path(analyses_dir, "artifacts")
  )
  invisible(lapply(nat_dirs, dir_create_safe))

  # ------------------------------------------------------------------
  # Helper: re-read phase contrast xlsx with full CI/SE/d columns
  # ------------------------------------------------------------------
  read_phase_contrasts_full <- function(path, window_label, metric_name, model_label) {
    if (!file.exists(path)) return(tibble::tibble())
    df <- tryCatch(openxlsx::read.xlsx(path), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) return(tibble::tibble())
    if (!all(c("contrast", "estimate", "p.adjust") %in% names(df))) return(tibble::tibble())

    df %>%
      dplyr::mutate(
        metric   = metric_name,
        model    = model_label,
        window   = window_label,
        pair     = normalize_pair(contrast),
        Change   = as.character(Change),
        Sex      = as.character(Sex),
        Phase    = as.character(Phase),
        SE       = if ("SE"    %in% names(df)) as.numeric(SE)    else NA_real_,
        lwr      = if ("lwr"   %in% names(df)) as.numeric(lwr)   else NA_real_,
        upr      = if ("upr"   %in% names(df)) as.numeric(upr)   else NA_real_,
        d        = if ("d"     %in% names(df)) as.numeric(d)     else NA_real_,
        lwr_d    = if ("lwr_d" %in% names(df)) as.numeric(lwr_d) else NA_real_,
        upr_d    = if ("upr_d" %in% names(df)) as.numeric(upr_d) else NA_real_,
        p.adjust = as.numeric(p.adjust),
        estimate = as.numeric(estimate),
        sig      = p.adjust < 0.05
      ) %>%
      dplyr::filter(!is.na(pair)) %>%
      dplyr::select(metric, model, window, Change, Sex, Phase, pair,
                    estimate, SE, lwr, upr, d, lwr_d, upr_d, p.adjust, sig)
  }

  collect_metric_contrasts_full <- function(metric_name) {
    out     <- list()
    ldir    <- analysis_results[[metric_name]]$dirs$tables_phase
    suffixes <- c(all = "", first2h = "_first2h", last2h = "_last2h")

    for (wname in names(suffixes)) {
      wsuf <- suffixes[[wname]]
      out[[length(out) + 1]] <- read_phase_contrasts_full(
        file.path(ldir, paste0("phaseAvg_contrasts_", run_scope_current, wsuf, ".xlsx")),
        wname, metric_name, "lmer"
      )
      if (RUN_AR1_SENSITIVITY && !is.null(ar1_results[[metric_name]])) {
        adir <- ar1_results[[metric_name]]$dirs$tables_phase
        out[[length(out) + 1]] <- read_phase_contrasts_full(
          file.path(adir, paste0("phaseAvg_contrasts_ar1_", run_scope_current, wsuf, ".xlsx")),
          wname, metric_name, "ar1"
        )
      }
    }
    dplyr::bind_rows(out)
  }

  integrated_full <- dplyr::bind_rows(lapply(metrics_to_analyze, collect_metric_contrasts_full))

  # ------------------------------------------------------------------
  # Helper: fit a fast representative lmer per metric/change/sex/phase
  # (intercept-only random effect; used for heterogeneity, ICC, LOO, power)
  # ------------------------------------------------------------------
  fit_representative_models <- function() {
    cat("  Fitting representative models...\n")
    ctrl_rep <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
    models_list <- list()

    changes_r <- if (includeChange) unique(data_filtered_agg$Change) else "allChanges"
    sexes_r   <- if (includeSex)    unique(data_filtered_agg$Sex)    else "allSexes"
    phases_r  <- if (includePhase)  unique(data_filtered_agg$Phase)  else "allPhases"

    for (metric in metrics_to_analyze) {
      models_list[[metric]] <- list()
      for (ch in changes_r) {
        for (sx in sexes_r) {
          for (ph in phases_r) {
            key <- paste(ch, sx, ph, sep = "_")
            dsub <- data_filtered_agg %>%
              dplyr::filter(
                (!includeChange | Change == ch),
                (!includeSex    | Sex    == sx),
                (!includePhase  | Phase  == ph),
                !is.na(.data[[metric]]), !is.na(AnimalNum), !is.na(Group)
              )
            if (nrow(dsub) < min_obs ||
                dplyr::n_distinct(dsub$AnimalNum) < min_animals ||
                dplyr::n_distinct(dsub$Group) < min_groups) next

            t0 <- min(dsub$HalfHourElapsed, na.rm = TRUE)
            dsub <- dsub %>%
              dplyr::mutate(
                TH        = (HalfHourElapsed - t0) / 2,
                TH_scaled = as.vector(scale(TH))
              ) %>%
              dplyr::filter(!is.na(TH_scaled))

            use_int   <- use_group_time_interaction && dplyr::n_distinct(dsub$TH_scaled) >= min_time_lv
            fixed_rhs <- if (use_int) "Group * TH_scaled" else "Group + TH_scaled"

            f <- as.formula(paste0(metric, " ~ ", fixed_rhs, " + (1 | AnimalNum)"))
            m <- tryCatch(
              suppressMessages(suppressWarnings(
                lmerTest::lmer(f, data = dsub, control = ctrl_rep)
              )),
              error = function(e) NULL
            )
            if (!is.null(m)) {
              models_list[[metric]][[key]] <- list(
                model  = m,
                data   = dsub,
                Change = as.character(ch),
                Sex    = as.character(sx),
                Phase  = as.character(ph)
              )
            }
          }
        }
      }
    }
    models_list
  }

  rep_models <- fit_representative_models()
  saveRDS(rep_models, file.path(nat_dirs$artifacts, paste0("rep_models_", run_scope_current, ".rds")), compress = "xz")

  # ================================================================
  # #1: Extended lmer-vs-AR1 comparison table (estimate + SE + CI + AIC/BIC)
  # ================================================================
  cat("  [1/10] lmer vs AR1 full comparison table...\n")

  if (nrow(integrated_full) > 0 && RUN_AR1_SENSITIVITY) {
    full_comp <- integrated_full %>%
      tidyr::pivot_wider(
        names_from  = model,
        values_from = c(estimate, SE, lwr, upr, d, lwr_d, upr_d, p.adjust, sig)
      ) %>%
      dplyr::filter(!is.na(estimate_lmer), !is.na(estimate_ar1)) %>%
      dplyr::mutate(
        delta_estimate = estimate_ar1 - estimate_lmer,
        sign_match     = sign(estimate_ar1) == sign(estimate_lmer),
        sig_match      = sig_ar1 == sig_lmer
      )

    aic_bic_rows <- dplyr::bind_rows(lapply(metrics_to_analyze, function(metric) {
      dplyr::bind_rows(lapply(names(rep_models[[metric]]), function(key) {
        m_obj <- rep_models[[metric]][[key]]
        if (is.null(m_obj)) return(tibble::tibble())
        tibble::tibble(
          metric   = metric,
          Change   = m_obj$Change,
          Sex      = m_obj$Sex,
          Phase    = m_obj$Phase,
          AIC_lmer = tryCatch(AIC(m_obj$model), error = function(e) NA_real_),
          BIC_lmer = tryCatch(BIC(m_obj$model), error = function(e) NA_real_)
        )
      }))
    }))

    if (nrow(aic_bic_rows) > 0) {
      full_comp <- full_comp %>%
        dplyr::left_join(aic_bic_rows, by = c("metric", "Change", "Sex", "Phase"))
    }

    write_pub_xlsx(full_comp,
      file.path(nat_dirs$splines_lmer_tables, paste0("lmer_vs_ar1_full_comparison.xlsx")),
      sheet_name = "Results", p_cols = character(0))
    cat("     \u2713 Saved\n")
  } else {
    cat("     - Skipped (AR1 not run or no integrated data)\n")
  }

  # ================================================================
  # #2: SUS-CON / SUS-RES contrast stability heatmap
  # ================================================================
  cat("  [2/10] Contrast stability heatmap (SUS-CON / SUS-RES)...\n")

  if (nrow(integrated_full) > 0) {
    stab_tbl <- integrated_full %>%
      dplyr::filter(pair %in% c("SUS-CON", "SUS-RES"), model == "lmer") %>%
      dplyr::mutate(
        row_label = paste(metric, pair, sep = "\n"),
        col_label = paste0(Sex, "/", Phase, "/", window),
        sig_star  = dplyr::case_when(
          p.adjust < 0.001 ~ "***",
          p.adjust < 0.01  ~ "**",
          p.adjust < 0.05  ~ "*",
          TRUE             ~ ""
        )
      )

    if (nrow(stab_tbl) > 0) {
      p_stab_heat <- ggplot(stab_tbl, aes(x = col_label, y = row_label, fill = d)) +
        geom_tile(color = "white", linewidth = 0.35) +
        geom_text(aes(label = sig_star), size = 2.5, vjust = 0.75) +
        scale_fill_gradient2(
          low = "#2c7bb6", mid = "white", high = "#d7191c",
          midpoint = 0, na.value = "grey90", name = "Cohen's d"
        ) +
        theme_minimal(base_size = 9) +
        theme(
          axis.text.x   = element_text(angle = 55, hjust = 1, size = 7),
          axis.text.y   = element_text(size = 8),
          panel.grid    = element_blank(),
          plot.title    = element_text(face = "bold", hjust = 0.5),
          legend.position = "right"
        ) +
        labs(
          title    = "Contrast stability: SUS-CON and SUS-RES",
          subtitle = "Rows = metric x contrast | Columns = sex x phase x window | * p<0.05 ** p<0.01 *** p<0.001",
          x = NULL, y = NULL
        )

      ggsave(
        file.path(nat_dirs$plots, paste0("contrast_stability_heatmap_", run_scope_current, ".svg")),
        p_stab_heat,
        width  = max(8, dplyr::n_distinct(stab_tbl$col_label) * 0.65),
        height = max(4, dplyr::n_distinct(stab_tbl$row_label) * 0.6)
      )
      write_pub_xlsx(stab_tbl,
        file.path(nat_dirs$tables, paste0("contrast_stability_data_", run_scope_current, ".xlsx")),
        sheet_name = "Results", p_cols = character(0))
      cat("     \u2713 Saved\n")
    }
  }

  # ================================================================
  # #3: Sex-specific DID decomposition
  # ================================================================
  cat("  [3/10] Sex-specific DID decomposition...\n")

  if (nrow(integrated_full) > 0 && includeSex && dplyr::n_distinct(integrated_full$Sex) >= 2) {
    did_base <- integrated_full %>%
      dplyr::filter(model == "lmer", !is.na(estimate), !is.na(SE)) %>%
      dplyr::select(metric, window, Change, Phase, pair, Sex, estimate, SE)

    sex_vals <- sort(unique(as.character(did_base$Sex)))
    s1 <- sex_vals[1]; s2 <- sex_vals[2]

    did_wide <- did_base %>%
      tidyr::pivot_wider(
        names_from  = Sex,
        values_from = c(estimate, SE),
        names_sep   = "__"
      )

    est_s1 <- paste0("estimate__", s1); se_s1 <- paste0("SE__", s1)
    est_s2 <- paste0("estimate__", s2); se_s2 <- paste0("SE__", s2)

    if (all(c(est_s1, est_s2, se_s1, se_s2) %in% names(did_wide))) {
      did_tbl <- did_wide %>%
        dplyr::filter(!is.na(.data[[est_s1]]), !is.na(.data[[est_s2]])) %>%
        dplyr::mutate(
          DID            = .data[[est_s1]] - .data[[est_s2]],
          SE_DID         = sqrt(.data[[se_s1]]^2 + .data[[se_s2]]^2),
          z_DID          = DID / SE_DID,
          p_DID          = 2 * pnorm(-abs(z_DID)),
          p_DID_BH       = p.adjust(p_DID, method = "BH"),
          sig_DID        = p_DID_BH < 0.05,
          sex_comparison = paste0(s1, " - ", s2)
        )

      write_pub_xlsx(did_tbl,
        file.path(nat_dirs$tables, paste0("sex_did_decomposition_", run_scope_current, ".xlsx")),
        sheet_name = "Results", p_cols = character(0))

      did_plot_tbl <- did_tbl %>%
        dplyr::filter(window == "all", !is.na(DID)) %>%
        dplyr::mutate(
          lwr      = DID - 1.96 * SE_DID,
          upr      = DID + 1.96 * SE_DID,
          pair     = factor(pair, levels = c("RES-CON", "SUS-CON", "SUS-RES")),
          sig_star = dplyr::case_when(
            p_DID_BH < 0.001 ~ "***", p_DID_BH < 0.01 ~ "**",
            p_DID_BH < 0.05  ~ "*",   TRUE             ~ ""
          )
        )

      if (nrow(did_plot_tbl) > 0) {
        p_did <- ggplot(
          did_plot_tbl,
          aes(x = DID, xmin = lwr, xmax = upr, y = pair, color = pair)
        ) +
          geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
          geom_pointrange(size = 0.5, linewidth = 0.9) +
          geom_text(aes(x = upr, label = sig_star), hjust = -0.2, size = 3.2, color = "black", show.legend = FALSE) +
          scale_color_manual(values = pair_colors) +
          facet_grid(metric ~ Phase) +
          theme_classic(base_size = 9) +
          theme(
            axis.line.y = element_blank(), axis.ticks.y = element_blank(),
            axis.text.y = element_blank(), legend.position = "bottom",
            strip.text = element_text(face = "bold"),
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.margin = margin(5, 45, 5, 5)
          ) +
          scale_x_continuous(expand = expansion(mult = c(0.1, 0.4))) +
          labs(
            title    = paste0("Sex-specific DID: ", s1, " \u2212 ", s2),
            subtitle = "Positive = larger effect in first sex group",
            x        = paste0("\u0394 (", s1, " \u2212 ", s2, ") estimate"),
            y        = NULL, color = "Contrast"
          )

        ggsave(
          file.path(nat_dirs$plots, paste0("sex_did_forest_", run_scope_current, ".svg")),
          p_did, width = 7, height = 5
        )
      }
      cat("     \u2713 Saved\n")
    }
  } else {
    cat("     - Skipped (sex stratification not active or < 2 sex levels)\n")
  }

  # ================================================================
  # #4: Spline trajectory analysis
  # ================================================================
  cat("  [4/10] Spline trajectory analysis...\n")

  library(splines)  # base R package, always available
  RUN_SPLINE_AUC_BOOTSTRAP <- TRUE
  N_SPLINE_AUC_BOOT <- 1000L

  trapz_auc <- function(x, y) {
    if (length(x) < 2 || length(y) < 2) return(NA_real_)
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]; y <- y[keep]
    if (length(x) < 2) return(NA_real_)
    ord <- order(x)
    x <- x[ord]; y <- y[ord]
    sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  }

  # ---- Helper: Extract and format fixed effects with interpretation ----
  extract_fixed_effects_with_metadata <- function(model_obj, model_type = "lmer",
                                                   metric, Change, Sex, Phase,
                                                   data_subset, group_levels) {
    # Extract fixed effects
    if (model_type == "lmer") {
      fixef_tbl <- lmerTest::as_lmerModLmerTest(model_obj) %>%
        summary() %>%
        .$coefficients %>%
        as.data.frame() %>%
        tibble::rownames_to_column("Coefficient") %>%
        tibble::as_tibble()
      names(fixef_tbl) <- c("Coefficient", "Estimate", "SE", "df", "t_value", "p_value")
    } else if (model_type == "ar1") {
      # For nlme::lme AR1 models
      coef_tbl <- summary(model_obj)$tTable
      fixef_tbl <- coef_tbl %>%
        as.data.frame() %>%
        tibble::rownames_to_column("Coefficient") %>%
        tibble::as_tibble()
      names(fixef_tbl) <- c("Coefficient", "Value", "SE", "df", "t_value", "p_value")
      fixef_tbl <- fixef_tbl %>%
        dplyr::rename(Estimate = Value)
    } else {
      return(list(fixed_effects = tibble::tibble(), metadata = tibble::tibble()))
    }

    # Add context columns
    fixef_tbl <- fixef_tbl %>%
      dplyr::mutate(
        metric = metric,
        Change = as.character(Change),
        Sex = as.character(Sex),
        Phase = as.character(Phase),
        model_type = model_type,
        .before = Coefficient
      )

    # Create interpretation guide
    group_ref <- group_levels[1]  # First group is reference
    group_other <- setdiff(group_levels, group_ref)

    # Create metadata table
    metadata <- tibble::tibble(
      Information_Type = c(
        "Model Type",
        "Formula",
        "Reference Group",
        "Data Subset - Metric",
        "Data Subset - Change",
        "Data Subset - Sex",
        "Data Subset - Phase",
        "Data Subset - N Animals",
        "Data Subset - N Observations",
        "Groups in Data",
        "Coefficient Interpretation - Intercept",
        "Coefficient Interpretation - Group Effects",
        "Coefficient Interpretation - Interaction Terms",
        "Statistical Inference",
        "Publication Note"
      ),
      Value = c(
        if_else(model_type == "lmer", "Mixed effects model (lmerTest::lmer)", "AR1 mixed model (nlme::lme with corAR1)"),
        if_else(model_type == "lmer",
                paste0(metric, " ~ Group * ns(TH_scaled, df) + (1|AnimalNum)"),
                paste0(metric, " ~ Group * ns(TH_scaled, df)")),
        paste0(group_ref, " (CON=Control)"),
        metric,
        as.character(Change),
        as.character(Sex),
        as.character(Phase),
        as.character(dplyr::n_distinct(data_subset$AnimalNum)),
        as.character(nrow(data_subset)),
        paste0(sort(unique(as.character(data_subset$Group))), collapse = ", "),
        paste0("Baseline response for ", group_ref, " at TH_scaled=0 (natural cubic spline intercept)"),
        paste0("Difference in intercept between each group and reference (", group_ref, "). ",
               "E.g., GroupRES = effect of RES vs ", group_ref, " group on baseline response"),
        "Spline coefficients × group interaction: how the time trajectory differs by group",
        "All p-values from lmerTest/nlme. 95% CI = Estimate ± 1.96*SE",
        "For interpretation: use AUC bootstrap contrasts table for group comparisons"
      )
    )

    list(fixed_effects = fixef_tbl, metadata = metadata)
  }

  bootstrap_auc_contrasts <- function(model_obj, th_min, th_max, group_levels,
                                      metric, Change, Sex, Phase,
                                      model_label = "lmer", n_boot = 1000L) {
    if (length(group_levels) < 2 || !is.finite(th_min) || !is.finite(th_max) || th_max <= th_min) {
      return(tibble::tibble())
    }

    fixed_formula <- if (inherits(model_obj, "merMod")) {
      lme4::nobars(formula(model_obj))
    } else {
      formula(model_obj)
    }

    beta_hat <- tryCatch(as.numeric(nlme::fixef(model_obj)), error = function(e) NULL)
    beta_nms <- tryCatch(names(nlme::fixef(model_obj)), error = function(e) NULL)
    if (is.null(beta_hat) || is.null(beta_nms)) {
      beta_hat <- tryCatch(as.numeric(lme4::fixef(model_obj)), error = function(e) NULL)
      beta_nms <- tryCatch(names(lme4::fixef(model_obj)), error = function(e) NULL)
    }
    if (is.null(beta_hat) || is.null(beta_nms)) return(tibble::tibble())
    names(beta_hat) <- beta_nms

    V <- tryCatch(as.matrix(vcov(model_obj)), error = function(e) NULL)
    if (is.null(V) || nrow(V) == 0) return(tibble::tibble())

    group_levels <- unique(as.character(group_levels))
    group_pref <- c("CON", "RES", "SUS")
    group_levels <- c(group_pref[group_pref %in% group_levels], setdiff(group_levels, group_pref))

    pred_grid <- expand.grid(
      TH_scaled = seq(th_min, th_max, length.out = 50),
      Group = group_levels,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::mutate(Group = factor(Group, levels = group_levels))

    X <- tryCatch(
      model.matrix(stats::delete.response(stats::terms(fixed_formula)), data = pred_grid),
      error = function(e) NULL
    )
    if (is.null(X) || ncol(X) == 0) return(tibble::tibble())

    common <- intersect(colnames(X), names(beta_hat))
    if (length(common) < 2) return(tibble::tibble())

    X <- X[, common, drop = FALSE]
    beta_hat <- beta_hat[common]
    V <- V[common, common, drop = FALSE]

    draws <- tryCatch(MASS::mvrnorm(n = n_boot, mu = beta_hat, Sigma = V), error = function(e) NULL)
    if (is.null(draws)) return(tibble::tibble())
    if (!is.matrix(draws)) draws <- matrix(draws, nrow = 1)
    n_draws <- nrow(draws)

    pred_mat <- X %*% t(draws)
    auc_draws <- matrix(NA_real_, nrow = n_draws, ncol = length(group_levels),
                        dimnames = list(NULL, group_levels))

    for (g in group_levels) {
      idx <- which(as.character(pred_grid$Group) == g)
      xg <- pred_grid$TH_scaled[idx]
      auc_draws[, g] <- apply(pred_mat[idx, , drop = FALSE], 2, function(y) trapz_auc(xg, y))
    }

    pred_obs <- as.vector(X %*% beta_hat)
    auc_obs <- sapply(group_levels, function(g) {
      idx <- which(as.character(pred_grid$Group) == g)
      trapz_auc(pred_grid$TH_scaled[idx], pred_obs[idx])
    })

    pairs <- list(
      c("RES", "CON"),
      c("SUS", "CON"),
      c("SUS", "RES")
    )

    out <- lapply(pairs, function(pp) {
      g1 <- pp[1]; g0 <- pp[2]
      if (!(g1 %in% colnames(auc_draws) && g0 %in% colnames(auc_draws))) return(NULL)

      dboot <- auc_draws[, g1] - auc_draws[, g0]
      dboot <- dboot[is.finite(dboot)]
      if (length(dboot) < 10) return(NULL)

      p_lo <- mean(dboot <= 0)
      p_hi <- mean(dboot >= 0)
      p_two <- min(1, 2 * min(p_lo, p_hi))

      tibble::tibble(
        metric = metric,
        Change = as.character(Change),
        Sex = as.character(Sex),
        Phase = as.character(Phase),
        model = model_label,
        contrast = paste0(g1, "-", g0),
        auc_diff = as.numeric(auc_obs[g1] - auc_obs[g0]),
        ci_low = as.numeric(stats::quantile(dboot, 0.025, na.rm = TRUE)),
        ci_high = as.numeric(stats::quantile(dboot, 0.975, na.rm = TRUE)),
        p_boot = as.numeric(p_two),
        n_boot = as.integer(length(dboot))
      )
    })

    dplyr::bind_rows(out)
  }

  # ---- Helper: Extract individual animal AUC values for correlation ----
  extract_individual_animal_auc <- function(model_obj, th_min, th_max,
                                            metric, Change, Sex, Phase,
                                            data_subset, model_label = "lmer") {
    # Get unique animals and their groups
    animal_info <- data_subset %>%
      dplyr::distinct(AnimalNum, Group, Batch) %>%
      dplyr::arrange(AnimalNum)

    if (nrow(animal_info) == 0) return(tibble::tibble())

    # Create prediction grid for AUC calculation
    pred_grid <- expand.grid(
      TH_scaled = seq(th_min, th_max, length.out = 100),
      AnimalNum = animal_info$AnimalNum,
      stringsAsFactors = FALSE
    )

    # Get predictions with random effects
    pred_with_re <- tryCatch({
      pred_list <- lapply(animal_info$AnimalNum, function(anim) {
        td <- pred_grid[pred_grid$AnimalNum == anim, ]
        p <- predict(model_obj, newdata = td, re.form = NULL)
        data.frame(td, predicted = p)
      })
      dplyr::bind_rows(pred_list)
    }, error = function(e) NULL)

    if (is.null(pred_with_re) || nrow(pred_with_re) == 0) return(tibble::tibble())

    # Calculate AUC for each animal
    auc_individual <- pred_with_re %>%
      dplyr::group_by(AnimalNum) %>%
      dplyr::summarise(
        AUC = trapz_auc(TH_scaled, predicted),
        .groups = "drop"
      ) %>%
      dplyr::left_join(animal_info, by = "AnimalNum") %>%
      dplyr::mutate(
        metric = metric,
        Change = as.character(Change),
        Sex = as.character(Sex),
        Phase = as.character(Phase),
        model = model_label,
        .before = AnimalNum
      )

    return(auc_individual)
  }

  spline_results <- list()
  spline_auc_boot_rows <- list()
  spline_auc_individual_rows <- list()  # NEW: collect individual animal AUC values
  ctrl_sp <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))

  changes_sp <- if (includeChange) unique(data_filtered_agg$Change) else "allChanges"
  sexes_sp   <- if (includeSex)    unique(data_filtered_agg$Sex)    else "allSexes"
  phases_sp  <- if (includePhase)  unique(data_filtered_agg$Phase)  else "allPhases"

  for (metric in metrics_to_analyze) {
    for (ch in changes_sp) {
      for (sx in sexes_sp) {
        for (ph in phases_sp) {
          dsub <- data_filtered_agg %>%
            dplyr::filter(
              (!includeChange | Change == ch),
              (!includeSex    | Sex    == sx),
              (!includePhase  | Phase  == ph),
              !is.na(.data[[metric]]), !is.na(AnimalNum), !is.na(Group)
            )
          if (nrow(dsub) < max(min_obs, 15) ||
              dplyr::n_distinct(dsub$AnimalNum) < min_animals ||
              dplyr::n_distinct(dsub$Group) < min_groups) next

          t0 <- min(dsub$HalfHourElapsed, na.rm = TRUE)
          dsub <- dsub %>%
            dplyr::mutate(
              TH        = (HalfHourElapsed - t0) / 2,
              TH_scaled = as.vector(scale(TH))
            ) %>%
            dplyr::filter(!is.na(TH_scaled))

          n_uniq_t <- dplyr::n_distinct(dsub$TH_scaled)
          if (n_uniq_t < 6) next

          df_sp <- min(4L, floor(n_uniq_t / 3L))
          f_sp <- as.formula(paste0(
            metric, " ~ Group * ns(TH_scaled, df = ", df_sp, ") + (1 | AnimalNum)"
          ))

          m_sp <- tryCatch(
            suppressMessages(suppressWarnings(
              lmerTest::lmer(f_sp, data = dsub, control = ctrl_sp)
            )),
            error = function(e) NULL
          )
          if (is.null(m_sp)) next

          if (RUN_SPLINE_AUC_BOOTSTRAP) {
            boot_tbl <- bootstrap_auc_contrasts(
              model_obj = m_sp,
              th_min = min(dsub$TH_scaled, na.rm = TRUE),
              th_max = max(dsub$TH_scaled, na.rm = TRUE),
              group_levels = unique(as.character(dsub$Group)),
              metric = metric,
              Change = as.character(ch),
              Sex = as.character(sx),
              Phase = as.character(ph),
              model_label = "lmer",
              n_boot = N_SPLINE_AUC_BOOT
            )
            if (nrow(boot_tbl) > 0) {
              spline_auc_boot_rows[[length(spline_auc_boot_rows) + 1]] <- boot_tbl
            }

            # NEW: Extract individual animal AUC values
            ind_auc_tbl <- extract_individual_animal_auc(
              model_obj = m_sp,
              th_min = min(dsub$TH_scaled, na.rm = TRUE),
              th_max = max(dsub$TH_scaled, na.rm = TRUE),
              metric = metric,
              Change = as.character(ch),
              Sex = as.character(sx),
              Phase = as.character(ph),
              data_subset = dsub,
              model_label = "lmer"
            )
            if (nrow(ind_auc_tbl) > 0) {
              spline_auc_individual_rows[[length(spline_auc_individual_rows) + 1]] <- ind_auc_tbl
            }
          }

          t_grid <- data.frame(
            TH_scaled = seq(min(dsub$TH_scaled), max(dsub$TH_scaled), length.out = 50),
            AnimalNum = dsub$AnimalNum[1]
          )
          pred_list <- lapply(unique(dsub$Group), function(g) {
            td <- t_grid; td$Group <- g
            p  <- tryCatch(predict(m_sp, newdata = td, re.form = NA), error = function(e) rep(NA_real_, nrow(td)))
            data.frame(td, predicted = p, metric = metric,
                       Change = as.character(ch), Sex = as.character(sx),
                       Phase = as.character(ph), df_spline = df_sp)
          })

          spline_results[[length(spline_results) + 1]] <- list(
            metric     = metric,
            Change     = as.character(ch),
            Sex        = as.character(sx),
            Phase      = as.character(ph),
            AIC        = AIC(m_sp),
            BIC        = BIC(m_sp),
            df_spline  = df_sp,
            predictions = dplyr::bind_rows(pred_list),
            model      = m_sp,                              # Store model for fixed effects
            data_subset = dsub,                             # Store data for metadata
            model_type = "lmer"                             # Model type identifier
          )
        }
      }
    }
  }

  if (length(spline_results) > 0) {
    spline_summary_tbl <- dplyr::bind_rows(lapply(spline_results, function(r) {
      tibble::tibble(metric = r$metric, Change = r$Change, Sex = r$Sex, Phase = r$Phase,
                     AIC = r$AIC, BIC = r$BIC, df_spline = r$df_spline)
    }))
    write_pub_xlsx(
      spline_summary_tbl,
      file.path(nat_dirs$splines_lmer_tables, paste0("model_summary_all.xlsx")),
      sheet_name = "spline_model_summary"
    )

    # Extract and export fixed effects with metadata
    spline_fixed_effects_all <- list()
    spline_metadata_all <- list()
    for (i in seq_along(spline_results)) {
      r <- spline_results[[i]]
      if (!is.null(r$model) && !is.null(r$data_subset)) {
        fe_result <- extract_fixed_effects_with_metadata(
          model_obj = r$model,
          model_type = r$model_type,
          metric = r$metric,
          Change = r$Change,
          Sex = r$Sex,
          Phase = r$Phase,
          data_subset = r$data_subset,
          group_levels = sort(unique(as.character(r$data_subset$Group)))
        )
        spline_fixed_effects_all[[i]] <- fe_result$fixed_effects
        # Store metadata with context
        spline_metadata_all[[i]] <- fe_result$metadata %>%
          dplyr::mutate(
            metric_var = r$metric,
            Change_var = r$Change,
            Sex_var = r$Sex,
            Phase_var = r$Phase
          )
      }
    }

    saveRDS(
      list(
        spline_results = spline_results,
        spline_auc_boot = if (exists("spline_auc_boot_tbl")) spline_auc_boot_tbl else tibble::tibble(),
        spline_auc_individual = if (length(spline_auc_individual_rows) > 0) dplyr::bind_rows(spline_auc_individual_rows) else tibble::tibble(),
        ar1_spline_results = if (exists("ar1_spline_results")) ar1_spline_results else list(),
        ar1_spline_auc_boot = if (exists("ar1_spline_auc_boot_tbl")) ar1_spline_auc_boot_tbl else tibble::tibble()
      ),
      file.path(nat_dirs$artifacts, paste0("spline_snapshot_", run_scope_current, ".rds")),
      compress = "xz"
    )

    if (length(spline_fixed_effects_all) > 0) {
      spline_fixef_tbl <- dplyr::bind_rows(spline_fixed_effects_all)
      spline_metadata_tbl <- dplyr::bind_rows(spline_metadata_all)

      # Write fixed effects with metadata sheets
      wb_fe <- openxlsx::createWorkbook()
      FN_fe <- "Arial"
      hdr_fe <- openxlsx::createStyle(
        fontName = FN_fe, fontSize = 10, textDecoration = "bold",
        fgFill = "#1F4E79", fontColour = "#FFFFFF", halign = "CENTER", valign = "CENTER",
        border = "Bottom", borderColour = "#173A5E", borderStyle = "medium"
      )
      s_fe_txt <- openxlsx::createStyle(fontName = FN_fe, fontSize = 9, halign = "LEFT",  valign = "CENTER")
      s_fe_num <- openxlsx::createStyle(fontName = FN_fe, fontSize = 9, halign = "RIGHT", valign = "CENTER", numFmt = "0.0000")
      s_fe_odd <- openxlsx::createStyle(fgFill = "#FFFFFF")
      s_fe_evn <- openxlsx::createStyle(fgFill = "#F7F7F7")
      s_fe_row_sig <- openxlsx::createStyle(fgFill = "#EAF6EC")
      s_fe_row_trd <- openxlsx::createStyle(fgFill = "#FFFDEB")
      apply_sheet <- function(wb, sht, tbl, p_col_name = "p.value") {
        openxlsx::addWorksheet(wb, sht, gridLines = TRUE)
        openxlsx::writeData(wb, sht, tbl, headerStyle = hdr_fe, rowNames = FALSE)
        openxlsx::setRowHeights(wb, sht, rows = 1, heights = 20)
        nr <- nrow(tbl); nc <- ncol(tbl)
        for (i in seq_len(nr)) {
          openxlsx::addStyle(wb, sht, style = if (i %% 2 == 1L) s_fe_odd else s_fe_evn,
                             rows = i + 1L, cols = 1:nc, gridExpand = TRUE, stack = FALSE)
        }
        if (p_col_name %in% names(tbl)) {
          for (i in seq_len(nr)) {
            pv_row <- suppressWarnings(as.numeric(tbl[[p_col_name]][i]))
            if (!is.finite(pv_row)) next
            s_row <- if (pv_row < 0.05) s_fe_row_sig else if (pv_row < 0.10) s_fe_row_trd else NULL
            if (!is.null(s_row)) {
              openxlsx::addStyle(wb, sht, style = s_row,
                                 rows = i + 1L, cols = 1:nc, gridExpand = TRUE, stack = TRUE)
            }
          }
        }
        for (j in seq_len(nc)) {
          s_a <- if (is.numeric(tbl[[j]])) s_fe_num else s_fe_txt
          openxlsx::addStyle(wb, sht, style = s_a,
                             rows = 2:(nr + 1), cols = j, gridExpand = FALSE, stack = TRUE)
        }
        # no cell-level p highlighting; row-level tint is the only significance styling
        openxlsx::freezePane(wb, sht, firstRow = TRUE)
        openxlsx::setColWidths(wb, sht, cols = 1:nc, widths = "auto")
        # Significance legend
        sig_r <- nr + 3L
        openxlsx::writeData(wb, sht,
          data.frame(x = "*** p\u202f<\u202f0.001   ** p\u202f<\u202f0.01   * p\u202f<\u202f0.05   \u2020 p\u202f<\u202f0.10   ns = not significant | row tint: green p<0.05, yellow p<0.10"),
          startRow = sig_r, startCol = 1, colNames = FALSE)
        openxlsx::addStyle(wb, sht,
          style = openxlsx::createStyle(fontName = FN_fe, fontSize = 8,
                                         fontColour = "#666666", textDecoration = "italic"),
          rows = sig_r, cols = 1, stack = FALSE)
        if (nc >= 2) openxlsx::mergeCells(wb, sht, cols = 1:min(nc, 10), rows = sig_r)
        invisible(wb)
      }

      # Fixed effects sheet
      apply_sheet(wb_fe, "Fixed Effects", spline_fixef_tbl)
      # Metadata/Interpretation sheet (no p-col)
      apply_sheet(wb_fe, "Metadata & Interpretation", spline_metadata_tbl, p_col_name = "NONE")
      openxlsx::setColWidths(wb_fe, "Metadata & Interpretation", cols = 1:2, widths = c(35, 65))
      
      openxlsx::saveWorkbook(
        wb_fe,
        file.path(nat_dirs$splines_lmer_tables, paste0("fixed_effects_all.xlsx")),
        overwrite = TRUE
      )
      cat("     ✓ Fixed effects with metadata exported\n")
    }

    spline_preds <- dplyr::bind_rows(lapply(spline_results, function(r) r$predictions))

    spline_auc_tbl <- spline_preds %>%
      dplyr::filter(!is.na(predicted), !is.na(TH_scaled), !is.na(Group)) %>%
      dplyr::group_by(metric, Change, Sex, Phase, Group, df_spline) %>%
      dplyr::summarise(
        th_min = min(TH_scaled, na.rm = TRUE),
        th_max = max(TH_scaled, na.rm = TRUE),
        th_range = th_max - th_min,
        AUC = trapz_auc(TH_scaled, predicted),
        AUC_norm = dplyr::if_else(th_range > 0, AUC / th_range, NA_real_),
        .groups = "drop"
      ) %>%
      dplyr::mutate(Group = factor(Group, levels = c("CON", "RES", "SUS")))

    if (nrow(spline_auc_tbl) > 0) {
      write_pub_xlsx(
        spline_auc_tbl,
        file.path(nat_dirs$splines_lmer_tables, paste0("auc_summary_all.xlsx")),
        sheet_name = "spline_auc"
      )

      auc_counts <- spline_auc_tbl %>%
        dplyr::count(metric, Sex, Phase, Group, name = "n")
      single_point_mode <- nrow(auc_counts) > 0 && max(auc_counts$n, na.rm = TRUE) <= 1

      p_auc <- ggplot(
        spline_auc_tbl,
        aes(x = Group, y = AUC, color = Group, fill = Group)
      ) +
        scale_color_manual(values = group_colors) +
        scale_fill_manual(values = group_colors) +
        labs(
          title = "Spline AUC by group",
          subtitle = ifelse(
            single_point_mode,
            "AUC from model-predicted trajectories (single estimate per group; point mode)",
            "AUC from model-predicted trajectories over TH_scaled range"
          ),
          x = "Group",
          y = "AUC (predicted response × TH_scaled)"
        ) +
        theme_classic(base_size = 9) +
        theme(
          legend.position = "none",
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)
        )

      if (single_point_mode) {
        p_auc <- p_auc +
          geom_line(aes(group = 1), color = "grey70", linewidth = 0.35, alpha = 0.75) +
          geom_point(size = 2.4, alpha = 0.95)
      } else {
        p_auc <- p_auc +
          geom_violin(alpha = 0.28, trim = TRUE, linewidth = 0.3) +
          geom_jitter(width = 0.12, height = 0, size = 1.8, alpha = 0.75) +
          stat_summary(fun = mean, geom = "crossbar", width = 0.35, linewidth = 0.5, color = "black")
      }

      if (includeSex && dplyr::n_distinct(spline_auc_tbl$Sex) > 1 &&
          includePhase && dplyr::n_distinct(spline_auc_tbl$Phase) > 1) {
        p_auc <- p_auc + facet_grid(metric + Phase ~ Sex, scales = "free_y")
      } else if (includeSex && dplyr::n_distinct(spline_auc_tbl$Sex) > 1) {
        p_auc <- p_auc + facet_grid(metric ~ Sex, scales = "free_y")
      } else if (includePhase && dplyr::n_distinct(spline_auc_tbl$Phase) > 1) {
        p_auc <- p_auc + facet_grid(metric ~ Phase, scales = "free_y")
      } else {
        p_auc <- p_auc + facet_wrap(~ metric, scales = "free_y")
      }

      ggsave(
        file.path(nat_dirs$splines_lmer_plots, paste0("auc_all.svg")),
        p_auc,
        width = 9,
        height = 6
      )

      if (RUN_SPLINE_AUC_BOOTSTRAP && length(spline_auc_boot_rows) > 0) {
        spline_auc_boot_tbl <- dplyr::bind_rows(spline_auc_boot_rows) %>%
          dplyr::mutate(p_boot_BH = p.adjust(p_boot, method = "BH"))

        write_pub_xlsx(
          spline_auc_boot_tbl,
          file.path(nat_dirs$splines_lmer_tables, paste0("auc_bootstrap_contrasts_all.xlsx")),
          sheet_name = "spline_auc_boot",
          p_cols = c("p_boot", "p_boot_BH")
        )

        # NEW: Save individual animal AUC values for proteomics correlation
        if (length(spline_auc_individual_rows) > 0) {
          spline_auc_individual_tbl <- dplyr::bind_rows(spline_auc_individual_rows) %>%
            dplyr::arrange(metric, Change, Sex, Phase, Group, AnimalNum)

          write_csv(
            spline_auc_individual_tbl,
            file.path(nat_dirs$splines_lmer_tables, paste0("auc_individual_animals_all.csv"))
          )

          cat("     [INFO] Saved individual animal AUC values to: auc_individual_animals_all.csv\n")
          cat("           Use x-axis: AUC values | y-axis: proteomics measurements (e.g., hippocampal oxphos log2)\n")
        }
      } else {
        spline_auc_boot_tbl <- tibble::tibble()
      }
    }

    # AR1-spline variant for all combinations (sensitivity check)
    ar1_spline_results <- list()
    ar1_spline_auc_boot_rows <- list()
    ctrl_ar1_sp <- nlme::lmeControl(
      opt = "optim", msMaxIter = 200, msMaxEval = 400,
      niterEM = 40, returnObject = TRUE
    )

    for (metric in metrics_to_analyze) {
      for (ch in changes_sp) {
        for (sx in sexes_sp) {
          for (ph in phases_sp) {
            dsub_ar1 <- data_filtered_agg %>%
              dplyr::filter(
                (!includeChange | Change == ch),
                (!includeSex    | Sex    == sx),
                (!includePhase  | Phase  == ph),
                !is.na(.data[[metric]]), !is.na(AnimalNum), !is.na(Group)
              )
            if (nrow(dsub_ar1) < max(min_obs, 15) ||
                dplyr::n_distinct(dsub_ar1$AnimalNum) < min_animals ||
                dplyr::n_distinct(dsub_ar1$Group) < min_groups) next

            t0_ar1 <- min(dsub_ar1$HalfHourElapsed, na.rm = TRUE)
            dsub_ar1 <- dsub_ar1 %>%
              dplyr::mutate(
                TH        = (HalfHourElapsed - t0_ar1) / 2,
                TH_scaled = as.vector(scale(TH))
              ) %>%
              dplyr::filter(!is.na(TH_scaled))

            n_uniq_t_ar1 <- dplyr::n_distinct(dsub_ar1$TH_scaled)
            if (n_uniq_t_ar1 < 6) next

            df_sp_ar1 <- min(4L, floor(n_uniq_t_ar1 / 3L))
            f_sp_ar1 <- as.formula(paste0(
              metric, " ~ Group * ns(TH_scaled, df = ", df_sp_ar1, ")"
            ))

            m_sp_ar1 <- tryCatch(
              nlme::lme(
                fixed = f_sp_ar1,
                random = ~ 1 | AnimalNum,
                correlation = nlme::corAR1(form = ~ TH_scaled | AnimalNum),
                data = dsub_ar1,
                control = ctrl_ar1_sp,
                na.action = na.omit
              ),
              error = function(e) NULL
            )
            if (is.null(m_sp_ar1)) next

            if (RUN_SPLINE_AUC_BOOTSTRAP) {
              boot_tbl_ar1 <- bootstrap_auc_contrasts(
                model_obj = m_sp_ar1,
                th_min = min(dsub_ar1$TH_scaled, na.rm = TRUE),
                th_max = max(dsub_ar1$TH_scaled, na.rm = TRUE),
                group_levels = unique(as.character(dsub_ar1$Group)),
                metric = metric,
                Change = as.character(ch),
                Sex = as.character(sx),
                Phase = as.character(ph),
                model_label = "ar1",
                n_boot = N_SPLINE_AUC_BOOT
              )
              if (nrow(boot_tbl_ar1) > 0) {
                ar1_spline_auc_boot_rows[[length(ar1_spline_auc_boot_rows) + 1]] <- boot_tbl_ar1
              }

              # NEW: Extract individual animal AUC values for AR1 model
              ind_auc_tbl_ar1 <- extract_individual_animal_auc(
                model_obj = m_sp_ar1,
                th_min = min(dsub_ar1$TH_scaled, na.rm = TRUE),
                th_max = max(dsub_ar1$TH_scaled, na.rm = TRUE),
                metric = metric,
                Change = as.character(ch),
                Sex = as.character(sx),
                Phase = as.character(ph),
                data_subset = dsub_ar1,
                model_label = "ar1"
              )
              if (nrow(ind_auc_tbl_ar1) > 0) {
                spline_auc_individual_rows[[length(spline_auc_individual_rows) + 1]] <- ind_auc_tbl_ar1
              }
            }

            # Predictions on grid
            t_grid_ar1 <- data.frame(
              TH_scaled = seq(min(dsub_ar1$TH_scaled), max(dsub_ar1$TH_scaled), length.out = 50),
              AnimalNum = dsub_ar1$AnimalNum[1]
            )

            pred_list_ar1 <- lapply(unique(dsub_ar1$Group), function(g) {
              td <- t_grid_ar1; td$Group <- g
              p <- tryCatch(
                predict(m_sp_ar1, newdata = td, re.form = NA),
                error = function(e) rep(NA_real_, nrow(td))
              )
              data.frame(
                td, predicted = p, metric = metric,
                Change = as.character(ch), Sex = as.character(sx),
                Phase = as.character(ph), df_spline = df_sp_ar1
              )
            })

            ar1_spline_results[[length(ar1_spline_results) + 1]] <- list(
              metric     = metric,
              Change     = as.character(ch),
              Sex        = as.character(sx),
              Phase      = as.character(ph),
              AIC        = tryCatch(AIC(m_sp_ar1), error = function(e) NA_real_),
              BIC        = tryCatch(BIC(m_sp_ar1), error = function(e) NA_real_),
              df_spline  = df_sp_ar1,
              predictions = dplyr::bind_rows(pred_list_ar1),
              model      = m_sp_ar1,                        # Store model for fixed effects
              data_subset = dsub_ar1,                       # Store data for metadata
              model_type = "ar1"                            # Model type identifier
            )
          }
        }
      }
    }

    if (length(ar1_spline_results) > 0) {
      # AR1 spline model summary table
      ar1_spline_summary_tbl <- dplyr::bind_rows(lapply(ar1_spline_results, function(r) {
        tibble::tibble(metric = r$metric, Change = r$Change, Sex = r$Sex, Phase = r$Phase,
                       AIC = r$AIC, BIC = r$BIC, df_spline = r$df_spline)
      }))
      write_pub_xlsx(
        ar1_spline_summary_tbl,
        file.path(nat_dirs$splines_ar1_tables, paste0("model_summary_all.xlsx")),
        sheet_name = "spline_model_ar1"
      )

      # Extract and export fixed effects with metadata for AR1 models
      ar1_fixed_effects_all <- list()
      ar1_metadata_all <- list()
      for (i in seq_along(ar1_spline_results)) {
        r <- ar1_spline_results[[i]]
        if (!is.null(r$model) && !is.null(r$data_subset)) {
          fe_result <- extract_fixed_effects_with_metadata(
            model_obj = r$model,
            model_type = r$model_type,
            metric = r$metric,
            Change = r$Change,
            Sex = r$Sex,
            Phase = r$Phase,
            data_subset = r$data_subset,
            group_levels = sort(unique(as.character(r$data_subset$Group)))
          )
          ar1_fixed_effects_all[[i]] <- fe_result$fixed_effects
          ar1_metadata_all[[i]] <- fe_result$metadata %>%
            dplyr::mutate(
              metric_var = r$metric,
              Change_var = r$Change,
              Sex_var = r$Sex,
              Phase_var = r$Phase
            )
        }
      }

      if (length(ar1_fixed_effects_all) > 0) {
        ar1_fixef_tbl <- dplyr::bind_rows(ar1_fixed_effects_all)
        ar1_metadata_tbl <- dplyr::bind_rows(ar1_metadata_all)

        wb_fe_ar1 <- openxlsx::createWorkbook()
        apply_sheet(wb_fe_ar1, "Fixed Effects (AR1)", ar1_fixef_tbl)
        apply_sheet(wb_fe_ar1, "Metadata & Interpretation", ar1_metadata_tbl, p_col_name = "NONE")
        openxlsx::setColWidths(wb_fe_ar1, "Metadata & Interpretation", cols = 1:2, widths = c(35, 65))
        
        openxlsx::saveWorkbook(
          wb_fe_ar1,
          file.path(nat_dirs$splines_ar1_tables, paste0("fixed_effects_all.xlsx")),
          overwrite = TRUE
        )
        cat("     ✓ AR1 fixed effects with metadata exported\n")
      }

      ar1_spline_preds <- dplyr::bind_rows(lapply(ar1_spline_results, function(r) r$predictions))

      ar1_spline_auc_tbl <- ar1_spline_preds %>%
        dplyr::filter(!is.na(predicted), !is.na(TH_scaled), !is.na(Group)) %>%
        dplyr::group_by(metric, Change, Sex, Phase, Group, df_spline) %>%
        dplyr::summarise(
          th_min = min(TH_scaled, na.rm = TRUE),
          th_max = max(TH_scaled, na.rm = TRUE),
          th_range = th_max - th_min,
          AUC = trapz_auc(TH_scaled, predicted),
          AUC_norm = dplyr::if_else(th_range > 0, AUC / th_range, NA_real_),
          .groups = "drop"
        ) %>%
        dplyr::mutate(Group = factor(Group, levels = c("CON", "RES", "SUS")))

      if (nrow(ar1_spline_auc_tbl) > 0) {
        write_pub_xlsx(
          ar1_spline_auc_tbl,
          file.path(nat_dirs$splines_ar1_tables, paste0("auc_summary_all.xlsx")),
          sheet_name = "spline_auc_ar1"
        )

        ar1_auc_counts <- ar1_spline_auc_tbl %>%
          dplyr::count(metric, Sex, Phase, Group, name = "n")
        single_point_mode_ar1 <- nrow(ar1_auc_counts) > 0 && max(ar1_auc_counts$n, na.rm = TRUE) <= 1

        # AR1 AUC plot
        p_ar1_auc <- ggplot(
          ar1_spline_auc_tbl,
          aes(x = Group, y = AUC, color = Group, fill = Group)
        ) +
          scale_color_manual(values = group_colors) +
          scale_fill_manual(values = group_colors) +
          labs(
            title = "AR1-spline AUC by group",
            subtitle = ifelse(
              single_point_mode_ar1,
              "AUC from AR1 model-predicted trajectories (single estimate per group; point mode)",
              "AUC from AR1-nlme model-predicted trajectories (sensitivity check)"
            ),
            x = "Group",
            y = "AUC (predicted response × TH_scaled)"
          ) +
          theme_classic(base_size = 9) +
          theme(
            legend.position = "none",
            strip.text = element_text(face = "bold"),
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, size = 8)
          )

        if (single_point_mode_ar1) {
          p_ar1_auc <- p_ar1_auc +
            geom_line(aes(group = 1), color = "grey70", linewidth = 0.35, alpha = 0.75) +
            geom_point(size = 2.4, alpha = 0.95)
        } else {
          p_ar1_auc <- p_ar1_auc +
            geom_violin(alpha = 0.28, trim = TRUE, linewidth = 0.3) +
            geom_jitter(width = 0.12, height = 0, size = 1.8, alpha = 0.75) +
            stat_summary(fun = mean, geom = "crossbar", width = 0.35, linewidth = 0.5, color = "black")
        }

        if (includeSex && dplyr::n_distinct(ar1_spline_auc_tbl$Sex) > 1 &&
            includePhase && dplyr::n_distinct(ar1_spline_auc_tbl$Phase) > 1) {
          p_ar1_auc <- p_ar1_auc + facet_grid(metric + Phase ~ Sex, scales = "free_y")
        } else if (includeSex && dplyr::n_distinct(ar1_spline_auc_tbl$Sex) > 1) {
          p_ar1_auc <- p_ar1_auc + facet_grid(metric ~ Sex, scales = "free_y")
        } else if (includePhase && dplyr::n_distinct(ar1_spline_auc_tbl$Phase) > 1) {
          p_ar1_auc <- p_ar1_auc + facet_grid(metric ~ Phase, scales = "free_y")
        } else {
          p_ar1_auc <- p_ar1_auc + facet_wrap(~ metric, scales = "free_y")
        }

        ggsave(
          file.path(nat_dirs$splines_ar1_plots, paste0("auc_all.svg")),
          p_ar1_auc,
          width = 9,
          height = 6
        )

        # AR1 trajectory plot
        ar1_spline_preds_plot <- ar1_spline_preds %>%
          dplyr::filter(!is.na(predicted)) %>%
          dplyr::mutate(Group = factor(Group, levels = c("CON", "RES", "SUS")))

        p_ar1_sp <- ggplot(ar1_spline_preds_plot, aes(x = TH_scaled, y = predicted, color = Group, group = Group)) +
          geom_line(linewidth = 1.0, alpha = 0.85) +
          scale_color_manual(values = group_colors) +
          labs(
            title = "AR1-spline trajectories (all combinations)",
            subtitle = "nlme::lme + corAR1 | ns(TH_scaled)",
            x = "Time (scaled)",
            y = "Predicted response",
            color = "Group"
          ) +
          theme_classic(base_size = 9) +
          theme(
            legend.position = "top",
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, size = 8),
            strip.text = element_text(face = "bold")
          )

        if (includeSex && dplyr::n_distinct(ar1_spline_preds_plot$Sex) > 1 &&
            includePhase && dplyr::n_distinct(ar1_spline_preds_plot$Phase) > 1) {
          p_ar1_sp <- p_ar1_sp + facet_grid(metric + Phase ~ Sex, scales = "free_y")
        } else if (includeSex && dplyr::n_distinct(ar1_spline_preds_plot$Sex) > 1) {
          p_ar1_sp <- p_ar1_sp + facet_grid(metric ~ Sex, scales = "free_y")
        } else if (includePhase && dplyr::n_distinct(ar1_spline_preds_plot$Phase) > 1) {
          p_ar1_sp <- p_ar1_sp + facet_grid(metric ~ Phase, scales = "free_y")
        } else {
          p_ar1_sp <- p_ar1_sp + facet_wrap(~ metric, scales = "free_y")
        }

        ggsave(
          file.path(nat_dirs$splines_ar1_plots, paste0("trajectory_all.svg")),
          p_ar1_sp,
          width = 9,
          height = 6
        )

        if (RUN_SPLINE_AUC_BOOTSTRAP && length(ar1_spline_auc_boot_rows) > 0) {
          ar1_spline_auc_boot_tbl <- dplyr::bind_rows(ar1_spline_auc_boot_rows) %>%
            dplyr::mutate(p_boot_BH = p.adjust(p_boot, method = "BH"))

          write_pub_xlsx(
            ar1_spline_auc_boot_tbl,
            file.path(nat_dirs$splines_ar1_tables, paste0("auc_bootstrap_contrasts_all.xlsx")),
            sheet_name = "spline_auc_boot_ar1",
            p_cols = c("p_boot", "p_boot_BH")
          )
        } else {
          ar1_spline_auc_boot_tbl <- tibble::tibble()
        }

        cat("     \u2713 Full AR1-spline results saved\n")
      }
    }

    # Additional focused export: first Change + first Active phase only
    change_levels_all <- levels(droplevels(data_filtered_agg$Change))
    if (length(change_levels_all) == 0) {
      change_levels_all <- sort(unique(as.character(data_filtered_agg$Change)))
    }
    change_ord_num <- suppressWarnings(as.numeric(stringr::str_extract(change_levels_all, "\\d+")))
    first_change_level <- if (length(change_levels_all) == 0) {
      NA_character_
    } else if (all(!is.na(change_ord_num))) {
      change_levels_all[order(change_ord_num)][1]
    } else {
      sort(change_levels_all)[1]
    }

    if (!is.na(first_change_level) && exists("spline_auc_boot_tbl") && nrow(spline_auc_boot_tbl) > 0) {
      spline_auc_boot_first <- spline_auc_boot_tbl %>%
        dplyr::filter(as.character(Change) == first_change_level, as.character(Phase) == "Active")
      if (nrow(spline_auc_boot_first) > 0) {
        write_pub_xlsx(
          spline_auc_boot_first,
          file.path(nat_dirs$splines_lmer_tables, paste0("auc_bootstrap_firstChangeActive.xlsx")),
          sheet_name = "spline_auc_boot_C1A",
          p_cols = c("p_boot", "p_boot_BH")
        )
      }
    }

    if (!is.na(first_change_level) && exists("ar1_spline_auc_boot_tbl") && nrow(ar1_spline_auc_boot_tbl) > 0) {
      ar1_spline_auc_boot_first <- ar1_spline_auc_boot_tbl %>%
        dplyr::filter(as.character(Change) == first_change_level, as.character(Phase) == "Active")
      if (nrow(ar1_spline_auc_boot_first) > 0) {
        write_pub_xlsx(
          ar1_spline_auc_boot_first,
          file.path(nat_dirs$splines_ar1_tables, paste0("auc_bootstrap_firstChangeActive.xlsx")),
          sheet_name = "spline_auc_boot_ar1_C1A",
          p_cols = c("p_boot", "p_boot_BH")
        )
      }
    }

    sp_first_active <- spline_preds %>%
      dplyr::filter(
        as.character(Change) == first_change_level,
        as.character(Phase) == "Active",
        !is.na(predicted)
      ) %>%
      dplyr::mutate(Group = factor(Group, levels = c("CON", "RES", "SUS")))

    if (!is.na(first_change_level) && nrow(sp_first_active) > 0) {
      # Focused model summary table (first change, Active)
      spline_summary_first_active <- spline_summary_tbl %>%
        dplyr::filter(as.character(Change) == first_change_level, as.character(Phase) == "Active")

      if (nrow(spline_summary_first_active) > 0) {
        write_pub_xlsx(
          spline_summary_first_active,
          file.path(nat_dirs$splines_lmer_tables, paste0("model_summary_firstChangeActive.xlsx")),
          sheet_name = "spline_model_C1A"
        )
      }

      # Focused AUC table (first change, Active)
      spline_auc_first_active <- sp_first_active %>%
        dplyr::group_by(metric, Change, Sex, Phase, Group, df_spline) %>%
        dplyr::summarise(
          th_min = min(TH_scaled, na.rm = TRUE),
          th_max = max(TH_scaled, na.rm = TRUE),
          th_range = th_max - th_min,
          AUC = trapz_auc(TH_scaled, predicted),
          AUC_norm = dplyr::if_else(th_range > 0, AUC / th_range, NA_real_),
          .groups = "drop"
        )

      if (nrow(spline_auc_first_active) > 0) {
        write_pub_xlsx(
          spline_auc_first_active,
          file.path(nat_dirs$splines_lmer_tables, paste0("auc_firstChangeActive.xlsx")),
          sheet_name = "spline_auc_C1A"
        )

        auc_first_counts <- spline_auc_first_active %>%
          dplyr::count(metric, Sex, Group, name = "n")
        single_point_mode_first <- nrow(auc_first_counts) > 0 && max(auc_first_counts$n, na.rm = TRUE) <= 1

        p_auc_first <- ggplot(
          spline_auc_first_active,
          aes(x = Group, y = AUC, color = Group, fill = Group)
        ) +
          scale_color_manual(values = group_colors) +
          scale_fill_manual(values = group_colors) +
          labs(
            title = paste0("Spline AUC (", first_change_level, ", Active phase)"),
            subtitle = ifelse(
              single_point_mode_first,
              "AUC from model-predicted trajectories (single estimate per group; point mode)",
              "AUC from model-predicted trajectories"
            ),
            x = "Group",
            y = "AUC (predicted response × TH_scaled)"
          ) +
          theme_classic(base_size = 9) +
          theme(
            legend.position = "none",
            strip.text = element_text(face = "bold"),
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5)
          )

        if (single_point_mode_first) {
          p_auc_first <- p_auc_first +
            geom_line(aes(group = 1), color = "grey70", linewidth = 0.35, alpha = 0.75) +
            geom_point(size = 2.4, alpha = 0.95)
        } else {
          p_auc_first <- p_auc_first +
            geom_violin(alpha = 0.30, trim = TRUE, linewidth = 0.3) +
            geom_jitter(width = 0.12, size = 1.8, alpha = 0.75) +
            stat_summary(fun = mean, geom = "crossbar", width = 0.35, linewidth = 0.5, color = "black")
        }

        if (includeSex && dplyr::n_distinct(spline_auc_first_active$Sex) > 1) {
          p_auc_first <- p_auc_first + facet_grid(metric ~ Sex, scales = "free_y")
        } else {
          p_auc_first <- p_auc_first + facet_wrap(~ metric, scales = "free_y")
        }

        ggsave(
          file.path(nat_dirs$splines_lmer_plots, paste0("auc_firstChangeActive.svg")),
          p_auc_first,
          width = 8,
          height = 5
        )
      }

      # Focused trajectory plot (first change, Active)
      p_sp_first <- ggplot(sp_first_active, aes(x = TH_scaled, y = predicted, color = Group, group = Group)) +
        geom_line(linewidth = 1.0) +
        scale_color_manual(values = group_colors) +
        labs(
          title = paste0("Spline trajectories (", first_change_level, ", Active phase)"),
          subtitle = paste0("ns(TH, df up to ", max(sp_first_active$df_spline, na.rm = TRUE), ")"),
          x = "Time (scaled)",
          y = "Predicted response",
          color = "Group"
        ) +
        theme_classic(base_size = 9) +
        theme(
          legend.position = "top",
          plot.title = element_text(face = "bold", hjust = 0.5),
          strip.text = element_text(face = "bold")
        )

      if (includeSex && dplyr::n_distinct(sp_first_active$Sex) > 1) {
        p_sp_first <- p_sp_first + facet_grid(metric ~ Sex, scales = "free_y")
      } else {
        p_sp_first <- p_sp_first + facet_wrap(~ metric, scales = "free_y")
      }

      ggsave(
        file.path(nat_dirs$splines_lmer_plots, paste0("trajectory_firstChangeActive.svg")),
        p_sp_first,
        width = 8,
        height = 5
      )
    } else {
      cat("     - Skipped focused first-change active spline (no data)\n")
    }

    # AR1-spline variant for first change + Active phase (sensitivity check)
    if (!is.na(first_change_level)) {
      cat("     [AR1-spline variant] Fitting first-change active AR1 splines...\n")

      ar1_spline_results_focused <- list()
      ctrl_ar1_sp <- nlme::lmeControl(
        opt = "optim", msMaxIter = 200, msMaxEval = 400,
        niterEM = 40, returnObject = TRUE
      )

      for (metric in metrics_to_analyze) {
        for (sx in sexes_sp) {
          dsub_ar1 <- data_filtered_agg %>%
            dplyr::filter(
              as.character(Change) == first_change_level,
              as.character(Phase) == "Active",
              !includeSex | Sex == sx,
              !is.na(.data[[metric]]), !is.na(AnimalNum), !is.na(Group)
            )
          if (nrow(dsub_ar1) < max(min_obs, 15) ||
              dplyr::n_distinct(dsub_ar1$AnimalNum) < min_animals ||
              dplyr::n_distinct(dsub_ar1$Group) < min_groups) next

          t0_ar1 <- min(dsub_ar1$HalfHourElapsed, na.rm = TRUE)
          dsub_ar1 <- dsub_ar1 %>%
            dplyr::mutate(
              TH = (HalfHourElapsed - t0_ar1) / 2,
              TH_scaled = as.vector(scale(TH)),
              AnimalNum = factor(AnimalNum)
            ) %>%
            dplyr::filter(!is.na(TH_scaled))

          n_uniq_t_ar1 <- dplyr::n_distinct(dsub_ar1$TH_scaled)
          if (n_uniq_t_ar1 < 6) next

          df_sp_ar1 <- min(4L, floor(n_uniq_t_ar1 / 3L))
          f_sp_ar1 <- as.formula(paste0(
            metric, " ~ Group * ns(TH_scaled, df = ", df_sp_ar1, ")"
          ))

          m_sp_ar1 <- tryCatch(
            nlme::lme(
              fixed = f_sp_ar1,
              random = ~ 1 | AnimalNum,
              correlation = nlme::corAR1(form = ~ TH_scaled | AnimalNum),
              data = dsub_ar1,
              control = ctrl_ar1_sp,
              na.action = na.omit
            ),
            error = function(e) NULL
          )
          if (is.null(m_sp_ar1)) next

          # Predictions on grid
          t_grid_ar1 <- data.frame(
            TH_scaled = seq(min(dsub_ar1$TH_scaled), max(dsub_ar1$TH_scaled), length.out = 50),
            AnimalNum = dsub_ar1$AnimalNum[1]
          )

          pred_list_ar1 <- lapply(unique(dsub_ar1$Group), function(g) {
            td <- t_grid_ar1; td$Group <- g
            p <- tryCatch(
              predict(m_sp_ar1, newdata = td, re.form = NA),
              error = function(e) rep(NA_real_, nrow(td))
            )
            data.frame(
              td, predicted = p, metric = metric, Sex = as.character(sx),
              df_spline = df_sp_ar1, model = "ar1"
            )
          })

          ar1_spline_results_focused[[length(ar1_spline_results_focused) + 1]] <- list(
            metric     = metric,
            Sex        = as.character(sx),
            AIC        = tryCatch(AIC(m_sp_ar1), error = function(e) NA_real_),
            BIC        = tryCatch(BIC(m_sp_ar1), error = function(e) NA_real_),
            df_spline  = df_sp_ar1,
            predictions = dplyr::bind_rows(pred_list_ar1)
          )
        }
      }

      if (length(ar1_spline_results_focused) > 0) {
        # AR1 spline model summary table
        ar1_sp_summary <- dplyr::bind_rows(lapply(ar1_spline_results_focused, function(r) {
          tibble::tibble(
            metric = r$metric, Sex = r$Sex, model = "ar1",
            AIC = r$AIC, BIC = r$BIC, df_spline = r$df_spline
          )
        }))

        write_pub_xlsx(
          ar1_sp_summary,
          file.path(nat_dirs$splines_ar1_tables, paste0("model_summary_firstChangeActive.xlsx")),
          sheet_name = "spline_ar1_C1A"
        )

        # AR1 spline predictions + AUC
        ar1_sp_preds <- dplyr::bind_rows(lapply(ar1_spline_results_focused, function(r) r$predictions))

        ar1_sp_auc <- ar1_sp_preds %>%
          dplyr::filter(!is.na(predicted), !is.na(TH_scaled), !is.na(Group)) %>%
          dplyr::group_by(metric, Sex, Group, df_spline, model) %>%
          dplyr::summarise(
            th_min = min(TH_scaled, na.rm = TRUE),
            th_max = max(TH_scaled, na.rm = TRUE),
            th_range = th_max - th_min,
            AUC = trapz_auc(TH_scaled, predicted),
            AUC_norm = dplyr::if_else(th_range > 0, AUC / th_range, NA_real_),
            .groups = "drop"
          )

        if (nrow(ar1_sp_auc) > 0) {
          write_pub_xlsx(
            ar1_sp_auc,
            file.path(nat_dirs$splines_ar1_tables, paste0("auc_firstChangeActive.xlsx")),
            sheet_name = "spline_auc_ar1_C1A"
          )

          ar1_first_counts <- ar1_sp_auc %>%
            dplyr::count(metric, Sex, Group, name = "n")
          single_point_mode_ar1_first <- nrow(ar1_first_counts) > 0 && max(ar1_first_counts$n, na.rm = TRUE) <= 1

          p_ar1_auc <- ggplot(
            ar1_sp_auc %>% dplyr::mutate(Group = factor(Group, levels = c("CON", "RES", "SUS"))),
            aes(x = Group, y = AUC, color = Group, fill = Group)
          ) +
            scale_color_manual(values = group_colors) +
            scale_fill_manual(values = group_colors) +
            labs(
              title = paste0("AR1-spline AUC (", first_change_level, ", Active)"),
              subtitle = ifelse(
                single_point_mode_ar1_first,
                "AR1-nlme model with corAR1 structure (single estimate per group; point mode)",
                "AR1-nlme model with corAR1 structure"
              ),
              x = "Group",
              y = "AUC"
            ) +
            theme_classic(base_size = 9) +
            theme(
              legend.position = "none",
              strip.text = element_text(face = "bold"),
              plot.title = element_text(face = "bold", hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5, size = 8, face = "italic")
            )

          if (single_point_mode_ar1_first) {
            p_ar1_auc <- p_ar1_auc +
              geom_line(aes(group = 1), color = "grey70", linewidth = 0.35, alpha = 0.75) +
              geom_point(size = 2.4, alpha = 0.95)
          } else {
            p_ar1_auc <- p_ar1_auc +
              geom_violin(alpha = 0.30, trim = TRUE, linewidth = 0.3) +
              geom_jitter(width = 0.12, size = 1.8, alpha = 0.75) +
              stat_summary(fun = mean, geom = "crossbar", width = 0.35, linewidth = 0.5, color = "black")
          }

          if (includeSex && dplyr::n_distinct(ar1_sp_auc$Sex) > 1) {
            p_ar1_auc <- p_ar1_auc + facet_grid(metric ~ Sex, scales = "free_y")
          } else {
            p_ar1_auc <- p_ar1_auc + facet_wrap(~ metric, scales = "free_y")
          }

          ggsave(
            file.path(nat_dirs$splines_ar1_plots, paste0("auc_firstChangeActive.svg")),
            p_ar1_auc,
            width = 8,
            height = 5
          )

          # AR1 trajectory plot
          ar1_sp_preds_plot <- ar1_sp_preds %>%
            dplyr::filter(!is.na(predicted)) %>%
            dplyr::mutate(Group = factor(Group, levels = c("CON", "RES", "SUS")))

          p_ar1_sp <- ggplot(ar1_sp_preds_plot, aes(x = TH_scaled, y = predicted, color = Group, group = Group)) +
            geom_line(linewidth = 1.0, alpha = 0.85) +
            scale_color_manual(values = group_colors) +
            labs(
              title = paste0("AR1-spline trajectories (", first_change_level, ", Active)"),
              subtitle = paste0("nlme::lme + corAR1 | ns(TH, df up to ", max(ar1_sp_preds_plot$df_spline, na.rm = TRUE), ")"),
              x = "Time (scaled)",
              y = "Predicted response",
              color = "Group"
            ) +
            theme_classic(base_size = 9) +
            theme(
              legend.position = "top",
              plot.title = element_text(face = "bold", hjust = 0.5),
              plot.subtitle = element_text(hjust = 0.5, size = 8, face = "italic"),
              strip.text = element_text(face = "bold")
            )

          if (includeSex && dplyr::n_distinct(ar1_sp_preds_plot$Sex) > 1) {
            p_ar1_sp <- p_ar1_sp + facet_grid(metric ~ Sex, scales = "free_y")
          } else {
            p_ar1_sp <- p_ar1_sp + facet_wrap(~ metric, scales = "free_y")
          }

          ggsave(
            file.path(nat_dirs$splines_ar1_plots, paste0("trajectory_firstChangeActive.svg")),
            p_ar1_sp,
            width = 8,
            height = 5
          )
          cat("     \u2713 AR1-spline saved\n")
        }
      } else {
        cat("     - AR1-spline skipped (no convergence)\n")
      }
    }

    for (metric in metrics_to_analyze) {
      sp_sub <- spline_preds %>%
        dplyr::filter(metric == !!metric, !is.na(predicted)) %>%
        dplyr::mutate(Group = factor(Group, levels = c("CON", "RES", "SUS")))
      if (nrow(sp_sub) == 0) next

      p_sp <- ggplot(sp_sub, aes(x = TH_scaled, y = predicted, color = Group, group = Group)) +
        geom_line(linewidth = 1.0) +
        scale_color_manual(values = group_colors) +
        labs(
          title    = paste0(metric, ": Spline-smoothed group trajectories"),
          subtitle = paste0("ns(TH, df = ", max(sp_sub$df_spline, na.rm = TRUE), ")"),
          x        = "Time (scaled)", y = paste0("Predicted ", metric), color = "Group"
        ) +
        theme_classic(base_size = 9) +
        theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5),
              strip.text = element_text(face = "bold"))

      if (includeSex && dplyr::n_distinct(sp_sub$Sex) > 1 && includePhase && dplyr::n_distinct(sp_sub$Phase) > 1) {
        p_sp <- p_sp + facet_grid(Phase ~ Sex, scales = "free_y")
      } else if (includeSex && dplyr::n_distinct(sp_sub$Sex) > 1) {
        p_sp <- p_sp + facet_wrap(~ Sex, scales = "free_y")
      } else if (includePhase && dplyr::n_distinct(sp_sub$Phase) > 1) {
        p_sp <- p_sp + facet_wrap(~ Phase, scales = "free_y")
      }

      ggsave(
        file.path(nat_dirs$plots, paste0("spline_trajectory_", metric, "_", run_scope_current, ".svg")),
        p_sp, width = 7, height = 5
      )
    }
    cat("     \u2713 Saved\n")
  } else {
    cat("     - No spline models converged\n")
  }

  # ================================================================
  # #5: Individual heterogeneity (random-effect dispersion)
  # ================================================================
  cat("  [5/10] Individual heterogeneity analysis...\n")

  hetero_rows <- list()

  for (metric in metrics_to_analyze) {
    for (key in names(rep_models[[metric]])) {
      m_obj <- rep_models[[metric]][[key]]
      if (is.null(m_obj)) next

      re <- tryCatch(as.data.frame(lme4::ranef(m_obj$model)$AnimalNum), error = function(e) NULL)
      if (is.null(re) || nrow(re) == 0) next

      re$AnimalNum <- as.character(rownames(re))
      if ("(Intercept)" %in% names(re)) re <- dplyr::rename(re, intercept = `(Intercept)`)

      animal_meta <- m_obj$data %>%
        dplyr::select(AnimalNum, Group, Sex) %>%
        dplyr::distinct(AnimalNum, .keep_all = TRUE) %>%
        dplyr::mutate(AnimalNum = as.character(AnimalNum))

      re_joined <- dplyr::left_join(re, animal_meta, by = "AnimalNum") %>%
        dplyr::mutate(metric = metric, Change = m_obj$Change, Sex = m_obj$Sex, Phase = m_obj$Phase)

      hetero_rows[[length(hetero_rows) + 1]] <- re_joined
    }
  }

  if (length(hetero_rows) > 0) {
    hetero_tbl <- dplyr::bind_rows(hetero_rows)

    hetero_var <- hetero_tbl %>%
      dplyr::filter(!is.na(intercept), !is.na(Group)) %>%
      dplyr::group_by(metric, Change, Sex, Phase, Group) %>%
      dplyr::summarise(
        n_animals     = dplyr::n(),
        mean_intercept = mean(intercept, na.rm = TRUE),
        sd_intercept   = sd(intercept,   na.rm = TRUE),
        var_intercept  = var(intercept,  na.rm = TRUE),
        .groups = "drop"
      )

    write_pub_xlsx(hetero_tbl,
      file.path(nat_dirs$tables, paste0("individual_random_effects_", run_scope_current, ".xlsx")),
      sheet_name = "Results", p_cols = character(0))
    write_pub_xlsx(hetero_var,
      file.path(nat_dirs$tables, paste0("individual_heterogeneity_variance_", run_scope_current, ".xlsx")),
      sheet_name = "Results", p_cols = character(0))

    p_hetero <- ggplot(
      hetero_tbl %>%
        dplyr::filter(!is.na(intercept)) %>%
        dplyr::mutate(Group = factor(Group, levels = c("CON", "RES", "SUS"))),
      aes(x = Group, y = intercept, fill = Group, color = Group)
    ) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
      geom_violin(alpha = 0.45, trim = TRUE) +
      geom_jitter(width = 0.12, size = 1.5, alpha = 0.6) +
      stat_summary(fun = mean, geom = "crossbar", width = 0.35, linewidth = 0.5, color = "black") +
      scale_fill_manual(values  = group_colors) +
      scale_color_manual(values = group_colors) +
      facet_wrap(~ metric + Phase, scales = "free_y", ncol = 3) +
      theme_classic(base_size = 9) +
      theme(
        legend.position = "none",
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5)
      ) +
      labs(
        title = "Individual heterogeneity: random intercepts by group",
        x = "Group", y = "Random intercept (deviation from grand mean)"
      )

    ggsave(
      file.path(nat_dirs$plots, paste0("individual_heterogeneity_violin_", run_scope_current, ".svg")),
      p_hetero, width = 9, height = 6
    )
    cat("     \u2713 Saved\n")
  } else {
    cat("     - No random-effect data available\n")
  }

  # ================================================================
  # #6: Within-animal ICC / repeatability
  # ================================================================
  cat("  [6/10] ICC / repeatability...\n")

  icc_rows <- list()

  for (metric in metrics_to_analyze) {
    for (key in names(rep_models[[metric]])) {
      m_obj <- rep_models[[metric]][[key]]
      if (is.null(m_obj)) next

      vc <- tryCatch(as.data.frame(lme4::VarCorr(m_obj$model)), error = function(e) NULL)
      if (is.null(vc)) next

      var_animal   <- sum(vc$vcov[vc$grp != "Residual"], na.rm = TRUE)
      var_residual <- sum(vc$vcov[vc$grp == "Residual"], na.rm = TRUE)
      total_var    <- var_animal + var_residual
      icc_val      <- if (total_var > 0) var_animal / total_var else NA_real_

      icc_rows[[length(icc_rows) + 1]] <- tibble::tibble(
        metric       = metric,
        Change       = m_obj$Change,
        Sex          = m_obj$Sex,
        Phase        = m_obj$Phase,
        var_animal   = var_animal,
        var_residual = var_residual,
        total_var    = total_var,
        ICC          = icc_val
      )
    }
  }

  if (length(icc_rows) > 0) {
    icc_tbl <- dplyr::bind_rows(icc_rows)
    write_pub_xlsx(icc_tbl,
      file.path(nat_dirs$tables, paste0("icc_reliability_", run_scope_current, ".xlsx")),
      sheet_name = "Results", p_cols = character(0))

    if (!all(is.na(icc_tbl$ICC))) {
      p_icc <- ggplot(
        icc_tbl %>%
          dplyr::filter(!is.na(ICC)) %>%
          dplyr::mutate(Phase = factor(Phase, levels = c("Active", "Inactive")), Sex = factor(Sex)),
        aes(x = interaction(Phase, Sex), y = metric, fill = ICC)
      ) +
        geom_tile(color = "white", linewidth = 0.4) +
        geom_text(aes(label = sprintf("%.2f", ICC)), size = 3) +
        scale_fill_gradient(low = "#f3f3f3", high = "#264653", limits = c(0, 1), name = "ICC") +
        theme_minimal(base_size = 10) +
        theme(
          axis.text.x = element_text(angle = 40, hjust = 1),
          panel.grid  = element_blank(),
          plot.title  = element_text(face = "bold", hjust = 0.5)
        ) +
        labs(
          title = "Within-animal ICC (random intercept / total variance)",
          x = "Phase \u00d7 Sex", y = "Metric"
        )

      ggsave(
        file.path(nat_dirs$plots, paste0("icc_heatmap_", run_scope_current, ".svg")),
        p_icc, width = 6, height = 4
      )
    }
    cat("     \u2713 Saved\n")
  } else {
    cat("     - No ICC data available\n")
  }

  # ================================================================
  # #7: Partial correlations + latent sociability factor (PCA)
  # ================================================================
  cat("  [7/10] Partial correlations + latent factor (PCA)...\n")

  anim_avgs <- data_filtered_agg %>%
    dplyr::filter(!is.na(Movement), !is.na(Proximity), !is.na(ActivityIndex)) %>%
    dplyr::group_by(AnimalNum, Group, Sex, Change, Phase) %>%
    dplyr::summarise(
      Movement      = mean(Movement, na.rm = TRUE),
      Proximity     = mean(Proximity, na.rm = TRUE),
      ActivityIndex = mean(ActivityIndex, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(Movement), is.finite(Proximity), is.finite(ActivityIndex))

  # Partial correlations via residualisation
  partial_corr_rows <- list()
  vars3 <- c("Movement", "Proximity", "ActivityIndex")

  for (g in unique(anim_avgs$Group)) {
    for (ph in unique(anim_avgs$Phase)) {
      dsub_pc <- anim_avgs %>% dplyr::filter(Group == g, Phase == ph)
      if (nrow(dsub_pc) < 5) next

      for (i in seq_along(vars3)) {
        for (j in seq_along(vars3)) {
          if (j <= i) next
          ctrl_v <- vars3[-c(i, j)]
          x <- dsub_pc[[vars3[i]]]; y <- dsub_pc[[vars3[j]]]; z <- dsub_pc[[ctrl_v]]
          if (sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0 || sd(z, na.rm = TRUE) == 0) next

          res_x <- tryCatch(residuals(lm(x ~ z)), error = function(e) NULL)
          res_y <- tryCatch(residuals(lm(y ~ z)), error = function(e) NULL)
          if (is.null(res_x) || is.null(res_y)) next

          ct <- tryCatch(cor.test(res_x, res_y), error = function(e) NULL)
          if (is.null(ct)) next

          partial_corr_rows[[length(partial_corr_rows) + 1]] <- tibble::tibble(
            Group     = g,
            Phase     = as.character(ph),
            var_x     = vars3[i],
            var_y     = vars3[j],
            ctrl      = ctrl_v,
            partial_r = as.numeric(ct$estimate),
            p_value   = ct$p.value,
            n         = nrow(dsub_pc)
          )
        }
      }
    }
  }

  if (length(partial_corr_rows) > 0) {
    partial_tbl <- dplyr::bind_rows(partial_corr_rows) %>%
      dplyr::mutate(p_BH = p.adjust(p_value, method = "BH"))
    write_pub_xlsx(partial_tbl,
      file.path(nat_dirs$tables, paste0("partial_correlations_", run_scope_current, ".xlsx")),
      sheet_name = "Results", p_cols = character(0))
  }

  # PCA — latent sociability factor
  pca_results_nat   <- list()
  pca_scores_all_nat <- list()

  for (ph in unique(anim_avgs$Phase)) {
    dsub_pca <- anim_avgs %>%
      dplyr::filter(Phase == ph) %>%
      dplyr::select(AnimalNum, Group, Sex, Movement, Proximity, ActivityIndex)
    if (nrow(dsub_pca) < 5) next

    mat <- scale(as.matrix(dsub_pca[, c("Movement", "Proximity", "ActivityIndex")]))
    pca <- tryCatch(prcomp(mat, center = FALSE, scale. = FALSE), error = function(e) NULL)
    if (is.null(pca)) next

    pca_results_nat[[as.character(ph)]] <- list(
      phase         = as.character(ph),
      sdev          = pca$sdev,
      rotation      = as.data.frame(pca$rotation),
      var_explained = (pca$sdev^2) / sum(pca$sdev^2)
    )

    scores_df <- as.data.frame(pca$x)
    names(scores_df) <- paste0("PC", seq_len(ncol(scores_df)))
    scores_df$AnimalNum <- as.character(dsub_pca$AnimalNum)
    scores_df$Group     <- dsub_pca$Group
    scores_df$Sex       <- dsub_pca$Sex
    scores_df$Phase     <- as.character(ph)
    pca_scores_all_nat[[as.character(ph)]] <- scores_df
  }

  if (length(pca_results_nat) > 0) {
    loadings_tbl <- dplyr::bind_rows(lapply(names(pca_results_nat), function(ph) {
      r   <- pca_results_nat[[ph]]
      rot <- r$rotation; rot$variable <- rownames(rot)
      rot$Phase             <- r$phase
      rot$var_explained_PC1 <- r$var_explained[1]
      rot$var_explained_PC2 <- r$var_explained[2]
      rot
    }))
    write_pub_xlsx(loadings_tbl,
      file.path(nat_dirs$tables, paste0("pca_loadings_", run_scope_current, ".xlsx")),
      sheet_name = "Results", p_cols = character(0))

    pca_scores_combined <- dplyr::bind_rows(pca_scores_all_nat)
    write_pub_xlsx(pca_scores_combined,
      file.path(nat_dirs$tables, paste0("pca_scores_", run_scope_current, ".xlsx")),
      sheet_name = "Results", p_cols = character(0))

    p_pca <- ggplot(
      pca_scores_combined %>%
        dplyr::filter(!is.na(PC1)) %>%
        dplyr::mutate(Group = factor(Group, levels = c("CON", "RES", "SUS"))),
      aes(x = Group, y = PC1, fill = Group, color = Group)
    ) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
      geom_violin(alpha = 0.45, trim = TRUE) +
      geom_jitter(width = 0.11, size = 1.8, alpha = 0.65) +
      stat_summary(fun = mean, geom = "crossbar", width = 0.3, linewidth = 0.5, color = "black") +
      scale_fill_manual(values  = group_colors) +
      scale_color_manual(values = group_colors) +
      facet_wrap(~ Phase, scales = "free_y") +
      theme_classic(base_size = 10) +
      theme(
        legend.position = "none",
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5)
      ) +
      labs(
        title = "Latent sociability factor (PC1)",
        x = "Group", y = "PC1 score"
      )

    ggsave(
      file.path(nat_dirs$plots, paste0("latent_factor_PC1_", run_scope_current, ".svg")),
      p_pca, width = 7, height = 4.5
    )
  }
  cat("     \u2713 Saved\n")

  # ================================================================
  # #8: TOST equivalence testing for RES-CON
  # ================================================================
  cat("  [8/10] TOST equivalence for RES-CON...\n")

  tost_rows <- list()
  metric_sds_nat <- sapply(metrics_to_analyze, function(m) sd(data_filtered_agg[[m]], na.rm = TRUE))

  rescon_rows <- integrated_full %>%
    dplyr::filter(pair == "RES-CON", model == "lmer", !is.na(estimate), !is.na(SE))

  if (nrow(rescon_rows) > 0) {
    for (i in seq_len(nrow(rescon_rows))) {
      row   <- rescon_rows[i, ]
      m_sd  <- metric_sds_nat[row$metric]
      if (is.na(m_sd) || m_sd == 0) next

      bound   <- 0.5 * m_sd
      est_val <- as.numeric(row$estimate)
      se_val  <- as.numeric(row$SE)

      # Approximate df from SE ratio (conservative floor at 5)
      df_t <- max(5L, as.integer(abs(est_val / se_val) * 20))

      t_lower  <- (est_val - (-bound)) / se_val   # test vs lower bound
      t_upper  <- (est_val -   bound)  / se_val   # test vs upper bound
      p_lower  <- pt(t_lower, df = df_t, lower.tail = FALSE)
      p_upper  <- pt(t_upper, df = df_t, lower.tail = TRUE)
      p_equiv  <- max(p_lower, p_upper)

      tost_rows[[length(tost_rows) + 1]] <- tibble::tibble(
        metric     = row$metric,
        window     = row$window,
        Change     = row$Change,
        Sex        = row$Sex,
        Phase      = row$Phase,
        pair       = row$pair,
        estimate   = est_val,
        SE         = se_val,
        bound      = bound,
        t_lower    = t_lower,
        t_upper    = t_upper,
        p_lower    = p_lower,
        p_upper    = p_upper,
        p_TOST     = p_equiv,
        equivalent = p_equiv < 0.05,
        conclusion = ifelse(
          p_equiv < 0.05,
          paste0("Equivalent (within \u00b1", sprintf("%.2f", bound), ")"),
          paste0("Not equivalent (\u00b1", sprintf("%.2f", bound), ")")
        )
      )
    }
  }

  if (length(tost_rows) > 0) {
    tost_tbl <- dplyr::bind_rows(tost_rows)
    write_pub_xlsx(tost_tbl,
      file.path(nat_dirs$tables, paste0("tost_equivalence_RESCON_", run_scope_current, ".xlsx")),
      sheet_name = "Results", p_cols = character(0))
    n_equiv <- sum(tost_tbl$equivalent, na.rm = TRUE)
    cat(sprintf("     \u2713 TOST: %d / %d RES-CON comparisons equivalent\n", n_equiv, nrow(tost_tbl)))
  } else {
    cat("     - No RES-CON data for TOST\n")
  }

  
  # ================================================================
  # #11: GAMM trajectory analysis (mgcv::bam)
  # ================================================================
  cat("  [11] GAMM trajectory analysis...\n")

  if (!requireNamespace("mgcv", quietly = TRUE)) {
    cat("     - Skipped: mgcv not installed. Install with install.packages('mgcv')\n")
  } else {
    library(mgcv)

    RUN_GAMM <- TRUE
    GAMM_K   <- 6L    # Basis dimension target for shared and group-deviation smooths
    GAMM_SELECT <- TRUE
    GAMM_GAMMA  <- 1.5
    GAMM_N_GRID <- 100L
    GAMM_AUC_NSIM <- 1000L
    GAMM_AR1 <- TRUE
    GAMM_AR1_MIN_RHO <- 0.05
    GAMM_BH_EXCLUDE_RE_SMOOTH <- TRUE
    set.seed(123L)     # Reproducibility for AUC bootstrap (GAMM_AUC_NSIM draws)

    add_gamm_halfhour_within_cc <- function(df, phase_local = FALSE) {
      if (isTRUE(phase_local) && all(c("Change", "Phase") %in% names(df))) {
        time_lookup <- df %>%
          dplyr::distinct(Change, HalfHourElapsed, Phase) %>%
          dplyr::arrange(Change, HalfHourElapsed) %>%
          dplyr::group_by(Change) %>%
          dplyr::mutate(
            .phase_block = cumsum(
              dplyr::row_number() == 1L |
                Phase != dplyr::lag(Phase) |
                (HalfHourElapsed - dplyr::lag(HalfHourElapsed)) > 1
            )
          ) %>%
          dplyr::group_by(Change, .phase_block) %>%
          dplyr::mutate(
            HalfHourWithinCC0 = as.numeric(HalfHourElapsed - min(HalfHourElapsed, na.rm = TRUE))
          ) %>%
          dplyr::ungroup() %>%
          dplyr::transmute(
            Change, HalfHourElapsed, Phase,
            GammTimeBlock = .phase_block,
            HalfHourWithinCC0
          )

        return(
          df %>%
            dplyr::select(-dplyr::any_of(c("HalfHourWithinCC0", "GammTimeBlock"))) %>%
            dplyr::left_join(time_lookup, by = c("Change", "HalfHourElapsed", "Phase"))
        )
      }

      if ("Change" %in% names(df)) {
        df %>%
          dplyr::group_by(Change) %>%
          dplyr::mutate(
            GammTimeBlock = 1L,
            HalfHourWithinCC0 = as.numeric(HalfHourElapsed - min(HalfHourElapsed, na.rm = TRUE))
          ) %>%
          dplyr::ungroup()
      } else {
        df %>%
          dplyr::mutate(
            GammTimeBlock = 1L,
            HalfHourWithinCC0 = as.numeric(HalfHourElapsed - min(HalfHourElapsed, na.rm = TRUE))
          )
      }
    }

    compute_gamm_resid_acf <- function(model_obj, data_obj, metric, Change, Sex, Phase,
                                       window = "all", max_lag = 6L) {
      resid_vec <- tryCatch(as.numeric(stats::residuals(model_obj, type = "pearson")), error = function(e) NULL)
      if (is.null(resid_vec) || length(resid_vec) != nrow(data_obj)) {
        resid_vec <- tryCatch(as.numeric(stats::residuals(model_obj, type = "response")), error = function(e) NULL)
      }
      if (is.null(resid_vec) || length(resid_vec) != nrow(data_obj)) {
        return(list(by_lag = tibble::tibble(), summary = tibble::tibble()))
      }

      d_acf <- data_obj %>%
        dplyr::mutate(.resid = resid_vec) %>%
        dplyr::filter(is.finite(.resid), !is.na(AnimalNum), !is.na(Change), !is.na(GammTimeBlock), !is.na(HalfHourWithinCC0)) %>%
        dplyr::arrange(AnimalNum, Change, GammTimeBlock, HalfHourWithinCC0) %>%
        dplyr::mutate(.acf_unit = interaction(AnimalNum, Change, GammTimeBlock, drop = TRUE))
      if (nrow(d_acf) < 2) {
        return(list(by_lag = tibble::tibble(), summary = tibble::tibble()))
      }

      acf_rows <- lapply(split(d_acf, d_acf$.acf_unit), function(ds) {
        if (nrow(ds) < 2 || stats::sd(ds$.resid, na.rm = TRUE) == 0) return(NULL)
        lag_max_i <- min(as.integer(max_lag), nrow(ds) - 1L)
        if (lag_max_i < 1L) return(NULL)
        acf_obj <- tryCatch(
          stats::acf(ds$.resid, lag.max = lag_max_i, plot = FALSE, na.action = stats::na.pass),
          error = function(e) NULL
        )
        if (is.null(acf_obj)) return(NULL)
        acf_vals <- as.numeric(acf_obj$acf)[-1]
        if (length(acf_vals) == 0) return(NULL)
        tibble::tibble(
          AnimalNum = as.character(ds$AnimalNum[1]),
          ChangeUnit = as.character(ds$Change[1]),
          GammTimeBlock = as.integer(ds$GammTimeBlock[1]),
          lag = seq_along(acf_vals),
          acf = acf_vals,
          n_obs = nrow(ds)
        )
      })
      acf_rows <- dplyr::bind_rows(acf_rows)
      if (nrow(acf_rows) == 0) {
        return(list(by_lag = tibble::tibble(), summary = tibble::tibble()))
      }

      by_lag <- acf_rows %>%
        dplyr::group_by(lag) %>%
        dplyr::summarise(
          mean_acf = mean(acf, na.rm = TRUE),
          median_acf = stats::median(acf, na.rm = TRUE),
          sd_acf = stats::sd(acf, na.rm = TRUE),
          mean_abs_acf = mean(abs(acf), na.rm = TRUE),
          n_acf_units = dplyr::n(),
          n_animals = dplyr::n_distinct(AnimalNum),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          metric = metric,
          Change = as.character(Change),
          Sex = as.character(Sex),
          Phase = as.character(Phase),
          window = window
        ) %>%
        dplyr::select(metric, Change, Sex, Phase, window, lag,
                      mean_acf, median_acf, sd_acf, mean_abs_acf, n_acf_units, n_animals)

      lag1_summary <- by_lag %>%
        dplyr::filter(lag == 1L) %>%
        dplyr::mutate(
          lag1_flag = dplyr::case_when(
            mean_abs_acf < 0.20 ~ "low",
            mean_abs_acf < 0.30 ~ "moderate",
            TRUE ~ "high"
          ),
          ar1_recommendation = dplyr::case_when(
            mean_abs_acf < 0.20 ~ "Residual lag-1 autocorrelation is low (<0.2): current GAMM is likely adequate.",
            mean_abs_acf < 0.30 ~ "Residual lag-1 autocorrelation is modest: AR1 optional as sensitivity check.",
            TRUE ~ "Residual lag-1 autocorrelation is elevated: consider AR1 sensitivity analysis."
          )
        ) %>%
        dplyr::rename(lag1_mean_acf = mean_acf, lag1_abs_acf = mean_abs_acf) %>%
        dplyr::select(metric, Change, Sex, Phase, window,
                      lag1_mean_acf, lag1_abs_acf, median_acf, sd_acf, n_animals,
                      lag1_flag, ar1_recommendation)

      list(by_lag = by_lag, summary = lag1_summary)
    }

    decode_group_term_label <- function(term, grp_factor) {
      if (!is.character(term) || length(term) != 1L || is.na(term)) return(NA_character_)
      if (term == "(Intercept)") return("(Intercept) [model intercept]")
      if (!grepl("^Group", term)) return(term)

      cmat <- tryCatch(stats::contrasts(grp_factor), error = function(e) NULL)
      lvls <- levels(grp_factor)
      if (is.null(cmat) || ncol(cmat) == 0 || length(lvls) == 0) {
        return(paste0(term, " [Group effect; contrast matrix unavailable]"))
      }

      key <- sub("^Group", "", term)
      c_names <- colnames(cmat)
      j <- match(key, c_names)
      if (is.na(j)) {
        key_num <- suppressWarnings(as.integer(key))
        if (is.finite(key_num) && key_num >= 1L && key_num <= ncol(cmat)) j <- key_num
      }
      if (is.na(j)) return(paste0(term, " [Group effect; contrast column not resolved]"))

      coefs <- as.numeric(cmat[, j])
      parts <- paste0(lvls, "=", format(round(coefs, 3), trim = TRUE, nsmall = 0))
      paste0(term, " [", paste(parts, collapse = ", "), "]")
    }

    extract_gamm_individual_auc <- function(model_obj, data_obj, th_seq,
                                            metric, Change, Sex, Phase,
                                            window = "all") {
      animal_cols <- intersect(c("AnimalNum", "Group", "Batch"), names(data_obj))
      animal_info <- data_obj %>%
        dplyr::distinct(dplyr::across(dplyr::all_of(animal_cols))) %>%
        dplyr::arrange(AnimalNum)

      if (!"Batch" %in% names(animal_info)) {
        animal_info <- animal_info %>% dplyr::mutate(Batch = NA_character_)
      }
      if (nrow(animal_info) == 0 || length(th_seq) < 2) {
        return(tibble::tibble())
      }

      pred_rows <- lapply(seq_len(nrow(animal_info)), function(i) {
        animal_id <- animal_info$AnimalNum[i]
        group_id <- as.character(animal_info$Group[i])

        # Proteomics-facing AUC: use each animal's own observed time support
        # (not the shared global grid), to preserve subject-level trajectory differences.
        x_i <- data_obj %>%
          dplyr::filter(AnimalNum == animal_id, Group == group_id) %>%
          dplyr::pull(HalfHourWithinCC0)
        x_i <- sort(unique(as.numeric(x_i[is.finite(x_i)])))
        if (length(x_i) < 2) {
          # Fallback only when an animal has <2 valid observed time points
          x_i <- sort(unique(as.numeric(th_seq[is.finite(th_seq)])))
        }
        if (length(x_i) < 2) return(NULL)

        nd <- data.frame(
          HalfHourWithinCC0 = x_i,
          Group = factor(group_id, levels = levels(data_obj$Group)),
          AnimalNum = factor(animal_id, levels = levels(data_obj$AnimalNum))
        )
        pred_fit <- tryCatch(
          as.numeric(mgcv::predict.bam(model_obj, newdata = nd, type = "response")),
          error = function(e) NULL
        )
        if (is.null(pred_fit) || length(pred_fit) != nrow(nd)) return(NULL)

        th_range_i <- max(x_i, na.rm = TRUE) - min(x_i, na.rm = TRUE)
        auc_i <- trapz_auc(x_i, pred_fit)

        tibble::tibble(
          AnimalNum = as.character(animal_id),
          Batch = as.character(animal_info$Batch[i]),
          Group = group_id,
          metric = metric,
          Change = as.character(Change),
          Sex = as.character(Sex),
          Phase = as.character(Phase),
          window = window,
          AUC = auc_i,
          AUC_norm = ifelse(is.finite(th_range_i) && th_range_i > 0, auc_i / th_range_i, NA_real_),
          prediction_type = "subject_specific_gamm_observed_grid",
          n_grid = length(x_i),
          halfhour_min = min(x_i, na.rm = TRUE),
          halfhour_max = max(x_i, na.rm = TRUE)
        )
      })

      dplyr::bind_rows(pred_rows)
    }

    gamm_results     <- list()
    gamm_smooth_rows <- list()
    gamm_pairwise_rows <- list()
    data_gamm_base <- add_gamm_halfhour_within_cc(data_filtered_agg, phase_local = isTRUE(includePhase))
    gamm_time_axis_label <- if (isTRUE(includePhase)) {
      "Half-hour within phase block (0-based)"
    } else {
      "Half-hour within cage change (0-based)"
    }

    changes_gam <- if (includeChange) unique(data_filtered_agg$Change) else "allChanges"
    sexes_gam   <- if (includeSex)    unique(data_filtered_agg$Sex)    else "allSexes"
    phases_gam  <- if (includePhase)  unique(data_filtered_agg$Phase)  else "allPhases"

    for (metric in metrics_to_analyze) {
      for (ch in changes_gam) {
        for (sx in sexes_gam) {
          for (ph in phases_gam) {
            dsub_gam <- data_gamm_base %>%
              dplyr::filter(
                (!includeChange | Change == ch),
                (!includeSex    | Sex    == sx),
                (!includePhase  | Phase  == ph),
                !is.na(.data[[metric]]), !is.na(AnimalNum), !is.na(Group)
              )
            if (nrow(dsub_gam) < max(min_obs, 20) ||
                dplyr::n_distinct(dsub_gam$AnimalNum) < min_animals ||
                dplyr::n_distinct(dsub_gam$Group) < min_groups) next

            dsub_gam <- dsub_gam %>%
              dplyr::mutate(
                AnimalNum = factor(AnimalNum),
                Group     = factor(Group, levels = c("CON", "RES", "SUS"))
              ) %>%
              dplyr::filter(!is.na(HalfHourWithinCC0))

            # Required ordering for AR.start boundaries (independent AR1 process per animal x cage change)
            dsub_gam <- dplyr::arrange(dsub_gam, AnimalNum, Change, GammTimeBlock, HalfHourWithinCC0)
            ar_start_vec <- seq_len(nrow(dsub_gam)) == 1L |
              dsub_gam$AnimalNum != dplyr::lag(dsub_gam$AnimalNum) |
              dsub_gam$Change != dplyr::lag(dsub_gam$Change) |
              dsub_gam$GammTimeBlock != dplyr::lag(dsub_gam$GammTimeBlock)

            n_uniq_t_gam <- dplyr::n_distinct(dsub_gam$HalfHourWithinCC0)
            if (n_uniq_t_gam < 8) next
            k_use <- min(GAMM_K, floor(n_uniq_t_gam / 2L))

            gam_formula <- as.formula(paste0(
              metric, " ~ Group + s(HalfHourWithinCC0, k = ", k_use, ", bs = 'tp') + ",
              "s(HalfHourWithinCC0, by = Group, k = ", k_use, ", bs = 'tp') + s(AnimalNum, bs = 're')"
            ))

            # Fit GAMM: shared time smooth + group-specific deviations + random intercept per animal
            # bs="tp" = thin-plate regression spline; re = random effect
            m_gam0 <- tryCatch(
              suppressMessages(suppressWarnings(
                mgcv::bam(
                  formula = gam_formula,
                  data    = dsub_gam,
                  method  = "fREML",
                  family  = gaussian(),
                  discrete = TRUE,
                  select   = GAMM_SELECT,
                  gamma    = GAMM_GAMMA
                )
              )),
              error = function(e) NULL
            )
            if (is.null(m_gam0)) next

            rho_est <- tryCatch({
              resid0 <- as.numeric(stats::residuals(m_gam0, type = "pearson"))
              if (length(resid0) == nrow(dsub_gam)) {
                animal_acfs <- tapply(seq_along(resid0), interaction(dsub_gam$AnimalNum, dsub_gam$Change, dsub_gam$GammTimeBlock, drop = TRUE), function(idx) {
                  r <- resid0[idx]
                  r <- r[is.finite(r)]
                  if (length(r) < 3L) return(NA_real_)
                  as.numeric(stats::acf(r, lag.max = 1L, plot = FALSE)$acf[2L])
                })
                mean(unlist(animal_acfs), na.rm = TRUE)
              } else 0
            }, error = function(e) 0)
            if (!is.finite(rho_est)) rho_est <- 0
            rho_est <- max(-0.99, min(0.99, rho_est))

            ar1_applied <- FALSE
            rho_used <- 0
            m_gam <- m_gam0
            if (GAMM_AR1 && abs(rho_est) > GAMM_AR1_MIN_RHO) {
              m_gam1 <- tryCatch(
                suppressMessages(suppressWarnings(
                  mgcv::bam(
                    formula = gam_formula,
                    data    = dsub_gam,
                    method  = "fREML",
                    family  = gaussian(),
                    rho = rho_est,
                    AR.start = ar_start_vec,
                    discrete = FALSE,
                    select   = GAMM_SELECT,
                    gamma    = GAMM_GAMMA
                  )
                )),
                error = function(e) NULL
              )
              if (!is.null(m_gam1)) {
                m_gam <- m_gam1
                ar1_applied <- TRUE
                rho_used <- rho_est
              }
            }

            # Parametric summary (Group intercept differences)
            sm <- summary(m_gam)
            ptbl <- as.data.frame(sm$p.table) %>%
              tibble::rownames_to_column("Term") %>%
              tibble::as_tibble()
            if (ncol(ptbl) >= 5) {
              names(ptbl)[2:5] <- c("Estimate", "SE", "t_value", "p_value")
            }
            ptbl <- ptbl %>%
              dplyr::mutate(
                Contrast = vapply(Term, decode_group_term_label, character(1), grp_factor = dsub_gam$Group),
                metric = metric, Change = as.character(ch),
                Sex = as.character(sx), Phase = as.character(ph),
                AIC = AIC(m_gam), BIC = BIC(m_gam),
                dev_expl = sm$dev.expl
              ) %>%
              dplyr::relocate(Contrast, .after = Term)

            # Smooth term significance (edf + F/Chi-sq)
            stbl <- as.data.frame(sm$s.table) %>%
              tibble::rownames_to_column("Smooth") %>%
              tibble::as_tibble()
            if ("Ref.df" %in% names(stbl)) names(stbl)[names(stbl) == "Ref.df"] <- "Ref_df"
            if ("F" %in% names(stbl)) names(stbl)[names(stbl) == "F"] <- "F_stat"
            if ("p-value" %in% names(stbl)) names(stbl)[names(stbl) == "p-value"] <- "p_value"
            stbl <- stbl %>%
              dplyr::mutate(
                metric = metric, Change = as.character(ch),
                Sex = as.character(sx), Phase = as.character(ph)
              )

            # Basis-dimension diagnostics (k-check)
            kchk <- tryCatch(mgcv::k.check(m_gam), error = function(e) NULL)
            if (!is.null(kchk)) {
              kdiag_tbl <- as.data.frame(kchk) %>%
                tibble::rownames_to_column("Smooth") %>%
                tibble::as_tibble()
              if ("k-index" %in% names(kdiag_tbl)) names(kdiag_tbl)[names(kdiag_tbl) == "k-index"] <- "k_index"
              if ("p-value" %in% names(kdiag_tbl)) names(kdiag_tbl)[names(kdiag_tbl) == "p-value"] <- "kcheck_p_value"
            } else {
              kdiag_tbl <- tibble::tibble(Smooth = NA_character_, k_index = NA_real_, kcheck_p_value = NA_real_)
            }
            kdiag_tbl <- kdiag_tbl %>%
              dplyr::mutate(
                metric = metric,
                Change = as.character(ch),
                Sex = as.character(sx),
                Phase = as.character(ph),
                window = "all",
                k_target = GAMM_K,
                k_used = k_use,
                select_enabled = GAMM_SELECT,
                gamma_value = GAMM_GAMMA,
                dev_expl = sm$dev.expl,
                AIC = AIC(m_gam),
                BIC = BIC(m_gam),
                rho_estimated = rho_est,
                rho_used = rho_used,
                ar1_applied = ar1_applied
              )

            acf_diag_tbl <- compute_gamm_resid_acf(
              model_obj = m_gam,
              data_obj = dsub_gam,
              metric = metric,
              Change = ch,
              Sex = sx,
              Phase = ph,
              window = "all"
            )

            # Prediction grid for plotting and AUC
            th_seq <- seq(min(dsub_gam$HalfHourWithinCC0), max(dsub_gam$HalfHourWithinCC0), length.out = GAMM_N_GRID)
            grp_levels <- levels(dsub_gam$Group)

            pred_list_gam <- lapply(grp_levels, function(g) {
              nd <- data.frame(
                HalfHourWithinCC0 = th_seq,
                Group     = factor(g, levels = grp_levels),
                AnimalNum = factor(dsub_gam$AnimalNum[1], levels = levels(dsub_gam$AnimalNum))
              )
              pred_obj <- tryCatch(
                mgcv::predict.bam(m_gam, newdata = nd, type = "response",
                                  exclude = "s(AnimalNum)", se.fit = TRUE),
                error = function(e) NULL
              )
              if (is.null(pred_obj)) return(NULL)
              data.frame(
                HalfHourWithinCC0 = th_seq, Group = g,
                predicted = as.numeric(pred_obj$fit),
                se        = as.numeric(pred_obj$se.fit),
                lwr       = as.numeric(pred_obj$fit) - 1.96 * as.numeric(pred_obj$se.fit),
                upr       = as.numeric(pred_obj$fit) + 1.96 * as.numeric(pred_obj$se.fit),
                metric    = metric,
                Change    = as.character(ch), Sex = as.character(sx), Phase = as.character(ph)
              )
            })
            pred_df_gam <- dplyr::bind_rows(pred_list_gam)

            # AUC per group (trapezoidal)
            auc_gam <- lapply(grp_levels, function(g) {
              sub <- pred_df_gam %>% dplyr::filter(Group == g)
              tibble::tibble(
                Group  = g,
                AUC    = trapz_auc(sub$HalfHourWithinCC0, sub$predicted),
                metric = metric,
                Change = as.character(ch), Sex = as.character(sx), Phase = as.character(ph)
              )
            })
            auc_gam_tbl <- dplyr::bind_rows(auc_gam)
            auc_gam_individual_tbl <- extract_gamm_individual_auc(
              model_obj = m_gam,
              data_obj = dsub_gam,
              th_seq = th_seq,
              metric = metric,
              Change = ch,
              Sex = sx,
              Phase = ph,
              window = "all"
            )

            # Pairwise differences at each grid point + integrated AUC statistics
            pw_pairs <- list(c("RES", "CON"), c("SUS", "CON"), c("SUS", "RES"))
            pw_rows <- list()
            pw_auc_rows <- list()

            beta_hat <- tryCatch(stats::coef(m_gam), error = function(e) NULL)
            Vp_hat <- tryCatch(mgcv::vcov.gam(m_gam, unconditional = TRUE), error = function(e) {
              tryCatch(vcov(m_gam), error = function(e2) NULL)
            })
            beta_draws <- NULL
            if (!is.null(beta_hat) && !is.null(Vp_hat) &&
                is.matrix(Vp_hat) && nrow(Vp_hat) == length(beta_hat) &&
                requireNamespace("MASS", quietly = TRUE)) {
              beta_draws <- tryCatch(
                MASS::mvrnorm(n = GAMM_AUC_NSIM, mu = beta_hat, Sigma = Vp_hat),
                error = function(e) NULL
              )
              if (!is.null(beta_draws) && !is.matrix(beta_draws)) {
                beta_draws <- matrix(beta_draws, nrow = 1)
              }
            }

            for (pp in pw_pairs) {
              g1 <- pp[1]; g0 <- pp[2]
              if (!(g1 %in% grp_levels && g0 %in% grp_levels)) next

              sub1 <- pred_df_gam %>% dplyr::filter(Group == g1)
              sub0 <- pred_df_gam %>% dplyr::filter(Group == g0)
              diff_vec <- sub1$predicted - sub0$predicted
              # NOTE: this initial diff_se assumes independence of group predictions (ignores
              # covariance). It is overwritten by the lpmatrix delta-method below when
              # predict.bam succeeds, which is the statistically correct approach.
              diff_se <- sqrt(sub1$se^2 + sub0$se^2)

              nd1 <- data.frame(
                HalfHourWithinCC0 = th_seq,
                Group = factor(g1, levels = grp_levels),
                AnimalNum = factor(dsub_gam$AnimalNum[1], levels = levels(dsub_gam$AnimalNum))
              )
              nd0 <- data.frame(
                HalfHourWithinCC0 = th_seq,
                Group = factor(g0, levels = grp_levels),
                AnimalNum = factor(dsub_gam$AnimalNum[1], levels = levels(dsub_gam$AnimalNum))
              )
              X1 <- tryCatch(mgcv::predict.bam(m_gam, newdata = nd1, type = "lpmatrix"), error = function(e) NULL)
              X0 <- tryCatch(mgcv::predict.bam(m_gam, newdata = nd0, type = "lpmatrix"), error = function(e) NULL)
              Xdiff <- NULL
              if (!is.null(X1) && !is.null(X0) && ncol(X1) == ncol(X0)) {
                Xdiff <- X1 - X0
                re_cols <- grepl("^s\\(AnimalNum\\)", colnames(Xdiff))
                if (any(re_cols)) Xdiff[, re_cols] <- 0
                if (!is.null(beta_hat) && ncol(Xdiff) == length(beta_hat)) {
                  diff_vec <- as.numeric(Xdiff %*% beta_hat)
                }
                if (!is.null(Vp_hat) && ncol(Xdiff) == nrow(Vp_hat)) {
                  diff_var <- rowSums((Xdiff %*% Vp_hat) * Xdiff)
                  diff_se <- sqrt(pmax(diff_var, 0))
                }
              }

              diff_lwr <- diff_vec - 1.96 * diff_se
              diff_upr <- diff_vec + 1.96 * diff_se

              pw_rows[[length(pw_rows) + 1]] <- tibble::tibble(
                HalfHourWithinCC0 = th_seq,
                diff = diff_vec,
                diff_lwr = diff_lwr,
                diff_upr = diff_upr,
                pair = paste0(g1, "-", g0),
                metric = metric,
                Change = as.character(ch),
                Sex = as.character(sx),
                Phase = as.character(ph),
                window = "all"
              )

              auc_diff <- trapz_auc(th_seq, diff_vec)
              th_range <- max(th_seq, na.rm = TRUE) - min(th_seq, na.rm = TRUE)
              auc_diff_norm <- ifelse(is.finite(th_range) && th_range > 0, auc_diff / th_range, NA_real_)
              frac_sig <- mean(diff_lwr > 0 | diff_upr < 0, na.rm = TRUE)

              auc_ci_low <- NA_real_
              auc_ci_high <- NA_real_
              p_AUC_raw <- NA_real_
              if (!is.null(Xdiff) && !is.null(beta_draws) && ncol(Xdiff) == ncol(beta_draws)) {
                diff_draw_mat <- Xdiff %*% t(beta_draws)
                auc_draws <- apply(diff_draw_mat, 2, function(y) trapz_auc(th_seq, y))
                auc_draws <- auc_draws[is.finite(auc_draws)]
                if (length(auc_draws) >= 50) {
                  auc_ci_low <- as.numeric(stats::quantile(auc_draws, 0.025, na.rm = TRUE))
                  auc_ci_high <- as.numeric(stats::quantile(auc_draws, 0.975, na.rm = TRUE))
                  p_low <- mean(auc_draws <= 0)
                  p_high <- mean(auc_draws >= 0)
                  p_AUC_raw <- min(1, 2 * min(p_low, p_high))
                }
              }

              pw_auc_rows[[length(pw_auc_rows) + 1]] <- tibble::tibble(
                metric = metric,
                Change = as.character(ch),
                Sex = as.character(sx),
                Phase = as.character(ph),
                window = "all",
                pair = paste0(g1, "-", g0),
                AUC_diff = auc_diff,
                AUC_diff_norm = auc_diff_norm,
                AUC_ci_low = auc_ci_low,
                AUC_ci_high = auc_ci_high,
                p_AUC_raw = p_AUC_raw,
                frac_sig = frac_sig
              )
            }

            gamm_results[[length(gamm_results) + 1]] <- list(
              metric      = metric, Change = as.character(ch),
              Sex         = as.character(sx), Phase = as.character(ph),
              model       = m_gam,
              predictions = pred_df_gam,
              auc         = auc_gam_tbl,
              auc_individual = auc_gam_individual_tbl,
              parametric  = ptbl,
              smooths     = stbl,
              pairwise    = dplyr::bind_rows(pw_rows),
              pairwise_auc = dplyr::bind_rows(pw_auc_rows)
              ,diag       = kdiag_tbl,
              acf         = acf_diag_tbl$by_lag,
              acf_summary = acf_diag_tbl$summary
            )
          }
        }
      }
    }

    # BH correction for GAMM p-values by clearly defined families
    # Family A: parametric Group terms, per metric x window
    # Family B: smooth terms, per metric x window
    apply_gamm_bh_families <- function(param_tbl, smooth_tbl) {
      # BH families are stratified by metric x Change x Sex x Phase x window so that
      # corrections are applied within each experimental stratum, not across all strata
      # pooled together (which would be overly conservative).
      if (nrow(param_tbl) > 0) {
        param_tbl <- param_tbl %>%
          dplyr::mutate(
            metric    = if ("metric" %in% names(.)) as.character(metric) else "unknown_metric",
            window    = if ("window" %in% names(.)) as.character(window) else "all",
            Change    = if ("Change" %in% names(.)) as.character(Change) else "all",
            Sex       = if ("Sex"    %in% names(.)) as.character(Sex)    else "all",
            Phase     = if ("Phase"  %in% names(.)) as.character(Phase)  else "all",
            Term      = if ("Term" %in% names(.)) as.character(Term) else NA_character_,
            p_value   = if ("p_value" %in% names(.)) p_value else NA_real_,
            p_value_raw = suppressWarnings(as.numeric(p_value)),
            p_BH = NA_real_
          ) %>%
          dplyr::group_by(metric, Change, Sex, Phase, window) %>%
          dplyr::group_modify(~ {
            d <- .x
            idx <- which(grepl("^Group", d$Term) & is.finite(d$p_value_raw))
            if (length(idx) > 0) d$p_BH[idx] <- p.adjust(d$p_value_raw[idx], method = "BH")
            d
          }) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(
            p_for_table = dplyr::coalesce(p_BH, p_value_raw)
          )
      }

      if (nrow(smooth_tbl) > 0) {
        smooth_tbl <- smooth_tbl %>%
          dplyr::mutate(
            metric    = if ("metric" %in% names(.)) as.character(metric) else "unknown_metric",
            window    = if ("window" %in% names(.)) as.character(window) else "all",
            Change    = if ("Change" %in% names(.)) as.character(Change) else "all",
            Sex       = if ("Sex"    %in% names(.)) as.character(Sex)    else "all",
            Phase     = if ("Phase"  %in% names(.)) as.character(Phase)  else "all",
            Smooth    = if ("Smooth" %in% names(.)) as.character(Smooth) else NA_character_,
            p_value   = if ("p_value" %in% names(.)) p_value else NA_real_,
            p_value_raw = suppressWarnings(as.numeric(p_value)),
            p_BH = NA_real_
          ) %>%
          dplyr::group_by(metric, Change, Sex, Phase, window) %>%
          dplyr::group_modify(~ {
            d <- .x
            idx <- which(
              is.finite(d$p_value_raw) &
                (!GAMM_BH_EXCLUDE_RE_SMOOTH | !grepl("^s\\(AnimalNum\\)", d$Smooth))
            )
            if (length(idx) > 0) d$p_BH[idx] <- p.adjust(d$p_value_raw[idx], method = "BH")
            d
          }) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(
            p_for_table = dplyr::coalesce(p_BH, p_value_raw)
          )
      }

      list(param = param_tbl, smooth = smooth_tbl)
    }

    plot_gamm_auc_contrasts_nature <- function(df, title, subtitle = NULL) {
      if (is.null(df) || nrow(df) == 0 ||
          !all(c("pair", "AUC_diff", "AUC_ci_low", "AUC_ci_high") %in% names(df))) {
        return(NULL)
      }

      plot_df <- df %>%
        dplyr::filter(
          pair %in% c("RES-CON", "SUS-CON", "SUS-RES"),
          is.finite(AUC_diff)
        ) %>%
        dplyr::mutate(
          AUC_ci_low = dplyr::if_else(is.finite(AUC_ci_low), AUC_ci_low, AUC_diff),
          AUC_ci_high = dplyr::if_else(is.finite(AUC_ci_high), AUC_ci_high, AUC_diff),
          pair = factor(pair, levels = c("RES-CON", "SUS-CON", "SUS-RES")),
          metric = if ("metric" %in% names(.)) as.character(metric) else "metric",
          Change = if ("Change" %in% names(.)) as.character(Change) else "all",
          Sex = if ("Sex" %in% names(.)) as.character(Sex) else "all",
          Phase = if ("Phase" %in% names(.)) as.character(Phase) else "all",
          p_plot = dplyr::coalesce(
            if ("p_AUC_BH" %in% names(.)) suppressWarnings(as.numeric(p_AUC_BH)) else NA_real_,
            if ("p_AUC_raw" %in% names(.)) suppressWarnings(as.numeric(p_AUC_raw)) else NA_real_
          ),
          sig_label = dplyr::case_when(
            is.na(p_plot) ~ "",
            p_plot < 0.001 ~ "***",
            p_plot < 0.01 ~ "**",
            p_plot < 0.05 ~ "*",
            p_plot < 0.10 ~ "\u2020",
            TRUE ~ ""
          ),
          facet_label = paste(metric, Change, Sex, Phase, sep = " | ")
        )

      if (nrow(plot_df) == 0) return(NULL)

      plot_df <- plot_df %>%
        dplyr::group_by(facet_label) %>%
        dplyr::mutate(
          x_abs_max = max(abs(c(AUC_ci_low, AUC_ci_high)), na.rm = TRUE),
          x_abs_max = dplyr::if_else(is.finite(x_abs_max) & x_abs_max > 0, x_abs_max, 1),
          x_pad = x_abs_max * 0.08
        ) %>%
        dplyr::ungroup()

      ggplot(plot_df, aes(x = AUC_diff, y = pair, color = pair)) +
        geom_vline(xintercept = 0, linewidth = 0.28, color = "grey55") +
        geom_segment(
          aes(x = AUC_ci_low, xend = AUC_ci_high, yend = pair),
          linewidth = 0.65,
          lineend = "round"
        ) +
        geom_point(size = 1.9, stroke = 0.25) +
        geom_text(
          aes(x = AUC_ci_high + x_pad * 0.25, label = sig_label),
          color = "black",
          size = 2.4,
          hjust = 0,
          vjust = 0.45,
          show.legend = FALSE
        ) +
        scale_color_manual(values = pair_colors, guide = "none") +
        scale_x_continuous(expand = expansion(mult = c(0.08, 0.22))) +
        facet_wrap(~ facet_label, scales = "free_x", ncol = 1) +
        labs(
          title = title,
          subtitle = subtitle,
          x = "\u0394AUC (Group A - Group B)",
          y = NULL
        ) +
        theme_classic(base_size = 7, base_family = "Arial") +
        theme(
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(color = "black", size = 6.5),
          axis.text.x = element_text(color = "black", size = 6),
          axis.title.x = element_text(size = 7),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 6.5, hjust = 0),
          plot.title = element_text(face = "bold", hjust = 0, size = 8),
          plot.subtitle = element_text(hjust = 0, size = 6.2),
          plot.margin = margin(3, 12, 3, 3)
        )
    }

    if (length(gamm_results) > 0) {
      cat(sprintf("     Fitted %d GAMM model(s)\n", length(gamm_results)))

      # ---- Combine tables ----
      gamm_param_tbl  <- dplyr::bind_rows(lapply(gamm_results, `[[`, "parametric"))
      gamm_smooth_tbl <- dplyr::bind_rows(lapply(gamm_results, `[[`, "smooths"))
      bh_main <- apply_gamm_bh_families(gamm_param_tbl, gamm_smooth_tbl)
      gamm_param_tbl <- bh_main$param
      gamm_smooth_tbl <- bh_main$smooth

      # ---- EDF console readout ----
      edf_console <- gamm_smooth_tbl %>%
        dplyr::filter(!grepl("AnimalNum", Smooth)) %>%
        dplyr::select(metric, Change, Sex, Phase, Smooth, edf, Ref_df) %>%
        dplyr::arrange(metric, Change, Sex, Phase, Smooth)
      if (nrow(edf_console) > 0) {
        cat("     --- EDF Summary (shared + group-deviation smooths; k=", GAMM_K, ", select=", GAMM_SELECT, ") ---\n", sep = "")
        for (i in seq_len(nrow(edf_console))) {
          cat(sprintf("       [%s | %s | %s | %s] %s: edf=%.3f (Ref_df=%.3f)\n",
            edf_console$metric[i], edf_console$Change[i],
            edf_console$Sex[i],    edf_console$Phase[i],
            edf_console$Smooth[i], edf_console$edf[i], edf_console$Ref_df[i]))
        }
        cat("     (edf=1 → linear; edf>1 → non-linear; edf≈0 → shrunk to zero by select)\n")
        cat("     ---\n")
      }

      gamm_auc_tbl    <- dplyr::bind_rows(lapply(gamm_results, `[[`, "auc"))
      gamm_auc_individual_tbl <- dplyr::bind_rows(lapply(gamm_results, `[[`, "auc_individual"))
      gamm_pw_tbl     <- dplyr::bind_rows(lapply(gamm_results, `[[`, "pairwise"))
      gamm_pw_auc_tbl <- dplyr::bind_rows(lapply(gamm_results, `[[`, "pairwise_auc"))
      gamm_diag_tbl   <- dplyr::bind_rows(lapply(gamm_results, `[[`, "diag"))
      gamm_acf_tbl    <- dplyr::bind_rows(lapply(gamm_results, `[[`, "acf"))
      gamm_acf_summary_tbl <- dplyr::bind_rows(lapply(gamm_results, `[[`, "acf_summary"))
      gamm_pred_tbl   <- dplyr::bind_rows(lapply(gamm_results, `[[`, "predictions")) %>%
        dplyr::mutate(Group = factor(Group, levels = c("CON", "RES", "SUS")))

      if (nrow(gamm_acf_summary_tbl) > 0) {
        cat("     --- Residual lag-1 ACF summary ---\n")
        for (i in seq_len(nrow(gamm_acf_summary_tbl))) {
          cat(sprintf(
            "       [%s | %s | %s | %s] lag1 mean=%.3f | abs=%.3f | %s\n",
            gamm_acf_summary_tbl$metric[i], gamm_acf_summary_tbl$Change[i],
            gamm_acf_summary_tbl$Sex[i], gamm_acf_summary_tbl$Phase[i],
            gamm_acf_summary_tbl$lag1_mean_acf[i], gamm_acf_summary_tbl$lag1_abs_acf[i],
            gamm_acf_summary_tbl$lag1_flag[i]
          ))
        }
        cat("     (rule of thumb: |lag-1 ACF| < 0.2 = likely safe; higher values justify AR1 sensitivity)\n")
        cat("     ---\n")
      }

      if (nrow(gamm_pw_auc_tbl) > 0 && "p_AUC_raw" %in% names(gamm_pw_auc_tbl)) {
        gamm_pw_auc_tbl <- gamm_pw_auc_tbl %>%
          dplyr::mutate(
            metric = if ("metric" %in% names(.)) as.character(metric) else "unknown_metric",
            Change = if ("Change" %in% names(.)) as.character(Change) else "all",
            Sex = if ("Sex" %in% names(.)) as.character(Sex) else "all",
            Phase = if ("Phase" %in% names(.)) as.character(Phase) else "all",
            window = if ("window" %in% names(.)) as.character(window) else "all",
            p_AUC_BH = NA_real_,
            bh_family_n_AUC = NA_integer_,
            bh_family_scope_AUC = "metric|Change|Sex|Phase|window"
          ) %>%
          dplyr::group_by(metric, Change, Sex, Phase, window) %>%
          dplyr::group_modify(~ {
            d <- .x
            idx <- which(is.finite(d$p_AUC_raw))
            d$bh_family_n_AUC <- length(idx)
            if (length(idx) > 0) d$p_AUC_BH[idx] <- p.adjust(d$p_AUC_raw[idx], method = "BH")
            d
          }) %>%
          dplyr::ungroup()
      }

      # ---- Export tables ----
      # Parametric terms (Group intercept effects)
      wb_gp <- openxlsx::createWorkbook()
      FN_gp <- "Arial"
      hdr_gp <- openxlsx::createStyle(
        fontName = FN_gp, fontSize = 10, textDecoration = "bold",
        fgFill = "#1F4E79", fontColour = "#FFFFFF", halign = "CENTER", valign = "CENTER",
        border = "Bottom", borderColour = "#173A5E", borderStyle = "medium"
      )
      s_gp_txt <- openxlsx::createStyle(fontName = FN_gp, fontSize = 9, halign = "LEFT",  valign = "CENTER")
      s_gp_num <- openxlsx::createStyle(fontName = FN_gp, fontSize = 9, halign = "RIGHT", valign = "CENTER", numFmt = "0.000")
      s_gp_4dp <- openxlsx::createStyle(fontName = FN_gp, fontSize = 9, halign = "RIGHT", valign = "CENTER", numFmt = "0.0000")
      s_gp_odd <- openxlsx::createStyle(fgFill = "#FFFFFF")
      s_gp_evn <- openxlsx::createStyle(fgFill = "#F7F7F7")
      s_gp_row_sig <- openxlsx::createStyle(fgFill = "#EAF6EC")
      s_gp_row_trd <- openxlsx::createStyle(fgFill = "#FFFDEB")
      apply_gamm_sheet <- function(wb, sht, tbl, p_col_name = "p_value") {
        openxlsx::addWorksheet(wb, sht, gridLines = TRUE)
        openxlsx::writeData(wb, sht, tbl, headerStyle = hdr_gp, rowNames = FALSE)
        openxlsx::setRowHeights(wb, sht, rows = 1, heights = 20)
        nr <- nrow(tbl); nc <- ncol(tbl)
        all_data_rows <- seq_len(nr) + 1L
        odd_rows <- all_data_rows[seq(1L, nr, by = 2L)]
        evn_rows <- if (nr >= 2L) all_data_rows[seq(2L, nr, by = 2L)] else integer(0)
        if (length(odd_rows) > 0)
          openxlsx::addStyle(wb, sht, style = s_gp_odd, rows = odd_rows, cols = 1:nc, gridExpand = TRUE, stack = FALSE)
        if (length(evn_rows) > 0)
          openxlsx::addStyle(wb, sht, style = s_gp_evn, rows = evn_rows, cols = 1:nc, gridExpand = TRUE, stack = FALSE)
        if (p_col_name %in% names(tbl)) {
          pv <- suppressWarnings(as.numeric(tbl[[p_col_name]]))
          sig_rows <- all_data_rows[is.finite(pv) & pv < 0.05]
          trd_rows <- all_data_rows[is.finite(pv) & pv >= 0.05 & pv < 0.10]
          if (length(sig_rows) > 0)
            openxlsx::addStyle(wb, sht, style = s_gp_row_sig, rows = sig_rows, cols = 1:nc, gridExpand = TRUE, stack = TRUE)
          if (length(trd_rows) > 0)
            openxlsx::addStyle(wb, sht, style = s_gp_row_trd, rows = trd_rows, cols = 1:nc, gridExpand = TRUE, stack = TRUE)
        }
        for (j in seq_len(nc)) {
          cn_g <- names(tbl)[j]
          is_p <- grepl("p_value|p\\.value|p_adj|p_BH|p_for_table", cn_g)
          s_a  <- if (!is.numeric(tbl[[j]])) s_gp_txt else if (is_p) s_gp_4dp else s_gp_num
          openxlsx::addStyle(wb, sht, style = s_a,
                             rows = 2:(nr + 1), cols = j, gridExpand = FALSE, stack = TRUE)
        }
        # no cell-level p highlighting; row-level tint is the only significance styling
        openxlsx::freezePane(wb, sht, firstRow = TRUE)
        openxlsx::setColWidths(wb, sht, cols = 1:nc, widths = "auto")
        sig_r <- nr + 3L
        openxlsx::writeData(wb, sht,
          data.frame(x = "*** p\u202f<\u202f0.001   ** p\u202f<\u202f0.01   * p\u202f<\u202f0.05   \u2020 p\u202f<\u202f0.10   ns = not significant | row tint: green p<0.05, yellow p<0.10"),
          startRow = sig_r, startCol = 1, colNames = FALSE)
        openxlsx::addStyle(wb, sht,
          style = openxlsx::createStyle(fontName = FN_gp, fontSize = 8,
                                         fontColour = "#666666", textDecoration = "italic"),
          rows = sig_r, cols = 1, stack = FALSE)
        if (nc >= 2) openxlsx::mergeCells(wb, sht, cols = 1:min(nc, 10), rows = sig_r)
        invisible(wb)
      }

      apply_gamm_sheet(wb_gp, "Parametric_Terms", gamm_param_tbl, p_col_name = "p_for_table")

      # Smooth terms (non-linear effects)
      apply_gamm_sheet(wb_gp, "Smooth_Terms", gamm_smooth_tbl, p_col_name = "p_for_table")

      # EDF summary (wide format: one row per model stratum, one col per smooth)
      edf_wide <- gamm_smooth_tbl %>%
        dplyr::filter(!grepl("AnimalNum", Smooth)) %>%
        dplyr::mutate(Smooth_short = sub(".*:Group", "", Smooth)) %>%
        dplyr::select(metric, Change, Sex, Phase, Smooth_short, edf, Ref_df) %>%
        tidyr::pivot_wider(
          names_from  = Smooth_short,
          values_from = c(edf, Ref_df),
          names_glue  = "{Smooth_short}_{.value}"
        ) %>%
        dplyr::arrange(metric, Change, Sex, Phase)
      apply_gamm_sheet(wb_gp, "EDF_Summary", edf_wide, p_col_name = "NONE")

      # AUC per group (no p column)
      apply_gamm_sheet(wb_gp, "AUC_per_Group", gamm_auc_tbl, p_col_name = "NONE")
      if (nrow(gamm_auc_individual_tbl) > 0) {
        apply_gamm_sheet(wb_gp, "AUC_per_Animal", gamm_auc_individual_tbl, p_col_name = "NONE")
      }

      if (nrow(gamm_diag_tbl) > 0) {
        apply_gamm_sheet(wb_gp, "Diagnostics", gamm_diag_tbl, p_col_name = "kcheck_p_value")
      }
      if (nrow(gamm_acf_summary_tbl) > 0) {
        apply_gamm_sheet(wb_gp, "Residual_ACF_Summary", gamm_acf_summary_tbl, p_col_name = "NONE")
      }
      if (nrow(gamm_acf_tbl) > 0) {
        apply_gamm_sheet(wb_gp, "Residual_ACF_ByLag", gamm_acf_tbl, p_col_name = "NONE")
      }

      # Metadata sheet
      gamm_meta <- tibble::tibble(
        Information_Type = c(
          "Analysis Type", "Model Engine", "Formula Template", "Time Axis",
          "Smooth Type", "Basis Dimension (k)", "Random Effect",
          "Reference Group", "Intercept Interpretation",
          "Smooth Term Interpretation",
          "Residual ACF Diagnostic", "AUC Column", "Pairwise Differences",
          "Statistical Testing", "Publication Notes"
        ),
        Value = c(
          "Generalized Additive Mixed Model (GAMM)",
          paste0("mgcv::bam two-stage AR1: Stage 1 discrete=TRUE (estimate rho), ",
              "Stage 2 AR1 refit with discrete=FALSE if |rho| > ", GAMM_AR1_MIN_RHO,
                 "; select=TRUE; gamma=", GAMM_GAMMA),
          "response ~ Group + s(HalfHourWithinCC0, k=k, bs='tp') + s(HalfHourWithinCC0, by=Group, k=k, bs='tp') + s(AnimalNum, bs='re')",
          if (isTRUE(includePhase)) "HalfHourWithinCC0 is reset within each contiguous Phase block inside each Change." else "HalfHourWithinCC0 is reset within each Change.",
          "Thin-plate regression spline ('tp') — allows data-adaptive non-linear trajectory",
          as.character(GAMM_K),
          "Random intercept per animal via s(AnimalNum, bs='re')",
          "CON (Control group); Group term captures CON vs RES/SUS baseline shifts",
          "Parametric Group terms = mean level difference vs. CON (like lmer intercept)",
          paste0("Smooth terms: s(HalfHourWithinCC0) = shared non-linear time trajectory; ",
                 "s(HalfHourWithinCC0, by=Group) = group-specific deviations from that shared trajectory. ",
                 "edf (effective df) > 1 indicates non-linearity. p-value tests H0: smooth = 0."),
             "Residual autocorrelation is summarized from Pearson residual ACF within AnimalNum x Change x GAMM time block. Rule of thumb: |lag-1 ACF| < 0.2 suggests AR1 is not urgently needed.",
          "Trapezoidal integration of the predicted group trajectory over HalfHourWithinCC0 range",
          paste0("Pairwise difference curves (e.g., SUS-CON predicted difference at each time point). ",
                 "Shaded band = pointwise 95% CI. Where band excludes 0, groups differ at that timepoint."),
             paste0("fREML smoothness selection (fast REML). Raw p-values are retained (p_value_raw). ",
               "BH-corrected p-values (p_BH) are computed within metric\u00d7Change\u00d7Sex\u00d7Phase\u00d7window strata: ",
               "(A) parametric Group terms, (B) smooth terms, (C) pairwise AUC differences."),
          paste0("GAMM is complementary to lmer spline analysis. ",
                 "Key advantage: data-adaptive smoothness (no fixed df required). ",
                 "Compare edf across groups to assess shape complexity differences.")
        )
      )
      apply_gamm_sheet(wb_gp, "Metadata", gamm_meta, p_col_name = "NONE")
      openxlsx::setColWidths(wb_gp, "Metadata", cols = 1:2, widths = c(35, 75))

      openxlsx::saveWorkbook(
        wb_gp,
        file.path(nat_dirs$gamm_tables, paste0("summary_all.xlsx")),
        overwrite = TRUE
      )

      if (nrow(gamm_auc_individual_tbl) > 0) {
        readr::write_csv(
          gamm_auc_individual_tbl,
          file.path(nat_dirs$gamm_tables, paste0("auc_individual_animals_all.csv"))
        )
      }

      # ---- Plot 1: GAMM trajectory (shared smooth + group deviations + 95% CI ribbon) ----
      p_gamm_traj <- ggplot(
        gamm_pred_tbl %>% dplyr::filter(!is.na(predicted)),
        aes(x = HalfHourWithinCC0, y = predicted, color = Group, fill = Group, group = Group)
      ) +
        geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.18, color = NA) +
        geom_line(linewidth = 1.0, alpha = 0.9) +
        scale_color_manual(values = group_colors) +
        scale_fill_manual(values  = group_colors) +
        labs(
          title    = "GAMM trajectories by group",
          subtitle = paste0("mgcv::bam | shared s(HalfHourWithinCC0) + by-Group deviations, bs='tp', k=", GAMM_K, " | ribbon = 95% CI"),
          x        = gamm_time_axis_label, y = "Predicted response", color = "Group", fill = "Group"
        ) +
        theme_classic(base_size = 9) +
        theme(
          legend.position = "top",
          plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          strip.text    = element_text(face = "bold")
        )

      if (includeSex && dplyr::n_distinct(gamm_pred_tbl$Sex) > 1 &&
          includePhase && dplyr::n_distinct(gamm_pred_tbl$Phase) > 1) {
        p_gamm_traj <- p_gamm_traj + facet_grid(metric + Phase ~ Sex, scales = "free_y")
      } else if (includeSex && dplyr::n_distinct(gamm_pred_tbl$Sex) > 1) {
        p_gamm_traj <- p_gamm_traj + facet_grid(metric ~ Sex, scales = "free_y")
      } else if (includePhase && dplyr::n_distinct(gamm_pred_tbl$Phase) > 1) {
        p_gamm_traj <- p_gamm_traj + facet_grid(metric ~ Phase, scales = "free_y")
      } else {
        p_gamm_traj <- p_gamm_traj + facet_wrap(~ metric, scales = "free_y")
      }

      ggsave(
        file.path(nat_dirs$gamm_plots, paste0("trajectory_all.svg")),
        p_gamm_traj, width = 5, height = 6
      )

      # ---- Plot 2: Pairwise difference curves ----
      if (nrow(gamm_pw_tbl) > 0) {
        gamm_pw_plot <- gamm_pw_tbl %>%
          dplyr::filter(!is.na(diff)) %>%
          dplyr::mutate(pair = factor(pair, levels = c("RES-CON", "SUS-CON", "SUS-RES")))

        pair_cols <- pair_colors

        p_gamm_pw <- ggplot(
          gamm_pw_plot,
          aes(x = HalfHourWithinCC0, y = diff, color = pair, fill = pair, group = pair)
        ) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
          geom_ribbon(aes(ymin = diff_lwr, ymax = diff_upr), alpha = 0.16, color = NA) +
          geom_line(linewidth = 0.85) +
          scale_color_manual(values = pair_cols) +
          scale_fill_manual(values  = pair_cols) +
          labs(
            title    = "GAMM pairwise difference curves",
            subtitle = "Predicted Group A − Group B ± pointwise 95% CI | Shading excludes 0 = significant region",
            x        = gamm_time_axis_label, y = "\u0394 Predicted response",
            color = "Pair", fill = "Pair"
          ) +
          theme_classic(base_size = 9) +
          theme(
            legend.position = "top",
            plot.title    = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, size = 8),
            strip.text    = element_text(face = "bold")
          )

        if (includeSex && dplyr::n_distinct(gamm_pw_plot$Sex) > 1 &&
            includePhase && dplyr::n_distinct(gamm_pw_plot$Phase) > 1) {
          p_gamm_pw <- p_gamm_pw + facet_grid(metric + Phase ~ Sex, scales = "free_y")
        } else if (includeSex && dplyr::n_distinct(gamm_pw_plot$Sex) > 1) {
          p_gamm_pw <- p_gamm_pw + facet_grid(metric ~ Sex, scales = "free_y")
        } else if (includePhase && dplyr::n_distinct(gamm_pw_plot$Phase) > 1) {
          p_gamm_pw <- p_gamm_pw + facet_grid(metric ~ Phase, scales = "free_y")
        } else {
          p_gamm_pw <- p_gamm_pw + facet_wrap(pair ~ metric, scales = "free_y")
        }

        ggsave(
          file.path(nat_dirs$gamm_plots, paste0("pairwise_diff_all.svg")),
          p_gamm_pw, width = 4, height = 6  
        )

        # Export pairwise diff table
        if (nrow(gamm_pw_auc_tbl) > 0) {
          write_pub_xlsx(
            gamm_pw_auc_tbl,
            file.path(nat_dirs$gamm_tables, paste0("pairwise_auc_diff_all.xlsx")),
            sheet_name = "gamm_pairwise_auc",
            p_cols = c("p_AUC_BH", "p_AUC_raw")
          )

          p_gamm_auc_contr <- plot_gamm_auc_contrasts_nature(
            gamm_pw_auc_tbl,
            title = "GAMM AUC contrasts",
            subtitle = "Points show integrated group differences; bars show 95% simulation intervals"
          )
          if (!is.null(p_gamm_auc_contr)) {
            n_facets_auc <- dplyr::n_distinct(paste(
              gamm_pw_auc_tbl$metric, gamm_pw_auc_tbl$Change,
              gamm_pw_auc_tbl$Sex, gamm_pw_auc_tbl$Phase,
              sep = "|"
            ))
            ggsave(
              file.path(nat_dirs$gamm_plots, paste0("auc_contrasts_all.svg")),
              p_gamm_auc_contr,
              width = 1.5,
              height = max(2.2, min(9, 1.5 * n_facets_auc)),
              device = "svg"
            )
          }
        }
      }

      # ---- Plot 3: AUC per group ----
      gamm_auc_plot <- gamm_auc_tbl %>%
        dplyr::mutate(Group = factor(Group, levels = c("CON", "RES", "SUS")))

      auc_gam_counts <- gamm_auc_plot %>%
        dplyr::count(metric, Sex, Phase, Group, name = "n")
      single_pt_gam <- nrow(auc_gam_counts) > 0 && max(auc_gam_counts$n, na.rm = TRUE) <= 1

      p_gamm_auc <- ggplot(
        gamm_auc_plot,
        aes(x = Group, y = AUC, color = Group, fill = Group)
      ) +
        scale_color_manual(values = group_colors) +
        scale_fill_manual(values  = group_colors) +
        labs(
          title    = "GAMM AUC by group",
          subtitle = if (single_pt_gam)
            "AUC from GAMM-predicted trajectory (single estimate; point mode)"
          else
            "AUC from GAMM-predicted trajectory over HalfHourWithinCC0 range",
          x = "Group", y = "AUC (predicted response × half-hour bin)"
        ) +
        theme_classic(base_size = 9) +
        theme(
          legend.position = "none",
          strip.text    = element_text(face = "bold"),
          plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 8)
        )

      if (single_pt_gam) {
        p_gamm_auc <- p_gamm_auc +
          geom_line(aes(group = 1), color = "grey70", linewidth = 0.35, alpha = 0.75) +
          geom_point(size = 2.4, alpha = 0.95)
      } else {
        p_gamm_auc <- p_gamm_auc +
          geom_violin(alpha = 0.28, trim = TRUE, linewidth = 0.3) +
          geom_jitter(width = 0.12, size = 1.8, alpha = 0.75) +
          stat_summary(fun = mean, geom = "crossbar", width = 0.35, linewidth = 0.5, color = "black")
      }

      if (includeSex && dplyr::n_distinct(gamm_auc_plot$Sex) > 1 &&
          includePhase && dplyr::n_distinct(gamm_auc_plot$Phase) > 1) {
        p_gamm_auc <- p_gamm_auc + facet_grid(metric + Phase ~ Sex, scales = "free_y")
      } else if (includeSex && dplyr::n_distinct(gamm_auc_plot$Sex) > 1) {
        p_gamm_auc <- p_gamm_auc + facet_grid(metric ~ Sex, scales = "free_y")
      } else if (includePhase && dplyr::n_distinct(gamm_auc_plot$Phase) > 1) {
        p_gamm_auc <- p_gamm_auc + facet_grid(metric ~ Phase, scales = "free_y")
      } else {
        p_gamm_auc <- p_gamm_auc + facet_wrap(~ metric, scales = "free_y")
      }

      ggsave(
        file.path(nat_dirs$gamm_plots, paste0("auc_all.svg")),
        p_gamm_auc, width = 5, height = 3
      )

      if (nrow(gamm_acf_tbl) > 0) {
        p_gamm_acf <- ggplot(
          gamm_acf_tbl,
          aes(x = lag, y = mean_acf, group = 1)
        ) +
          geom_hline(yintercept = 0, color = "grey55", linewidth = 0.4) +
          geom_hline(yintercept = c(-0.2, 0.2), linetype = "dotted", color = "#c08c5c", linewidth = 0.4) +
          geom_line(color = "#1F4E79", linewidth = 0.8) +
          geom_point(color = "#1F4E79", size = 1.8) +
          labs(
            title = "GAMM residual autocorrelation by lag",
            subtitle = "Mean within-animal Pearson residual ACF; dotted lines mark |ACF| = 0.2",
            x = "Lag",
            y = "Mean residual ACF"
          ) +
          theme_classic(base_size = 9) +
          theme(
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, size = 8),
            strip.text = element_text(face = "bold")
          )

        if (includeSex && dplyr::n_distinct(gamm_acf_tbl$Sex) > 1 &&
            includePhase && dplyr::n_distinct(gamm_acf_tbl$Phase) > 1) {
          p_gamm_acf <- p_gamm_acf + facet_grid(metric + Phase ~ Sex + Change, scales = "free_y")
        } else {
          p_gamm_acf <- p_gamm_acf + facet_wrap(~ metric + Change + Sex + Phase, scales = "free_y")
        }

        ggsave(
          file.path(nat_dirs$gamm_plots, paste0("residual_acf_all.svg")),
          p_gamm_acf, width = 5, height = 3
        )
      }

      cat(sprintf("     \u2713 GAMM: %d models | residual ACF diagnostics + 4 plots saved\n", length(gamm_results)))
    } else {
      cat("     - No GAMM models converged\n")
    }

    # Additional GAMM runs: first Change+Active and first/last 2h windows
    run_gamm_variant <- function(data_source, suffix_tag, analysis_label, window_label = "all") {
      if (nrow(data_source) == 0) {
        cat(sprintf("     - GAMM %s skipped: no input rows\n", analysis_label))
        return(invisible(NULL))
      }

      variant_results <- list()
      data_source_gamm <- add_gamm_halfhour_within_cc(data_source, phase_local = isTRUE(includePhase))
      changes_var <- if (includeChange) unique(as.character(data_source$Change)) else "allChanges"
      sexes_var   <- if (includeSex)    unique(as.character(data_source$Sex))    else "allSexes"
      phases_var  <- if (includePhase)  unique(as.character(data_source$Phase))  else "allPhases"

      for (metric in metrics_to_analyze) {
        for (ch in changes_var) {
          for (sx in sexes_var) {
            for (ph in phases_var) {
                dsub_gam <- data_source_gamm %>%
                dplyr::filter(
                  (!includeChange | as.character(Change) == ch),
                  (!includeSex    | as.character(Sex)    == sx),
                  (!includePhase  | as.character(Phase)  == ph),
                  !is.na(.data[[metric]]), !is.na(AnimalNum), !is.na(Group)
                )

              if (window_label %in% c("first2h", "last2h") && nrow(dsub_gam) > 0) {
                t0_phase <- min(dsub_gam$HalfHourElapsed, na.rm = TRUE)
                t1_phase <- max(dsub_gam$HalfHourElapsed, na.rm = TRUE)

                if (window_label == "first2h") {
                  dsub_gam <- dsub_gam %>%
                    dplyr::mutate(rel_phase_h = (HalfHourElapsed - t0_phase) / 2) %>%
                    dplyr::filter(!is.na(rel_phase_h), rel_phase_h >= 0, rel_phase_h <= 2) %>%
                    dplyr::select(-rel_phase_h)
                }
                if (window_label == "last2h") {
                  dsub_gam <- dsub_gam %>%
                    dplyr::mutate(rel_phase_h_from_end = (t1_phase - HalfHourElapsed) / 2) %>%
                    dplyr::filter(!is.na(rel_phase_h_from_end), rel_phase_h_from_end >= 0, rel_phase_h_from_end <= 2) %>%
                    dplyr::select(-rel_phase_h_from_end)
                }
              }

              if (nrow(dsub_gam) < max(min_obs, 20) ||
                  dplyr::n_distinct(dsub_gam$AnimalNum) < min_animals ||
                  dplyr::n_distinct(dsub_gam$Group) < min_groups) next

              dsub_gam <- dsub_gam %>%
                dplyr::mutate(
                  AnimalNum = factor(AnimalNum),
                  Group     = factor(Group, levels = c("CON", "RES", "SUS"))
                ) %>%
                dplyr::filter(!is.na(HalfHourWithinCC0))

              # Required ordering for AR.start boundaries (independent AR1 process per animal x cage change)
              dsub_gam <- dplyr::arrange(dsub_gam, AnimalNum, Change, GammTimeBlock, HalfHourWithinCC0)
              ar_start_vec <- seq_len(nrow(dsub_gam)) == 1L |
                dsub_gam$AnimalNum != dplyr::lag(dsub_gam$AnimalNum) |
                dsub_gam$Change != dplyr::lag(dsub_gam$Change) |
                dsub_gam$GammTimeBlock != dplyr::lag(dsub_gam$GammTimeBlock)

              n_uniq_t_gam <- dplyr::n_distinct(dsub_gam$HalfHourWithinCC0)
              if (n_uniq_t_gam < 8) next
              k_use <- min(GAMM_K, floor(n_uniq_t_gam / 2L))

              gam_formula <- as.formula(paste0(
                metric, " ~ Group + s(HalfHourWithinCC0, k = ", k_use, ", bs = 'tp') + ",
                "s(HalfHourWithinCC0, by = Group, k = ", k_use, ", bs = 'tp') + s(AnimalNum, bs = 're')"
              ))

              m_gam0 <- tryCatch(
                suppressMessages(suppressWarnings(
                  mgcv::bam(
                    formula = gam_formula,
                    data    = dsub_gam,
                    method  = "fREML",
                    family  = gaussian(),
                    discrete = TRUE,
                    select   = GAMM_SELECT,
                    gamma    = GAMM_GAMMA
                  )
                )),
                error = function(e) NULL
              )
              if (is.null(m_gam0)) next

              rho_est <- tryCatch({
                resid0 <- as.numeric(stats::residuals(m_gam0, type = "pearson"))
                if (length(resid0) == nrow(dsub_gam)) {
                  animal_acfs <- tapply(seq_along(resid0), interaction(dsub_gam$AnimalNum, dsub_gam$Change, dsub_gam$GammTimeBlock, drop = TRUE), function(idx) {
                    r <- resid0[idx]
                    r <- r[is.finite(r)]
                    if (length(r) < 3L) return(NA_real_)
                    as.numeric(stats::acf(r, lag.max = 1L, plot = FALSE)$acf[2L])
                  })
                  mean(unlist(animal_acfs), na.rm = TRUE)
                } else 0
              }, error = function(e) 0)
              if (!is.finite(rho_est)) rho_est <- 0
              rho_est <- max(-0.99, min(0.99, rho_est))

              ar1_applied <- FALSE
              rho_used <- 0
              m_gam <- m_gam0
              if (GAMM_AR1 && abs(rho_est) > GAMM_AR1_MIN_RHO) {
                m_gam1 <- tryCatch(
                  suppressMessages(suppressWarnings(
                    mgcv::bam(
                      formula = gam_formula,
                      data    = dsub_gam,
                      method  = "fREML",
                      family  = gaussian(),
                      rho = rho_est,
                      AR.start = ar_start_vec,
                      discrete = FALSE,
                      select   = GAMM_SELECT,
                      gamma    = GAMM_GAMMA
                    )
                  )),
                  error = function(e) NULL
                )
                if (!is.null(m_gam1)) {
                  m_gam <- m_gam1
                  ar1_applied <- TRUE
                  rho_used <- rho_est
                }
              }

              sm <- summary(m_gam)

              ptbl <- as.data.frame(sm$p.table) %>%
                tibble::rownames_to_column("Term") %>%
                tibble::as_tibble()
              if (ncol(ptbl) >= 5) {
                names(ptbl)[2:5] <- c("Estimate", "SE", "t_value", "p_value")
              }
              ptbl <- ptbl %>%
                dplyr::mutate(
                  Contrast = vapply(Term, decode_group_term_label, character(1), grp_factor = dsub_gam$Group),
                  metric = metric,
                  Change = as.character(ch),
                  Sex = as.character(sx),
                  Phase = as.character(ph),
                  window = window_label,
                  AIC = AIC(m_gam),
                  BIC = BIC(m_gam),
                  dev_expl = sm$dev.expl
                ) %>%
                dplyr::relocate(Contrast, .after = Term)

              stbl <- as.data.frame(sm$s.table) %>%
                tibble::rownames_to_column("Smooth") %>%
                tibble::as_tibble()
              if ("Ref.df" %in% names(stbl)) names(stbl)[names(stbl) == "Ref.df"] <- "Ref_df"
              if ("F" %in% names(stbl))      names(stbl)[names(stbl) == "F"] <- "F_stat"
              if ("p-value" %in% names(stbl)) names(stbl)[names(stbl) == "p-value"] <- "p_value"
              stbl <- stbl %>%
                dplyr::mutate(
                  metric = metric,
                  Change = as.character(ch),
                  Sex = as.character(sx),
                  Phase = as.character(ph),
                  window = window_label
                )

              kchk <- tryCatch(mgcv::k.check(m_gam), error = function(e) NULL)
              if (!is.null(kchk)) {
                kdiag_tbl <- as.data.frame(kchk) %>%
                  tibble::rownames_to_column("Smooth") %>%
                  tibble::as_tibble()
                if ("k-index" %in% names(kdiag_tbl)) names(kdiag_tbl)[names(kdiag_tbl) == "k-index"] <- "k_index"
                if ("p-value" %in% names(kdiag_tbl)) names(kdiag_tbl)[names(kdiag_tbl) == "p-value"] <- "kcheck_p_value"
              } else {
                kdiag_tbl <- tibble::tibble(Smooth = NA_character_, k_index = NA_real_, kcheck_p_value = NA_real_)
              }
              kdiag_tbl <- kdiag_tbl %>%
                dplyr::mutate(
                  metric = metric,
                  Change = as.character(ch),
                  Sex = as.character(sx),
                  Phase = as.character(ph),
                  window = window_label,
                  k_target = GAMM_K,
                  k_used = k_use,
                  select_enabled = GAMM_SELECT,
                  gamma_value = GAMM_GAMMA,
                  dev_expl = sm$dev.expl,
                  AIC = AIC(m_gam),
                  BIC = BIC(m_gam),
                  rho_estimated = rho_est,
                  rho_used = rho_used,
                  ar1_applied = ar1_applied
                )

              acf_diag_tbl <- compute_gamm_resid_acf(
                model_obj = m_gam,
                data_obj = dsub_gam,
                metric = metric,
                Change = ch,
                Sex = sx,
                Phase = ph,
                window = window_label
              )

              th_seq <- seq(min(dsub_gam$HalfHourWithinCC0), max(dsub_gam$HalfHourWithinCC0), length.out = GAMM_N_GRID)
              grp_levels <- levels(dsub_gam$Group)

              pred_list_gam <- lapply(grp_levels, function(g) {
                nd <- data.frame(
                  HalfHourWithinCC0 = th_seq,
                  Group     = factor(g, levels = grp_levels),
                  AnimalNum = factor(dsub_gam$AnimalNum[1], levels = levels(dsub_gam$AnimalNum))
                )
                pred_obj <- tryCatch(
                  mgcv::predict.bam(m_gam, newdata = nd, type = "response", exclude = "s(AnimalNum)", se.fit = TRUE),
                  error = function(e) NULL
                )
                if (is.null(pred_obj)) return(NULL)
                data.frame(
                  HalfHourWithinCC0 = th_seq,
                  Group = g,
                  predicted = as.numeric(pred_obj$fit),
                  se = as.numeric(pred_obj$se.fit),
                  lwr = as.numeric(pred_obj$fit) - 1.96 * as.numeric(pred_obj$se.fit),
                  upr = as.numeric(pred_obj$fit) + 1.96 * as.numeric(pred_obj$se.fit),
                  metric = metric,
                  Change = as.character(ch),
                  Sex = as.character(sx),
                  Phase = as.character(ph),
                  window = window_label
                )
              })
              pred_df_gam <- dplyr::bind_rows(pred_list_gam)

              auc_gam <- lapply(grp_levels, function(g) {
                sub <- pred_df_gam %>% dplyr::filter(Group == g)
                tibble::tibble(
                  Group = g,
                  AUC = trapz_auc(sub$HalfHourWithinCC0, sub$predicted),
                  metric = metric,
                  Change = as.character(ch),
                  Sex = as.character(sx),
                  Phase = as.character(ph),
                  window = window_label
                )
              })
              auc_gam_tbl <- dplyr::bind_rows(auc_gam)
              auc_gam_individual_tbl <- extract_gamm_individual_auc(
                model_obj = m_gam,
                data_obj = dsub_gam,
                th_seq = th_seq,
                metric = metric,
                Change = ch,
                Sex = sx,
                Phase = ph,
                window = window_label
              )

              pw_pairs <- list(c("RES", "CON"), c("SUS", "CON"), c("SUS", "RES"))
              pw_rows <- list()
              pw_auc_rows <- list()

              beta_hat <- tryCatch(stats::coef(m_gam), error = function(e) NULL)
              Vp_hat <- tryCatch(mgcv::vcov.gam(m_gam, unconditional = TRUE), error = function(e) {
                tryCatch(vcov(m_gam), error = function(e2) NULL)
              })
              beta_draws <- NULL
              if (!is.null(beta_hat) && !is.null(Vp_hat) &&
                  is.matrix(Vp_hat) && nrow(Vp_hat) == length(beta_hat) &&
                  requireNamespace("MASS", quietly = TRUE)) {
                beta_draws <- tryCatch(
                  MASS::mvrnorm(n = GAMM_AUC_NSIM, mu = beta_hat, Sigma = Vp_hat),
                  error = function(e) NULL
                )
                if (!is.null(beta_draws) && !is.matrix(beta_draws)) {
                  beta_draws <- matrix(beta_draws, nrow = 1)
                }
              }

              for (pp in pw_pairs) {
                g1 <- pp[1]; g0 <- pp[2]
                if (!(g1 %in% grp_levels && g0 %in% grp_levels)) next

                sub1 <- pred_df_gam %>% dplyr::filter(Group == g1)
                sub0 <- pred_df_gam %>% dplyr::filter(Group == g0)
                diff_vec <- sub1$predicted - sub0$predicted
                # NOTE: this initial diff_se assumes independence of group predictions (ignores
                # covariance). It is overwritten by the lpmatrix delta-method below when
                # predict.bam succeeds, which is the statistically correct approach.
                diff_se  <- sqrt(sub1$se^2 + sub0$se^2)

                nd1 <- data.frame(
                  HalfHourWithinCC0 = th_seq,
                  Group = factor(g1, levels = grp_levels),
                  AnimalNum = factor(dsub_gam$AnimalNum[1], levels = levels(dsub_gam$AnimalNum))
                )
                nd0 <- data.frame(
                  HalfHourWithinCC0 = th_seq,
                  Group = factor(g0, levels = grp_levels),
                  AnimalNum = factor(dsub_gam$AnimalNum[1], levels = levels(dsub_gam$AnimalNum))
                )
                X1 <- tryCatch(mgcv::predict.bam(m_gam, newdata = nd1, type = "lpmatrix"), error = function(e) NULL)
                X0 <- tryCatch(mgcv::predict.bam(m_gam, newdata = nd0, type = "lpmatrix"), error = function(e) NULL)
                Xdiff <- NULL
                if (!is.null(X1) && !is.null(X0) && ncol(X1) == ncol(X0)) {
                  Xdiff <- X1 - X0
                  re_cols <- grepl("^s\\(AnimalNum\\)", colnames(Xdiff))
                  if (any(re_cols)) Xdiff[, re_cols] <- 0
                  if (!is.null(beta_hat) && ncol(Xdiff) == length(beta_hat)) {
                    diff_vec <- as.numeric(Xdiff %*% beta_hat)
                  }
                  if (!is.null(Vp_hat) && ncol(Xdiff) == nrow(Vp_hat)) {
                    diff_var <- rowSums((Xdiff %*% Vp_hat) * Xdiff)
                    diff_se <- sqrt(pmax(diff_var, 0))
                  }
                }

                diff_lwr <- diff_vec - 1.96 * diff_se
                diff_upr <- diff_vec + 1.96 * diff_se

                pw_rows[[length(pw_rows) + 1]] <- tibble::tibble(
                  HalfHourWithinCC0 = th_seq,
                  diff = diff_vec,
                  diff_lwr = diff_lwr,
                  diff_upr = diff_upr,
                  pair = paste0(g1, "-", g0),
                  metric = metric,
                  Change = as.character(ch),
                  Sex = as.character(sx),
                  Phase = as.character(ph),
                  window = window_label
                )

                auc_diff <- trapz_auc(th_seq, diff_vec)
                th_range <- max(th_seq, na.rm = TRUE) - min(th_seq, na.rm = TRUE)
                auc_diff_norm <- ifelse(is.finite(th_range) && th_range > 0, auc_diff / th_range, NA_real_)
                frac_sig <- mean(diff_lwr > 0 | diff_upr < 0, na.rm = TRUE)

                auc_ci_low <- NA_real_
                auc_ci_high <- NA_real_
                p_AUC_raw <- NA_real_
                if (!is.null(Xdiff) && !is.null(beta_draws) && ncol(Xdiff) == ncol(beta_draws)) {
                  diff_draw_mat <- Xdiff %*% t(beta_draws)
                  auc_draws <- apply(diff_draw_mat, 2, function(y) trapz_auc(th_seq, y))
                  auc_draws <- auc_draws[is.finite(auc_draws)]
                  if (length(auc_draws) >= 50) {
                    auc_ci_low <- as.numeric(stats::quantile(auc_draws, 0.025, na.rm = TRUE))
                    auc_ci_high <- as.numeric(stats::quantile(auc_draws, 0.975, na.rm = TRUE))
                    p_low <- mean(auc_draws <= 0)
                    p_high <- mean(auc_draws >= 0)
                    p_AUC_raw <- min(1, 2 * min(p_low, p_high))
                  }
                }

                pw_auc_rows[[length(pw_auc_rows) + 1]] <- tibble::tibble(
                  metric = metric,
                  Change = as.character(ch),
                  Sex = as.character(sx),
                  Phase = as.character(ph),
                  window = window_label,
                  pair = paste0(g1, "-", g0),
                  AUC_diff = auc_diff,
                  AUC_diff_norm = auc_diff_norm,
                  AUC_ci_low = auc_ci_low,
                  AUC_ci_high = auc_ci_high,
                  p_AUC_raw = p_AUC_raw,
                  frac_sig = frac_sig
                )
              }

              variant_results[[length(variant_results) + 1]] <- list(
                parametric = ptbl,
                smooths = stbl,
                auc = auc_gam_tbl,
                auc_individual = auc_gam_individual_tbl,
                pairwise = dplyr::bind_rows(pw_rows),
                pairwise_auc = dplyr::bind_rows(pw_auc_rows),
                predictions = pred_df_gam,
                diag = kdiag_tbl,
                acf = acf_diag_tbl$by_lag,
                acf_summary = acf_diag_tbl$summary
              )
            }
          }
        }
      }

      if (length(variant_results) == 0) {
        cat(sprintf("     - GAMM %s skipped: no converged models\n", analysis_label))
        return(invisible(NULL))
      }

      var_param <- dplyr::bind_rows(lapply(variant_results, `[[`, "parametric"))
      var_smooth <- dplyr::bind_rows(lapply(variant_results, `[[`, "smooths"))
      bh_var <- apply_gamm_bh_families(var_param, var_smooth)
      var_param <- bh_var$param
      var_smooth <- bh_var$smooth
      var_auc <- dplyr::bind_rows(lapply(variant_results, `[[`, "auc"))
      var_auc_individual <- dplyr::bind_rows(lapply(variant_results, `[[`, "auc_individual"))
      var_pw <- dplyr::bind_rows(lapply(variant_results, `[[`, "pairwise"))
      var_pw_auc <- dplyr::bind_rows(lapply(variant_results, `[[`, "pairwise_auc"))
      var_diag <- dplyr::bind_rows(lapply(variant_results, `[[`, "diag"))
      var_acf <- dplyr::bind_rows(lapply(variant_results, `[[`, "acf"))
      var_acf_summary <- dplyr::bind_rows(lapply(variant_results, `[[`, "acf_summary"))
      var_pred <- dplyr::bind_rows(lapply(variant_results, `[[`, "predictions")) %>%
        dplyr::mutate(Group = factor(Group, levels = c("CON", "RES", "SUS")))

      if (nrow(var_pw_auc) > 0 && "p_AUC_raw" %in% names(var_pw_auc)) {
        var_pw_auc <- var_pw_auc %>%
          dplyr::mutate(
            metric = if ("metric" %in% names(.)) as.character(metric) else "unknown_metric",
            Change = if ("Change" %in% names(.)) as.character(Change) else "all",
            Sex = if ("Sex" %in% names(.)) as.character(Sex) else "all",
            Phase = if ("Phase" %in% names(.)) as.character(Phase) else "all",
            window = if ("window" %in% names(.)) as.character(window) else "all",
            p_AUC_BH = NA_real_,
            bh_family_n_AUC = NA_integer_,
            bh_family_scope_AUC = "metric|Change|Sex|Phase|window"
          ) %>%
          dplyr::group_by(metric, Change, Sex, Phase, window) %>%
          dplyr::group_modify(~ {
            d <- .x
            idx <- which(is.finite(d$p_AUC_raw))
            d$bh_family_n_AUC <- length(idx)
            if (length(idx) > 0) d$p_AUC_BH[idx] <- p.adjust(d$p_AUC_raw[idx], method = "BH")
            d
          }) %>%
          dplyr::ungroup()
      }

      wb_var <- openxlsx::createWorkbook()
      apply_gamm_sheet(wb_var, "Parametric_Terms", var_param, p_col_name = "p_for_table")
      apply_gamm_sheet(wb_var, "Smooth_Terms", var_smooth, p_col_name = "p_for_table")
      apply_gamm_sheet(wb_var, "AUC_per_Group", var_auc, p_col_name = "NONE")
      if (nrow(var_auc_individual) > 0) {
        apply_gamm_sheet(wb_var, "AUC_per_Animal", var_auc_individual, p_col_name = "NONE")
      }
      if (nrow(var_diag) > 0) {
        apply_gamm_sheet(wb_var, "Diagnostics", var_diag, p_col_name = "kcheck_p_value")
      }
      if (nrow(var_acf_summary) > 0) {
        apply_gamm_sheet(wb_var, "Residual_ACF_Summary", var_acf_summary, p_col_name = "NONE")
      }
      if (nrow(var_acf) > 0) {
        apply_gamm_sheet(wb_var, "Residual_ACF_ByLag", var_acf, p_col_name = "NONE")
      }
      apply_gamm_sheet(wb_var, "Metadata",
        tibble::tibble(
          Information_Type = c("Analysis Type", "Subset", "Window", "Model Engine", "Formula", "Time Axis", "Residual ACF Diagnostic"),
          Value = c("Generalized Additive Mixed Model (GAMM)", analysis_label, window_label, "mgcv::bam",
                    "response ~ Group + s(HalfHourWithinCC0, k=k, bs='tp') + s(HalfHourWithinCC0, by=Group, k=k, bs='tp') + s(AnimalNum, bs='re')",
                    if (isTRUE(includePhase)) "HalfHourWithinCC0 is reset within each contiguous Phase block inside each Change." else "HalfHourWithinCC0 is reset within each Change.",
                    "Mean Pearson residual lag-1 ACF is reported within AnimalNum x Change x GAMM time block; |lag-1 ACF| < 0.2 is a practical threshold for deeming AR1 optional.")
        ),
        p_col_name = "NONE"
      )
      openxlsx::setColWidths(wb_var, "Metadata", cols = 1:2, widths = c(25, 75))

      openxlsx::saveWorkbook(
        wb_var,
        file.path(nat_dirs$gamm_tables, paste0("summary", suffix_tag, ".xlsx")),
        overwrite = TRUE
      )

      if (nrow(var_auc_individual) > 0) {
        readr::write_csv(
          var_auc_individual,
          file.path(nat_dirs$gamm_tables, paste0("auc_individual_animals", suffix_tag, ".csv"))
        )
      }

      if (nrow(var_pw_auc) > 0) {
        write_pub_xlsx(
          var_pw_auc,
          file.path(nat_dirs$gamm_tables, paste0("pairwise_auc_diff", suffix_tag, ".xlsx")),
          sheet_name = "gamm_pairwise_auc",
          p_cols = c("p_AUC_BH", "p_AUC_raw")
        )

        p_var_auc_contr <- plot_gamm_auc_contrasts_nature(
          var_pw_auc,
          title = paste0("GAMM AUC contrasts (", analysis_label, ")"),
          subtitle = paste0("window=", window_label, " | points = \u0394AUC; bars = 95% simulation intervals")
        )
        if (!is.null(p_var_auc_contr)) {
          n_facets_var_auc <- dplyr::n_distinct(paste(
            var_pw_auc$metric, var_pw_auc$Change,
            var_pw_auc$Sex, var_pw_auc$Phase,
            sep = "|"
          ))
          ggsave(
            file.path(nat_dirs$gamm_plots, paste0("auc_contrasts", suffix_tag, ".svg")),
            p_var_auc_contr,
            width = 1.5,
            height = max(2.2, min(9, 1.05 * n_facets_var_auc)),
            device = "svg"
          )
        }
      }

      p_var_traj <- ggplot(
        var_pred %>% dplyr::filter(!is.na(predicted)),
        aes(x = HalfHourWithinCC0, y = predicted, color = Group, fill = Group, group = Group)
      ) +
        geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.18, color = NA) +
        geom_line(linewidth = 1.0, alpha = 0.9) +
        scale_color_manual(values = group_colors) +
        scale_fill_manual(values  = group_colors) +
        labs(
          title = paste0("GAMM trajectories by group (", analysis_label, ")"),
          subtitle = paste0("window=", window_label),
          x = gamm_time_axis_label, y = "Predicted response", color = "Group", fill = "Group"
        ) +
        theme_classic(base_size = 9) +
        theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5))

      if (includeSex && dplyr::n_distinct(var_pred$Sex) > 1 &&
          includePhase && dplyr::n_distinct(var_pred$Phase) > 1) {
        p_var_traj <- p_var_traj + facet_grid(metric + Phase ~ Sex, scales = "free_y")
      } else if (includeSex && dplyr::n_distinct(var_pred$Sex) > 1) {
        p_var_traj <- p_var_traj + facet_grid(metric ~ Sex, scales = "free_y")
      } else if (includePhase && dplyr::n_distinct(var_pred$Phase) > 1) {
        p_var_traj <- p_var_traj + facet_grid(metric ~ Phase, scales = "free_y")
      } else {
        p_var_traj <- p_var_traj + facet_wrap(~ metric, scales = "free_y")
      }

      ggsave(
        file.path(nat_dirs$gamm_plots, paste0("trajectory", suffix_tag, ".svg")),
        p_var_traj, width = 4, height = 4
      )

      if (nrow(var_pw) > 0) {
        pair_cols <- pair_colors
        p_var_pw <- ggplot(
          var_pw %>% dplyr::mutate(pair = factor(pair, levels = c("RES-CON", "SUS-CON", "SUS-RES"))),
          aes(x = HalfHourWithinCC0, y = diff, color = pair, fill = pair, group = pair)
        ) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
          geom_ribbon(aes(ymin = diff_lwr, ymax = diff_upr), alpha = 0.16, color = NA) +
          geom_line(linewidth = 0.85) +
          scale_color_manual(values = pair_cols) +
          scale_fill_manual(values = pair_cols) +
          labs(
            title = paste0("GAMM pairwise differences (", analysis_label, ")"),
            subtitle = paste0("window=", window_label),
            x = gamm_time_axis_label, y = "\u0394 Predicted response", color = "Pair", fill = "Pair"
          ) +
          theme_classic(base_size = 9) +
          theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5))

        if (includeSex && dplyr::n_distinct(var_pw$Sex) > 1 &&
            includePhase && dplyr::n_distinct(var_pw$Phase) > 1) {
          p_var_pw <- p_var_pw + facet_grid(metric + Phase ~ Sex, scales = "free_y")
        } else if (includeSex && dplyr::n_distinct(var_pw$Sex) > 1) {
          p_var_pw <- p_var_pw + facet_grid(metric ~ Sex, scales = "free_y")
        } else if (includePhase && dplyr::n_distinct(var_pw$Phase) > 1) {
          p_var_pw <- p_var_pw + facet_grid(metric ~ Phase, scales = "free_y")
        } else {
          p_var_pw <- p_var_pw + facet_wrap(pair ~ metric, scales = "free_y")
        }

        ggsave(
          file.path(nat_dirs$gamm_plots, paste0("pairwise_diff", suffix_tag, ".svg")),
          p_var_pw, width = 5, height = 3
        )
      }

      p_var_auc <- ggplot(
        var_auc %>% dplyr::mutate(Group = factor(Group, levels = c("CON", "RES", "SUS"))),
        aes(x = Group, y = AUC, color = Group, fill = Group)
      ) +
        scale_color_manual(values = group_colors) +
        scale_fill_manual(values  = group_colors) +
        geom_violin(alpha = 0.28, trim = TRUE, linewidth = 0.3) +
        geom_jitter(width = 0.12, size = 1.8, alpha = 0.75) +
        stat_summary(fun = mean, geom = "crossbar", width = 0.35, linewidth = 0.5, color = "black") +
        labs(
          title = paste0("GAMM AUC by group (", analysis_label, ")"),
          subtitle = paste0("window=", window_label),
          x = "Group", y = "AUC (predicted response \u00d7 half-hour bin)"
        ) +
        theme_classic(base_size = 9) +
        theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5))

      if (includeSex && dplyr::n_distinct(var_auc$Sex) > 1 &&
          includePhase && dplyr::n_distinct(var_auc$Phase) > 1) {
        p_var_auc <- p_var_auc + facet_grid(metric + Phase ~ Sex, scales = "free_y")
      } else if (includeSex && dplyr::n_distinct(var_auc$Sex) > 1) {
        p_var_auc <- p_var_auc + facet_grid(metric ~ Sex, scales = "free_y")
      } else if (includePhase && dplyr::n_distinct(var_auc$Phase) > 1) {
        p_var_auc <- p_var_auc + facet_grid(metric ~ Phase, scales = "free_y")
      } else {
        p_var_auc <- p_var_auc + facet_wrap(~ metric, scales = "free_y")
      }

      ggsave(
        file.path(nat_dirs$gamm_plots, paste0("auc", suffix_tag, ".svg")),
        p_var_auc, width = 5, height = 3
      )

      if (nrow(var_acf) > 0) {
        p_var_acf <- ggplot(
          var_acf,
          aes(x = lag, y = mean_acf, group = 1)
        ) +
          geom_hline(yintercept = 0, color = "grey55", linewidth = 0.4) +
          geom_hline(yintercept = c(-0.2, 0.2), linetype = "dotted", color = "#c08c5c", linewidth = 0.4) +
          geom_line(color = "#1F4E79", linewidth = 0.8) +
          geom_point(color = "#1F4E79", size = 1.8) +
          labs(
            title = paste0("GAMM residual autocorrelation (", analysis_label, ")"),
            subtitle = paste0("window=", window_label, " | mean within-animal Pearson residual ACF"),
            x = "Lag",
            y = "Mean residual ACF"
          ) +
          theme_classic(base_size = 9) +
          theme(
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, size = 8),
            strip.text = element_text(face = "bold")
          )

        if (includeSex && dplyr::n_distinct(var_acf$Sex) > 1 &&
            includePhase && dplyr::n_distinct(var_acf$Phase) > 1) {
          p_var_acf <- p_var_acf + facet_grid(metric + Phase ~ Sex + Change, scales = "free_y")
        } else {
          p_var_acf <- p_var_acf + facet_wrap(~ metric + Change + Sex + Phase, scales = "free_y")
        }

        ggsave(
          file.path(nat_dirs$gamm_plots, paste0("residual_acf", suffix_tag, ".svg")),
          p_var_acf, width = 5, height = 3
        )
      }

      cat(sprintf("     \u2713 GAMM %s: %d models | residual ACF diagnostics + 4 plots + tables saved\n", analysis_label, length(variant_results)))
      invisible(NULL)
    }

    change_levels_gam <- levels(droplevels(data_filtered_agg$Change))
    if (length(change_levels_gam) == 0) {
      change_levels_gam <- sort(unique(as.character(data_filtered_agg$Change)))
    }
    change_ord_num_gam <- suppressWarnings(as.numeric(stringr::str_extract(change_levels_gam, "\\d+")))
    first_change_level_gam <- if (length(change_levels_gam) == 0) {
      NA_character_
    } else if (all(!is.na(change_ord_num_gam))) {
      change_levels_gam[order(change_ord_num_gam)][1]
    } else {
      sort(change_levels_gam)[1]
    }

    if (!is.na(first_change_level_gam)) {
      data_first_change_active_gam <- data_filtered_agg %>%
        dplyr::filter(
          as.character(Change) == first_change_level_gam,
          as.character(Phase) == "Active"
        )
      run_gamm_variant(
        data_source = data_first_change_active_gam,
        suffix_tag = "_firstChangeActive",
        analysis_label = paste0("first Change (", first_change_level_gam, ") + Active phase"),
        window_label = "all"
      )
    } else {
      cat("     - GAMM firstChangeActive skipped: could not determine first change level\n")
    }

    run_gamm_variant(
      data_source = data_filtered_agg,
      suffix_tag = "_first2h",
      analysis_label = "first 2h within each phase",
      window_label = "first2h"
    )
    run_gamm_variant(
      data_source = data_filtered_agg,
      suffix_tag = "_last2h",
      analysis_label = "last 2h within each phase",
      window_label = "last2h"
    )

    saveRDS(
      list(
        gamm_results = if (exists("gamm_results")) gamm_results else list(),
        gamm_param = if (exists("gamm_param_tbl")) gamm_param_tbl else tibble::tibble(),
        gamm_smooth = if (exists("gamm_smooth_tbl")) gamm_smooth_tbl else tibble::tibble(),
        gamm_auc = if (exists("gamm_auc_tbl")) gamm_auc_tbl else tibble::tibble(),
        gamm_pairwise = if (exists("gamm_pw_tbl")) gamm_pw_tbl else tibble::tibble(),
        gamm_pairwise_auc = if (exists("gamm_pw_auc_tbl")) gamm_pw_auc_tbl else tibble::tibble(),
        gamm_diag = if (exists("gamm_diag_tbl")) gamm_diag_tbl else tibble::tibble(),
        gamm_acf = if (exists("gamm_acf_tbl")) gamm_acf_tbl else tibble::tibble(),
        gamm_acf_summary = if (exists("gamm_acf_summary_tbl")) gamm_acf_summary_tbl else tibble::tibble()
      ),
      file.path(nat_dirs$artifacts, paste0("gamm_snapshot_", run_scope_current, ".rds")),
      compress = "xz"
    )
  }

  saveRDS(
    list(
      run_scope = run_scope_current,
      integrated_full = integrated_full,
      rep_models = rep_models
    ),
    file.path(nat_dirs$artifacts, paste0("analyses_snapshot_", run_scope_current, ".rds")),
    compress = "xz"
  )

  cat(sprintf("\n=== Analyses complete. Results saved to: %s ===\n", analyses_dir))
}

# -------------------------------------------------
# Optional: Data Quality Diagnostics
# -------------------------------------------------
# Set RUN_DIAGNOSTICS to TRUE to verify data merge integrity
RUN_DIAGNOSTICS <- FALSE

if (RUN_DIAGNOSTICS) {
  library(dplyr)
  library(purrr)
  
  cat("\n === DATA QUALITY DIAGNOSTICS ===\n")
  
  # 1. Check Animal ID consistency across primary and activity datasets
  id_check <- df_primary %>%
    group_by(Batch) %>%
    summarise(
      Primary_IDs = paste(sort(unique(AnimalNum))[1:3], collapse = ", "),
      Activity_IDs = paste(sort(unique(df_activity$AnimalNum[df_activity$Batch == first(Batch)]))[1:3], collapse = ", "),
      Count_Primary = n_distinct(AnimalNum),
      Count_Activity = n_distinct(df_activity$AnimalNum[df_activity$Batch == first(Batch)]),
      All_Match = all(unique(AnimalNum) %in% unique(df_activity$AnimalNum[df_activity$Batch == first(Batch)]))
    )
  
  cat("\nAnimal ID Summary by Batch:\n")
  print(id_check)
  
  # 2. Identify any animals in primary data but not in activity data
  missing_animals <- setdiff(unique(df_primary$AnimalNum), unique(df_activity$AnimalNum))
  if (length(missing_animals) > 0) {
    cat("\n⚠️  WARNING: The following IDs exist in Primary data but NOT in Activity data:\n")
    print(missing_animals)
  } else {
    cat("\n✅ PASS: All Animal IDs match perfectly between datasets.\n")
  }
  
  # 3. Find animals with missing ActivityIndex after join
  ghost_animals <- data_filtered_agg %>%
    filter(is.na(ActivityIndex)) %>%
    group_by(Batch, AnimalNum) %>%
    summarise(Missing_Rows = n(), .groups = 'drop')
  
  if (nrow(ghost_animals) > 0) {
    cat("\n⚠️  WARNING: Animals with missing ActivityIndex values:\n")
    print(ghost_animals)
  } else {
    cat("\n✅ PASS: No missing ActivityIndex values found.\n")
  }
  
  cat("\n === END DIAGNOSTICS ===\n")
}

# Get all unique animals from the activity file
activity_registry <- df_activity %>% 
  select(Batch, AnimalNum) %>% 
  distinct() %>% 
  mutate(In_Activity_File = TRUE)

# Join this registry to your ghost list
diagnostic_report <- ghost_animals %>%
  left_join(activity_registry, by = c("Batch", "AnimalNum"))

cat("\nDiagnostic Report:\n")
print(diagnostic_report)

# Replace '3' with a known missing ID from your ghost list
problem_id <- "3" 
problem_batch <- "B1"

cat("\nChecking alignment for Animal:", problem_id, "in Batch:", problem_batch, "\n")

# Check time range in Primary
primary_range <- df_primary %>% 
  filter(Batch == problem_batch, AnimalNum == problem_id) %>% 
  pull(HalfHourElapsed) %>% range()

# Check time range in Activity
activity_range <- df_activity %>% 
  filter(Batch == problem_batch, AnimalNum == problem_id) %>% 
  pull(HalfHourElapsed) %>% range()

cat("Primary Time Range: ", primary_range[1], "to", primary_range[2], "\n")
cat("Activity Time Range:", activity_range[1], "to", activity_range[2], "\n")

# Check Phase naming
cat("Primary Phases: ", unique(df_primary$Phase), "\n")
cat("Activity Phases:", unique(df_activity$Phase), "\n")

# Define the problem animals found in your report
problem_ids <- c("4", "OR426")

# Check metadata in Primary
meta_primary <- df_primary %>%
  filter(AnimalNum %in% problem_ids) %>%
  group_by(Batch, AnimalNum, Change, Sex, Group) %>%
  summarise(Rows_Primary = n(), .groups = 'drop')

# Check metadata in Activity
meta_activity <- df_activity %>%
  filter(AnimalNum %in% problem_ids) %>%
  group_by(Batch, AnimalNum, Change, Sex, Group) %>%
  summarise(Rows_Activity = n(), .groups = 'drop')

# View them together
cat("Metadata in Primary Dataset:\n")
print(meta_primary)

cat("\nMetadata in Activity Dataset:\n")
print(meta_activity)
