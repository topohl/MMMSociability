# -------------------------------------------------
# correlate_gamm_auc_CombZ_two_auc_files_nature_FIXED_SCALES_PADDED.R
# -------------------------------------------------

library(dplyr)
library(readr)
library(readxl)
library(ggplot2)
library(broom)
library(tidyr)
library(stringr)
library(purrr)
library(forcats)
library(patchwork)

# -------------------------------------------------
# 1) PATHS
# -------------------------------------------------

auc_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/statistics/gamm_new/analyses/gamm/tables"

auc_files <- list(
  firstChangeActive = file.path(auc_dir, "auc_individual_animals_firstChangeActive.csv"),
  allPhases = file.path(auc_dir, "auc_individual_animals_all.csv")
)

classic_data <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/SIS_Analysis/E9_Behavior_Data.xlsx"

base_saving_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability/statistics/correlation_gamm_auc_CombZ_sex_specific_nature/"

dir.create(base_saving_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------
# 2) STYLE
# -------------------------------------------------

theme_nature <- function(base_size = 7) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = "black", family = "sans"),
      axis.text = element_text(size = base_size - 1, color = "black"),
      axis.title = element_text(size = base_size, face = "bold", color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.25),
      axis.ticks = element_line(color = "black", linewidth = 0.25),
      axis.ticks.length = grid::unit(1.0, "mm"),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.justification = "left",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 1),
      legend.key.size = grid::unit(3.0, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 1, hjust = 0, color = "grey25"),
      plot.margin = margin(5, 7, 5, 8)
    )
}

group_colors <- c(
  "CON" = "#3E3C6F",
  "RES" = "#9E9A92",
  "SUS" = "#D7303F"
)

group_fills <- scales::alpha(group_colors, 0.42)

sex_labels <- c("f" = "Females", "m" = "Males")

# -------------------------------------------------
# 3) HELPERS
# -------------------------------------------------

save_svg <- function(plot, filename, width = 2.5, height = 2.5) {
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300
  )
}

safe_cor <- function(data, x, y) {
  d <- data %>%
    filter(is.finite(.data[[x]]), is.finite(.data[[y]]))

  if (
    nrow(d) < 4 ||
    sd(d[[x]], na.rm = TRUE) == 0 ||
    sd(d[[y]], na.rm = TRUE) == 0
  ) {
    return(tibble(
      estimate = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      p.value = NA_real_,
      n = nrow(d)
    ))
  }

  ct <- suppressWarnings(cor.test(d[[x]], d[[y]], method = "pearson"))

  tibble(
    estimate = unname(ct$estimate),
    conf.low = unname(ct$conf.int[1]),
    conf.high = unname(ct$conf.int[2]),
    p.value = ct$p.value,
    n = nrow(d)
  )
}

make_stat_text <- function(data, x, y) {
  s <- safe_cor(data, x, y)

  if (!is.finite(s$estimate)) {
    return(sprintf("n = %d; correlation not estimated", s$n))
  }

  sprintf(
    "r = %.2f, 95%% CI %.2f to %.2f, p = %.3g, n = %d",
    s$estimate,
    s$conf.low,
    s$conf.high,
    s$p.value,
    s$n
  )
}

get_limits_padded <- function(data, variable, q = c(0.02, 0.98), pad_frac = 0.25, min_pad = 0.25) {
  v <- data[[variable]]
  v <- v[is.finite(v)]

  if (length(v) == 0) return(NULL)

  if (length(v) < 5) {
    r <- range(v, na.rm = TRUE)
  } else {
    r <- as.numeric(quantile(v, probs = q, na.rm = TRUE))
  }

  d <- diff(r)

  if (!is.finite(d) || d == 0) {
    d <- sd(v, na.rm = TRUE)
  }

  if (!is.finite(d) || d == 0) {
    d <- 1
  }

  pad <- max(d * pad_frac, min_pad)

  c(r[1] - pad, r[2] + pad)
}

get_limits_robust <- function(data, variable, pad = 0.08, q = c(0.02, 0.98)) {
  v <- data[[variable]]
  v <- v[is.finite(v)]

  if (length(v) == 0) return(NULL)
  if (length(v) < 5) return(range(v, na.rm = TRUE))

  qs <- quantile(v, probs = q, na.rm = TRUE)
  d <- diff(qs)

  if (!is.finite(d) || d == 0) {
    d <- sd(v, na.rm = TRUE)
  }

  if (!is.finite(d) || d == 0) {
    d <- 1
  }

  c(qs[1] - d * pad, qs[2] + d * pad)
}

format_p <- function(p) {
  case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

fit_lm_safe <- function(data, formula, model_name) {
  vars <- all.vars(formula)

  d <- data %>%
    drop_na(any_of(vars)) %>%
    mutate(across(where(is.factor), droplevels))

  if (nrow(d) < 6) {
    return(tibble(
      model = model_name,
      formula_used = NA_character_,
      term = NA_character_,
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      n = nrow(d),
      r.squared = NA_real_,
      adj.r.squared = NA_real_
    ))
  }

  one_level_factors <- vars[
    vars %in% names(d) &
      vapply(
        vars,
        function(v) is.factor(d[[v]]) && nlevels(d[[v]]) < 2,
        logical(1)
      )
  ]

  rhs_terms <- attr(terms(formula), "term.labels")
  response <- all.vars(formula)[1]

  if (length(one_level_factors) > 0) {
    rhs_terms <- rhs_terms[
      !vapply(
        rhs_terms,
        function(term) {
          any(vapply(
            one_level_factors,
            function(fac) str_detect(term, paste0("(^|[:*])", fac, "($|[:*])")),
            logical(1)
          ))
        },
        logical(1)
      )
    ]
  }

  if (length(rhs_terms) == 0) {
    return(tibble(
      model = model_name,
      formula_used = NA_character_,
      term = NA_character_,
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      n = nrow(d),
      r.squared = NA_real_,
      adj.r.squared = NA_real_
    ))
  }

  formula_fixed <- as.formula(paste(response, "~", paste(rhs_terms, collapse = " + ")))

  fit <- tryCatch(
    lm(formula_fixed, data = d),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    return(tibble(
      model = model_name,
      formula_used = deparse(formula_fixed),
      term = NA_character_,
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_,
      n = nrow(d),
      r.squared = NA_real_,
      adj.r.squared = NA_real_
    ))
  }

  gl <- broom::glance(fit)

  broom::tidy(fit, conf.int = TRUE) %>%
    mutate(
      model = model_name,
      formula_used = deparse(formula_fixed),
      n = nrow(d),
      r.squared = gl$r.squared,
      adj.r.squared = gl$adj.r.squared,
      .before = 1
    )
}

# -------------------------------------------------
# 4) LOAD PHYSIOLOGY ONCE
# -------------------------------------------------

if (!file.exists(classic_data)) {
  stop(sprintf("Classic data file not found: %s", classic_data))
}

zscore_data <- read_excel(classic_data, sheet = "zScore")

df_physio <- zscore_data %>%
  rename(
    AnimalNum = ID,
    PhysioGroup = Group
  ) %>%
  mutate(
    AnimalNum = as.character(AnimalNum),
    Sex = as.character(Sex),
    Sex = factor(Sex, levels = c("f", "m")),
    delta_cort = as.numeric(delta_cort),
    sucrose_pref = as.numeric(sucrose_pref),
    CombZ = as.numeric(CombZ)
  ) %>%
  select(AnimalNum, Sex, delta_cort, sucrose_pref, CombZ)

# -------------------------------------------------
# 5) MAIN ANALYSIS FUNCTION
# -------------------------------------------------

run_auc_analysis <- function(auc_file, analysis_name) {

  cat("\n-------------------------------------------------\n")
  cat("Running analysis:", analysis_name, "\n")
  cat("AUC file:", auc_file, "\n")
  cat("-------------------------------------------------\n")

  if (!file.exists(auc_file)) {
    stop(sprintf("AUC file not found: %s", auc_file))
  }

  saving_dir <- file.path(base_saving_dir, analysis_name)
  table_dir <- file.path(saving_dir, "tables")
  fig_dir <- file.path(saving_dir, "figures")
  supp_dir <- file.path(saving_dir, "supplementary")
  qc_dir <- file.path(saving_dir, "qc")

  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

  auc_data <- read_csv(auc_file, show_col_types = FALSE)

  auc_movement <- auc_data %>%
    mutate(
      AnimalNum = as.character(AnimalNum),
      Batch = as.character(Batch),
      Group = toupper(as.character(Group)),
      Group = factor(Group, levels = c("CON", "RES", "SUS")),
      Sex = as.character(Sex),
      Sex = factor(Sex, levels = c("f", "m")),
      Phase = as.character(Phase),
      Phase = factor(Phase, levels = c("Active", "Inactive")),
      Change = as.character(Change),
      window = as.character(window),
      metric = as.character(metric),
      MeanPredictedMovement = as.numeric(AUC_norm)
    ) %>%
    filter(
      metric == "Movement",
      is.finite(MeanPredictedMovement),
      !is.na(Group),
      !is.na(Sex),
      !is.na(Phase)
    ) %>%
    mutate(across(where(is.factor), droplevels))

  control_stats_batch <- auc_movement %>%
    filter(Group == "CON") %>%
    group_by(Batch, Sex, Phase, Change, window) %>%
    summarise(
      con_mean_batch = mean(MeanPredictedMovement, na.rm = TRUE),
      con_sd_batch = sd(MeanPredictedMovement, na.rm = TRUE),
      con_n_batch = sum(is.finite(MeanPredictedMovement)),
      .groups = "drop"
    )

  control_stats_global <- auc_movement %>%
    filter(Group == "CON") %>%
    group_by(Sex, Phase, Change, window) %>%
    summarise(
      con_mean_global = mean(MeanPredictedMovement, na.rm = TRUE),
      con_sd_global = sd(MeanPredictedMovement, na.rm = TRUE),
      con_n_global = sum(is.finite(MeanPredictedMovement)),
      .groups = "drop"
    )

  auc_z <- auc_movement %>%
    left_join(
      control_stats_batch,
      by = c("Batch", "Sex", "Phase", "Change", "window")
    ) %>%
    left_join(
      control_stats_global,
      by = c("Sex", "Phase", "Change", "window")
    ) %>%
    mutate(
      use_batch_reference = !is.na(con_sd_batch) &
        con_sd_batch > 0 &
        con_n_batch >= 3,
      CON_mean_used = if_else(use_batch_reference, con_mean_batch, con_mean_global),
      CON_sd_used = if_else(use_batch_reference, con_sd_batch, con_sd_global),
      Movement_AUC_z_vs_CON = if_else(
        is.na(CON_sd_used) | CON_sd_used == 0,
        NA_real_,
        (MeanPredictedMovement - CON_mean_used) / CON_sd_used
      )
    ) %>%
    filter(is.finite(Movement_AUC_z_vs_CON)) %>%
    mutate(across(where(is.factor), droplevels))

  df_trait <- auc_z %>%
    left_join(df_physio, by = c("AnimalNum", "Sex")) %>%
    filter(
      is.finite(Movement_AUC_z_vs_CON),
      is.finite(CombZ)
    ) %>%
    mutate(across(where(is.factor), droplevels))

  write_csv(auc_z, file.path(table_dir, "gamm_auc_z_scored_vs_control.csv"))
  write_csv(df_trait, file.path(table_dir, "merged_full_table.csv"))

  qc_overview <- df_trait %>%
    count(Sex, Phase, Group, Change, window, name = "n")

  write_csv(qc_overview, file.path(qc_dir, "qc_counts_by_sex_phase_group_change_window.csv"))

  scale_summary <- df_trait %>%
    summarise(
      Movement_min = min(Movement_AUC_z_vs_CON, na.rm = TRUE),
      Movement_q02 = quantile(Movement_AUC_z_vs_CON, 0.02, na.rm = TRUE),
      Movement_q05 = quantile(Movement_AUC_z_vs_CON, 0.05, na.rm = TRUE),
      Movement_median = median(Movement_AUC_z_vs_CON, na.rm = TRUE),
      Movement_q95 = quantile(Movement_AUC_z_vs_CON, 0.95, na.rm = TRUE),
      Movement_q98 = quantile(Movement_AUC_z_vs_CON, 0.98, na.rm = TRUE),
      Movement_max = max(Movement_AUC_z_vs_CON, na.rm = TRUE),
      CombZ_min = min(CombZ, na.rm = TRUE),
      CombZ_q02 = quantile(CombZ, 0.02, na.rm = TRUE),
      CombZ_q98 = quantile(CombZ, 0.98, na.rm = TRUE),
      CombZ_max = max(CombZ, na.rm = TRUE)
    )

  write_csv(scale_summary, file.path(qc_dir, "qc_axis_scale_summary.csv"))

  summary_by_sex_phase_group <- df_trait %>%
    group_by(Sex, Phase, Group) %>%
    summarise(
      n = n(),
      Movement_mean = mean(Movement_AUC_z_vs_CON, na.rm = TRUE),
      Movement_sd = sd(Movement_AUC_z_vs_CON, na.rm = TRUE),
      CombZ_mean = mean(CombZ, na.rm = TRUE),
      CombZ_sd = sd(CombZ, na.rm = TRUE),
      delta_cort_mean = mean(delta_cort, na.rm = TRUE),
      delta_cort_sd = sd(delta_cort, na.rm = TRUE),
      sucrose_pref_mean = mean(sucrose_pref, na.rm = TRUE),
      sucrose_pref_sd = sd(sucrose_pref, na.rm = TRUE),
      .groups = "drop"
    )

  write_csv(summary_by_sex_phase_group, file.path(table_dir, "summary_by_sex_phase_group.csv"))

  # -------------------------------------------------
  # Models
  # -------------------------------------------------

  models_primary <- df_trait %>%
    group_by(Sex, Phase) %>%
    group_split() %>%
    map_dfr(function(d) {
      sex_i <- as.character(unique(d$Sex))
      phase_i <- as.character(unique(d$Phase))

      bind_rows(
        fit_lm_safe(
          d,
          CombZ ~ Movement_AUC_z_vs_CON + Group,
          paste0(sex_i, "_", phase_i, "_primary_CombZ_movement_plus_group")
        ),
        fit_lm_safe(
          d,
          CombZ ~ Movement_AUC_z_vs_CON * Group,
          paste0(sex_i, "_", phase_i, "_group_interaction")
        )
      )
    })

  pooled_models <- bind_rows(
    fit_lm_safe(
      df_trait,
      CombZ ~ Movement_AUC_z_vs_CON * Sex + Group,
      "pooled_sex_interaction"
    ),
    fit_lm_safe(
      df_trait,
      CombZ ~ Movement_AUC_z_vs_CON * Phase + Group + Sex,
      "pooled_phase_interaction_sex_adjusted"
    ),
    fit_lm_safe(
      df_trait,
      CombZ ~ Movement_AUC_z_vs_CON + Sex + Group + Batch,
      "pooled_batch_sensitivity"
    )
  )

  all_models <- bind_rows(models_primary, pooled_models)

  write_csv(
    all_models,
    file.path(table_dir, "PRIMARY_CombZ_models.csv")
  )

  # -------------------------------------------------
  # Correlations
  # -------------------------------------------------

  cor_sex_phase <- df_trait %>%
    group_by(Sex, Phase) %>%
    group_modify(~ safe_cor(.x, "Movement_AUC_z_vs_CON", "CombZ")) %>%
    ungroup() %>%
    mutate(
      SexLabel = recode(as.character(Sex), !!!sex_labels),
      Label = paste0(SexLabel, ", ", Phase),
      p.adj = p.adjust(p.value, method = "BH")
    )

  write_csv(
    cor_sex_phase,
    file.path(table_dir, "PRIMARY_CombZ_correlations_by_sex_phase.csv")
  )

  # -------------------------------------------------
  # Distribution scale only
  # -------------------------------------------------

  x_lim_dist <- get_limits_padded(
    df_trait,
    "Movement_AUC_z_vs_CON",
    q = c(0.01, 0.99),
    pad_frac = 0.18,
    min_pad = 0.35
  )

  # -------------------------------------------------
  # Plot functions
  # -------------------------------------------------

  plot_distribution <- function(data, sex_label, phase_label, filename) {
    d <- data %>%
      filter(is.finite(Movement_AUC_z_vs_CON), !is.na(Group))

    if (nrow(d) < 4) return(NULL)

    p <- ggplot(d, aes(x = Group, y = Movement_AUC_z_vs_CON)) +
      geom_hline(yintercept = 0, color = "grey55", linewidth = 0.25) +
      geom_violin(
        aes(fill = Group),
        width = 0.72,
        alpha = 0.50,
        linewidth = 0.25,
        trim = FALSE,
        color = NA
      ) +
      geom_boxplot(
        aes(color = Group),
        width = 0.16,
        outlier.shape = NA,
        alpha = 0.85,
        linewidth = 0.28,
        fill = "white"
      ) +
      geom_point(
        aes(color = Group),
        position = position_jitter(width = 0.075, height = 0),
        size = 1.65,
        alpha = 0.82
      ) +
      stat_summary(
        fun = mean,
        geom = "point",
        shape = 95,
        size = 5,
        linewidth = 0.7,
        color = "black"
      ) +
      scale_fill_manual(values = group_fills, drop = FALSE) +
      scale_color_manual(values = group_colors, drop = FALSE) +
      coord_cartesian(ylim = x_lim_dist, clip = "on") +
      labs(
        title = paste0(sex_label, ": ", phase_label),
        subtitle = "GAMM-derived movement AUC",
        x = NULL,
        y = "Movement AUC, Z vs CON"
      ) +
      theme_nature() +
      theme(legend.position = "none")

    save_svg(p, file.path(fig_dir, filename), width = 2.25, height = 2.35)
    p
  }

  plot_combz <- function(data, sex_label, phase_label, filename) {
    d <- data %>%
      filter(
        is.finite(Movement_AUC_z_vs_CON),
        is.finite(CombZ),
        !is.na(Group)
      )

    if (nrow(d) < 4) return(NULL)

    x_lim_local <- get_limits_padded(
      d,
      "Movement_AUC_z_vs_CON",
      q = c(0.02, 0.98),
      pad_frac = 0.25,
      min_pad = 0.35
    )

    y_lim_local <- get_limits_padded(
      d,
      "CombZ",
      q = c(0.02, 0.98),
      pad_frac = 0.22,
      min_pad = 0.25
    )

    stat_text <- make_stat_text(d, "Movement_AUC_z_vs_CON", "CombZ")

    p <- ggplot(d, aes(x = Movement_AUC_z_vs_CON, y = CombZ)) +
      geom_smooth(
        method = "lm",
        se = TRUE,
        color = "black",
        fill = "grey88",
        linewidth = 0.42,
        alpha = 0.40,
        fullrange = FALSE
      ) +
      geom_point(
        aes(color = Group),
        size = 2.00,
        alpha = 0.90
      ) +
      scale_color_manual(values = group_colors, drop = FALSE) +
      coord_cartesian(
        xlim = x_lim_local,
        ylim = y_lim_local,
        clip = "on"
      ) +
      labs(
        title = paste0(sex_label, ": movement tracks stress burden"),
        subtitle = paste0(phase_label, "; ", stat_text),
        x = "Movement AUC, Z vs CON",
        y = "Composite stress burden (CombZ)"
      ) +
      theme_nature() +
      theme(
        plot.margin = margin(5, 7, 5, 8)
      )

    save_svg(p, file.path(fig_dir, filename), width = 2.75, height = 2.45)
    p
  }

  plot_group_slope <- function(data, sex_label, phase_label, filename) {
    d <- data %>%
      filter(
        is.finite(Movement_AUC_z_vs_CON),
        is.finite(CombZ),
        !is.na(Group)
      )

    if (nrow(d) < 4) return(NULL)

    x_lim_local <- get_limits_padded(
      d,
      "Movement_AUC_z_vs_CON",
      q = c(0.02, 0.98),
      pad_frac = 0.25,
      min_pad = 0.35
    )

    y_lim_local <- get_limits_padded(
      d,
      "CombZ",
      q = c(0.02, 0.98),
      pad_frac = 0.22,
      min_pad = 0.25
    )

    p <- ggplot(d, aes(x = Movement_AUC_z_vs_CON, y = CombZ)) +
      geom_smooth(
        method = "lm",
        se = TRUE,
        color = "black",
        fill = "grey90",
        linewidth = 0.38,
        alpha = 0.28,
        fullrange = FALSE
      ) +
      geom_smooth(
        aes(color = Group),
        method = "lm",
        se = FALSE,
        linewidth = 0.35,
        alpha = 0.80,
        fullrange = FALSE
      ) +
      geom_point(
        aes(color = Group),
        size = 1.95,
        alpha = 0.86
      ) +
      scale_color_manual(values = group_colors, drop = FALSE) +
      coord_cartesian(
        xlim = x_lim_local,
        ylim = y_lim_local,
        clip = "on"
      ) +
      labs(
        title = paste0(sex_label, ": group-slope diagnostic"),
        subtitle = phase_label,
        x = "Movement AUC, Z vs CON",
        y = "Composite stress burden (CombZ)"
      ) +
      theme_nature() +
      theme(
        plot.margin = margin(5, 7, 5, 8)
      )

    save_svg(p, file.path(supp_dir, filename), width = 2.75, height = 2.45)
    p
  }

  plot_component <- function(data, sex_label, phase_label, y, y_label, filename) {
    d <- data %>%
      filter(
        is.finite(Movement_AUC_z_vs_CON),
        is.finite(.data[[y]]),
        !is.na(Group)
      )

    if (nrow(d) < 4) return(NULL)

    x_lim_local <- get_limits_padded(
      d,
      "Movement_AUC_z_vs_CON",
      q = c(0.02, 0.98),
      pad_frac = 0.25,
      min_pad = 0.35
    )

    y_lim_component <- get_limits_padded(
      d,
      y,
      q = c(0.02, 0.98),
      pad_frac = 0.22,
      min_pad = 0.25
    )

    stat_text <- make_stat_text(d, "Movement_AUC_z_vs_CON", y)

    p <- ggplot(d, aes(x = Movement_AUC_z_vs_CON, y = .data[[y]])) +
      geom_smooth(
        method = "lm",
        se = TRUE,
        color = "black",
        fill = "grey90",
        linewidth = 0.36,
        alpha = 0.32,
        fullrange = FALSE
      ) +
      geom_point(
        aes(color = Group),
        size = 1.95,
        alpha = 0.84
      ) +
      scale_color_manual(values = group_colors, drop = FALSE) +
      coord_cartesian(
        xlim = x_lim_local,
        ylim = y_lim_component,
        clip = "on"
      ) +
      labs(
        title = paste0(sex_label, ": ", y_label),
        subtitle = paste0(phase_label, "; ", stat_text),
        x = "Movement AUC, Z vs CON",
        y = y_label
      ) +
      theme_nature() +
      theme(
        plot.margin = margin(5, 7, 5, 8)
      )

    save_svg(p, file.path(supp_dir, filename), width = 2.55, height = 2.35)
    p
  }

  # -------------------------------------------------
  # Generate sex x phase plots
  # -------------------------------------------------

  analysis_levels <- df_trait %>%
    distinct(Phase) %>%
    filter(!is.na(Phase)) %>%
    pull(Phase) %>%
    as.character()

  plot_objects <- list()

  for (phase_i in analysis_levels) {
    for (sex_i in c("f", "m")) {

      d_i <- df_trait %>%
        filter(Sex == sex_i, Phase == phase_i)

      if (nrow(d_i) < 4) next

      sex_label <- sex_labels[[sex_i]]
      phase_label <- phase_i
      prefix <- paste0(sex_i, "_", phase_i)

      plot_objects[[paste0(prefix, "_dist")]] <- plot_distribution(
        d_i,
        sex_label,
        phase_label,
        paste0("Distribution_", prefix, "_Movement_AUC_Z.svg")
      )

      plot_objects[[paste0(prefix, "_combz")]] <- plot_combz(
        d_i,
        sex_label,
        phase_label,
        paste0("PRIMARY_", prefix, "_Movement_AUC_Z_vs_CombZ.svg")
      )

      plot_objects[[paste0(prefix, "_groupslope")]] <- plot_group_slope(
        d_i,
        sex_label,
        phase_label,
        paste0("GroupSlopeDiagnostic_", prefix, "_Movement_CombZ.svg")
      )

      plot_component(
        d_i,
        sex_label,
        phase_label,
        "delta_cort",
        "Corticosterone Z-score",
        paste0("Component_", prefix, "_Movement_vs_Corticosterone.svg")
      )

      plot_component(
        d_i,
        sex_label,
        phase_label,
        "sucrose_pref",
        "Sucrose preference (%)",
        paste0("Component_", prefix, "_Movement_vs_Sucrose.svg")
      )
    }
  }

  # -------------------------------------------------
  # Forest plot
  # -------------------------------------------------

  p_forest <- cor_sex_phase %>%
    filter(is.finite(estimate)) %>%
    mutate(
      Label = fct_rev(factor(Label)),
      p_label = paste0("p=", format_p(p.value))
    ) %>%
    ggplot(aes(x = estimate, y = Label)) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.28) +
    geom_errorbarh(
      aes(xmin = conf.low, xmax = conf.high),
      height = 0,
      linewidth = 0.34,
      color = "black"
    ) +
    geom_point(size = 2.15, color = "black") +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
    labs(
      title = "Movement–CombZ association",
      subtitle = "Pearson r with 95% CI",
      x = "Correlation coefficient",
      y = NULL
    ) +
    theme_nature() +
    theme(legend.position = "none")

  save_svg(
    p_forest,
    file.path(fig_dir, "PRIMARY_Forest_Movement_CombZ_by_Sex_Phase.svg"),
    width = 2.75,
    height = 2.0
  )

  # -------------------------------------------------
  # Composite figures
  # -------------------------------------------------

  if (all(c(
    "f_Active_dist", "m_Active_dist",
    "f_Active_combz", "m_Active_combz"
  ) %in% names(plot_objects))) {

    fig_active <- (
      plot_objects[["f_Active_dist"]] |
        plot_objects[["m_Active_dist"]] |
        plot_objects[["f_Active_combz"]] |
        plot_objects[["m_Active_combz"]] |
        p_forest
    ) +
      plot_annotation(tag_levels = "A") &
      theme(plot.tag = element_text(face = "bold", size = 9))

    save_svg(
      fig_active,
      file.path(fig_dir, "MAIN_FIGURE_Active_SexSpecific_Movement_CombZ.svg"),
      width = 10.8,
      height = 2.55
    )
  }

  if (all(c(
    "f_Inactive_dist", "m_Inactive_dist",
    "f_Inactive_combz", "m_Inactive_combz"
  ) %in% names(plot_objects))) {

    fig_inactive <- (
      plot_objects[["f_Inactive_dist"]] |
        plot_objects[["m_Inactive_dist"]] |
        plot_objects[["f_Inactive_combz"]] |
        plot_objects[["m_Inactive_combz"]] |
        p_forest
    ) +
      plot_annotation(tag_levels = "A") &
      theme(plot.tag = element_text(face = "bold", size = 9))

    save_svg(
      fig_inactive,
      file.path(fig_dir, "MAIN_FIGURE_Inactive_SexSpecific_Movement_CombZ.svg"),
      width = 10.8,
      height = 2.55
    )
  }

  if (all(c(
    "f_Active_combz", "m_Active_combz",
    "f_Inactive_combz", "m_Inactive_combz"
  ) %in% names(plot_objects))) {

    fig_all_phase_grid <- (
      (plot_objects[["f_Active_combz"]] | plot_objects[["m_Active_combz"]]) /
        (plot_objects[["f_Inactive_combz"]] | plot_objects[["m_Inactive_combz"]])
    ) +
      plot_annotation(tag_levels = "A") &
      theme(plot.tag = element_text(face = "bold", size = 9))

    save_svg(
      fig_all_phase_grid,
      file.path(fig_dir, "MAIN_FIGURE_CombZ_2x2_Sex_by_Phase.svg"),
      width = 5.8,
      height = 5.0
    )
  }

  capture.output(sessionInfo(), file = file.path(qc_dir, "sessionInfo.txt"))

  cat("Done:", analysis_name, "\n")
  cat("Saved to:", saving_dir, "\n")

  invisible(list(
    df_trait = df_trait,
    correlations = cor_sex_phase,
    models = all_models,
    saving_dir = saving_dir
  ))
}

# -------------------------------------------------
# 6) RUN BOTH ANALYSES
# -------------------------------------------------

results <- imap(
  auc_files,
  ~ run_auc_analysis(auc_file = .x, analysis_name = .y)
)

cat("\nAll analyses complete.\n")
cat("Base output folder:\n", base_saving_dir, "\n")