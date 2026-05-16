# ================================================================
# Behavioral dynamics helper functions
# MMMSociability
# ================================================================
# Purpose:
#   Shared utilities for temporal instability, behavioral state-space,
#   and early prediction analyses based on movement / entropy / proximity
#   time-series data.
#
# Design:
#   - Self-contained and tolerant to common column-name variants.
#   - Does not modify existing pipeline objects.
#   - Writes analysis-ready tables and Nature-style ggplot figures.
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

mmm_group_levels <- c("CON", "RES", "SUS")
mmm_group_colors <- c("CON" = "#3d3b6e", "RES" = "#C6C3BB", "SUS" = "#e63947", "All" = "grey55")
mmm_pair_colors <- c("RES-CON" = "#3d3b6e", "SUS-CON" = "#e63947", "SUS-RES" = "#8A817C")
mmm_state_colors <- c("#2F4858", "#4D908E", "#7E9F35", "#F2A65A", "#B23A48", "#6D597A")
mmm_diverging_colors <- c(low = "#3d3b6e", mid = "white", high = "#e63947")

safe_name <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "") %>%
    str_to_lower()
}

pretty_metric_label <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("__", " | ") %>%
    str_replace_all("_", " ") %>%
    str_squish() %>%
    str_to_sentence()
}

analysis_output_dirs <- function(output_dir) {
  dirs <- list(
    root = output_dir,
    tables = file.path(output_dir, "tables"),
    stats = file.path(output_dir, "stats_tables"),
    qc = file.path(output_dir, "qc"),
    table_qc = file.path(output_dir, "tables", "qc"),
    table_sensitivity = file.path(output_dir, "tables", "duration_sensitivity"),
    figure_root = file.path(output_dir, "figures"),
    figure_publication = file.path(output_dir, "figures", "publication_panels"),
    figure_supplementary = file.path(output_dir, "figures", "supplementary"),
    figure_qc = file.path(output_dir, "figures", "qc"),
    figure_exploratory = file.path(output_dir, "figures", "exploratory"),
    figure_interactive = file.path(output_dir, "figures", "interactive")
  )
  purrr::walk(unlist(dirs), ensure_dir)
  dirs
}

write_output_manifest <- function(output_dir,
                                  script_name,
                                  analysis_name,
                                  primary_tables = character(),
                                  primary_figures = character(),
                                  notes = character()) {
  analysis_output_dirs(output_dir)
  manifest <- tibble(
    script = script_name,
    analysis = analysis_name,
    output_type = c(rep("table", length(primary_tables)), rep("figure", length(primary_figures)), rep("note", length(notes))),
    output = c(primary_tables, primary_figures, notes),
    recommended_use = c(
      rep("analysis_ready_table", length(primary_tables)),
      rep("publication_or_supplementary_figure", length(primary_figures)),
      rep("interpretation_note", length(notes))
    )
  )
  readr::write_csv(manifest, file.path(output_dir, "output_manifest.csv"))
  invisible(manifest)
}

first_existing_col <- function(dat, candidates, required = TRUE, label = "column") {
  hit <- candidates[candidates %in% names(dat)][1]
  if (is.na(hit) && required) {
    stop(
      "Could not find ", label, ". Tried: ", paste(candidates, collapse = ", "),
      "\nAvailable columns: ", paste(names(dat), collapse = ", "),
      call. = FALSE
    )
  }
  hit %||% NA_character_
}

standardize_behavior_columns <- function(dat,
                                         animal_col = NULL,
                                         time_col = NULL,
                                         group_col = NULL,
                                         sex_col = NULL,
                                         phase_col = NULL,
                                         cage_col = NULL,
                                         movement_col = NULL,
                                         entropy_col = NULL,
                                         proximity_col = NULL) {
  stopifnot(is.data.frame(dat))

  animal_col <- animal_col %||% first_existing_col(dat, c("AnimalNum", "Animal", "MouseID", "Mouse", "ID", "RFID", "animal_id"), TRUE, "animal ID column")
  time_col <- time_col %||% first_existing_col(dat, c("TimeIndex", "BinStart", "HalfHourElapsed", "HalfHourWithinCC0", "HalfHour", "Time", "TimeBin", "ZeitgeberTime", "ZT", "datetime", "DateTime"), TRUE, "time column")
  group_col <- group_col %||% first_existing_col(dat, c("Group", "Phenotype", "Condition", "Treatment", "StressGroup"), FALSE, "group column")
  sex_col <- sex_col %||% first_existing_col(dat, c("Sex", "sex"), FALSE, "sex column")
  phase_col <- phase_col %||% first_existing_col(dat, c("Phase", "phase", "LightDark", "DayNight", "CircadianPhase"), FALSE, "phase column")
  cage_col <- cage_col %||% first_existing_col(dat, c("CageChange", "CC", "CageChangeNum", "Regrouping", "Batch", "Cage"), FALSE, "cage-change column")
  movement_col <- movement_col %||% first_existing_col(dat, c("Movement", "movement", "Distance", "distance", "Activity", "activity"), TRUE, "movement column")
  entropy_col <- entropy_col %||% first_existing_col(dat, c("Entropy", "entropy", "ShannonEntropy", "shannon_entropy", "PositionEntropy"), TRUE, "entropy column")
  proximity_col <- proximity_col %||% first_existing_col(dat, c("Proximity", "proximity", "MeanProximity", "SocialProximity", "CloseProximity"), TRUE, "proximity column")

  out <- dat %>%
    mutate(
      AnimalNum = .data[[animal_col]],
      TimeIndex = .data[[time_col]],
      Movement = suppressWarnings(as.numeric(.data[[movement_col]])),
      Entropy = suppressWarnings(as.numeric(.data[[entropy_col]])),
      Proximity = suppressWarnings(as.numeric(.data[[proximity_col]]))
    )

  out$Group <- if (!is.na(group_col)) as.factor(dat[[group_col]]) else factor("All")
  out$Sex <- if (!is.na(sex_col)) as.factor(dat[[sex_col]]) else factor("All")
  out$Phase <- if (!is.na(phase_col)) as.factor(dat[[phase_col]]) else factor("All")
  out$CageChange <- if (!is.na(cage_col)) as.factor(dat[[cage_col]]) else factor("All")

  out %>%
    filter(!is.na(AnimalNum), !is.na(TimeIndex)) %>%
    arrange(AnimalNum, CageChange, Phase, TimeIndex)
}

read_behavior_table <- function(input_file) {
  if (is.null(input_file) || !file.exists(input_file)) {
    stop("input_file does not exist: ", input_file, call. = FALSE)
  }
  ext <- tools::file_ext(input_file) %>% tolower()
  if (ext %in% c("csv")) return(readr::read_csv(input_file, show_col_types = FALSE))
  if (ext %in% c("tsv", "txt")) return(readr::read_tsv(input_file, show_col_types = FALSE))
  if (ext %in% c("rds")) return(readRDS(input_file))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) stop("Install readxl to read Excel files.")
    return(readxl::read_excel(input_file))
  }
  stop("Unsupported file extension: ", ext, call. = FALSE)
}

z_within_metric <- function(x) {
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

calc_rmssd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  sqrt(mean(diff(x)^2, na.rm = TRUE))
}

calc_acf1 <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4 || sd(x, na.rm = TRUE) == 0) return(NA_real_)
  suppressWarnings(stats::acf(x, plot = FALSE, lag.max = 1, na.action = na.pass)$acf[2])
}

calc_instability_metrics <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0) {
    return(tibble(n_bins = 0, mean = NA_real_, sd = NA_real_, cv = NA_real_,
                  fano = NA_real_, rmssd = NA_real_, acf1 = NA_real_, p95 = NA_real_,
                  max = NA_real_))
  }
  mu <- mean(x, na.rm = TRUE)
  sig <- sd(x, na.rm = TRUE)
  tibble(
    n_bins = n,
    mean = mu,
    sd = sig,
    cv = ifelse(is.finite(mu) && abs(mu) > .Machine$double.eps, sig / abs(mu), NA_real_),
    fano = ifelse(is.finite(mu) && abs(mu) > .Machine$double.eps, stats::var(x, na.rm = TRUE) / abs(mu), NA_real_),
    rmssd = calc_rmssd(x),
    acf1 = calc_acf1(x),
    p95 = suppressWarnings(quantile(x, 0.95, na.rm = TRUE, names = FALSE)),
    max = suppressWarnings(max(x, na.rm = TRUE))
  )
}

calc_burst_events <- function(x, time_index = seq_along(x), threshold_quantile = 0.90) {
  x <- suppressWarnings(as.numeric(x))
  threshold <- suppressWarnings(stats::quantile(x, threshold_quantile, na.rm = TRUE, names = FALSE))
  if (!is.finite(threshold)) {
    return(tibble(
      burst_threshold = NA_real_, burst_count = NA_integer_, burst_fraction = NA_real_,
      mean_burst_duration_bins = NA_real_, median_burst_duration_bins = NA_real_,
      mean_interburst_interval_bins = NA_real_, mean_burst_amplitude = NA_real_,
      max_burst_amplitude = NA_real_
    ))
  }
  is_burst <- is.finite(x) & x >= threshold
  r <- rle(is_burst)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  burst_runs <- tibble(value = r$values, start = starts, end = ends, duration = r$lengths) %>%
    filter(value)
  if (nrow(burst_runs) == 0) {
    return(tibble(
      burst_threshold = threshold, burst_count = 0L, burst_fraction = mean(is_burst, na.rm = TRUE),
      mean_burst_duration_bins = NA_real_, median_burst_duration_bins = NA_real_,
      mean_interburst_interval_bins = NA_real_, mean_burst_amplitude = NA_real_,
      max_burst_amplitude = NA_real_
    ))
  }
  burst_values <- unlist(map2(burst_runs$start, burst_runs$end, ~x[.x:.y]), use.names = FALSE)
  interburst <- if (nrow(burst_runs) > 1) burst_runs$start[-1] - burst_runs$end[-nrow(burst_runs)] else NA_real_
  tibble(
    burst_threshold = threshold,
    burst_count = nrow(burst_runs),
    burst_fraction = mean(is_burst, na.rm = TRUE),
    mean_burst_duration_bins = mean(burst_runs$duration, na.rm = TRUE),
    median_burst_duration_bins = median(burst_runs$duration, na.rm = TRUE),
    mean_interburst_interval_bins = mean(interburst, na.rm = TRUE),
    mean_burst_amplitude = mean(burst_values, na.rm = TRUE),
    max_burst_amplitude = max(burst_values, na.rm = TRUE)
  )
}

assign_movement_state <- function(x,
                                  probs = c(0.25, 0.60, 0.85),
                                  labels = c("inactive", "low", "active", "burst")) {
  x <- suppressWarnings(as.numeric(x))
  q <- suppressWarnings(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
  if (any(!is.finite(q)) || length(unique(q)) < length(q)) {
    return(factor(rep(NA_character_, length(x)), levels = labels))
  }
  cut(x, breaks = c(-Inf, q, Inf), labels = labels, include.lowest = TRUE, ordered_result = TRUE)
}

calc_transition_metrics <- function(state) {
  state <- as.character(state)
  state <- state[!is.na(state)]
  if (length(state) < 3) {
    return(tibble(
      n_transitions = length(state) - 1L,
      transition_entropy = NA_real_, switch_probability = NA_real_,
      mean_dwell_bins = NA_real_, max_dwell_bins = NA_real_
    ))
  }
  from <- state[-length(state)]
  to <- state[-1]
  tab <- table(from, to)
  p <- as.numeric(tab) / sum(tab)
  p <- p[p > 0]
  r <- rle(state)
  tibble(
    n_transitions = length(from),
    transition_entropy = -sum(p * log2(p)),
    switch_probability = mean(from != to, na.rm = TRUE),
    mean_dwell_bins = mean(r$lengths, na.rm = TRUE),
    max_dwell_bins = max(r$lengths, na.rm = TRUE)
  )
}

make_transition_matrix_long <- function(state) {
  state <- as.character(state)
  state <- state[!is.na(state)]
  if (length(state) < 2) return(tibble(From = character(), To = character(), n = integer(), probability = numeric()))
  tibble(From = state[-length(state)], To = state[-1]) %>%
    count(From, To, name = "n") %>%
    group_by(From) %>%
    mutate(probability = n / sum(n)) %>%
    ungroup()
}

calc_circadian_metrics <- function(value, phase) {
  value <- suppressWarnings(as.numeric(value))
  phase <- as.character(phase)
  ok <- is.finite(value) & !is.na(phase)
  value <- value[ok]
  phase <- phase[ok]
  if (length(value) < 6 || length(unique(phase)) < 2) {
    return(tibble(
      intradaily_variability = NA_real_, phase_stability = NA_real_,
      relative_amplitude = NA_real_, active_inactive_ratio = NA_real_
    ))
  }
  mu <- mean(value, na.rm = TRUE)
  iv <- ifelse(stats::var(value, na.rm = TRUE) > 0,
               mean(diff(value)^2, na.rm = TRUE) / stats::var(value, na.rm = TRUE),
               NA_real_)
  phase_means <- tapply(value, phase, mean, na.rm = TRUE)
  phase_ns <- tapply(value, phase, length)
  ps <- ifelse(stats::var(value, na.rm = TRUE) > 0,
               sum(phase_ns * (phase_means - mu)^2, na.rm = TRUE) / ((length(value) - 1) * stats::var(value, na.rm = TRUE)),
               NA_real_)
  high <- max(phase_means, na.rm = TRUE)
  low <- min(phase_means, na.rm = TRUE)
  ra <- ifelse(is.finite(high + low) && abs(high + low) > .Machine$double.eps, (high - low) / (high + low), NA_real_)
  tibble(
    intradaily_variability = iv,
    phase_stability = ps,
    relative_amplitude = ra,
    active_inactive_ratio = ifelse(is.finite(low) && abs(low) > .Machine$double.eps, high / abs(low), NA_real_)
  )
}

calc_dfa_alpha <- function(x,
                           window_sizes = c(4, 6, 8, 12, 16, 24, 32),
                           min_windows = 3) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < max(window_sizes, na.rm = TRUE) * min_windows) {
    return(tibble(dfa_alpha = NA_real_, dfa_n_scales = 0L))
  }
  y <- cumsum(x - mean(x, na.rm = TRUE))
  fluct <- map_dfr(window_sizes, function(w) {
    n_seg <- floor(n / w)
    if (n_seg < min_windows) return(tibble(window = w, fluctuation = NA_real_))
    f <- map_dbl(seq_len(n_seg), function(i) {
      idx <- ((i - 1) * w + 1):(i * w)
      fit <- lm(y[idx] ~ idx)
      sqrt(mean(stats::resid(fit)^2, na.rm = TRUE))
    })
    tibble(window = w, fluctuation = sqrt(mean(f^2, na.rm = TRUE)))
  }) %>%
    filter(is.finite(fluctuation), fluctuation > 0, is.finite(window), window > 1)
  if (nrow(fluct) < 3) return(tibble(dfa_alpha = NA_real_, dfa_n_scales = nrow(fluct)))
  fit <- lm(log10(fluctuation) ~ log10(window), data = fluct)
  tibble(dfa_alpha = unname(coef(fit)[2]), dfa_n_scales = nrow(fluct))
}

make_nature_theme <- function(base_size = 7, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.28, colour = "black"),
      axis.ticks = element_line(linewidth = 0.22, colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = "black"),
      legend.title = element_blank(),
      legend.position = "top",
      legend.key.height = unit(3.4, "mm"),
      legend.key.width = unit(6, "mm"),
      panel.spacing = unit(1.0, "lines"),
      plot.title = element_text(face = "bold", hjust = 0, margin = margin(b = 2)),
      plot.subtitle = element_text(hjust = 0, colour = "grey25", margin = margin(b = 3)),
      plot.caption = element_text(hjust = 0, colour = "grey35", size = rel(0.82)),
      plot.margin = margin(4, 4, 4, 4)
    )
}

make_publication_theme <- function(base_size = 7) {
  make_nature_theme(base_size = base_size) +
    theme(
      panel.grid.major.y = element_line(linewidth = 0.13, colour = "grey92"),
      panel.grid.major.x = element_blank(),
      legend.position = "top"
    )
}

save_plot_svg_pdf <- function(plot, filename_base, width = 85, height = 65, units = "mm", dpi = 600) {
  ensure_dir(dirname(filename_base))
  ggplot2::ggsave(paste0(filename_base, ".svg"), plot, width = width, height = height, units = units)
  pdf_device <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else "pdf"
  ggplot2::ggsave(paste0(filename_base, ".pdf"), plot, width = width, height = height, units = units, device = pdf_device)
  ggplot2::ggsave(paste0(filename_base, ".png"), plot, width = width, height = height, units = units, dpi = dpi)
  invisible(filename_base)
}

write_table <- function(x, path) {
  ensure_dir(dirname(path))
  readr::write_csv(x, path)
  invisible(path)
}

safe_cor <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}
