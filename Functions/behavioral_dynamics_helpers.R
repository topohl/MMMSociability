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

safe_name <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_|_$", "") %>%
    str_to_lower()
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
  time_col <- time_col %||% first_existing_col(dat, c("HalfHourElapsed", "HalfHourWithinCC0", "HalfHour", "Time", "TimeBin", "ZeitgeberTime", "ZT", "datetime", "DateTime"), TRUE, "time column")
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

make_nature_theme <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = element_line(linewidth = 0.25, colour = "black"),
      axis.text = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.title = element_blank(),
      legend.position = "top",
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(hjust = 0)
    )
}

save_plot_svg_pdf <- function(plot, filename_base, width = 85, height = 65, units = "mm") {
  ensure_dir(dirname(filename_base))
  ggplot2::ggsave(paste0(filename_base, ".svg"), plot, width = width, height = height, units = units)
  ggplot2::ggsave(paste0(filename_base, ".pdf"), plot, width = width, height = height, units = units)
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
