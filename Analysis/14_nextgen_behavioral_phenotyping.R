# ================================================================
# Next-Generation Behavioral Phenotyping
# MMMSociability
# ================================================================
# Adds a systems-level phenotyping layer for long-term RFID/home-cage
# stress-resilience data. This script is designed to run after:
#   Analysis/03_build_multiscale_behavior_metrics.R
#
# Main analyses:
#   1. Multiscale complexity and entropy fingerprints
#   2. Early-warning signals of later stress burden
#   3. HMM-style latent behavioral states with robust k-means fallback
#   4. Behavioral energy landscapes and attractor depth
#   5. Dynamic social coupling / network-state summaries
#   6. Integrated next-generation phenotyping index
#
# Output:
#   analysis_ready/14_nextgen_behavioral_phenotyping/<bin_level>/
#     tables/
#     stats_tables/
#     figures/publication_panels/
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(purrr)
  library(tibble)
  library(stringr)
})

# -----------------------------
# User settings
# -----------------------------
project_root <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID"
repo_root <- "C:/Users/topohl/Documents/GitHub/MMMSociability"
bin_level <- "5min_based"

input_file <- file.path(project_root, "analysis_ready/03_derived_metrics", bin_level, "all_behavior_metrics.csv")
output_dir <- file.path(project_root, "analysis_ready/14_nextgen_behavioral_phenotyping", bin_level)

endpoint_file <- NULL
endpoint_cols <- c("CombZ", "stress_z_score", "SucrosePreference", "Corticosterone", "DeltaCorticosterone")
primary_outcome <- "CombZ"

early_phase_pattern <- "active|dark|night"
early_window_hours <- 12
rolling_window_hours <- 2
n_behavioral_states <- 4
state_labels <- c("Inactive", "Explore", "Social", "Burst")
random_seed <- 123

group_levels <- c("CON", "RES", "SUS")
group_colors <- c("CON" = "#3d3b6e", "RES" = "#C6C3BB", "SUS" = "#e63947")
sex_levels <- c("Female", "Male")

# -----------------------------
# Shared helpers
# -----------------------------
source_if_exists <- function(path) if (file.exists(path)) source(path)
source_if_exists(file.path(repo_root, "Functions/behavioral_dynamics_helpers.R"))
source_if_exists(file.path(repo_root, "Functions/behavioral_dynamics_stats_helpers.R"))

if (!exists("ensure_dir")) ensure_dir <- function(path) { if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE); invisible(path) }
if (!exists("write_table")) write_table <- function(x, path) { ensure_dir(dirname(path)); readr::write_csv(x, path); invisible(path) }
if (!exists("save_plot_svg_pdf")) save_plot_svg_pdf <- function(plot, filename_base, width = 85, height = 65, units = "mm") { ensure_dir(dirname(filename_base)); ggsave(paste0(filename_base, ".svg"), plot, width = width, height = height, units = units); ggsave(paste0(filename_base, ".pdf"), plot, width = width, height = height, units = units); invisible(filename_base) }
if (!exists("first_existing_col")) first_existing_col <- function(dat, candidates, required = TRUE, label = "column") { hit <- candidates[candidates %in% names(dat)][1]; if (is.na(hit) && required) stop("Missing ", label, call. = FALSE); hit }
if (!exists("make_nature_theme")) make_nature_theme <- function(base_size = 7) theme_classic(base_size = base_size) + theme(axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.2), strip.background = element_blank(), strip.text = element_text(face = "bold"), legend.title = element_blank(), legend.position = "top", plot.title = element_text(face = "bold", hjust = 0), plot.subtitle = element_text(hjust = 0, colour = "grey35"))

panel_theme <- function(base_size = 7) make_nature_theme(base_size) + theme(panel.grid = element_blank(), legend.position = "right")
safe_num <- function(x) suppressWarnings(as.numeric(x))
safe_scale <- function(x) { s <- sd(x, na.rm = TRUE); m <- mean(x, na.rm = TRUE); if (!is.finite(s) || s == 0) return(rep(0, length(x))); (x - m) / s }
safe_cor <- function(x, y, method = "spearman") { ok <- is.finite(x) & is.finite(y); if (sum(ok) < 4 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_); suppressWarnings(cor(x[ok], y[ok], method = method)) }
calc_acf1 <- function(x) { x <- x[is.finite(x)]; if (length(x) < 4 || sd(x) == 0) return(NA_real_); suppressWarnings(acf(x, plot = FALSE, lag.max = 1)$acf[2]) }
calc_rmssd <- function(x) { x <- x[is.finite(x)]; if (length(x) < 3) return(NA_real_); sqrt(mean(diff(x)^2, na.rm = TRUE)) }
shannon_entropy <- function(p) { p <- p[is.finite(p) & p > 0]; if (length(p) == 0) return(NA_real_); -sum(p * log2(p)) }

read_any <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") return(read_csv(path, show_col_types = FALSE))
  if (ext %in% c("tsv", "txt")) return(read_tsv(path, show_col_types = FALSE))
  if (ext == "rds") return(readRDS(path))
  if (ext %in% c("xlsx", "xls")) { if (!requireNamespace("readxl", quietly = TRUE)) stop("Install readxl"); return(readxl::read_excel(path)) }
  NULL
}

standardize_input <- function(dat) {
  animal_col <- first_existing_col(dat, c("AnimalNum", "Animal", "AnimalID", "MouseID", "Mouse", "ID", "RFID"), TRUE, "animal ID")
  time_col <- first_existing_col(dat, c("TimeIndex", "BinStart", "HalfHourElapsed", "HalfHourWithinCC0", "Time", "DateTime"), TRUE, "time")
  group_col <- first_existing_col(dat, c("Group", "Phenotype", "Condition", "Treatment", "StressGroup"), FALSE, "group")
  sex_col <- first_existing_col(dat, c("Sex", "sex"), FALSE, "sex")
  phase_col <- first_existing_col(dat, c("Phase", "phase", "LightDark", "DayNight"), FALSE, "phase")
  cage_col <- first_existing_col(dat, c("CageChange", "CC", "CageChangeNum", "Regrouping", "Batch", "Cage"), FALSE, "cage")
  batch_col <- first_existing_col(dat, c("Batch", "batch"), FALSE, "batch")
  system_col <- first_existing_col(dat, c("System", "system", "Cage", "CageID"), FALSE, "system")
  movement_col <- first_existing_col(dat, c("Movement", "movement", "MovementDistance", "Distance", "Activity"), TRUE, "movement")
  entropy_col <- first_existing_col(dat, c("Entropy", "entropy", "ShannonEntropy", "PositionEntropy"), TRUE, "entropy")
  proximity_col <- first_existing_col(dat, c("ProximityFraction", "Proximity", "proximity", "MeanProximity", "SocialProximity"), TRUE, "proximity")

  dat %>%
    mutate(
      AnimalNum = .data[[animal_col]],
      TimeIndex = .data[[time_col]],
      Group = if (!is.na(group_col)) as.character(.data[[group_col]]) else "All",
      Sex = if (!is.na(sex_col)) as.character(.data[[sex_col]]) else "All",
      Phase = if (!is.na(phase_col)) as.character(.data[[phase_col]]) else "All",
      CageChange = if (!is.na(cage_col)) as.character(.data[[cage_col]]) else "All",
      Batch = if (!is.na(batch_col)) as.character(.data[[batch_col]]) else "All",
      System = if (!is.na(system_col)) as.character(.data[[system_col]]) else "All",
      Movement = safe_num(.data[[movement_col]]),
      Entropy = safe_num(.data[[entropy_col]]),
      Proximity = safe_num(.data[[proximity_col]])
    ) %>%
    mutate(
      Group = factor(Group, levels = unique(c(group_levels, sort(unique(Group))))),
      Sex = factor(Sex, levels = unique(c(sex_levels, sort(unique(Sex)))))
    ) %>%
    filter(!is.na(AnimalNum), !is.na(TimeIndex)) %>%
    arrange(AnimalNum, CageChange, Phase, TimeIndex)
}

ordinal_time_index <- function(x) {
  if (inherits(x, c("POSIXct", "POSIXt", "Date"))) return(as.numeric(x))
  y <- suppressWarnings(as.numeric(x))
  if (all(is.na(y))) return(as.numeric(factor(x, levels = unique(x))))
  y
}

coarse_grain <- function(x, scale) {
  x <- x[is.finite(x)]
  if (length(x) < scale * 3) return(numeric(0))
  n <- floor(length(x) / scale) * scale
  x <- x[seq_len(n)]
  rowMeans(matrix(x, ncol = scale, byrow = TRUE), na.rm = TRUE)
}

sample_entropy <- function(x, m = 2, r = 0.2) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 20 || sd(x) == 0) return(NA_real_)
  tol <- r * sd(x)
  count_matches <- function(mm) {
    emb <- embed(x, mm)
    cnt <- 0
    for (i in seq_len(nrow(emb) - 1)) {
      d <- apply(abs(t(t(emb[(i + 1):nrow(emb), , drop = FALSE]) - emb[i, ])), 1, max)
      cnt <- cnt + sum(d <= tol, na.rm = TRUE)
    }
    cnt
  }
  a <- count_matches(m + 1)
  b <- count_matches(m)
  if (!is.finite(a) || !is.finite(b) || a == 0 || b == 0) return(NA_real_)
  -log(a / b)
}

permutation_entropy <- function(x, order = 3) {
  x <- x[is.finite(x)]
  if (length(x) < order + 2) return(NA_real_)
  patterns <- map_chr(seq_len(length(x) - order + 1), function(i) paste(order(x[i:(i + order - 1)]), collapse = "-"))
  p <- as.numeric(table(patterns)) / length(patterns)
  shannon_entropy(p) / log2(factorial(order))
}

calc_multiscale_entropy <- function(x, scales = c(1, 2, 3, 6, 12)) {
  map_dfr(scales, function(sc) {
    z <- coarse_grain(x, sc)
    tibble(scale = sc, sample_entropy = sample_entropy(z), permutation_entropy = permutation_entropy(z), acf1 = calc_acf1(z), rmssd = calc_rmssd(z))
  })
}

state_from_centres <- function(centres) {
  label_tbl <- centres %>% mutate(State = NA_character_)
  label_tbl$State[which.min(label_tbl$Movement_z)] <- "Inactive"
  label_tbl$State[which.max(label_tbl$Proximity_z)] <- "Social"
  label_tbl$State[which.max(label_tbl$Movement_z)] <- "Burst"
  label_tbl$State[is.na(label_tbl$State)] <- setdiff(state_labels, label_tbl$State[!is.na(label_tbl$State)])[seq_len(sum(is.na(label_tbl$State)))]
  label_tbl
}

# -----------------------------
# Load and prepare data
# -----------------------------
ensure_dir(output_dir)
ensure_dir(file.path(output_dir, "tables"))
ensure_dir(file.path(output_dir, "stats_tables"))
ensure_dir(file.path(output_dir, "figures/publication_panels"))

raw <- read_any(input_file)
if (is.null(raw)) stop("Missing input file. Run 03_build_multiscale_behavior_metrics.R first.", call. = FALSE)
behav <- standardize_input(raw)

endpoint_dat <- read_any(endpoint_file)
if (!is.null(endpoint_dat)) {
  animal_col <- first_existing_col(endpoint_dat, c("AnimalNum", "Animal", "AnimalID", "MouseID", "Mouse", "ID", "RFID"), TRUE, "endpoint animal ID")
  keep <- intersect(endpoint_cols, names(endpoint_dat))
  endpoint_dat <- endpoint_dat %>% transmute(AnimalNum = .data[[animal_col]], across(all_of(keep), safe_num))
  behav <- left_join(behav, endpoint_dat, by = "AnimalNum")
}

available_outcomes <- intersect(endpoint_cols, names(behav))
available_outcomes <- available_outcomes[map_lgl(available_outcomes, ~ any(is.finite(safe_num(behav[[.x]]))))]
outcome_to_use <- if (primary_outcome %in% available_outcomes) primary_outcome else ifelse(length(available_outcomes) > 0, available_outcomes[1], NA_character_)

bins_per_hour <- case_when(str_detect(bin_level, "10sec") ~ 360, str_detect(bin_level, "5min") ~ 12, str_detect(bin_level, "10min") ~ 6, str_detect(bin_level, "30min") ~ 2, TRUE ~ 12)
max_early_bins <- early_window_hours * bins_per_hour
rolling_n <- max(3, rolling_window_hours * bins_per_hour)

behav <- behav %>% mutate(TimeNumeric = ordinal_time_index(TimeIndex))

early <- behav %>%
  filter(str_detect(str_to_lower(as.character(Phase)), early_phase_pattern)) %>%
  group_by(AnimalNum, CageChange, Phase) %>%
  arrange(TimeNumeric, .by_group = TRUE) %>%
  mutate(EarlyBin = row_number(), EarlyHour = (EarlyBin - 1) / bins_per_hour) %>%
  filter(EarlyBin <= max_early_bins) %>%
  ungroup()
if (nrow(early) == 0) {
  early <- behav %>% group_by(AnimalNum, CageChange, Phase) %>% arrange(TimeNumeric, .by_group = TRUE) %>% mutate(EarlyBin = row_number(), EarlyHour = (EarlyBin - 1) / bins_per_hour) %>% filter(EarlyBin <= max_early_bins) %>% ungroup()
}

# ================================================================
# 1. Multiscale complexity and entropy fingerprints
# ================================================================
complexity_long <- early %>%
  group_by(AnimalNum, Group, Sex) %>%
  arrange(EarlyHour, .by_group = TRUE) %>%
  group_modify(~ {
    bind_rows(
      calc_multiscale_entropy(.x$Movement) %>% mutate(metric = "Movement"),
      calc_multiscale_entropy(.x$Entropy) %>% mutate(metric = "Entropy"),
      calc_multiscale_entropy(.x$Proximity) %>% mutate(metric = "Proximity")
    )
  }) %>%
  ungroup()

complexity_features <- complexity_long %>%
  group_by(AnimalNum, Group, Sex, metric) %>%
  summarise(
    mse_auc = sum(sample_entropy, na.rm = TRUE),
    perm_entropy_mean = mean(permutation_entropy, na.rm = TRUE),
    mean_acf1 = mean(acf1, na.rm = TRUE),
    mean_rmssd = mean(rmssd, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = metric, values_from = c(mse_auc, perm_entropy_mean, mean_acf1, mean_rmssd), names_sep = "_")

write_table(complexity_long, file.path(output_dir, "tables/multiscale_entropy_by_animal.csv"))
write_table(complexity_features, file.path(output_dir, "tables/multiscale_complexity_features.csv"))

p_complexity <- complexity_long %>%
  ggplot(aes(scale, sample_entropy, colour = Group, group = interaction(Group, AnimalNum))) +
  geom_line(alpha = 0.20, linewidth = 0.25) +
  stat_summary(aes(group = Group), fun = mean, geom = "line", linewidth = 0.95) +
  facet_grid(Sex ~ metric, scales = "free_y") +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  labs(title = "A  Multiscale behavioral complexity", subtitle = "Sample entropy across temporal coarse-graining scales in the early active window", x = "Coarse-graining scale", y = "Sample entropy") +
  panel_theme(6)
save_plot_svg_pdf(p_complexity, file.path(output_dir, "figures/publication_panels/Fig14A_multiscale_complexity"), width = 170, height = 95)

# ================================================================
# 2. Early-warning signals: critical slowing, variance and flickering
# ================================================================
rolling_features <- early %>%
  group_by(AnimalNum, Group, Sex) %>%
  arrange(EarlyHour, .by_group = TRUE) %>%
  group_modify(~ {
    x <- .x
    starts <- seq(1, max(1, nrow(x) - rolling_n + 1), by = max(1, floor(rolling_n / 2)))
    map_dfr(starts, function(st) {
      en <- min(nrow(x), st + rolling_n - 1)
      w <- x[st:en, , drop = FALSE]
      tibble(
        WindowStartHour = min(w$EarlyHour, na.rm = TRUE),
        WindowMidHour = mean(range(w$EarlyHour, na.rm = TRUE)),
        Movement_var = var(w$Movement, na.rm = TRUE),
        Movement_acf1 = calc_acf1(w$Movement),
        Movement_rmssd = calc_rmssd(w$Movement),
        Entropy_var = var(w$Entropy, na.rm = TRUE),
        Entropy_acf1 = calc_acf1(w$Entropy),
        Entropy_rmssd = calc_rmssd(w$Entropy),
        Proximity_var = var(w$Proximity, na.rm = TRUE),
        Proximity_acf1 = calc_acf1(w$Proximity),
        Proximity_rmssd = calc_rmssd(w$Proximity)
      )
    })
  }) %>%
  ungroup()

early_warning_summary <- rolling_features %>%
  group_by(AnimalNum, Group, Sex) %>%
  summarise(
    across(ends_with("_var"), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
    across(ends_with("_acf1"), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
    across(ends_with("_rmssd"), ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
    acf1_rise_Movement = safe_cor(WindowMidHour, Movement_acf1, "spearman"),
    acf1_rise_Entropy = safe_cor(WindowMidHour, Entropy_acf1, "spearman"),
    acf1_rise_Proximity = safe_cor(WindowMidHour, Proximity_acf1, "spearman"),
    variance_rise_Movement = safe_cor(WindowMidHour, Movement_var, "spearman"),
    variance_rise_Entropy = safe_cor(WindowMidHour, Entropy_var, "spearman"),
    variance_rise_Proximity = safe_cor(WindowMidHour, Proximity_var, "spearman"),
    .groups = "drop"
  ) %>%
  mutate(
    CriticalSlowingIndex = rowMeans(cbind(safe_scale(mean_Movement_acf1), safe_scale(mean_Entropy_acf1), safe_scale(mean_Proximity_acf1)), na.rm = TRUE),
    FlickeringIndex = rowMeans(cbind(safe_scale(mean_Movement_var), safe_scale(mean_Entropy_var), safe_scale(mean_Proximity_var)), na.rm = TRUE),
    InstabilityRiseIndex = rowMeans(cbind(safe_scale(variance_rise_Movement), safe_scale(variance_rise_Entropy), safe_scale(variance_rise_Proximity)), na.rm = TRUE)
  )

write_table(rolling_features, file.path(output_dir, "tables/rolling_early_warning_features.csv"))
write_table(early_warning_summary, file.path(output_dir, "tables/early_warning_summary_by_animal.csv"))

p_ews <- rolling_features %>%
  select(AnimalNum, Group, Sex, WindowMidHour, Movement_acf1, Entropy_acf1, Proximity_acf1) %>%
  pivot_longer(ends_with("_acf1"), names_to = "metric", values_to = "ACF1") %>%
  mutate(metric = str_remove(metric, "_acf1")) %>%
  ggplot(aes(WindowMidHour, ACF1, colour = Group, group = interaction(Group, AnimalNum))) +
  geom_line(alpha = 0.18, linewidth = 0.25) +
  stat_summary(aes(group = Group), fun = mean, geom = "line", linewidth = 0.95) +
  facet_grid(Sex ~ metric) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  labs(title = "B  Early-warning dynamics", subtitle = "Rolling lag-1 autocorrelation: rising ACF1 is consistent with critical slowing down", x = "Early active hour", y = "Rolling ACF1") +
  panel_theme(6)
save_plot_svg_pdf(p_ews, file.path(output_dir, "figures/publication_panels/Fig14B_early_warning_acf1"), width = 170, height = 95)

# ================================================================
# 3. HMM-style latent states with robust fallback
# ================================================================
state_input <- early %>%
  transmute(AnimalNum, Group, Sex, CageChange, Phase, EarlyBin, EarlyHour,
            Movement_z = safe_scale(Movement), Entropy_z = safe_scale(Entropy), Proximity_z = safe_scale(Proximity)) %>%
  filter(if_all(c(Movement_z, Entropy_z, Proximity_z), is.finite))

set.seed(random_seed)
state_method <- "kmeans_fallback"
if (requireNamespace("depmixS4", quietly = TRUE) && nrow(state_input) >= n_behavioral_states * 30) {
  # depmixS4 fits one pooled Gaussian-emission HMM. The posterior states are then
  # biologically labelled by their emission centroids.
  hmm_dat <- state_input %>% select(Movement_z, Entropy_z, Proximity_z)
  suppressWarnings({
    hmm_mod <- depmixS4::depmix(list(Movement_z ~ 1, Entropy_z ~ 1, Proximity_z ~ 1),
                                data = hmm_dat, nstates = n_behavioral_states,
                                family = list(gaussian(), gaussian(), gaussian()))
    hmm_fit <- try(depmixS4::fit(hmm_mod, verbose = FALSE), silent = TRUE)
  })
  if (!inherits(hmm_fit, "try-error")) {
    post <- depmixS4::posterior(hmm_fit)
    state_input$Cluster <- post$state
    state_method <- "pooled_gaussian_hmm_depmixS4"
  }
}
if (!"Cluster" %in% names(state_input)) {
  km <- kmeans(state_input %>% select(Movement_z, Entropy_z, Proximity_z), centers = n_behavioral_states, nstart = 50)
  state_input$Cluster <- km$cluster
}

centres <- state_input %>% group_by(Cluster) %>% summarise(across(c(Movement_z, Entropy_z, Proximity_z), mean, na.rm = TRUE), .groups = "drop")
label_tbl <- state_from_centres(centres)
state_tbl <- state_input %>% left_join(label_tbl %>% select(Cluster, State), by = "Cluster") %>% mutate(State = factor(State, levels = state_labels), StateMethod = state_method)

state_occupancy <- state_tbl %>% count(AnimalNum, Group, Sex, State, name = "n") %>% group_by(AnimalNum, Group, Sex) %>% mutate(Occupancy = n / sum(n)) %>% ungroup()
state_switching <- state_tbl %>% arrange(AnimalNum, CageChange, Phase, EarlyBin) %>% group_by(AnimalNum, Group, Sex, CageChange, Phase) %>% mutate(NextState = lead(State), Switch = as.numeric(State != NextState)) %>% ungroup() %>% group_by(AnimalNum, Group, Sex) %>% summarise(SwitchingRate = mean(Switch, na.rm = TRUE), StateEntropy = shannon_entropy(as.numeric(table(State)) / sum(table(State))), .groups = "drop")
state_transitions <- state_tbl %>% arrange(AnimalNum, CageChange, Phase, EarlyBin) %>% group_by(AnimalNum, Group, Sex, CageChange, Phase) %>% mutate(From = State, To = lead(State)) %>% ungroup() %>% filter(!is.na(To)) %>% count(Group, Sex, From, To, name = "n") %>% group_by(Group, Sex, From) %>% mutate(Probability = n / sum(n)) %>% ungroup()

write_table(label_tbl %>% mutate(method = state_method), file.path(output_dir, "tables/hmm_style_state_labels.csv"))
write_table(state_tbl, file.path(output_dir, "tables/hmm_style_latent_states_by_bin.csv"))
write_table(state_occupancy, file.path(output_dir, "tables/hmm_style_state_occupancy_by_animal.csv"))
write_table(state_switching, file.path(output_dir, "tables/hmm_style_state_switching_by_animal.csv"))
write_table(state_transitions, file.path(output_dir, "tables/hmm_style_transition_probabilities.csv"))

p_switch <- state_switching %>%
  ggplot(aes(Group, SwitchingRate, fill = Group, colour = Group)) +
  geom_violin(alpha = 0.32, colour = NA, trim = FALSE) +
  geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.25, alpha = 0.75) +
  geom_jitter(width = 0.07, size = 0.9, alpha = 0.75, show.legend = FALSE) +
  facet_grid(. ~ Sex) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  labs(title = "C  Latent-state flexibility", subtitle = paste0("State model: ", state_method, "; higher switching suggests greater behavioral flexibility"), x = NULL, y = "State switching rate") +
  panel_theme(7) + theme(legend.position = "none")
save_plot_svg_pdf(p_switch, file.path(output_dir, "figures/publication_panels/Fig14C_latent_state_flexibility"), width = 105, height = 75)

p_hmm_trans <- state_transitions %>%
  mutate(From = factor(From, levels = state_labels), To = factor(To, levels = state_labels)) %>%
  ggplot(aes(To, From, fill = Probability)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  facet_grid(Sex ~ Group) +
  scale_fill_gradient(low = "white", high = "#333333") +
  labs(title = "D  Latent-state transition architecture", subtitle = "Transition matrix from HMM-style or fallback unsupervised states", x = "To state", y = "From state", fill = "P") +
  panel_theme(6) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot_svg_pdf(p_hmm_trans, file.path(output_dir, "figures/publication_panels/Fig14D_latent_state_transition_heatmap"), width = 165, height = 95)

# ================================================================
# 4. Behavioral energy landscapes and attractor depth
# ================================================================
latent_mat <- state_tbl %>% select(Movement_z, Entropy_z, Proximity_z) %>% as.matrix()
pca <- prcomp(latent_mat, center = TRUE, scale. = TRUE)
landscape_points <- state_tbl %>% bind_cols(as_tibble(pca$x[, 1:2]) %>% setNames(c("Latent1", "Latent2")))

dens_to_energy <- function(dat, bins = 45) {
  x <- dat$Latent1; y <- dat$Latent2
  xr <- range(x, na.rm = TRUE); yr <- range(y, na.rm = TRUE)
  xb <- seq(xr[1], xr[2], length.out = bins + 1)
  yb <- seq(yr[1], yr[2], length.out = bins + 1)
  xi <- pmax(1, pmin(bins, findInterval(x, xb, all.inside = TRUE)))
  yi <- pmax(1, pmin(bins, findInterval(y, yb, all.inside = TRUE)))
  grid <- tibble(xi = xi, yi = yi) %>% count(xi, yi, name = "n") %>% complete(xi = 1:bins, yi = 1:bins, fill = list(n = 0)) %>% mutate(P = (n + 1) / sum(n + 1), Energy = -log(P), Latent1 = (xb[xi] + xb[xi + 1]) / 2, Latent2 = (yb[yi] + yb[yi + 1]) / 2)
  grid
}

energy_maps <- landscape_points %>% group_by(Group, Sex) %>% group_split() %>% map_dfr(function(dat) dens_to_energy(dat) %>% mutate(Group = unique(dat$Group), Sex = unique(dat$Sex)))
energy_summary <- energy_maps %>% group_by(Group, Sex) %>% summarise(MinEnergy = min(Energy, na.rm = TRUE), MedianEnergy = median(Energy, na.rm = TRUE), AttractorDepth = MedianEnergy - MinEnergy, LowEnergyArea = mean(Energy <= quantile(Energy, 0.15, na.rm = TRUE), na.rm = TRUE), .groups = "drop")

write_table(landscape_points, file.path(output_dir, "tables/behavioral_landscape_latent_points.csv"))
write_table(energy_maps, file.path(output_dir, "tables/behavioral_energy_landscape_grid.csv"))
write_table(energy_summary, file.path(output_dir, "tables/behavioral_energy_landscape_summary.csv"))

p_energy <- energy_maps %>%
  ggplot(aes(Latent1, Latent2, fill = Energy)) +
  geom_raster() +
  facet_grid(Sex ~ Group) +
  scale_fill_viridis_c(option = "magma", direction = -1) +
  labs(title = "E  Behavioral energy landscapes", subtitle = "Energy = -log occupancy in latent behavior space; deeper wells indicate attractor-like trapping", x = "Latent behavior 1", y = "Latent behavior 2", fill = "Energy") +
  panel_theme(6)
save_plot_svg_pdf(p_energy, file.path(output_dir, "figures/publication_panels/Fig14E_behavioral_energy_landscape"), width = 165, height = 95)

p_depth <- energy_summary %>%
  ggplot(aes(Group, AttractorDepth, fill = Group, colour = Group)) +
  geom_col(alpha = 0.75, linewidth = 0.25) +
  facet_grid(. ~ Sex) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  labs(title = "F  Attractor depth", subtitle = "Higher values indicate more concentrated low-energy occupancy", x = NULL, y = "Median energy - minimum energy") +
  panel_theme(7) + theme(legend.position = "none")
save_plot_svg_pdf(p_depth, file.path(output_dir, "figures/publication_panels/Fig14F_attractor_depth"), width = 105, height = 75)

# ================================================================
# 5. Dynamic social coupling and cage-level network states
# ================================================================
# Uses synchronized animal-level time bins. Pairwise coupling is computed within
# Batch/System/CageChange/Phase blocks. This is not a dyadic contact graph unless
# dyadic raw data are supplied, but it captures shared temporal coupling of
# movement, entropy and proximity dynamics inside each social unit.

pairwise_coupling <- behav %>%
  group_by(Batch, System, CageChange, Phase, TimeIndex) %>%
  filter(n_distinct(AnimalNum) >= 2) %>%
  ungroup() %>%
  group_by(Batch, System, CageChange, Phase) %>%
  group_modify(~ {
    x <- .x
    animals <- sort(unique(x$AnimalNum))
    if (length(animals) < 2) return(tibble())
    pairs <- t(combn(animals, 2))
    map_dfr(seq_len(nrow(pairs)), function(i) {
      a <- pairs[i, 1]; b <- pairs[i, 2]
      xa <- x %>% filter(AnimalNum == a) %>% arrange(TimeIndex)
      xb <- x %>% filter(AnimalNum == b) %>% arrange(TimeIndex)
      joined <- inner_join(xa %>% select(TimeIndex, Group_a = Group, Sex_a = Sex, Movement_a = Movement, Entropy_a = Entropy, Proximity_a = Proximity), xb %>% select(TimeIndex, Group_b = Group, Sex_b = Sex, Movement_b = Movement, Entropy_b = Entropy, Proximity_b = Proximity), by = "TimeIndex")
      tibble(
        AnimalA = a,
        AnimalB = b,
        GroupA = first(joined$Group_a),
        GroupB = first(joined$Group_b),
        SexA = first(joined$Sex_a),
        SexB = first(joined$Sex_b),
        n_bins = nrow(joined),
        MovementSynchrony = safe_cor(joined$Movement_a, joined$Movement_b, "spearman"),
        EntropySynchrony = safe_cor(joined$Entropy_a, joined$Entropy_b, "spearman"),
        ProximitySynchrony = safe_cor(joined$Proximity_a, joined$Proximity_b, "spearman"),
        CouplingStrength = rowMeans(cbind(MovementSynchrony, EntropySynchrony, ProximitySynchrony), na.rm = TRUE)
      )
    })
  }) %>%
  ungroup() %>%
  filter(is.finite(CouplingStrength))

animal_coupling <- pairwise_coupling %>%
  pivot_longer(c(AnimalA, AnimalB), names_to = "PairRole", values_to = "AnimalNum") %>%
  mutate(Group = if_else(PairRole == "AnimalA", as.character(GroupA), as.character(GroupB)), Sex = if_else(PairRole == "AnimalA", as.character(SexA), as.character(SexB))) %>%
  group_by(AnimalNum, Group, Sex) %>%
  summarise(
    MeanCouplingStrength = mean(CouplingStrength, na.rm = TRUE),
    MeanMovementSynchrony = mean(MovementSynchrony, na.rm = TRUE),
    MeanEntropySynchrony = mean(EntropySynchrony, na.rm = TRUE),
    MeanProximitySynchrony = mean(ProximitySynchrony, na.rm = TRUE),
    CouplingVolatility = sd(CouplingStrength, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Group = factor(Group, levels = unique(c(group_levels, sort(unique(Group))))), Sex = factor(Sex, levels = unique(c(sex_levels, sort(unique(Sex))))))

write_table(pairwise_coupling, file.path(output_dir, "tables/pairwise_dynamic_coupling.csv"))
write_table(animal_coupling, file.path(output_dir, "tables/animal_social_coupling_summary.csv"))

p_coupling <- animal_coupling %>%
  ggplot(aes(Group, MeanCouplingStrength, fill = Group, colour = Group)) +
  geom_violin(alpha = 0.32, colour = NA, trim = FALSE) +
  geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.25, alpha = 0.75) +
  geom_jitter(width = 0.07, size = 0.9, alpha = 0.75, show.legend = FALSE) +
  facet_grid(. ~ Sex) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  labs(title = "G  Dynamic social coupling", subtitle = "Pairwise synchrony of movement, entropy and proximity within social housing blocks", x = NULL, y = "Mean coupling strength") +
  panel_theme(7) + theme(legend.position = "none")
save_plot_svg_pdf(p_coupling, file.path(output_dir, "figures/publication_panels/Fig14G_dynamic_social_coupling"), width = 105, height = 75)

# ================================================================
# 6. Integrated next-generation phenotype index and endpoint links
# ================================================================
feature_blocks <- list(complexity_features, early_warning_summary, state_switching, animal_coupling)
nextgen_features <- reduce(feature_blocks, full_join, by = c("AnimalNum", "Group", "Sex")) %>%
  left_join(energy_summary %>% group_by(Group, Sex) %>% summarise(GroupSexAttractorDepth = first(AttractorDepth), .groups = "drop"), by = c("Group", "Sex"))

num_features <- nextgen_features %>% select(where(is.numeric)) %>% names() %>% setdiff("AnimalNum")
nextgen_features <- nextgen_features %>%
  mutate(
    ComplexityIndex = rowMeans(across(any_of(names(.)[str_detect(names(.), "mse_auc|perm_entropy")])), na.rm = TRUE),
    EarlyWarningIndex = rowMeans(across(any_of(c("CriticalSlowingIndex", "FlickeringIndex", "InstabilityRiseIndex"))), na.rm = TRUE),
    FlexibilityIndex = rowMeans(cbind(safe_scale(SwitchingRate), safe_scale(StateEntropy)), na.rm = TRUE),
    SocialCouplingIndex = safe_scale(MeanCouplingStrength),
    NextGenPhenotypeIndex = rowMeans(cbind(
      safe_scale(ComplexityIndex),
      safe_scale(EarlyWarningIndex),
      safe_scale(FlexibilityIndex),
      safe_scale(SocialCouplingIndex),
      safe_scale(GroupSexAttractorDepth)
    ), na.rm = TRUE)
  )

if (!is.na(outcome_to_use)) {
  outcome_tbl <- behav %>% group_by(AnimalNum) %>% summarise(outcome = first(na.omit(safe_num(.data[[outcome_to_use]]))), .groups = "drop")
  nextgen_features <- left_join(nextgen_features, outcome_tbl, by = "AnimalNum")
}

nextgen_correlations <- if ("outcome" %in% names(nextgen_features)) {
  nextgen_features %>%
    select(AnimalNum, Group, Sex, outcome, where(is.numeric)) %>%
    select(-AnimalNum) %>%
    pivot_longer(-c(Group, Sex, outcome), names_to = "feature", values_to = "value") %>%
    group_by(feature) %>%
    summarise(n = sum(is.finite(value) & is.finite(outcome)), spearman_rho = safe_cor(value, outcome, "spearman"), pearson_r = safe_cor(value, outcome, "pearson"), .groups = "drop") %>%
    arrange(desc(abs(spearman_rho)))
} else tibble(feature = character(), n = integer(), spearman_rho = numeric(), pearson_r = numeric())

write_table(nextgen_features, file.path(output_dir, "tables/nextgen_behavioral_phenotype_matrix.csv"))
write_table(nextgen_correlations, file.path(output_dir, "stats_tables/nextgen_feature_endpoint_correlations.csv"))

p_index <- nextgen_features %>%
  ggplot(aes(Group, NextGenPhenotypeIndex, fill = Group, colour = Group)) +
  geom_hline(yintercept = 0, linewidth = 0.25, linetype = "dashed", colour = "grey55") +
  geom_violin(alpha = 0.32, colour = NA, trim = FALSE) +
  geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.25, alpha = 0.75) +
  geom_jitter(width = 0.07, size = 0.9, alpha = 0.75, show.legend = FALSE) +
  facet_grid(. ~ Sex) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  scale_colour_manual(values = group_colors, drop = FALSE) +
  labs(title = "H  Integrated next-generation behavioral phenotype", subtitle = "Composite of complexity, early-warning, flexibility, coupling and attractor-depth features", x = NULL, y = "Next-generation phenotype index") +
  panel_theme(7) + theme(legend.position = "none")
save_plot_svg_pdf(p_index, file.path(output_dir, "figures/publication_panels/Fig14H_nextgen_phenotype_index"), width = 105, height = 75)

if ("outcome" %in% names(nextgen_features)) {
  p_endpoint <- nextgen_features %>%
    filter(is.finite(outcome), is.finite(NextGenPhenotypeIndex)) %>%
    ggplot(aes(NextGenPhenotypeIndex, outcome, colour = Group, fill = Group, shape = Sex)) +
    geom_point(size = 2.1, alpha = 0.92, stroke = 0.25) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.55, colour = "grey20") +
    scale_colour_manual(values = group_colors, drop = FALSE) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    labs(title = "I  Systems phenotype-to-endpoint coupling", subtitle = paste0("Outcome: ", outcome_to_use, "; association is predictive/screening, not causal inference"), x = "Next-generation phenotype index", y = outcome_to_use) +
    panel_theme(7)
  save_plot_svg_pdf(p_endpoint, file.path(output_dir, "figures/publication_panels/Fig14I_nextgen_index_endpoint_coupling"), width = 95, height = 78)
}

# -----------------------------
# Analysis manifest
# -----------------------------
manifest <- tibble(
  analysis = c("multiscale_complexity", "early_warning_signals", "hmm_style_latent_states", "energy_landscapes", "dynamic_social_coupling", "nextgen_index"),
  output_prefix = c("Fig14A", "Fig14B", "Fig14C-D", "Fig14E-F", "Fig14G", "Fig14H-I"),
  interpretation = c(
    "Temporal complexity across scales; reduced complexity may indicate stereotypy, rigidity or collapse of adaptive behavioral repertoire.",
    "Rolling variance and autocorrelation; rising ACF1 is consistent with critical slowing but requires validation against endpoints.",
    "Latent state flexibility and transition architecture; HMM is used when depmixS4 is installed, otherwise k-means fallback is used.",
    "Occupancy-derived energy landscapes; deep attractor wells suggest trapping in restricted behavioral regions.",
    "Within-cage synchrony of animal-level dynamics; proxy for social coupling, not raw dyadic contact causality.",
    "Composite systems phenotype; intended for hypothesis generation and endpoint association, not as a finalized biomarker without cross-validation."
  )
)
write_table(manifest, file.path(output_dir, "analysis_manifest.csv"))

message("Next-generation behavioral phenotyping complete: ", output_dir)
