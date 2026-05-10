# ================================================================
# Dynamic Social Network Analysis
# MMMSociability
# ================================================================
# Goal:
#   Convert proximity / interaction data into time-resolved social-network
#   summaries suitable for stress-resilience analyses.
#
# Two supported input formats:
#   A) animal-level table with Proximity per animal x time bin
#      -> estimates social connectedness / fragmentation per animal-period
#   B) optional dyadic table with focal/partner/proximity columns
#      -> computes true graph/network metrics if available
#
# Outputs:
#   analysis_ready/06_behavioral_dynamics/social_networks/
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

source("Functions/behavioral_dynamics_helpers.R")

# ------------------------------------------------
# USER INPUT
# ------------------------------------------------

input_file <- "analysis_ready/03_derived_metrics/phase_based/all_behavior_metrics.csv"
output_dir <- "analysis_ready/06_behavioral_dynamics/social_networks"

# Optional true dyadic proximity file. Leave NULL if unavailable.
dyad_file <- NULL

# Threshold for defining social contact from proximity. If your proximity is
# distance-like, lower values mean closer. If it is closeness-like, higher means closer.
# Set proximity_is_distance accordingly.
proximity_is_distance <- FALSE
contact_quantile <- 0.75

# ------------------------------------------------
# LOAD ANIMAL-LEVEL DATA
# ------------------------------------------------

raw_dat <- read_behavior_table(input_file)
behav <- standardize_behavior_columns(raw_dat)

ensure_dir(output_dir)
ensure_dir(file.path(output_dir, "tables"))
ensure_dir(file.path(output_dir, "figures"))

# ------------------------------------------------
# ANIMAL-LEVEL SOCIAL CONNECTEDNESS / FRAGMENTATION
# ------------------------------------------------

contact_threshold <- if (proximity_is_distance) {
  quantile(behav$Proximity, probs = 1 - contact_quantile, na.rm = TRUE)
} else {
  quantile(behav$Proximity, probs = contact_quantile, na.rm = TRUE)
}

animal_social <- behav %>%
  mutate(
    social_contact_bin = if (proximity_is_distance) Proximity <= contact_threshold else Proximity >= contact_threshold,
    Proximity_z = z_within_metric(Proximity),
    Movement_z = z_within_metric(Movement),
    Entropy_z = z_within_metric(Entropy),
    ActiveIsolation = -Proximity_z + Movement_z,
    PassiveIsolation = -Proximity_z - Movement_z,
    SocialEngagement = Proximity_z + Movement_z
  ) %>%
  group_by(Group, Sex, Phase, CageChange, AnimalNum) %>%
  arrange(TimeIndex, .by_group = TRUE) %>%
  summarise(
    mean_proximity = mean(Proximity, na.rm = TRUE),
    proximity_rmssd = calc_rmssd(Proximity),
    proximity_acf1 = calc_acf1(Proximity),
    contact_fraction = mean(social_contact_bin, na.rm = TRUE),
    contact_switch_rate = mean(social_contact_bin != lag(social_contact_bin), na.rm = TRUE),
    active_isolation = mean(ActiveIsolation, na.rm = TRUE),
    passive_isolation = mean(PassiveIsolation, na.rm = TRUE),
    social_engagement = mean(SocialEngagement, na.rm = TRUE),
    n_bins = n(),
    .groups = "drop"
  )

write_table(animal_social, file.path(output_dir, "tables", "animal_level_social_dynamics.csv"))

p_social <- animal_social %>%
  pivot_longer(
    cols = c(contact_fraction, contact_switch_rate, active_isolation, passive_isolation, social_engagement),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  ggplot(aes(Group, Value, fill = Group)) +
  geom_violin(alpha = 0.5, linewidth = 0.2, trim = FALSE) +
  geom_jitter(width = 0.1, size = 0.8, alpha = 0.7) +
  facet_grid(Metric ~ Phase, scales = "free_y") +
  labs(title = "Animal-level social dynamics", x = NULL, y = NULL) +
  make_nature_theme()

save_plot_svg_pdf(p_social, file.path(output_dir, "figures", "animal_level_social_dynamics"), width = 180, height = 140)

# ------------------------------------------------
# TRUE DYADIC GRAPH ANALYSIS IF DYAD FILE EXISTS
# ------------------------------------------------

if (!is.null(dyad_file) && file.exists(dyad_file)) {

  dyad_raw <- read_behavior_table(dyad_file)

  focal_col <- first_existing_col(dyad_raw, c("Focal", "AnimalNum", "Animal", "MouseID", "animal_1", "id1"), TRUE, "dyad focal column")
  partner_col <- first_existing_col(dyad_raw, c("Partner", "OtherAnimal", "PartnerID", "animal_2", "id2"), TRUE, "dyad partner column")
  time_col <- first_existing_col(dyad_raw, c("HalfHourElapsed", "HalfHourWithinCC0", "HalfHour", "Time", "TimeBin", "datetime", "DateTime"), TRUE, "dyad time column")
  prox_col <- first_existing_col(dyad_raw, c("Proximity", "proximity", "Distance", "distance", "Interaction", "interaction_weight"), TRUE, "dyad proximity column")
  group_col <- first_existing_col(dyad_raw, c("Group", "Phenotype", "Condition", "Treatment", "StressGroup"), FALSE, "dyad group column")
  sex_col <- first_existing_col(dyad_raw, c("Sex", "sex"), FALSE, "dyad sex column")
  phase_col <- first_existing_col(dyad_raw, c("Phase", "phase", "LightDark", "DayNight"), FALSE, "dyad phase column")
  cage_col <- first_existing_col(dyad_raw, c("CageChange", "CC", "CageChangeNum", "Regrouping", "Batch", "Cage"), FALSE, "dyad cage-change column")

  dyad <- dyad_raw %>%
    transmute(
      Focal = .data[[focal_col]],
      Partner = .data[[partner_col]],
      TimeIndex = .data[[time_col]],
      Weight = suppressWarnings(as.numeric(.data[[prox_col]])),
      Group = if (!is.na(group_col)) as.factor(.data[[group_col]]) else factor("All"),
      Sex = if (!is.na(sex_col)) as.factor(.data[[sex_col]]) else factor("All"),
      Phase = if (!is.na(phase_col)) as.factor(.data[[phase_col]]) else factor("All"),
      CageChange = if (!is.na(cage_col)) as.factor(.data[[cage_col]]) else factor("All")
    ) %>%
    filter(!is.na(Focal), !is.na(Partner), Focal != Partner, is.finite(Weight))

  threshold <- if (proximity_is_distance) {
    quantile(dyad$Weight, probs = 1 - contact_quantile, na.rm = TRUE)
  } else {
    quantile(dyad$Weight, probs = contact_quantile, na.rm = TRUE)
  }

  edge_tbl <- dyad %>%
    mutate(Contact = if (proximity_is_distance) Weight <= threshold else Weight >= threshold)

  if (requireNamespace("igraph", quietly = TRUE)) {

    graph_metrics <- edge_tbl %>%
      group_by(Group, Sex, Phase, CageChange, TimeIndex) %>%
      group_modify(~{
        edges <- .x %>% filter(Contact) %>% select(Focal, Partner, Weight)
        nodes <- unique(c(.x$Focal, .x$Partner))
        if (length(nodes) < 2) {
          return(tibble(n_nodes = length(nodes), n_edges = 0, density = NA_real_, modularity = NA_real_, mean_degree = NA_real_))
        }
        g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = tibble(name = nodes))
        tibble(
          n_nodes = igraph::vcount(g),
          n_edges = igraph::ecount(g),
          density = igraph::edge_density(g, loops = FALSE),
          modularity = ifelse(igraph::ecount(g) > 0, igraph::modularity(igraph::cluster_louvain(g)), NA_real_),
          mean_degree = mean(igraph::degree(g), na.rm = TRUE)
        )
      }) %>%
      ungroup()

    write_table(graph_metrics, file.path(output_dir, "tables", "dyadic_graph_metrics_by_time.csv"))

    node_metrics <- edge_tbl %>%
      group_by(Group, Sex, Phase, CageChange, TimeIndex) %>%
      group_modify(~{
        edges <- .x %>% filter(Contact) %>% select(Focal, Partner, Weight)
        nodes <- unique(c(.x$Focal, .x$Partner))
        if (length(nodes) < 2) return(tibble(AnimalNum = nodes, degree = NA_real_, strength = NA_real_, betweenness = NA_real_))
        g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = tibble(name = nodes))
        tibble(
          AnimalNum = igraph::V(g)$name,
          degree = igraph::degree(g),
          strength = igraph::strength(g, weights = igraph::E(g)$Weight),
          betweenness = igraph::betweenness(g, normalized = TRUE)
        )
      }) %>%
      ungroup()

    write_table(node_metrics, file.path(output_dir, "tables", "dyadic_node_metrics_by_time.csv"))
  } else {
    message("Package igraph not installed. Dyadic table detected, but true graph metrics skipped.")
  }
} else {
  message("No dyadic file provided. Wrote animal-level social dynamics only.")
}

message("Dynamic social network analysis complete.")
