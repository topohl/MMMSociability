# ================================================================
# Behavior ↔ Proteomics Integration
# MMMSociability
# ================================================================
# Goal:
#   Integrate behavioral dynamics with hippocampal proteomic modules.
#
# Compatible with:
#   - module scores
#   - region/layer proteomics
#   - neuropil enrichment signatures
#
# Outputs:
#   - correlation matrices
#   - behavior/proteomics latent relationships
#   - clustered heatmaps
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(purrr)
})

source("Functions/behavioral_dynamics_helpers.R")

behavior_file <- "analysis_ready/06_behavioral_dynamics/early_prediction/tables/early_behavior_features.csv"
proteomics_file <- "analysis_ready/proteomics/module_scores.csv"
output_dir <- "analysis_ready/06_behavioral_dynamics/proteomics_integration"

ensure_dir(output_dir)
ensure_dir(file.path(output_dir, "tables"))
ensure_dir(file.path(output_dir, "figures"))

if (!file.exists(behavior_file)) {
  stop("Behavior feature file not found: ", behavior_file)
}

if (!file.exists(proteomics_file)) {
  stop("Proteomics module file not found: ", proteomics_file)
}

beh <- read_behavior_table(behavior_file)
prot <- read_behavior_table(proteomics_file)

beh_animal_col <- first_existing_col(beh, c("AnimalNum", "Animal", "MouseID", "Mouse", "ID"), TRUE, "behavior animal column")
prot_animal_col <- first_existing_col(prot, c("AnimalNum", "Animal", "MouseID", "Mouse", "ID"), TRUE, "proteomics animal column")

beh <- beh %>% rename(AnimalNum = all_of(beh_animal_col))
prot <- prot %>% rename(AnimalNum = all_of(prot_animal_col))

merged <- beh %>% inner_join(prot, by = "AnimalNum")

write_table(merged, file.path(output_dir, "tables", "behavior_proteomics_merged.csv"))

beh_num <- merged %>%
  select(where(is.numeric))

cor_mat <- suppressWarnings(cor(beh_num, use = "pairwise.complete.obs", method = "spearman"))

cor_tbl <- as.data.frame(as.table(cor_mat)) %>%
  rename(Feature1 = Var1, Feature2 = Var2, SpearmanRho = Freq) %>%
  filter(Feature1 != Feature2) %>%
  arrange(desc(abs(SpearmanRho)))

write_table(cor_tbl, file.path(output_dir, "tables", "behavior_proteomics_correlations.csv"))

# crude separation of behavior vs proteomics variables
beh_vars <- names(beh %>% select(where(is.numeric)))
prot_vars <- setdiff(names(prot %>% select(where(is.numeric))), beh_vars)

cross_tbl <- cor_tbl %>%
  filter(
    (Feature1 %in% beh_vars & Feature2 %in% prot_vars) |
    (Feature2 %in% beh_vars & Feature1 %in% prot_vars)
  )

write_table(cross_tbl, file.path(output_dir, "tables", "crossmodal_behavior_proteomics_correlations.csv"))

plot_tbl <- cross_tbl %>%
  slice_max(abs(SpearmanRho), n = 80)

p_heat <- plot_tbl %>%
  ggplot(aes(Feature1, Feature2, fill = SpearmanRho)) +
  geom_tile() +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0) +
  labs(
    title = "Behavior ↔ proteomics associations",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 7) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid = element_blank()
  )

save_plot_svg_pdf(p_heat, file.path(output_dir, "figures", "behavior_proteomics_heatmap"), width = 180, height = 140)

if (requireNamespace("mixOmics", quietly = TRUE)) {

  x <- merged %>%
    select(any_of(beh_vars)) %>%
    select(where(is.numeric)) %>%
    as.matrix()

  y <- merged %>%
    select(any_of(prot_vars)) %>%
    select(where(is.numeric)) %>%
    as.matrix()

  if (ncol(x) > 1 && ncol(y) > 1) {

    pls_fit <- mixOmics::pls(x, y, ncomp = 2)

    latent_tbl <- tibble(
      AnimalNum = merged$AnimalNum,
      comp1_behavior = pls_fit$variates$X[, 1],
      comp1_proteomics = pls_fit$variates$Y[, 1]
    )

    write_table(latent_tbl, file.path(output_dir, "tables", "pls_latent_scores.csv"))

    p_pls <- latent_tbl %>%
      ggplot(aes(comp1_behavior, comp1_proteomics)) +
      geom_point(size = 1.5, alpha = 0.8) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 0.5) +
      labs(
        title = "Latent behavior ↔ proteomics relationship",
        x = "Behavior latent component",
        y = "Proteomics latent component"
      ) +
      make_nature_theme()

    save_plot_svg_pdf(p_pls, file.path(output_dir, "figures", "latent_behavior_proteomics_relationship"), width = 85, height = 75)
  }
}

message("Behavior-proteomics integration complete.")
