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
analysis_output_dirs(output_dir)
write_output_manifest(
  output_dir,
  script_name = "12_behavior_proteomics_integration.R",
  analysis_name = "low-dimensional behavior-proteomics systems alignment",
  primary_tables = c(
    "tables/behavioral_latent_axes.csv",
    "tables/proteomic_latent_axes.csv",
    "tables/latent_axis_associations.csv",
    "tables/curated_behavior_proteomics_models.csv"
  ),
  primary_figures = c("figures/latent_behavior_proteomics_relationship.svg"),
  notes = c("Use curated low-dimensional axes; feature-explosion correlations are retained only as exploratory supplements.")
)

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

make_axis <- function(dat, regex_pattern, axis_name) {
  cols <- names(dat)[str_detect(names(dat), regex(regex_pattern, ignore_case = TRUE))]
  cols <- cols[sapply(dat[cols], is.numeric)]
  if (length(cols) == 0) return(tibble(AnimalNum = dat$AnimalNum, !!axis_name := NA_real_, n_features = 0L))
  tibble(
    AnimalNum = dat$AnimalNum,
    !!axis_name := rowMeans(as.data.frame(lapply(dat[cols], z_within_metric)), na.rm = TRUE),
    n_features = length(cols)
  )
}

behavior_axis_specs <- tibble(
  Axis = c("locomotor_adaptation_axis", "temporal_organization_axis", "phase_organization_axis", "social_organization_axis", "inactivity_quiescence_axis"),
  Regex = c("movement|locomotor|adaptation|recovery|slope|half_life", "rmssd|acf1|entropy|instability|volatility|persistence", "phase|active_minus_inactive|phase_contrast|timing", "social|proximity|partner|network|degree|strength|fragmentation", "inactive|inactivity|quiescence|rest_like|bout")
)

proteomic_axis_specs <- tibble(
  Axis = c("RNA_RNP_splicing_axis", "translation_ribosome_axis", "mitochondrial_OXPHOS_axis", "proteostasis_endolysosomal_axis", "synaptic_plasticity_axis"),
  Regex = c("rna|rnp|splice|splicing", "translation|ribosome|rpl|rps|eif", "mitochond|oxphos|respirat|electron", "proteostasis|lysosom|endosom|ubiquitin|proteasom", "synap|plastic|vesicle|glutamate|gaba")
)

behavioral_latent_axes <- reduce(
  map2(behavior_axis_specs$Regex, behavior_axis_specs$Axis, ~make_axis(beh, .x, .y) %>% select(-n_features)),
  full_join,
  by = "AnimalNum"
)

proteomic_latent_axes <- reduce(
  map2(proteomic_axis_specs$Regex, proteomic_axis_specs$Axis, ~make_axis(prot, .x, .y) %>% select(-n_features)),
  full_join,
  by = "AnimalNum"
)

axis_feature_inventory <- bind_rows(
  map2_dfr(behavior_axis_specs$Regex, behavior_axis_specs$Axis, ~tibble(Modality = "behavior", Axis = .y, Feature = names(beh)[str_detect(names(beh), regex(.x, ignore_case = TRUE))])),
  map2_dfr(proteomic_axis_specs$Regex, proteomic_axis_specs$Axis, ~tibble(Modality = "proteomics", Axis = .y, Feature = names(prot)[str_detect(names(prot), regex(.x, ignore_case = TRUE))]))
)

axis_merged <- behavioral_latent_axes %>%
  inner_join(proteomic_latent_axes, by = "AnimalNum")

behavior_axis_cols <- setdiff(names(behavioral_latent_axes), "AnimalNum")
proteomic_axis_cols <- setdiff(names(proteomic_latent_axes), "AnimalNum")

latent_axis_associations <- expand_grid(BehaviorAxis = behavior_axis_cols, ProteomicAxis = proteomic_axis_cols) %>%
  mutate(
    n = map2_int(BehaviorAxis, ProteomicAxis, ~sum(is.finite(axis_merged[[.x]]) & is.finite(axis_merged[[.y]]))),
    spearman_rho = map2_dbl(BehaviorAxis, ProteomicAxis, ~safe_cor(axis_merged[[.x]], axis_merged[[.y]], "spearman")),
    pearson_r = map2_dbl(BehaviorAxis, ProteomicAxis, ~safe_cor(axis_merged[[.x]], axis_merged[[.y]], "pearson")),
    EvidenceUse = if_else(n < 12, "exploratory_effect_size_small_n", "exploratory_association")
  ) %>%
  arrange(desc(abs(spearman_rho)))

curated_behavior_proteomics_models <- latent_axis_associations %>%
  mutate(
    BiologicalModel = paste(BehaviorAxis, "aligned_with", ProteomicAxis),
    ClaimType = "associative",
    AllowedInterpretation = "Behavioral adaptation axis covaries with a curated hippocampal proteomic module axis.",
    ReviewerRisk = if_else(n < 12, "high_small_n", "medium_exploratory"),
    StableForMainText = FALSE
  )

write_table(behavioral_latent_axes, file.path(output_dir, "tables", "behavioral_latent_axes.csv"))
write_table(proteomic_latent_axes, file.path(output_dir, "tables", "proteomic_latent_axes.csv"))
write_table(latent_axis_associations, file.path(output_dir, "tables", "latent_axis_associations.csv"))
write_table(curated_behavior_proteomics_models, file.path(output_dir, "tables", "curated_behavior_proteomics_models.csv"))
write_table(axis_feature_inventory, file.path(output_dir, "tables", "latent_axis_feature_inventory.csv"))

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
