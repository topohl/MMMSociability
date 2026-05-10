suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(mgcv)
  library(pracma)
})

source("Functions/behavioral_dynamics_helpers.R")

input_file <- "analysis_ready/03_derived_metrics/phase_based/all_behavior_metrics.csv"
output_dir <- "analysis_ready/06_behavioral_dynamics/gamm_features"

ensure_dir(output_dir)
ensure_dir(file.path(output_dir, "tables"))
ensure_dir(file.path(output_dir, "figures"))

raw_dat <- read_behavior_table(input_file)
behav <- standardize_behavior_columns(raw_dat)

metric_names <- c("Movement", "Entropy", "Proximity")
all_features <- list()

for (metric in metric_names) {

  dat <- behav %>%
    transmute(
      AnimalNum,
      Group,
      Sex,
      Phase,
      CageChange,
      TimeIndex,
      Value = .data[[metric]]
    ) %>%
    filter(is.finite(Value))

  if (nrow(dat) < 20) next

  dat <- dat %>% mutate(TimeScaled = as.numeric(scale(TimeIndex)))

  fit <- mgcv::bam(
    Value ~ Group + s(TimeScaled, by = Group, k = 6, bs = 'tp') + s(AnimalNum, bs = 're'),
    data = dat,
    discrete = TRUE,
    method = "fREML"
  )

  dat$Pred <- predict(fit, newdata = dat)

  feature_tbl <- dat %>%
    group_by(Group, Sex, Phase, CageChange, AnimalNum) %>%
    arrange(TimeIndex, .by_group = TRUE) %>%
    summarise(
      Metric = metric,
      auc = pracma::trapz(TimeIndex, Pred),
      mean_pred = mean(Pred, na.rm = TRUE),
      peak = max(Pred, na.rm = TRUE),
      trough = min(Pred, na.rm = TRUE),
      dynamic_range = peak - trough,
      time_to_peak = TimeIndex[which.max(Pred)][1],
      rmssd_pred = calc_rmssd(Pred),
      acf1_pred = calc_acf1(Pred),
      .groups = "drop"
    )

  all_features[[metric]] <- feature_tbl

  write_table(feature_tbl, file.path(output_dir, "tables", paste0(safe_name(metric), "_gamm_features.csv")))
}

combined_features <- bind_rows(all_features)
write_table(combined_features, file.path(output_dir, "tables", "combined_gamm_features.csv"))

message("GAMM trajectory-feature extraction complete.")
