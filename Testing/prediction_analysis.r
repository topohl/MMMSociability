# ==============================================================================
# Resilience Prediction: Sex-Specific Analysis (Nature Style)
# ==============================================================================

library(tidyverse)
library(pROC)
library(patchwork)

# 1. Load Data
data_file <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability/processed_data/data_lme_format/data_filtered_agg.csv"
output_dir <- "S:/Lab_Member/Tobi/Experiments/Exp9_Social-Stress/Analysis/Behavior/RFID/MMMSociability/statistics/prediction/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

df <- read_csv(data_file, show_col_types = FALSE)

# 2. Nature-Style Theme Helper
theme_nature <- function() {
    theme_classic(base_size = 7) +
        theme(
            text = element_text(family = "sans", color = "black"),
            axis.text = element_text(size = 7, color = "black"),
            axis.title = element_text(size = 8, face = "bold"),
            plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
            panel.grid = element_blank(),
            axis.line = element_line(linewidth = 0.3)
        )
}

# 3. Generalized Processing Function
run_prediction_comparison <- function(target_sex, plot_color, compare_groups, compare_label) {

    d_sex <- df %>%
        filter(
            toupper(Sex) == toupper(target_sex),
            Change == "CC1",
            Phase == "Active",
            HalfHourElapsed <= 23,
            Group %in% c("RES", compare_groups),
            !is.na(Movement)
        ) %>%
        mutate(Group2 = if_else(Group == "RES", "RES", compare_label)) %>%
        group_by(AnimalNum, Group2) %>%
        summarise(mean_act = mean(Movement, na.rm = TRUE), .groups = "drop") %>%
        mutate(label = factor(Group2, levels = c(compare_label, "RES")))

    if (nrow(d_sex) < 5 || length(unique(d_sex$label)) < 2) return(NULL)

    roc_obj <- roc(d_sex$label, d_sex$mean_act, direction = ">")
    best_stats <- coords(roc_obj, "best", ret = c("threshold", "specificity", "sensitivity"))

    other_color <- case_when(
        compare_label == "SUS" ~ "#e63947",
        compare_label == "CON" ~ "#7f7f7f",
        compare_label == "OTHER" ~ "#555555",
        TRUE ~ "#7f7f7f"
    )

    p_dist <- ggplot(d_sex, aes(x = label, y = mean_act, color = label)) +
        geom_boxplot(aes(fill = label), alpha = 0.1, outlier.shape = NA, linewidth = 0.3) +
        geom_jitter(width = 0.15, size = 1.2, alpha = 0.7) +
        geom_hline(yintercept = best_stats$threshold, linetype = "dashed", linewidth = 0.3) +
        scale_color_manual(values = c(compare_label = other_color, "RES" = plot_color)) +
        scale_fill_manual(values = c(compare_label = other_color, "RES" = plot_color)) +
        labs(
            title = paste(target_sex, "RES vs", compare_label),
            x = "Outcome",
            y = "Mean Movement"
        ) +
        theme_nature() +
        theme(legend.position = "none")

    p_roc <- ggroc(roc_obj, color = plot_color, linewidth = 0.8) +
        geom_abline(slope = 1, intercept = 1, linetype = "dotted", color = "grey") +
        annotate(
            "text",
            x = 0.3, y = 0.15,
            label = paste("AUC =", round(auc(roc_obj), 3)),
            size = 3, fontface = "bold"
        ) +
        labs(
            title = paste(target_sex, "Prediction"),
            x = "Specificity",
            y = "Sensitivity"
        ) +
        theme_nature()

    p_dist + p_roc + plot_annotation(tag_levels = "a")
}

# 4. Helper to Save Figures
save_prediction_set <- function(target_sex, sex_color, sex_name) {

    comparisons <- list(
        list(groups = c("SUS"), label = "SUS"),
        list(groups = c("CON"), label = "CON"),
        list(groups = c("SUS", "CON"), label = "OTHER")
    )

    for (cmp in comparisons) {
        fig <- run_prediction_comparison(
            target_sex = target_sex,
            plot_color = sex_color,
            compare_groups = cmp$groups,
            compare_label = cmp$label
        )

        if (!is.null(fig)) {
            base_name <- paste0("Prediction_", sex_name, "_RES_vs_", cmp$label)

            ggsave(
                file.path(output_dir, paste0(base_name, ".svg")),
                fig, width = 100, height = 50, units = "mm"
            )
        }
    }
}

# 5. Generate and Save Figures
save_prediction_set("f", "#3d3b6e", "Females")
save_prediction_set("m", "#2a9d8f", "Males")

cat("\n--- Sex-Specific Prediction Complete ---\n")