# ==============================================================================
# Script Name: Fig3F_GluN_EFeatures_ADORA1_Low_vs_High.R
# Description: Generates Figure 3F: Comparison of electrophysiological features 
#              (UpStroke, Up/Down Ratio, AP HalfWidth, AP Duration) between 
#              ADORA1-Low and ADORA1-High human epileptic Glu.N neurons.
# Author: Auto-refactored for EpiAdo STM Revision
# Date: 2026-09-08
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(rstatix)
})

# Define paths
ps_seurat_path <- "Z:/2023_EpiAdo/data/PatchSeq_CM/merged/ps_itg_seuratobj"
ps_cpm_path <- "Z:/2023_EpiAdo/data/PatchSeq_CM/merged/ps_itg_logcpm"
out_dir_data <- "../sourcedata/Figure3"
out_dir_plot <- "../sourcedata/check_plots"

dir.create(out_dir_data, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_plot, recursive = TRUE, showWarnings = FALSE)

# 1. Load data
cat("Loading Patch-seq data...\n")
ps_obj <- qread(ps_seurat_path)
ps_cpm <- qread(ps_cpm_path)

properties <- c("AP_UpStroke", "AP_UpDownRatio", "AP_HalfWidth", "AP_Duration")
prop_labels <- c(
  "AP_UpStroke" = "UpStroke",
  "AP_UpDownRatio" = "Up/Down Ratio",
  "AP_HalfWidth" = "AP HalfWidth",
  "AP_Duration" = "AP Duration (mV)"
)

# Extract filtered cohort: Epilepsy, Glu.N, HS cohort
meta_df <- ps_obj@meta.data
selected_cells <- meta_df %>%
  tibble::rownames_to_column("Cell_ID") %>%
  filter(Group == "Epilepsy", predicted.id == "Glu.N", Cohort %in% c("HS")) %>%
  mutate(ADORA1_Exp = ps_cpm["ADORA1", Cell_ID]) %>%
  select(Cell_ID, Group, Cohort, ADORA1_Exp, all_of(properties))

exp_threshold <- mean(selected_cells$ADORA1_Exp)
selected_cells$Ado_Group <- ifelse(selected_cells$ADORA1_Exp > exp_threshold, "ADORA1High", "ADORA1Low")

# Export cell-level SourceData
csv_file_data <- file.path(out_dir_data, "Fig3F_GluN_EFeatures_ADORA1_Low_vs_High.csv")
write.csv(selected_cells, csv_file_data, row.names = FALSE)
cat("Cell-level SourceData exported to:", csv_file_data, "\n")
cat("Sample counts:", table(selected_cells$Ado_Group), "\n")

# Calculate statistics and build subplots
color_map <- c("ADORA1High" = "#386CB0", "ADORA1Low" = "#7FBC41")
stat_records <- list()
plot_list <- list()

for (prop in properties) {
  sub_df <- selected_cells %>%
    select(Cell_ID, Ado_Group, Feature_Val = all_of(prop)) %>%
    filter(!is.na(Feature_Val))
  
  t_res <- t_test(sub_df, Feature_Val ~ Ado_Group, paired = FALSE)
  
  p_val <- t_res$p
  p_label <- paste0("P = ", signif(p_val, 2))
  
  # Summary stats for export
  sum_stats <- sub_df %>%
    group_by(Ado_Group) %>%
    summarise(
      N = n(),
      Mean = mean(Feature_Val),
      SD = sd(Feature_Val),
      Median = median(Feature_Val),
      IQR = IQR(Feature_Val),
      .groups = "drop"
    )
  
  stat_records[[prop]] <- data.frame(
    Feature = prop_labels[prop],
    N_High = sum_stats$N[sum_stats$Ado_Group == "ADORA1High"],
    Mean_High = round(sum_stats$Mean[sum_stats$Ado_Group == "ADORA1High"], 4),
    SD_High = round(sum_stats$SD[sum_stats$Ado_Group == "ADORA1High"], 4),
    N_Low = sum_stats$N[sum_stats$Ado_Group == "ADORA1Low"],
    Mean_Low = round(sum_stats$Mean[sum_stats$Ado_Group == "ADORA1Low"], 4),
    SD_Low = round(sum_stats$SD[sum_stats$Ado_Group == "ADORA1Low"], 4),
    T_statistic = round(t_res$statistic, 4),
    p_value = signif(p_val, 4),
    stringsAsFactors = FALSE
  )
  
  set.seed(42)
  sub_df$Jitter_X <- ifelse(sub_df$Ado_Group == "ADORA1High", 1.25, 2.25) + runif(nrow(sub_df), -0.07, 0.07)
  
  max_y <- max(sub_df$Feature_Val, na.rm = TRUE)
  min_y <- min(sub_df$Feature_Val, na.rm = TRUE)
  y_range <- max_y - min_y
  bracket_y <- max_y + y_range * 0.08
  text_y <- max_y + y_range * 0.15
  
  p <- ggplot(sub_df, aes(x = Ado_Group, y = Feature_Val)) +
    geom_violin(aes(fill = Ado_Group, color = Ado_Group), width = 0.55, alpha = 0.35, position = position_nudge(x = -0.15)) +
    geom_boxplot(aes(color = Ado_Group), width = 0.18, fill = "white", outlier.shape = NA, 
                 linewidth = 0.5, position = position_nudge(x = -0.15)) +
    geom_point(aes(x = Jitter_X, y = Feature_Val, fill = Ado_Group), shape = 21, color = "black", 
               size = 1.8, stroke = 0.3, alpha = 0.85) +
    annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y, linewidth = 0.4) +
    annotate("segment", x = 1, xend = 1, y = bracket_y - y_range * 0.03, yend = bracket_y, linewidth = 0.4) +
    annotate("segment", x = 2, xend = 2, y = bracket_y - y_range * 0.03, yend = bracket_y, linewidth = 0.4) +
    annotate("text", x = 1.5, y = text_y, label = p_label, size = 3, family = "sans") +
    scale_fill_manual(values = color_map) +
    scale_color_manual(values = color_map) +
    scale_x_discrete(labels = c("ADORA1High", "ADORA1Low")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
    labs(title = prop_labels[prop], x = NULL, y = NULL) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5, color = "black"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color = "black", size = 9),
      axis.line = element_line(colour = "black", linewidth = 0.5),
      axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
      axis.ticks.length = unit(1.5, "mm"),
      legend.position = "none"
    )
  
  plot_list[[prop]] <- p
}

# Export statistics
stat_df <- bind_rows(stat_records)
csv_file_stat <- file.path(out_dir_data, "Fig3F_GluN_EFeatures_ADORA1_Low_vs_High_Statistics.csv")
write.csv(stat_df, csv_file_stat, row.names = FALSE)
cat("Statistics SourceData exported to:", csv_file_stat, "\n")
print(stat_df)

# Assemble 4 panels horizontally
title_header <- ggdraw() + 
  draw_label(expression(bold("Comparison of E-Features in ") * bolditalic("ADORA1")^bold("Low") * bold(" V.S. ") * bolditalic("ADORA1")^bold("High") * bold(" Epilepsy Glu.N")), 
             size = 11, x = 0.5, hjust = 0.5)

legend_plot <- ggplot(data.frame(x = c(1, 2), y = c(1, 1), Ado_Group = c("ADORA1High", "ADORA1Low")), 
                      aes(x = x, y = y, fill = Ado_Group)) +
  geom_point(shape = 22, size = 4) +
  scale_fill_manual(values = color_map, labels = c(
    expression(italic(ADORA1)^High),
    expression(italic(ADORA1)^Low)
  )) +
  theme_void() +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 9))

plots_grid <- plot_grid(plotlist = plot_list, nrow = 1, align = "h")
legend_box <- get_legend(legend_plot)

final_p3f <- plot_grid(
  title_header,
  plots_grid,
  legend_box,
  ncol = 1,
  rel_heights = c(0.12, 1, 0.12)
)

check_plot_file <- file.path(out_dir_plot, "check_Fig3F_GluN_EFeatures_ADORA1_Low_vs_High.png")
ggsave(check_plot_file, final_p3f, width = 7.5, height = 3.6, dpi = 300)
cat("Check plot saved to:", check_plot_file, "\n")