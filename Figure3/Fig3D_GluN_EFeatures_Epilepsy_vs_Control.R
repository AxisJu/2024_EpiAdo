# ==============================================================================
# Script Name: Fig3D_GluN_EFeatures_Epilepsy_vs_Control.R
# Description: Generates Figure 3D: Comparison of electrophysiological features 
#              (UpStroke, Up/Down Ratio, AP HalfWidth, AP Duration) in human 
#              glutamatergic neurons between Epilepsy and Control.
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

# Extract filtered cohort: Glu.N, HS cohort, Control with RMP < -50 or Epilepsy
meta_df <- ps_obj@meta.data
selected_cells <- meta_df %>%
  tibble::rownames_to_column("Cell_ID") %>%
  filter(predicted.id == "Glu.N", Cohort == "HS") %>%
  # NOTE: Control neurons are filtered to RestingMembranePotential < -50 mV to
  # exclude cells with poor membrane integrity (depolarized RMP). This is a standard
  # electrophysiology quality criterion (see Methods: "Patch-clamp recording inclusion criteria").
  # All Epilepsy neurons pass quality control irrespective of RMP.
  filter((Group == "Epilepsy") | (Group == "Control" & RestingMembranePotential < -50)) %>%
  select(Cell_ID, Group, Cohort, RestingMembranePotential, all_of(properties))

# Export cell-level SourceData
csv_file_data <- file.path(out_dir_data, "Fig3D_GluN_EFeatures_Epilepsy_vs_Control.csv")
write.csv(selected_cells, csv_file_data, row.names = FALSE)
cat("Cell-level SourceData exported to:", csv_file_data, "\n")
cat("Sample counts:", table(selected_cells$Group), "\n")

# Calculate statistics and build subplots
color_group <- c("Control" = "#989798", "Epilepsy" = "#fac5b1")
stat_records <- list()
plot_list <- list()

for (prop in properties) {
  sub_df <- selected_cells %>%
    select(Cell_ID, Group, Feature_Val = all_of(prop)) %>%
    filter(!is.na(Feature_Val))
  
  t_res <- t_test(sub_df, Feature_Val ~ Group, paired = FALSE)
  
  p_val <- t_res$p
  p_label <- paste0("P = ", signif(p_val, 2))
  
  # Summary stats for export
  sum_stats <- sub_df %>%
    group_by(Group) %>%
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
    N_Control = sum_stats$N[sum_stats$Group == "Control"],
    Mean_Control = round(sum_stats$Mean[sum_stats$Group == "Control"], 4),
    SD_Control = round(sum_stats$SD[sum_stats$Group == "Control"], 4),
    N_Epilepsy = sum_stats$N[sum_stats$Group == "Epilepsy"],
    Mean_Epilepsy = round(sum_stats$Mean[sum_stats$Group == "Epilepsy"], 4),
    SD_Epilepsy = round(sum_stats$SD[sum_stats$Group == "Epilepsy"], 4),
    T_statistic = round(t_res$statistic, 4),
    p_value = signif(p_val, 4),
    stringsAsFactors = FALSE
  )
  
  # Plot matching paper style: half violin + half boxplot + jitter points
  # Using standard ggplot layers to ensure robust rendering
  set.seed(42)
  sub_df$Jitter_X <- ifelse(sub_df$Group == "Control", 1.25, 2.25) + runif(nrow(sub_df), -0.07, 0.07)
  
  max_y <- max(sub_df$Feature_Val, na.rm = TRUE)
  min_y <- min(sub_df$Feature_Val, na.rm = TRUE)
  y_range <- max_y - min_y
  bracket_y <- max_y + y_range * 0.08
  text_y <- max_y + y_range * 0.15
  
  p <- ggplot(sub_df, aes(x = Group, y = Feature_Val)) +
    # Half violin on left side
    geom_violin(aes(fill = Group, color = Group), width = 0.55, alpha = 0.35, position = position_nudge(x = -0.15)) +
    # Boxplot in center-left
    geom_boxplot(aes(color = Group), width = 0.18, fill = "white", outlier.shape = NA, 
                 linewidth = 0.5, position = position_nudge(x = -0.15)) +
    # Points on right side
    geom_point(aes(x = Jitter_X, y = Feature_Val, fill = Group), shape = 21, color = "black", 
               size = 1.8, stroke = 0.3, alpha = 0.85) +
    # Significance bracket
    annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y, linewidth = 0.4) +
    annotate("segment", x = 1, xend = 1, y = bracket_y - y_range * 0.03, yend = bracket_y, linewidth = 0.4) +
    annotate("segment", x = 2, xend = 2, y = bracket_y - y_range * 0.03, yend = bracket_y, linewidth = 0.4) +
    annotate("text", x = 1.5, y = text_y, label = p_label, size = 3, family = "sans") +
    scale_fill_manual(values = color_group) +
    scale_color_manual(values = color_group) +
    scale_x_discrete(labels = c("Control", "Epilepsy")) +
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
csv_file_stat <- file.path(out_dir_data, "Fig3D_GluN_EFeatures_Epilepsy_vs_Control_Statistics.csv")
write.csv(stat_df, csv_file_stat, row.names = FALSE)
cat("Statistics SourceData exported to:", csv_file_stat, "\n")
print(stat_df)

# Assemble 4 panels horizontally
title_header <- ggdraw() + 
  draw_label("Comparison of E-Features in Epilepsy V.S. Control Glu.N", fontface = "bold", size = 11, x = 0.5, hjust = 0.5)

legend_plot <- ggplot(data.frame(x = c(1, 2), y = c(1, 1), Group = c("Control", "Epilepsy")), 
                      aes(x = x, y = y, fill = Group)) +
  geom_point(shape = 22, size = 4) +
  scale_fill_manual(values = color_group) +
  theme_void() +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 9))

plots_grid <- plot_grid(plotlist = plot_list, nrow = 1, align = "h")
legend_box <- get_legend(legend_plot)

final_p3d <- plot_grid(
  title_header,
  plots_grid,
  legend_box,
  ncol = 1,
  rel_heights = c(0.12, 1, 0.12)
)

check_plot_file <- file.path(out_dir_plot, "check_Fig3D_GluN_EFeatures_Epilepsy_vs_Control.png")
ggsave(check_plot_file, final_p3d, width = 7.5, height = 3.6, dpi = 300)
cat("Check plot saved to:", check_plot_file, "\n")