# ==============================================================================
# Script Name: Fig3G_Clinical_SeizureBurden_ADORA1_Correlation.R
# Description: Generates Figure 3G: Pearson correlation between pseudo-bulked 
#              ADORA1 expression in Glu.N or GABA.N neurons and clinical features 
#              (OnsetToSurgery Time and Calculated Total Seizure Burden) in epileptic patients.
# Author: Auto-refactored for EpiAdo STM Revision
# Date: 2026-09-08
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(Seurat)
  library(ggplot2)
  library(RColorBrewer)
  library(cowplot)
})

# Define paths
ps_seurat_path <- "Z:/2023_EpiAdo/data/PatchSeq_CM/merged/ps_itg_seuratobj"
out_dir_data <- "../sourcedata/Figure3"
out_dir_plot <- "../sourcedata/check_plots"

dir.create(out_dir_data, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_plot, recursive = TRUE, showWarnings = FALSE)

# 1. Load Patch-seq Seurat object
cat("Loading Patch-seq Seurat object...\n")
ps_obj <- qread(ps_seurat_path)

# Calculate Total Seizure Burden formula
ps_obj[["TotalSeizureBurden.calc"]] <- log(ps_obj$OTS_Time * ps_obj$GS_Freq * 365 + 
                                          ps_obj$OTS_Time * ps_obj$FS_Freq * 365 + 1)

# SCTransform normalization for scale.data ADORA1
cat("Running SCTransform normalization...\n")
if (!"SCT" %in% names(ps_obj@assays)) {
  ps_obj <- SCTransform(ps_obj, return.only.var.genes = FALSE, vars.to.regress = c("nCount_RNA"), verbose = FALSE)
}

adora1_scale <- GetAssayData(ps_obj, assay = "SCT", slot = "scale.data")["ADORA1", ]

# 2. Extract and summarize data for Glu.N and GABA.N
get_clinical_df <- function(target_celltype) {
  df <- data.frame(
    Sample = ps_obj$Sample_Patchseq,
    Celltype = ps_obj$predicted.id,
    Group = ps_obj$Group,
    Disease = ps_obj$Disease,
    Cohort = ps_obj$Cohort,
    OTS_Time = ps_obj$OTS_Time,
    TotalSeizureBurden_calc = ps_obj$TotalSeizureBurden.calc,
    ADORA1_Scale = adora1_scale
  ) %>%
    filter(Celltype == target_celltype, Group == "Epilepsy", Cohort %in% c("HS")) %>%
    filter(!is.na(OTS_Time) & !is.na(TotalSeizureBurden_calc)) %>%
    group_by(Sample, Disease, OTS_Time, TotalSeizureBurden_calc) %>%
    summarise(Mean_ADORA1_Scale = mean(ADORA1_Scale), .groups = "drop") %>%
    mutate(Celltype = target_celltype)
  return(df)
}

glu_clin <- get_clinical_df("Glu.N")
gaba_clin <- get_clinical_df("GABA.N")

all_clin <- bind_rows(glu_clin, gaba_clin)

csv_file_data <- file.path(out_dir_data, "Fig3G_Clinical_SeizureBurden_ADORA1_Correlation.csv")
write.csv(all_clin, csv_file_data, row.names = FALSE)
cat("SourceData exported to:", csv_file_data, "\n")

# 3. Calculate correlations and build subplots
cor_stats <- list()

make_plot <- function(data_sub, clin_col, x_label, title_celltype, show_y_title = TRUE) {
  x_vals <- data_sub[[clin_col]]
  y_vals <- data_sub$Mean_ADORA1_Scale
  
  res <- cor.test(x_vals, y_vals, method = "pearson")
  r_est <- round(res$estimate, 2)
  p_val <- round(res$p.value, 2)
  
  label_text <- paste0("Cor = ", r_est, "  P = ", p_val)
  
  cor_stats[[paste(title_celltype, clin_col, sep = "_")]] <<- data.frame(
    Celltype = title_celltype,
    Clinical_Feature = x_label,
    N_Samples = nrow(data_sub),
    Pearson_R = r_est,
    p_value = p_val,
    stringsAsFactors = FALSE
  )
  
  # Colors
  line_col <- ifelse(title_celltype == "Glu.N", "#d6604d", "#7FBC41")
  fill_col <- ifelse(title_celltype == "Glu.N", "#f4a582", "#b8e186")
  
  p <- ggplot(data_sub, aes(x = .data[[clin_col]], y = Mean_ADORA1_Scale)) +
    geom_smooth(method = "lm", se = TRUE, color = line_col, fill = fill_col, 
                linetype = "dashed", linewidth = 0.8, alpha = 0.25) +
    geom_point(aes(fill = Mean_ADORA1_Scale), shape = 21, color = "black", size = 3, stroke = 0.4) +
    scale_fill_gradientn(
      colours = brewer.pal(9, "Blues")[2:7],
      limits = c(-0.5, 0.5),
      guide = guide_colorbar(
        title = "ADORA1 Scaled Expression",
        title.position = "top",
        barwidth = unit(2.5, "cm"),
        barheight = unit(0.3, "cm"),
        ticks.colour = "black",
        frame.colour = "black",
        frame.linewidth = 0.4
      )
    ) +
    annotate("text", x = min(x_vals) + (max(x_vals) - min(x_vals)) * 0.5, 
             y = max(y_vals) + 0.15, label = label_text, size = 3, family = "sans") +
    labs(
      x = x_label,
      y = if (show_y_title) "Expression of ADORA1 (LogCPM)" else NULL
    ) +
    theme_classic() +
    theme(
      axis.title = element_text(size = 9.5, color = "black"),
      axis.text = element_text(size = 8.5, color = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.5),
      axis.ticks = element_line(colour = "black", linewidth = 0.5),
      axis.ticks.length = unit(1.5, "mm"),
      legend.position = "none"
    )
  return(p)
}

p_glu_ots <- make_plot(glu_clin, "OTS_Time", "OnsetToSurgery Time (yr)", "Glu.N", TRUE)
p_glu_bur <- make_plot(glu_clin, "TotalSeizureBurden_calc", "Calculated Total Seizure Burden", "Glu.N", FALSE)

p_gaba_ots <- make_plot(gaba_clin, "OTS_Time", "OnsetToSurgery Time (yr)", "GABA.N", TRUE)
p_gaba_bur <- make_plot(gaba_clin, "TotalSeizureBurden_calc", "Calculated Total Seizure Burden", "GABA.N", FALSE)

# Export correlation statistics
stat_df <- bind_rows(cor_stats)
csv_file_stat <- file.path(out_dir_data, "Fig3G_Clinical_Correlation_Statistics.csv")
write.csv(stat_df, csv_file_stat, row.names = FALSE)
cat("Statistics SourceData exported to:", csv_file_stat, "\n")
print(stat_df)

# Assemble 4 panels (Glu.N on left, GABA.N on right)
top_header <- ggdraw() + 
  draw_label("Pearson's Correlation of Clinical Features and ADORA1 Expression in Pseudo-bulked Epilepsy Neurons\nTotal Seizure Burden = log2(Generalized Seizure Burden + Focal Seizure Burden + 1)", 
             fontface = "bold", size = 10, hjust = 0.5, x = 0.5)

grid_glu <- plot_grid(p_glu_ots, p_glu_bur, nrow = 1, align = "h")
grid_gaba <- plot_grid(p_gaba_ots, p_gaba_bur, nrow = 1, align = "h")

title_glu <- ggdraw() + draw_label("Glu.N", fontface = "bold", color = "#c51b7d", size = 11, hjust = 0, x = 0.05)
title_gaba <- ggdraw() + draw_label("GABA.N", fontface = "bold", color = "#4d9221", size = 11, hjust = 0, x = 0.05)

col_glu <- plot_grid(title_glu, grid_glu, ncol = 1, rel_heights = c(0.12, 1))
col_gaba <- plot_grid(title_gaba, grid_gaba, ncol = 1, rel_heights = c(0.12, 1))

final_grid <- plot_grid(col_glu, col_gaba, nrow = 1, rel_widths = c(1, 1))

final_p3g <- plot_grid(
  top_header,
  final_grid,
  ncol = 1,
  rel_heights = c(0.18, 1)
)

check_plot_file <- file.path(out_dir_plot, "check_Fig3G_Clinical_SeizureBurden_ADORA1_Correlation.png")
ggsave(check_plot_file, final_p3g, width = 9.5, height = 4.2, dpi = 300)
cat("Check plot saved to:", check_plot_file, "\n")