# ==============================================================================
# Script Name: Fig3C_PatchSeq_Marker_BubblePlot.R
# Description: Generates Figure 3C: Bubble plot of marker genes across human 
#              Patch-seq Glu.N and GABA.N neurons (SLC17A7, RBFOX3, CAMK2A, SYN3, GAD1, GAD2).
# Author: Auto-refactored for EpiAdo STM Revision
# Date: 2026-09-08
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(ggplot2)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
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

markers <- c("SLC17A7", "RBFOX3", "CAMK2A", "SYN3", "GAD1", "GAD2")

# Calculate dotplot statistics manually for exact reproducibility and clean SourceData
exp_mat <- GetAssayData(ps_obj, assay = "RNA", slot = "data")[markers, ]
cell_groups <- ps_obj$predicted.id

dot_data_list <- list()
for (gene in markers) {
  for (grp in c("Glu.N", "GABA.N")) {
    cells_in_grp <- names(cell_groups)[cell_groups == grp]
    vals <- exp_mat[gene, cells_in_grp]
    pct <- mean(vals > 0) * 100
    mean_exp <- mean(expm1(vals)) # natural average expression before log
    mean_log_exp <- mean(vals)
    dot_data_list[[length(dot_data_list) + 1]] <- data.frame(
      Gene = gene,
      Celltype = grp,
      Mean_Expression = round(mean_log_exp, 4),
      Percent_Expressed = round(pct, 2),
      stringsAsFactors = FALSE
    )
  }
}

dot_df <- bind_rows(dot_data_list)

# Scale mean expression per gene across cell types
dot_df <- dot_df %>%
  group_by(Gene) %>%
  mutate(Scaled_Expression = round(as.numeric(scale(Mean_Expression)), 4)) %>%
  ungroup()

# Export SourceData
csv_file <- file.path(out_dir_data, "Fig3C_PatchSeq_Marker_BubblePlot.csv")
write.csv(dot_df, csv_file, row.names = FALSE)
cat("SourceData exported to:", csv_file, "\n")

# Order factors
dot_df$Gene <- factor(dot_df$Gene, levels = markers)
dot_df$Celltype <- factor(dot_df$Celltype, levels = c("GABA.N", "Glu.N"))

piyg_palette <- rev(brewer.pal(11, "PiYG")[2:7])

p3c <- ggplot(dot_df, aes(x = Gene, y = Celltype)) +
  geom_point(aes(size = Percent_Expressed, fill = Scaled_Expression), 
             shape = 21, color = "black", stroke = 0.4) +
  scale_fill_gradientn(
    colours = piyg_palette,
    name = "Average Expression",
    limits = c(-1, 1),
    breaks = c(-0.4, 0, 0.4),
    guide = guide_colorbar(
      title.position = "top",
      barwidth = unit(2.5, "cm"),
      barheight = unit(0.3, "cm"),
      ticks.colour = "black",
      frame.colour = "black",
      frame.linewidth = 0.4
    )
  ) +
  scale_size(
    range = c(2.5, 6.5),
    breaks = c(40, 60, 80),
    name = "Percent Expressed",
    guide = guide_legend(title.position = "top")
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 9, color = "black"),
    legend.text = element_text(size = 8, color = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length = unit(1.5, "mm"),
    aspect.ratio = 0.45
  )

check_plot_file <- file.path(out_dir_plot, "check_Fig3C_PatchSeq_Marker_BubblePlot.png")
ggsave(check_plot_file, p3c, width = 5, height = 3.8, dpi = 300)
cat("Check plot saved to:", check_plot_file, "\n")