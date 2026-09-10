# ==============================================================================
# Script Name: Fig3B_PatchSeq_UMAP_Mapping.R
# Description: Generates Figure 3B: UMAP visualization of 281 human Patch-seq 
#              neurons with cell-type prediction mapping score (Glu.N vs GABA.N).
# Author: Auto-refactored for EpiAdo STM Revision
# Date: 2026-09-08
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(ggplot2)
  library(RColorBrewer)
  library(dplyr)
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

# 2. Extract UMAP coordinates and cell-type prediction scores
umap_coords <- as.data.frame(Embeddings(ps_obj, "umap"))
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

sourcedata_fig3b <- data.frame(
  Cell_ID = Cells(ps_obj),
  UMAP_1 = round(umap_coords$UMAP_1, 4),
  UMAP_2 = round(umap_coords$UMAP_2, 4),
  Predicted_Celltype = ps_obj$predicted.id,
  Prediction_Score_GluN = round(ps_obj$prediction.score.Glu.N, 4),
  Prediction_Score_GABAN = round(ps_obj$prediction.score.GABA.N, 4),
  stringsAsFactors = FALSE
)

# Export SourceData
csv_file <- file.path(out_dir_data, "Fig3B_PatchSeq_Cell_Mapping_Coordinates.csv")
write.csv(sourcedata_fig3b, csv_file, row.names = FALSE)
cat("SourceData exported to:", csv_file, "\n")
cat("Total Patch-seq cells:", nrow(sourcedata_fig3b), "\n")
print(table(sourcedata_fig3b$Predicted_Celltype))

# 3. Generate Plot matching Figure 3B
# Color gradient from PiYG (green to pink)
piyg_cols <- rev(brewer.pal(11, "PiYG")[3:9])

p3b <- ggplot(sourcedata_fig3b, aes(x = UMAP_1, y = UMAP_2, color = Prediction_Score_GluN)) +
  geom_point(size = 2.2, alpha = 0.9) +
  scale_color_gradientn(
    colours = piyg_cols,
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = c("GABA.N\n0", "0.5", "Glu.N\n1"),
    guide = guide_colorbar(
      title = "Mapping Cell-type Score",
      title.position = "top",
      title.hjust = 0.5,
      barwidth = unit(4, "cm"),
      barheight = unit(0.4, "cm"),
      ticks.colour = "black",
      frame.colour = "black",
      frame.linewidth = 0.5
    )
  ) +
  labs(
    title = "Patch-seq: Total 281 Cells",
    x = "UMAP1",
    y = "UMAP2"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold", color = "black"),
    legend.text = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 11, color = "black"),
    axis.text = element_text(size = 9, color = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length = unit(1.5, "mm"),
    aspect.ratio = 0.85
  )

check_plot_file <- file.path(out_dir_plot, "check_Fig3B_PatchSeq_UMAP_Mapping.png")
ggsave(check_plot_file, p3b, width = 4.5, height = 4.5, dpi = 300)
cat("Check plot saved to:", check_plot_file, "\n")