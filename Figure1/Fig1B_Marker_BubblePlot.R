# ==============================================================================
# Script: Fig1B_Marker_BubblePlot.R
# Purpose: Bubble plot of scaled logCPM for canonical marker genes across main cell types
# Article: EpiAdo (Science Translational Medicine) - Figure 1B
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
})

out_dir_sd <- "../sourcedata/Figure1"
out_dir_plot <- "../sourcedata/check_plots"
dir.create(out_dir_sd, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_plot, showWarnings = FALSE, recursive = TRUE)

cat("Loading Seurat object for Figure 1B...\n")
# [DATA REQUIRED] The line below originally loaded 'merged_ps_log_pca_harmony_umap_final_seuratobj' from a private path.
# Please provide 'merged_seuratobj' by loading the appropriate data before running this script.
# Example: merged_seuratobj <- qread('path/to/merged_ps_log_pca_harmony_umap_final_seuratobj')
merged_seuratobj <- subset(merged_seuratobj, subset = celltype_coarse == "LowQuality", invert = TRUE)
Idents(merged_seuratobj) <- "celltype_coarse"

target_levels <- c(
  "Glu.N", "GABA.N", "OPCs", "Oligo.", "Mo/Mφ", "Microglia", "Astro.", "Endo.", "Epi."
)
Idents(merged_seuratobj) <- factor(Idents(merged_seuratobj), levels = rev(target_levels))

marker_genes <- c(
  "SLC17A7", "SLC17A6", "CAMK2A", "RBFOX3", "SYN3",
  "GAD1", "GAD2",
  "VCAN", "PDGFRA",
  "MOBP", "MOG",
  "CD96", "PTPRC", "DOCK8", "APBB1IP", "P2RY12",
  "FGFR3", "AQP4", "GJA1", "ALDH1L1", "GFAP",
  "FLT1", "VWF", "PECAM1", "CDH5",
  "HTR2C"
)
marker_genes <- intersect(marker_genes, rownames(merged_seuratobj))

cat("Generating DotPlot...\n")
p_dot <- DotPlot(
  merged_seuratobj,
  scale = TRUE,
  scale.by = "radius",
  features = marker_genes
) +
  scale_color_gradientn(colours = rev(brewer.pal(11, "PiYG")[c(2:10)])) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, colour = "black"),
    axis.text.y = element_text(size = 11, colour = "black"),
    axis.title = element_blank(),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5)
  ) +
  labs(title = "Marker Expression Across Cell Types")

# Extract SourceData directly from ggplot DotPlot data
cat("Exporting SourceData for Figure 1B...\n")
dot_data <- p_dot$data %>%
  dplyr::select(
    Gene = features.plot,
    CellType = id,
    AverageExpression = avg.exp,
    ScaledAverageExpression = avg.exp.scaled,
    PercentExpressed = pct.exp
  )
write.csv(dot_data, file.path(out_dir_sd, "Fig1B_Marker_BubblePlot.csv"), row.names = FALSE)

ggsave(file.path(out_dir_plot, "check_Fig1B_Marker_BubblePlot.png"), plot = p_dot, width = 10, height = 4.5, dpi = 300)
cat("Figure 1B finished successfully.\n")
