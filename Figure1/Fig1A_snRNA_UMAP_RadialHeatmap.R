# ==============================================================================
# Script: Fig1A_snRNA_UMAP_RadialHeatmap.R
# Purpose: Epileptic human brain atlas (UMAP 421,714 cells) & Radial Heatmap
# Article: EpiAdo (Science Translational Medicine) - Figure 1A
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(circlize)
})

out_dir_sd <- "../sourcedata/Figure1"
out_dir_plot <- "../sourcedata/check_plots"
dir.create(out_dir_sd, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_plot, showWarnings = FALSE, recursive = TRUE)

cat("Loading Seurat object for Figure 1A...\n")
# [DATA REQUIRED] The line below originally loaded 'merged_ps_log_pca_harmony_umap_final_seuratobj' from a private path.
# Please provide 'merged_seuratobj' by loading the appropriate data before running this script.
# Example: merged_seuratobj <- qread('path/to/merged_ps_log_pca_harmony_umap_final_seuratobj')
merged_seuratobj <- subset(merged_seuratobj, subset = celltype_coarse == "LowQuality", invert = TRUE)

# Color palettes
color_res.coarse <- c(
  "#d0a9cb", "#9dcb84", "#cea27b", "#f3d490",
  "#d25f61", "#dc8384", "#adaed3", "#adc1d3", "#c5d3e0"
)
names(color_res.coarse) <- c(
  "Glu.N", "GABA.N", "OPCs", "Oligo.",
  "Mo/Mφ", "Microglia", "Astro.", "Endo.", "Epi."
)

# 1. Export SourceData for 1A
cat("Exporting SourceData for Figure 1A...\n")
umap_coords <- Embeddings(merged_seuratobj, "umap")
meta <- merged_seuratobj@meta.data

sd_umap <- data.frame(
  CellID = rownames(meta),
  UMAP_1 = round(umap_coords[, 1], 3),
  UMAP_2 = round(umap_coords[, 2], 3),
  CellType_Coarse = meta$celltype_coarse,
  CellType_Fine = meta$celltype_fine,
  Group = meta$Group,
  Gender = meta$Gender,
  BrainRegion = meta$BrainRegion,
  Cohort = meta$Cohort,
  Patient = meta$Patient
)
write.csv(sd_umap[seq(1, nrow(sd_umap), by = 5), ], 
          file.path(out_dir_sd, "Fig1A_snRNA_UMAP_Coordinates_Subsample.csv"), row.names = FALSE)

# Composition summary
comp_summary <- meta %>%
  dplyr::group_by(celltype_coarse, Group, BrainRegion, Cohort) %>%
  dplyr::summarise(CellCount = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(celltype_coarse) %>%
  dplyr::mutate(Percentage = round(CellCount / sum(CellCount) * 100, 2))
write.csv(comp_summary, file.path(out_dir_sd, "Fig1A_Sample_Composition_Summary.csv"), row.names = FALSE)

# 2. Generate Check Plots
cat("Generating Check Plots for Figure 1A...\n")
p_umap <- DimPlot(
  merged_seuratobj,
  group.by = "celltype_coarse",
  reduction = "umap",
  pt.size = 0.5,
  cols = color_res.coarse,
  raster = TRUE, raster.dpi = c(1024, 1024)
) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "right"
  ) +
  labs(title = "Epileptic human brain atlas (snRNA-seq: 421,714 cells)")

ggsave(file.path(out_dir_plot, "check_Fig1A_UMAP.png"), plot = p_umap, width = 8, height = 7, dpi = 300)
cat("Figure 1A finished successfully.\n")
