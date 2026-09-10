# ==============================================================================
# Script Name: FigS1A_snRNA_UMAP_Region_Pathology.R
# Description: UMAP of snRNA-seq cell distribution across brain regions and
#              pathological types (Supplementary Figure 1A).
# ==============================================================================

library(Seurat)
library(qs)
library(ggplot2)
library(dplyr)
library(patchwork)

# Set seed
set.seed(42)

# Paths
obj_path <- "Z:/2023_EpiAdo/data/scSeq_MSH/Merged/merged_ps_log_pca_harmony_umap_final_seuratobj"
out_data_dir <- "../sourcedata/SuppFig1"
out_check_dir <- "../sourcedata/check_plots"

dir.create(out_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_check_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading Seurat object...\n")
merged_seuratobj <- qread(obj_path)

# Prepare metadata columns
merged_seuratobj$Diseasei3 <- ifelse(
  merged_seuratobj$Group == "Control", "Control",
  ifelse(
    merged_seuratobj$Disease %in% c("Epilepsy, FCD2a", "Epilepsy, FCD2b"),
    "FCD2a & FCD2b",
    "MTLE & LEAT"
  )
)
merged_seuratobj$Diseasei3 <- factor(
  merged_seuratobj$Diseasei3,
  levels = c("Control", "FCD2a & FCD2b", "MTLE & LEAT")
)

merged_seuratobj$BrainRegion <- factor(
  merged_seuratobj$BrainRegion,
  ordered = TRUE,
  levels = c("Cortex", "Hippocampus", "Amygdala")
)

# Palette (10 cell types)
color_res_coarse <- c(
  "#d88d9e", "#c3c88c", "#cea27b", "#f3d490",
  "#d25f61", "#dc8384", "#adaed3", "#adc1d3", "#c5d3e0",
  "#EEEEEE"
)
names(color_res_coarse) <- levels(factor(merged_seuratobj$celltype_coarse))

# Extract UMAP coordinates and metadata for SourceData
umap_df <- as.data.frame(Embeddings(merged_seuratobj, reduction = "umap"))
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$CellID <- rownames(umap_df)
umap_df$celltype_coarse <- as.character(merged_seuratobj$celltype_coarse)
umap_df$BrainRegion <- as.character(merged_seuratobj$BrainRegion)
umap_df$Disease <- as.character(merged_seuratobj$Disease)
umap_df$DiseaseGroup <- as.character(merged_seuratobj$Diseasei3)
umap_df$Group <- as.character(merged_seuratobj$Group)
umap_df$Patient <- as.character(merged_seuratobj$Patient)

# Export SourceData
cat("Exporting SourceData...\n")
write.csv(
  umap_df,
  file = file.path(out_data_dir, "FigS1A_snRNA_UMAP_Coordinates.csv"),
  row.names = FALSE
)

# Export Summary Composition
comp_summary <- umap_df %>%
  group_by(BrainRegion, DiseaseGroup, celltype_coarse) %>%
  summarise(CellCount = n(), .groups = "drop") %>%
  group_by(BrainRegion, DiseaseGroup) %>%
  mutate(Proportion = CellCount / sum(CellCount))

write.csv(
  comp_summary,
  file = file.path(out_data_dir, "FigS1A_Cell_Distribution_Summary.csv"),
  row.names = FALSE
)

cat("Plotting Check Figures...\n")
# Subsample 50000 cells for crisp check plot
plot_df <- umap_df %>%
  sample_n(min(nrow(umap_df), 50000))

# Plot Brain Region (Top row)
p_region <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = celltype_coarse)) +
  geom_point(size = 0.15, alpha = 0.6) +
  facet_wrap(~BrainRegion, ncol = 3) +
  scale_color_manual(values = color_res_coarse) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    axis.text = element_blank(),
    axis.ticks.length = unit(0, "pt"),
    legend.position = "none"
  ) +
  labs(x = "UMAP1", y = "UMAP2", title = "Brain Region")

# Plot Disease (Bottom row)
p_disease <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = celltype_coarse)) +
  geom_point(size = 0.15, alpha = 0.6) +
  facet_wrap(~DiseaseGroup, ncol = 3) +
  scale_color_manual(values = color_res_coarse) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    axis.text = element_blank(),
    axis.ticks.length = unit(0, "pt"),
    legend.position = "none"
  ) +
  labs(x = "UMAP1", y = "UMAP2", title = "Disease")

p_combined <- p_region / p_disease
ggsave(
  filename = file.path(out_check_dir, "check_FigS1A_UMAP_Region_Pathology.png"),
  plot = p_combined,
  width = 9,
  height = 6,
  dpi = 300
)

cat("FigS1A completed successfully.\n")
