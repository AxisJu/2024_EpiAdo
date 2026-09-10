# ==============================================================================
# Script Name: Fig7G_snRNA_UMAP_Mouse_DBS.R
# Description: Generates Figure 7G: UMAP visualization of 31,001 single nuclei 
#              from mouse anterior thalamic nucleus (ANT) across Sham, DBS Contra., and DBS groups.
# Author: Auto-refactored for EpiAdo STM Revision
# Date: 2026-09-08
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

# Define paths
dbs_seurat_path <- "Z:/2023_EpiAdo/data/snSeq_JSH/processed/merged_ps_log_pca_hrm_final_seuratobj"
out_dir_data <- "../sourcedata/Figure7"
out_dir_plot <- "../sourcedata/check_plots"

dir.create(out_dir_data, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_plot, recursive = TRUE, showWarnings = FALSE)

# 1. Load Seurat object
cat("Loading mouse DBS Seurat object...\n")
obj <- qread(dbs_seurat_path)

umap_coords <- as.data.frame(Embeddings(obj, "umap"))
colnames(umap_coords) <- c("UMAP_1", "UMAP_2")

# Harmonize celltype labels with manuscript
celltype_map <- c(
  "VGLUT2+ Glu.N" = "Vglut2+ Glu.N",
  "Oth Glu.N" = "Vglut1+ Glu.N",
  "GABA.N" = "GABA.N",
  "OPCs" = "OPCs",
  "Oligo." = "Oligo.",
  "Microglia" = "Microglia",
  "Astro." = "Astro.",
  "Endo." = "Endo.",
  "Peri." = "Peri.",
  "Fb." = "Fb.",
  "Epi." = "Epi.",
  "LowQuality" = "LowQuality"
)

group_map <- c(
  "Ctrl_2" = "Sham",
  "Ctrl_1" = "DBS Contra.",
  "DBS_1" = "DBS"
)

cell_types <- celltype_map[as.character(obj$celltype_coarse)]
conditions <- group_map[as.character(obj$Sample)]

# 2. Export cell-level SourceData
sourcedata_fig7g <- data.frame(
  Cell_ID = Cells(obj),
  UMAP_1 = round(umap_coords$UMAP_1, 4),
  UMAP_2 = round(umap_coords$UMAP_2, 4),
  Celltype = as.character(cell_types),
  Sample = as.character(obj$Sample),
  Condition = as.character(conditions),
  stringsAsFactors = FALSE
)

csv_file_full <- file.path(out_dir_data, "Fig7G_snRNA_UMAP_Coordinates.csv")
write.csv(sourcedata_fig7g, csv_file_full, row.names = FALSE)
cat("Full UMAP SourceData exported to:", csv_file_full, "\n")
cat("Total cells:", nrow(sourcedata_fig7g), "\n")

# Summary composition table
comp_summary <- sourcedata_fig7g %>%
  filter(Celltype != "LowQuality") %>%
  group_by(Condition, Celltype) %>%
  summarise(Cell_Count = n(), .groups = "drop") %>%
  group_by(Condition) %>%
  mutate(Proportion_Percent = round(Cell_Count / sum(Cell_Count) * 100, 2)) %>%
  ungroup()

csv_file_sum <- file.path(out_dir_data, "Fig7G_Cell_Composition_Summary.csv")
write.csv(comp_summary, csv_file_sum, row.names = FALSE)
cat("Composition SourceData exported to:", csv_file_sum, "\n")

# 3. Generate Plot matching Figure 7G
colors_7g <- c(
  "Vglut2+ Glu.N" = "#bd86b6",
  "Vglut1+ Glu.N" = "#e2cbdf",
  "GABA.N" = "#9dcb84",
  "OPCs" = "#cea27b",
  "Oligo." = "#f3d490",
  "Microglia" = "#d25f61",
  "Astro." = "#adaed3",
  "Endo." = "#adc1d3",
  "Peri." = "#9bbdbd",
  "Fb." = "#add3d2",
  "Epi." = "#c5d3e0",
  "LowQuality" = "#EEEEEE"
)

plot_df <- sourcedata_fig7g %>%
  filter(Celltype != "LowQuality")

plot_df$Celltype <- factor(plot_df$Celltype, levels = c(
  "Vglut2+ Glu.N", "Vglut1+ Glu.N", "GABA.N",
  "OPCs", "Oligo.", "Microglia",
  "Astro.", "Endo.", "Peri.", "Fb.", "Epi."
))

p7g <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = Celltype)) +
  geom_point(size = 0.7, alpha = 0.8, stroke = 0) +
  scale_color_manual(values = colors_7g) +
  labs(
    title = "snRNA-Seq: Total 31,001 cells",
    x = "UMAP1",
    y = "UMAP2"
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3)) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5, color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    aspect.ratio = 0.9
  )

check_plot_file <- file.path(out_dir_plot, "check_Fig7G_snRNA_UMAP_Mouse_DBS.png")
ggsave(check_plot_file, p7g, width = 5.2, height = 5.5, dpi = 300)
cat("Check plot saved to:", check_plot_file, "\n")