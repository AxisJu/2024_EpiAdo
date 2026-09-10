# ==============================================================================
# Script: Fig2D_Candidate_Gene_Heatmap_Correlation_CohenD.R
# Purpose: Heatmap of Candidate Genes Correlation with Adenosine Metabolism & Cohen's d (Epi vs Ctrl)
# Article: EpiAdo (Science Translational Medicine) - Figure 2D
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
  library(gridExtra)
  library(ggplot2)
})

out_dir_sd <- "../sourcedata/Figure2"
out_dir_plot <- "../sourcedata/check_plots"
dir.create(out_dir_sd, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_plot, showWarnings = FALSE, recursive = TRUE)

cat("Loading Fig 2D data tables...\n")
# [DATA REQUIRED] The line below originally loaded 'final.csv' from a private path.
# Please provide 'final_df' by loading the appropriate data before running this script.
# Example: final_df <- read.csv('path/to/final.csv')
# [DATA REQUIRED] The line below originally loaded 'gene_catagory.csv' from a private path.
# Please provide 'gene_category' by loading the appropriate data before running this script.
# Example: gene_category <- read.csv('path/to/gene_catagory.csv')

category_order <- c(
  "Synaptic Structure",
  "Synaptic Transmission",
  "Neuronal Differentiation",
  "Epigenetic Regulation",
  "Cell Signaling",
  "Mitochondrial Function",
  "Membrane Function",
  "Protein Process",
  "Metabolic Process",
  "Others"
)

# Correlation matrix
heat_df <- final_df %>%
  dplyr::select(Gene, dplyr::starts_with("Correlation_")) %>%
  tibble::column_to_rownames("Gene")

mean_cor <- rowMeans(abs(heat_df), na.rm = TRUE)

heat_order_df <- data.frame(
  Gene = rownames(heat_df),
  mean_abs_cor = mean_cor,
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(gene_category, by = "Gene") %>%
  dplyr::mutate(
    Category = ifelse(is.na(Category) | Category == "", "Others", Category),
    Category = factor(Category, levels = category_order)
  ) %>%
  dplyr::arrange(Category, dplyr::desc(mean_abs_cor))

# Map SuperCategory
super_cat_map <- c(
  "Synaptic Structure" = "Synaptic Function",
  "Synaptic Transmission" = "Synaptic Function",
  "Neuronal Differentiation" = "Synaptic Function",
  "Epigenetic Regulation" = "Cellular Function",
  "Cell Signaling" = "Cellular Function",
  "Mitochondrial Function" = "Cellular Function",
  "Membrane Function" = "Cellular Function",
  "Protein Process" = "Cellular Function",
  "Metabolic Process" = "Cellular Function",
  "Others" = "Others"
)
heat_order_df$SuperCategory <- super_cat_map[as.character(heat_order_df$Category)]

ordered_genes <- heat_order_df$Gene
cor_mat <- as.matrix(heat_df[ordered_genes, , drop = FALSE])
cor_mat_t <- t(cor_mat)
rownames(cor_mat_t) <- gsub("Correlation_", "", rownames(cor_mat_t))

# Cohen's d matrix
cohen_df <- final_df %>%
  dplyr::select(Gene, dplyr::starts_with("Cohen.d_")) %>%
  tibble::column_to_rownames("Gene")
cohen_mat <- as.matrix(cohen_df[ordered_genes, , drop = FALSE])
cohen_mat_t <- t(cohen_mat)
rownames(cohen_mat_t) <- gsub("Cohen.d_", "", rownames(cohen_mat_t))

# 1. Export SourceData
cat("Exporting SourceData for Figure 2D...\n")
sd_cor <- data.frame(
  Gene = heat_order_df$Gene,
  SuperCategory = heat_order_df$SuperCategory,
  Category = heat_order_df$Category,
  Mean_Abs_Correlation = round(heat_order_df$mean_abs_cor, 4),
  Correlation_Amygdala = round(final_df$Correlation_Amygdala[match(ordered_genes, final_df$Gene)], 4),
  Correlation_Cortex = round(final_df$Correlation_Cortex[match(ordered_genes, final_df$Gene)], 4),
  Correlation_Hippocampus = round(final_df$Correlation_Hippocampus[match(ordered_genes, final_df$Gene)], 4)
)
write.csv(sd_cor, file.path(out_dir_sd, "Fig2D_Candidate_Gene_Correlation_Adenosine.csv"), row.names = FALSE)

sd_cohen <- data.frame(
  Gene = heat_order_df$Gene,
  SuperCategory = heat_order_df$SuperCategory,
  Category = heat_order_df$Category,
  Cohen.d_Amygdala = round(final_df$Cohen.d_Amygdala[match(ordered_genes, final_df$Gene)], 4),
  Cohen.d_Cortex = round(final_df$Cohen.d_Cortex[match(ordered_genes, final_df$Gene)], 4),
  Cohen.d_Hippocampus = round(final_df$Cohen.d_Hippocampus[match(ordered_genes, final_df$Gene)], 4)
)
write.csv(sd_cohen, file.path(out_dir_sd, "Fig2D_Candidate_Gene_CohenD_Epi_vs_Ctrl.csv"), row.names = FALSE)

# 2. Generate Check Plots
cat("Generating Check Plot for Figure 2D...\n")
annotation_col <- heat_order_df %>%
  dplyr::select(Gene, Category) %>%
  tibble::column_to_rownames("Gene")

category_colors <- RColorBrewer::brewer.pal(n = length(category_order), name = "Set3")
names(category_colors) <- category_order
ann_colors <- list(Category = category_colors)

color_cor <- colorRampPalette(c("#2166AC", "white"))(100)
color_cohen <- colorRampPalette(c("white", "#B2182B"))(100)

p_cor <- pheatmap(
  cor_mat_t,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = color_cor,
  cellwidth = 11, cellheight = 11,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = TRUE, show_colnames = FALSE,
  fontsize_row = 8, fontsize_col = 8,
  silent = TRUE
)

p_cohen <- pheatmap(
  cohen_mat_t,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = color_cohen,
  cellwidth = 11, cellheight = 11,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = TRUE, show_colnames = TRUE,
  fontsize_row = 8, fontsize_col = 8, angle_col = 45,
  silent = TRUE
)

png(file.path(out_dir_plot, "check_Fig2D_Candidate_Gene_Heatmap.png"), width = 3600, height = 1800, res = 300)
grid.arrange(p_cor$gtable, p_cohen$gtable, ncol = 1)
dev.off()

cat("Figure 2D finished successfully.\n")
