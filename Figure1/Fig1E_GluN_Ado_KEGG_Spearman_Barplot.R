# ==============================================================================
# Script: Fig1E_GluN_Ado_KEGG_Spearman_Barplot.R
# Purpose: Mean Spearman correlation between adenosine metabolic genes and 85 KEGG metabolic pathways in Glu.N
# Article: EpiAdo (Science Translational Medicine) - Figure 1E
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(ggplot2)
  library(viridis)
  library(dplyr)
})

out_dir_sd <- "../sourcedata/Figure1"
out_dir_plot <- "../sourcedata/check_plots"
dir.create(out_dir_sd, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_plot, showWarnings = FALSE, recursive = TRUE)

cat("Loading data for Figure 1E...\n")
# [DATA REQUIRED] The line below originally loaded 'merged_ps_log_pca_harmony_umap_final_seuratobj' from a private path.
# Please provide 'merged_seuratobj' by loading the appropriate data before running this script.
# Example: merged_seuratobj <- qread('path/to/merged_ps_log_pca_harmony_umap_final_seuratobj')
# [DATA REQUIRED] The line below originally loaded 'merged_ps_log_pca_harmony_umap_final_logcpm' from a private path.
# Please provide 'cpm' by loading the appropriate data before running this script.
# Example: cpm <- qread('path/to/merged_ps_log_pca_harmony_umap_final_logcpm')
# [DATA REQUIRED] The line below originally loaded 'signature_exp' from a private path.
# Please provide 'signature_exp' by loading the appropriate data before running this script.
# Example: signature_exp <- qread('path/to/signature_exp')

gene.list <- c(
  "ENTPD1","ENTPD2","ENTPD3","ENTPD8","ENPP1","ENPP3","NT5E","ALPL",
  "CD38","BST1","AK1","NME1","NME2","NT5C1A","NT5C2", 
  "SLC29A1","SLC29A2","SLC29A3","SLC29A4","SLC28A1","SLC28A2","SLC28A3",
  "ADK","ADA","PNP","AHCY","DNMT1","DNMT3A","DNMT3B","DNMT3L"
)

cells <- Cells(merged_seuratobj)[merged_seuratobj$celltype_coarse == "Glu.N"]
cpm.cells <- cpm[, cells]
sig_cells <- gsub("-", "\\.", cells)
sig_cells <- intersect(sig_cells, rownames(signature_exp))
valid_cells_original <- gsub("\\.", "-", sig_cells)
cpm.cells <- cpm.cells[, valid_cells_original]

valid_genes <- intersect(gene.list, rownames(cpm.cells))
mat_gene <- t(as.matrix(cpm.cells[valid_genes, ]))
mat_ptw <- as.matrix(signature_exp[sig_cells, ])

cat("Calculating Spearman correlation matrix...\n")
rho_mat <- cor(mat_ptw, mat_gene, method = "spearman")

n <- nrow(mat_ptw)
t_stat <- rho_mat * sqrt((n - 2) / (1 - rho_mat^2))
p_mat <- 2 * pt(abs(t_stat), df = n - 2, lower.tail = FALSE)

cor.mat <- data.frame(
  Pathway = rownames(rho_mat),
  Celltype = "Glu.N",
  SigNum = rowSums(p_mat < 0.01, na.rm = TRUE),
  SigPercent = rowMeans(p_mat < 0.01, na.rm = TRUE),
  Mean_Rho = rowMeans(rho_mat, na.rm = TRUE),
  Median_Rho = apply(rho_mat, 1, median, na.rm = TRUE),
  Mean_Abs = rowMeans(abs(rho_mat), na.rm = TRUE),
  Median_Abs = apply(abs(rho_mat), 1, median, na.rm = TRUE),
  stringsAsFactors = FALSE
)

cor.mat <- cor.mat[order(cor.mat$Mean_Rho, decreasing = TRUE), ]
cor.mat$Pathway <- factor(cor.mat$Pathway, levels = cor.mat$Pathway)
cor.mat$Category <- ifelse(cor.mat$Pathway == "Purine metabolism", "Purine metabolism", "Others")

# Save SourceData
sd_file <- file.path(out_dir_sd, "Fig1E_GluN_Ado_KEGG_Spearman_Correlation.csv")
write.csv(cor.mat, sd_file, row.names = FALSE)
cat("Saved SourceData to:", sd_file, "\n")

# Plot Check Figure
p <- ggplot(cor.mat, aes(x = Pathway, y = Mean_Rho)) +
  geom_col(aes(fill = Category), width = 0.35) +
  geom_point(aes(size = Mean_Rho, color = Mean_Rho), shape = 16) +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title = element_text(colour = "black", size = 11, family = "sans"),
    axis.text.x = element_blank(), 
    axis.text.y = element_text(colour = "black", size = 11, family = "sans"),
    axis.line.x = element_blank(),
    axis.line.y = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
    legend.position = "right"
  ) +
  scale_size_continuous(range = c(2, 5)) +
  scale_color_gradientn(colours = viridis(50)) +
  scale_fill_manual(values = c("Purine metabolism" = "#d66f70", "Others" = "#CCCCCC")) +
  labs(x = NULL, y = "Mean Correlation (rho)", title = "Metabolism dysfunction correlated with adenosine in Glu.N")

plot_file <- file.path(out_dir_plot, "check_Fig1E_Spearman_Barplot.png")
ggsave(plot_file, plot = p, width = 7, height = 4.5, dpi = 300)
cat("Figure 1E finished successfully.\n")
