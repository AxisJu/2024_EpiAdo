# ==============================================================================
# Script: Fig1D_Adenosine_Metabolism_Gene_Heatmap.R
# Purpose: Pseudo-bulked gene expression levels and Cohen's d effect sizes for key adenosine genes
# Article: EpiAdo (Science Translational Medicine) - Figure 1D
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(dplyr)
  library(Seurat)
  library(pheatmap)
  library(RColorBrewer)
})

out_dir_sd <- "../sourcedata/Figure1"
out_dir_plot <- "../sourcedata/check_plots"
dir.create(out_dir_sd, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir_plot, showWarnings = FALSE, recursive = TRUE)

cat("Loading data for Figure 1D...\n")
# [DATA REQUIRED] The line below originally loaded 'merged_ps_log_pca_harmony_umap_final_seuratobj' from a private path.
# Please provide 'merged_seuratobj' by loading the appropriate data before running this script.
# Example: merged_seuratobj <- qread('path/to/merged_ps_log_pca_harmony_umap_final_seuratobj')
# [DATA REQUIRED] The line below originally loaded 'merged_ps_log_pca_harmony_umap_final_logcpm' from a private path.
# Please provide 'cpm' by loading the appropriate data before running this script.
# Example: cpm <- qread('path/to/merged_ps_log_pca_harmony_umap_final_logcpm')

ado_genes_ordered <- c(
  "ENTPD1", "ENPP1", "NT5E",
  "NT5C2",
  "SLC29A2", "SLC29A1", "SLC28A3", "SLC28A2",
  "ADK", "ADA",
  "ADORA1", "ADORA2A"
)

gene_categories <- c(
  "ENTPD1" = "Extracellular Production", "ENPP1" = "Extracellular Production", "NT5E" = "Extracellular Production",
  "NT5C2" = "Intracellular Production",
  "SLC29A2" = "Adenosine Transporter", "SLC29A1" = "Adenosine Transporter",
  "SLC28A3" = "Adenosine Transporter", "SLC28A2" = "Adenosine Transporter",
  "ADK" = "Intracellular Degradation", "ADA" = "Intracellular Degradation",
  "ADORA1" = "Adenosine Receptor", "ADORA2A" = "Adenosine Receptor"
)

# Calculate pseudobulk expression and Cohen's d in Glu.N
cells_glun <- Cells(merged_seuratobj)[merged_seuratobj$celltype_coarse == "Glu.N"]
meta_glun <- merged_seuratobj@meta.data[cells_glun, ]

res_df <- data.frame()
for (g in ado_genes_ordered) {
  if (!g %in% rownames(cpm)) next
  exp_vals <- cpm[g, cells_glun]
  df_g <- data.frame(Sample = meta_glun$Sample_snSeq.processed, Group = meta_glun$Group, Exp = exp_vals)
  pb <- df_g %>%
    dplyr::group_by(Sample, Group) %>%
    dplyr::summarise(PBExp = mean(Exp), .groups = "drop")
  
  pb_ctrl <- pb$PBExp[pb$Group == "Control"]
  pb_epi  <- pb$PBExp[pb$Group == "Epilepsy"]
  
  m_c <- mean(pb_ctrl)
  m_e <- mean(pb_epi)
  s_p <- sqrt((var(pb_ctrl) + var(pb_epi)) / 2)
  cohen_d <- ifelse(s_p > 0, (m_e - m_c) / s_p, 0)
  p_val <- wilcox.test(pb_epi, pb_ctrl)$p.value
  
  res_df <- rbind(res_df, data.frame(
    Gene = g,
    FunctionalCategory = gene_categories[g],
    logCPM_Control = round(m_c, 3),
    logCPM_Epilepsy = round(m_e, 3),
    logCPM_Epi_vs_Ctrl = round(m_e - m_c, 3),
    CohenD = round(cohen_d, 3),
    PValue = p_val
  ))
}

# Export SourceData
write.csv(res_df, file.path(out_dir_sd, "Fig1D_Adenosine_Metabolism_Gene_Heatmap.csv"), row.names = FALSE)
cat("SourceData saved for Figure 1D.\n")

# Heatmap matrices
mat_exp <- as.matrix(res_df[, c("logCPM_Control", "logCPM_Epilepsy")])
rownames(mat_exp) <- res_df$Gene
colnames(mat_exp) <- c("Ctrl.", "Epi.")

mat_d <- as.matrix(res_df[, c("CohenD"), drop = FALSE])
rownames(mat_d) <- res_df$Gene
colnames(mat_d) <- c("Cohen's D")

# Plot Check Heatmap
plot_file <- file.path(out_dir_plot, "check_Fig1D_Adenosine_Metabolism_Gene_Heatmap.png")
png(plot_file, width = 1800, height = 2200, res = 300)
pheatmap(
  mat_exp,
  cellwidth = 25, cellheight = 14,
  cluster_rows = FALSE, cluster_cols = FALSE,
  color = colorRampPalette(brewer.pal(9, "Purples"))(50),
  main = "Adenosine Pathway Key Genes logCPM"
)
dev.off()
cat("Figure 1D finished successfully.\n")
