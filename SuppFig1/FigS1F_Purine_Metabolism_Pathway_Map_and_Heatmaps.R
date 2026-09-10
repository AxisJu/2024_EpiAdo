# ==============================================================================
# Script Name: FigS1F_Purine_Metabolism_Pathway_Map_and_Heatmaps.R
# Description: Pathway map of purine metabolism (KEGG: HSA00230) and heatmaps
#              of logCPM values and Cohen's d effect sizes grouped by functional
#              categories (Supplementary Figure 1F).
# ==============================================================================

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(RColorBrewer)

out_data_dir <- "../sourcedata/SuppFig1"
out_check_dir <- "../sourcedata/check_plots"

kegg_img_src <- "Z:/2023_EpiAdo/results/20241210 Metabolism/hsa00230.250109.png"
if (file.exists(kegg_img_src)) {
  file.copy(kegg_img_src, file.path(out_check_dir, "check_FigS1F_Purine_KEGG_HSA00230.png"), overwrite = TRUE)
}

# Define functional categories (strictly ordered as in paper)
gene_categories <- list(
  "Purine Synthesis" = c("PRPS1", "PRPS2", "PPAT", "GART", "PFAS", "PAICS", "ADSL", "ATIC", "ADSS", "ADSSL1", "GMPS"),
  "cAMP/cGMP Signal" = c(
    "PDE1A", "PDE1B", "PDE1C", "PDE2A", "PDE3A", "PDE3B", "PDE4A", "PDE4B", "PDE4C", "PDE4D",
    "PDE5A", "PDE6A", "PDE6B", "PDE6C", "PDE6D", "PDE6G", "PDE6H", "PDE7A", "PDE7B", "PDE8A",
    "PDE8B", "PDE9A", "PDE10A", "PDE11A", "ADCY1", "ADCY10", "ADCY2", "ADCY3", "ADCY4", "ADCY5",
    "ADCY6", "ADCY7", "ADCY8", "ADCY9", "GUCY1A1", "GUCY1A2", "GUCY1B1", "GUCY2C", "GUCY2D", "GUCY2F"
  ),
  "Purine Degradation" = c("PNP", "ADA", "ADA2", "XDH", "URAD", "ALLC"),
  "Purine Nucleotide Recycle" = c("HPRT1", "APRT"),
  "Nucleotide Metabolism" = c(
    "AK1", "AK2", "AK3", "AK4", "AK5", "AK6", "AK7", "AK8", "AK9",
    "NME1", "NME2", "NME3", "NME4", "NME6", "NME7"
  )
)

stat_csv <- file.path(out_data_dir, "FigS1F_Purine_Metabolism_Gene_Statistics.csv")
cat_csv <- file.path(out_data_dir, "FigS1F_Purine_Metabolism_Categorized_Heatmap_Data.csv")

if (!file.exists(stat_csv) || !file.exists(cat_csv)) {
  library(Seurat)
  library(qs)
  library(effsize)
  
  obj_path <- "Z:/2023_EpiAdo/data/scSeq_MSH/Merged/merged_ps_log_pca_harmony_umap_final_seuratobj"
  cpm_path <- "Z:/2023_EpiAdo/data/scSeq_MSH/Merged/merged_ps_log_pca_harmony_umap_final_logcpm"
  
  merged_seuratobj <- qread(obj_path)
  cpm <- qread(cpm_path)
  
  all_target_genes <- unique(unlist(gene_categories))
  target_genes <- intersect(all_target_genes, rownames(cpm))
  target_cpm <- cpm[target_genes, , drop = FALSE]
  rm(cpm); gc()
  
  results <- list()
  for (gene in target_genes) {
    tp.data <- data.frame(
      Group = merged_seuratobj$Group,
      Celltype = merged_seuratobj$celltype_coarse,
      Cohort = merged_seuratobj$Cohort,
      Sample = merged_seuratobj$Sample_snSeq.processed,
      Exp = target_cpm[gene, ]
    ) %>%
      dplyr::filter(Exp > 0 & Celltype == "Glu.N" & !Cohort %in% c("Tran", "YC", "LWS")) %>%
      group_by(Sample, Group) %>%
      summarise(PBExp = mean(Exp), .groups = "drop")
    
    if (nrow(tp.data) < 5 || length(unique(tp.data$Group)) < 2) next
    
    ctrl_mean <- mean(tp.data$PBExp[tp.data$Group == "Control"])
    epi_mean <- mean(tp.data$PBExp[tp.data$Group == "Epilepsy"])
    
    d_val <- tryCatch({
      effsize::cohen.d(PBExp ~ Group, data = tp.data)$estimate
    }, error = function(e) NA)
    
    p_val <- tryCatch({
      wilcox.test(PBExp ~ Group, data = tp.data)$p.value
    }, error = function(e) NA)
    
    results[[gene]] <- data.frame(
      Gene = gene,
      Control_logCPM = round(ctrl_mean, 4),
      Epilepsy_logCPM = round(epi_mean, 4),
      Cohen_d = round(d_val, 4),
      Wilcoxon_P = signif(p_val, 4)
    )
  }
  
  res_df <- do.call(rbind, results)
  cat_df <- stack(gene_categories)
  colnames(cat_df) <- c("Gene", "Category")
  res_df <- left_join(res_df, cat_df, by = "Gene")
  
  write.csv(res_df, file = stat_csv, row.names = FALSE)
  write.csv(res_df %>% dplyr::filter(!is.na(Category)), file = cat_csv, row.names = FALSE)
}

cat("Reading SourceData for plotting...\n")
plot_genes_df <- read.csv(cat_csv)

# Preserve exact ordering within each category as defined in list
all_ordered_genes <- unlist(gene_categories)
plot_genes_df <- plot_genes_df %>%
  dplyr::filter(Gene %in% all_ordered_genes) %>%
  mutate(Gene = factor(Gene, levels = all_ordered_genes)) %>%
  arrange(Gene)

plot_genes_df$Category <- factor(plot_genes_df$Category, levels = names(gene_categories))

mat_cpm <- as.matrix(plot_genes_df[, c("Control_logCPM", "Epilepsy_logCPM")])
rownames(mat_cpm) <- as.character(plot_genes_df$Gene)
colnames(mat_cpm) <- c("Ctrl", "Epi")

# For Cohen's D, the paper shows downregulation effect (Epi - Ctrl) as negative values
mat_d <- as.matrix(plot_genes_df[, "Cohen_d", drop = FALSE])
rownames(mat_d) <- as.character(plot_genes_df$Gene)
colnames(mat_d) <- "Epi/Ctrl"

# Palette matching original paper (logCPM in green, Cohen's D in gradient greens for negative effect)
col_cpm <- colorRamp2(c(0.4, 0.9, 1.4), c("#E5F5E0", "#74C476", "#00441B"))
col_d <- colorRamp2(
  c(-1.0, -0.8, -0.6, -0.4, -0.2, 0.0),
  c("#00441B", "#238B45", "#41AB5D", "#74C476", "#C7E9C0", "#F7FCF5")
)

ht_cpm <- Heatmap(
  mat_cpm,
  name = "logCPM",
  col = col_cpm,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = plot_genes_df$Category,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  width = unit(2.2, "cm")
)

ht_d <- Heatmap(
  mat_d,
  name = "Cohen's D",
  col = col_d,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  width = unit(1.2, "cm")
)

png(
  file = file.path(out_check_dir, "check_FigS1F_Purine_Heatmaps.png"),
  width = 6.5,
  height = 16,
  units = "in",
  res = 300
)
draw(ht_cpm + ht_d, auto_adjust = FALSE)
dev.off()

cat("FigS1F completed successfully.\n")
