# ==============================================================================
# Script Name: FigS2E_Astrocyte_Adenosine_Metabolism_Heatmap.R
# Description: Heatmap of Cohen's d for adenosine metabolism genes in astrocytes
#              across 3 brain regions (Supplementary Figure 2E).
# ==============================================================================

library(dplyr)
library(tidyr)
library(pheatmap)
library(circlize)
library(effsize)

out_data_dir <- "../sourcedata/SuppFig2"
out_check_dir <- "../sourcedata/check_plots"

cache_data <- file.path(out_data_dir, "ado_genes_glun_astro_meta_exp.rds")

if (!file.exists(cache_data)) {
  stop("Please run FigS2A first to generate ado_genes_glun_astro_meta_exp.rds cache.")
}

cat("Loading targeted data from cache...\n")
save_payload <- readRDS(cache_data)
meta_df <- save_payload$meta
sub_cpm <- save_payload$cpm

astro_cells <- rownames(meta_df)[meta_df$celltype_coarse == "Astro."]
astro_cells <- intersect(astro_cells, colnames(sub_cpm))

# 12 genes displayed in Figure S2E (matching author script 260205_revSTM.R and paper panel E)
gene_list <- c(
  "ADORA1", "ADORA2A",
  "ENTPD1", "ENPP1", "NT5E", "NT5C2",
  "SLC29A1", "SLC29A2",
  "SLC28A2", "SLC28A3",
  "ADK", "ADA"
)
gene_list <- intersect(gene_list, rownames(sub_cpm))

# Region order in paper: Amygdala, Cortex, Hippocampus
regions <- c("Amygdala", "Cortex", "Hippocampus")

cat("Computing effect sizes in astrocytes across regions...\n")
# If Seurat data slot was cached or direct from obj
if (!"astro_data" %in% names(save_payload)) {
  library(qs)
  library(Seurat)
  obj_path <- "Z:/2023_EpiAdo/data/scSeq_MSH/Merged/merged_ps_log_pca_harmony_umap_final_seuratobj"
  cat("Loading Seurat Astrocyte data layer...\n")
  obj <- qread(obj_path)
  astro <- subset(obj, subset = celltype_coarse == "Astro.")
  plot_data <- FetchData(astro, vars = c(gene_list, "Group", "BrainRegion"), slot = "data")
  rm(obj, astro); gc()
  save_payload$astro_data <- plot_data
  saveRDS(save_payload, file = cache_data)
} else {
  plot_data <- save_payload$astro_data
}

stats_results <- expand.grid(Region = regions, Gene = gene_list) %>%
  group_by(Region, Gene) %>%
  do({
    region_data <- plot_data[plot_data$BrainRegion == .$Region, ]
    gene_name <- as.character(.$Gene)
    if (length(unique(region_data$Group)) == 2) {
      d_res <- effsize::cohen.d(region_data[[gene_name]] ~ region_data$Group)
      d_val <- -d_res$estimate
      w_res <- wilcox.test(region_data[[gene_name]] ~ region_data$Group)
      p_val <- w_res$p.value
      sig_label <- case_when(
        abs(d_val) > 0.6 ~ "***",
        abs(d_val) > 0.4 ~ "**",
        abs(d_val) > 0.2 ~ "*",
        TRUE ~ ""
      )
      data.frame(CohenD = round(d_val, 4), PValue = signif(p_val, 4), Signif = sig_label)
    } else {
      data.frame(CohenD = NA, PValue = NA, Signif = "")
    }
  }) %>%
  ungroup()

stats_df <- stats_results

# Pivot to matrix
d_mat <- stats_df %>%
  dplyr::select(Region, Gene, CohenD) %>%
  pivot_wider(names_from = Gene, values_from = CohenD) %>%
  tibble::column_to_rownames("Region") %>%
  as.matrix()

sig_mat <- stats_df %>%
  dplyr::select(Region, Gene, Signif) %>%
  pivot_wider(names_from = Gene, values_from = Signif) %>%
  tibble::column_to_rownames("Region") %>%
  as.matrix()

d_mat <- d_mat[regions, gene_list]
sig_mat <- sig_mat[regions, gene_list]

# Export SourceData
cat("Exporting SourceData...\n")
write.csv(
  stats_df,
  file = file.path(out_data_dir, "FigS2E_Astrocyte_Ado_Genes_Statistics.csv"),
  row.names = FALSE
)

write.csv(
  d_mat,
  file = file.path(out_data_dir, "FigS2E_Astrocyte_Ado_Genes_CohenD_Matrix.csv"),
  row.names = TRUE
)

cat("Plotting Check Figures...\n")
png(
  file = file.path(out_check_dir, "check_FigS2E_Astrocyte_Ado_Heatmap.png"),
  width = 7.5,
  height = 3.5,
  units = "in",
  res = 300
)

pheatmap(
  as.matrix(d_mat),
  display_numbers = as.matrix(sig_mat),
  color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
  breaks = seq(-1, 1, length.out = 101),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cellwidth = 24,
  cellheight = 20,
  main = "Astrocyte Adenosine Metabolism",
  fontsize_number = 14,
  number_color = "black",
  na_col = "grey90",
  border_color = "white",
  angle_col = 45
)
dev.off()

cat("FigS2E completed successfully.\n")
