# ==============================================================================
# Script Name: Fig7J_Ado_Pathways_Expression_Violin.R
# Description: Generates Figure 7J: Expression violins of core genes across 
#              Purine Metabolism, Adenosine Metabolism, and Adenosine Transport pathways.
# Author: Auto-refactored for EpiAdo STM Revision
# Date: 2026-09-08
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(rstatix)
})

# Define paths
dbs_seurat_path <- "Z:/2023_EpiAdo/data/snSeq_JSH/processed/merged_ps_log_pca_hrm_final_seuratobj"
dbs_cpm_path <- "Z:/2023_EpiAdo/data/snSeq_JSH/processed/merged_ps_log_pca_hrm_final_logcpm"
out_dir_data <- "../sourcedata/Figure7"
out_dir_plot <- "../sourcedata/check_plots"

dir.create(out_dir_data, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_plot, recursive = TRUE, showWarnings = FALSE)

# 1. Load data
obj <- qread(dbs_seurat_path)
obj_vglut2 <- subset(obj, subset = celltype_coarse == "VGLUT2+ Glu.N")
cpm <- qread(dbs_cpm_path)

cells_vglut2 <- Cells(obj_vglut2)
cpm_vglut2 <- cpm[, colnames(cpm) %in% cells_vglut2]

# Define gene sets matching paper Figure 7J
genes_purine <- c("Ppat", "Paics", "Pfas", "Gart", "Adsl", "Adss1", "Adss2", "Ampd1", "Ampd2", "Ampd3", "Gmpr", "Gmpr2", "Impdh1", "Impdh2", "Hprt", "Aprt")
genes_ado_metab <- c("Nt5e", "Nt5c2", "Nt5c1a", "Adk", "Ada", "Enpp1", "Entpd1")
genes_ado_trans <- c("Slc29a1", "Slc29a2", "Slc28a2", "Slc28a3")

genes_purine <- intersect(genes_purine, rownames(cpm_vglut2))
genes_ado_metab <- intersect(genes_ado_metab, rownames(cpm_vglut2))
genes_ado_trans <- intersect(genes_ado_trans, rownames(cpm_vglut2))

# Compute mean expression of each gene in Sham (Ctrl_2), DBS Contra. (Ctrl_1), and DBS (DBS_1)
calc_pathway_genes <- function(gene_list, pathway_name) {
  res_list <- list()
  for (gene in gene_list) {
    exp_vec <- cpm_vglut2[gene, ]
    for (smp in c("Ctrl_2", "Ctrl_1", "DBS_1")) {
      grp_cells <- Cells(obj_vglut2)[obj_vglut2$Sample == smp]
      vals <- exp_vec[grp_cells]
      mean_val <- mean(vals[vals > 0], na.rm = TRUE)
      if (is.nan(mean_val)) mean_val <- 0
      
      grp_label <- case_when(
        smp == "Ctrl_2" ~ "Sham",
        smp == "Ctrl_1" ~ "DBS Contra.",
        smp == "DBS_1" ~ "DBS"
      )
      
      res_list[[length(res_list) + 1]] <- data.frame(
        Pathway = pathway_name,
        Gene = gene,
        Sample = smp,
        Group = grp_label,
        Mean_Expression = round(mean_val, 4),
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows(res_list)
}

df_purine <- calc_pathway_genes(genes_purine, "Purine Metabolism")
df_metab <- calc_pathway_genes(genes_ado_metab, "Adenosine Metabolism")
df_trans <- calc_pathway_genes(genes_ado_trans, "Adenosine Transport")

all_pathways_df <- bind_rows(df_purine, df_metab, df_trans)
all_pathways_df$Group <- factor(all_pathways_df$Group, levels = c("Sham", "DBS Contra.", "DBS"))

csv_file_data <- file.path(out_dir_data, "Fig7J_Ado_Pathways_Gene_Expression.csv")
write.csv(all_pathways_df, csv_file_data, row.names = FALSE)
cat("SourceData exported to:", csv_file_data, "\n")

# Compute t-test statistics
pathway_names <- c("Purine Metabolism", "Adenosine Metabolism", "Adenosine Transport")
stat_list <- list()
plot_list <- list()

color_palette <- c("Sham" = "#989798", "DBS Contra." = "#a1d99b", "DBS" = "#31a354")

for (p_name in pathway_names) {
  sub_df <- all_pathways_df %>% filter(Pathway == p_name)
  
  t_res <- sub_df %>% t_test(Mean_Expression ~ Group, paired = FALSE)
  
  for (r in 1:nrow(t_res)) {
    stat_list[[length(stat_list) + 1]] <- data.frame(
      Pathway = p_name,
      Comparison = paste(t_res$group1[r], "vs", t_res$group2[r]),
      p_value = signif(t_res$p[r], 3),
      stringsAsFactors = FALSE
    )
  }
  
  set.seed(42)
  sub_df$Jitter_X <- as.numeric(sub_df$Group) + 0.18 + runif(nrow(sub_df), -0.06, 0.06)
  
  max_y <- max(sub_df$Mean_Expression)
  min_y <- min(sub_df$Mean_Expression)
  y_rng <- max_y - min_y
  
  # Format p labels from paper
  p_sham_dbs <- t_res$p[t_res$group1 == "Sham" & t_res$group2 == "DBS"]
  p_sham_contra <- t_res$p[t_res$group1 == "Sham" & t_res$group2 == "DBS Contra."]
  p_contra_dbs <- t_res$p[t_res$group1 == "DBS Contra." & t_res$group2 == "DBS"]
  
  lbl1 <- if (p_sham_contra < 0.001) "P < 10^-3" else paste0("P = ", round(p_sham_contra, 2))
  lbl2 <- if (p_contra_dbs < 0.001) "P < 10^-3" else paste0("P = ", round(p_contra_dbs, 3))
  lbl3 <- if (p_sham_dbs < 0.001) "P < 10^-3" else paste0("P = ", round(p_sham_dbs, 3))
  
  p <- ggplot(sub_df, aes(x = Group, y = Mean_Expression)) +
    geom_violin(aes(fill = Group, color = Group), width = 0.55, alpha = 0.35, position = position_nudge(x = -0.15)) +
    geom_boxplot(aes(color = Group), width = 0.18, fill = "white", outlier.shape = NA, linewidth = 0.5, position = position_nudge(x = -0.15)) +
    geom_point(aes(x = Jitter_X, y = Mean_Expression, fill = Group), shape = 21, color = "black", size = 1.8, stroke = 0.3, alpha = 0.85) +
    annotate("segment", x = 1, xend = 2, y = max_y + y_rng * 0.08, yend = max_y + y_rng * 0.08, linewidth = 0.35) +
    annotate("text", x = 1.5, y = max_y + y_rng * 0.13, label = lbl1, size = 2.6) +
    annotate("segment", x = 2, xend = 3, y = max_y + y_rng * 0.22, yend = max_y + y_rng * 0.22, linewidth = 0.35) +
    annotate("text", x = 2.5, y = max_y + y_rng * 0.27, label = lbl2, size = 2.6) +
    annotate("segment", x = 1, xend = 3, y = max_y + y_rng * 0.36, yend = max_y + y_rng * 0.36, linewidth = 0.35) +
    annotate("text", x = 2.0, y = max_y + y_rng * 0.41, label = lbl3, size = 2.6) +
    scale_fill_manual(values = color_palette) +
    scale_color_manual(values = color_palette) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    labs(title = p_name, x = NULL, y = if (p_name == "Purine Metabolism") "Expression of Genes in Pathway\n(LogCPM)" else NULL) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 9.5, face = "bold", hjust = 0.5, color = "black"),
      axis.title = element_text(size = 9, color = "black"),
      axis.text.x = element_text(size = 8.5, color = "black"),
      axis.text.y = element_text(size = 8, color = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.5),
      axis.ticks = element_line(colour = "black", linewidth = 0.5),
      legend.position = "none"
    )
  
  plot_list[[p_name]] <- p
}

stat_df <- bind_rows(stat_list)
csv_file_stat <- file.path(out_dir_data, "Fig7J_Ado_Pathways_Statistics.csv")
write.csv(stat_df, csv_file_stat, row.names = FALSE)
cat("Statistics SourceData exported to:", csv_file_stat, "\n")
print(stat_df)

final_p7j <- plot_grid(plotlist = plot_list, nrow = 1, align = "h")
check_plot_file <- file.path(out_dir_plot, "check_Fig7J_Ado_Pathways_Expression_Violin.png")
ggsave(check_plot_file, final_p7j, width = 7.5, height = 3.6, dpi = 300)
cat("Check plot saved to:", check_plot_file, "\n")