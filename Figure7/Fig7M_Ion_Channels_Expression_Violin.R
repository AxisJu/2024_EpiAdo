# ==============================================================================
# Script Name: Fig7M_Ion_Channels_Expression_Violin.R
# Description: Generates Figure 7M: Expression of 6 key ion channels altered 
#              in ADORA1-low Glu.N (Scn1b, Kcna2, Kcnh2, Kcnc1, Cacna2d2, Clcn2) 
#              in mouse ANT across Sham, DBS Contra., and DBS groups.
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
sample_vec <- obj_vglut2$Sample

channels <- c("Scn1b", "Kcna2", "Kcnh2", "Kcnc1", "Cacna2d2", "Clcn2")

data_list <- list()
stat_list <- list()
plot_list <- list()

color_palette <- c("Sham" = "#989798", "DBS Contra." = "#a1d99b", "DBS" = "#31a354")

for (gene in channels) {
  exp_vec <- cpm[gene, cells_vglut2]
  df_gene <- data.frame(
    Cell_ID = cells_vglut2,
    Gene = gene,
    Sample = sample_vec,
    Exp = exp_vec,
    stringsAsFactors = FALSE
  ) %>%
    filter(Exp > 0)
  
  df_gene$Group <- factor(
    case_when(
      df_gene$Sample == "Ctrl_2" ~ "Sham",
      df_gene$Sample == "Ctrl_1" ~ "DBS Contra.",
      df_gene$Sample == "DBS_1" ~ "DBS"
    ),
    levels = c("Sham", "DBS Contra.", "DBS")
  )
  
  data_list[[gene]] <- df_gene
  
  # Holm adjusted t-test matching original script
  t_res <- df_gene %>% t_test(Exp ~ Group, p.adjust.method = "holm")
  
  for (r in 1:nrow(t_res)) {
    stat_list[[length(stat_list) + 1]] <- data.frame(
      Gene = gene,
      Comparison = paste(t_res$group1[r], "vs", t_res$group2[r]),
      p_unadj = signif(t_res$p[r], 3),
      p_adj = signif(t_res$p.adj[r], 3),
      Test = "Holm-adjusted t-test",
      stringsAsFactors = FALSE
    )
  }
  
  # Format p-labels matching paper Figure 7M
  p1 <- t_res$p.adj[t_res$group1 == "Sham" & t_res$group2 == "DBS Contra."]
  p2 <- t_res$p.adj[t_res$group1 == "DBS Contra." & t_res$group2 == "DBS"]
  
  format_p <- function(p) {
    if (p < 1e-10) return("P < 1 %*% 10^-10")
    if (p < 1e-6) return("P < 1 %*% 10^-6")
    if (p < 1e-4) return("P < 1 %*% 10^-4")
    if (p < 0.05) return(paste0("P = ", signif(p, 2)))
    return(paste0("P = ", round(p, 2)))
  }
  
  lbl1 <- format_p(p1)
  lbl2 <- format_p(p2)
  
  set.seed(42)
  df_gene$Jitter_X <- as.numeric(df_gene$Group) + 0.18 + runif(nrow(df_gene), -0.06, 0.06)
  
  max_y <- max(df_gene$Exp)
  y_rng <- max_y - min(df_gene$Exp)
  
  p <- ggplot(df_gene, aes(x = Group, y = Exp)) +
    geom_violin(aes(fill = Group, color = Group), width = 0.55, alpha = 0.35, position = position_nudge(x = -0.15)) +
    geom_boxplot(aes(color = Group), width = 0.18, fill = "white", outlier.shape = NA, linewidth = 0.5, position = position_nudge(x = -0.15)) +
    geom_point(aes(x = Jitter_X, y = Exp, fill = Group), shape = 21, color = "black", size = 1.4, stroke = 0.25, alpha = 0.75) +
    annotate("segment", x = 1, xend = 2, y = max_y + y_rng * 0.08, yend = max_y + y_rng * 0.08, linewidth = 0.35) +
    annotate("text", x = 1.5, y = max_y + y_rng * 0.14, label = lbl1, parse = TRUE, size = 2.4) +
    annotate("segment", x = 2, xend = 3, y = max_y + y_rng * 0.24, yend = max_y + y_rng * 0.24, linewidth = 0.35) +
    annotate("text", x = 2.5, y = max_y + y_rng * 0.30, label = lbl2, parse = TRUE, size = 2.4) +
    scale_fill_manual(values = color_palette) +
    scale_color_manual(values = color_palette) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
    labs(title = gene, x = NULL, y = if (gene %in% c("Scn1b", "Kcnc1")) "Expression (LogCPM)" else NULL) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 9.5, face = "bold.italic", hjust = 0.5, color = "black"),
      axis.title = element_text(size = 8.5, color = "black"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 8, color = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.5),
      axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
      legend.position = "none"
    )
  
  plot_list[[gene]] <- p
}

# Export SourceData
all_channels_df <- bind_rows(data_list)
csv_file_data <- file.path(out_dir_data, "Fig7M_Ion_Channels_SingleCell_Expression.csv")
write.csv(all_channels_df, csv_file_data, row.names = FALSE)
cat("Cell-level SourceData exported to:", csv_file_data, "\n")

all_stats_df <- bind_rows(stat_list)
csv_file_stat <- file.path(out_dir_data, "Fig7M_Ion_Channels_Statistics.csv")
write.csv(all_stats_df, csv_file_stat, row.names = FALSE)
cat("Statistics SourceData exported to:", csv_file_stat, "\n")
print(all_stats_df)

# Assemble 2x3 grid
top_title <- ggdraw() + 
  draw_label(expression(bold("Expression of ion channels altered in ") * bolditalic("ADORA1")^bold("Low") * bold(" Glu.N")^bold("Epi")), 
             size = 10.5, hjust = 0.5)

legend_plot <- ggplot(data.frame(x = c(1, 2, 3), y = c(1, 1, 1), Group = c("Sham", "DBS Contra.", "DBS")), 
                      aes(x = x, y = y, fill = Group)) +
  geom_point(shape = 22, size = 3.5) +
  scale_fill_manual(values = color_palette) +
  theme_void() +
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 8.5))

grid_plots <- plot_grid(plotlist = plot_list, ncol = 3, align = "hv")
leg_box <- get_legend(legend_plot)

final_p7m <- plot_grid(
  top_title,
  grid_plots,
  leg_box,
  ncol = 1,
  rel_heights = c(0.08, 1, 0.08)
)

check_plot_file <- file.path(out_dir_plot, "check_Fig7M_Ion_Channels_Expression_Violin.png")
ggsave(check_plot_file, final_p7m, width = 7.8, height = 5.2, dpi = 300)
cat("Check plot saved to:", check_plot_file, "\n")