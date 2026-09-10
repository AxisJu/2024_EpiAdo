# ==============================================================================
# Script Name: FigS2B_ADORA1_Control_FCD_TLE_Violin.R
# Description: Violin plots showing ADORA1 expression in Glu.N across Control,
#              FCD, and TLE samples (Supplementary Figure 2B).
# ==============================================================================

library(dplyr)
library(ggplot2)
library(gghalves)
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

glun_cells <- rownames(meta_df)[meta_df$celltype_coarse == "Glu.N"]
glun_cells <- intersect(glun_cells, colnames(sub_cpm))

# Author filtering from 260205_revSTM.R: exclude outlier samples
meta_glun <- meta_df[glun_cells, ] %>%
  dplyr::filter(!Sample_snSeq.processed %in% c("CTX_Epi_002", "AMY_Epi_002")) %>%
  dplyr::mutate(
    Pathology_Group = dplyr::case_when(
      grepl("Control", Disease) ~ "Control",
      Disease %in% c("Epilepsy, TLE", "Epilepsy, MTLE", "Epilepsy, LEAT") ~ "TLE",
      Disease %in% c("Epilepsy, FCD2a", "Epilepsy, FCD2b") ~ "FCD",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(Pathology_Group))

meta_glun$Pathology_Group <- factor(meta_glun$Pathology_Group, levels = c("Control", "FCD", "TLE"))

tp_data <- data.frame(
  CellID = rownames(meta_glun),
  Pathology_Group = meta_glun$Pathology_Group,
  Sample = meta_glun$Sample_snSeq.processed,
  Exp = as.numeric(sub_cpm["ADORA1", rownames(meta_glun)])
) %>%
  dplyr::filter(Exp > 0)

# Statistical tests: FCD vs Control, TLE vs Control
comparisons <- list(c("FCD", "Control"), c("TLE", "Control"))
stats_list <- list()

for (comp in comparisons) {
  sub_d <- tp_data %>% dplyr::filter(Pathology_Group %in% comp)
  w_res <- wilcox.test(Exp ~ Pathology_Group, data = sub_d)
  d_res <- effsize::cohen.d(sub_d$Exp ~ sub_d$Pathology_Group)$estimate
  
  p_val <- w_res$p.value
  if (p_val < 1e-100) {
    p_lbl <- "P < 1 x 10^-100"
  } else {
    p_lbl <- sprintf("P = %.2e", p_val)
  }
  
  stats_list[[paste(comp, collapse = "_vs_")]] <- data.frame(
    Comparison = paste(comp[1], "vs", comp[2]),
    Control_Cells = sum(sub_d$Pathology_Group == "Control"),
    Disease_Cells = sum(sub_d$Pathology_Group == comp[1]),
    PValue = p_val,
    CohenD = d_res,
    Label = p_lbl
  )
}
stats_df <- do.call(rbind, stats_list)

# Export SourceData
cat("Exporting SourceData...\n")
write.csv(
  tp_data,
  file = file.path(out_data_dir, "FigS2B_ADORA1_Pathology_Expression.csv"),
  row.names = FALSE
)

write.csv(
  stats_df,
  file = file.path(out_data_dir, "FigS2B_ADORA1_Pathology_Statistics.csv"),
  row.names = FALSE
)

cat("Plotting Check Figures...\n")
# Palette from 260205_revSTM.R: Control="#fac5b1", FCD="#7fb7be", TLE="#9fa1cb"
color_pathology <- c("Control" = "#fac5b1", "FCD" = "#7fb7be", "TLE" = "#9fa1cb")

max_exp <- max(tp_data$Exp, na.rm = TRUE)

p_violin <- ggplot(tp_data, aes(x = Pathology_Group, y = Exp, fill = Pathology_Group, color = Pathology_Group)) +
  geom_half_violin(
    linewidth = 0.5, alpha = 0.3,
    position = position_dodge(width = 0.15), width = 0.75, side = "l"
  ) +
  geom_jitter(shape = 16, size = 0.01, alpha = 0.15, width = 0.15) +
  # Brackets for comparisons
  # FCD vs Control (x=1 to x=2)
  annotate("segment", x = 1, xend = 2, y = max_exp * 1.06, yend = max_exp * 1.06, color = "black", linewidth = 0.3) +
  annotate("segment", x = 1, xend = 1, y = max_exp * 1.03, yend = max_exp * 1.06, color = "black", linewidth = 0.3) +
  annotate("segment", x = 2, xend = 2, y = max_exp * 1.03, yend = max_exp * 1.06, color = "black", linewidth = 0.3) +
  annotate("text", x = 1.5, y = max_exp * 1.10, label = stats_df$Label[1], size = 3, color = "black", fontface = "italic") +
  # TLE vs Control (x=1 to x=3)
  annotate("segment", x = 1, xend = 3, y = max_exp * 1.18, yend = max_exp * 1.18, color = "black", linewidth = 0.3) +
  annotate("segment", x = 1, xend = 1, y = max_exp * 1.15, yend = max_exp * 1.18, color = "black", linewidth = 0.3) +
  annotate("segment", x = 3, xend = 3, y = max_exp * 1.15, yend = max_exp * 1.18, color = "black", linewidth = 0.3) +
  annotate("text", x = 2.2, y = max_exp * 1.22, label = stats_df$Label[2], size = 3, color = "black", fontface = "italic") +
  scale_fill_manual(values = color_pathology) +
  scale_color_manual(values = color_pathology) +
  theme_classic() +
  theme(
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = "cm"),
    title = element_text(colour = "black", size = 12, face = "bold"),
    axis.text.x = element_text(colour = "black", size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 10),
    axis.title.y = element_text(size = 11),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  labs(
    title = "ADORA1",
    y = "log2CPM"
  ) +
  ylim(0, max_exp * 1.28)

ggsave(
  filename = file.path(out_check_dir, "check_FigS2B_ADORA1_Pathology_Violin.png"),
  plot = p_violin,
  width = 4.5,
  height = 5,
  dpi = 300
)

cat("FigS2B completed successfully.\n")
