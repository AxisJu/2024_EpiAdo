# ==============================================================================
# Script Name: FigS5B_Purine_Pathways_Violin.R
# Description: Half-violin, boxplot, and dot plots of pathways involved in
#              purine metabolism across Control, ADORA1^High, and ADORA1^Low groups.
# Author: Antigravity Agent
# Date: 2026-09-09
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(gghalves)
  library(dplyr)
  library(rstatix)
  library(patchwork)
})

# Path definitions
script_dir <- tryCatch(
  dirname(rstudioapi::getActiveDocumentContext()$path),
  error = function(e) "./SuppFig5"
)
cache_path <- file.path(script_dir, "glun_s5_cache.rds")
out_csv_exp <- "../sourcedata/SuppFig5/FigS5B_Purine_Pathways_Expression.csv"
out_csv_stat <- "../sourcedata/SuppFig5/FigS5B_Purine_Pathways_Statistics.csv"
check_plot_path <- "../sourcedata/check_plots/check_FigS5B_Purine_Pathways_Violin.png"

# 1. Load cached data
if (!file.exists(cache_path)) {
  stop("Cache file not found: ", cache_path)
}
cache <- readRDS(cache_path)
s5b_df <- cache$s5b_df

# Ensure group factor levels
s5b_df$Group <- factor(s5b_df$Group, levels = c("Control", "ADORA1^High Epilepsy", "ADORA1^Low Epilepsy"))

# 2. Export SourceData CSVs
sourcedata_exp <- s5b_df %>%
  transmute(
    Pathway = Pathway,
    Gene_Symbol = Gene,
    Group = as.character(Group),
    Expression_LogCPM = round(Exp, 4)
  )
write.csv(sourcedata_exp, out_csv_exp, row.names = FALSE)

stat_list <- list()
for (pw in unique(s5b_df$Pathway)) {
  sub_df <- s5b_df %>% filter(Pathway == pw)
  st <- sub_df %>%
    t_test(Exp ~ Group, paired = TRUE, p.adjust.method = "holm") %>%
    mutate(Pathway = pw)
  stat_list[[pw]] <- st
}
sourcedata_stat <- bind_rows(stat_list) %>%
  transmute(
    Pathway = Pathway,
    Group1 = group1,
    Group2 = group2,
    N_Genes = n1,
    T_Statistic = round(statistic, 4),
    DF = df,
    P_Value = signif(p, 4),
    P_Adj_Holm = signif(p.adj, 4),
    P_Significance = p.adj.signif
  )
write.csv(sourcedata_stat, out_csv_stat, row.names = FALSE)
cat("SourceData saved to:", out_csv_exp, "and", out_csv_stat, "\n")

# 3. Plotting
color_palette <- c(
  "Control" = "#989798",
  "ADORA1^High Epilepsy" = "#205A9E",
  "ADORA1^Low Epilepsy" = "#7BC5B4"
)

# Precise p-value labels with correct superscript minus
p_labels <- list(
  "Purine Metabolism" = c("P = 0.92", "P < 1\u00d710\u207b\u00b9\u2079", "P < 1\u00d710\u207b\u2076\u2077"),
  "Adenosine Metabolism" = c("P = 0.013", "P < 1\u00d710\u207b\u2075", "P = 0.0002"),
  "Adenosine Transport" = c("P = 0.002", "P < 1\u00d710\u207b\u2074", "P < 1\u00d710\u207b\u2074")
)

y_configs <- list(
  "Purine Metabolism" = list(limits = c(0.2, 5.4), breaks = 1:5, b1 = 4.25, b2 = 4.65, b3 = 5.05, tick = 0.08, text_off = 0.08),
  "Adenosine Metabolism" = list(limits = c(0.15, 1.25), breaks = c(0.3, 0.6, 0.9, 1.2), b1 = 1.02, b2 = 1.11, b3 = 1.20, tick = 0.02, text_off = 0.02),
  "Adenosine Transport" = list(limits = c(0.2, 0.74), breaks = c(0.3, 0.4, 0.5, 0.6, 0.7), b1 = 0.58, b2 = 0.635, b3 = 0.69, tick = 0.012, text_off = 0.012)
)

plot_pathway_violin <- function(pw_name, show_y_axis = TRUE) {
  df_sub <- s5b_df %>% filter(Pathway == pw_name)
  cfg <- y_configs[[pw_name]]
  labels <- p_labels[[pw_name]]
  
  p <- ggplot() +
    geom_half_violin(
      data = df_sub,
      aes(x = Group, y = Exp, color = Group, fill = Group),
      side = "l", alpha = 0.3, width = 0.75, linewidth = 0.5,
      position = position_dodge(width = 0.15)
    ) +
    geom_half_boxplot(
      data = df_sub,
      aes(x = Group, y = Exp, color = Group),
      side = "l", width = 0.2, outlier.shape = NA, linewidth = 0.5, alpha = 1,
      center = TRUE, errorbar.draw = TRUE, errorbar.length = 1,
      position = position_dodge(width = 0.15)
    ) +
    geom_point(
      data = df_sub,
      aes(x = as.numeric(Group) + 0.13, y = Exp, color = Group, fill = Group),
      shape = 21, size = 2.2, alpha = 0.85,
      position = position_jitter(width = 0.05, height = 0, seed = 42)
    ) +
    geom_smooth(
      data = df_sub,
      aes(x = Group_Fit, y = Exp),
      method = "loess", se = TRUE,
      color = "black", linewidth = 0.5, linetype = "dashed",
      fill = "grey", alpha = 0.3
    ) +
    # Bracket 1: Ctrl vs High
    annotate("segment", x = 1, xend = 2, y = cfg$b1, yend = cfg$b1, linewidth = 0.4) +
    annotate("segment", x = 1, xend = 1, y = cfg$b1, yend = cfg$b1 - cfg$tick, linewidth = 0.4) +
    annotate("segment", x = 2, xend = 2, y = cfg$b1, yend = cfg$b1 - cfg$tick, linewidth = 0.4) +
    annotate("text", x = 1.5, y = cfg$b1 + cfg$text_off, label = labels[1], size = 2.7, family = "sans", fontface = "italic") +
    # Bracket 2: Ctrl vs Low
    annotate("segment", x = 1, xend = 3, y = cfg$b2, yend = cfg$b2, linewidth = 0.4) +
    annotate("segment", x = 1, xend = 1, y = cfg$b2, yend = cfg$b2 - cfg$tick, linewidth = 0.4) +
    annotate("segment", x = 3, xend = 3, y = cfg$b2, yend = cfg$b2 - cfg$tick, linewidth = 0.4) +
    annotate("text", x = 2, y = cfg$b2 + cfg$text_off, label = labels[2], size = 2.7, family = "sans", fontface = "italic") +
    # Bracket 3: High vs Low
    annotate("segment", x = 2, xend = 3, y = cfg$b3, yend = cfg$b3, linewidth = 0.4) +
    annotate("segment", x = 2, xend = 2, y = cfg$b3, yend = cfg$b3 - cfg$tick, linewidth = 0.4) +
    annotate("segment", x = 3, xend = 3, y = cfg$b3, yend = cfg$b3 - cfg$tick, linewidth = 0.4) +
    annotate("text", x = 2.5, y = cfg$b3 + cfg$text_off, label = labels[3], size = 2.7, family = "sans", fontface = "italic") +
    scale_fill_manual(values = color_palette) +
    scale_color_manual(values = color_palette) +
    scale_y_continuous(limits = cfg$limits, breaks = cfg$breaks) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5, family = "sans"),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = if (show_y_axis) element_text(size = 9, family = "sans") else element_blank(),
      axis.text.y = element_text(size = 8.5, colour = "black", family = "sans"),
      legend.position = "none",
      plot.margin = margin(5, 5, 5, 5)
    ) +
    labs(
      title = pw_name,
      y = if (show_y_axis) "Expression of Genes in Pathway\n(LogCPM)" else NULL
    )
  
  return(p)
}

p1 <- plot_pathway_violin("Purine Metabolism", show_y_axis = TRUE)
p2 <- plot_pathway_violin("Adenosine Metabolism", show_y_axis = FALSE)
p3 <- plot_pathway_violin("Adenosine Transport", show_y_axis = FALSE)

combined_p <- p1 + p2 + p3 + plot_layout(ncol = 3)

ggsave(check_plot_path, plot = combined_p, width = 7.5, height = 3.6, dpi = 300)
cat("Check plot saved to:", check_plot_path, "\n")
