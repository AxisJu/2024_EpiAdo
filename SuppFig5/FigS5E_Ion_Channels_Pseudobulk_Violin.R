# ==============================================================================
# Script Name: FigS5E_Ion_Channels_Pseudobulk_Violin.R
# Description: Half-violin, boxplot, and dot plots of key ion channels
#              (Sodium, Potassium, Calcium, Chloride) across Control,
#              ADORA1^High, and ADORA1^Low Epilepsy at pseudobulk level.
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
out_csv_exp <- "../sourcedata/SuppFig5/FigS5E_Ion_Channels_Pseudobulk_Expression.csv"
out_csv_stat <- "../sourcedata/SuppFig5/FigS5E_Ion_Channels_Statistics.csv"
check_plot_path <- "../sourcedata/check_plots/check_FigS5E_Ion_Channels_Pseudobulk_Violin.png"

# 1. Load cached data
if (!file.exists(cache_path)) {
  stop("Cache file not found: ", cache_path)
}
cache <- readRDS(cache_path)
s5e_df <- cache$s5e_df

# Ensure factor order
s5e_df$Group <- factor(s5e_df$Group, levels = c("Control", "ADORA1^High Epilepsy", "ADORA1^Low Epilepsy"))

# 2. Export SourceData CSVs
sourcedata_exp <- s5e_df %>%
  transmute(
    Gene = Gene,
    Sample_ID = Sample,
    Group = as.character(Group),
    Expression_LogCPM = round(PB, 4)
  )
write.csv(sourcedata_exp, out_csv_exp, row.names = FALSE)

stat_list <- list()
for (g in unique(s5e_df$Gene)) {
  sub_d <- s5e_df %>% filter(Gene == g)
  st <- sub_d %>%
    t_test(PB ~ Group, p.adjust.method = "holm") %>%
    mutate(Gene = g)
  stat_list[[g]] <- st
}
sourcedata_stat <- bind_rows(stat_list) %>%
  transmute(
    Gene = Gene,
    Group1 = group1,
    Group2 = group2,
    N_Samples_Group1 = n1,
    N_Samples_Group2 = n2,
    T_Statistic = round(statistic, 4),
    DF = df,
    P_Value = signif(p, 4),
    P_Adj_Holm = signif(p.adj, 4),
    P_Significance = p.adj.signif
  )
write.csv(sourcedata_stat, out_csv_stat, row.names = FALSE)
cat("SourceData saved to:", out_csv_exp, "and", out_csv_stat, "\n")

# 3. Plotting configurations matching original publication
color_palette <- c(
  "Control" = "#989798",
  "ADORA1^High Epilepsy" = "#205A9E",
  "ADORA1^Low Epilepsy" = "#7BC5B4"
)

gene_meta <- list(
  "SCN1B" = list(
    channel_class = "Sodium Channels",
    display_title = "SCN1B (Na\u1d65\u03b21)",
    p_labels = c("P = 0.24", "P = 0.06", "P = 0.001"),
    limits = c(0.15, 1.25), breaks = c(0.3, 0.6, 0.9, 1.2),
    b1 = 1.02, b2 = 1.10, b3 = 1.18, tick = 0.02, text_off = 0.02,
    y_title = "Expression of Gene\n(LogCPM)"
  ),
  "KCNA2" = list(
    channel_class = "Potassium Channels (Voltage-gated)",
    display_title = "KCNA2 (K\u1d651.2)",
    p_labels = c("P = 0.79", "P = 0.02", "P = 0.02"),
    limits = c(0.2, 1.65), breaks = c(0.5, 1.0, 1.5),
    b1 = 1.34, b2 = 1.45, b3 = 1.56, tick = 0.025, text_off = 0.025,
    y_title = "Expression of Genes in Pathway\n(LogCPM)"
  ),
  "KCNC1" = list(
    channel_class = "Potassium Channels (Voltage-gated)",
    display_title = "KCNC1 (K\u1d6511.1)",
    p_labels = c("P = 0.69", "P = 0.28", "P = 0.03"),
    limits = c(0.2, 1.65), breaks = c(0.5, 1.0, 1.5),
    b1 = 1.34, b2 = 1.45, b3 = 1.56, tick = 0.025, text_off = 0.025,
    y_title = NULL
  ),
  "KCNH2" = list(
    channel_class = "Potassium Channels (Voltage-gated)",
    display_title = "KCNH2 (K\u1d653.1)",
    p_labels = c("P = 0.53", "P = 0.013", "P = 0.009"),
    limits = c(0.1, 0.85), breaks = c(0.25, 0.50, 0.75),
    b1 = 0.69, b2 = 0.75, b3 = 0.81, tick = 0.015, text_off = 0.015,
    y_title = NULL
  ),
  "CACNA2D2" = list(
    channel_class = "Calcium Channels",
    display_title = "CACNA2D2 (Ca\u1d65\u03b12\u03b4-2)",
    p_labels = c("P = 0.81", "P = 0.13", "P = 0.014"),
    limits = c(0.2, 2.2), breaks = c(0.5, 1.0, 1.5, 2.0),
    b1 = 1.80, b2 = 1.94, b3 = 2.08, tick = 0.035, text_off = 0.035,
    y_title = "Expression of Genes in Pathway\n(LogCPM)"
  ),
  "CLCN2" = list(
    channel_class = "Chloride Channels",
    display_title = "CLCN2 (CLC-2)",
    p_labels = c("P = 0.64", "P = 0.023", "P = 0.023"),
    limits = c(0.1, 1.1), breaks = c(0.25, 0.50, 0.75, 1.00),
    b1 = 0.88, b2 = 0.95, b3 = 1.02, tick = 0.018, text_off = 0.018,
    y_title = "Expression of Genes in Pathway\n(LogCPM)"
  )
)

plot_single_ion_gene <- function(gene_name) {
  df_sub <- s5e_df %>% filter(Gene == gene_name)
  meta <- gene_meta[[gene_name]]
  labels <- meta$p_labels
  
  p <- ggplot() +
    geom_half_violin(
      data = df_sub,
      aes(x = Group, y = PB, color = Group, fill = Group),
      side = "l", alpha = 0.3, width = 0.75, linewidth = 0.5,
      position = position_dodge(width = 0.15)
    ) +
    geom_half_boxplot(
      data = df_sub,
      aes(x = Group, y = PB, color = Group),
      side = "l", width = 0.2, outlier.shape = NA, linewidth = 0.5, alpha = 1,
      center = TRUE, errorbar.draw = TRUE, errorbar.length = 1,
      position = position_dodge(width = 0.15)
    ) +
    geom_point(
      data = df_sub,
      aes(x = as.numeric(Group) + 0.13, y = PB, color = Group, fill = Group),
      shape = 21, size = 2.0, alpha = 0.85,
      position = position_jitter(width = 0.05, height = 0, seed = 42)
    ) +
    geom_smooth(
      data = df_sub,
      aes(x = Group_Fit, y = PB),
      method = "loess", se = TRUE,
      color = "black", linewidth = 0.5, linetype = "dashed",
      fill = "grey", alpha = 0.3
    ) +
    # Bracket 1: Ctrl vs High
    annotate("segment", x = 1, xend = 2, y = meta$b1, yend = meta$b1, linewidth = 0.4) +
    annotate("segment", x = 1, xend = 1, y = meta$b1, yend = meta$b1 - meta$tick, linewidth = 0.4) +
    annotate("segment", x = 2, xend = 2, y = meta$b1, yend = meta$b1 - meta$tick, linewidth = 0.4) +
    annotate("text", x = 1.5, y = meta$b1 + meta$text_off, label = labels[1], size = 2.5, family = "sans", fontface = "italic") +
    # Bracket 2: Ctrl vs Low
    annotate("segment", x = 1, xend = 3, y = meta$b2, yend = meta$b2, linewidth = 0.4) +
    annotate("segment", x = 1, xend = 1, y = meta$b2, yend = meta$b2 - meta$tick, linewidth = 0.4) +
    annotate("segment", x = 3, xend = 3, y = meta$b2, yend = meta$b2 - meta$tick, linewidth = 0.4) +
    annotate("text", x = 2, y = meta$b2 + meta$text_off, label = labels[2], size = 2.5, family = "sans", fontface = "italic") +
    # Bracket 3: High vs Low
    annotate("segment", x = 2, xend = 3, y = meta$b3, yend = meta$b3, linewidth = 0.4) +
    annotate("segment", x = 2, xend = 2, y = meta$b3, yend = meta$b3 - meta$tick, linewidth = 0.4) +
    annotate("segment", x = 3, xend = 3, y = meta$b3, yend = meta$b3 - meta$tick, linewidth = 0.4) +
    annotate("text", x = 2.5, y = meta$b3 + meta$text_off, label = labels[3], size = 2.5, family = "sans", fontface = "italic") +
    scale_fill_manual(values = color_palette) +
    scale_color_manual(values = color_palette) +
    scale_y_continuous(limits = meta$limits, breaks = meta$breaks) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(size = 9.5, face = "bold.italic", hjust = 0.5, family = "sans"),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = if (!is.null(meta$y_title)) element_text(size = 8.5, family = "sans") else element_blank(),
      axis.text.y = element_text(size = 8.5, colour = "black", family = "sans"),
      legend.position = "none",
      plot.margin = margin(3, 4, 3, 4)
    ) +
    labs(
      title = meta$display_title,
      y = meta$y_title
    )
  
  return(p)
}

p_scn1b    <- plot_single_ion_gene("SCN1B")
p_kcna2    <- plot_single_ion_gene("KCNA2")
p_kcnc1    <- plot_single_ion_gene("KCNC1")
p_kcnh2    <- plot_single_ion_gene("KCNH2")
p_cacna2d2 <- plot_single_ion_gene("CACNA2D2")
p_clcn2    <- plot_single_ion_gene("CLCN2")

# Super-title grobs
make_super_title <- function(title_text) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.25, label = title_text, size = 3.2, fontface = "bold", family = "sans") +
    annotate("segment", x = 0.05, xend = 0.95, y = 0.05, yend = 0.05, linewidth = 0.5) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 0.5), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(0, 2, 2, 2))
}

st_na <- make_super_title("Sodium Channels")
st_k  <- make_super_title("Potassium Channels (Voltage-gated)")
st_ca <- make_super_title("Calcium Channels")
st_cl <- make_super_title("Chloride Channels")

# Combine super titles with plots
col_na <- (st_na / p_scn1b) + plot_layout(heights = c(0.08, 1))
col_k  <- (st_k / (p_kcna2 | p_kcnc1 | p_kcnh2)) + plot_layout(heights = c(0.08, 1))
col_ca <- (st_ca / p_cacna2d2) + plot_layout(heights = c(0.08, 1))
col_cl <- (st_cl / p_clcn2) + plot_layout(heights = c(0.08, 1))

combined_e <- (col_na | col_k | col_ca | col_cl) +
  plot_layout(widths = c(1.1, 2.9, 1.1, 1.1))

# Add bottom legend with properly spaced items
legend_grob <- ggplot() +
  annotate("rect", xmin = 0.52, xmax = 0.54, ymin = 0.35, ymax = 0.65, fill = "#989798") +
  annotate("text", x = 0.55, y = 0.5, label = "Control", size = 2.8, hjust = 0, family = "sans") +
  annotate("rect", xmin = 0.65, xmax = 0.67, ymin = 0.35, ymax = 0.65, fill = "#205A9E") +
  annotate("text", x = 0.68, y = 0.5, label = "ADORA1^High Epilepsy", size = 2.8, hjust = 0, family = "sans") +
  annotate("rect", xmin = 0.82, xmax = 0.84, ymin = 0.35, ymax = 0.65, fill = "#7BC5B4") +
  annotate("text", x = 0.85, y = 0.5, label = "ADORA1^Low Epilepsy", size = 2.8, hjust = 0, family = "sans") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_void()

final_e_plot <- (combined_e / legend_grob) + plot_layout(heights = c(1, 0.05))

ggsave(check_plot_path, plot = final_e_plot, width = 12.5, height = 3.8, dpi = 300)
cat("Check plot saved to:", check_plot_path, "\n")
