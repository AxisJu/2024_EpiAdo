# ==============================================================================
# Script Name: FigS5D_Electrophysiology_Pathways_ADORA1_Correlation.R
# Description: Correlation between ADORA1 expression and 8 action potential/membrane
#              potential pathways in ADORA1^High and ADORA1^Low Epilepsy glutamatergic neurons.
# Author: Antigravity Agent
# Date: 2026-09-09
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# Path definitions
script_dir <- tryCatch(
  dirname(rstudioapi::getActiveDocumentContext()$path),
  error = function(e) "./SuppFig5"
)
cache_path <- file.path(script_dir, "glun_s5_cache.rds")
out_csv_data <- "../sourcedata/SuppFig5/FigS5D_Electrophysiology_Pathways_ADORA1_Correlation.csv"
out_csv_stat <- "../sourcedata/SuppFig5/FigS5D_Electrophysiology_Pathways_Statistics.csv"
check_plot_path <- "../sourcedata/check_plots/check_FigS5D_Electrophysiology_Pathways_ADORA1_Correlation.png"

# 1. Load cached data
if (!file.exists(cache_path)) {
  stop("Cache file not found: ", cache_path)
}
cache <- readRDS(cache_path)
s5d_df <- cache$s5d_df

s5d_names <- c(
  "Neuronal action potential",
  "Membrane depolarization during action potential",
  "Membrane repolarization during action potential",
  "Regulation of presynaptic membrane potential",
  "Regulation of neuronal action potential",
  "Regulation of membrane depolarization during action potential",
  "Regulation of membrane repolarization during action potential",
  "Regulation of postsynaptic membrane potential"
)

# Filter and set factor levels
s5d_df <- s5d_df %>%
  filter(Pathway %in% s5d_names) %>%
  mutate(
    Pathway = factor(Pathway, levels = s5d_names),
    Group = factor(Group, levels = c("ADORA1^High Epilepsy", "ADORA1^Low Epilepsy"))
  )

# 2. Export SourceData CSVs
sourcedata_df <- s5d_df %>%
  transmute(
    Pathway = as.character(Pathway),
    Sample_ID = Sample,
    Group = as.character(Group),
    ADORA1_Expression_LogCPM = round(ADORA1_PB, 4),
    Pathway_Mean_Expression_LogCPM = round(Pathway_PB, 4)
  )
write.csv(sourcedata_df, out_csv_data, row.names = FALSE)

stat_list <- list()
for (pw in s5d_names) {
  for (grp in c("ADORA1^High Epilepsy", "ADORA1^Low Epilepsy")) {
    sub_d <- s5d_df %>% filter(Pathway == pw, Group == grp)
    ct <- cor.test(sub_d$ADORA1_PB, sub_d$Pathway_PB, method = "pearson")
    stat_list[[paste(pw, grp, sep = "_")]] <- data.frame(
      Pathway = pw,
      Group = grp,
      N_Samples = nrow(sub_d),
      Pearson_R = round(ct$estimate, 4),
      T_Statistic = round(ct$statistic, 4),
      DF = ct$parameter,
      P_Value = signif(ct$p.value, 4)
    )
  }
}
sourcedata_stat <- bind_rows(stat_list)
write.csv(sourcedata_stat, out_csv_stat, row.names = FALSE)
cat("SourceData saved to:", out_csv_data, "and", out_csv_stat, "\n")

# 3. Label formatting matching publication
format_stat_label <- function(r, p) {
  # Format r with trailing zeros stripped if integer tenths (e.g. 0.3, 0.1, 0.8)
  r_val <- round(r, 2)
  if (abs(r_val * 10 - round(r_val * 10)) < 1e-4) {
    r_str <- sprintf("%.1f", r_val)
  } else {
    r_str <- sprintf("%.2f", r_val)
  }
  
  if (p < 0.005) {
    p_str <- "P < 1\u00d710\u207b\u2075"
  } else if (p < 0.05) {
    p_str <- sprintf("P = %.2f", p)
  } else {
    p_val <- round(p, 2)
    if (abs(p_val * 10 - round(p_val * 10)) < 1e-4) {
      p_str <- sprintf("P = %.1f", p_val)
    } else {
      p_str <- sprintf("P = %.2f", p_val)
    }
  }
  paste0("\u03c1 = ", r_str, ", ", p_str)
}

plot_single_pathway <- function(pw_name) {
  d_pw <- s5d_df %>% filter(Pathway == pw_name)
  d_high <- d_pw %>% filter(Group == "ADORA1^High Epilepsy")
  d_low  <- d_pw %>% filter(Group == "ADORA1^Low Epilepsy")
  
  ct_high <- cor.test(d_high$ADORA1_PB, d_high$Pathway_PB, method = "pearson")
  ct_low  <- cor.test(d_low$ADORA1_PB, d_low$Pathway_PB, method = "pearson")
  
  lbl_high <- format_stat_label(ct_high$estimate, ct_high$p.value)
  lbl_low  <- format_stat_label(ct_low$estimate, ct_low$p.value)
  
  # Coordinate boundaries
  y_min <- min(d_pw$Pathway_PB, na.rm = TRUE)
  y_max <- max(d_pw$Pathway_PB, na.rm = TRUE)
  y_span <- y_max - y_min
  y_top <- y_max + y_span * 0.32
  y_lbl <- y_max + y_span * 0.16
  y_bot <- y_min - y_span * 0.08
  
  # Left sub-panel (High)
  p_left <- ggplot(d_high, aes(x = ADORA1_PB, y = Pathway_PB)) +
    geom_point(color = "#205A9E", fill = "#205A9E", size = 1.6, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, color = "#205A9E", linetype = "dashed",
                linewidth = 0.5, fill = "#c5d7ed", alpha = 0.45) +
    annotate("text", x = min(d_high$ADORA1_PB), y = y_lbl, label = lbl_high,
             size = 2.4, hjust = 0, family = "sans") +
    scale_y_continuous(limits = c(y_bot, y_top)) +
    theme_classic(base_size = 9) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.margin = margin(0, 0, 0, 0)
    )
  
  # Right sub-panel (Low)
  p_right <- ggplot(d_low, aes(x = ADORA1_PB, y = Pathway_PB)) +
    geom_point(color = "#7BC5B4", fill = "#7BC5B4", size = 1.6, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE, color = "#5aa695", linetype = "dashed",
                linewidth = 0.5, fill = "#d5ebe5", alpha = 0.45) +
    annotate("text", x = min(d_low$ADORA1_PB), y = y_lbl, label = lbl_low,
             size = 2.4, hjust = 0, family = "sans") +
    scale_y_continuous(limits = c(y_bot, y_top)) +
    theme_classic(base_size = 9) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.margin = margin(0, 0, 0, 0)
    )
  
  # Center separator: vertical dashed line
  sep_line <- ggplot() +
    annotate("segment", x = 0, xend = 0, y = y_bot, yend = y_top,
             linetype = "dashed", color = "grey40", linewidth = 0.45) +
    scale_x_continuous(limits = c(-0.01, 0.01), expand = c(0, 0)) +
    scale_y_continuous(limits = c(y_bot, y_top), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
  
  combined_box <- (p_left | sep_line | p_right) +
    plot_layout(widths = c(1, 0.005, 1))
  
  # Title on top
  title_grob <- ggplot() +
    annotate("text", x = 0, y = 0.5, label = pw_name, size = 2.8, fontface = "plain", hjust = 0, family = "sans") +
    scale_x_continuous(limits = c(0, 1)) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 1, 0))
  
  (title_grob / combined_box) + plot_layout(heights = c(0.18, 1))
}

plot_list <- lapply(s5d_names, plot_single_pathway)
combined_grid <- wrap_plots(plot_list, ncol = 4, nrow = 2)

# Global layout with axis titles and legend
y_axis_label <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Mean Expression of Pathway (LogCPM)",
           angle = 90, size = 3.6, family = "sans") +
  theme_void()

x_axis_label <- ggplot() +
  annotate("text", x = 0.4, y = 0.5, label = "ADORA1 Expression (LogCPM)",
           size = 3.6, family = "sans", fontface = "italic") +
  theme_void()

legend_grob <- ggplot() +
  annotate("rect", xmin = 0.05, xmax = 0.12, ymin = 0.35, ymax = 0.65, fill = "#205A9E") +
  annotate("text", x = 0.15, y = 0.5, label = "ADORA1^High Epilepsy", size = 2.8, hjust = 0, family = "sans") +
  annotate("rect", xmin = 0.55, xmax = 0.62, ymin = 0.35, ymax = 0.65, fill = "#7BC5B4") +
  annotate("text", x = 0.65, y = 0.5, label = "ADORA1^Low Epilepsy", size = 2.8, hjust = 0, family = "sans") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_void()

bottom_bar <- (x_axis_label | legend_grob) + plot_layout(widths = c(1, 0.45))

final_d_plot <- (y_axis_label | combined_grid) + plot_layout(widths = c(0.025, 1))
final_d_plot <- (final_d_plot / bottom_bar) + plot_layout(heights = c(1, 0.08))

ggsave(check_plot_path, plot = final_d_plot, width = 11.5, height = 3.8, dpi = 300)
cat("Check plot saved to:", check_plot_path, "\n")
