# ==============================================================================
# Script: FigS8DE_GABA_EFeatures_Violin.R
# Description: Action potential features in GABA.N neurons across disease groups
#              (Epilepsy vs Control) and ADORA1 expression groups (High vs Low).
# Panel: Supplementary Figure 8D and 8E
# Source Data:
#   - FigS8D_GABA_Epilepsy_vs_Control_EFeatures.csv
#   - FigS8E_GABA_ADORA1High_vs_Low_EFeatures.csv
# Output:
#   - ../../sourcedata-260907/check_plots/check_FigS8DE_GABA_EFeatures_Violin.png
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(gghalves)
  library(patchwork)
})

base_dir <- here::here()  # Project root; use here::here() or setwd() as appropriate
code_dir <- "."  # SuppFig8
data_dir <- file.path("..", "sourcedata", "SuppFig8")
plot_dir <- file.path("..", "sourcedata", "check_plots")

if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

cache_path <- file.path(code_dir, "ps_s8_cache.rds")
s8_cache <- readRDS(cache_path)
gaba_all <- s8_cache$gaba_data

features <- c("AP_UpStroke", "AP_UpDownRatio", "AP_HalfWidth", "AP_Duration")
titles <- c("UpStroke", "Up/Down Ratio", "AP HalfWidth", "AP Duration (mV)")

# ==============================================================================
# 1. Prepare and export Panel D SourceData (Epilepsy vs Control)
# ==============================================================================
df_d <- gaba_all %>%
  filter(Group %in% c("Control", "Epilepsy")) %>%
  mutate(Group_Label = factor(ifelse(Group == "Control", "Ctrl.", "Epi."), levels = c("Ctrl.", "Epi.")))

s8d_export <- df_d %>%
  dplyr::select(Cell_ID, Sample, Group = Group_Label, all_of(features))

write.csv(s8d_export, file.path(data_dir, "FigS8D_GABA_Epilepsy_vs_Control_EFeatures.csv"), row.names = FALSE)

# ==============================================================================
# 2. Prepare and export Panel E SourceData (ADORA1 High vs Low in Epilepsy)
# ==============================================================================
df_e <- gaba_all %>%
  filter(Group == "Epilepsy")

thresh_e <- mean(df_e$ADORA1_Exp)
df_e$ADORA1_Group <- factor(
  ifelse(df_e$ADORA1_Exp > thresh_e, "ADORA1high", "ADORA1low"),
  levels = c("ADORA1high", "ADORA1low")
)

s8e_export <- df_e %>%
  dplyr::select(Cell_ID, Sample, ADORA1_Exp, ADORA1_Group, all_of(features))

write.csv(s8e_export, file.path(data_dir, "FigS8E_GABA_ADORA1High_vs_Low_EFeatures.csv"), row.names = FALSE)
cat("SourceData exported to:", data_dir, "\n")

# ==============================================================================
# 3. Plotting function definitions
# ==============================================================================
cols_d <- c("Ctrl." = "#989798", "Epi." = "#FAC5B1")
cols_e <- c("ADORA1high" = "#205A9E", "ADORA1low" = "#7BC5B4")

plot_d_list <- list()
plot_e_list <- list()

# Generate Panel D sub-plots
for (i in seq_along(features)) {
  f <- features[i]
  sub <- df_d %>% filter(!is.na(.data[[f]]))
  tt <- t.test(sub[[f]] ~ sub$Group_Label)
  p_val <- tt$p.value
  p_label <- sprintf("P = %.3f", p_val)
  
  y_max <- max(sub[[f]], na.rm = TRUE)
  y_min <- min(sub[[f]], na.rm = TRUE)
  y_bracket <- y_max + (y_max - y_min) * 0.12
  y_text <- y_bracket + (y_max - y_min) * 0.08
  
  p <- ggplot(sub, aes(x = as.numeric(Group_Label), y = .data[[f]], fill = Group_Label, color = Group_Label)) +
    geom_half_violin(aes(group = Group_Label), side = "l", alpha = 0.35, linewidth = 0.5, trim = FALSE) +
    geom_half_boxplot(aes(group = Group_Label), side = "l", width = 0.25, alpha = 0.8, fill = "white", color = "black",
                      outlier.shape = NA, errorbar.length = 0.5) +
    geom_point(aes(x = as.numeric(Group_Label) + 0.15), shape = 21, size = 1.8, alpha = 0.85, stroke = 0.4,
               position = position_jitter(width = 0.04, seed = 42)) +
    annotate("segment", x = 1, xend = 2, y = y_bracket, yend = y_bracket, color = "black", linewidth = 0.4) +
    annotate("text", x = 1.5, y = y_text, label = p_label, size = 2.8, fontface = "bold", color = "black") +
    scale_fill_manual(values = cols_d) +
    scale_color_manual(values = cols_d) +
    scale_x_continuous(breaks = 1:2, labels = levels(df_d$Group_Label), limits = c(0.5, 2.5)) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.25))) +
    theme_classic(base_size = 11) +
    labs(title = titles[i], x = NULL, y = NULL) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9.5),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color = "black", size = 8.5),
      legend.position = "none"
    )
  
  plot_d_list[[i]] <- p
}

# Generate Panel E sub-plots
for (i in seq_along(features)) {
  f <- features[i]
  sub <- df_e %>% filter(!is.na(.data[[f]]))
  tt <- t.test(sub[[f]] ~ sub$ADORA1_Group)
  p_val <- tt$p.value
  p_label <- sprintf("P = %.3f", p_val)
  
  y_max <- max(sub[[f]], na.rm = TRUE)
  y_min <- min(sub[[f]], na.rm = TRUE)
  y_bracket <- y_max + (y_max - y_min) * 0.12
  y_text <- y_bracket + (y_max - y_min) * 0.08
  
  p <- ggplot(sub, aes(x = as.numeric(ADORA1_Group), y = .data[[f]], fill = ADORA1_Group, color = ADORA1_Group)) +
    geom_half_violin(aes(group = ADORA1_Group), side = "l", alpha = 0.35, linewidth = 0.5, trim = FALSE) +
    geom_half_boxplot(aes(group = ADORA1_Group), side = "l", width = 0.25, alpha = 0.8, fill = "white", color = "black",
                      outlier.shape = NA, errorbar.length = 0.5) +
    geom_point(aes(x = as.numeric(ADORA1_Group) + 0.15), shape = 21, size = 1.8, alpha = 0.85, stroke = 0.4,
               position = position_jitter(width = 0.04, seed = 42)) +
    annotate("segment", x = 1, xend = 2, y = y_bracket, yend = y_bracket, color = "black", linewidth = 0.4) +
    annotate("text", x = 1.5, y = y_text, label = p_label, size = 2.8, fontface = "bold", color = "black") +
    scale_fill_manual(values = cols_e) +
    scale_color_manual(values = cols_e) +
    scale_x_continuous(breaks = 1:2, labels = levels(df_e$ADORA1_Group), limits = c(0.5, 2.5)) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.25))) +
    theme_classic(base_size = 11) +
    labs(title = titles[i], x = NULL, y = NULL) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 9.5),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color = "black", size = 8.5),
      legend.position = "none"
    )
  
  plot_e_list[[i]] <- p
}

# Assemble Panel D row
p_d_row <- wrap_plots(plot_d_list, nrow = 1)
# Assemble Panel E row
p_e_row <- wrap_plots(plot_e_list, nrow = 1)

# Add panel title labels for D and E
title_d <- ggplot() + 
  annotate("text", x = 0, y = 0.5, label = "Comparison of E-Features in Epilepsy V.S. Control GABA.N", 
           hjust = 0, fontface = "bold", size = 4.2) + 
  theme_void()

title_e <- ggplot() + 
  annotate("text", x = 0, y = 0.5, label = "Comparison of E-Features in ADORA1Low V.S. ADORA1High Epilepsy GABA.N", 
           hjust = 0, fontface = "bold", size = 4.2) + 
  theme_void()

# Assemble final composite figure
p_final <- (title_d / p_d_row / plot_spacer() / title_e / p_e_row) +
  plot_layout(heights = c(0.4, 4, 0.3, 0.4, 4))

check_png <- file.path(plot_dir, "check_FigS8DE_GABA_EFeatures_Violin.png")
ggsave(check_png, plot = p_final, width = 12, height = 7.5, dpi = 300)
cat("Check plot saved to:", check_png, "\n")
