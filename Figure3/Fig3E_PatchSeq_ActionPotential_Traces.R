# ==============================================================================
# Script Name: Fig3E_PatchSeq_ActionPotential_Traces.R
# Description: Generates Figure 3E: Action potential waveforms overlay 
#              from ADORA1-High and ADORA1-Low neurons (Glu.N and GABA.N).
# Author: Auto-refactored for EpiAdo STM Revision
# Date: 2026-09-08
# ==============================================================================

suppressPackageStartupMessages({
  library(readABF)
  library(pracma)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})

# Define paths
out_dir_data <- "../sourcedata/Figure3"
out_dir_plot <- "../sourcedata/check_plots"

dir.create(out_dir_data, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_plot, recursive = TRUE, showWarnings = FALSE)

get_traces <- function(dir_path, indices = NULL, min_h = 0) {
  files <- list.files(dir_path, full.names = TRUE, pattern = "\\.abf$")
  if (!is.null(indices)) {
    files <- files[indices]
  }
  traces_list <- list()
  for (f in files) {
    cname <- tools::file_path_sans_ext(basename(f))
    tryCatch({
      abf <- readABF(f)
      for (j in seq_along(abf$data)) {
        df <- as.data.frame(abf, sweep = j)
        peaks <- pracma::findpeaks(x = df[, 2], minpeakheight = min_h)
        if (!is.null(peaks) && nrow(peaks) > 0 && nrow(peaks) <= 3) {
          for (k in 1:nrow(peaks)) {
            peak_df <- df[peaks[k, 3]:peaks[k, 4], ]
            names(peak_df) <- c("Time", "mV")
            peak_df$Time_ms <- round((peak_df$Time - df[peaks[k, 2], 1]) * 1000, 4)
            peak_df$mV <- round(peak_df$mV, 2)
            peak_df$Cell_Sweep <- paste(cname, j, k, sep = "_")
            traces_list[[length(traces_list) + 1]] <- peak_df[, c("Time_ms", "mV", "Cell_Sweep")]
          }
        }
      }
    }, error = function(e) {})
  }
  bind_rows(traces_list)
}

cat("Extracting Glu.N AP traces...\n")
glu_low <- get_traces("Z:/2023_EpiAdo/data/PatchClampTraces/processed/ADORA1 Low GluN/", c(6, 8, 9, 12, 13), min_h = 0)
glu_low$Ado_Group <- "ADORA1Low"
glu_low$Celltype <- "Glu.N"

glu_high <- get_traces("Z:/2023_EpiAdo/data/PatchClampTraces/processed/ADORA1 High GluN/", c(2, 5, 8), min_h = 30)
glu_high$Ado_Group <- "ADORA1High"
glu_high$Celltype <- "Glu.N"

cat("Extracting GABA.N AP traces...\n")
gaba_low_files <- list.files("Z:/2023_EpiAdo/data/PatchClampTraces/processed/ADORA1 Low GABAN/", pattern = "\\.abf$")
gaba_high_files <- list.files("Z:/2023_EpiAdo/data/PatchClampTraces/processed/ADORA1 High GABAN/", pattern = "\\.abf$")

gaba_low <- get_traces("Z:/2023_EpiAdo/data/PatchClampTraces/processed/ADORA1 Low GABAN/", 1:length(gaba_low_files), min_h = 0)
gaba_low$Ado_Group <- "ADORA1Low"
gaba_low$Celltype <- "GABA.N"

gaba_high <- get_traces("Z:/2023_EpiAdo/data/PatchClampTraces/processed/ADORA1 High GABAN/", 1:length(gaba_high_files), min_h = 30)
gaba_high$Ado_Group <- "ADORA1High"
gaba_high$Celltype <- "GABA.N"

all_trace_df <- bind_rows(glu_low, glu_high, gaba_low, gaba_high)

# Filter within visual window (-1.5 to 3.0 ms) to keep SourceData clean
all_trace_filtered <- all_trace_df %>%
  filter(Time_ms >= -1.5 & Time_ms <= 3.0)

csv_file <- file.path(out_dir_data, "Fig3E_PatchSeq_ActionPotential_Traces.csv")
write.csv(all_trace_filtered, csv_file, row.names = FALSE)
cat("SourceData exported to:", csv_file, "\n")
cat("Summary of traces:\n")
print(table(all_trace_filtered$Celltype, all_trace_filtered$Ado_Group))

# Colors matching paper Figure 3E:
# ADORA1High: dark teal #2B5B6E
# ADORA1Low: green #7FBC41
color_map <- c("ADORA1High" = "#2B5B6E", "ADORA1Low" = "#7FBC41")

p_glu <- ggplot(all_trace_filtered %>% filter(Celltype == "Glu.N"), 
                aes(x = Time_ms, y = mV, group = Cell_Sweep, color = Ado_Group)) +
  geom_line(alpha = 0.35, linewidth = 0.4) +
  scale_color_manual(values = color_map) +
  scale_x_continuous(limits = c(-1.5, 3), breaks = -1:3) +
  scale_y_continuous(limits = c(-60, 70), breaks = seq(-50, 70, 20)) +
  labs(title = "Glu.N", x = "Time from Peak (ms)", y = "mV") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5, color = "black"),
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 9, color = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length = unit(1.5, "mm"),
    legend.position = "none"
  )

p_gaba <- ggplot(all_trace_filtered %>% filter(Celltype == "GABA.N"), 
                 aes(x = Time_ms, y = mV, group = Cell_Sweep, color = Ado_Group)) +
  geom_line(alpha = 0.35, linewidth = 0.4) +
  scale_color_manual(values = color_map) +
  scale_x_continuous(limits = c(-1.5, 3), breaks = -1:3) +
  scale_y_continuous(limits = c(-60, 70), breaks = seq(-50, 70, 20)) +
  labs(title = "GABA.N", x = "Time from Peak (ms)", y = NULL) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5, color = "black"),
    axis.title = element_text(size = 10, color = "black"),
    axis.text = element_text(size = 9, color = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length = unit(1.5, "mm"),
    legend.position = "none"
  )

# Legend dummy
p_legend <- ggplot(data.frame(x = c(1, 2), y = c(1, 1), Ado_Group = c("ADORA1High", "ADORA1Low")), 
                   aes(x = x, y = y, color = Ado_Group)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = color_map, labels = c(
    expression(italic(ADORA1)^High),
    expression(italic(ADORA1)^Low)
  )) +
  theme_void() +
  theme(legend.position = "top", legend.title = element_blank(), legend.text = element_text(size = 9))

leg <- get_legend(p_legend)
trace_grid <- plot_grid(p_glu, p_gaba, nrow = 1, align = "h")

final_p3e <- plot_grid(
  leg,
  trace_grid,
  ncol = 1,
  rel_heights = c(0.12, 1)
)

check_plot_file <- file.path(out_dir_plot, "check_Fig3E_PatchSeq_ActionPotential_Traces.png")
ggsave(check_plot_file, final_p3e, width = 6.2, height = 3.6, dpi = 300)
cat("Check plot saved to:", check_plot_file, "\n")