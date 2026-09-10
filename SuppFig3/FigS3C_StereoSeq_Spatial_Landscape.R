# ==============================================================================
# Script Name: FigS3C_StereoSeq_Spatial_Landscape.R
# Description: Spatial mapping of cell-type distributions, purine metabolism
#              activity, and ADORA1 expression levels in cell-binned Stereo-seq
#              sections from control and epileptic human hippocampus
#              (Supplementary Figure 3C).
# ==============================================================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(ggrastr)

out_data_dir <- "../sourcedata/SuppFig3"
out_check_dir <- "../sourcedata/check_plots"

cache_data <- file.path(out_data_dir, "stereoseq_s3_cache.rds")

if (!file.exists(cache_data)) {
  stop("Please run FigS3A first to generate stereoseq_s3_cache.rds cache.")
}

cat("Loading Stereo-seq cache...\n")
combined_data <- readRDS(cache_data)

celltype_levels <- c(
  "Glu.N", "GABA.N", "OPCs", "Oligo.", 
  "Mo/Mφ", "Microglia", "Astro.", "Endo.", "Epi.", "LowQuality"
)

color_res.coarse <- c(
  "Glu.N" = "#d0a9cb",
  "GABA.N" = "#9dcb84",
  "OPCs" = "#cea27b",
  "Oligo." = "#f3d490",
  "Mo/Mφ" = "#adaed3",
  "Microglia" = "#adc1d3",
  "Astro." = "#d25f61",
  "Endo." = "#dc8384",
  "Epi." = "#c5d3e0",
  "LowQuality" = "#EEEEEE"
)

# Prepare slice data with anatomical alignment:
# T1035 (Control-2) was anatomically cropped to hippocampus proper (Y > -21500)
# and rotated 90 degrees clockwise (X' = Y, Y' = -X) in the publication figure.
df_list <- list()

for (sm in c("T1036", "T1035", "T940", "T941")) {
  sub <- combined_data %>% filter(SampleID == sm)
  if (sm == "T1035") {
    sub <- sub %>% filter(Y > -21500)
    # 90 deg clockwise rotation
    orig_x <- sub$X
    orig_y <- sub$Y
    sub$X <- orig_y
    sub$Y <- -orig_x
  }
  df_list[[sm]] <- sub
}

# Coordinate parameters for each aligned slice
# 4000 units = 2 mm
samples_params <- list(
  T1036 = list(
    name = "Slice: Control-1",
    sm = "T1036",
    x_lim = c(2000, 27000), y_lim = c(-26000, -1000),
    sb_x = 18000, sb_xend = 22000, sb_y = -24500,
    txt_x = 20000, txt_y = -23500
  ),
  T1035 = list(
    name = "Slice: Control-2",
    sm = "T1035",
    x_lim = c(-21000, -2000), y_lim = c(-23000, -2000),
    sb_x = -7000, sb_xend = -3000, sb_y = -21500,
    txt_x = -5000, txt_y = -20500
  ),
  T940 = list(
    name = "Slice: Epilepsy-1",
    sm = "T940",
    x_lim = c(10000, 38000), y_lim = c(-49000, -21000),
    sb_x = 28000, sb_xend = 32000, sb_y = -47000,
    txt_x = 30000, txt_y = -46000
  ),
  T941 = list(
    name = "Slice: Epilepsy-2",
    sm = "T941",
    x_lim = c(9000, 37000), y_lim = c(-48000, -20000),
    sb_x = 28000, sb_xend = 32000, sb_y = -46000,
    txt_x = 30000, txt_y = -45000
  )
)

# Export SourceData
cat("Exporting SourceData...\n")
set.seed(42)
spatial_sub <- bind_rows(df_list) %>%
  group_by(Slice) %>%
  sample_n(min(n(), 5000)) %>%
  ungroup() %>%
  dplyr::select(
    Slice, SampleID, SpotID,
    Spatial_X = X, Spatial_Y = Y,
    Celltype,
    Purine_Metabolism = PurineMetabolism,
    ADORA1_Expression = ADORA1
  )

write.csv(
  spatial_sub,
  file = file.path(out_data_dir, "FigS3C_StereoSeq_Spatial_Coordinates_and_Scores.csv"),
  row.names = FALSE
)

cat("Generating Spatial Landscape Plots...\n")
plot_list_celltype <- list()
plot_list_purine <- list()
plot_list_adora1 <- list()

for (sm in names(samples_params)) {
  params <- samples_params[[sm]]
  df_sm <- df_list[[sm]]
  df_sm$Celltype <- factor(df_sm$Celltype, levels = celltype_levels)
  
  # 1. Cell-type plot
  p_cell <- ggplot() +
    theme_void() +
    coord_fixed(xlim = params$x_lim, ylim = params$y_lim) +
    rasterise(
      geom_point(data = df_sm[df_sm$Celltype == "LowQuality", ], aes(x = X, y = Y), size = 0.15, color = "#EEEEEE"),
      dpi = 200
    ) +
    rasterise(
      geom_point(data = df_sm[df_sm$Celltype == "Oligo.", ], aes(x = X, y = Y), size = 0.15, color = "#f3d490", alpha = 0.5),
      dpi = 200
    ) +
    rasterise(
      geom_point(data = df_sm[!df_sm$Celltype %in% c("LowQuality", "Oligo."), ], aes(x = X, y = Y, color = Celltype), size = 0.25),
      dpi = 200
    ) +
    scale_color_manual(values = color_res.coarse) +
    geom_segment(aes(x = params$sb_x, xend = params$sb_xend, y = params$sb_y, yend = params$sb_y), linewidth = 0.8, color = "black") +
    annotate("text", x = params$txt_x, y = params$txt_y, label = "2 mm", size = 2.5, fontface = "bold", color = "black") +
    theme(
      legend.position = if (sm == "T941") "right" else "none",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold")
    ) +
    labs(title = params$name)
  
  # 2. Purine metabolism plot
  df_purine <- df_sm[order(df_sm$PurineMetabolism), ]
  p_purine <- ggplot() +
    theme_void() +
    coord_fixed(xlim = params$x_lim, ylim = params$y_lim) +
    rasterise(
      geom_point(data = df_purine[df_purine$Celltype != "Glu.N", ], aes(x = X, y = Y), size = 0.15, color = "#DDDDDD"),
      dpi = 200
    ) +
    rasterise(
      geom_point(data = df_purine[df_purine$Celltype == "Glu.N", ], aes(x = X, y = Y, color = PurineMetabolism), size = 0.4),
      dpi = 200
    ) +
    scale_color_gradientn(
      colours = c(brewer.pal(9, "PuBu")[2:9]),
      limits = c(0, 0.4),
      name = "Metabolism\nActivity"
    ) +
    theme(
      legend.position = if (sm == "T941") "right" else "none",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8)
    )
  
  # 3. ADORA1 expression plot
  df_ado <- df_sm[order(df_sm$ADORA1), ]
  p_ado <- ggplot() +
    theme_void() +
    coord_fixed(xlim = params$x_lim, ylim = params$y_lim) +
    rasterise(
      geom_point(data = df_ado[df_ado$Celltype != "Glu.N", ], aes(x = X, y = Y), size = 0.15, color = "#DDDDDD"),
      dpi = 200
    ) +
    rasterise(
      geom_point(data = df_ado[df_ado$Celltype == "Glu.N", ], aes(x = X, y = Y, color = ADORA1), size = 0.4),
      dpi = 200
    ) +
    scale_color_gradientn(
      colours = c(brewer.pal(9, "PuRd")[2:9]),
      limits = c(0, 5),
      name = "Gene\nExpression"
    ) +
    theme(
      legend.position = if (sm == "T941") "right" else "none",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8)
    )
  
  plot_list_celltype[[sm]] <- p_cell
  plot_list_purine[[sm]]   <- p_purine
  plot_list_adora1[[sm]]   <- p_ado
}

cat("Combining into 4x3 Spatial Landscape figure...\n")
grid_plots <- c(
  plot_list_celltype,
  plot_list_purine,
  plot_list_adora1
)

p_spatial_landscape <- wrap_plots(grid_plots, ncol = 4, nrow = 3)

ggsave(
  filename = file.path(out_check_dir, "check_FigS3C_Spatial_Landscape.png"),
  plot = p_spatial_landscape,
  width = 16,
  height = 12,
  dpi = 300
)

cat("FigS3C completed successfully.\n")
