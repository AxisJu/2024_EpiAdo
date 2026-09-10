# ==============================================================================
# Script: Fig1C_KEGG_Metabolic_Pathway_Heatmap.R
# Description: Reproduces the published cell-type-specific metabolic dysfunction
#              heatmap across 3 epileptogenic regions (Figure 1C).
# Publication: Science Translational Medicine
# Panel: Figure 1C
# ==============================================================================

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(readr)
  library(dplyr)
  library(grid)
})

# Define paths
work_dir <- 'E:/Transfer/Projects/EpiAdo/7th_Science Translational Medicine Rev3'
data_file <- file.path(work_dir, 'sourcedata-260907/Figure1/Fig1C_KEGG_Metabolic_Pathway_Heatmap.csv')
check_dir <- file.path(work_dir, 'sourcedata-260907/check_plots')
dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load Data
df <- as.data.frame(read_csv(data_file, show_col_types = FALSE))

target_cols <- c(
  'Glu.N_Cortex', 'Glu.N_Amygdala', 'Glu.N_Hippocampus',
  'GABA.N_Cortex', 'GABA.N_Amygdala', 'GABA.N_Hippocampus',
  'OPCs_Cortex', 'OPCs_Amygdala', 'OPCs_Hippocampus',
  'Oligo._Cortex', 'Oligo._Amygdala', 'Oligo._Hippocampus',
  'Microglia_Cortex', 'Microglia_Amygdala', 'Microglia_Hippocampus',
  'Astro._Cortex', 'Astro._Amygdala', 'Astro._Hippocampus'
)

mat <- as.matrix(df[, target_cols, drop = FALSE])
rownames(mat) <- make.unique(df[['Pathway']])

# 2. Setup Annotations
# Column annotation
col_celltypes <- rep(c('Glu.N', 'GABA.N', 'OPCs', 'Oligo.', 'Microglia', 'Astro.'), each = 3)
col_regions <- rep(c('Cortex', 'Amygdala', 'Hippocampus'), 6)

celltype_colors <- c(
  'Glu.N' = '#d0a9cb',
  'GABA.N' = '#9dcb84',
  'OPCs' = '#cea27b',
  'Oligo.' = '#f3d490',
  'Microglia' = '#dc8384',
  'Astro.' = '#adaed3'
)

region_colors <- c(
  'Cortex' = '#e5a4a6',
  'Amygdala' = '#f5dba3',
  'Hippocampus' = '#adc1d3'
)

ha_col <- HeatmapAnnotation(
  Celltype = col_celltypes,
  BrainRegion = col_regions,
  col = list(Celltype = celltype_colors, BrainRegion = region_colors),
  show_legend = TRUE,
  annotation_name_side = 'right',
  annotation_name_gp = gpar(fontsize = 9, fontface = 'bold'),
  simple_anno_size = unit(3.5, 'mm'),
  gap = unit(1, 'mm')
)

# Row annotation (left side)
row_celltypes <- df[['Screened_Celltype']]
ha_row <- rowAnnotation(
  Screened_In = row_celltypes,
  col = list(Screened_In = celltype_colors),
  show_legend = FALSE,
  show_annotation_name = FALSE,
  simple_anno_size = unit(3.5, 'mm')
)

# 3. Setup Color Mapping
# PiYG: Green for down (negative), Magenta for up (positive), White for 0
col_fun <- colorRamp2(
  seq(-1, 1, length.out = 11),
  c(brewer.pal(11, 'PiYG')[11:7], '#FFFFFF', rev(brewer.pal(11, 'PiYG')[1:5]))
)

# Custom row label color: Highlight 'Purine metabolism' in pink/magenta
row_colors <- ifelse(seq_len(nrow(mat)) == 1, '#C75D89', 'black')
row_fontface <- ifelse(seq_len(nrow(mat)) == 1, 'bold', 'plain')

# 4. Draw Heatmap with Dashed Highlight Boxes
# Dashed box 1: Glu.N (Cols 1-3, Rows 1-8)
# Dashed box 2: GABA.N (Cols 4-6, Rows 9-16)
ht <- Heatmap(
  mat,
  name = 'Cohen d',
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  top_annotation = ha_col,
  left_annotation = ha_row,
  rect_gp = gpar(col = 'white', lwd = 1),
  row_names_side = 'right',
  row_names_gp = gpar(fontsize = 8, col = row_colors, fontface = row_fontface),
  row_labels = df[['Pathway']],
  show_column_names = FALSE,
  column_title = 'Celltype specific metabolism dysfunction in 3 epileptogenic-regions',
  column_title_gp = gpar(fontsize = 11, fontface = 'bold'),
  heatmap_legend_param = list(
    title = 'Cohen d',
    at = c(-1, -0.5, 0, 0.5, 1),
    labels = c('-1', '-0.5', '0', '0.5', '1')
  )
)

# 5. Save Output Files and decorate with dashed boxes
png_path <- file.path(check_dir, 'check_Fig1C_KEGG_Metabolic_Pathway_Heatmap.png')
pdf_path <- file.path(check_dir, 'check_Fig1C_KEGG_Metabolic_Pathway_Heatmap.pdf')

png(png_path, width = 9.5, height = 7, units = 'in', res = 300)
ht_drawn <- draw(ht)

# Decorate with dashed bounding boxes
decorate_heatmap_body('Cohen d', {
  n_rows <- 29
  n_cols <- 18
  
  # Box 1: Glu.N (Cols 1 to 3, Rows 1 to 8)
  x_left_1 <- 0 / n_cols
  x_right_1 <- 3 / n_cols
  y_top_1 <- (n_rows - 0) / n_rows
  y_bottom_1 <- (n_rows - 8) / n_rows
  grid.rect(
    x = unit((x_left_1 + x_right_1)/2, 'npc'),
    y = unit((y_top_1 + y_bottom_1)/2, 'npc'),
    width = unit(x_right_1 - x_left_1, 'npc'),
    height = unit(y_top_1 - y_bottom_1, 'npc'),
    gp = gpar(lty = 2, lwd = 1.5, col = 'black', fill = NA)
  )
  
  # Box 2: GABA.N (Cols 4 to 6, Rows 9 to 16)
  x_left_2 <- 3 / n_cols
  x_right_2 <- 6 / n_cols
  y_top_2 <- (n_rows - 8) / n_rows
  y_bottom_2 <- (n_rows - 16) / n_rows
  grid.rect(
    x = unit((x_left_2 + x_right_2)/2, 'npc'),
    y = unit((y_top_2 + y_bottom_2)/2, 'npc'),
    width = unit(x_right_2 - x_left_2, 'npc'),
    height = unit(y_top_2 - y_bottom_2, 'npc'),
    gp = gpar(lty = 2, lwd = 1.5, col = 'black', fill = NA)
  )
})
dev.off()

pdf(pdf_path, width = 9.5, height = 7)
ht_drawn <- draw(ht)
decorate_heatmap_body('Cohen d', {
  n_rows <- 29
  n_cols <- 18
  
  # Box 1: Glu.N (Cols 1 to 3, Rows 1 to 8)
  x_left_1 <- 0 / n_cols
  x_right_1 <- 3 / n_cols
  y_top_1 <- (n_rows - 0) / n_rows
  y_bottom_1 <- (n_rows - 8) / n_rows
  grid.rect(
    x = unit((x_left_1 + x_right_1)/2, 'npc'),
    y = unit((y_top_1 + y_bottom_1)/2, 'npc'),
    width = unit(x_right_1 - x_left_1, 'npc'),
    height = unit(y_top_1 - y_bottom_1, 'npc'),
    gp = gpar(lty = 2, lwd = 1.5, col = 'black', fill = NA)
  )
  
  # Box 2: GABA.N (Cols 4 to 6, Rows 9 to 16)
  x_left_2 <- 3 / n_cols
  x_right_2 <- 6 / n_cols
  y_top_2 <- (n_rows - 8) / n_rows
  y_bottom_2 <- (n_rows - 16) / n_rows
  grid.rect(
    x = unit((x_left_2 + x_right_2)/2, 'npc'),
    y = unit((y_top_2 + y_bottom_2)/2, 'npc'),
    width = unit(x_right_2 - x_left_2, 'npc'),
    height = unit(y_top_2 - y_bottom_2, 'npc'),
    gp = gpar(lty = 2, lwd = 1.5, col = 'black', fill = NA)
  )
})
dev.off()

cat('Fig 1C Heatmap generated successfully!\n')
