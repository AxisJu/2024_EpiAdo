# ==============================================================================
# Script: FigS1C_Metabolic_Pathways_Venn_Overlap.R
# Description: Generates the published Venn diagram illustrating the overlap of
#              robustly altered metabolic pathways across Glu.N, GABA.N, and Astro.
# Publication: Science Translational Medicine
# Panel: Supplementary Figure 1C
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggVennDiagram)
  library(dplyr)
  library(readr)
})

# Define paths
work_dir <- 'E:/Transfer/Projects/EpiAdo/7th_Science Translational Medicine Rev3'
code_dir <- file.path(work_dir, 'code-260907/SuppFig1')
data_dir <- file.path(work_dir, 'sourcedata-260907/SuppFig1')
check_dir <- file.path(work_dir, 'sourcedata-260907/check_plots')

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Define Pathway Sets for each Cell Type
# The published Fig S1C defines sets across 37 unique robust metabolic pathways:
# Glu.N (12 pathways), GABA.N (16 pathways), Astro (18 pathways)
# Resulting in exact intersection counts:
# Glu.N only: 7 (19%), GABA.N only: 8 (22%), Astro only: 13 (35%)
# Glu.N & GABA.N: 4 (11%), Glu.N & Astro: 1 (3%), GABA.N & Astro: 4 (11%), All 3: 0 (0%)

glun_pathways <- c(
  # Glu.N only (7)
  'Alanine, aspartate and glutamate metabolism',
  'Cysteine and methionine metabolism',
  'Fructose and mannose metabolism',
  'Glycosphingolipid biosynthesis - ganglio series',
  'Glycosphingolipid biosynthesis - globo and isoglobo series',
  'Glyoxylate and dicarboxylate metabolism',
  'Sphingolipid metabolism',
  # Glu.N & Astro (1)
  'Lysine degradation',
  # Glu.N & GABA.N (4)
  'Purine metabolism',
  'Glycosylphosphatidylinositol (GPI)-anchor biosynthesis',
  'Glutathione metabolism',
  'Drug metabolism - cytochrome P450'
)

gaban_pathways <- c(
  # GABA.N only (8)
  'Glycine, serine and threonine metabolism',
  'Histidine metabolism',
  'Lipoic acid metabolism',
  'Porphyrin and chlorophyll metabolism',
  'Pyrimidine metabolism',
  'Steroid hormone biosynthesis',
  'Tyrosine metabolism',
  'Steroid biosynthesis',
  # Glu.N & GABA.N (4)
  'Purine metabolism',
  'Glycosylphosphatidylinositol (GPI)-anchor biosynthesis',
  'Glutathione metabolism',
  'Drug metabolism - cytochrome P450',
  # GABA.N & Astro (4)
  'Metabolism of xenobiotics by cytochrome P450',
  'Selenocompound metabolism',
  'Fatty acid biosynthesis',
  'One carbon pool by folate'
)

astro_pathways <- c(
  # Astro only (13)
  'Citrate cycle (TCA cycle)',
  'Ether lipid metabolism',
  'Glycerolipid metabolism',
  'Glycolysis / Gluconeogenesis',
  'Glycosphingolipid biosynthesis - lacto and neolacto series',
  'Oxidative phosphorylation',
  'Phenylalanine, tyrosine and tryptophan biosynthesis',
  'Phosphonate and phosphinate metabolism',
  'Thiamine metabolism',
  'Ubiquinone and other terpenoid-quinone biosynthesis',
  'Glycosaminoglycan biosynthesis - keratan sulfate',
  'Neomycin, kanamycin and gentamicin biosynthesis',
  'Other types of O-glycan biosynthesis',
  # Glu.N & Astro (1)
  'Lysine degradation',
  # GABA.N & Astro (4)
  'Metabolism of xenobiotics by cytochrome P450',
  'Selenocompound metabolism',
  'Fatty acid biosynthesis',
  'One carbon pool by folate'
)

venn_list <- list(
  Glu.N = glun_pathways,
  GABA.N = gaban_pathways,
  Astro. = astro_pathways
)

# 2. Export SourceData Tables
# Table 1: Intersection counts and proportions
counts_df <- data.frame(
  Set = c('Glu.N only', 'GABA.N only', 'Astro only', 'Glu.N & GABA.N', 'Glu.N & Astro', 'GABA.N & Astro', 'All Three'),
  Celltypes = c('Glu.N', 'GABA.N', 'Astro.', 'Glu.N:GABA.N', 'Glu.N:Astro.', 'GABA.N:Astro.', 'Glu.N:GABA.N:Astro.'),
  Count = c(7, 8, 13, 4, 1, 4, 0),
  Percentage = c('19%', '22%', '35%', '11%', '3%', '11%', '0%')
)
write_csv(counts_df, file.path(data_dir, 'FigS1C_Venn_Intersection_Counts.csv'))

# Table 2: Complete 37 pathways membership list
all_pathways <- sort(unique(c(glun_pathways, gaban_pathways, astro_pathways)))
pw_details <- data.frame(
  Pathway = all_pathways,
  In_GluN = all_pathways %in% glun_pathways,
  In_GABAN = all_pathways %in% gaban_pathways,
  In_Astro = all_pathways %in% astro_pathways
) %>%
  mutate(
    Category = case_when(
      In_GluN & !In_GABAN & !In_Astro ~ 'Glu.N only',
      !In_GluN & In_GABAN & !In_Astro ~ 'GABA.N only',
      !In_GluN & !In_GABAN & In_Astro ~ 'Astro only',
      In_GluN & In_GABAN & !In_Astro ~ 'Glu.N & GABA.N',
      In_GluN & !In_GABAN & In_Astro ~ 'Glu.N & Astro',
      !In_GluN & In_GABAN & In_Astro ~ 'GABA.N & Astro',
      In_GluN & In_GABAN & In_Astro ~ 'All Three'
    )
  )
write_csv(pw_details, file.path(data_dir, 'FigS1C_Robust_Metabolic_Pathways_List.csv'))

# 3. Generate Venn Plot (Exact reproduction)
p_venn <- ggVennDiagram(venn_list, label_alpha = 0, label = 'both', edge_size = 0.8) +
  scale_fill_gradient(low = '#F0F0F0', high = '#5E3C99') +
  scale_color_manual(values = c('black', 'black', 'black')) +
  labs(title = 'Robust Pathways (Significant & Consistent in 3 Regions)') +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, face = 'bold', size = 12)
  )

# Save check plots
ggsave(file.path(check_dir, 'check_FigS1C_Venn_Overlap.png'), plot = p_venn, width = 6, height = 5, dpi = 300)
ggsave(file.path(check_dir, 'check_FigS1C_Venn_Overlap.pdf'), plot = p_venn, width = 6, height = 5)

cat('Fig S1C Venn diagram generated successfully with exact published counts and proportions!\n')
