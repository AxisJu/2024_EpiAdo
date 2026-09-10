# ==============================================================================
# install_packages.R
# EpiAdo: Glutamatergic neuron-associated adenosine deficits in drug-resistant epilepsy
# ==============================================================================
# Run this script once to install all required R packages before executing
# any figure reproduction scripts.
# ==============================================================================

# 1. CRAN Packages
cran_packages <- c(
  # Data manipulation
  "tidyverse", "dplyr", "tidyr", "tibble", "readr", "magrittr",
  # Visualization
  "ggplot2", "ggrepel", "ggpubr", "gghalves", "cowplot", "patchwork",
  "ggsci", "scales", "RColorBrewer", "ggrastr",
  # Heatmaps and grids
  "pheatmap", "corrplot", "grid", "gridExtra",
  # Statistics
  "effsize", "car", "rstatix", "PMCMRplus",
  # Mixed models
  "lme4", "lmerTest",
  # I/O utilities
  "qs", "writexl"
)

# 2. Bioconductor Packages
bioc_packages <- c(
  "Seurat",       # Single-cell analysis
  "DESeq2",       # Differential expression
  "limma",        # Linear models for omics
  "clusterProfiler", # Pathway enrichment (KEGG/GO)
  "enrichplot",   # Enrichment visualization
  "fgsea",        # GSEA
  "ComplexHeatmap", # Advanced heatmaps
  "org.Hs.eg.db", # Human gene annotation
  "org.Mm.eg.db"  # Mouse gene annotation
)

message("=== Installing CRAN packages ===")
install.packages(setdiff(cran_packages, rownames(installed.packages())),
                 dependencies = TRUE)

message("=== Installing Bioconductor packages ===")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(
  setdiff(bioc_packages, rownames(installed.packages())),
  update = FALSE, ask = FALSE
)

message("=== All packages installed successfully! ===")
message("You can now run any figure reproduction script.")
