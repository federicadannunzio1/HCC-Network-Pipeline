# =============================================================================
# Title:       Find Switch Gene Outliers from RWR Results
# Description: Identifies switch genes that are also IQR-based outliers in the
#              degree-preserving RWR z-score distribution. Compares with a
#              previously saved outlier list to find new additions.
# Inputs:      - data/switch_module.xlsx (sheet 2): switch module gene list
#              - results/degree_preserving/dp_iqr_outliers.txt
#              - results/degree_preserving/dp_targetseed_probability_z_score.txt
# Outputs:     - results/switch_outliers/all_selected_switch_outliers_iqr.xlsx
# Author:      Federica D'Annunzio
# =============================================================================

rm(list = ls())

library(openxlsx)

# --- Configuration -----------------------------------------------------------
base_path <- "cell_lines/04_random_walk"   # <-- update before running

# --- Load reference switch module ---------------------------------------------
switch_module <- read.xlsx(
  file.path(base_path, "data/switch_module.xlsx"),
  sheet = 2
)

# --- Load degree-preserving IQR outliers --------------------------------------
high_outliers <- read.table(
  file.path(base_path, "results/degree_preserving/dp_iqr_outliers.txt"),
  sep = "\t", header = TRUE, quote = "", check.names = FALSE
)

# --- Filter outliers to switch genes ------------------------------------------
switch_outliers <- high_outliers[which(high_outliers$gene_symbols %in% switch_module$hgnc_symbol), ]

# --- Export current switch outliers -------------------------------------------
out_dir <- file.path(base_path, "results/switch_outliers")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write.xlsx(
  switch_outliers,
  file.path(out_dir, "all_selected_switch_outliers_iqr.xlsx"),
  rowNames = FALSE
)

cat("Switch gene outliers:", nrow(switch_outliers), "\n")
cat("Genes:", switch_outliers$gene_symbols, "\n")
