# =============================================================================
# Title:       Find Switch Gene Outliers from RWR Results
# Description: Identifies switch genes that are also IQR-based outliers in the
#              RWR z-score distribution. Compares newly identified outliers
#              against a previously saved outlier list to find additions.
#              Assumes `high_outliers` is already present in the R environment
#              (e.g., computed upstream in the same session).
# Inputs:      - list_switch_combo_cluster2_and_nn_neg_cluster3.xlsx
#                  (sheet 2): reference list of switch module genes
#              - iqr_extreme_switch_genes_z_score.txt: previously identified
#                  IQR outliers
#              - high_outliers (in-memory): data frame of high-scoring outliers
#                  with a `gene_symbols` column
# Outputs:     - all_selected_switch_ouliers_iqr.xlsx: switch genes among the
#                  current IQR outliers
# Author:      Federica D'Annunzio
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(xlsx)

# --- Load reference switch module ---------------------------------------------
switch_module <- read.xlsx(
  "~/Documents/projects/HCC_data/cell_lines/hepatosome_Paola/list_switch_combo_cluster2_and_nn_neg_cluster3.xlsx",
  sheetIndex = 2
)

# --- Filter outliers to switch genes ------------------------------------------
# Subset the in-memory high_outliers to retain only genes in the switch module
switch_outliers <- high_outliers[which(high_outliers$gene_symbols %in% switch_module$hgnc_symbol), ]

# --- Export current switch outliers -------------------------------------------
write.xlsx(
  switch_outliers,
  "~/Documents/projects/HCC_data/cell_lines/random_walk/res/targets_seed/all_selected_switch_ouliers_iqr.xlsx",
  row.names = FALSE
)

# --- Compare with previously saved outlier list -------------------------------
old_iqr_outliers <- read.table(
  "~/Documents/projects/HCC_data/cell_lines/random_walk/res/targets_seed/iqr_extreme_switch_genes_z_score.txt",
  sep    = "\t",
  quote  = "",
  header = TRUE
)

# Genes present in both the old and new outlier lists
x <- intersect(old_iqr_outliers$gene_symbols, switch_outliers$gene_symbols)

# Genes newly added (present in current run but not in the old list)
added_switch <- setdiff(switch_outliers$gene_symbols, x)
ind          <- which(switch_outliers$gene_symbols %in% added_switch)
