# =============================================================================
# Hypergeometric Enrichment Test
# =============================================================================
# Description: Defines and applies a hypergeometric test to assess whether
#              switch genes are significantly enriched in DisGeNET
#              hepato-biliary disease gene sets.
#
# Input:
#   - data/hypergeom_data/universo_EZE_DEGs.txt        : Universe gene list
#   - data/hypergeom_data/DisGeNET_proprietà_hepato.txt : DisGeNET gene set
#   - data/hypergeom_data/*_switch_converted.txt        : Switch gene lists
#
# Output:
#   - results/output/*_output.txt : Results table (p-value, OR, gene count)
#
# Author: Federica D'Annunzio
# =============================================================================

rm(list = ls())

# =============================================================================
# Function: hypergeom_test
# =============================================================================
hypergeom_test <- function(universe_genes, property_genes, selection_genes, filename_output) {

  universe  <- unique(universe_genes)
  selection <- unique(selection_genes)
  property  <- unique(property_genes)

  u <- length(universe)
  k <- length(selection)
  m <- length(property)
  n <- u - m
  x <- length(intersect(property, selection))

  pval <- sum(dhyper(x:k, m, n, k))
  OR   <- (m / u) / (x / k)

  result <- data.frame(
    pval       = pval,
    OR         = OR,
    num_gene   = x,
    name_genes = paste(intersect(property, selection), collapse = ", ")
  )

  write.table(result, filename_output, quote = FALSE, sep = "\t", row.names = FALSE)
  return(result)
}

# =============================================================================
# Application
# =============================================================================

path <- "cell_lines/02_hepatosome"   # <-- update before running

mat_universo    <- read.table(file.path(path, "data/hypergeom_data/universo_EZE_DEGs.txt"), check.names = FALSE, quote = "")
disgenet_hepato <- read.table(file.path(path, "data/hypergeom_data/DisGeNET_proprietà_hepato.txt"), check.names = FALSE, quote = "")

common_genes    <- intersect(disgenet_hepato$V1, mat_universo$V1)
disgenet_hepato <- as.data.frame(disgenet_hepato[disgenet_hepato$V1 %in% common_genes, ])
colnames(disgenet_hepato) <- "V1"

switch_EZE25     <- read.table(file.path(path, "data/hypergeom_data/EZE25_vs_CTRL_switch_converted.txt"), quote = "")
switch_ATO5      <- read.table(file.path(path, "data/hypergeom_data/ATO5_vs_CTRL_switch_converted.txt"), quote = "")
switch_EZE25_ATO <- read.table(file.path(path, "data/hypergeom_data/EZE25_ATO_vs_CTRL_switch_converted.txt"), quote = "")

switch_list <- list(ATO5 = switch_ATO5, EZE25_ATO = switch_EZE25_ATO, EZE25 = switch_EZE25)

out_dir <- file.path(path, "results/output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

results <- lapply(names(switch_list), function(cond) {
  hypergeom_test(
    universe_genes  = mat_universo$V1,
    property_genes  = disgenet_hepato$V1,
    selection_genes = switch_list[[cond]]$V1,
    filename_output = file.path(out_dir, paste0(cond, "_output.txt"))
  )
})
names(results) <- names(switch_list)
