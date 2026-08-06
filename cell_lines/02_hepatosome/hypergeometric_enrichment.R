# =============================================================================
# Hypergeometric Enrichment Test
# =============================================================================
# Description: Defines and applies a hypergeometric test to assess whether
#              switch genes are significantly enriched in DisGeNET
#              hepato-biliary disease gene sets.
#
# Input:
#   - mat_universo.txt               : Universe gene list (HGNC symbols)
#   - DisGeNET_proprietà_hepato.txt  : DisGeNET hepato-biliary gene set
#   - *_switch_converted.txt         : Switch gene lists per condition
#
# Output:
#   - output/*_output.txt : Results table (p-value, OR, gene count, gene names)
#
# Author: Federica D'Annunzio
# =============================================================================

# -----------------------------------------------------------------------------
# Function: hypergeom_test
# -----------------------------------------------------------------------------
# Performs a one-tailed hypergeometric test to evaluate enrichment of
# selection genes within a property gene set, relative to a universe.
#
# Args:
#   universe_genes  : character vector — full gene universe
#   property_genes  : character vector — disease-associated genes (DisGeNET)
#   selection_genes : character vector — switch genes for one condition
#   filename_output : string — path to save results
#
# Returns:
#   data.frame with columns: pval, OR, num_gene, name_genes
# -----------------------------------------------------------------------------
hypergeom_test <- function(universe_genes, property_genes, selection_genes, filename_output) {

  universe  <- unique(universe_genes)
  selection <- unique(selection_genes)
  property  <- unique(property_genes)

  u <- length(universe)   # total universe size
  k <- length(selection)  # number of selected genes (switch genes)
  m <- length(property)   # number of genes in property (disease set)
  n <- u - m              # genes NOT in property
  x <- length(intersect(property, selection))  # overlap

  # One-tailed hypergeometric test (enrichment)
  pval <- sum(dhyper(x:k, m, n, k))
  OR   <- (m / u) / (x / k)

  result <- data.frame(
    pval      = pval,
    OR        = OR,
    num_gene  = x,
    name_genes = paste(intersect(property, selection), collapse = ", ")
  )

  write.table(result, filename_output, quote = FALSE, sep = "\t", row.names = FALSE)

  return(result)
}

# =============================================================================
# Application
# =============================================================================

rm(list = ls())
source("hypergeometric_enrichment.R")   # loads hypergeom_test()

path <- "path/to/02_hepatosome"   # <-- update before running

# Load universe and property gene sets
mat_universo     <- read.table(paste0(path, "/hypergeom_data/universo_EZE_DEGs.txt"), check.names = FALSE, quote = "")
disgenet_hepato  <- read.table(paste0(path, "/hypergeom_data/DisGeNET_proprietà_hepato.txt"), check.names = FALSE, quote = "")

# Subset property to genes present in universe
common_genes    <- intersect(disgenet_hepato$V1, mat_universo$V1)
disgenet_hepato <- as.data.frame(disgenet_hepato[disgenet_hepato$V1 %in% common_genes, ])
colnames(disgenet_hepato) <- "V1"

# Load switch gene lists
switch_EZE25 <- read.table(paste0(path, "/hypergeom_data/EZE25_vs_CTRL_switch_converted.txt"), quote = "")

switch_list        <- list(EZE25 = switch_EZE25)
# To run all conditions, use:
# switch_ATO5     <- read.table(paste0(path, "/hypergeom_data/ATO5_vs_CTRL_switch_converted.txt"), quote = "")
# switch_EZE25_ATO <- read.table(paste0(path, "/hypergeom_data/EZE25_ATO_vs_CTRL_switch_converted.txt"), quote = "")
# switch_list <- list(ATO5 = switch_ATO5, EZE25_ATO = switch_EZE25_ATO, EZE25 = switch_EZE25)

# Run hypergeometric test for each condition
results <- lapply(names(switch_list), function(cond) {
  hypergeom_test(
    universe_genes  = mat_universo$V1,
    property_genes  = disgenet_hepato$V1,
    selection_genes = switch_list[[cond]]$V1,
    filename_output = paste0(path, "/output/", cond, "_output.txt")
  )
})
names(results) <- names(switch_list)
