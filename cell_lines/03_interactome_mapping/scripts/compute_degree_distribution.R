# =============================================================================
# Function: computeDegreeDistribution
# =============================================================================
# Description: Computes the degree distribution for a list of genes within
#              the global interactome graph.
#
# Args:
#   gene_ids : character vector of node IDs present in the interactome
#
# Returns:
#   data.frame with columns: degree, frequency
#
# Requires (from global env): graph, node
#
# Author: Federica D'Annunzio
# =============================================================================

computeDegreeDistribution <- function(gene_ids) {

  idx    <- which(node %in% gene_ids)
  deg    <- degree(graph, v = V(graph)[idx])
  deg_t  <- table(deg)

  data.frame(
    degree    = as.numeric(names(deg_t)),
    frequency = as.numeric(deg_t)
  )
}
