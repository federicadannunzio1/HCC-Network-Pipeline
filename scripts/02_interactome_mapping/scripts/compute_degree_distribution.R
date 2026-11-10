# =============================================================================
# Function: computeDegreeDistribution
# =============================================================================
# Description: Computes the degree distribution for a list of genes within
#              the interactome graph.
#
# Args:
#   gene_ids : character vector of node IDs present in the interactome
#   graph    : igraph object (the interactome)
#   node     : character vector of all node names in graph
#
# Returns:
#   data.frame with columns: degree, frequency
#
# Author: Federica D'Annunzio
# =============================================================================

computeDegreeDistribution <- function(gene_ids, graph = graph, node = node) {

  idx    <- which(node %in% gene_ids)
  deg    <- igraph::degree(graph, v = igraph::V(graph)[idx])
  deg_t  <- table(deg)

  data.frame(
    degree    = as.numeric(names(deg_t)),
    frequency = as.numeric(deg_t)
  )
}
