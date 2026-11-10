# =============================================================================
# Function: extractLargestConnectedComponent
# =============================================================================
# Description: Extracts the Largest Connected Component (LCC) from an igraph
#              graph object and returns its node and edge counts.
#
# Args:
#   graph : igraph object
#
# Returns:
#   data.frame with columns: num_node, num_edge
#
# Author: Federica D'Annunzio
# =============================================================================

extractLargestConnectedComponent <- function(graph) {

  component <- igraph::components(graph)
  lcc_idx   <- which(component$membership == which.max(component$csize))
  lcc       <- igraph::induced_subgraph(graph, igraph::V(graph)[lcc_idx])

  data.frame(
    num_node = length(igraph::V(lcc)),
    num_edge = length(igraph::E(lcc))
  )
}
