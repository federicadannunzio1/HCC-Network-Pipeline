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

  component <- components(graph)
  lcc_idx   <- which(component$membership == which.max(component$csize))
  lcc       <- induced_subgraph(graph, V(graph)[lcc_idx])

  data.frame(
    num_node = length(V(lcc)),
    num_edge = length(E(lcc))
  )
}
