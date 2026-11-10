# =============================================================================
# Function: performIteration
# =============================================================================
# Description: Runs the random permutation loop to build the null distribution
#              for LCC size, LCC edge count, and total edge count.
#
# Args:
#   iter                : integer — number of permutation iterations (default: 1000)
#   graph               : igraph object (the interactome)
#   degree_distribution : data.frame from computeDegreeDistribution()
#   node                : character vector of all node names
#   degree              : named numeric vector of node degrees
#
# Returns:
#   data.frame with columns: num_node_LCC, num_edge_LCC, num_edge_tot
#
# Author: Federica D'Annunzio
# =============================================================================

performIteration <- function(iter = 1000, graph, degree_distribution, node, degree) {

  v_lcc  <- numeric(iter)
  e_lcc  <- numeric(iter)
  e_tot  <- numeric(iter)

  for (i in seq_len(iter)) {
    random_nodes    <- selectRandomNodes(degree_distribution, node, degree)
    random_subgraph <- igraph::induced_subgraph(graph, random_nodes)
    random_lcc      <- extractLargestConnectedComponent(random_subgraph)

    v_lcc[i] <- random_lcc$num_node
    e_lcc[i] <- random_lcc$num_edge
    e_tot[i] <- length(igraph::E(random_subgraph))
  }

  data.frame(num_node_LCC = v_lcc, num_edge_LCC = e_lcc, num_edge_tot = e_tot)
}
