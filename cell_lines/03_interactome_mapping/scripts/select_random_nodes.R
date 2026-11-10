# =============================================================================
# Function: selectRandomNodes
# =============================================================================
# Description: Randomly samples nodes from the interactome matching the degree
#              distribution of the query gene set (degree-preserving null model).
#
# Args:
#   degree_distribution : data.frame with columns "degree" and "frequency"
#   node                : character vector of all node names
#   degree              : named numeric vector of node degrees
#
# Returns:
#   character vector of randomly selected node IDs
#
# Author: Federica D'Annunzio
# =============================================================================

selectRandomNodes <- function(degree_distribution, node, degree) {

  random_nodes <- unlist(lapply(degree_distribution$degree, function(d) {
    freq <- degree_distribution$frequency[degree_distribution$degree == d]
    sample(node[degree == d], freq)
  }))

  return(random_nodes)
}
