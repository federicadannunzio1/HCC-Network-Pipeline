# =============================================================================
# Function: extractSubgraph
# =============================================================================
# Description: Extracts the subgraph induced by a given gene list from the
#              global interactome. Nodes are relabelled with gene symbols.
#              Writes the subgraph edge list to file.
#
# Args:
#   gene_list : data.frame with columns "ID" (Entrez) and "Symbol" (HGNC)
#
# Returns:
#   igraph subgraph object with vertex names as gene symbols
#
# Requires (from global env): graph, node, disease
#
# Author: Federica D'Annunzio
# =============================================================================

extractSubgraph <- function(gene_list) {

  ids     <- as.character(gene_list$ID)
  symbols <- gene_list$Symbol

  # Keep only genes present in the interactome
  in_interactome  <- which(ids %in% node)
  ids_found       <- ids[in_interactome]
  symbols_found   <- symbols[in_interactome]

  node_df          <- data.frame(ID = ids_found, symbol = symbols_found)
  rownames(node_df) <- ids_found

  # Induce subgraph on matching nodes
  subgraph <- induced_subgraph(graph, ids_found)
  node_df  <- node_df[names(V(subgraph)), ]

  # Relabel vertices with gene symbols
  subgraph <- set_vertex_attr(subgraph, "name", index = V(subgraph), node_df$symbol)

  # Save edge list
  write_graph(subgraph, paste0("res/", disease, "_subnetwork.txt"), format = "ncol")

  return(subgraph)
}
