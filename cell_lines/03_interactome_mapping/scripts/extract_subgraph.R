# =============================================================================
# Function: extractSubgraph
# =============================================================================
# Description: Extracts the subgraph induced by a given gene list from the
#              interactome. Nodes are relabelled with gene symbols.
#              Writes the subgraph edge list to file.
#
# Args:
#   gene_list  : data.frame with columns "ID" (Entrez) and "Symbol" (HGNC)
#   graph      : igraph object (the interactome)
#   node       : character vector of all node names in graph
#   disease    : string label used in the output filename
#   output_dir : directory path for the output file
#
# Returns:
#   igraph subgraph object with vertex names as gene symbols
#
# Author: Federica D'Annunzio
# =============================================================================

extractSubgraph <- function(gene_list, graph, node, disease, output_dir = "results") {

  ids     <- as.character(gene_list$ID)
  symbols <- gene_list$Symbol

  in_interactome  <- which(ids %in% node)
  ids_found       <- ids[in_interactome]
  symbols_found   <- symbols[in_interactome]

  node_df          <- data.frame(ID = ids_found, symbol = symbols_found)
  rownames(node_df) <- ids_found

  subgraph <- igraph::induced_subgraph(graph, ids_found)
  node_df  <- node_df[names(igraph::V(subgraph)), ]

  subgraph <- igraph::set_vertex_attr(subgraph, "name",
                                       index = igraph::V(subgraph), node_df$symbol)

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  igraph::write_graph(subgraph,
                      file.path(output_dir, paste0(disease, "_subnetwork.txt")),
                      format = "ncol")

  return(subgraph)
}
