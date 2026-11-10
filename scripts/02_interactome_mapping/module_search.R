# =============================================================================
# Interactome Module Analysis
# =============================================================================
# Description: Maps switch genes onto the human PPI interactome, extracts the
#              induced subgraph and its LCC, and evaluates statistical
#              significance against a degree-preserving random null model
#              (1000 permutations).
#
# Input:
#   - data/selected_switch_eze_ato.txt  : Switch gene list (Ensembl, Symbol, Entrez)
#   - data/interactome_PNAS2021.txt     : Human interactome edge list
#
# Output:
#   - results/*_subnetwork.txt             : Subgraph edge list
#   - results/eze_combined_plots_high_res.png : Distribution plots
#
# Author: Federica D'Annunzio
# =============================================================================

rm(list = ls())

library(igraph)

options(stringsAsFactors = FALSE)

# --- Configuration -----------------------------------------------------------
base_path <- "cell_lines/03_interactome_mapping"   # <-- update before running

# Load helper functions
source(file.path(base_path, "scripts/compute_degree_distribution.R"))
source(file.path(base_path, "scripts/extract_subgraph.R"))
source(file.path(base_path, "scripts/extract_lcc.R"))
source(file.path(base_path, "scripts/select_random_nodes.R"))
source(file.path(base_path, "scripts/perform_iteration.R"))
source(file.path(base_path, "scripts/compute_statistics.R"))

disease <- "HCC_eze_ato"

# --- Load data ----------------------------------------------------------------
my_genes <- read.table(file.path(base_path, "data/selected_switch_eze_ato.txt"),
                       header = TRUE, sep = "\t", check.names = FALSE, quote = "")
colnames(my_genes) <- c("ensembl", "Symbol", "ID")

interactome <- read.table(file.path(base_path, "data/interactome_PNAS2021.txt"),
                          header = TRUE, sep = "\t", check.names = FALSE, quote = "")

# --- Build interactome graph (largest connected component only) ---------------
graph <- graph_from_data_frame(interactome, directed = FALSE)
graph <- simplify(graph, remove.multiple = TRUE, remove.loops = TRUE)

components_all   <- components(graph)
largest_comp_idx <- which.max(components_all$csize)
graph <- induced_subgraph(graph, which(components_all$membership == largest_comp_idx))

node   <- names(V(graph))
degree <- degree(graph, v = V(graph))

rm(interactome)

# --- Extract disease module and compute LCC -----------------------------------
degree_distribution <- computeDegreeDistribution(my_genes$ID)
subgraph_result     <- extractSubgraph(my_genes, graph, node, disease,
                                       file.path(base_path, "results"))
subgraph_lcc        <- extractLargestConnectedComponent(subgraph_result)

# --- Permutation test (1000 iterations) ---------------------------------------
random_distribution <- performIteration(iter = 1000, graph, degree_distribution,
                                        node, degree)

# --- Plot and compute statistics ----------------------------------------------
res_dir <- file.path(base_path, "results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

png(file.path(res_dir, "eze_combined_plots_high_res.png"),
    width = 4800, height = 1600, res = 600)

par(mfrow = c(1, 3),
    mar   = c(5, 2, 4, 2) + 0.1,
    oma   = c(0, 0, 2, 0))

computeStatistics(random_distribution$num_node_LCC, subgraph_lcc$num_node, "Size of LCC")
computeStatistics(random_distribution$num_edge_LCC, subgraph_lcc$num_edge, "Number of interactions in LCC")
computeStatistics(random_distribution$num_edge_tot, length(E(subgraph_result)), "Number of total interactions")

dev.off()
