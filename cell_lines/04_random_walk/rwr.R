# =============================================================================
# Title:       Random Walk with Restart on PPI Multiplex Network
# Description: Runs Random Walk with Restart (RWR) on a protein-protein
#              interaction (PPI) multiplex network using drug-target genes
#              as seeds. Converts Entrez gene IDs in the results to HGNC
#              symbols and exports the scored gene list.
# Inputs:      - interactome.rds         : igraph object of the PPI network
#              - adj_normalized_matrix.rds: precomputed normalized adjacency matrix
#              - eze_ato_targets.txt      : seed gene list (Entrez IDs)
# Outputs:     - res_targetseed_rwr.txt  : RWR scores with HGNC gene symbols
# Author:      Federica D'Annunzio
# =============================================================================

rm(list = ls())

# --- Libraries ----------------------------------------------------------------
library(RandomWalkRestartMH)
library(igraph)
library(biomaRt)
library(org.Hs.eg.db)

# --- Paths --------------------------------------------------------------------
path_files <- "~/Documents/HCC_data/cell_lines/interactome_mapping/code"
path_res   <- "~/Documents/HCC_data/cell_lines/random_walk/res"

# --- Load network -------------------------------------------------------------
interactome_graph <- readRDS(paste0(path_files, "/interactome.rds"))

# Create multiplex object from the interactome graph
PPI_MultiplexObject <- create.multiplex(list(interactome = interactome_graph))

# Load precomputed normalized adjacency matrix
# (commented-out lines show how to recompute it from scratch if needed)
# AdjMatrix_PPI     <- compute.adjacency.matrix(PPI_MultiplexObject)
# AdjMatrixNorm_PPI <- normalize.multiplex.adjacency(AdjMatrix_PPI)
AdjMatrixNorm_PPI <- readRDS("~/Documents/HCC_data/cell_lines/random_walk/res/adj_normalized_matrix.rds")

# --- Load seed genes ----------------------------------------------------------
# Seed genes are the drug targets for the ezetimibe + atorvastatin combination
rwr_seed_targets <- read.table("~/Documents/HCC_data/cell_lines/random_walk/file/eze_ato_targets.txt")
SeedGene <- as.character(unlist(rwr_seed_targets))

# --- Run RWR ------------------------------------------------------------------
RWR_PPI_Results <- Random.Walk.Restart.Multiplex(
  x              = AdjMatrixNorm_PPI,
  MultiplexObject = PPI_MultiplexObject,
  Seeds          = SeedGene,
  r              = 0.7
)

# --- Convert Entrez IDs to HGNC symbols ---------------------------------------
res_df <- RWR_PPI_Results$RWRM_Results

# Approach 1: biomaRt conversion (kept for reference)
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
res_converted_neg_nn <- getBM(
  attributes = c('entrezgene_id', "hgnc_symbol"),
  filters    = "hgnc_symbol",
  values     = res_df$NodeNames,
  mart       = ensembl
)

# Approach 2: org.Hs.eg.db mapping (used for the exported table)
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys     = res_df$NodeNames,
  column   = "SYMBOL",
  keytype  = "ENTREZID",
  multiVals = "first"
)

# --- Export results -----------------------------------------------------------
gene_symbols_df <- data.frame(
  gene_symbols,
  entrez_geneid = res_df$NodeNames,
  score         = res_df$Score
)

write.table(
  gene_symbols_df,
  paste0(path_res, "/res_targetseed_rwr.txt"),
  sep       = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote     = FALSE
)
