# =============================================================================
# Title:       Random Walk with Restart on PPI Multiplex Network
# Description: Runs Random Walk with Restart (RWR) on a protein-protein
#              interaction (PPI) multiplex network using drug-target genes
#              as seeds. Converts Entrez gene IDs in the results to HGNC
#              symbols and exports the scored gene list.
# Inputs:      - data/interactome.rds          : igraph object of the PPI network
#              - data/adj_normalized_matrix.rds : precomputed normalized adjacency
#              - data/eze_ato_targets.txt       : seed gene list (Entrez IDs)
# Outputs:     - results/rwr_scores/res_targetseed_rwr.txt
# Author:      Federica D'Annunzio
# =============================================================================

rm(list = ls())

# --- Libraries ----------------------------------------------------------------
library(RandomWalkRestartMH)
library(igraph)
library(org.Hs.eg.db)

# --- Paths (relative to repo root) -------------------------------------------
base_path <- "cell_lines/04_random_walk"   # <-- update before running

# --- Load network -------------------------------------------------------------
interactome_graph <- readRDS(file.path(base_path, "data/interactome.rds"))

PPI_MultiplexObject <- create.multiplex(list(interactome = interactome_graph))

AdjMatrixNorm_PPI <- readRDS(file.path(base_path, "data/adj_normalized_matrix.rds"))

# --- Load seed genes ----------------------------------------------------------
rwr_seed_targets <- read.table(file.path(base_path, "data/eze_ato_targets.txt"))
SeedGene <- as.character(unlist(rwr_seed_targets))

cat("Number of seed genes:", length(SeedGene), "\n")
cat("Seed genes:", SeedGene, "\n")

# --- Run RWR ------------------------------------------------------------------
RWR_PPI_Results <- Random.Walk.Restart.Multiplex(
  x               = AdjMatrixNorm_PPI,
  MultiplexObject = PPI_MultiplexObject,
  Seeds           = SeedGene,
  r               = 0.7
)

# --- Convert Entrez IDs to HGNC symbols ---------------------------------------
res_df <- RWR_PPI_Results$RWRM_Results

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys      = res_df$NodeNames,
  column    = "SYMBOL",
  keytype   = "ENTREZID",
  multiVals = "first"
)

# --- Export results -----------------------------------------------------------
gene_symbols_df <- data.frame(
  gene_symbols  = gene_symbols,
  entrez_geneid = res_df$NodeNames,
  score         = res_df$Score
)

res_dir <- file.path(base_path, "results/rwr_scores")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

write.table(
  gene_symbols_df,
  file.path(res_dir, "res_targetseed_rwr.txt"),
  sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE
)
