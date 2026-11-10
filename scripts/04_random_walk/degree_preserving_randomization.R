# =============================================================================
# Degree-Preserving Null Model for RWR
# =============================================================================
# Description: For each of 1,000 iterations, generates a random seed set by
#              replacing each real drug-target seed with a node drawn uniformly
#              from all interactome nodes sharing the same degree (with +/-1
#              or +/-2 tolerance as fallback). Produces degree-preserving (DP)
#              z-scores (standard and MAD-based modified) and identifies
#              outliers.
#
# Input:
#   - data/interactome.rds          : igraph PPI network
#   - data/adj_normalized_matrix.rds: precomputed adjacency matrix
#   - data/eze_ato_targets.txt      : seed gene list (Entrez IDs)
#   - data/rwr_seed.txt             : switch gene list for outlier detection
#
# Output (in results/degree_preserving/):
#   - dp_targetseed_probability_z_score.txt
#   - dp_iqr_outliers.txt
#   - dp_extreme_z3_outliers.txt
#
# Author: Federica D'Annunzio
# =============================================================================

rm(list = ls())

library(RandomWalkRestartMH)
library(igraph)
library(org.Hs.eg.db)

set.seed(42)

# --- Paths (relative to repo root) -------------------------------------------
base_path <- "cell_lines/04_random_walk"   # <-- update before running
res_path  <- file.path(base_path, "results/degree_preserving")
dir.create(res_path, showWarnings = FALSE, recursive = TRUE)

# --- Load interactome & adjacency matrix --------------------------------------
interactome_graph <- readRDS(file.path(base_path, "data/interactome.rds"))

PPI_MultiplexObject <- create.multiplex(list(interactome = interactome_graph))
AdjMatrixNorm_PPI   <- readRDS(file.path(base_path, "data/adj_normalized_matrix.rds"))

# --- Load drug-target seeds ---------------------------------------------------
rwr_seed_targets <- read.table(file.path(base_path, "data/eze_ato_targets.txt"))
SeedGene         <- as.character(unlist(rwr_seed_targets))

cat("Number of seed genes:", length(SeedGene), "\n")
cat("Seed genes:", SeedGene, "\n\n")

# --- 1. Observed RWR ----------------------------------------------------------
cat("Running observed RWR...\n")
RWR_observed <- Random.Walk.Restart.Multiplex(
  x               = AdjMatrixNorm_PPI,
  MultiplexObject = PPI_MultiplexObject,
  Seeds           = SeedGene,
  r               = 0.7
)

observed_scores <- RWR_observed$RWRM_Results
rownames(observed_scores) <- observed_scores$NodeNames

# --- 2. Degree-preserving null model ------------------------------------------
all_nodes   <- names(V(interactome_graph))
all_degrees <- degree(interactome_graph)
names(all_degrees) <- all_nodes

seed_degrees <- all_degrees[SeedGene]
cat("Seed gene degrees:", seed_degrees, "\n\n")

select_degree_preserving_seeds <- function(seed_degrees, all_degrees, all_nodes) {
  random_seeds <- character(length(seed_degrees))
  used <- character(0)

  for (i in seq_along(seed_degrees)) {
    target_deg <- seed_degrees[i]
    candidates <- all_nodes[all_degrees == target_deg]
    candidates <- setdiff(candidates, used)

    if (length(candidates) == 0) {
      candidates <- all_nodes[abs(all_degrees - target_deg) <= 1]
      candidates <- setdiff(candidates, used)
    }
    if (length(candidates) == 0) {
      candidates <- all_nodes[abs(all_degrees - target_deg) <= 2]
      candidates <- setdiff(candidates, used)
    }

    random_seeds[i] <- sample(candidates, 1)
    used <- c(used, random_seeds[i])
  }

  return(random_seeds)
}

n_iter <- 1000
cat("Running", n_iter, "degree-preserving random RWR iterations...\n")

null_scores_matrix <- matrix(NA, nrow = nrow(observed_scores), ncol = n_iter)
rownames(null_scores_matrix) <- observed_scores$NodeNames

pb <- txtProgressBar(min = 0, max = n_iter, style = 3)

for (it in 1:n_iter) {
  random_seeds <- select_degree_preserving_seeds(seed_degrees, all_degrees, all_nodes)

  rwr_random <- Random.Walk.Restart.Multiplex(
    x               = AdjMatrixNorm_PPI,
    MultiplexObject = PPI_MultiplexObject,
    Seeds           = random_seeds,
    r               = 0.7
  )

  rand_res <- rwr_random$RWRM_Results
  common   <- intersect(rand_res$NodeNames, rownames(null_scores_matrix))
  null_scores_matrix[common, it] <- rand_res$Score[match(common, rand_res$NodeNames)]

  setTxtProgressBar(pb, it)
}
close(pb)

# --- 3. Compute degree-preserving z-scores ------------------------------------
cat("\nComputing degree-preserving z-scores...\n")

null_mean   <- rowMeans(null_scores_matrix, na.rm = TRUE)
null_sd     <- apply(null_scores_matrix, 1, sd, na.rm = TRUE)
null_median <- apply(null_scores_matrix, 1, median, na.rm = TRUE)
null_mad    <- apply(null_scores_matrix, 1, mad, na.rm = TRUE)

dp_zscore     <- (observed_scores$Score - null_mean[observed_scores$NodeNames]) /
                  null_sd[observed_scores$NodeNames]

dp_zscore_mod <- (observed_scores$Score - null_median[observed_scores$NodeNames]) /
                  null_mad[observed_scores$NodeNames]
dp_zscore_mod[is.infinite(dp_zscore_mod)] <- NA

# --- 4. Convert to gene symbols -----------------------------------------------
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = observed_scores$NodeNames,
                       column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")

# --- 5. Assemble and save results ---------------------------------------------
results <- data.frame(
  gene_symbols  = gene_symbols,
  entrez_geneid = observed_scores$NodeNames,
  score         = observed_scores$Score,
  null_mean     = null_mean[observed_scores$NodeNames],
  null_sd       = null_sd[observed_scores$NodeNames],
  dp_zscore     = dp_zscore,
  dp_zscore_mod = dp_zscore_mod,
  stringsAsFactors = FALSE
)

results <- results[order(-results$dp_zscore_mod, na.last = TRUE), ]

write.table(results, file.path(res_path, "dp_targetseed_probability_z_score.txt"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# --- 6. Outlier detection -----------------------------------------------------
Q1_dp <- quantile(results$dp_zscore_mod, 0.25, na.rm = TRUE)
Q3_dp <- quantile(results$dp_zscore_mod, 0.75, na.rm = TRUE)
IQR_dp <- Q3_dp - Q1_dp
upper_bound_dp <- Q3_dp + 1.5 * IQR_dp

high_outliers_dp <- results[!is.na(results$dp_zscore_mod) &
                            results$dp_zscore_mod > upper_bound_dp, ]

extreme_dp <- results[!is.na(results$dp_zscore_mod) &
                      results$dp_zscore_mod > 3, ]

write.table(high_outliers_dp, file.path(res_path, "dp_iqr_outliers.txt"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

write.table(extreme_dp, file.path(res_path, "dp_extreme_z3_outliers.txt"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# --- 7. Switch genes in outliers ----------------------------------------------
switch_genes <- read.table(file.path(base_path, "data/rwr_seed.txt"))
switch_in_dp_iqr <- high_outliers_dp[high_outliers_dp$entrez_geneid %in% switch_genes$V1, ]

cat("\nSwitch genes in degree-preserving IQR outliers:", nrow(switch_in_dp_iqr), "\n")
cat("Switch gene outliers:", switch_in_dp_iqr$gene_symbols, "\n")

cat("\n--- Summary ---\n")
cat("Total genes scored:", nrow(results), "\n")
cat("IQR outliers (dp):", nrow(high_outliers_dp), "\n")
cat("Z > 3 outliers (dp):", nrow(extreme_dp), "\n")
cat("Results saved in:", res_path, "\n")
