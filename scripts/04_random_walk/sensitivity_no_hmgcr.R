# =============================================================================
# Sensitivity Analysis: Leave-HMGCR-Out RWR
# =============================================================================
# Description: Two complementary tests addressing reviewer concern about
#              UBIAD1's ranking being driven by HMGCR proximity:
#
#   Test A: Among all direct PPI neighbors of HMGCR, where does UBIAD1 rank
#           in the full (8-seed) degree-preserving analysis?
#
#   Test B: Re-run RWR with 7 seeds (excluding HMGCR) + 1,000 degree-preserving
#           iterations to confirm UBIAD1 remains an outlier.
#
# Input:
#   - data/interactome.rds
#   - data/adj_normalized_matrix.rds
#   - data/eze_ato_targets.txt
#   - data/rwr_seed.txt
#   - results/degree_preserving/dp_targetseed_probability_z_score.txt
#
# Output (in results/sensitivity/):
#   - hmgcr_neighbors_dp_scores.txt
#   - no_hmgcr_dp_results.txt
#   - no_hmgcr_dp_iqr_outliers.txt
#   - sensitivity_summary.txt
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
out_path  <- file.path(base_path, "results/sensitivity")
dir.create(out_path, showWarnings = FALSE, recursive = TRUE)

# --- Load interactome ---------------------------------------------------------
interactome_graph <- readRDS(file.path(base_path, "data/interactome.rds"))

PPI_MultiplexObject <- create.multiplex(list(interactome = interactome_graph))
AdjMatrixNorm_PPI   <- readRDS(file.path(base_path, "data/adj_normalized_matrix.rds"))

# --- Identifiers --------------------------------------------------------------
HMGCR_ENTREZ  <- "3156"
UBIAD1_ENTREZ <- "29914"

# --- Load seeds ---------------------------------------------------------------
all_seeds <- as.character(unlist(read.table(file.path(base_path, "data/eze_ato_targets.txt"))))
cat("All seeds:", all_seeds, "\n")
cat("HMGCR (", HMGCR_ENTREZ, ") is seed:", HMGCR_ENTREZ %in% all_seeds, "\n\n")

switch_genes <- read.table(file.path(base_path, "data/rwr_seed.txt"))

# ==========================================================================
# TEST A: HMGCR-neighbor comparison (uses existing DP results)
# ==========================================================================
cat("=== TEST A: UBIAD1 vs other direct HMGCR neighbors ===\n\n")

hmgcr_neighbors <- names(neighbors(interactome_graph, HMGCR_ENTREZ))
cat("HMGCR has", length(hmgcr_neighbors), "direct PPI neighbors\n")
cat("UBIAD1 among them:", UBIAD1_ENTREZ %in% hmgcr_neighbors, "\n\n")

dp_full <- read.table(
  file.path(base_path, "results/degree_preserving/dp_targetseed_probability_z_score.txt"),
  header = TRUE, sep = "\t", check.names = FALSE, quote = ""
)

dp_neighbors <- dp_full[dp_full$entrez_geneid %in% hmgcr_neighbors, ]
dp_neighbors <- dp_neighbors[order(-dp_neighbors$dp_zscore_mod, na.last = TRUE), ]
dp_neighbors$rank_among_neighbors <- seq_len(nrow(dp_neighbors))

ubiad1_row <- dp_neighbors[dp_neighbors$entrez_geneid == UBIAD1_ENTREZ, ]

cat("UBIAD1 rank among", nrow(dp_neighbors), "HMGCR neighbors:",
    ubiad1_row$rank_among_neighbors, "\n")
cat("UBIAD1 dp_zscore_mod:", ubiad1_row$dp_zscore_mod, "\n\n")

write.table(dp_neighbors, file.path(out_path, "hmgcr_neighbors_dp_scores.txt"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# ==========================================================================
# TEST B: Leave-HMGCR-out RWR + degree-preserving randomisation
# ==========================================================================
cat("=== TEST B: RWR without HMGCR as seed ===\n\n")

seeds_no_hmgcr <- setdiff(all_seeds, HMGCR_ENTREZ)
cat("Seeds without HMGCR (n =", length(seeds_no_hmgcr), "):", seeds_no_hmgcr, "\n\n")

cat("Running observed RWR without HMGCR...\n")
RWR_no_hmgcr <- Random.Walk.Restart.Multiplex(
  x               = AdjMatrixNorm_PPI,
  MultiplexObject = PPI_MultiplexObject,
  Seeds           = seeds_no_hmgcr,
  r               = 0.7
)

observed_no_hmgcr <- RWR_no_hmgcr$RWRM_Results
rownames(observed_no_hmgcr) <- observed_no_hmgcr$NodeNames

all_nodes   <- names(V(interactome_graph))
all_degrees <- degree(interactome_graph)
names(all_degrees) <- all_nodes

seed_degrees_no_hmgcr <- all_degrees[seeds_no_hmgcr]

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
cat("Running", n_iter, "degree-preserving iterations (no HMGCR)...\n")

null_matrix <- matrix(NA, nrow = nrow(observed_no_hmgcr), ncol = n_iter)
rownames(null_matrix) <- observed_no_hmgcr$NodeNames

pb <- txtProgressBar(min = 0, max = n_iter, style = 3)

for (it in 1:n_iter) {
  rand_seeds <- select_degree_preserving_seeds(seed_degrees_no_hmgcr, all_degrees, all_nodes)

  rwr_rand <- Random.Walk.Restart.Multiplex(
    x               = AdjMatrixNorm_PPI,
    MultiplexObject = PPI_MultiplexObject,
    Seeds           = rand_seeds,
    r               = 0.7
  )

  rand_res <- rwr_rand$RWRM_Results
  common   <- intersect(rand_res$NodeNames, rownames(null_matrix))
  null_matrix[common, it] <- rand_res$Score[match(common, rand_res$NodeNames)]

  setTxtProgressBar(pb, it)
}
close(pb)

cat("\nComputing degree-preserving z-scores (no HMGCR)...\n")

null_median <- apply(null_matrix, 1, median, na.rm = TRUE)
null_mad    <- apply(null_matrix, 1, mad, na.rm = TRUE)
null_mean   <- rowMeans(null_matrix, na.rm = TRUE)
null_sd     <- apply(null_matrix, 1, sd, na.rm = TRUE)

dp_zscore_mod_no_hmgcr <- (observed_no_hmgcr$Score -
                            null_median[observed_no_hmgcr$NodeNames]) /
                           null_mad[observed_no_hmgcr$NodeNames]
dp_zscore_mod_no_hmgcr[is.infinite(dp_zscore_mod_no_hmgcr)] <- NA

gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = observed_no_hmgcr$NodeNames,
                       column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")

results_no_hmgcr <- data.frame(
  gene_symbols  = gene_symbols,
  entrez_geneid = observed_no_hmgcr$NodeNames,
  score         = observed_no_hmgcr$Score,
  null_mean     = null_mean[observed_no_hmgcr$NodeNames],
  null_sd       = null_sd[observed_no_hmgcr$NodeNames],
  dp_zscore_mod = dp_zscore_mod_no_hmgcr,
  stringsAsFactors = FALSE
)

results_no_hmgcr <- results_no_hmgcr[order(-results_no_hmgcr$dp_zscore_mod, na.last = TRUE), ]

write.table(results_no_hmgcr, file.path(out_path, "no_hmgcr_dp_results.txt"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

Q3_nh <- quantile(results_no_hmgcr$dp_zscore_mod, 0.75, na.rm = TRUE)
IQR_nh <- IQR(results_no_hmgcr$dp_zscore_mod, na.rm = TRUE)
upper_nh <- Q3_nh + 1.5 * IQR_nh

outliers_no_hmgcr <- results_no_hmgcr[!is.na(results_no_hmgcr$dp_zscore_mod) &
                                       results_no_hmgcr$dp_zscore_mod > upper_nh, ]

write.table(outliers_no_hmgcr, file.path(out_path, "no_hmgcr_dp_iqr_outliers.txt"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

ubiad1_no_hmgcr <- results_no_hmgcr[results_no_hmgcr$entrez_geneid == UBIAD1_ENTREZ, ]
ubiad1_rank_no_hmgcr <- which(results_no_hmgcr$entrez_geneid == UBIAD1_ENTREZ)
ubiad1_is_outlier_no_hmgcr <- UBIAD1_ENTREZ %in% outliers_no_hmgcr$entrez_geneid

switch_outliers_no_hmgcr <- outliers_no_hmgcr[outliers_no_hmgcr$entrez_geneid %in% switch_genes$V1, ]

# --- Save summary ------------------------------------------------------------
sink(file.path(out_path, "sensitivity_summary.txt"))
cat("SENSITIVITY ANALYSIS: UBIAD1 RWR RANKING AND HMGCR DEPENDENCY\n")
cat("==============================================================\n\n")
cat("Test A: UBIAD1 vs other direct HMGCR neighbors\n")
cat("  Total HMGCR neighbors:", nrow(dp_neighbors), "\n")
cat("  UBIAD1 rank among neighbors:", ubiad1_row$rank_among_neighbors, "/", nrow(dp_neighbors), "\n")
cat("  UBIAD1 dp_zscore_mod:", round(ubiad1_row$dp_zscore_mod, 2), "\n\n")
cat("Test B: RWR without HMGCR seed (7 seeds, 1000 DP iterations)\n")
cat("  UBIAD1 overall rank:", ubiad1_rank_no_hmgcr, "/", nrow(results_no_hmgcr), "\n")
cat("  UBIAD1 dp_zscore_mod:", round(ubiad1_no_hmgcr$dp_zscore_mod, 2), "\n")
cat("  UBIAD1 is IQR outlier:", ubiad1_is_outlier_no_hmgcr, "\n")
cat("  Total IQR outliers:", nrow(outliers_no_hmgcr), "\n")
cat("  Switch genes among outliers:", nrow(switch_outliers_no_hmgcr), "\n")
if (nrow(switch_outliers_no_hmgcr) > 0) {
  cat("  Switch outlier names:", switch_outliers_no_hmgcr$gene_symbols, "\n")
}
sink()

cat("\nResults saved in:", out_path, "\n")
