# =============================================================================
# Title:       Density Plots of RWR Score Distributions
# Description: Visualises the distribution of modified z-scores and raw RWR
#              scores from a Random Walk with Restart analysis. Produces two
#              density plots on a log2 scale:
#                1. Modified z-score distribution with quartile and IQR-outlier
#                   threshold lines, and reference gene markers (BRCA1, RELB).
#                2. Raw RWR score distribution with quartile lines and the
#                   UBIAD1 gene marker.
#              Also classifies UBIAD1 into its z-score quartile.
# Inputs:      - targetseed_probability_z_score.txt: RWR results table with
#                  columns `zscore_mod`, `score`, and `gene_symbols`
# Outputs:     - Two base-R density plots rendered to the active graphics device
# Author:      Federica D'Annunzio
# =============================================================================

rm(list = ls())

# --- Load data ----------------------------------------------------------------
path    <- "~/Documents/HCC_data/cell_lines/random_walk/res/targets_seed"
rwr_res <- read.table(
  paste0(path, "/targetseed_probability_z_score.txt"),
  header      = TRUE,
  check.names = FALSE,
  quote       = ""
)

# --- Quartile classification of UBIAD1 (modified z-score) --------------------
quartiles_z_score_modified <- quantile(rwr_res$zscore_mod, probs = c(0.25, 0.5, 0.75))

z_score_modified_ubiad1 <- rwr_res$zscore_mod[which(rwr_res$gene_symbols == "UBIAD1")]
z_score_modified_brca1  <- rwr_res$zscore_mod[which(rwr_res$gene_symbols == "BRCA1")]
z_score_modified_relb   <- rwr_res$zscore_mod[which(rwr_res$gene_symbols == "RELB")]

# Determine which quartile UBIAD1 falls into
if (z_score_modified_ubiad1 <= quartiles_z_score_modified[1]) {
  quartile <- "1° quartile"
} else if (z_score_modified_ubiad1 <= quartiles_z_score_modified[2]) {
  quartile <- "2° quartile"
} else if (z_score_modified_ubiad1 <= quartiles_z_score_modified[3]) {
  quartile <- "3° quartile"
} else {
  quartile <- "4° quartile"
}

quartile

# --- IQR outlier threshold (modified z-score) ---------------------------------
rwr_res$log_zscore_mod <- log2(rwr_res$zscore_mod + 1)

Q3              <- quantile(rwr_res$zscore_mod, 0.75)
iqr             <- IQR(rwr_res$zscore_mod)
threshold_upper <- Q3 + 1.5 * iqr

# --- Plot 1: Modified z-score distribution (log2 scale) ----------------------
plot(
  density(rwr_res$log_zscore_mod),
  main = "Modified z-score distribution",
  xlab = "Modified z-score (log scale)",
  col  = "blue"
)

# Reference gene markers (UBIAD1 line commented out)
# abline(v = log2(z_score_modified_ubiad1+1), col = "red", lwd = 2, lty = 4)
abline(v = log2(z_score_modified_brca1 + 1), col = "red",      lwd = 2, lty = 4)
abline(v = log2(z_score_modified_relb  + 1), col = "red",      lwd = 2, lty = 4)

# Quartile and outlier threshold lines
abline(v = log2(quartiles_z_score_modified[1] + 1), col = "lightblue", lwd = 2, lty = "dashed")
abline(v = log2(quartiles_z_score_modified[2] + 1), col = "lightblue", lwd = 2, lty = "dashed")
abline(v = log2(quartiles_z_score_modified[3] + 1), col = "lightblue", lwd = 2, lty = "dashed")
abline(v = log2(threshold_upper              + 1), col = "orange",    lwd = 2, lty = "dashed")

# Text labels for reference genes
# text(log2(z_score_modified_ubiad1+1), 0.05, labels = "UBIAD1", col = "red", pos = 4)
text(log2(z_score_modified_brca1 + 1), 0.08, labels = "BRCA1", col = "red", pos = 4)
text(log2(z_score_modified_relb  + 1), 0.06, labels = "RELB",  col = "red", pos = 2)

# Text labels for quartiles and outlier threshold
text(log2(quartiles_z_score_modified[1] + 1), 0.1, labels = "Q1",              col = "darkblue", pos = 3)
text(log2(quartiles_z_score_modified[2] + 1), 0.1, labels = "Q2",              col = "darkblue", pos = 3)
text(log2(quartiles_z_score_modified[3] + 1), 0.1, labels = "Q3",              col = "darkblue", pos = 3)
text(log2(threshold_upper              + 1), 0.2, labels = "outliers \nthreshold", col = "darkblue", pos = 4)

# --- Plot 2: Raw RWR score distribution (log2 scale) -------------------------
quartili            <- quantile(rwr_res$score, probs = c(0.25, 0.5, 0.75))
valore_gene         <- rwr_res$score[rwr_res$gene_symbols == "UBIAD1"]
rwr_res$log_score   <- log2(rwr_res$score + 1e-6)

plot(
  density(rwr_res$log_score),
  main = "RWR score distribution",
  xlab = "RWR score (log scale)",
  col  = "green"
)

# UBIAD1 marker and quartile lines
abline(v = log2(valore_gene  + 1e-6), col = "red",       lwd = 2, lty = 4)       # UBIAD1
abline(v = log2(quartili[1]  + 1e-6), col = "lightblue", lwd = 2, lty = "dashed") # Q1
abline(v = log2(quartili[2]  + 1e-6), col = "lightblue", lwd = 2, lty = "dashed") # Q2 (median)
abline(v = log2(quartili[3]  + 1e-6), col = "lightblue", lwd = 2, lty = "dashed") # Q3

# Text labels
text(log2(valore_gene + 1e-6), 0.05, labels = "UBIAD1", col = "red",      pos = 4)
text(log2(quartili[1] + 1e-6), 0.1,  labels = "Q1",     col = "darkblue", pos = 3)
text(log2(quartili[2] + 1e-6), 0.1,  labels = "Q2",     col = "darkblue", pos = 3)
text(log2(quartili[3] + 1e-6), 0.1,  labels = "Q3",     col = "darkblue", pos = 3)
