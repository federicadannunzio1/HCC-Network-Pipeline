# =============================================================================
# Sensitivity Analysis Plots
# =============================================================================
# Description: Generates publication-quality plots for the HMGCR sensitivity
#              analysis:
#   1. Bar plot: DP z-scores of HMGCR direct PPI neighbors (Test A)
#   2. Grouped bar: switch-gene DP z-scores 8 seeds vs 7 seeds (Test B)
#   3. Scatter: rank comparison full vs no-HMGCR (Test B)
#
# Input:
#   - results/sensitivity/hmgcr_neighbors_dp_scores.txt
#   - results/sensitivity/no_hmgcr_dp_results.txt
#   - results/degree_preserving/dp_targetseed_probability_z_score.txt
#   - data/rwr_seed.txt
#
# Output (in results/plots/):
#   - test_a_hmgcr_neighbors_barplot.png
#   - test_b_switch_genes_comparison.png
#   - test_b_rank_comparison_scatter.png
#
# Author: Federica D'Annunzio
# =============================================================================

library(ggplot2)
library(tidyr)

# --- Paths (relative to repo root) -------------------------------------------
base_path  <- "cell_lines/04_random_walk"   # <-- update before running
sens_path  <- file.path(base_path, "results/sensitivity")
plot_path  <- file.path(base_path, "results/plots")
dir.create(plot_path, showWarnings = FALSE, recursive = TRUE)

# ==========================================================================
# PLOT 1: Bar plot Test A - HMGCR neighbors DP z-scores
# ==========================================================================
nb <- read.table(file.path(sens_path, "hmgcr_neighbors_dp_scores.txt"),
                 header = TRUE, sep = "\t", check.names = FALSE, quote = "")

nb <- nb[order(nb$dp_zscore_mod, decreasing = TRUE), ]
nb$gene_symbols <- factor(nb$gene_symbols, levels = nb$gene_symbols)
nb$is_ubiad1 <- ifelse(nb$gene_symbols == "UBIAD1", "UBIAD1", "Other")

p1 <- ggplot(nb, aes(x = gene_symbols, y = dp_zscore_mod, fill = is_ubiad1)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("UBIAD1" = "#E63946", "Other" = "#457B9D"), guide = "none") +
  geom_hline(yintercept = median(nb$dp_zscore_mod, na.rm = TRUE),
             linetype = "dashed", color = "grey40", linewidth = 0.6) +
  labs(title = "Modified z-scores of HMGCR direct PPI neighbors",
       x = NULL, y = "Modified z-score") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        plot.title = element_text(face = "bold", size = 13))

ggsave(file.path(plot_path, "test_a_hmgcr_neighbors_barplot.png"), p1,
       width = 10, height = 5.5, dpi = 300)
cat("Plot 1 saved.\n")

# ==========================================================================
# PLOT 2: Test B - switch genes: 8 seeds vs 7 seeds
# ==========================================================================
dp_full <- read.table(
  file.path(base_path, "results/degree_preserving/dp_targetseed_probability_z_score.txt"),
  header = TRUE, sep = "\t", check.names = FALSE, quote = ""
)

dp_no <- read.table(file.path(sens_path, "no_hmgcr_dp_results.txt"),
                    header = TRUE, sep = "\t", check.names = FALSE, quote = "")

switch_genes <- read.table(file.path(base_path, "data/rwr_seed.txt"))

Q3_full <- quantile(dp_full$dp_zscore_mod, 0.75, na.rm = TRUE)
iqr_full <- IQR(dp_full$dp_zscore_mod, na.rm = TRUE)
thresh_full <- Q3_full + 1.5 * iqr_full

sw_full <- dp_full[dp_full$entrez_geneid %in% switch_genes$V1 &
                   !is.na(dp_full$dp_zscore_mod) &
                   dp_full$dp_zscore_mod > thresh_full, ]
sw_full <- sw_full[order(-sw_full$dp_zscore_mod), ]

sw_merged <- merge(
  sw_full[, c("gene_symbols", "entrez_geneid", "dp_zscore_mod")],
  dp_no[, c("entrez_geneid", "dp_zscore_mod")],
  by = "entrez_geneid", suffixes = c("_8seed", "_7seed")
)
sw_merged <- sw_merged[order(-sw_merged$dp_zscore_mod_8seed), ]
sw_top <- head(sw_merged, 10)

sw_long <- pivot_longer(sw_top,
                        cols = c(dp_zscore_mod_8seed, dp_zscore_mod_7seed),
                        names_to = "condition", values_to = "dp_zscore_mod")
sw_long$condition <- ifelse(sw_long$condition == "dp_zscore_mod_8seed",
                            "Full (8 seeds)", "No HMGCR (7 seeds)")
sw_long$gene_symbols <- factor(sw_long$gene_symbols, levels = sw_top$gene_symbols)

p2 <- ggplot(sw_long, aes(x = gene_symbols, y = dp_zscore_mod, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.65) +
  scale_fill_manual(values = c("Full (8 seeds)" = "#457B9D", "No HMGCR (7 seeds)" = "#E9C46A")) +
  labs(title = "Switch-gene modified z-scores with and without HMGCR seed",
       x = NULL, y = "Modified z-score", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        plot.title = element_text(face = "bold", size = 13),
        legend.position = "top")

ggsave(file.path(plot_path, "test_b_switch_genes_comparison.png"), p2,
       width = 10, height = 5.5, dpi = 300)
cat("Plot 2 saved.\n")

# ==========================================================================
# PLOT 3: Rank scatter - full vs no-HMGCR (top 500)
# ==========================================================================
dp_full_ranked <- dp_full[order(-dp_full$dp_zscore_mod), ]
dp_full_ranked$rank_full <- seq_len(nrow(dp_full_ranked))

dp_no_ranked <- dp_no[order(-dp_no$dp_zscore_mod), ]
dp_no_ranked$rank_no <- seq_len(nrow(dp_no_ranked))

rank_merged <- merge(
  dp_full_ranked[, c("gene_symbols", "entrez_geneid", "rank_full")],
  dp_no_ranked[, c("entrez_geneid", "rank_no")],
  by = "entrez_geneid"
)

rank_merged$group <- ifelse(rank_merged$entrez_geneid %in% switch_genes$V1, "Switch gene", "Other")
rank_merged$group[rank_merged$gene_symbols == "UBIAD1"] <- "UBIAD1"
rank_merged$group <- factor(rank_merged$group, levels = c("Other", "Switch gene", "UBIAD1"))

rank_sub <- rank_merged[rank_merged$rank_full <= 500 | rank_merged$rank_no <= 500, ]

p3 <- ggplot(rank_sub, aes(x = rank_full, y = rank_no, color = group, size = group, alpha = group)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Other" = "grey70", "Switch gene" = "#2A9D8F", "UBIAD1" = "#E63946")) +
  scale_size_manual(values = c("Other" = 0.8, "Switch gene" = 2, "UBIAD1" = 4)) +
  scale_alpha_manual(values = c("Other" = 0.3, "Switch gene" = 0.8, "UBIAD1" = 1)) +
  annotate("text",
           x = rank_sub$rank_full[rank_sub$group == "UBIAD1"] + 10,
           y = rank_sub$rank_no[rank_sub$group == "UBIAD1"] + 15,
           label = "UBIAD1", color = "#E63946", fontface = "bold", size = 4, hjust = 0) +
  labs(title = "Rank comparison: full (8 seeds) vs no-HMGCR (7 seeds)",
       x = "Rank (8 seeds, full)", y = "Rank (7 seeds, no HMGCR)",
       color = NULL, size = NULL, alpha = NULL) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13), legend.position = "top") +
  guides(size = "none", alpha = "none")

ggsave(file.path(plot_path, "test_b_rank_comparison_scatter.png"), p3,
       width = 7, height = 7, dpi = 300)
cat("Plot 3 saved.\n")

cat("All plots saved in:", plot_path, "\n")
