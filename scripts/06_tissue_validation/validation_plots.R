# =============================================================================
# Title:       UBIAD1 and HMGCR Expression Validation Plots (Human Tissues)
# Description: Generates expression boxplots and scatter/correlation plots for
#              UBIAD1 and HMGCR in human tissue samples (Simvastatin+Ezetimibe
#              vs. Placebo) from the TSS clinical trial.
# Inputs:      - data/matrix.txt: normalised expression matrix
# Outputs:     - results/UBIAD1_boxplot_tissues.png
#              - results/HMGCR_boxplot_tissues.png
#              - results/Correlation_UBIAD1_HMGCR_*.png
# Author:      Federica D'Annunzio
# =============================================================================

rm(list = ls())

library(dplyr)

# --- Paths (relative to repo root) -------------------------------------------
base_path <- "human_tissues"   # <-- update before running
res_dir   <- file.path(base_path, "results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load expression matrix ---------------------------------------------------
matrix_combo <- read.table(
  file.path(base_path, "data/matrix.txt"),
  header = TRUE, sep = "\t", quote = "", check.names = FALSE, row.names = 1
)

# --- Select and order columns by condition ------------------------------------
cols_simv <- grep("Simvastatin\\+Ezetimibe", colnames(matrix_combo), value = TRUE)
cols_plac <- grep("Placebo",                  colnames(matrix_combo), value = TRUE)
cols_ordered        <- c(cols_simv, cols_plac)
matrix_combo_subset <- matrix_combo %>% select(all_of(cols_ordered))

# --- Gene IDs -----------------------------------------------------------------
ubiad1 <- "ENSG00000120942"
hmgcr  <- "ENSG00000113161"

# =============================================================================
# Function: plot_gene_expression_boxplot
# =============================================================================
plot_gene_expression_boxplot <- function(gene_id, gene_name, matrix_data, cols_simv, cols_plac, box_colors) {
  matrix_data$gene <- gsub("\\..*$", "", rownames(matrix_data))

  gene_mat <- matrix_data[which(matrix_data$gene %in% gene_id), ]
  if (nrow(gene_mat) == 0) {
    warning(paste("Gene", gene_name, "not found in matrix."))
    return(NULL)
  }
  gene_mat <- gene_mat %>% select(-gene)

  M <- length(cols_simv)
  N <- length(cols_plac)
  expression_values <- as.numeric(gene_mat[1, ])

  df <- data.frame(
    expression = expression_values,
    condition  = factor(
      c(rep("Simvastatin + Ezetimibe", M), rep("Placebo", N)),
      levels = c("Simvastatin + Ezetimibe", "Placebo")
    )
  )

  pval <- t.test(expression_values[1:M], expression_values[(M + 1):(M + N)], paired = FALSE)$p.value

  boxplot(
    expression ~ condition, data = df,
    main  = paste0(gene_name, ", p-value = ", format(pval, digits = 2)),
    ylab  = "Gene expression value", xlab  = "Condition",
    col   = box_colors, notch = FALSE,
    pars  = list(boxwex = 0.3, staplewex = 0.6)
  )
}

# --- Boxplots -----------------------------------------------------------------
png(file.path(res_dir, "UBIAD1_boxplot_tissues.png"), width = 1800, height = 1500, res = 300)
plot_gene_expression_boxplot(ubiad1, "UBIAD1", matrix_combo_subset, cols_simv, cols_plac,
                             box_colors = c("purple", "green"))
dev.off()

png(file.path(res_dir, "HMGCR_boxplot_tissues.png"), width = 1800, height = 1500, res = 300)
plot_gene_expression_boxplot(hmgcr, "HMGCR", matrix_combo_subset, cols_simv, cols_plac,
                             box_colors = c("purple", "green"))
dev.off()

# =============================================================================
# Function: plot_gene_correlation_grouped
# =============================================================================
plot_gene_correlation_grouped <- function(gene_id1, gene_name1, gene_id2, gene_name2,
                                          matrix_data, cols_simv, cols_plac, output_path) {
  matrix_data$gene <- gsub("\\..*$", "", rownames(matrix_data))

  exp1 <- as.numeric(matrix_data[matrix_data$gene == gene_id1, !(colnames(matrix_data) %in% "gene")][1, ])
  exp2 <- as.numeric(matrix_data[matrix_data$gene == gene_id2, !(colnames(matrix_data) %in% "gene")][1, ])

  df_corr <- data.frame(
    Gene1 = exp1, Gene2 = exp2,
    Group = ifelse(colnames(matrix_data[, !(colnames(matrix_data) %in% "gene")]) %in% cols_simv,
                   "Simv+Eze", "Placebo")
  )

  res_tot  <- cor.test(df_corr$Gene1, df_corr$Gene2)
  dat_simv <- subset(df_corr, Group == "Simv+Eze")
  res_simv <- cor.test(dat_simv$Gene1, dat_simv$Gene2)
  dat_plac <- subset(df_corr, Group == "Placebo")
  res_plac <- cor.test(dat_plac$Gene1, dat_plac$Gene2)

  png(output_path, width = 2000, height = 2000, res = 300)
  plot(df_corr$Gene1, df_corr$Gene2, type = "n",
       xlab = paste(gene_name1, "expression"), ylab = paste(gene_name2, "expression"),
       main = paste("Correlation:", gene_name1, "vs", gene_name2))
  points(dat_simv$Gene1, dat_simv$Gene2, col = "purple", pch = 19)
  points(dat_plac$Gene1, dat_plac$Gene2, col = "green",  pch = 19)
  abline(lm(Gene2 ~ Gene1, data = dat_simv), col = "purple", lwd = 2, lty = 1)
  abline(lm(Gene2 ~ Gene1, data = dat_plac), col = "green",  lwd = 2, lty = 1)
  abline(lm(Gene2 ~ Gene1, data = df_corr),  col = "black",  lwd = 1, lty = 3)
  legend_text <- c(
    paste0("Simv+Eze: r = ", round(res_simv$estimate, 2), " (p = ", format(res_simv$p.value, digits = 2), ")"),
    paste0("Placebo: r = ",  round(res_plac$estimate, 2), " (p = ", format(res_plac$p.value, digits = 2), ")"),
    paste0("Total: r = ",    round(res_tot$estimate,  2), " (p = ", format(res_tot$p.value,  digits = 2), ")")
  )
  legend("topleft", legend = legend_text, col = c("purple", "green", "black"),
         pch = c(19, 19, NA), lty = c(1, 1, 3), bty = "n", cex = 0.8)
  dev.off()
}

plot_gene_correlation_grouped(ubiad1, "UBIAD1", hmgcr, "HMGCR",
                              matrix_combo_subset, cols_simv, cols_plac,
                              file.path(res_dir, "Correlation_UBIAD1_HMGCR_Grouped.png"))

# =============================================================================
# Function: plot_gene_correlation_treated
# =============================================================================
plot_gene_correlation_treated <- function(gene_id1, gene_name1, gene_id2, gene_name2,
                                          matrix_data, cols_simv, output_path) {
  matrix_data$gene <- gsub("\\..*$", "", rownames(matrix_data))

  exp1 <- as.numeric(matrix_data[matrix_data$gene == gene_id1, !(colnames(matrix_data) %in% "gene")][1, ])
  exp2 <- as.numeric(matrix_data[matrix_data$gene == gene_id2, !(colnames(matrix_data) %in% "gene")][1, ])

  df_all <- data.frame(Gene1 = exp1, Gene2 = exp2,
                       Sample = colnames(matrix_data[, !(colnames(matrix_data) %in% "gene")]))
  df_treated  <- df_all[df_all$Sample %in% cols_simv, ]
  res_treated <- cor.test(df_treated$Gene1, df_treated$Gene2)

  png(output_path, width = 2000, height = 2000, res = 300)
  plot(df_treated$Gene1, df_treated$Gene2, pch = 19, col = "purple",
       xlab = paste(gene_name1, "expression"), ylab = paste(gene_name2, "expression"),
       main = paste("Correlation in Treated Group:", gene_name1, "vs", gene_name2))
  abline(lm(Gene2 ~ Gene1, data = df_treated), col = "purple", lwd = 2)
  stat_text <- paste0("Treated (Simv+Eze):\nr = ", round(res_treated$estimate, 2),
                      "\np = ", format(res_treated$p.value, digits = 2))
  legend("topleft", legend = stat_text, bty = "n", text.col = "purple", cex = 1.0)
  dev.off()
}

plot_gene_correlation_treated(ubiad1, "UBIAD1", hmgcr, "HMGCR",
                              matrix_combo_subset, cols_simv,
                              file.path(res_dir, "Correlation_UBIAD1_HMGCR_TreatedOnly.png"))

# =============================================================================
# Function: plot_gene_correlation_total
# =============================================================================
plot_gene_correlation_total <- function(gene_id1, gene_name1, gene_id2, gene_name2,
                                        matrix_data, output_path) {
  matrix_data$gene <- gsub("\\..*$", "", rownames(matrix_data))

  exp1 <- as.numeric(matrix_data[matrix_data$gene == gene_id1, !(colnames(matrix_data) %in% "gene")][1, ])
  exp2 <- as.numeric(matrix_data[matrix_data$gene == gene_id2, !(colnames(matrix_data) %in% "gene")][1, ])

  df_total <- na.omit(data.frame(Gene1 = exp1, Gene2 = exp2))
  res_total <- cor.test(df_total$Gene1, df_total$Gene2)

  png(output_path, width = 2000, height = 2000, res = 300)
  plot(df_total$Gene1, df_total$Gene2, pch = 19, col = "darkblue",
       xlab = paste(gene_name1, "expression"), ylab = paste(gene_name2, "expression"),
       main = paste("Global Correlation:", gene_name1, "vs", gene_name2))
  abline(lm(Gene2 ~ Gene1, data = df_total), col = "darkblue", lwd = 2)
  stat_text <- paste0("Total Samples (N=", nrow(df_total), "):\nr = ",
                      round(res_total$estimate, 2), "\np = ", format(res_total$p.value, digits = 2))
  legend("topleft", legend = stat_text, bty = "n", text.col = "darkblue", cex = 1.0)
  dev.off()
}

plot_gene_correlation_total(ubiad1, "UBIAD1", hmgcr, "HMGCR",
                            matrix_combo_subset,
                            file.path(res_dir, "Correlation_UBIAD1_HMGCR_Global.png"))
