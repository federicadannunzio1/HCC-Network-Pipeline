# =============================================================================
# Title:       UBIAD1 and HMGCR Expression Validation Plots (Human Tissues)
# Description: Generates expression boxplots and scatter/correlation plots for
#              UBIAD1 and HMGCR in human tissue samples (Simvastatin+Ezetimibe
#              vs. Placebo). Includes three correlation plot functions:
#                - Grouped (both conditions, colour-coded)
#                - Treated-only (Simvastatin+Ezetimibe samples)
#                - Global (all samples combined)
# Inputs:      - matrix.txt: normalised expression matrix (Ensembl gene IDs as
#                  row names, samples as columns) from the Combo_placebo
#                  human tissues comparison
# Outputs:     - UBIAD1_boxplot_tissues.png
#              - HMGCR_boxplot_tissues.png
#              - Correlation_UBIAD1_HMGCR_Grouped.png
#              - Correlation_UBIAD1_HMGCR_TreatedOnly.png
#              - Correlation_UBIAD1_HMGCR_Global.png
# Author:      Federica D'Annunzio
# =============================================================================

rm(list = ls())

# --- Libraries ----------------------------------------------------------------
library(dplyr)

# --- Paths --------------------------------------------------------------------
res_dir <- "~/Documents/projects/HCC_data/human_tissues/res_gene_overview"
if (!dir.exists(res_dir)) {
  dir.create(res_dir)
}

# --- Load expression matrix ---------------------------------------------------
matrix_combo <- read.table(
  "~/Documents/projects/HCC_data/human_tissues/Combo_placebo/matrix/matrix.txt",
  header      = TRUE,
  sep         = "\t",
  quote       = "",
  check.names = FALSE,
  row.names   = 1
)

# --- Select and order columns by condition ------------------------------------
cols_to_select <- grep("Placebo|Simvastatin\\+Ezetimibe", colnames(matrix_combo), value = TRUE)
matrix_combo_subset <- matrix_combo %>% select(all_of(cols_to_select))

cols_simv <- grep("Simvastatin\\+Ezetimibe", colnames(matrix_combo), value = TRUE)
cols_plac <- grep("Placebo",                  colnames(matrix_combo), value = TRUE)

# Reorder columns: treated group first, then placebo
cols_ordered        <- c(cols_simv, cols_plac)
matrix_combo_subset <- matrix_combo %>% select(all_of(cols_ordered))

# --- Gene IDs -----------------------------------------------------------------
ubiad1 <- "ENSG00000120942"
hmgcr  <- "ENSG00000113161"
# brca1  <- "ENSG00000012048"  # reserved for future use
# relb   <- "ENSG00000104856"  # reserved for future use
# srebf2 <- "ENSG00000198911"  # reserved for future use

# =============================================================================
# Function: plot_gene_expression_boxplot
# Description: Draws a boxplot of expression values for a single gene across
#              treated and placebo samples, with a Welch t-test p-value in the
#              title.
# =============================================================================
plot_gene_expression_boxplot <- function(gene_id, gene_name, matrix_data, cols_simv, cols_plac, box_colors) {
  matrix_data$gene <- rownames(matrix_data)
  matrix_data$gene <- gsub("\\..*$", "", matrix_data$gene)

  gene_mat <- matrix_data[which(matrix_data$gene %in% gene_id), ]
  if (nrow(gene_mat) == 0) {
    warning(paste("Gene", gene_name, "not found in matrix."))
    return(NULL)
  }

  gene_mat <- gene_mat %>% select(-gene)

  M                 <- length(cols_simv)
  N                 <- length(cols_plac)
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
    expression ~ condition,
    data  = df,
    main  = paste0(gene_name, ", p-value = ", format(pval, digits = 2)),
    ylab  = "Gene expression value",
    xlab  = "Condition",
    col   = box_colors,
    notch = FALSE,
    pars  = list(boxwex = 0.3, staplewex = 0.6)
  )
}

# --- Boxplot: UBIAD1 ----------------------------------------------------------
png(
  "/Users/federicadannunzio/Documents/projects/HCC_data/human_tissues/res_gene_overview/UBIAD1_boxplot_tissues.png",
  width = 1800, height = 1500, res = 300
)
plot_gene_expression_boxplot(ubiad1, "UBIAD1", matrix_combo_subset, cols_simv, cols_plac,
                             box_colors = c("purple", "green"))
dev.off()

# --- Boxplot: HMGCR -----------------------------------------------------------
png(
  "/Users/federicadannunzio/Documents/projects/HCC_data/human_tissues/res_gene_overview/HMGCR_boxplot_tissues.png",
  width = 1800, height = 1500, res = 300
)
plot_gene_expression_boxplot(hmgcr, "HMGCR", matrix_combo_subset, cols_simv, cols_plac,
                             box_colors = c("purple", "green"))
dev.off()
# plot_gene_expression_boxplot(relb, "RELB", matrix_combo_subset, cols_simv, cols_plac,
#                              box_colors = c("magenta", "lightblue"))

# =============================================================================
# Function: plot_gene_correlation_grouped
# Description: Scatter plot of expression for two genes across all samples,
#              colour-coded by treatment group. Computes and displays Pearson
#              correlations for each group and overall, with regression lines.
# =============================================================================
plot_gene_correlation_grouped <- function(gene_id1, gene_name1, gene_id2, gene_name2,
                                          matrix_data, cols_simv, cols_plac, output_path) {
  # Clean gene IDs
  matrix_data$gene <- gsub("\\..*$", "", rownames(matrix_data))

  exp1 <- as.numeric(matrix_data[matrix_data$gene == gene_id1, !(colnames(matrix_data) %in% "gene")][1, ])
  exp2 <- as.numeric(matrix_data[matrix_data$gene == gene_id2, !(colnames(matrix_data) %in% "gene")][1, ])

  df_corr <- data.frame(
    Gene1 = exp1,
    Gene2 = exp2,
    Group = ifelse(
      colnames(matrix_data[, !(colnames(matrix_data) %in% "gene")]) %in% cols_simv,
      "Simv+Eze", "Placebo"
    )
  )

  # Compute correlations: overall, treated, placebo
  res_tot  <- cor.test(df_corr$Gene1, df_corr$Gene2)
  dat_simv <- subset(df_corr, Group == "Simv+Eze")
  res_simv <- cor.test(dat_simv$Gene1, dat_simv$Gene2)
  dat_plac <- subset(df_corr, Group == "Placebo")
  res_plac <- cor.test(dat_plac$Gene1, dat_plac$Gene2)

  png(output_path, width = 2000, height = 2000, res = 300)

  plot(df_corr$Gene1, df_corr$Gene2, type = "n",
       xlab = paste(gene_name1, "expression"),
       ylab = paste(gene_name2, "expression"),
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

  legend("topleft",
         legend = legend_text,
         col    = c("purple", "green", "black"),
         pch    = c(19, 19, NA),
         lty    = c(1, 1, 3),
         bty    = "n",
         cex    = 0.8)

  dev.off()
}

plot_gene_correlation_grouped(
  gene_id1    = ubiad1,
  gene_name1  = "UBIAD1",
  gene_id2    = hmgcr,
  gene_name2  = "HMGCR",
  matrix_data = matrix_combo_subset,
  cols_simv   = cols_simv,
  cols_plac   = cols_plac,
  output_path = paste0(res_dir, "/Correlation_UBIAD1_HMGCR_Grouped.png")
)

# =============================================================================
# Function: plot_gene_correlation_treated
# Description: Scatter plot restricted to the treated (Simvastatin+Ezetimibe)
#              group, with Pearson correlation shown in a top-left legend.
# =============================================================================
plot_gene_correlation_treated <- function(gene_id1, gene_name1, gene_id2, gene_name2,
                                          matrix_data, cols_simv, output_path) {
  # Clean gene IDs
  matrix_data$gene <- gsub("\\..*$", "", rownames(matrix_data))

  # Extract expression for all samples, then filter to treated
  exp1 <- as.numeric(matrix_data[matrix_data$gene == gene_id1, !(colnames(matrix_data) %in% "gene")][1, ])
  exp2 <- as.numeric(matrix_data[matrix_data$gene == gene_id2, !(colnames(matrix_data) %in% "gene")][1, ])

  df_all <- data.frame(
    Gene1  = exp1,
    Gene2  = exp2,
    Sample = colnames(matrix_data[, !(colnames(matrix_data) %in% "gene")])
  )

  df_treated  <- df_all[df_all$Sample %in% cols_simv, ]
  res_treated <- cor.test(df_treated$Gene1, df_treated$Gene2)

  png(output_path, width = 2000, height = 2000, res = 300)

  plot(df_treated$Gene1, df_treated$Gene2,
       pch  = 19,
       col  = "purple",
       xlab = paste(gene_name1, "expression"),
       ylab = paste(gene_name2, "expression"),
       main = paste("Correlation in Treated Group:", gene_name1, "vs", gene_name2))

  abline(lm(Gene2 ~ Gene1, data = df_treated), col = "purple", lwd = 2)

  stat_text <- paste0(
    "Treated (Simv+Eze):\n",
    "r = ", round(res_treated$estimate, 2),
    "\np = ", format(res_treated$p.value, digits = 2)
  )

  legend("topleft", legend = stat_text, bty = "n", text.col = "purple", cex = 1.0)

  dev.off()
}

plot_gene_correlation_treated(
  gene_id1    = ubiad1,
  gene_name1  = "UBIAD1",
  gene_id2    = hmgcr,
  gene_name2  = "HMGCR",
  matrix_data = matrix_combo_subset,
  cols_simv   = cols_simv,
  output_path = paste0(res_dir, "/Correlation_UBIAD1_HMGCR_TreatedOnly.png")
)

# =============================================================================
# Function: plot_gene_correlation_total
# Description: Scatter plot using all samples (treated + placebo) to compute a
#              global Pearson correlation between two genes.
# =============================================================================
plot_gene_correlation_total <- function(gene_id1, gene_name1, gene_id2, gene_name2,
                                        matrix_data, output_path) {
  # Clean gene IDs (strip Ensembl version suffixes)
  matrix_data$gene <- gsub("\\..*$", "", rownames(matrix_data))

  # Extract expression across all available samples
  exp1 <- as.numeric(matrix_data[matrix_data$gene == gene_id1, !(colnames(matrix_data) %in% "gene")][1, ])
  exp2 <- as.numeric(matrix_data[matrix_data$gene == gene_id2, !(colnames(matrix_data) %in% "gene")][1, ])

  df_total <- data.frame(Gene1 = exp1, Gene2 = exp2)
  df_total <- na.omit(df_total)

  res_total <- cor.test(df_total$Gene1, df_total$Gene2)

  png(output_path, width = 2000, height = 2000, res = 300)

  plot(df_total$Gene1, df_total$Gene2,
       pch  = 19,
       col  = "darkblue",
       xlab = paste(gene_name1, "expression"),
       ylab = paste(gene_name2, "expression"),
       main = paste("Global Correlation:", gene_name1, "vs", gene_name2))

  abline(lm(Gene2 ~ Gene1, data = df_total), col = "darkblue", lwd = 2)

  stat_text <- paste0(
    "Total Samples (N=", nrow(df_total), "):\n",
    "r = ", round(res_total$estimate, 2),
    "\np = ", format(res_total$p.value, digits = 2)
  )

  legend("topleft", legend = stat_text, bty = "n", text.col = "darkblue", cex = 1.0)

  dev.off()
}

plot_gene_correlation_total(
  gene_id1    = ubiad1,
  gene_name1  = "UBIAD1",
  gene_id2    = hmgcr,
  gene_name2  = "HMGCR",
  matrix_data = matrix_combo_subset,
  output_path = paste0(res_dir, "/Correlation_UBIAD1_HMGCR_Global.png")
)
