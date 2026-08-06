# =============================================================================
# Title:       UBIAD1 and HMGCR Expression Boxplots and Correlation Plots
#              (HCC Cell Lines)
# Description: Generates expression boxplots and scatter/correlation plots for
#              UBIAD1 and HMGCR in HCC cell lines treated with EZE25+ATO vs.
#              CTRL. Includes three correlation plot functions:
#                - Grouped (treated + placebo, colour-coded)
#                - Treated-only (two implementations)
#                - Global (all samples)
# Inputs:      - matrix.txt: normalised expression matrix (Ensembl gene IDs as
#                  row names, samples as columns) from EZE25_ATO_vs_CTRL
# Outputs:     - UBIAD1_boxplot_cells.png
#              - HMGCR_boxplot_cells.png
#              - Correlation_UBIAD1_HMGCR_Grouped.png
#              - Correlation_UBIAD1_HMGCR_TREATED_ONLY.png
#              - Correlation_UBIAD1_HMGCR_TreatedOnly.png
#              - Correlation_UBIAD1_HMGCR_Global.png
# Author:      Federica D'Annunzio
# =============================================================================

rm(list = ls())

# --- Libraries ----------------------------------------------------------------
library(dplyr)

# --- Paths --------------------------------------------------------------------
dir_res <- "~/Documents/projects/HCC_data/cell_lines/genes_overview/ubiad1/ubiad1_vs_hmgcr/"
if (!dir.exists(res_dir)) {
  dir.create(res_dir)
}

# --- Load expression matrix ---------------------------------------------------
matrix_combo <- read.table(
  "~/Documents/projects/HCC_data/cell_lines/cells/EZE25_ATO_vs_CTRL/matrix/matrix.txt",
  header      = TRUE,
  sep         = "\t",
  quote       = "",
  check.names = FALSE,
  row.names   = 1
)

# --- Select and order columns by condition ------------------------------------
cols_to_select <- grep("CTRL|EZE25_ATO", colnames(matrix_combo), value = TRUE)
matrix_combo_subset <- matrix_combo %>% select(all_of(cols_to_select))

cols_simv <- grep("EZE25_ATO", colnames(matrix_combo), value = TRUE)
cols_plac <- grep("CTRL",      colnames(matrix_combo), value = TRUE)

# Reorder columns: treated group first, then placebo
cols_ordered        <- c(cols_simv, cols_plac)
matrix_combo_subset <- matrix_combo %>% select(all_of(cols_ordered))

# --- Gene IDs -----------------------------------------------------------------
ubiad1 <- "ENSG00000120942"
hmgcr  <- "ENSG00000113161"
# brca1 <- "ENSG00000012048"   # reserved for future use
# relb  <- "ENSG00000104856"   # reserved for future use
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
png(paste0(dir_res, "/UBIAD1_boxplot_cells.png"), width = 1800, height = 1500, res = 300)
plot_gene_expression_boxplot(ubiad1, "UBIAD1", matrix_combo_subset, cols_simv, cols_plac,
                             box_colors = c("orange", "turquoise"))
dev.off()

# --- Boxplot: HMGCR -----------------------------------------------------------
png(paste0(dir_res, "/HMGCR_boxplot_cells.png"), width = 1800, height = 1500, res = 300)
plot_gene_expression_boxplot(hmgcr, "HMGCR", matrix_combo_subset, cols_simv, cols_plac,
                             box_colors = c("orange", "turquoise"))
dev.off()
# plot_gene_expression_boxplot(relb, "RELB", matrix_combo_subset, cols_simv, cols_plac,
#                              box_colors = c("magenta", "lightblue"))

# =============================================================================
# Function: plot_gene_correlation_grouped
# Description: Scatter plot of expression for two genes across all samples,
#              colour-coded by treatment group, with per-group and overall
#              Pearson correlation reported in the legend.
# =============================================================================
plot_gene_correlation_grouped <- function(gene_id1, gene_name1, gene_id2, gene_name2,
                                          matrix_data, cols_simv, cols_plac, output_path) {
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

  cor_gen <- cor.test(df_corr$Gene1, df_corr$Gene2)

  png(output_path, width = 1800, height = 1800, res = 300)

  plot(df_corr$Gene1, df_corr$Gene2, type = "n",
       xlab = paste(gene_name1, "expression"),
       ylab = paste(gene_name2, "expression"),
       main = paste("Correlation:", gene_name1, "vs", gene_name2),
       sub  = paste0("Overall Pearson r = ", round(cor_gen$estimate, 2),
                     " (p = ", format(cor_gen$p.value, digits = 2), ")"))

  points(df_corr$Gene1[df_corr$Group == "Simv+Eze"], df_corr$Gene2[df_corr$Group == "Simv+Eze"],
         col = "orange",    pch = 19)
  points(df_corr$Gene1[df_corr$Group == "Placebo"],  df_corr$Gene2[df_corr$Group == "Placebo"],
         col = "turquoise", pch = 19)

  legend("topleft", legend = c("Ato+Eze", "Placebo"),
         col = c("orange", "turquoise"), pch = 19)

  abline(lm(Gene2 ~ Gene1, data = subset(df_corr, Group == "Simv+Eze")),
         col = "orange",    lwd = 1, lty = 2)
  abline(lm(Gene2 ~ Gene1, data = subset(df_corr, Group == "Placebo")),
         col = "turquoise", lwd = 1, lty = 2)

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
# Function: plot_gene_correlation_treated (version 1)
# Description: Scatter plot restricted to treated (Simv+Eze) samples only,
#              with Pearson correlation displayed as a subtitle.
# =============================================================================
plot_gene_correlation_treated <- function(gene_id1, gene_name1, gene_id2, gene_name2,
                                          matrix_data, cols_simv, output_path) {
  # Clean gene IDs (strip Ensembl version suffixes)
  matrix_data$gene <- gsub("\\..*$", "", rownames(matrix_data))

  # Extract expression for treated samples only
  exp1 <- as.numeric(matrix_data[matrix_data$gene == gene_id1, cols_simv][1, ])
  exp2 <- as.numeric(matrix_data[matrix_data$gene == gene_id2, cols_simv][1, ])

  df_corr  <- data.frame(Gene1 = exp1, Gene2 = exp2)
  cor_res  <- cor.test(df_corr$Gene1, df_corr$Gene2, method = "pearson")

  png(output_path, width = 1800, height = 1800, res = 300)

  plot(df_corr$Gene1, df_corr$Gene2,
       pch  = 19,
       col  = "orange",
       xlab = paste(gene_name1, "(Normalized Expression)"),
       ylab = paste(gene_name2, "(Normalized Expression)"),
       main = paste("Treated Group Correlation:", gene_name1, "vs", gene_name2),
       sub  = paste0("Pearson r = ", round(cor_res$estimate, 2),
                     " | p-value = ", format(cor_res$p.value, digits = 3)))

  abline(lm(Gene2 ~ Gene1, data = df_corr), col = "darkorange", lwd = 2)
  grid(lty = "dotted", col = "lightgray")

  dev.off()
}

plot_gene_correlation_treated(
  gene_id1    = ubiad1,
  gene_name1  = "UBIAD1",
  gene_id2    = hmgcr,
  gene_name2  = "HMGCR",
  matrix_data = matrix_combo_subset,
  cols_simv   = cols_simv,
  output_path = paste0(res_dir, "/Correlation_UBIAD1_HMGCR_TREATED_ONLY.png")
)

# =============================================================================
# Function: plot_gene_correlation_treated (version 2 — redefines the function)
# Description: Alternate treated-only scatter plot that builds the full sample
#              data frame first and then filters to the treated group. Pearson
#              correlation statistics are shown in a top-left legend.
# =============================================================================
plot_gene_correlation_treated <- function(gene_id1, gene_name1, gene_id2, gene_name2,
                                          matrix_data, cols_simv, output_path) {
  # Clean gene IDs
  matrix_data$gene <- gsub("\\..*$", "", rownames(matrix_data))

  # Extract expression for all samples
  exp1 <- as.numeric(matrix_data[matrix_data$gene == gene_id1, !(colnames(matrix_data) %in% "gene")][1, ])
  exp2 <- as.numeric(matrix_data[matrix_data$gene == gene_id2, !(colnames(matrix_data) %in% "gene")][1, ])

  df_all <- data.frame(
    Gene1  = exp1,
    Gene2  = exp2,
    Sample = colnames(matrix_data[, !(colnames(matrix_data) %in% "gene")])
  )

  # Filter to treated group
  df_treated  <- df_all[df_all$Sample %in% cols_simv, ]
  res_treated <- cor.test(df_treated$Gene1, df_treated$Gene2)

  png(output_path, width = 2000, height = 2000, res = 300)

  plot(df_treated$Gene1, df_treated$Gene2,
       pch  = 19,
       col  = "orange",
       xlab = paste(gene_name1, "expression"),
       ylab = paste(gene_name2, "expression"),
       main = paste("Correlation in Treated Group:", gene_name1, "vs", gene_name2))

  abline(lm(Gene2 ~ Gene1, data = df_treated), col = "orange", lwd = 2)

  stat_text <- paste0(
    "Treated (Ato+Eze):\n",
    "r = ", round(res_treated$estimate, 2),
    "\np = ", format(res_treated$p.value, digits = 2)
  )

  legend("topleft", legend = stat_text, bty = "n", text.col = "orange", cex = 1.0)

  dev.off()
}

plot_gene_correlation_treated(
  gene_id1    = ubiad1,
  gene_name1  = "UBIAD1",
  gene_id2    = hmgcr,
  gene_name2  = "HMGCR",
  matrix_data = matrix_combo_subset,
  cols_simv   = cols_simv,
  output_path = paste0(dir_res, "/Correlation_UBIAD1_HMGCR_TreatedOnly.png")
)

# =============================================================================
# Function: plot_gene_correlation_total
# Description: Scatter plot using all samples (treated + placebo combined) to
#              compute a global Pearson correlation between two genes.
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
