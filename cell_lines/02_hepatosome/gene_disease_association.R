# =============================================================================
# Title:       Gene-Disease Association Table Builder
# Description: Reads DisGeNET enrichment result files, extracts the Term
#              (disease name) and Genes columns, splits the semicolon-delimited
#              gene lists, and reshapes each file into a long-format
#              Disease-Gene table.
# Inputs:      *.txt files in data/DisGeNET_hepato_cardio/
# Outputs:     One file per input in results/gene_disease_association/
# Author:      Federica D'Annunzio
# =============================================================================

rm(list = ls())

# --- Configuration -----------------------------------------------------------
path       <- "cell_lines/02_hepatosome"   # <-- update before running
input_dir  <- file.path(path, "data/DisGeNET_hepato_cardio")
output_dir <- file.path(path, "results/gene_disease_association")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Process each enrichment result file --------------------------------------
file_list <- list.files(path = input_dir, pattern = "*.txt", full.names = TRUE)

for (file_path in file_list) {

  data <- read.table(file_path, header = TRUE, quote = "", sep = "\t", check.names = FALSE)
  data <- data[, c(1, 9)]

  data$Genes <- gsub('^"|"$', '', data$Genes)
  gene_list  <- strsplit(data$Genes, ";")

  disease_gene_df <- data.frame(Disease = character(0), Gene = character(0))

  for (i in 1:nrow(data)) {
    disease <- data$Term[i]
    genes   <- gene_list[[i]]

    disease_gene_df <- rbind(
      disease_gene_df,
      data.frame(Disease = rep(disease, length(genes)), Gene = genes)
    )
  }

  output_file <- file.path(
    output_dir,
    paste0("disease_gene_df_", tools::file_path_sans_ext(basename(file_path)), ".txt")
  )

  write.table(disease_gene_df, output_file,
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}
