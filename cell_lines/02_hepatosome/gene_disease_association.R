# =============================================================================
# Title:       Gene–Disease Association Table Builder
# Description: Reads DisGeNET enrichment result files (one per analysis), extracts
#              the Term (disease name) and Genes columns, splits the
#              semicolon-delimited gene lists, and reshapes each file into a
#              long-format Disease–Gene table that is saved as a tab-separated
#              text file.
# Inputs:      *.txt files in:
#              ~/Documents/HCC_data/cell_lines/functional_enrichment/enrichR/
#              DisGeNET_hepato_cardio/
#              Each file must be tab-separated with a header row and contain at
#              least columns "Term" (col 1) and "Genes" (col 9).
# Outputs:     One file per input, written to:
#              <path>/Res_gene_disease_association_hepato_cardio/
#              Named: disease_gene_df_<input_basename>.txt
# Author:      Federica D'Annunzio
# =============================================================================

rm(list = ls())

# ---------------------------------------------------------------------------
# Paths and output directory setup
# ---------------------------------------------------------------------------

path <- "~/Documents/HCC_data/cell_lines/functional_enrichment/enrichR/DisGeNET_hepato_cardio"
setwd(path)

# Create the results subdirectory if it does not already exist
if (!dir.exists(paste0(path, "/Res_gene_disease_association_hepato_cardio"))) {
  dir.create(paste0(path, "/Res_gene_disease_association_hepato_cardio"))
} else {
  print("The directory already exists")
}

output_dir <- paste0(path, "/Res_gene_disease_association_hepato_cardio")

# ---------------------------------------------------------------------------
# Process each enrichment result file
# ---------------------------------------------------------------------------

file_list <- list.files(path = path, pattern = "*.txt", full.names = TRUE)

for (file_path in file_list) {

  # Read the enrichment result table
  data <- read.table(file_path, header = TRUE, quote = "", sep = "\t", check.names = FALSE)

  # Retain only the Term (disease name) and Genes columns
  data <- data[, c(1, 9)]

  # Strip surrounding quotes from the Genes column if present
  data$Genes <- gsub('^"|"$', '', data$Genes)

  # Split each semicolon-delimited gene string into a list of gene symbols
  gene_list <- strsplit(data$Genes, ";")

  # Initialise an empty long-format data frame
  disease_gene_df <- data.frame(Disease = character(0), Gene = character(0))

  # Expand each disease row into one row per associated gene
  for (i in 1:nrow(data)) {
    disease <- data$Term[i]   # Disease term for this row
    genes   <- gene_list[[i]] # Vector of associated gene symbols

    disease_gene_df <- rbind(
      disease_gene_df,
      data.frame(
        Disease = rep(disease, length(genes)),
        Gene    = genes
      )
    )
  }

  # Build output filename from the input basename (without extension)
  output_file <- paste0(
    output_dir, "/disease_gene_df_",
    tools::file_path_sans_ext(basename(file_path)),
    ".txt"
  )

  # Save the long-format table as a tab-separated file
  write.table(
    disease_gene_df,
    output_file,
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}
