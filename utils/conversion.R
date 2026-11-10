# =============================================================================
# Title:       ENSEMBL to Gene Symbol Conversion
# Description: Reads lists of ENSEMBL gene IDs from input files, converts them
#              to HGNC gene symbols using Biomart, removes empty/duplicate
#              entries, and writes the converted lists to output text files.
# Inputs:      - Tab-delimited files with ENSEMBL IDs in column V1
# Outputs:     - One converted text file per input, with suffix "_converted.txt"
# Author:      Federica D'Annunzio
# =============================================================================

rm(list = ls())

library(biomaRt)

# --- Configuration -----------------------------------------------------------
# Set these before running:
dir_input_name  <- "data/input/"          # directory containing input files
dir_output_name <- "data/output/"         # directory for converted output
filenames_input <- list.files(dir_input_name, pattern = "\\.txt$")

# --- Connect to Ensembl BioMart ----------------------------------------------
ensembl    <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl)
filters    <- listFilters(ensembl)

# --- Process each input file --------------------------------------------------
lapply(filenames_input, function(x) {

  input_list <- read.table(
    paste0(dir_input_name, x),
    header = FALSE, sep = ".", check.names = FALSE, quote = ""
  )$V1

  input_list <- getBM(
    attributes = c("hgnc_symbol", "ensembl_gene_id"),
    filters    = "ensembl_gene_id",
    values     = input_list,
    mart       = ensembl
  )$hgnc_symbol

  input_list <- unique(input_list)
  ind <- which(input_list == "")
  if (length(ind) > 0) input_list <- input_list[-ind]

  output_file <- paste0(dir_output_name, sub("\\..*", "", x), "_converted.txt")

  write.table(input_list, file = output_file,
              quote = FALSE, row.names = FALSE, col.names = FALSE)
})
