# =============================================================================
# Title:       ENSEMBL to Gene Symbol Conversion
# Description: Reads lists of ENSEMBL gene IDs from input files, converts them
#              to HGNC gene symbols using Biomart, removes empty/duplicate
#              entries, and writes the converted lists to output text files.
# Inputs:      - filenames_input: character vector of input filenames
#              - dir_input_name:  path to the directory containing input files
#              Each input file is tab-delimited with ENSEMBL IDs in column V1.
# Outputs:     One converted text file per input file, written to:
#              ~/Documents/HCC_data/cell_lines/specific_selected_switch/
#              with the suffix "_converted.txt".
# Author:      Federica D'Annunzio
# =============================================================================

# Connect to Ensembl BioMart and retrieve available attributes and filters
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl)
filters <- listFilters(ensembl)

# Process each input file: convert ENSEMBL IDs to HGNC symbols and save output
lapply(filenames_input, function(x) {

  # Read ENSEMBL gene IDs from the input file
  input_list <- read.table(
    paste0(dir_input_name, x),
    header = F,
    sep = ".",
    check.names = F,
    quote = ""
  )$V1

  # Convert ENSEMBL IDs to HGNC gene symbols via Biomart
  input_list <- getBM(
    attributes = c("hgnc_symbol", "ensembl_gene_id"),
    filters    = "ensembl_gene_id",
    values     = input_list,
    mart       = ensembl
  )$hgnc_symbol

  # Remove duplicates and entries with empty gene symbols
  input_list <- unique(input_list)
  ind <- which(input_list == "")
  if (length(ind) > 0) input_list <- input_list[-ind]

  # Build output file path, stripping the original file extension
  dir_output_name <- "~/Documents/HCC_data/cell_lines/specific_selected_switch"
  output_file <- paste0(dir_output_name, sub("\\..*", "", x), "_converted.txt")

  # Write the converted gene symbols to a plain text file
  write.table(
    input_list,
    file      = output_file,
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

})
