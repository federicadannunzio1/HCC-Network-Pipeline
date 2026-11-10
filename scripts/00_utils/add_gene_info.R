# =============================================================================
# Title:       Add Gene Information via Biomart
# Description: Reads a list of HGNC gene symbols from an Excel file, queries
#              Ensembl Biomart to retrieve Ensembl IDs, synonyms, and
#              descriptions for each gene, and writes the annotated table to
#              an output Excel file.
# Inputs:      - genes_in_manuscript.xlsx (Sheet 1, no header; symbols in col X1)
# Outputs:     - Table_S4.xlsx (gene annotation table)
# Author:      Federica D'Annunzio
# =============================================================================

rm(list = ls())

library(biomaRt)
library(openxlsx)

# --- Configuration -----------------------------------------------------------
input_file  <- "data/genes_in_manuscript.xlsx"    # <-- update before running
output_file <- "results/Table_S4.xlsx"            # <-- update before running

# --- Load gene list -----------------------------------------------------------
file <- read.xlsx(input_file, sheet = 1, colNames = FALSE)

genes <- trimws(as.character(file$X1))
genes <- genes[genes != ""]

# --- Query Biomart for gene annotations ---------------------------------------
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

annot <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "external_synonym", "description"),
  filters    = "hgnc_symbol",
  values     = genes,
  mart       = mart
)

# --- Export annotated table ---------------------------------------------------
write.xlsx(annot, output_file, rowNames = FALSE)
