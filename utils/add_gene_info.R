# =============================================================================
# Title:       Add Gene Information via Biomart
# Description: Reads a list of HGNC gene symbols from an Excel file, queries
#              Ensembl Biomart to retrieve Ensembl IDs, synonyms, and
#              descriptions for each gene, and writes the annotated table to
#              an output Excel file.
# Inputs:      ~/Desktop/genes_in_manuscript.xlsx
#              - Sheet 1, no header; gene symbols are expected in column X1.
# Outputs:     ~/Documents/HCC_data/cell_lines/paper_draft/for_publication/
#              supplementary/supplementary_tables/Table_S2.xlsx
# Author:      Federica D'Annunzio
# =============================================================================

rm(list = ls())

library(biomaRt)
library(xlsx)

# ---------------------------------------------------------------------------
# Load gene list
# ---------------------------------------------------------------------------

# Read gene symbols from the first sheet of the input Excel file
file <- read.xlsx("~/Desktop/genes_in_manuscript.xlsx", sheetIndex = 1, header = F)

# Clean whitespace and drop any empty entries
genes <- trimws(as.character(file$X1))
genes <- genes[genes != ""]

# ---------------------------------------------------------------------------
# Query Biomart for gene annotations
# ---------------------------------------------------------------------------

# Connect to the human Ensembl dataset
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve Ensembl ID, HGNC symbol, synonyms, and gene description
annot <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "external_synonym", "description"),
  filters    = "hgnc_symbol",
  values     = genes,
  mart       = mart
)

# ---------------------------------------------------------------------------
# Export annotated table
# ---------------------------------------------------------------------------

# Write annotation table to the supplementary Excel output file
#write.table(annot, "~/Documents/HCC_data/cell_lines/attribute/selected_switch_combo_gene_names.txt", sep = "\t", quote = F, row.names = F)
write.xlsx(
  annot,
  "~/Documents/HCC_data/cell_lines/paper_draft/for_publication/supplementary/supplementary_tables/Table_S2.xlsx",
  row.names = FALSE
)
