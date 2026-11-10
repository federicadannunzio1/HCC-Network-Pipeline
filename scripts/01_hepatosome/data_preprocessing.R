# =============================================================================
# Data Preprocessing for Hypergeometric Enrichment
# =============================================================================
# Description: Converts Ensembl gene IDs to HGNC symbols using biomaRt and
#              prepares the universe gene list and DisGeNET property list
#              for downstream hypergeometric enrichment analysis.
#
# Input:
#   - data/matrix/matrix.txt                    : RNA-seq expression matrix
#   - data/hypergeom_data/*_switch_converted.txt : Switch gene lists per condition
#   - data/DisGeNET_cardio.txt                  : DisGeNET disease-gene associations
#
# Output:
#   - data/mat_universo.txt                     : Universe gene list (HGNC symbols)
#   - data/DisGeNET_proprietà_cardio.txt        : Unique DisGeNET genes
#
# Author: Federica D'Annunzio
# =============================================================================

rm(list = ls())

library(biomaRt)

# --- Configuration -----------------------------------------------------------
path <- "cell_lines/02_hepatosome"   # <-- update before running

# --- Load input files ---------------------------------------------------------
mat_universo <- read.table(
  file.path(path, "data/matrix/matrix.txt"),
  sep = "\t", quote = "", header = TRUE, check.names = FALSE, row.names = 1
)

switch_ATO5      <- read.table(file.path(path, "data/hypergeom_data/ATO5_vs_CTRL_switch_converted.txt"), quote = "")
switch_EZE25_ATO <- read.table(file.path(path, "data/hypergeom_data/EZE25_ATO_vs_CTRL_switch_converted.txt"), quote = "")
switch_EZE25     <- read.table(file.path(path, "data/hypergeom_data/EZE25_vs_CTRL_switch_converted.txt"), quote = "")

# --- Extract unique DisGeNET genes --------------------------------------------
disgenet <- read.table(
  file.path(path, "data/DisGeNET_cardio.txt"),
  header = TRUE, sep = "\t", fill = TRUE, quote = "", check.names = FALSE
)

all_genes    <- unlist(disgenet[, -1])
all_genes    <- all_genes[!is.na(all_genes) & all_genes != ""]
unique_genes <- unique(all_genes)
disgenet_out <- data.frame(Gene = unique_genes)

write.table(disgenet_out, file.path(path, "data/DisGeNET_proprietà_cardio.txt"),
            sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)

# --- Convert Ensembl IDs to HGNC symbols (universe) ---------------------------
ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

ensembl_ids <- gsub("\\.\\d+$", "", rownames(mat_universo))

hgnc_symbols <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters    = "ensembl_gene_id",
  values     = ensembl_ids,
  mart       = ensembl
)$hgnc_symbol

hgnc_symbols <- unique(hgnc_symbols)
hgnc_symbols <- hgnc_symbols[hgnc_symbols != ""]

mat_universo_out <- as.data.frame(hgnc_symbols)

write.table(mat_universo_out, file.path(path, "data/mat_universo.txt"),
            sep = "\t", col.names = TRUE, quote = FALSE, row.names = FALSE)
