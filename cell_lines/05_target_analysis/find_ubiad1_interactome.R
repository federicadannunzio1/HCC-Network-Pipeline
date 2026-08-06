# =============================================================================
# Title:       Extract UBIAD1 Interactome from PPI Network
# Description: Subsets the full protein-protein interaction (PPI) interactome
#              to retain only interactions involving UBIAD1. Converts Entrez
#              gene IDs to HGNC symbols using biomaRt and exports both the
#              Entrez-based and symbol-based edge lists.
# Inputs:      - interactome_PNAS2021.txt: full PPI edge list with `source`,
#                  `target`, and `databases` columns (Entrez gene IDs)
# Outputs:     - ubiad1_interactome.txt        : UBIAD1 edges (Entrez IDs)
#              - ubiad1_interactome_symbols.txt : UBIAD1 edges (HGNC symbols)
# Author:      Federica D'Annunzio
# =============================================================================

rm(list = ls())

# --- Libraries ----------------------------------------------------------------
library(dplyr)
library(biomaRt)

# --- Paths --------------------------------------------------------------------
interactome_path <- "~/Documents/HCC_data/cell_lines/interactome_mapping/code/files/"
output_path      <- "~/Documents/HCC_data/cell_lines/genes_overview/brca-relb/"

# --- Load interactome ---------------------------------------------------------
interactome <- read.table(
  paste0(interactome_path, "interactome_PNAS2021.txt"),
  sep         = "\t",
  header      = TRUE,
  quote       = "",
  check.names = FALSE
)

# --- Subset to UBIAD1 interactions --------------------------------------------
ubiad1_id <- "29914"
brca1_id  <- "672"   # retained for reference; not used in subsetting

# Keep all edges where UBIAD1 is either source or target
subset_ubiad1 <- interactome[
  which(interactome$source %in% ubiad1_id | interactome$target %in% ubiad1_id),
]

# Remove the databases column (not needed for downstream analysis)
subset_ubiad1 <- subset_ubiad1 %>%
  select(-databases)

# Export Entrez-ID-based edge list
write.table(
  subset_ubiad1,
  paste0(output_path, "ubiad1_interactome.txt"),
  sep       = "\t",
  col.names = TRUE,
  row.names = FALSE,
  quote     = FALSE
)

# --- Convert Entrez IDs to HGNC symbols ---------------------------------------
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

all_entrez_ids <- unique(c(subset_ubiad1$source, subset_ubiad1$target))

conversion <- getBM(
  attributes = c("entrezgene_id", "hgnc_symbol"),
  filters    = "entrezgene_id",
  values     = all_entrez_ids,
  mart       = ensembl
)

conversion <- unique(conversion)

# Join HGNC symbols onto source and target columns
subset_ubiad1_conv <- subset_ubiad1 %>%
  left_join(conversion, by = c("source" = "entrezgene_id")) %>%
  rename(source_symbol = hgnc_symbol) %>%
  left_join(conversion, by = c("target" = "entrezgene_id")) %>%
  rename(target_symbol = hgnc_symbol)

# Retain only the symbol columns
subset_ubiad1_symbols <- subset_ubiad1_conv %>%
  dplyr::select(source_symbol, target_symbol)

# --- Export symbol-based edge list --------------------------------------------
write.table(
  subset_ubiad1_symbols,
  paste0(output_path, "ubiad1_interactome_symbols.txt"),
  sep       = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote     = FALSE
)
