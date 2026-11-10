# =============================================================================
# Title:       Extract UBIAD1 Interactome from PPI Network
# Description: Subsets the full PPI interactome to retain only interactions
#              involving UBIAD1. Converts Entrez gene IDs to HGNC symbols
#              using biomaRt and exports both edge lists.
# Inputs:      - data/interactome_PNAS2021.txt: full PPI edge list
# Outputs:     - results/ubiad1_interactome.txt         (Entrez IDs)
#              - results/ubiad1_interactome_symbols.txt  (HGNC symbols)
# Author:      Federica D'Annunzio
# =============================================================================

rm(list = ls())

library(dplyr)
library(biomaRt)

# --- Paths (relative to repo root) -------------------------------------------
base_path <- "cell_lines/05_target_analysis"   # <-- update before running

# --- Load interactome ---------------------------------------------------------
interactome <- read.table(
  file.path(base_path, "data/interactome_PNAS2021.txt"),
  sep = "\t", header = TRUE, quote = "", check.names = FALSE
)

# --- Subset to UBIAD1 interactions --------------------------------------------
ubiad1_id <- "29914"

subset_ubiad1 <- interactome[
  which(interactome$source %in% ubiad1_id | interactome$target %in% ubiad1_id),
]

subset_ubiad1 <- subset_ubiad1 %>% select(-databases)

# --- Export Entrez-ID-based edge list -----------------------------------------
res_dir <- file.path(base_path, "results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

write.table(subset_ubiad1, file.path(res_dir, "ubiad1_interactome.txt"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

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

subset_ubiad1_conv <- subset_ubiad1 %>%
  left_join(conversion, by = c("source" = "entrezgene_id")) %>%
  rename(source_symbol = hgnc_symbol) %>%
  left_join(conversion, by = c("target" = "entrezgene_id")) %>%
  rename(target_symbol = hgnc_symbol)

subset_ubiad1_symbols <- subset_ubiad1_conv %>%
  dplyr::select(source_symbol, target_symbol)

write.table(subset_ubiad1_symbols, file.path(res_dir, "ubiad1_interactome_symbols.txt"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
