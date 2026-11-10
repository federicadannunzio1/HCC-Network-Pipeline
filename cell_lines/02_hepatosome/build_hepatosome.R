# =============================================================================
# Build Hepatosome Network
# =============================================================================
# Description: Constructs the hepatosome — a bipartite disease-gene network —
#              from DisGeNET data. Each edge connects a hepato-biliary disease
#              to an associated gene. Also computes disease node degrees.
#
# Input:
#   - data/DisGeNET_table.xlsx : DisGeNET table with "Term" and "Genes" columns
#
# Output:
#   - results/hepatosome.txt       : Edge list (source = disease, target = gene)
#   - results/degree_diseases.txt  : Disease node degree
#   - results/diseases.txt         : List of all disease nodes
#
# Author: Federica D'Annunzio
# =============================================================================

rm(list = ls())

library(openxlsx)
library(data.table)

# --- Configuration -----------------------------------------------------------
path <- "cell_lines/02_hepatosome"   # <-- update before running

# --- Load DisGeNET disease-gene table -----------------------------------------
data <- read.xlsx(file.path(path, "data/DisGeNET_table.xlsx"),
                  colNames = TRUE, check.names = FALSE, sheet = 1)

gene_list        <- strsplit(data$Genes, ";")
names(gene_list) <- data$Term

# --- Build edge list (disease -> gene) ----------------------------------------
edge_list <- rbindlist(lapply(names(gene_list), function(disease) {
  data.frame(source = disease, target = gene_list[[disease]])
}))

# --- Compute disease degree ---------------------------------------------------
disease_names  <- names(gene_list)
disease_degree <- sapply(gene_list, length)
degree_df      <- data.frame(disease = disease_names, degree = disease_degree)

# --- Save outputs -------------------------------------------------------------
res_dir <- file.path(path, "results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

write.table(edge_list,     file.path(res_dir, "hepatosome.txt"),      sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(degree_df,     file.path(res_dir, "degree_diseases.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(disease_names, file.path(res_dir, "diseases.txt"),        sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
