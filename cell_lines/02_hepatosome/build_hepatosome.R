# =============================================================================
# Build Hepatosome Network
# =============================================================================
# Description: Constructs the hepatosome — a bipartite disease-gene network —
#              from DisGeNET data. Each edge connects a hepato-biliary disease
#              to an associated gene. Also computes disease node degrees.
#
# Input:
#   - DisGeNET_table.xlsx : DisGeNET table with columns "Term" and "Genes"
#                           (genes separated by ";")
#
# Output:
#   - hepatosome.txt      : Edge list (source = disease, target = gene)
#   - degree_diseases.txt : Disease node degree (number of associated genes)
#   - diseases.txt        : List of all disease nodes
#
# Author: Federica D'Annunzio
# =============================================================================

rm(list = ls())

library(openxlsx)
library(data.table)

path <- "path/to/02_hepatosome"   # <-- update before running
setwd(path)

# -----------------------------------------------------------------------------
# Load DisGeNET disease-gene table
# -----------------------------------------------------------------------------
data <- read.xlsx("DisGeNET_table.xlsx",
                  colNames = TRUE, check.names = FALSE, sheet = 1)

# Split gene lists (semicolon-separated) into individual genes per disease
gene_list        <- strsplit(data$Genes, ";")
names(gene_list) <- data$Term

# -----------------------------------------------------------------------------
# Build edge list (disease → gene)
# -----------------------------------------------------------------------------
edge_list <- rbindlist(lapply(names(gene_list), function(disease) {
  data.frame(source = disease, target = gene_list[[disease]])
}))

# -----------------------------------------------------------------------------
# Compute disease degree (number of associated genes)
# -----------------------------------------------------------------------------
disease_names  <- names(gene_list)
disease_degree <- sapply(gene_list, length)
degree_df      <- data.frame(disease = disease_names, degree = disease_degree)

# -----------------------------------------------------------------------------
# Save outputs
# -----------------------------------------------------------------------------
write.table(edge_list,    "hepatosome.txt",      sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(degree_df,    "degree_diseases.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(disease_names, "diseases.txt",       sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
