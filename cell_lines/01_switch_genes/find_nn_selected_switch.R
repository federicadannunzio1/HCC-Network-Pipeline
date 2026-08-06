# =============================================================================
# Title:       Find Nearest-Neighbour Genes of Selected Switch Genes
# Description: For each dataset, reads the nearest-neighbour table of switch
#              genes, filters it to the selected switch genes, resolves the
#              union of their neighbours, removes neighbours that are
#              themselves switch genes, annotates them with expression direction
#              (UP/DOWN based on logFC), and appends the result as a new sheet
#              in an output Excel workbook.
# Inputs:      For each dataset directory under `path`:
#              - <path>/<dataset>/switch/txtFile/nn_<tag>_switch.txt
#                Nearest-neighbour table (rows = switch genes, columns = neighbours).
#              - <path>/<dataset>/switch/txtFile/attribute.txt
#                Attribute table with at least a logFC column.
#              - ../selected_switch/<dataset>.txt
#                List of selected switch gene IDs (one per row, no header).
# Outputs:     ../attribute/attribute_nn_<tag>_selected_switch.xlsx
#              One sheet per dataset containing the attribute table of
#              nearest-neighbour genes with an added "direction" column.
# Author:      Federica D'Annunzio
# =============================================================================

rm(list = ls())

library(openxlsx)

options(stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

setwd("G:/Drive condivisi/Giulia-Fede/Collaborazioni/PaoloParini/SWIMmeR/utils")

# Root path containing one subdirectory per dataset
path <- "G:/Drive condivisi/Giulia-Fede/SWIMmeR/code/project/PaoloParini/dataset/"

# Tag identifying the switch type to process (e.g. "pos" or "neg")
tag <- "pos"

# Discover all dataset directories
datasets <- dir(path, ".")

# ---------------------------------------------------------------------------
# Initialise output workbook
# ---------------------------------------------------------------------------

fileout_nn <- paste0("../attribute/attribute_nn_", tag, "_selected_switch.xlsx")
wb_nn <- createWorkbook()

# ---------------------------------------------------------------------------
# Process each dataset
# ---------------------------------------------------------------------------

lapply(datasets, function(d) {

  # Load nearest-neighbour table for switch genes
  nn <- read.table(
    paste0(path, d, "/switch/txtFile/nn_", tag, "_switch.txt"),
    row.names   = 1,
    check.names = FALSE,
    header      = TRUE,
    quote       = "",
    sep         = "\t"
  )

  # Load full attribute table for all genes in this dataset
  attribute <- read.table(
    paste0(path, d, "/switch/txtFile/attribute.txt"),
    row.names   = 1,
    check.names = FALSE,
    header      = TRUE,
    quote       = "",
    sep         = "\t"
  )

  # Load the list of selected switch genes for this dataset
  switch_selected <- read.table(
    paste0("../selected_switch/", d, ".txt"),
    check.names = FALSE,
    header      = FALSE,
    quote       = ""
  )$V1

  # Record all switch gene IDs (used later to exclude them from neighbours)
  all_switch <- rownames(nn)

  # Subset the nearest-neighbour table to the selected switch genes
  nn <- nn[switch_selected, ]

  # Flatten all neighbour entries into a unique vector and remove blank entries
  nn <- unique(unlist(nn))
  ind <- which(nn == " ")
  nn <- nn[-ind]

  # Exclude genes that are themselves switch genes
  nn <- setdiff(nn, all_switch)

  # Retrieve attribute rows for the identified nearest-neighbour genes
  attribute_nn <- attribute[nn, ]

  # Annotate expression direction based on logFC sign
  attribute_nn$direction <- ifelse(attribute_nn$logFC > 0, "UP", "DOWN")

  # Write results as a new sheet in the output workbook
  addWorksheet(wb = wb_nn, sheetName = d)
  writeData(wb = wb_nn, sheet = d, x = attribute_nn, rowNames = TRUE)
  saveWorkbook(wb_nn, fileout_nn, overwrite = TRUE)

})
