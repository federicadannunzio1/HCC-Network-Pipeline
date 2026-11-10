# =============================================================================
# Title:       Find Nearest-Neighbour Genes of Selected Switch Genes
# Description: For each dataset, reads the nearest-neighbour table of switch
#              genes, filters it to the selected switch genes, resolves the
#              union of their neighbours, removes neighbours that are
#              themselves switch genes, annotates them with expression direction
#              (UP/DOWN based on logFC), and appends the result as a new sheet
#              in an output Excel workbook.
# Inputs:      For each dataset directory under `data_path`:
#              - <data_path>/<dataset>/switch/txtFile/nn_<tag>_switch.txt
#              - <data_path>/<dataset>/switch/txtFile/attribute.txt
#              - data/selected_switch/<dataset>.txt
# Outputs:     results/attribute_nn_<tag>_selected_switch.xlsx
# Author:      Federica D'Annunzio
# =============================================================================

rm(list = ls())

library(openxlsx)

options(stringsAsFactors = FALSE)

# --- Configuration -----------------------------------------------------------
data_path <- "data/swim_output/"                # <-- update before running
tag       <- "pos"

datasets <- dir(data_path, full.names = FALSE)

# --- Initialise output workbook -----------------------------------------------
fileout_nn <- paste0("results/attribute_nn_", tag, "_selected_switch.xlsx")
wb_nn      <- createWorkbook()

# --- Process each dataset -----------------------------------------------------
lapply(datasets, function(d) {

  nn <- read.table(
    file.path(data_path, d, "switch/txtFile", paste0("nn_", tag, "_switch.txt")),
    row.names = 1, check.names = FALSE, header = TRUE, quote = "", sep = "\t"
  )

  attribute <- read.table(
    file.path(data_path, d, "switch/txtFile/attribute.txt"),
    row.names = 1, check.names = FALSE, header = TRUE, quote = "", sep = "\t"
  )

  switch_selected <- read.table(
    file.path("data/selected_switch", paste0(d, ".txt")),
    check.names = FALSE, header = FALSE, quote = ""
  )$V1

  all_switch <- rownames(nn)
  nn <- nn[switch_selected, ]

  nn <- unique(unlist(nn))
  ind <- which(nn == " ")
  nn <- nn[-ind]

  nn <- setdiff(nn, all_switch)

  attribute_nn <- attribute[nn, ]
  attribute_nn$direction <- ifelse(attribute_nn$logFC > 0, "UP", "DOWN")

  addWorksheet(wb = wb_nn, sheetName = d)
  writeData(wb = wb_nn, sheet = d, x = attribute_nn, rowNames = TRUE)
  saveWorkbook(wb_nn, fileout_nn, overwrite = TRUE)
})
