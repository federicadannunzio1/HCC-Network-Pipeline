# HCC Network Pipeline

Network medicine framework for the identification of novel drug targets in hepatocellular carcinoma (HCC).

---

## Pipeline

```
RNA-seq data
         │
         ▼
  Switch Gene Detection (SWIM — see SWIMmeR repository)
         SWIM-based co-expression network analysis
         Identifies genes whose network role shifts between conditions
         │
         ▼
  01 · Hepatosome Construction
         Hypergeometric enrichment of switch genes in DisGeNET
         Bipartite network: switch genes ↔ hepato-biliary diseases
         │
         ▼
  02 · Interactome Mapping
         Mapping switch genes onto the human PPI network
         LCC (Largest Connected Component) analysis
         Statistical validation against random gene sets
         │
         ▼
  03 · Random Walk with Restart (RWR)
         Propagation from drug-target seed genes
         Degree-preserving randomization (1 000 iterations)
         Sensitivity analysis: leave-one-seed-out test
         │
         ▼
  04 · Target Analysis
         Target characterisation in the interactome
         Expression analysis across treatment conditions
         │
         ▼
  05 · Tissue Validation
         mRNA expression in human liver biopsies
```

---

## Repository Structure

```
HCC-Network-Pipeline/
│
├── scripts/
│   ├── 01_hepatosome/             # Hypergeometric enrichment + disease network
│   ├── 02_interactome_mapping/    # LCC analysis in the human interactome
│   │   └── scripts/               #   Helper functions for module analysis
│   ├── 03_random_walk/            # RWR, degree-preserving randomization,
│   │                              #   sensitivity analysis & plots
│   ├── 04_target_analysis/        # Target interactome + expression plots
│   └── 05_tissue_validation/      # Expression validation in liver biopsies
│
├── .gitignore
├── LICENSE
└── README.md
```

---

## Methods & Key Tools

| Step | Method | Reference |
|------|--------|-----------|
| Switch gene detection | SWIM | Paci P. et al., *Sci Rep*, 2017 |
| Disease enrichment | Hypergeometric test (DisGeNET) | This work |
| Interactome module | LCC analysis | This work |
| Target prioritisation | Random Walk with Restart | Köhler S. et al., *Am J Hum Genet*, 2008 |
| Statistical validation | Degree-preserving randomization (1 000 iter.) | This work |
| Sensitivity analysis | Leave-one-seed-out | This work |

> **Note on SWIM:** Switch gene detection was performed using the SWIM method (Paci P. et al., *Sci Rep* 2017) implemented via the SWIMmeR R package (Paci P. & Fiscon G., *Bioinformatics* 2022). SWIM code is not included in this repository; please refer to the [SWIMmeR package](https://github.com/sportingCode/SWIMmeR).

---

## Requirements

```r
# R >= 4.0
install.packages(c(
  "igraph",       # network analysis
  "ggplot2",      # plotting
  "tidyr",        # data reshaping
  "RColorBrewer"  # colour palettes
))

# Bioconductor
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c(
  "org.Hs.eg.db",           # Entrez ↔ Symbol mapping
  "RandomWalkRestartMH"     # multiplex RWR
))
```
