# HCC Network Pipeline

**Network medicine framework for the identification of novel drug targets in hepatocellular carcinoma (HCC) through lipid-lowering therapy.**

> 📄 *Manuscript in preparation*
> D'Annunzio F., Pedrelli M., Fiscon G., Paci P., Parini P.
> *From hepatosome to interactome: network-based discovery of UBIAD1 as a lipid-metabolism target in hepatocellular carcinoma*

---

## Overview

This repository contains the R-based analysis pipeline developed to identify **UBIAD1** as a novel putative drug target in HCC, by integrating RNA-seq transcriptomics with a network medicine approach.

The pipeline was applied to two datasets:
- **Cell lines** — SOAT2-only HepG2 cells treated with atorvastatin (ATO5, 5 µM), ezetimibe (EZE25, 25 µM), their combination, or placebo
- **Human tissues** — liver biopsies from the TSS clinical trial (simvastatin + ezetimibe vs placebo)

---

## Pipeline

```
RNA-seq data (SOAT2-only HepG2 cells)
         │
         ▼
  01 · Switch Gene Detection
         SWIM-based co-expression network analysis
         Identifies genes whose network role shifts between conditions
         │
         ▼
  02 · Hepatosome Construction
         Hypergeometric enrichment of switch genes in DisGeNET
         Bipartite network: switch genes ↔ hepato-biliary diseases
         │
         ▼
  03 · Interactome Mapping
         Mapping switch genes onto the human PPI network
         LCC (Largest Connected Component) analysis
         Statistical validation against random gene sets
         │
         ▼
  04 · Random Walk with Restart (RWR)
         Propagation from 8 atorvastatin/ezetimibe drug targets
         Degree-preserving randomization (1 000 iterations)
         Sensitivity analysis: leave-HMGCR-out test
         │
         ▼
  05 · Target Analysis
         UBIAD1 characterisation in the interactome and hepatosome
         Expression analysis: UBIAD1 vs HMGCR across conditions
         │
         ▼
  Human tissue validation
         UBIAD1 & HMGCR mRNA in liver biopsies (TSS trial)
```

---

## Repository Structure

```
HCC-Network-Pipeline/
│
├── scripts/
│   ├── 02_hepatosome/             # Hypergeometric enrichment + disease network
│   ├── 03_interactome_mapping/    # LCC analysis in the human interactome
│   │   └── scripts/               #   Helper functions for module analysis
│   ├── 04_random_walk/            # RWR, degree-preserving randomization,
│   │                              #   sensitivity analysis & plots
│   ├── 05_target_analysis/        # UBIAD1 interactome + expression plots
│   └── 06_tissue_validation/      # UBIAD1/HMGCR in liver biopsies (TSS trial)
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
| Statistical validation | Degree-preserving randomization (1 000 iterations) | This work |
| Sensitivity analysis | Leave-one-seed-out (HMGCR) | This work |
| Gene ID mapping | org.Hs.eg.db | Carlson M., Bioconductor |

> **Note on SWIM:** Switch gene detection was performed using the SWIM method (Paci P. et al., *Sci Rep* 2017) implemented via the SWIMmeR R package (Paci P. & Fiscon G., *Bioinformatics* 2022). SWIM code is not included in this repository; please refer to the [SWIMmeR package](https://github.com/sportingCode/SWIMmeR) for installation and usage.

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

---

## Data Availability

RNA-seq data from the TSS clinical trial are available upon reasonable request to the corresponding authors. Cell line RNA-seq data will be made available upon publication.

---

## Contact

**Federica D'Annunzio**
Dept. of Translational and Precision Medicine, Sapienza University of Rome
✉ federica.dannunzio@uniroma1.it
