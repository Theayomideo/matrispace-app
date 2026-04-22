# Breast cancer case study

Worked example accompanying the MatriSpace manuscript: downstream analysis of a MatriSpace-profiled breast cancer Visium dataset (Janesick et al., 2023).

The script consumes MatriSpace output (per-spot matrisome/niche scores exported as CSV) rather than calling MatriSpace functions directly. It reproduces the figures in the case-study section of the manuscript.

## Inputs

| File | Source |
|------|--------|
| `visium.seurat.breast.rds` | Janesick et al. 2023 ([10x Visium preview dataset](https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast/)). Also available in MatriSpace online as *"Ductal carcinoma in situ, Invasive carcinoma (Breast) 1 10x"*. |
| `MatriSpace_Metadata_Breast_cancer_Ductal_Carcinoma_In_Situ_Invasive.csv` | Export from MatriSpace (Export → Metadata CSV) on the same sample. |
| `Curate.txt` | Curated TF–target table (same source as used in MatriCom). |

Place all three alongside `case_study.R` before running.

## Run

```r
setwd("analyses/breast-case-study")
source("case_study.R")
```

## Dependencies

`dplyr`, `data.table`, `ggplot2`, `Seurat`, `corrplot`, `scales`, `UCell`, `FNN`, `ggsci`, `pheatmap`, `effects`.
