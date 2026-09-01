# UCell score benchmark

This notebook compares MatriSpace UCell scores with simple raw-count and SCT-normalized expression sums for the basement membrane and interstitial ECM signatures. It also tests the effect of multiplying each of the three highest-contributing genes by 2, 5, or 10.

## Input

Download [`Breast_cancer_Ductal_Carcinoma_In_Situ_Invasive.RDS`](https://zenodo.org/records/22231956/files/Breast_cancer_Ductal_Carcinoma_In_Situ_Invasive.RDS?download=1) from [Zenodo (doi:10.5281/zenodo.22231956)](https://doi.org/10.5281/zenodo.22231956). Place it in `data/`, or provide its path through the notebook parameter `data_file`. The object is the 4,992-spot breast cancer Visium sample used in the MatriSpace manuscript.

The 1.10 GB RDS is not stored in Git. Its exact size and SHA-256 checksum are recorded in `data_manifest.csv`.

The basement and interstitial signatures are read from the `matrispace` package. Scores are generated with `score_ecm_niche_signatures()` and retrieved with `get_ecm_niche_signature_scores()`.

## Render

```r
knitr::knit("ucell_score_benchmark.Rmd", output = "ucell_score_benchmark.md")
```

Required packages: `ggplot2`, `matrispace`, `Matrix`, `Seurat`, `SeuratObject`, `UCell`, and `knitr`.
