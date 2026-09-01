# SHG–CosMx collagen validation

This notebook compares MatriSpace collagen scores from a CosMx-profiled lung adenocarcinoma section with second-harmonic generation (SHG) collagen imaging from the adjacent section described by Pentimalli et al. (2025).

## Input

Place the following files under `data/`, preserving the subdirectories shown below, or provide the input directory through the notebook parameter `data_dir`.

```text
data/
├── cosmx_section04/
│   └── Run1129_Section_04_exprMat_file.csv
└── shg/
    ├── Collagen.tif
    └── section4_aligned_metadata.csv
```

Sources and checksums are recorded in `data_manifest.csv`. The collagen signature is read from the `matrispace` package. Scores are generated with `score_matrisome()` and retrieved with `get_matrisome_scores()`.

Download [`section4_aligned_metadata.csv`](https://zenodo.org/records/22231956/files/section4_aligned_metadata.csv?download=1) from [Zenodo (doi:10.5281/zenodo.22231956)](https://doi.org/10.5281/zenodo.22231956). The raw CosMx matrix and SHG image are available from the source dataset ([doi:10.5281/zenodo.7899173](https://doi.org/10.5281/zenodo.7899173)).

## Render

```r
knitr::knit("shg_cosmx_validation.Rmd", output = "shg_cosmx_validation.md")
```

Required packages: `data.table`, `ggplot2`, `hexbin`, `matrispace`, `Matrix`, `Seurat`, `SeuratObject`, `tiff`, `UCell`, and `knitr`.
