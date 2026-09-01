# Reproducible analyses

This directory contains analyses accompanying the MatriSpace manuscript.

- [`breast-case-study`](breast-case-study): breast cancer case study.
- [`ucell-score-benchmark`](ucell-score-benchmark): comparison of UCell scores with raw and normalized expression sums.
- [`shg-cosmx-validation`](shg-cosmx-validation): adjacent-section comparison of CosMx collagen scores with SHG collagen imaging.

Each analysis directory contains an R Markdown notebook, its rendered Markdown output, individual figures, summary results, and an input manifest. Large source data are not stored in Git.

Large inputs generated for the UCell benchmark and SHG–CosMx validation are deposited at [Zenodo (doi:10.5281/zenodo.22231956)](https://doi.org/10.5281/zenodo.22231956). Original public inputs are linked in each analysis-specific `data_manifest.csv`.
