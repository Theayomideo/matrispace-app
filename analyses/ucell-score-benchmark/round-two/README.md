# Round-two UCell sensitivity analyses

`run_ucell_benchmark.R` is the source used for the extended benchmark, including maxRank sensitivity. `results/` preserves its compact tables and per-spot scores. `verify_maxrank.R` independently reconstructs Spearman correlations and top-decile Jaccard overlap from the saved scores (4,992 spots).

From this directory:

```sh
Rscript verify_maxrank.R results
Rscript run_ucell_benchmark.R --input=/path/to/Breast_cancer_Ductal_Carcinoma_In_Situ_Invasive.RDS --signatures=../../../inst/extdata/ecm_ucell_signatures.rds --hgnc=../../../inst/extdata/hgnc_long_data.rds --output=rerun --strength-tests
```

When running from a repository checkout, the bundled references are at `../../../inst/extdata` relative to `analyses/ucell-score-benchmark/round-two`. The Janesick input is deposited at DOI [10.5281/zenodo.22231956](https://doi.org/10.5281/zenodo.22231956). `input_checksums.csv` identifies the exact input and frozen reference files used for the saved run. The original reference snapshots have the same content as the app's corresponding bundled HGNC/niche files.

The script lists and checks its R dependencies. This extended runner uses UCell directly; the earlier notebook in the parent directory uses MatriSpace 0.3.0. The public MatriSpace package update is pending and is not replaced here. `results/run_metadata.csv` records the original run environment. `plot_strength_tests.R` rebuilds plots from `results/`.
