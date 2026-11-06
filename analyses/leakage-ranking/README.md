# Leakage-ranking comparison

Reproduces rebuttal Table R1 using the app's default correction (alpha 0.10, six neighbours, distance cutoff 1.6) and the Janesick breast cancer section. The six profiles sum SCT counts over canonical category genes. Spearman correlations are unchanged by the profile's monotonic scaling.

Run with the updated `matrispace.app` installed:

```sh
Rscript run_leakage_ranking.R /path/to/Breast_cancer_Ductal_Carcinoma_In_Situ_Invasive.RDS results
```

The input is deposited at DOI [10.5281/zenodo.22231956](https://doi.org/10.5281/zenodo.22231956); its SHA-256 is recorded in `../ucell-score-benchmark/data_manifest.csv`.

This is a reconstruction using the harmonized canonical lists, not the recovered original script. All six correlations reproduce Table R1 at its four-decimal precision: 0.9993, 0.9992, 0.9991, 0.9992, 0.9991 and 0.9986 in the saved table's order. The input has 4,992 spots and 990 detected canonical genes. F8 and USH2A remain available to functional signatures but do not contribute to these six profiles.

Saved results, input MD5, parameters and session versions are in `results/`. The comparison measures sensitivity to this correction; it does not independently estimate contamination.
