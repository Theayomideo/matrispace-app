# Round-two SHG comparisons

These sources extend the collagen validation with COL1A1, DCN, fibroblast density/fraction and fibrillar-collagen comparators, partial correlations, block-bootstrap intervals and the COL1A2 leave-one-out check. The saved results use 1,545 matched hexagons.

From this directory, where `work` is a writable analysis directory:

```sh
Rscript 01_prepare_cosmx_scores.R work
Rscript 06_hexagonal_validation.R work
Rscript compare_shg_baselines.R work results
Rscript check_col1a2_overlap.R work results
```

Required inputs under `work/input/` are listed with SHA-256 checksums in `input_checksums.csv`:

- `cosmx_section04/Run1129_Section_04_exprMat_file.csv` and `Run1129_Section_04_metadata_file.csv`: original [Zenodo dataset](https://doi.org/10.5281/zenodo.7899173).
- `shg/Collagen.tif`: the corresponding SHG image from that dataset.
- `shg/section4_aligned_metadata.csv`: deposited alignment from [the analysis archive](https://doi.org/10.5281/zenodo.22231956).
- `shg/coordinates.csv`: the `xc` and `yc` columns of the aligned metadata, in the same row order. The preparation script verifies agreement.
- `matrispace_0.3.0/matrisome_signatures.rda`: copy the frozen signature snapshot from this directory's `data/`.

The first two scripts generate the `derived/` and `results/` inputs consumed by the comparison scripts. Checksums of the original derived inputs are also recorded. Do not mix hexagon files from the older notebook with these area-matched results.

Full regeneration requires the analysis version **MatriSpace 0.3.0**, with `score_matrisome(signatures=..., min_genes=..., maxRank=...)` and `get_matrisome_scores()`. The currently published 0.1.0 package does not provide that API. Its release update is pending; this folder archives the analysis sources without replacing the package. The COL1A2 check additionally asserts UCell 2.14.0, as in the saved analysis.
