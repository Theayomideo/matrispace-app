# Matrisome references

The harmonized reference contains 1,027 canonical genes plus F8 and USH2A for functional signatures. Canonical counts, selectors and the six main categories exclude those two genes. Functional and niche memberships retain them. HSPC159 is mapped to LGALSL and is not counted as a second gene.

The pair resource uses the exact 16,876 unique non-self pairs from `Theayomideo/matrispace` commit `518b7d2`, preserving provenance columns. Beta stores its endpoints as `Ligand`/`Receptor`; the offline app uses `Gene1`/`Gene2`. These pairs cover 929 canonical genes. Pair scoring and functional-signature membership are independent of the canonical count exclusion.

The Janesick breast cancer input has 990 detected canonical genes. Manuscript Table 1, S1, S3A and the round-two rebuttal comments were checked for this update on 2026-11-06. The S3A display label `TDGF1 (CRIPTO)` corresponds to reference symbol `CRIPTO`.

Run `Rscript dev/check_core.R` from the offline repository root to verify membership, pair coverage, arithmetic means and score refreshes.
