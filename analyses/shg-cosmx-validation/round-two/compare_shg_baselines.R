suppressPackageStartupMessages({
  library(data.table)
  library(hexbin)
  library(Seurat)
})
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Usage: Rscript compare_shg_baselines.R analysis_root [output_dir]")
root <- args[1]
out <- if (length(args) > 1L) args[2] else "results"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cells <- fread(file.path(root, "derived/section04_cell_collagen_scores.csv.gz"))
hex <- fread(file.path(root, "results/primary_hex_data_area_matched.csv.gz"))
obj <- readRDS(file.path(root, "derived/section04_matrispace_collagen_scored.rds"))
sct <- SeuratObject::LayerData(obj, assay = "SCT", layer = "data")
stopifnot(!anyDuplicated(cells$barcode), setequal(cells$barcode, colnames(sct)),
          !anyNA(cells$published_celltype), all(c("COL1A1", "DCN") %in% rownames(sct)))
for (gene in c("COL1A1", "DCN"))
  cells[[gene]] <- as.numeric(sct[gene, cells$barcode])
bins <- hexbin(cells$x_um, cells$y_um, xbins = 37L, xbnds = c(0,4000),
               ybnds = c(0,4000), shape = 1, IDs = TRUE)
cells[, hex_id := bins@cID]
cell_means <- cells[, .(n = .N, full_check = mean(collagens_ucell),
                         COL1A1 = mean(COL1A1), DCN = mean(DCN),
                         fibroblasts = sum(published_celltype == "Fibroblasts")), by = hex_id]
h <- merge(hex, cell_means, by = "hex_id")
stopifnot(nrow(h) == 1545, all(h$n == h$n_cells),
          max(abs(h$full_check - h$collagens_ucell)) < 1e-12)
# Use Figure 6's saved within-ROI percentiles for rank calculations. They
# preserve the exact ordering of the original means, including close ties.
# Clipped hexagons at the ROI boundary use their actual imaged area.
h[, fibroblast_density := fibroblasts / (shg_pixels * (4000/1320)^2 / 1e6)]
h[, fibroblast_fraction := fibroblasts / n_cells]
measures <- c("collagens_ucell","COL1A1","DCN","fibroblast_density",
              "fibroblast_fraction","fibrillar_collagens_ucell")
rho <- function(a,b) cor(a,b,method = "spearman")
calc <- function(idx) vapply(measures, function(m) {
  x <- if (m=="collagens_ucell") h$cosmx_percentile else h[[m]]
  rho(x[idx],h$shg_percentile[idx])
}, numeric(1))
obs <- calc(seq_len(nrow(h)))
stopifnot(abs(obs[1] - 0.489245181950197) < 1e-10)
blocks <- split(seq_len(nrow(h)), paste(floor(h$x_um/400),floor(h$y_um/400),sep="_"))
set.seed(15240431L + 101L)
boot <- replicate(1000, {
  idx <- unlist(blocks[sample(names(blocks), length(blocks), replace=TRUE)], use.names=FALSE)
  calc(idx)
})
ci <- function(x) quantile(x, c(.025,.975), type=8, names=FALSE)
ci_each <- t(apply(boot,1,ci))
diff_ci <- t(apply(sweep(-boot,2,boot[1,],"+"),1,ci))
result <- data.frame(measure=measures, rho=obs, lo=ci_each[,1], hi=ci_each[,2],
                     full_minus_measure=obs[1]-obs, diff_lo=diff_ci[,1], diff_hi=diff_ci[,2])
fwrite(result,file.path(out,"shg_comparators.csv"))
rank_residual <- function(x, covars) resid(lm(rank(x) ~ ., data=covars))
adjust <- function(vars) {
  covars <- as.data.frame(lapply(h[, ..vars],rank))
  cor(rank_residual(h$shg_percentile,covars),rank_residual(h$cosmx_percentile,covars))
}
controls <- list(original=c("n_cells","nFeature_Spatial"),
                 fibroblast_density=c("n_cells","nFeature_Spatial","fibroblast_density"),
                 all_comparators=c("n_cells","nFeature_Spatial","fibroblast_density","COL1A1","DCN"))
partial <- vapply(controls,adjust,numeric(1))
stopifnot(abs(partial[1]-0.478714566335356) < 1e-10)
fwrite(data.frame(adjustment=names(partial),partial_rho=partial),file.path(out,"shg_partial.csv"))
fwrite(h[,c("hex_id","n_cells","nFeature_Spatial","shg_mean",measures),with=FALSE],
       file.path(out,"shg_comparator_bins.csv"))
print(result,digits=4); print(partial)
