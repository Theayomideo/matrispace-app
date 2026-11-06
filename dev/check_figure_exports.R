# check_figure_exports.R
#
# Focused checks for the figure export (helpers.R) on a synthetic Visium
# fixture. Run from the app root:
#
#   Rscript dev/check_figure_exports.R

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(patchwork); library(ggrepel); library(pheatmap)
})
pkgload::load_all(".", quiet = TRUE)

check <- function(ok, what) {
  cat(if (isTRUE(ok)) "PASS" else "FAIL", what, "\n")
  if (!isTRUE(ok)) failures <<- failures + 1
}
failures <- 0

# Hexagonal Visium-like lattice: 100 um pitch at 2 px/um, so the Space Ranger
# 65 um scale-factor spot is 130 px.
fixture <- function(spot_px = 130, n_row = 12, n_col = 12) {
  grid <- expand.grid(r = seq_len(n_row), c = seq_len(n_col))
  coords <- data.frame(
    tissue = 1L, row = grid$r, col = grid$c,
    imagerow = 400 + grid$r * 200 * sqrt(3) / 2,
    imagecol = 400 + (grid$c + (grid$r %% 2) / 2) * 200,
    row.names = sprintf("BC%03d", seq_len(nrow(grid)))
  )
  counts <- matrix(rpois(20 * nrow(coords), 5), 20, dimnames = list(paste0("G", 1:20), rownames(coords)))
  obj <- CreateSeuratObject(as(counts, "dgCMatrix"), assay = "Spatial")
  set.seed(1)
  img <- new("VisiumV1", assay = "Spatial", key = "slice1_", coordinates = coords,
             scale.factors = scalefactors(spot = spot_px, fiducial = 200, hires = 0.2, lowres = 0.05),
             image = array(runif(160 * 160 * 3), c(160, 160, 3)), spot.radius = 0.01)
  obj@images$slice1 <- img
  obj$cell_type <- rep(c("A_cells", "B_cells", "not.assigned"), length.out = ncol(obj))
  obj
}

obj <- fixture()
cells <- colnames(obj)
data <- data.frame(row.names = cells, cell_type = obj$cell_type,
                   feature1 = runif(length(cells)), feature2 = runif(length(cells)),
                   LISA = sample(c("High.High", "High.Low", "Low.High", "Low.Low"), length(cells), TRUE))
snap <- list(
  sample = "fixture", time = Sys.time(), preset = "publication", obj = obj,
  annotation_col = "cell_type", annotation = setNames(obj$cell_type, cells),
  colours = c(A_cells = "#1b9e77", B_cells = "#d95f02", not.assigned = "#7f7f7f"),
  feature = list(data = data, gene1 = "G1", gene2 = "G2", sel1 = "matrisome gene", sel2 = "matrisome gene",
                 annotation_col = "cell_type", expression = list(assay = "Spatial", layer = "counts"),
                 selected = c("B_cells", "A_cells"), stale = FALSE, contrast = list(), violin = TRUE,
                 independent_y = FALSE, coex = list(palette = "Classic", alpha = 0.7, pt_size = 3)),
  matrisome = NULL, lr = NULL, views = list(main = list(), ecm = list(), pt = list())
)
style <- export_style("publication")

# Calibration from scale factors, cross-checked against the lattice.
ctx <- resolve_export_image_context(obj)
check(abs(ctx$px_per_um - 2) < 1e-6, "Visium calibration recovers 2 px/um")
check(!is.null(ctx$image), "registered histology is kept")
check(is.na(resolve_export_image_context(fixture(spot_px = 60))$px_per_um), "inconsistent spot pitch leaves scale uncalibrated")
obj_shuffled <- obj
obj_shuffled@images$slice1@coordinates <- obj@images$slice1@coordinates[rev(cells), ]
check(identical(resolve_export_image_context(obj_shuffled)$coords[cells, ], ctx$coords[cells, ]),
      "coordinates align by barcode, not row order")

out <- tempfile("figures_")
dirs <- lapply(list(pdf = "PDF", png = "PNG", svg = "SVG", components = "SVG/Components"), function(d) file.path(out, d))
for (d in dirs) dir.create(d, recursive = TRUE)
registry <- build_export_plot_registry(snap, style)
devices <- grDevices::dev.list()
broken <- list(id = "broken", stem = "99_broken", build = function() stop("induced failure"))
rows <- do.call(rbind, lapply(c(list(broken), registry$entries), write_export_plot, dirs = dirs, mode = "plots", style = style))
check(identical(grDevices::dev.list(), devices), "no graphics devices left open")
check(rows$status[rows$output_id == "broken"] == "error" && grepl("induced", rows$reason[rows$output_id == "broken"]),
      "a failing figure is reported with its error")
check(all(rows$status[rows$output_id == "primary_feature"] == "success"), "later figures still export after a failure")
check(all(rows$status[grepl("^matrisome_", rows$output_id)] == "not_available"), "missing analyses are not_available")

svg <- readLines(file.path(dirs$svg, "02_primary_feature.svg"), warn = FALSE, encoding = "UTF-8")
n_selected <- sum(data$cell_type %in% snap$feature$selected)
check(sum(lengths(regmatches(svg, gregexpr("<circle", svg)))) == n_selected, "feature map draws only the selected spots")
check(any(grepl("<text", svg)) && !any(grepl("textLength", svg)), "SVG text is live and not length-locked")
check(any(grepl("data:image/png;base64", svg)) && !any(grepl("href=['\"](file|http)", svg)), "SVG is self-contained")
check(any(grepl("viewBox", svg)), "SVG declares its size")
check(file.exists(file.path(dirs$components, "02_primary_feature__legend.svg")), "legend component is written")
check(any(grepl("[0-9] (\u00b5m|mm)", svg)), "calibrated scale bar is drawn")

bad_dirs <- dirs
writeLines("", file.path(out, "not_a_directory"))
bad_dirs$png <- file.path(out, "not_a_directory", "PNG")
png_fail <- write_export_plot(Filter(function(e) e$id == "primary_feature", registry$entries)[[1]], bad_dirs, "plots", style)
status <- write_export_manifest(png_fail, snap, style, registry$ctx, out, plot_dirs = dirs)
check(any(png_fail$status == "error") && grepl("0 figures exported, 0 not available, 1 failed", status),
      "a failed format counts the figure as failed")

snap$feature$selected <- character(0)
empty <- write_export_plot(Filter(function(e) e$id == "primary_feature", build_export_plot_registry(snap, style)$entries)[[1]],
                           dirs, "plots", style)
check(empty$status == "not_available", "an empty selection is not treated as select all")


pair_names <- c("COL1A1-DCN", "COL1A2-DCN", "COL1A1-COL1A2")
snap$pairs <- list(
  stats = data.frame(feature = pair_names, cluster = c("Basement", "Interstitial", "Basement"),
                     avg_log2FC = c(1.5, 2, 1), perc_difference = c(0.3, 0.4, 0.2)),
  means = matrix(c(8, 1, 6, 2, 9, 3), 3,
                 dimnames = list(pair_names, c("Basement", "Interstitial")))
)
pair_entries <- Filter(function(e) grepl("^pair_", e$id), build_export_plot_registry(snap, style)$entries)
pair_rows <- do.call(rbind, lapply(pair_entries, write_export_plot, dirs = dirs, mode = "plots", style = style))
check(nrow(pair_rows) > 0 && all(pair_rows$status == "success"), "pair heatmap and enrichment plots export")

cat(if (failures) sprintf("%d check(s) failed\n", failures) else "All checks passed\n")
quit(status = if (failures) 1 else 0)
