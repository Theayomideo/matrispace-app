options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
analysis_root <- if (length(args) >= 1L) args[[1L]] else "orthogonal_validation"
input_root <- file.path(analysis_root, "input")
output_dir <- file.path(analysis_root, "derived")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required <- c(
  "data.table", "Matrix", "Seurat", "SeuratObject", "UCell",
  "matrispace", "tiff"
)
missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

expr_path <- file.path(input_root, "cosmx_section04", "Run1129_Section_04_exprMat_file.csv")
meta_path <- file.path(input_root, "cosmx_section04", "Run1129_Section_04_metadata_file.csv")
aligned_path <- file.path(input_root, "shg", "section4_aligned_metadata.csv")
coordinate_path <- file.path(input_root, "shg", "coordinates.csv")
collagen_path <- file.path(input_root, "shg", "Collagen.tif")
signature_path <- file.path(input_root, "matrispace_0.3.0", "matrisome_signatures.rda")

paths <- c(expr_path, meta_path, aligned_path, coordinate_path, collagen_path, signature_path)
if (any(!file.exists(paths))) stop("Missing input(s): ", paste(paths[!file.exists(paths)], collapse = ", "))

message("Reading raw Section 04 expression matrix...")
expr <- data.table::fread(expr_path, showProgress = TRUE)
expr <- expr[cell_ID != 0]
expr[, barcode := sprintf("section_4_fov_%s_ID_%s", fov, cell_ID)]
if (anyDuplicated(expr$barcode)) stop("Expression barcodes are duplicated")

message("Reading raw cell metadata and deposited alignment...")
raw_meta <- data.table::fread(meta_path)
raw_meta <- raw_meta[cell_ID != 0]
raw_meta[, barcode := sprintf("section_4_fov_%s_ID_%s", fov, cell_ID)]
if (anyDuplicated(raw_meta$barcode)) stop("Metadata barcodes are duplicated")

aligned <- data.table::fread(aligned_path)
coords <- data.table::fread(coordinate_path)
if (nrow(aligned) != nrow(coords)) stop("Aligned metadata and coordinates.csv have different row counts")
if (anyDuplicated(aligned$barcode)) stop("Aligned barcodes are duplicated")
if (!isTRUE(all.equal(aligned$xc, coords$xc)) || !isTRUE(all.equal(aligned$yc, coords$yc))) {
  stop("Aligned metadata and coordinates.csv do not describe the same coordinates")
}

expr_index <- match(aligned$barcode, expr$barcode)
meta_index <- match(aligned$barcode, raw_meta$barcode)
if (anyNA(expr_index) || anyNA(meta_index)) stop("Some aligned cells are absent from the raw files")
expr <- expr[expr_index]
raw_meta <- raw_meta[meta_index]
if (!identical(expr$barcode, aligned$barcode) || !identical(raw_meta$barcode, aligned$barcode)) {
  stop("Failed to preserve the deposited aligned-cell order")
}

id_columns <- c("fov", "cell_ID", "barcode")
all_features <- setdiff(names(expr), id_columns)
technical_features <- grep(
  "^(NegPrb|Negative|SystemControl|FalseCode|Control)",
  all_features, value = TRUE, ignore.case = TRUE
)
gene_features <- setdiff(all_features, technical_features)

message("Building sparse gene-by-cell count matrix...")
dense_counts <- t(as.matrix(expr[, ..gene_features]))
storage.mode(dense_counts) <- "double"
counts <- methods::as(Matrix::Matrix(dense_counts, sparse = TRUE), "dgCMatrix")
rownames(counts) <- gene_features
colnames(counts) <- aligned$barcode
rm(dense_counts, expr)
invisible(gc())

collagen_img <- tiff::readTIFF(collagen_path, native = FALSE, as.is = TRUE)
if (length(dim(collagen_img)) == 3L) collagen_img <- collagen_img[, , 1L]
storage.mode(collagen_img) <- "double"
image_height <- nrow(collagen_img)
image_width <- ncol(collagen_img)

# The deposited coordinates are already aligned to the Section 03 SHG crop.
# Recompute the published linear map explicitly with 1-based pixel coordinates.
x_fraction <- (aligned$xc - min(aligned$xc)) / diff(range(aligned$xc))
y_fraction <- (aligned$yc - min(aligned$yc)) / diff(range(aligned$yc))
aligned[, image_x_px := 1 + (1 - x_fraction) * (image_width - 1)]
aligned[, image_y_px := 1 + y_fraction * (image_height - 1)]

# The deposited workflow documents a 4000-um field across the 1320-pixel crop.
field_width_um <- 4000
pixel_size_um <- field_width_um / image_width
aligned[, x_um := (image_x_px - 1) * pixel_size_um]
aligned[, y_um := (image_y_px - 1) * pixel_size_um]

keep_meta <- intersect(
  c(
    "fov", "cell_ID", "Area", "AspectRatio", "CenterX_local_px",
    "CenterY_local_px", "CenterX_global_px", "CenterY_global_px",
    "Mean.MembraneStain", "Mean.PanCK", "Mean.CD45", "Mean.CD3", "Mean.DAPI"
  ),
  names(raw_meta)
)
seurat_meta <- as.data.frame(raw_meta[, ..keep_meta])
rownames(seurat_meta) <- aligned$barcode
seurat_meta$image_x_px <- aligned$image_x_px
seurat_meta$image_y_px <- aligned$image_y_px
seurat_meta$x_um <- aligned$x_um
seurat_meta$y_um <- aligned$y_um
if ("celltypes" %in% names(aligned)) seurat_meta$published_celltype <- aligned$celltypes

message("Creating and normalizing the Section 04 Seurat object...")
obj <- Seurat::CreateSeuratObject(
  counts = counts,
  assay = "Spatial",
  meta.data = seurat_meta,
  min.cells = 1,
  min.features = 0,
  project = "LUAD_SHG_CosMx_validation"
)
rm(counts)
invisible(gc())

# Match the published CosMx processing while retaining every assayed gene.
obj <- Seurat::SCTransform(
  obj,
  assay = "Spatial",
  new.assay.name = "SCT",
  do.correct.umi = TRUE,
  return.only.var.genes = FALSE,
  conserve.memory = TRUE,
  verbose = TRUE
)

sig_env <- new.env(parent = emptyenv())
load(signature_path, envir = sig_env)
bundled_signatures <- sig_env[[ls(sig_env)[1L]]]
official_collagens <- bundled_signatures[["collagens"]]

# Fibril-forming collagen types I, II, III, V, XI, XXIV, and XXVII.
fibrillar_collagens <- c(
  "COL1A1", "COL1A2", "COL2A1", "COL3A1",
  "COL5A1", "COL5A2", "COL5A3", "COL11A1", "COL11A2",
  "COL24A1", "COL27A1"
)

score_definitions <- list(
  collagens = official_collagens,
  fibrillar_collagens = fibrillar_collagens
)

message("Scoring collagen signatures with MatriSpace 0.3.0...")
obj <- matrispace::score_matrisome(
  obj,
  signatures = score_definitions,
  assay = "SCT",
  layer = "data",
  min_genes = 3L,
  maxRank = 1500L
)
score_matrix <- as.matrix(matrispace::get_matrisome_scores(obj))
# Seurat replaces underscores in assay feature names with hyphens. MatriSpace
# 0.3.0 currently restores canonical pair ids for LR results but not canonical
# signature keys for matrisome results, so reconcile against the definitions.
stored_keys <- rownames(score_matrix)
canonical_keys <- names(score_definitions)
key_index <- match(canonical_keys, stored_keys)
missing_key <- is.na(key_index)
key_index[missing_key] <- match(gsub("_", "-", canonical_keys[missing_key]), stored_keys)
if (anyNA(key_index)) {
  stop(
    "Could not reconcile stored MatriSpace score keys. Stored: ",
    paste(stored_keys, collapse = ", ")
  )
}
score_matrix <- score_matrix[key_index, , drop = FALSE]
rownames(score_matrix) <- canonical_keys
scores <- t(score_matrix)

sct_data <- SeuratObject::LayerData(obj, assay = "SCT", layer = "data")
matched <- lapply(score_definitions, intersect, y = rownames(sct_data))
normalized_sums <- vapply(
  matched,
  function(g) as.numeric(Matrix::colSums(sct_data[g, , drop = FALSE])),
  numeric(ncol(obj))
)
colnames(normalized_sums) <- paste0(colnames(normalized_sums), "_sct_sum")
rownames(normalized_sums) <- colnames(obj)

cell_scores <- data.frame(
  barcode = colnames(obj),
  image_x_px = obj$image_x_px,
  image_y_px = obj$image_y_px,
  x_um = obj$x_um,
  y_um = obj$y_um,
  fov = obj$fov,
  nCount_Spatial = obj$nCount_Spatial,
  nFeature_Spatial = obj$nFeature_Spatial,
  collagens_ucell = scores[colnames(obj), "collagens"],
  fibrillar_collagens_ucell = scores[colnames(obj), "fibrillar_collagens"],
  collagens_sct_sum = normalized_sums[colnames(obj), "collagens_sct_sum"],
  fibrillar_collagens_sct_sum = normalized_sums[colnames(obj), "fibrillar_collagens_sct_sum"],
  stringsAsFactors = FALSE
)
if ("published_celltype" %in% colnames(obj@meta.data)) {
  cell_scores$published_celltype <- obj$published_celltype
}

data.table::fwrite(cell_scores, file.path(output_dir, "section04_cell_collagen_scores.csv.gz"))
saveRDS(obj, file.path(output_dir, "section04_matrispace_collagen_scored.rds"), compress = "gzip")
saveRDS(collagen_img, file.path(output_dir, "section03_shg_collagen_image.rds"), compress = "gzip")

coverage <- data.frame(
  signature = names(score_definitions),
  requested_genes = lengths(score_definitions),
  measured_genes = lengths(matched),
  coverage_fraction = lengths(matched) / lengths(score_definitions),
  measured = vapply(matched, paste, character(1L), collapse = ";"),
  missing = vapply(
    names(score_definitions),
    function(k) paste(setdiff(score_definitions[[k]], rownames(sct_data)), collapse = ";"),
    character(1L)
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(coverage, file.path(output_dir, "score_signature_coverage.csv"), row.names = FALSE)

run_metadata <- data.frame(
  field = c(
    "timestamp", "r_version", "matrispace_version", "seurat_version",
    "seuratobject_version", "ucell_version", "cells", "assayed_genes",
    "technical_features_excluded", "normalization", "score_assay",
    "score_layer", "requested_maxRank", "effective_maxRank",
    "shg_image_width_px", "shg_image_height_px", "shg_field_width_um",
    "shg_pixel_size_um"
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    R.version.string,
    as.character(utils::packageVersion("matrispace")),
    as.character(utils::packageVersion("Seurat")),
    as.character(utils::packageVersion("SeuratObject")),
    as.character(utils::packageVersion("UCell")),
    ncol(obj),
    nrow(obj),
    length(technical_features),
    "SCTransform; return.only.var.genes=FALSE; conserve.memory=TRUE",
    "SCT", "data", 1500,
    min(1500L, nrow(sct_data)),
    image_width, image_height, field_width_um, pixel_size_um
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(run_metadata, file.path(output_dir, "score_run_metadata.csv"), row.names = FALSE)

cat("Prepared and scored", ncol(obj), "aligned cells\n")
print(coverage[, 1:4], row.names = FALSE)
