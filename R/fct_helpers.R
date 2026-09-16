# helpers.R
# This file contains all helper functions used in the server logic.

#' Drop SCT scale.data, reductions, and graphs from a Seurat object before
#' storing in rv$active_seurat_object. scale.data is the biggest RAM offender
#' post-ScType. Returns the original object unchanged if DietSeurat throws.
slim_seurat_for_app <- function(obj) {
  if (!inherits(obj, "Seurat")) return(obj)
  tryCatch({
    Seurat::DietSeurat(obj, layers = c("counts", "data"),
               dimreducs = character(0), graphs = character(0))
  }, error = function(e) {
    warning(sprintf("slim_seurat_for_app: DietSeurat failed (%s); keeping full object.",
                    conditionMessage(e)))
    obj
  })
}

#' Safe assay data accessor for Seurat v4/v5 compatibility
#'
#' @param object A Seurat object
#' @param assay Assay name (defaults to SCT > Spatial > first available)
#' @param slot Data slot: "counts", "data", or "scale.data"
#' @return Sparse matrix of expression data or NULL if not found
safe_get_assay_data <- function(object, assay = NULL, slot = NULL) {
  # Step 1: Determine which assay to use
  if (is.null(assay)) {
    if ("SCT" %in% SeuratObject::Assays(object)) {
      assay <- "SCT"
    } else if ("Spatial" %in% SeuratObject::Assays(object)) {
      assay <- "Spatial"
    } else if (length(SeuratObject::Assays(object)) >= 1) {
      assay <- SeuratObject::Assays(object)[1]
    } else {
      warning("No assays found in object")
      return(NULL)
    }
  }

  # Step 2: Get the assay object
  assay_obj <- object[[assay]]
  is_v5 <- inherits(assay_obj, "Assay5")

  # Step 3: Try to get data in preferred order
  tryCatch({
    if (is_v5) {
      # For Seurat v5
      if (is.null(slot)) {
        # Try slots in preferred order
        if (!is.null(assay_obj$counts)) return(assay_obj$counts)
        if (!is.null(assay_obj$data)) return(assay_obj$data)
        if (!is.null(assay_obj$scale.data)) return(assay_obj$scale.data)
        warning(sprintf("No data found in any slot for assay %s", assay))
        return(NULL)
      } else {
        # If slot is specified, try that specific slot
        if (slot == "counts") return(assay_obj$counts)
        if (slot == "data") return(assay_obj$data)
        if (slot == "scale.data") return(assay_obj$scale.data)
      }
    } else {
      # For older Seurat versions
      if (is.null(slot)) {
        # Try slots in preferred order
        tryCatch({
          return(GetAssayData(object = object, assay = assay, layer = "counts"))
        }, error = function(e) {
          tryCatch({
            return(GetAssayData(object = object, assay = assay, layer = "data"))
          }, error = function(e) {
            tryCatch({
              return(GetAssayData(object = object, assay = assay, layer = "scale.data"))
            }, error = function(e) {
              warning(sprintf("No data found in any slot for assay %s", assay))
              return(NULL)
            })
          })
        })
      } else {
        # If slot is specified, try that specific slot
        return(GetAssayData(object = object, assay = assay, layer = slot))
      }
    }
  }, error = function(e) {
    warning(sprintf("Error accessing %s slot in %s assay: %s", slot, assay, e$message))
    return(NULL)
  })
}

#' Get gene names from Seurat object safely
#'
#' @param object A Seurat object
#' @param assay Assay name (optional)
#' @return Character vector of gene names or NULL
safe_get_rownames <- function(object, assay = NULL) {
  data <- safe_get_assay_data(object, assay = assay)
  if (!is.null(data)) return(rownames(data))
  return(NULL)
}

#' Standardize gene symbols to current HGNC nomenclature
#'
#' Gets the HGNC mapping via scCustomize, then renames features at the slot
#' level across all assays (Assay5 and v4 Assay). Duplicates after renaming
#' are reverted to the original symbol to avoid collisions.
#'
#' @param seurat_obj A Seurat object
#' @return Seurat object with updated gene symbols
#' Update gene symbols against a bundled HGNC long-format table
#'
#' Offline replacement for scCustomize::Updated_HGNC_Symbols. Looks up each
#' input symbol in the HGNC table built by scripts/build_hgnc_table.R: if
#' it's already an approved current symbol it stays, if it matches a
#' previous symbol it gets renamed to the current one, otherwise it stays.
#'
#' @param genes Character vector of gene symbols.
#' @param hgnc_long_data Data frame with columns "symbol" and "prev_symbol".
#' @return Data frame with input_features and Output_Features (same shape as
#'   scCustomize::Updated_HGNC_Symbols), so callers can swap in directly.
update_hgnc_symbols <- function(genes, hgnc_long_data) {
  approved <- unique(hgnc_long_data$symbol)
  # Map prev_symbol -> current symbol; drop rows where prev is itself an
  # approved symbol so we never demote a current name.
  rename_tbl <- hgnc_long_data[!hgnc_long_data$prev_symbol %in% approved, , drop = FALSE]
  rename_map <- setNames(rename_tbl$symbol, rename_tbl$prev_symbol)
  rename_map <- rename_map[!duplicated(names(rename_map))]

  output <- ifelse(
    genes %in% approved,
    genes,
    ifelse(genes %in% names(rename_map), unname(rename_map[genes]), genes)
  )

  data.frame(
    input_features = genes,
    Output_Features = output,
    stringsAsFactors = FALSE
  )
}

standardize_gene_symbols <- function(
  seurat_obj,
  hgnc_long_data = readRDS(extdata_path("hgnc_long_data.rds"))
) {
  # Collect all unique features across every assay
  all_features <- unique(unlist(lapply(SeuratObject::Assays(seurat_obj), function(a) {
    rownames(seurat_obj[[a]])
  })))

  # Get HGNC mapping using the bundled offline table.
  map_df <- update_hgnc_symbols(all_features, hgnc_long_data)

  changed_idx <- which(map_df$input_features != map_df$Output_Features)
  if (length(changed_idx) == 0) {
    message("Gene symbols: all up to date.")
    return(seurat_obj)
  }

  # Build a safe rename map: only include renames where the target name
  # does not collide with any existing (unchanged) feature name
  occupied <- map_df$input_features  # all current names are "taken"
  safe_rename <- setNames(character(0), character(0))

  for (i in changed_idx) {
    target <- map_df$Output_Features[i]
    source <- map_df$input_features[i]
    if (!(target %in% occupied)) {
      safe_rename[source] <- target
      occupied[i] <- target  # mark this slot as the new name
    }
  }

  n_changed <- length(safe_rename)
  if (n_changed == 0) {
    message("Gene symbols: all up to date (all updates skipped due to collisions).")
    return(seurat_obj)
  }

  # Helper: apply safe_rename to a feature-name vector
  apply_rename <- function(feat_names) {
    idx <- match(feat_names, names(safe_rename))
    has_match <- !is.na(idx)
    feat_names[has_match] <- safe_rename[idx[has_match]]
    feat_names
  }

  # Rename at the slot level so feature metadata and matrices stay consistent
  for (assay_name in SeuratObject::Assays(seurat_obj)) {
    assay_obj <- seurat_obj[[assay_name]]

    if (inherits(assay_obj, "Assay5")) {
      # --- Seurat v5 Assay5 ---
      # 1. Rename master feature list
      rownames(assay_obj@features) <- apply_rename(rownames(assay_obj@features))
      # 2. Rename each layer matrix
      for (lname in names(assay_obj@layers)) {
        rn <- rownames(assay_obj@layers[[lname]])
        if (!is.null(rn)) rownames(assay_obj@layers[[lname]]) <- apply_rename(rn)
      }
    } else {
      # --- Seurat v4 Assay ---
      # Rename slot matrices
      if (nrow(assay_obj@counts) > 0)
        rownames(assay_obj@counts) <- apply_rename(rownames(assay_obj@counts))
      if (nrow(assay_obj@data) > 0)
        rownames(assay_obj@data) <- apply_rename(rownames(assay_obj@data))
      if (length(assay_obj@scale.data) > 0 && nrow(assay_obj@scale.data) > 0)
        rownames(assay_obj@scale.data) <- apply_rename(rownames(assay_obj@scale.data))
      # Rename feature metadata
      if (nrow(assay_obj@meta.features) > 0)
        rownames(assay_obj@meta.features) <- apply_rename(rownames(assay_obj@meta.features))
    }

    # Update variable features
    vf <- SeuratObject::VariableFeatures(assay_obj)
    if (length(vf) > 0) {
      SeuratObject::VariableFeatures(assay_obj) <- apply_rename(vf)
    }

    seurat_obj[[assay_name]] <- assay_obj
  }

  message(paste0("Gene symbols: updated ", n_changed, " to current HGNC nomenclature."))
  return(seurat_obj)
}

#' Get feature count from Seurat object safely
#'
#' @param object A Seurat object
#' @param assay Assay name (optional)
#' @return Integer count of features or 0
safe_get_nfeatures <- function(object, assay = NULL) {
  data <- safe_get_assay_data(object, assay = assay)
  if (!is.null(data)) return(nrow(data))
  return(0)
}

#' Combine spatial coordinates with metadata for plotting.
prepare_data_plot <- function(obj) {
  spot_coords <- align_spatial_coordinates(obj)
  metadata <- obj@meta.data
  metadata <- metadata[rownames(spot_coords), , drop = FALSE]
  data_plot <- cbind(spot_coords, metadata)
  data_plot
}

#' Align a named vector to Seurat cell order when names are available.
align_to_cells <- function(values, cells) {
  if (!is.null(names(values)) && all(cells %in% names(values))) {
    values <- values[cells]
  }
  values
}

#' Return spatial coordinates in the requested Seurat cell order.
align_spatial_coordinates <- function(obj, coords = NULL, cells = NULL) {
  if (is.null(coords)) coords <- GetTissueCoordinates(obj)
  coords <- as.data.frame(coords)
  if (is.null(cells)) cells <- rownames(obj@meta.data)

  coord_cells <- rownames(coords)
  if ((is.null(coord_cells) || anyNA(coord_cells) || any(!nzchar(coord_cells))) &&
      "cell" %in% colnames(coords)) {
    coord_cells <- as.character(coords$cell)
  }
  if ((is.null(coord_cells) || anyNA(coord_cells) || any(!nzchar(coord_cells))) &&
      "barcode" %in% colnames(coords)) {
    coord_cells <- as.character(coords$barcode)
  }

  if (!is.null(coord_cells) && length(coord_cells) == nrow(coords) &&
      !anyNA(coord_cells) && all(nzchar(coord_cells))) {
    rownames(coords) <- make.unique(coord_cells)
    if (all(cells %in% rownames(coords))) {
      coords <- coords[cells, , drop = FALSE]
    } else {
      common_cells <- intersect(cells, rownames(coords))
      if (length(common_cells) > 0) {
        coords <- coords[common_cells, , drop = FALSE]
      }
    }
  }

  if (!all(c("row", "col") %in% names(coords))) {
    if (length(obj@images) > 0 && "coordinates" %in% slotNames(obj@images[[1]])) {
      raw_coords <- obj@images[[1]]@coordinates
      raw_coords <- align_spatial_coordinates(obj, raw_coords, rownames(coords))
      coords$row <- raw_coords$row
      coords$col <- raw_coords$col
    } else if (all(c("x", "y") %in% names(coords))) {
      coords$row <- coords$y
      coords$col <- coords$x
    }
  }

  if (!all(c("imagerow", "imagecol") %in% names(coords)) &&
      all(c("x", "y") %in% names(coords))) {
    coords$imagerow <- coords$x
    coords$imagecol <- coords$y
  }

  coords
}

#' Convert VisiumV2 images to VisiumV1 format for compatibility
#'
#' Seurat 5.1+ uses VisiumV2 by default, which is incompatible with some
#' functions in Seurat 5.0.x. This function converts VisiumV2 images to
#' VisiumV1 format to enable spatial plotting.
#'
#' @param obj A Seurat object potentially containing VisiumV2 images
#' @param verbose Logical, whether to print conversion messages
#' @return Seurat object with VisiumV1 images (or unchanged if no VisiumV2)
convert_visiumv2_to_v1 <- function(obj, verbose = TRUE) {
  # Access @images slot directly - don't use Images() which fails on VisiumV2
  img_names <- names(obj@images)

  if (length(img_names) == 0) {
    return(obj)
  }

  for (img_name in img_names) {
    img <- obj@images[[img_name]]

    # Check if this is a VisiumV2 image
    if (!inherits(img, "VisiumV2")) {
      next
    }

    if (verbose) message(sprintf("Converting VisiumV2 image '%s' to VisiumV1...", img_name))

    tryCatch({
      # Extract components from VisiumV2 using slot access
      assay_name <- tryCatch(slot(img, "assay"), error = function(e) "Spatial")
      key_val <- tryCatch(slot(img, "key"), error = function(e) paste0(img_name, "_"))
      scale_factors <- tryCatch(slot(img, "scale.factors"), error = function(e) NULL)
      image_data <- tryCatch(slot(img, "image"), error = function(e) NULL)

      # Get coordinates - try multiple approaches
      coords <- NULL

      # Method 1: Try cells slot
      coords <- tryCatch({
        cells_data <- slot(img, "cells")
        if (!is.null(cells_data) && (is.data.frame(cells_data) || is.matrix(cells_data))) {
          as.data.frame(cells_data)
        } else NULL
      }, error = function(e) NULL)

      # Method 2: Try boundaries$centroids slot (VisiumV2 structure)
      if (is.null(coords)) {
        coords <- tryCatch({
          boundaries <- slot(img, "boundaries")
          if (!is.null(boundaries) && "centroids" %in% names(boundaries)) {
            centroids <- boundaries[["centroids"]]
            # GetTissueCoordinates returns x, y, cell columns
            tc <- GetTissueCoordinates(centroids)
            if (!is.null(tc) && nrow(tc) > 0) {
              # Convert to VisiumV1 format: x -> imagerow, y -> imagecol
              data.frame(
                tissue = rep(1L, nrow(tc)),
                row = seq_len(nrow(tc)),
                col = seq_len(nrow(tc)),
                imagerow = tc$x,
                imagecol = tc$y,
                row.names = tc$cell
              )
            } else NULL
          } else NULL
        }, error = function(e) NULL)
      }

      # Method 3: Build from object metadata
      if (is.null(coords)) {
        coords <- tryCatch({
          meta <- obj@meta.data
          coord_cols <- c("imagerow", "imagecol", "pxl_row_in_fullres", "pxl_col_in_fullres",
                          "x", "y", "spatial_x", "spatial_y")
          found_cols <- intersect(coord_cols, colnames(meta))
          if (length(found_cols) >= 2) {
            df <- data.frame(row.names = rownames(meta))
            # Try to find imagerow/imagecol equivalents
            # Note: For VisiumV2, x maps to imagerow and y maps to imagecol
            if ("imagerow" %in% found_cols) df$imagerow <- meta$imagerow
            else if ("pxl_row_in_fullres" %in% found_cols) df$imagerow <- meta$pxl_row_in_fullres
            else if ("x" %in% found_cols) df$imagerow <- meta$x
            else if ("spatial_x" %in% found_cols) df$imagerow <- meta$spatial_x

            if ("imagecol" %in% found_cols) df$imagecol <- meta$imagecol
            else if ("pxl_col_in_fullres" %in% found_cols) df$imagecol <- meta$pxl_col_in_fullres
            else if ("y" %in% found_cols) df$imagecol <- meta$y
            else if ("spatial_y" %in% found_cols) df$imagecol <- meta$spatial_y

            if (ncol(df) >= 2) df else NULL
          } else NULL
        }, error = function(e) NULL)
      }

      if (is.null(coords)) {
        warning(sprintf("Could not extract coordinates from VisiumV2 image '%s'. Skipping conversion.", img_name))
        next
      }

      # Ensure coordinates have required columns for VisiumV1
      if (!"tissue" %in% colnames(coords)) coords$tissue <- 1
      if (!"row" %in% colnames(coords)) {
        if ("array_row" %in% colnames(coords)) coords$row <- coords$array_row
        else coords$row <- seq_len(nrow(coords))
      }
      if (!"col" %in% colnames(coords)) {
        if ("array_col" %in% colnames(coords)) coords$col <- coords$array_col
        else coords$col <- seq_len(nrow(coords))
      }
      if (!"imagerow" %in% colnames(coords)) coords$imagerow <- coords$row * 100
      if (!"imagecol" %in% colnames(coords)) coords$imagecol <- coords$col * 100

      # Create scale factors if missing
      if (is.null(scale_factors)) {
        scale_factors <- tryCatch(
          scalefactors(spot = 0.02, fiducial = 0.1, hires = 0.17, lowres = 0.05),
          error = function(e) NULL
        )
      }

      # Create image matrix if missing
      if (is.null(image_data) || !is.array(image_data)) {
        image_data <- array(0.95, dim = c(100, 100, 3))
      }

      # Create VisiumV1 object
      visium_v1 <- new(
        Class = "VisiumV1",
        assay = assay_name,
        key = key_val,
        coordinates = coords,
        scale.factors = scale_factors,
        image = image_data
      )

      # Set spot radius - VisiumV2 scale.factors$spot is in pixels, need to convert
      # to normalized coordinates (fraction of image) that VisiumV1 expects
      spot_radius <- tryCatch({
        sf <- scale_factors
        # scale.factors$spot is spot diameter in full-res pixels (~70-90 typical)
        # VisiumV1 expects spot.radius as a small normalized value (~0.01)
        if (sf$spot > 1) {
          # Convert: (spot_diameter_pixels / 2) / fullres_image_size
          img_height <- dim(image_data)[1]
          fullres_height <- img_height / sf$lowres
          (sf$spot / 2) / fullres_height
        } else {
          # Already normalized
          sf$spot / 2
        }
      }, error = function(e) {
        0.01  # Safe default similar to working collection objects
      })

      # Ensure spot_radius is valid and in reasonable range
      if (length(spot_radius) == 0 || is.null(spot_radius) || is.na(spot_radius) || spot_radius > 0.1) {
        spot_radius <- 0.01
      }

      visium_v1@spot.radius <- spot_radius

      # Align cells with object
      v1_cells <- tryCatch(Cells(visium_v1), error = function(e) rownames(coords))
      obj_cells <- Cells(obj)
      common_cells <- intersect(v1_cells, obj_cells)

      if (length(common_cells) > 0 && length(common_cells) < length(v1_cells)) {
        visium_v1 <- visium_v1[common_cells]
      }

      # Replace in object
      obj@images[[img_name]] <- visium_v1

      if (verbose) message(sprintf("Successfully converted '%s' to VisiumV1", img_name))

    }, error = function(e) {
      warning(sprintf("Failed to convert VisiumV2 image '%s': %s", img_name, e$message))
    })
  }

  obj
}

#' Return detected canonical matrisome genes.
available_matrisome_genes <- function(obj, matrisome_genes) {
  object_genes <- safe_get_rownames(obj, assay = DefaultAssay(obj))
  sort(unique(intersect(object_genes, matrisome_genes)))
}

#' Detect custom UCell scores with limited colour contrast.
diagnose_ucell_contrast <- function(scores, upper_cutoff = 0.95,
                                    min_upper_fraction = 0.40, max_iqr = 0.075) {
  values <- as.numeric(scores)
  values <- values[is.finite(values)]
  upper_fraction <- if (length(values) > 0) mean(values >= upper_cutoff) else NA_real_
  score_iqr <- if (length(values) > 0) unname(stats::IQR(values)) else NA_real_
  has_variation <- length(unique(values)) > 1

  list(
    flagged = isTRUE(has_variation && upper_fraction >= min_upper_fraction && score_iqr <= max_iqr),
    upper_fraction = upper_fraction,
    iqr = score_iqr,
    upper_cutoff = upper_cutoff
  )
}

#' Convert UCell scores to within-tissue percentiles for display.
ucell_relative_contrast <- function(scores) {
  values <- as.numeric(scores)
  result <- rep(NA_real_, length(values))
  names(result) <- names(scores)
  keep <- is.finite(values)
  if (sum(keep) < 2 || length(unique(values[keep])) < 2) {
    result[keep] <- 0.5
    return(result)
  }
  result[keep] <- (base::rank(values[keep], ties.method = "average") - 1) / (sum(keep) - 1)
  result
}

#' Score custom matrisome gene lists with UCell.
score_custom_gene_sets <- function(obj, signatures) {
  assay <- DefaultAssay(obj)
  requested_genes <- unique(unlist(signatures))
  layer <- "data"
  expression <- safe_get_assay_data(obj, assay = assay, slot = layer)

  if (is.null(expression) || nrow(expression) == 0 ||
      !all(requested_genes %in% rownames(expression))) {
    layer <- "counts"
    expression <- safe_get_assay_data(obj, assay = assay, slot = layer)
  }
  if (is.null(expression) || nrow(expression) == 0) {
    stop("Could not retrieve expression data for custom matrisome scoring.")
  }

  missing_genes <- setdiff(requested_genes, rownames(expression))
  if (length(missing_genes) > 0) {
    stop("The custom matrisome list contains genes not found in the active assay.")
  }

  scores <- as.data.frame(
    UCell::ScoreSignatures_UCell(
      expression,
      features = signatures,
      BPPARAM = BiocParallel::SerialParam()
    ),
    check.names = FALSE
  )
  cells <- rownames(obj@meta.data)
  score_columns <- paste0(names(signatures), "_UCell")

  if (!all(score_columns %in% colnames(scores)) ||
      is.null(rownames(scores)) || !all(cells %in% rownames(scores))) {
    stop("Custom matrisome scores could not be aligned to all spots.")
  }

  scores <- scores[cells, score_columns, drop = FALSE]
  colnames(scores) <- names(signatures)
  if (any(!is.finite(as.matrix(scores)))) {
    stop("Custom matrisome scoring returned non-finite values.")
  }

  diagnostics <- lapply(names(signatures), function(id) {
    diagnose_ucell_contrast(scores[[id]])
  })
  names(diagnostics) <- names(signatures)

  contrast_scores <- lapply(names(signatures), function(id) {
    if (!diagnostics[[id]]$flagged) return(NULL)
    values <- scores[[id]]
    names(values) <- rownames(scores)
    ucell_relative_contrast(values)
  })
  names(contrast_scores) <- names(signatures)

  list(
    scores = scores,
    assay = assay,
    layer = layer,
    diagnostics = diagnostics,
    contrast_scores = contrast_scores
  )
}

#' Add feature analysis results and LISA statistics to Seurat object.
addfeat <- function(obj, feat1, sel1, feat2, sel2,
                    feature1_values = NULL, feature2_values = NULL) {
  # Coerce character(0) inputs (from un-selected radio buttons) to length-1
  # strings so downstream `||` chains don't blow up on logical(0).
  norm_str <- function(x, default = "") {
    if (is.null(x) || length(x) == 0) return(default)
    x <- as.character(x)[[1L]]
    if (is.na(x)) return(default)
    x
  }
  feat1 <- norm_str(feat1, "")
  sel1  <- norm_str(sel1,  "")
  feat2 <- norm_str(feat2, "")
  sel2  <- norm_str(sel2,  "none")

  m <- obj@meta.data

  # Handle primary feature (feature1)
  if (!is.null(feature1_values)) {
    m$feature1 <- align_to_cells(feature1_values, rownames(m))
    if (length(m$feature1) != nrow(m)) stop("Primary custom score length does not match the object.")
  } else if (!nzchar(feat1)) {
    m$feature1 <- 0
  } else if (identical(sel1, "matrisome gene") || identical(sel1, "any gene")) {
    feature1_data <- safe_get_assay_data(obj)
    if (!is.null(feature1_data) && feat1 %in% rownames(feature1_data)) {
      m$feature1 <- align_to_cells(feature1_data[feat1, ], rownames(m))
    } else {
      warning("Could not retrieve data for feature1")
      m$feature1 <- 0
    }
  } else if (identical(sel1, "matrisome signature")) {
    if (feat1 %in% colnames(obj@meta.data)) {
      m$feature1 <- obj@meta.data[, feat1]
    } else {
      m$feature1 <- 0
    }
  } else {
    m$feature1 <- 0
  }

  # Handle secondary feature (feature2)
  if (!is.null(feature2_values)) {
    m$feature2 <- align_to_cells(feature2_values, rownames(m))
    if (length(m$feature2) != nrow(m)) stop("Secondary custom score length does not match the object.")
  } else if (!nzchar(feat2) || identical(sel2, "none")) {
    m$feature2 <- 0
  } else if (identical(sel2, "matrisome gene") || identical(sel2, "any gene")) {
    feature2_data <- safe_get_assay_data(obj)
    if (!is.null(feature2_data) && feat2 %in% rownames(feature2_data)) {
      m$feature2 <- align_to_cells(feature2_data[feat2, ], rownames(m))
    } else {
      warning("Could not retrieve data for feature2")
      m$feature2 <- 0
    }
  } else if (identical(sel2, "matrisome signature")) {
    if (feat2 %in% colnames(obj@meta.data)) {
      m$feature2 <- obj@meta.data[, feat2]
    } else {
      m$feature2 <- 0
    }
  } else {
    m$feature2 <- 0
  }

  # Get spatial coordinates using GetTissueCoordinates for consistency
  coords <- align_spatial_coordinates(obj, cells = rownames(m))
  m <- m[rownames(coords), , drop = FALSE]
  coord_cols <- intersect(c("row", "col", "imagerow", "imagecol"), colnames(coords))
  m[coord_cols] <- coords[coord_cols]

  # Clean up negative values
  m$feature1[m$feature1 < 0] <- 0
  m$feature2[m$feature2 < 0] <- 0

  # Sparse LISA: dense n*n weight matrix OOMs the worker past ~5k spots.
  # The labels only depend on the sign of the spatial lag, so a sparse k-NN
  # adjacency gives identical High/Low quadrants at O(n*k) memory.
  run_lisa <- !identical(sel2, "none") && nzchar(feat2) &&
              isTRUE(var(m$feature2, na.rm = TRUE) > 0)

  if (run_lisa) {
    lisa_result <- tryCatch({
      # Pixel coords (imagerow/imagecol), not array indices: uploaded objects
      # can have character-typed or collinear row/col after VisiumV2->V1
      # fallback, which degenerates MERINGUE's Delaunay step.
      coords_df <- data.frame(x = as.numeric(m$imagerow), y = as.numeric(m$imagecol))
      Wij <- MERINGUE::getSpatialNeighbors(coords_df, filterDist = NA)

      rs <- Matrix::rowSums(Wij)
      rs[rs == 0] <- 1
      Wij <- Wij / rs

      x <- m$feature1
      y <- m$feature2
      x_scaled <- if (isTRUE(sd(x, na.rm = TRUE) > 0)) as.numeric(scale(x)) else rep(0, length(x))
      y_scaled <- if (isTRUE(sd(y, na.rm = TRUE) > 0)) as.numeric(scale(y)) else rep(0, length(y))
      x_scaled[is.na(x_scaled)] <- 0
      y_scaled[is.na(y_scaled)] <- 0

      lag_y <- as.numeric(Wij %*% y_scaled)

      # Dot-joined labels ("High.High", "Low.High") match the legacy
      # interaction() output that downstream plotting depends on.
      x_lab <- ifelse(x_scaled > 0, "High", "Low")
      y_lab <- ifelse(lag_y    > 0, "High", "Low")
      paste(x_lab, y_lab, sep = ".")
    }, error = function(e) {
      warning(sprintf("LISA calculation failed (%s); falling back to 'not.applicable'.",
                      conditionMessage(e)))
      NULL
    })

    m$LISA <- if (is.null(lisa_result)) "not.applicable" else lisa_result
  } else {
    m$LISA <- "not.applicable"
  }

  return(m)
}

#' Process matrisome expression for a gene category
#'
#' Calculates aggregate expression scores for a list of matrisome genes.
#' Handles gene synonym resolution when primary names are not found.
#'
#' @param gene_list Character vector of gene symbols
#' @param seurat_obj Seurat object with spatial data
#' @param name Name for the output score column
#' @param counts_matrix Pre-fetched expression matrix
#' @return Named list with robust_score, log_scaled_score, and gene_count
process_matrisome_expression <- function(
  gene_list, seurat_obj, name, counts_matrix, matrisome
) {
  counts <- counts_matrix
  if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0) {
    warning(paste("No counts matrix available for", name))
    return(NULL)
  }

  available_genes <- character(0)

  # For each gene, check primary name first, then synonyms if needed
  for(gene in gene_list) {
    if(gene %in% rownames(counts)) {
      available_genes <- c(available_genes, gene)
    } else {
      # Only check synonyms if primary gene not found
      gene_row <- matrisome[matrisome$gene == gene,]
      if(nrow(gene_row) == 1) {
        if(!is.na(gene_row$Synonyms)) {
          synonyms <- strsplit(gene_row$Synonyms, "\\||;|,|\\s+|/|-")[[1]]
          matching_synonym <- synonyms[synonyms %in% rownames(counts)]
          if(length(matching_synonym) > 0) {
            available_genes <- c(available_genes, matching_synonym[1])
          }
        }
      }
    }
  }

  if(length(available_genes) == 0) {
    warning(paste("No matching genes found for", name))
    return(NULL)
  }

  # Get spot coordinates
  coords <- align_spatial_coordinates(seurat_obj)

  # Matrix::colSums keeps the sparse path and avoids densifying dgCMatrix.
  base_counts <- Matrix::colSums(counts[available_genes, , drop = FALSE])
  base_counts <- align_to_cells(base_counts, rownames(coords))

  # Guard the degenerate min==max case so downstream colorRamp doesn't NaN.
  sums_log  <- log1p(base_counts)
  sl_min    <- min(sums_log, na.rm = TRUE)
  sl_max    <- max(sums_log, na.rm = TRUE)
  log_scaled <- if (isTRUE(sl_max > sl_min)) (sums_log - sl_min) / (sl_max - sl_min)
                else rep(0, length(sums_log))

  # IQR can be 0 when >50% of spots have zero counts for a small gene list.
  iqr_val <- IQR(base_counts, na.rm = TRUE)
  robust  <- if (isTRUE(iqr_val > 0)) (base_counts - median(base_counts, na.rm = TRUE)) / iqr_val
             else rep(0, length(base_counts))
  coord_rows <- if (!is.null(names(base_counts)) &&
                    all(names(base_counts) %in% rownames(coords))) {
    names(base_counts)
  } else {
    seq_len(min(length(base_counts), nrow(coords)))
  }

  # Create data frame - handle VisiumV1 vs VisiumV2 coordinate structures
  if (all(c("imagerow", "imagecol") %in% colnames(coords))) {
    df <- data.frame(
      spot = names(base_counts),
      base_counts = base_counts,
      log_scaled = log_scaled,
      robust = robust,
      imagerow = coords[coord_rows, "imagerow"],
      imagecol = coords[coord_rows, "imagecol"]
    )
  } else {
    # VisiumV2: use x/y as fallback
    df <- data.frame(
      spot = names(base_counts),
      base_counts = base_counts,
      log_scaled = log_scaled,
      robust = robust,
      imagerow = coords[coord_rows, "x"],
      imagecol = coords[coord_rows, "y"]
    )
  }
  return(df)
}


#' Compute per-category matrisome plot data (no plotly objects).
#' Plot construction is deferred to build_matrisome_plotly() at render time
#' so only the visible category materializes plotly objects.
#' Pass coords_cache to avoid re-running GetTissueCoordinates per call.
compute_matrisome_plot_data <- function(seurat_obj, display_name, internal_name, type,
                                        annotation_col, coords_cache = NULL) {

  robust_feature_name <- paste0(internal_name, "_robust_score")
  log_scaled_feature_name <- paste0(internal_name, "_log_scaled_score")

  # --- Robust Validation ---
  if (!robust_feature_name %in% colnames(seurat_obj@meta.data) ||
      !log_scaled_feature_name %in% colnames(seurat_obj@meta.data)) {
    return(NULL)
  }
  feature_values <- seurat_obj@meta.data[[robust_feature_name]]
  valid_indices <- !is.na(feature_values)
  non_na_count <- sum(valid_indices)
  total_count <- length(feature_values)
  if (non_na_count < 20 || non_na_count < (total_count * 0.01)) {
    return(NULL)
  }

  # Align coords + FetchData by spot identity before Plotly sees them. D3 does
  # the same via barcode merge; relying on cbind row position can desync hovers.
  fetched <- FetchData(seurat_obj, vars = c(robust_feature_name, log_scaled_feature_name, annotation_col))
  coords  <- align_spatial_coordinates(seurat_obj, coords_cache, cells = rownames(fetched))
  fetched <- fetched[rownames(coords), , drop = FALSE]
  df      <- cbind(coords, fetched)
  df[[annotation_col]] <- as.character(df[[annotation_col]])

  seurat_spectral_palette <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
  seurat_color_ramp <- colorRamp(seurat_spectral_palette)

  to_hex <- function(values) {
    vmin <- suppressWarnings(min(values, na.rm = TRUE))
    vmax <- suppressWarnings(max(values, na.rm = TRUE))
    if (!isTRUE(is.finite(vmin)) || !isTRUE(is.finite(vmax)) || vmax == vmin) {
      base <- seurat_color_ramp(0)
      return(rep(rgb(base[1, 1], base[1, 2], base[1, 3], maxColorValue = 255),
                 length(values)))
    }
    normalized <- (values - vmin) / (vmax - vmin)
    rgbm <- seurat_color_ramp(normalized)
    bad  <- is.na(rgbm[, 1]) | is.na(rgbm[, 2]) | is.na(rgbm[, 3])
    out  <- rgb(rgbm[, 1], rgbm[, 2], rgbm[, 3], maxColorValue = 255)
    out[bad] <- "#808080"
    out
  }

  robust_values  <- df[[robust_feature_name]]
  hotspot_values <- df[[log_scaled_feature_name]]

  list(
    robust_feature      = robust_feature_name,
    log_scaled_feature  = log_scaled_feature_name,
    annotation_col      = annotation_col,
    df                  = df,
    robust_colors       = to_hex(robust_values),
    hotspot_colors      = to_hex(hotspot_values),
    robust_min          = suppressWarnings(min(robust_values,  na.rm = TRUE)),
    robust_max          = suppressWarnings(max(robust_values,  na.rm = TRUE)),
    hotspot_min         = suppressWarnings(min(hotspot_values, na.rm = TRUE)),
    hotspot_max         = suppressWarnings(max(hotspot_values, na.rm = TRUE)),
    palette             = seurat_spectral_palette
  )
}

#' Build the plotly object for one matrisome category. Called lazily from
#' renderPlotly() so only the visible category materializes a plotly.
#' plot_data is an entry from rv$matrisome_results; NULL-safe.
build_matrisome_plotly <- function(plot_data, kind = c("spatial", "hotspot")) {
  if (is.null(plot_data)) return(NULL)
  kind <- match.arg(kind)

  if (kind == "spatial") {
    colors       <- plot_data$robust_colors
    vmin         <- plot_data$robust_min
    vmax         <- plot_data$robust_max
    feature_name <- plot_data$robust_feature
  } else {
    colors       <- plot_data$hotspot_colors
    vmin         <- plot_data$hotspot_min
    vmax         <- plot_data$hotspot_max
    feature_name <- plot_data$log_scaled_feature
  }

  ticks     <- pretty(c(vmin, vmax), n = 5)
  tick_text <- round(ticks, 1)
  df             <- plot_data$df
  annotation_col <- plot_data$annotation_col

  plot_ly(df, x = ~imagecol, y = ~imagerow, type = 'scatter', mode = 'markers',
          marker = list(
            color = colors,
            colorbar = list(
              title = "",
              tickvals = ticks,
              ticktext = tick_text,
              len = 0.5, thickness = 15
            ),
            cmin = vmin,
            cmax = vmax,
            colorscale = plot_data$palette,
            showscale = FALSE,
            size = 5
          ),
          text = as.formula(paste0(
            "~paste('Value:', round(`", feature_name, "`, 3), ",
            "'<br>Annotation:', `", annotation_col, "`)"
          )),
          hoverinfo = 'text') %>%
    layout(xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
           yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE,
                        scaleanchor = "x", scaleratio = 1, autorange = "reversed"),
           plot_bgcolor = 'rgba(0,0,0,0)', paper_bgcolor = 'rgba(0,0,0,0)') %>%
    config(responsive = TRUE)
}


#' Calculate spatial statistics for feature correlation.
calcspatstat <- function(obj,annots,minexp,sels){
  kz <- obj
  kz$id <- kz[,colnames(kz)%in%as.character(annots)]  # Added as.character()

  kz <- kz[order(-kz$feature1),]
  vvv <- kz$feature1[kz$feature1>0]
  q <- ntile(vvv,10)
  kz$q <- c(q,rep(0,length(kz$feature1)-length(q)))
  kz$feature1[kz$q < minexp/10] <- 0

  kz <- kz[order(-kz$feature2),]
  vvv <- kz$feature2[kz$feature2>0]
  q <- ntile(vvv,10)
  kz$q2 <- c(q,rep(0,length(kz$feature2)-length(q)))
  kz$feature2[kz$q < minexp/10] <- 0

  kz <- kz[kz$id %in% sels,]

  mt <- FindSpatiallyVariableFeatures(
    t(cbind(
      kz[,colnames(kz)%in%c("feature1","feature2")],
      ifelse(kz$feature1>0 & kz$feature2>0,1,0)
    )),
    spatial.location = kz[,colnames(kz)%in%c("row","col")],
    selection.method = "moransi"
  )

  # Simplified results creation - removed intermediate variables
  df <- data.frame(a=c(
    round(mt[1,1],3),
    round(mt[2,1],3),
    round(mt[3,1],3)
  ))

  return(df)
}

#' Corner colours for the co-expression blend, keyed by palette name.
coex_blend_colors <- function(palette) {
  if (palette == "Classic") {
    list(bottom_left = "#d3d3d3", bottom_right = "#FF0000", top_left = "#00FF00", top_right = "#FFFF00")
  } else if (palette == "Vibrant") {
    list(bottom_left = "white", bottom_right = "orange", top_left = "#0000FF", top_right = "#FF0000")
  } else {
    list(bottom_left = "#d3d3d3", bottom_right = "#FF00FF", top_left = "#00FF00", top_right = "#FFFFFF")
  }
}

#' Create blended spatial plot showing two features.
#' Credit: https://github.com/george-hall-ucl/SpatialFeaturePlotBlend
SpatialFeaturePlotBlend <- function(object, features, combine = TRUE,
                                    feature_1_alt_name = NULL,
                                    feature_2_alt_name = NULL, assay = NULL,
                                    bottom_left = "#000000",
                                    bottom_right = "#FF0000",
                                    top_left = "#00FF00",
                                    top_right = "#FFFF00",
                                    use_seurat_backend = FALSE,
                                    fp_extra_arguments = list(),
                                    sfp_extra_arguments = list())  {

  # Generate a grid of RGB color values given the requested corner colours.
  gen_color_grid <- function(side_length, bottom_left, bottom_right,
                             top_left, top_right) {

    grad_gen <- function(start, end, n = side_length) {
      colfunc <- colorRampPalette(c(start, end))
      return(colfunc(n))
    }

    # x_y = "x to y"; "bl" = "bottom left", etc
    bl_tl <- grad_gen(bottom_left, bottom_right)
    br_tr <- grad_gen(top_left, top_right)

    l <- lapply(seq_len(length(bl_tl)),
                function(i) {
                  start <- bl_tl[i]
                  end <- br_tr[i]
                  new_grad <- grad_gen(start, end)
                })

    return(t(matrix(unlist(l), ncol = side_length, nrow = side_length)))
  }
  custom_color_SpatialDimPlot <- function(cells_obj, image_name,
                                          new_md_column_name,
                                          colors_per_spot, ...) {
    cells_obj[[new_md_column_name]] <- colors_per_spot
    names(colors_per_spot) <- as.character(colors_per_spot)

    p <- SpatialDimPlot(cells_obj, new_md_column_name,
                        cols = colors_per_spot, images = image_name, ...) +
      ggtitle(new_md_column_name) +
      blend_plot_theme
    return(p)
  }

  extract_colors_from_ggplot <- function(p) {
    built <- ggplot_build(p)$data[[1]]
    if (!is.na(built[1, "fill"])) {
      col_to_use <- "fill"
    } else {
      col_to_use <- "colour"
    }
    return(built[, col_to_use])
  }

  if (length(features) != 2) {
    stop(paste(c("Incorrect number of features. ",
                 "Requires two features, received ",
                 length(features))))
  }

  if (!is.null(assay)) {
    DefaultAssay(object) <- assay
  }

  if (length(fp_extra_arguments) > 0) {
    use_seurat_backend <- TRUE
  }

  blend_plot_theme <- theme(legend.position = "none",
                            plot.title = element_text(hjust = 0.5))

  plot_list_outer <- list()

  for (i in Images(object)) {
    cell_barcodes <- Seurat:::CellsByImage(object, images = i,
                                           unlist = TRUE)
    cells_obj_sub <- object[, cell_barcodes]
    images_sub_list <- list(object[[i]])
    names(images_sub_list) <- i
    cells_obj_sub@images <- images_sub_list
    if (!use_seurat_backend) {
      plot_list <- lapply(features,
                          function(feature) {
                            max_color <- ifelse(feature == features[1],
                                                bottom_right, top_left)
                            SpatialFeaturePlot(object, feature,
                                               images = i,
                                               sfp_extra_arguments) +
                              scale_fill_gradient(low = bottom_left,
                                                  high = max_color) +
                              ggtitle(feature) +
                              blend_plot_theme
                          })
      colors_list <- lapply(plot_list, extract_colors_from_ggplot)

      # Now construct the blended plot
      dat <- FetchData(cells_obj_sub, features)
      side_length <- 100
      col_grid <- gen_color_grid(side_length, bottom_left, bottom_right,
                                 top_left, top_right)
      dat_norm <- apply(dat, 2,
                        function(x) {
                          round((side_length - 1) * x / max(x)) + 1
                        })
      colors_list[[3]] <- sapply(seq_len(nrow(dat_norm)),
                                 function(x) {
                                   col_grid[dat_norm[x, 1],
                                            dat_norm[x, 2]]
                                 })
      legend_grid <- expand.grid(seq(from = min(dat[, features[1]]),
                                     to = max(dat[, features[1]]),
                                     length.out = side_length),
                                 seq(from = min(dat[, features[2]]),
                                     to = max(dat[, features[2]]),
                                     length.out = side_length))
      colnames(legend_grid) <- features
      legend_colors <- c(col_grid)
      legend_grid$color <- legend_colors
      names(legend_colors) <- legend_colors

      legend <- ggplot(legend_grid,
                       aes(x = .data[[features[1]]],
                           y = .data[[features[2]]],
                           color = color)) +
        geom_point(shape = 15, size = 1.9) +
        scale_color_manual(values = legend_colors) +
        coord_cartesian(expand = FALSE) +
        theme(legend.position = "none", aspect.ratio = 1,
              panel.background = element_blank(),
              axis.text.x = element_text(angle = 45,
                                         hjust = 1)) +
        xlab(ifelse(is.null(feature_1_alt_name),
                    features[1], feature_1_alt_name)) +
        ylab(ifelse(is.null(feature_2_alt_name),
                    features[2], feature_2_alt_name))
    } else {
      if (top_right != "#FFFF00") {
        warning(paste("Cannot alter color in top right corner when",
                      "use_seurat_backend is TRUE"))
      }
      vis_reduc <- cells_obj_sub@images[[i]]@coordinates[, c(3, 2)]
      colnames(vis_reduc) <- c("vis_1", "vis_2")
      vis_reduc$vis_2 <- -1 * vis_reduc$vis_2
      vis_reduc_mat <- as.matrix(vis_reduc)
      vis_reduc_obj <- CreateDimReducObject(embeddings = vis_reduc_mat,
                                            key = "vis_")
      cells_obj_sub@reductions$vis <- vis_reduc_obj
      seurat_fp <- do.call(FeaturePlot, c(list(object = cells_obj_sub,
                                               features = features,
                                               reduction = "vis",
                                               blend = TRUE,
                                               cols = c(bottom_left,
                                                        bottom_right,
                                                        top_left),
                                               combine = FALSE),
                                          fp_extra_arguments))
      colors_list <- lapply(seurat_fp[1:3], extract_colors_from_ggplot)
      legend <- seurat_fp[[4]]
    }

    names(colors_list) <- c(features, paste0(features[1], "_", features[2]))
    plot_list <- lapply(names(colors_list),
                        function(x) {
                          do.call(custom_color_SpatialDimPlot,
                                  c(list(cells_obj = cells_obj_sub, i,
                                         x, colors_list[[x]]),
                                    sfp_extra_arguments))
                        })

    plot_list[[4]] <- wrap_plots(ggplot() + theme_void(), legend,
                                 ggplot() + theme_void(), ncol = 1,
                                 heights = c(0.2, 0.6, 0.2))

    plot_list_outer[[i]] <- plot_list
  }

  if (combine == FALSE) {
    return(plot_list_outer)
  } else {
    plot_list_outer <- lapply(plot_list_outer,
                              function(p) {
                                wrap_plots(p, nrow = 1,
                                           widths = c(0.28, 0.28,
                                                      0.28, 0.16))
                              })
    p <- wrap_plots(plot_list_outer, ncol = 1)

    return(p)
  }
}

#' Generate consistent color palette for annotations.
create_custom_color_map <- function(levels_vector) {
  # Ensure input is unique characters and sorted for consistency
  all_levels <- sort(unique(as.character(levels_vector)))

  # Define the neutral color (hex code for JavaScript/CSS compatibility)
  not_assigned_color <- "#7f7f7f"  # This is grey50

  # Check if "not.assigned" exists
  if ("not.assigned" %in% all_levels) {
    # Separate the "real" annotations from "not.assigned"
    real_levels <- all_levels[all_levels != "not.assigned"]
    num_real_levels <- length(real_levels)

    # Generate palette ONLY for the real levels
    if (num_real_levels > 0) {

      # Check if the number of levels exceeds the palette's limit
      if (num_real_levels > 32) {
        # If too many, use a continuous palette (viridis) to generate colors
        palette_func <- colorRampPalette(viridis::viridis(num_real_levels))
        real_palette <- palette_func(num_real_levels)
      } else {
        # If within the limit, use the highly distinct "glasbey" palette
        real_palette <- DiscretePalette_scCustomize(num_colors = num_real_levels, palette = "glasbey")
      }

      real_color_map <- setNames(real_palette, real_levels)

    } else {
      real_color_map <- c() # Handle case with only "not.assigned"
    }

    # Create the named vector for the "not.assigned" color
    not_assigned_map <- setNames(not_assigned_color, "not.assigned")

    # Combine them. The "not.assigned" part is added to the main map.
    final_color_map <- c(real_color_map, not_assigned_map)

  } else {
    # If "not.assigned" is not present, just generate the palette as usual
    num_all_levels <- length(all_levels)
    if (num_all_levels > 0) {

      if (num_all_levels > 32) {
        palette_func <- colorRampPalette(viridis::viridis(num_all_levels))
        palette <- palette_func(num_all_levels)
      } else {
        palette <- DiscretePalette_scCustomize(num_colors = num_all_levels, palette = "glasbey")
      }

      final_color_map <- setNames(palette, all_levels)

    } else {
      final_color_map <- c() # Handle case of no levels
    }
  }

  return(final_color_map)
}

#' Print hex codes for cell types to console.
#'
#' @param seurat_obj A Seurat object.
#' @param annotation_col Name of the metadata column containing cell type labels.
#' @return Invisibly returns the named vector of hex codes.
print_celltype_colors <- function(seurat_obj, annotation_col = "cell_type") {
  if (!annotation_col %in% colnames(seurat_obj@meta.data)) {
    stop("Column '", annotation_col, "' not found in metadata. Available columns:\n",
         paste(colnames(seurat_obj@meta.data), collapse = ", "))
  }

  levels_vector <- seurat_obj@meta.data[[annotation_col]]
  color_map <- create_custom_color_map(levels_vector)

  invisible(color_map)
}

# ============================================================================ #
#                  NEW: ON-THE-FLY PREPROCESSING ENGINE                        #
# ============================================================================ #

# Source the required ScType scripts safely
tryCatch({
  source(extdata_path("gene_sets_prepare.R"), local = TRUE)
  source(extdata_path("sctype_score_.R"), local = TRUE)
}, error = function(e) {
  warning("Could not source ScType helper scripts. Annotation will fail.")
})

#' Run ScType to generate ecm_domain_annotation
sctype_annotate_ecm <- function(seurat_obj) {
  if ("ecm_domain_annotation" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj@meta.data$ecm_domain_annotation <- NULL
  }
  if (!file.exists(extdata_path("ECM_domains_transformed4ScType.xlsx"))) {
    warning("ECM domain database not found. Skipping annotation.")
    return(seurat_obj)
  }
  if (!"SCT" %in% SeuratObject::Assays(seurat_obj)) {
    warning("SCT assay not found. Skipping ECM niche annotation.")
    return(seurat_obj)
  }

  gs_list <- gene_sets_prepare(extdata_path("ECM_domains_transformed4ScType.xlsx"), "ECM")
  ecm_classes <- c("Interstitial_ECM", "Basement_ECM")
  if (!all(ecm_classes %in% names(gs_list$gs_positive))) {
    warning("Required ECM classes were not found in the ScType marker database.")
    return(seurat_obj)
  }
  gs_list$gs_positive <- gs_list$gs_positive[ecm_classes]
  gs_list$gs_negative <- gs_list$gs_negative[ecm_classes]

  marker_genes <- unique(c(unlist(gs_list$gs_positive), unlist(gs_list$gs_negative)))
  marker_genes <- marker_genes[!is.na(marker_genes) & nzchar(marker_genes)]
  scale_genes <- rownames(seurat_obj[["SCT"]]@scale.data)
  missing_markers <- intersect(setdiff(marker_genes, scale_genes), rownames(seurat_obj[["SCT"]]))

  if (length(missing_markers) > 0 && length(seurat_obj[["SCT"]]@SCTModel.list) > 0) {
    umi_assay <- tryCatch(Seurat::SCTResults(seurat_obj[["SCT"]], slot = "umi.assay")[[1]],
                          error = function(e) NULL)
    if (!is.null(umi_assay) && umi_assay %in% SeuratObject::Assays(seurat_obj)) {
      seurat_obj <- tryCatch(
        Seurat::GetResidual(seurat_obj, features = missing_markers, assay = "SCT",
                            umi.assay = umi_assay, verbose = FALSE),
        error = function(e) seurat_obj
      )
    }
  }

  scale_data <- seurat_obj[["SCT"]]@scale.data
  marker_data <- scale_data[intersect(marker_genes, rownames(scale_data)), , drop = FALSE]
  usable_markers <- character(0)
  if (nrow(marker_data) > 0) {
    usable_markers <- rownames(marker_data)[apply(marker_data, 1, function(x) {
      x <- x[is.finite(x)]
      length(x) > 1 && max(x) > min(x)
    })]
  }
  if (any(lengths(lapply(gs_list$gs_positive, intersect, y = usable_markers)) == 0)) {
    warning("ECM niche annotation unavailable: a required class has no usable markers.")
    return(seurat_obj)
  }
  marker_data <- marker_data[usable_markers, , drop = FALSE]

  es.max <- sctype_score(scRNAseqData = marker_data, scaled = TRUE,
                         gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

  # Get top score for each spot
  spot_results <- do.call("rbind", lapply(rownames(seurat_obj@meta.data), function(coord) {
    es.max.subset <- es.max[, coord, drop = FALSE]
    es.max.coord <- sort(rowSums(as.matrix(es.max.subset)), decreasing = TRUE)
    head(data.frame(spot = coord, type = names(es.max.coord), scores = es.max.coord), 1)
  }))

  spot_results$type[as.numeric(as.character(spot_results$scores)) <= 0] <- "not.assigned"

  # Add to metadata, ensuring correct order
  seurat_obj$ecm_domain_annotation <- spot_results$type[match(rownames(seurat_obj@meta.data), spot_results$spot)]

  return(seurat_obj)
}

#' Translate gene symbols to current HGNC and filter to those present in data
translate_gene_signatures <- function(seurat_obj, signature_list) {
  seurat_genes <- rownames(seurat_obj)
  translated_signatures <- list()

  for (sig_name in names(signature_list)) {
    original_genes <- signature_list[[sig_name]]
    result <- update_hgnc_symbols(original_genes, hgnc_long_data)
    corrected_genes <- result$Output_Features
    valid_genes <- unique(corrected_genes[corrected_genes %in% seurat_genes])
    if (length(valid_genes) > 0) {
      translated_signatures[[sig_name]] <- valid_genes
    }
  }
  return(translated_signatures)
}

#' Add Matrisome Feature Scores (for Feature Selection dropdown)
add_matrisome_feature_scores <- function(seurat_obj, signature_list) {
  translated_sigs <- translate_gene_signatures(seurat_obj, signature_list)

  # UCell requires list of lists for multiple signatures
  seurat_obj <- AddModuleScore_UCell(
    seurat_obj,
    features = translated_sigs,
    name = "_score" # Appends this suffix to each signature name
  )

  # Clean up column names (e.g., "collagens_score1" -> "collagens")
  new_colnames <- gsub("_score1$", "", colnames(seurat_obj@meta.data))
  colnames(seurat_obj@meta.data) <- new_colnames

  # Create display-named columns to match UI dropdown names
  # Maps UCell column names (with _score suffix) to display names used in Feature Selection
  display_name_map <- c(
    # Main categories
    "ecm_glycoproteins_score" = "ECM Glycoproteins",
    "collagens_score" = "Collagens",
    "proteoglycans_score" = "Proteoglycans",
    "ecm_regulators_score" = "ECM Regulators",
    "secreted_factors_score" = "Secreted Factors",
    "ecm-affiliated_proteins_score" = "ECM-affiliated Proteins",
    # Subcategories
    "basement_membrane_score" = "Basement membrane",
    "hemostasis_score" = "Hemostasis",
    "elastic_fibers_score" = "Elastic fibers",
    "growth_factor-binding_score" = "Growth-factor binding",
    # Gene families
    "laminin_-_basement_membrane_score" = "Laminins",
    "matricellular_score" = "Matricellular proteins",
    "syndecan_score" = "Syndecans",
    "glypican_score" = "Glypicans",
    # Others
    "annexin_score" = "Annexins",
    "cathepsin_score" = "Cathepsins",
    "ccn_family_score" = "CCNs",
    "cystatin_score" = "Cystatins",
    "facit_score" = "FACITs",
    "fibulin_score" = "Fibulins",
    "galectin_score" = "Galectins",
    "mucin_score" = "Mucins",
    "plexin_score" = "Plexins",
    "semaphorin_score" = "Semaphorins"
  )

  # Copy UCell columns to display-named columns
  for (internal_name in names(display_name_map)) {
    if (internal_name %in% colnames(seurat_obj@meta.data)) {
      display_name <- display_name_map[[internal_name]]
      seurat_obj@meta.data[[display_name]] <- seurat_obj@meta.data[[internal_name]]
    }
  }

  return(seurat_obj)
}

#' Master orchestrator for preparing an uploaded Seurat object
prepare_uploaded_object <- function(seurat_obj, mat_feat_sigs, ecm_ucell_sigs) {

  # Fix metadata/count matrix cell name mismatch (common with subsetted objects)
  if (!identical(colnames(seurat_obj), rownames(seurat_obj@meta.data))) {
    common_cells <- intersect(colnames(seurat_obj), rownames(seurat_obj@meta.data))
    if (length(common_cells) == 0) stop("No matching cell names between metadata and count matrix.")
    if (length(common_cells) < ncol(seurat_obj)) {
      seurat_obj <- subset(seurat_obj, cells = common_cells)
    } else {
      seurat_obj@meta.data <- seurat_obj@meta.data[colnames(seurat_obj), , drop = FALSE]
    }
  }

  # SCTransform if ScType has neither residuals nor a fitted SCT model
  has_sct <- "SCT" %in% SeuratObject::Assays(seurat_obj)
  has_scale_data <- has_sct && nrow(seurat_obj[["SCT"]]@scale.data) > 0
  has_sct_model <- has_sct && length(seurat_obj[["SCT"]]@SCTModel.list) > 0

  if (!has_scale_data && !has_sct_model) {
    showNotification("SCT data not found. Running SCTransform...", type = "message", duration = 8)
    # Use default assay dynamically (handles both native Seurat "Spatial" and converted SPE "originalexp"/"RNA")
    seurat_obj <- SCTransform(seurat_obj, assay = DefaultAssay(seurat_obj), verbose = FALSE)
  }
  if ("SCT" %in% SeuratObject::Assays(seurat_obj)) {
    DefaultAssay(seurat_obj) <- "SCT"
  }

  # --- Step 2: Check & Add Matrisome Feature Signatures ---
  if (!"collagens" %in% colnames(seurat_obj@meta.data)) {
    showNotification("Matrisome feature scores not found. Calculating now...", type = "message", duration = 8)
    seurat_obj <- add_matrisome_feature_scores(seurat_obj, mat_feat_sigs)
  }

  # --- Step 3: Add ECM Domain Annotation ---
  seurat_obj <- sctype_annotate_ecm(seurat_obj)
  if (!"ecm_domain_annotation" %in% colnames(seurat_obj@meta.data)) {
    showNotification(
      "ECM niche annotations are unavailable because this dataset does not contain the required marker genes.",
      type = "warning", duration = 10
    )
  }

  # --- Step 4: Check & Add ECM Niche Signatures ---
  if (!"Interstitial_UCell" %in% colnames(seurat_obj@meta.data)) {
    showNotification("ECM niche scores not found. Running UCell...", type = "message", duration = 8)
    translated_sigs <- translate_gene_signatures(seurat_obj, ecm_ucell_sigs)
    if (length(translated_sigs) > 0) {
      seurat_obj <- AddModuleScore_UCell(seurat_obj, features = translated_sigs)
    }
  }

  showNotification("Preprocessing complete!", type = "message", duration = 5)
  return(seurat_obj)
}

# --------------------------------------------------------------------------- #
# Spatial Statistics & Visualization Helpers
# --------------------------------------------------------------------------- #

#' Permutation test for spatial cross-correlation (bivariate Moran's I).
#'
#' Sparse, allocation-free reformulation of MERINGUE::spatialCrossCor
#' (cv = dx'Wdy + dy'Wdx) so the label-permutation null is cheap to build
#' and avoids the dense N x N outer products of the original implementation.
#'
#' @param x,y Named feature vectors aligned to the spatial weight matrix.
#' @param weight Spatial neighbour weight matrix (from getSpatialNeighbors).
#' @param n Number of label permutations for the null distribution.
#' @param seed Seed for reproducible permutation p-values.
#' @return list(cor = observed statistic, pval = two-sided permutation p-value).
spatial_cross_cor_test <- function(x, y, weight, n = 999, seed = 1) {
  common <- intersect(names(x), rownames(weight))
  if (length(common) < 3) return(list(cor = NA_real_, pval = NA_real_))

  W <- methods::as(weight[common, common], "CsparseMatrix")
  x <- x[common]; y <- y[common]
  rs <- Matrix::rowSums(W); rs[rs == 0] <- 1
  W <- Matrix::Diagonal(x = 1 / rs) %*% W

  N <- length(x); Wsum <- sum(W)
  dy <- y - mean(y); v_y <- sum(dy^2); Wdy <- as.numeric(W %*% dy)

  scc <- function(dx) {
    v <- sqrt(sum(dx^2) * v_y)
    if (v == 0) return(0)
    cv <- sum(dx * Wdy) + sum(dy * as.numeric(W %*% dx))
    (N / Wsum) * (cv / v) / 2
  }

  observed <- scc(x - mean(x))

  set.seed(seed)
  perm <- vapply(seq_len(n), function(i) {
    xp <- x[sample.int(N)]
    scc(xp - mean(xp))
  }, numeric(1))

  list(cor = observed, pval = sum(abs(c(perm, observed)) >= abs(observed)) / (n + 1))
}

#' Calculate Spatial Autocorrelation (Moran's I) for Matrisome Scores
#'
#' @param seurat_obj A Seurat object containing matrisome scores in metadata.
#' @return A data frame with Moran's I and p-value for each score category.
calculate_autocorrelation <- function(seurat_obj) {
  tryCatch({
    coords <- GetTissueCoordinates(seurat_obj)
    if (nrow(coords) < 3) return(NULL)

    # Handle VisiumV2 coordinate column names (x/y instead of imagecol/imagerow)
    if(!all(c("imagerow","imagecol")%in%colnames(coords))){
      coords$imagerow <- coords$x
      coords$imagecol <- coords$y
    }
    coords <- coords[,colnames(coords)%in%c("imagerow","imagecol")]

    weight_matrix <- getSpatialNeighbors(coords)
    score_cols <- grep("_robust_score", colnames(seurat_obj@meta.data), value = TRUE)

    results_list <- lapply(score_cols, function(col_name) {
      scores <- seurat_obj@meta.data[[col_name]]
      scores <- scores[is.finite(scores)]

      if (length(scores) < 3 || var(scores) == 0) {
        return(NULL)
      }

      names(scores) <- rownames(seurat_obj@meta.data)[is.finite(seurat_obj@meta.data[[col_name]])]
      test_result <- moranTest(scores, weight_matrix)

      data.frame(
        category = col_name,
        morans_I = test_result["observed"],
        p_value = test_result["p.value"],
        stringsAsFactors = FALSE
      )
    })

    results_df <- do.call(rbind, results_list[!sapply(results_list, is.null)])
    return(results_df)

  }, error = function(e) {
    warning(paste("Error in autocorrelation calculation:", e$message))
    return(NULL)
  })
}

#' Create an Interactive Spatial Feature Plot using Plotly
#'
#' @param data A data frame with coordinates, feature values, and annotations.
#' @param feature_col The name of the column containing the feature values.
#' @param feature_name The display name for the feature (for tooltips).
#' @param annotation_vector A vector of annotations for hover text.
#' @return A plotly object.
create_feature_plot <- function(data, feature_col, feature_name, annotation_vector,
                                colour_values = NULL, relative_contrast = FALSE) {
  if (is.null(data) || !feature_col %in% colnames(data)) {
    return(plot_ly(type = 'scatter', mode = 'markers') %>%
             add_annotations(
               text = "No data available",
               showarrow = FALSE,
               font = list(size = 16)
             ))
  }

  data <- as.data.frame(data)
  if (!is.null(names(annotation_vector)) && all(rownames(data) %in% names(annotation_vector))) {
    annotation_vector <- annotation_vector[rownames(data)]
  }
  if (length(annotation_vector) == nrow(data)) {
    data$Annotation <- as.character(annotation_vector)
  } else if ("Annotation" %in% colnames(data)) {
    data$Annotation <- as.character(data$Annotation)
  } else {
    data$Annotation <- NA_character_
  }
  if (!is.null(colour_values) && !is.null(names(colour_values)) &&
      all(rownames(data) %in% names(colour_values))) {
    colour_values <- colour_values[rownames(data)]
  }
  if (is.null(colour_values) || length(colour_values) != nrow(data)) {
    colour_values <- data[[feature_col]]
    relative_contrast <- FALSE
  }

  data$.hover_text <- paste(
    feature_name, ":", round(data[[feature_col]], 3),
    "<br>Annotation:", data$Annotation
  )
  if (relative_contrast) {
    data$.hover_text <- paste0(
      data$.hover_text,
      "<br>Relative percentile: ", round(100 * colour_values, 1), "%"
    )
  }

  colors <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  max_val <- if (relative_contrast) 1 else max(colour_values, na.rm = TRUE)
  colorbar_title <- if (relative_contrast) {
    "UCell relative<br>contrast"
  } else if (identical(feature_name, "Custom matrisome list")) {
    "Custom matrisome<br>list"
  } else {
    feature_name
  }

  plot_ly(data,
          x = ~imagecol,
          y = ~imagerow,
          type = 'scatter',
          mode = 'markers',
          marker = list(
            size = 5,
            color = colour_values,
            colorscale = list(
              list(0, colors[1]),
              list(0.5, colors[6]),
              list(1, colors[11])
            ),
            showscale = TRUE,
            cmin = 0,
            cmax = max_val,
            colorbar = list(
              title = list(text = ""),
              thickness = 15, len = 0.5, x = 0.86, xanchor = "left", y = 0.5
            )
          ),
          text = ~.hover_text,
          hoverinfo = 'text') %>%
    layout(
      xaxis = list(title = "", domain = c(0, 0.78), constrain = "domain",
                   showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
      yaxis = list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE,
                   scaleanchor = "x", scaleratio = 1, constrain = "domain",
                   autorange = "reversed"),
      annotations = list(list(
        x = 0.89, y = 0.745, xref = "paper", yref = "paper",
        text = colorbar_title, showarrow = FALSE,
        xanchor = "center", yanchor = "bottom", align = "left",
        font = list(size = 14)
      )),
      margin = list(l = 20, r = 20, t = 20, b = 20),
      plot_bgcolor = 'rgba(0,0,0,0)', paper_bgcolor = 'rgba(0,0,0,0)'
    )
}

#' Create Primary Feature Expression Plot (Violin/Box)
#'
#' @param data_df The analysis results data frame.
#' @param feature_name The display name of the feature.
#' @param stable_color_map A named vector of colors for annotations.
#' @param selected_types A vector of the currently selected annotations to display.
#' @param cleaned_annot_vector The full, cleaned annotation vector for all spots.
#' @param use_violin Boolean, TRUE for violin plot, FALSE for box plot.
#' @param independent_y Boolean, TRUE for independent y-axis scaling.
#' @return A plotly object.
p3a <- function(data_df, feature_name, stable_color_map, selected_types, cleaned_annot_vector, use_violin, independent_y) {
  # Validate feature1 exists and has data
  if (!"feature1" %in% colnames(data_df) || all(is.na(data_df$feature1)) || nrow(data_df) == 0) {
    return(plot_ly(type = 'scatter', mode = 'markers') %>%
             layout(title = paste("No data available for", feature_name)))
  }

  kz <- data.frame(clust = cleaned_annot_vector, value = data_df$feature1)
  kz <- kz[kz$clust %in% selected_types, ]
  if (nrow(kz) == 0) return(plot_ly(type = 'scatter', mode = 'markers') %>% layout(title = "No data for selected annotations"))

  kz$clust <- factor(kz$clust, levels = selected_types)
  kz <- kz %>% group_by(clust) %>% filter(n() > 0) %>% ungroup()

  # Only calculate p-value if there are 2 or more groups to compare
  p_value <- if (nlevels(kz$clust) >= 2) {
    tryCatch(summary(aov(value ~ clust, data = kz))[[1]][["Pr(>F)"]][1], error = function(e) NA)
  } else {
    NA
  }

  plot_type <- if (use_violin) 'violin' else 'box'
  p <- plot_ly(kz, x = ~clust, y = ~value, type = plot_type, color = ~clust, colors = stable_color_map)

  p <- p %>% layout(
    title = list(text = feature_name),
    xaxis = list(title = '', showticklabels = FALSE),
    yaxis = list(title = 'Expression', autorange = independent_y),
    showlegend = FALSE,
    annotations = list(x = 0.5, y = 1.05, xref = 'paper', yref = 'paper',
                       text = if (is.na(p_value)) "Single group selected" else paste("ANOVA p-value:", format.pval(p_value, digits = 4)),
                       showarrow = FALSE)
  ) %>% config(displayModeBar = FALSE)

  if (!independent_y) {
    shared_range <- c(0, max(c(data_df$feature1, data_df$feature2), na.rm = TRUE))
    p <- p %>% layout(yaxis = list(title = 'Expression', range = shared_range))
  }
  return(p)
}

#' Create Secondary Feature Expression Plot (Violin/Box)
#'
#' @description This function is largely identical to p3a but targets 'feature2'.
p3b <- function(data_df, feature_name, stable_color_map, selected_types, cleaned_annot_vector, use_violin, independent_y) {
  # Validate feature2 exists and has data
  if (!"feature2" %in% colnames(data_df) || all(is.na(data_df$feature2)) || nrow(data_df) == 0) {
    return(plot_ly(type = 'scatter', mode = 'markers') %>%
             layout(title = paste("No data available for", feature_name)))
  }

  kz <- data.frame(clust = cleaned_annot_vector, value = data_df$feature2)
  kz <- kz[kz$clust %in% selected_types, ]
  if (nrow(kz) == 0) return(plot_ly(type = 'scatter', mode = 'markers') %>% layout(title = "No data for selected annotations"))

  kz$clust <- factor(kz$clust, levels = selected_types)
  kz <- kz %>% group_by(clust) %>% filter(n() > 0) %>% ungroup()

  # Only calculate p-value if there are 2 or more groups to compare
  p_value <- if (nlevels(kz$clust) >= 2) {
    tryCatch(summary(aov(value ~ clust, data = kz))[[1]][["Pr(>F)"]][1], error = function(e) NA)
  } else {
    NA
  }

  plot_type <- if (use_violin) 'violin' else 'box'
  p <- plot_ly(kz, x = ~clust, y = ~value, type = plot_type, color = ~clust, colors = stable_color_map)

  p <- p %>% layout(
    title = list(text = feature_name),
    xaxis = list(title = '', showticklabels = FALSE),
    yaxis = list(title = 'Expression', autorange = independent_y),
    showlegend = FALSE,
    annotations = list(x = 0.5, y = 1.05, xref = 'paper', yref = 'paper',
                       text = if (is.na(p_value)) "Single group selected" else paste("ANOVA p-value:", format.pval(p_value, digits = 4)),
                       showarrow = FALSE)
  ) %>% config(displayModeBar = FALSE)

  if (!independent_y) {
    shared_range <- c(0, max(c(data_df$feature1, data_df$feature2), na.rm = TRUE))
    p <- p %>% layout(yaxis = list(title = 'Expression', range = shared_range))
  }
  return(p)
}

#' Add scale bar to spatial plots based on spot spacing.
add_scale_bar <- function(p, seurat_obj) {
  # 1. Get pixel coordinates
  coords <- GetTissueCoordinates(seurat_obj)
  if (is.null(coords) || nrow(coords) < 10) {
    return(p) # Not enough data to calculate
  }

  # 2. Calculate pixel-to-micron ratio (based on 100um between adjacent spots)
  # Use a subset for efficiency
  sample_size <- min(100, nrow(coords))
  dist_matrix <- as.matrix(dist(coords[1:sample_size, c("imagecol", "imagerow")]))
  diag(dist_matrix) <- NA
  # The smallest non-zero distance is our best guess for 100um
  pixels_per_100um <- min(dist_matrix, na.rm = TRUE)

  # 3. Define a fixed length for the scale bar
  bar_length_um <- 1000
  bar_length_px <- (bar_length_um / 100) * pixels_per_100um

  # 4. Define position for the bar (bottom-right corner, OUTSIDE the panel)
  margin_x <- diff(range(coords$imagecol)) * 0.05
  x_end <- max(coords$imagecol) - margin_x
  x_start <- x_end - bar_length_px

  # Y-positioning for OUTSIDE the plot
  y_offset <- diff(range(coords$imagerow)) * 0.05
  y_pos <- min(coords$imagerow) - y_offset

  tick_height <- y_offset * 0.4
  bar_label <- paste(bar_length_um / 1000, "mm")

  # 5. Add the components to the plot
  # --- START of KEY CHANGES ---
  # A. Modify the existing coordinate system instead of adding a new one.
  # This prevents the warning and allows drawing outside the plot panel.
  p$coordinates$clip <- "off"

  p_with_bar <- p +
    # B. Add all scale bar components using annotate()
    annotate("segment", x = x_start, xend = x_end, y = y_pos, yend = y_pos,
             color = "black", linewidth = 1) +
    annotate("segment", x = x_start, xend = x_start, y = y_pos, yend = y_pos + tick_height,
             color = "black", linewidth = 1) +
    annotate("segment", x = x_end, xend = x_end, y = y_pos, yend = y_pos + tick_height,
             color = "black", linewidth = 1) +
    annotate("text", x = x_start + bar_length_px / 2, y = y_pos,
             label = bar_label, vjust = 2.0,
             color = "black", size = 4) +

    # C. Add a margin at the bottom of the plot to make space for the bar
    theme(plot.margin = margin(t = 5, r = 5, b = 30, l = 5, unit = "pt"))

  return(p_with_bar)
}


# ============================================================================ #
#                  SPATIAL RELATIONSHIP PLOTTING FUNCTIONS                     #
# ============================================================================ #

#' Create a Split-View LISA Clustering Plot with a Left-Side Legend and Sorted Bars
plot_lisa_clustering <- function(k, selectedCellTypes = NULL, ann_col = "id", f1_name = "F1", f2_name = "F2") {
  k$id <- k[[ann_col]]

  # Handle missing feature2
  if (all(k$LISA == "not.applicable")) {
    return(
      ggplot(k, aes(imagecol, imagerow)) +
        geom_point(color = "grey80", size = 1.5) +
        theme_void() + scale_y_reverse() + coord_fixed() +
        labs(title = "LISA Clustering", subtitle = "Secondary feature not selected")
    )
  }

  # Create dynamic labels
  k$LISA_dynamic <- k$LISA
  k$LISA_dynamic <- gsub("High\\.High", paste0("High ", f1_name, " | High ", f2_name, " Neighbors"), k$LISA_dynamic)
  k$LISA_dynamic <- gsub("High\\.Low",  paste0("High ", f1_name, " | Low ",  f2_name, " Neighbors"), k$LISA_dynamic)
  k$LISA_dynamic <- gsub("Low\\.High",  paste0("Low ",  f1_name, " | High ", f2_name, " Neighbors"), k$LISA_dynamic)
  k$LISA_dynamic <- gsub("Low\\.Low",   paste0("Low ",  f1_name, " | Low ",  f2_name, " Neighbors"), k$LISA_dynamic)

  dynamic_levels <- c(
    paste0("High ", f1_name, " | High ", f2_name, " Neighbors"),
    paste0("High ", f1_name, " | Low ",  f2_name, " Neighbors"),
    paste0("Low ",  f1_name, " | High ", f2_name, " Neighbors"),
    paste0("Low ",  f1_name, " | Low ",  f2_name, " Neighbors")
  )
  k$LISA_dynamic <- factor(k$LISA_dynamic, levels = dynamic_levels)
  lisa_colors <- c("#d62728", "#ff7f0e", "#9467bd", "#1f77b4")
  names(lisa_colors) <- dynamic_levels

  # Filter by selected cell types
  if (!is.null(selectedCellTypes)) {
    k <- k[k$id %in% selectedCellTypes, ]
  }
  if (nrow(k) == 0) return(ggplot() + theme_void() + labs(title="No data for selected annotations"))

  # PLOT 1: The Spatial Map (with legend for extraction)
  p_spatial_with_legend <- ggplot(k, aes(imagecol, imagerow, color = LISA_dynamic)) +
    geom_point(size = 1.5) +
    scale_color_manual(values = lisa_colors, name = "Spatial Association Pattern", drop = FALSE) +
    theme_void() +
    scale_y_reverse() +
    coord_fixed() +
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 1)) +
    theme(legend.position = "right", legend.box.margin = margin(0, 0, 0, 10))

  # EXTRACT THE LEGEND as a separate grob
  legend_grob <- ggpubr::get_legend(p_spatial_with_legend)
  p_legend <- as_ggplot(legend_grob)

  # RE-RENDER SPATIAL PLOT with no legend
  p_spatial_no_legend <- p_spatial_with_legend + theme(legend.position = "none")

  # PLOT 2: The Quantitative Bar Chart
  summary_data <- k %>%
    dplyr::count(LISA_dynamic, id, .drop = FALSE) %>%
    dplyr::group_by(LISA_dynamic) %>%
    dplyr::mutate(proportion = n / sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(id = reorder_within(id, proportion, LISA_dynamic))

  p_summary <- ggplot(summary_data, aes(x = proportion, y = id, fill = LISA_dynamic)) +
    geom_col() +
    scale_fill_manual(values = lisa_colors, drop = FALSE) +
    scale_y_reordered() +
    scale_x_continuous(labels = scales::percent_format(), expand = c(0, 0.01)) +
    facet_wrap(~LISA_dynamic, ncol = 1, scales = "free_y") +
    theme_light(base_size = 11) +
    theme(
      legend.position = "none",
      strip.text = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.y = element_text(size = 9),
      axis.title.y = element_blank()
    ) +
    labs(x = "Proportion", y = "Annotation")

  # Combine plots using patchwork in the desired [Legend | Spatial | Summary] order
  combined_plot <- p_legend + p_spatial_no_legend + p_summary +
    plot_layout(widths = c(1.5, 3, 2))

  return(combined_plot)
}

#' Reorder factor levels within facet groups
#'
#' Helper for ggplot2 faceted bar charts to sort bars independently per facet.
#'
#' @param x Factor to reorder
#' @param by Numeric values for ordering
#' @param within Grouping variable (facet)
#' @param fun Aggregation function (default: mean)
#' @param sep Separator string for internal encoding
#' @param ... Additional arguments passed to reorder()
#' @return Reordered factor with encoded levels
reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

#' Scale for reordered factors in faceted plots
#'
#' Companion to reorder_within() - strips the internal encoding from labels.
#'
#' @param ... Arguments passed to scale_y_discrete()
#' @param sep Separator used in reorder_within()
#' @return ggplot2 scale object
scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_y_discrete(labels = function(x) gsub(reg, "", x), ...)
}

#' Create color map for ECM domains
#'
#' Generates consistent colors for ECM domain annotations.
#' Basement: purple, Interstitial: green, others: grey.
#'
#' @param levels_vector Character vector of unique ECM domain levels
#' @return Named character vector of hex color codes
create_ecm_color_map <- function(levels_vector) {
  color_spec <- c(
    Basement = "#6a3d9aff",
    Interstitial = "#2b9e2bff",
    not.assigned = "#7f7f7f"
  )

  color_map <- character(length(levels_vector))
  names(color_map) <- levels_vector

  for (level in levels_vector) {
    if (grepl("Basement", level, ignore.case = TRUE)) {
      color_map[level] <- color_spec["Basement"]
    } else if (grepl("Interstitial", level, ignore.case = TRUE)) {
      color_map[level] <- color_spec["Interstitial"]
    } else if (level == "not.assigned") {
      color_map[level] <- color_spec["not.assigned"]
    } else {
      # Fallback color (#CCCCCC - light grey) for any unrecognized annotations
      color_map[level] <- "#CCCCCC"
    }
  }

  return(color_map)
}

#' Find Matrisome-Pair Enrichment
#'
#' This function calculates enrichment statistics (avg_log2FC and percentage
#' difference) for features in a given activity assay, comparing each group
#' (cluster) to all other groups.
#'
#' @param seurat_obj A Seurat object containing the matrisome-pair assay.
#' @param assay_name A string specifying the name of the assay to analyze
#'        (for example, the "MATRISOMEPAIR" assay).
#' @param group_by A string specifying the metadata column to use for grouping
#'        cells (e.g., "seurat_clusters", "ecm_domain_annotation").
#'
#' @return A tidy dataframe with enrichment statistics for each feature in each group.
#'         Columns: feature, cluster, avg_log2FC, perc_difference, pct.1, pct.2
find_matrisome_pair_enrichment <- function(seurat_obj, assay_name, group_by) {

  cat(sprintf("Calculating enrichment for assay '%s' grouped by '%s'...\n", assay_name, group_by))

  score_matrix <- GetAssayData(seurat_obj, assay = assay_name, layer = "data")
  groups <- seurat_obj[[group_by, drop = TRUE]]
  cell_groups <- split(colnames(score_matrix), groups)

  all_cluster_stats <- lapply(names(cell_groups), function(cluster_name) {
    # Split cells into in-group (current cluster) vs out-of-group (all other clusters)
    group_1_cells <- cell_groups[[cluster_name]]  # In-group
    group_2_cells <- unlist(cell_groups[names(cell_groups) != cluster_name], use.names = FALSE)  # Out-of-group

    # Small constant to prevent log(0) errors in fold change calculation
    # 1e-9 is negligible compared to typical expression values but prevents division by zero
    epsilon <- 1e-9
    group_1_avg <- rowMeans(score_matrix[, group_1_cells, drop = FALSE])
    group_2_avg <- rowMeans(score_matrix[, group_2_cells, drop = FALSE])
    avg_log2FC <- log2((group_1_avg + epsilon) / (group_2_avg + epsilon))

    pct_1 <- rowMeans(score_matrix[, group_1_cells, drop = FALSE] > 0)
    pct_2 <- rowMeans(score_matrix[, group_2_cells, drop = FALSE] > 0)
    perc_difference <- pct_1 - pct_2

    data.frame(
      feature = rownames(score_matrix),
      cluster = cluster_name,
      avg_log2FC = avg_log2FC,
      perc_difference = perc_difference,
      pct.1 = pct_1, # Keeping this for potential filtering
      pct.2 = pct_2
    )
  })

  return(bind_rows(all_cluster_stats))
}


#' Coalesce Reciprocal Rows into Unique Matrisome Pairs
#'
#' Identifies reciprocal rows (e.g., A-B and B-A) in statistics and mean-score
#' matrices and combines them into one unordered heterotypic pair (e.g., A-B).
#' Non-reciprocal pairs are retained under the same canonical naming scheme.
#'
#' @param stats_df A long-format dataframe of statistics from
#'   `find_matrisome_pair_enrichment`.
#'        Must contain 'feature' and 'cluster' columns.
#' @param means_matrix A wide-format matrix of mean scores, with features as
#'        rows and clusters as columns.
#'
#' @return A list containing the coalesced `stats` dataframe and `means` matrix.
#'
coalesce_reciprocal_pairs <- function(stats_df, means_matrix) {

  cat("Coalescing reciprocal rows into unique matrisome pairs...\n")

  # Early return if inputs are empty
  if (nrow(stats_df) == 0 || nrow(means_matrix) == 0) {
    warning("Empty stats_df or means_matrix - returning empty results")
    return(list(
      stats = data.frame(feature = character(), cluster = character(),
                         avg_log2FC = numeric(), perc_difference = numeric()),
      means = matrix(nrow = 0, ncol = ncol(means_matrix),
                     dimnames = list(NULL, colnames(means_matrix)))
    ))
  }

  # STEP 1: Standardize cluster names
  # Remove spaces and hyphens to ensure consistent matching across datasets
  # (e.g., "Interstitial ECM" -> "Interstitial_ECM")
  stats_df$cluster <- gsub(" |\\-", "_", stats_df$cluster)
  colnames(means_matrix) <- gsub(" |\\-", "_", colnames(means_matrix))

  # STEP 2: Create canonical feature names
  # Sort pair components alphabetically to identify reciprocal rows
  # Example: "TGFB1-TGFBR1" and "TGFBR1-TGFB1" both become "TGFB1-TGFBR1"
  original_features <- rownames(means_matrix)
  canonical_names <- sapply(original_features, function(name) {
    paste(sort(strsplit(name, "-")[[1]]), collapse = "-")
  })

  # STEP 3: Use canonical names for every unordered pair. MatriSpace does not
  # infer directionality or label these co-expression scores as communication.
  new_feature_names <- canonical_names

  # Create mapping dataframe for later joining
  name_map <- data.frame(
    original_feature = original_features,
    new_feature = new_feature_names
  )

  # STEP 5: Coalesce using rowsum
  # Combines reciprocal pairs by summing their scores (averaging happens implicitly
  # when divided by 2 in subsequent analysis)
  coalesced_means_matrix <- rowsum(as.matrix(means_matrix), group = new_feature_names)

  # STEP 6: Coalesce statistics dataframe
  # Join with mapping and average enrichment statistics for reciprocal pairs
  coalesced_stats <- stats_df %>%
    left_join(name_map, by = c("feature" = "original_feature"))

  # Handle case where join produced no matches (new_feature all NA)
  if (!"new_feature" %in% colnames(coalesced_stats)) {
    warning("'new_feature' column not found after join. Using original feature names.")
    coalesced_stats$new_feature <- coalesced_stats$feature
  } else if (all(is.na(coalesced_stats$new_feature))) {
    warning("'new_feature' column is all NA. Using original feature names.")
    coalesced_stats$new_feature <- coalesced_stats$feature
  }

  coalesced_stats <- coalesced_stats %>%
    group_by(.data$new_feature, .data$cluster) %>%
    summarise(
      avg_log2FC = mean(avg_log2FC, na.rm = TRUE),
      perc_difference = mean(perc_difference, na.rm = TRUE),
      .groups = 'drop'
    )

  # Rename new_feature to feature (use base R to avoid dplyr scoping issues)
  colnames(coalesced_stats)[colnames(coalesced_stats) == "new_feature"] <- "feature"

  return(list(
    stats = coalesced_stats,
    means = coalesced_means_matrix
  ))
}

#' Calculate Spatial Matrisome-Pair Co-expression Scores
#'
#' This function calculates pair scores on a single core, designed to be
#' called from within a Shiny observer with progress feedback (via shinybusy).
#' It avoids all parallel backend complexity to ensure stability in any environment.
#
#' @param seurat_obj A Seurat object with expression data.
#' @param pair_db A MatriComDB pair dataframe. The `Gene1` and `Gene2` columns
#'   identify the two members of each pair.
#' @param adj_matrix A sparse adjacency matrix defining spatial neighbors.
#' @param assay The assay to use for expression data (default: "SCT").
#' @param layer The layer to use within the assay (default: "data").
#' @param update_progress Optional progress update function (e.g., from shinybusy).
#'        Called at key steps with signature: update_progress(value, text)
#' @param chunk_size Number of pairs to process per batch (default: 5000).
#'        Lower values reduce peak memory usage but may be slightly slower.
#'
#' @return A Seurat object with the new "MATRISOMEPAIR" assay added.
#'
compute_spatial_matrisome_pair_scores_single_core <- function(seurat_obj, pair_db, adj_matrix, assay = "SCT", layer = "data", update_progress = NULL, chunk_size = 5000) {

  # 1. Prepare data
  expr <- LayerData(seurat_obj, assay = assay, layer = layer)

  # Ensure adjacency matrix matches expression matrix spot order
  if (!all(colnames(expr) == colnames(adj_matrix))) {
    adj_matrix <- adj_matrix[colnames(expr), colnames(expr)]
  }

  # Keep pairs whose two genes are present in this dataset
  genes_in_data <- rownames(expr)
  pair_db_filtered <- pair_db %>%
    filter(Gene1 %in% genes_in_data & Gene2 %in% genes_in_data)

  if (nrow(pair_db_filtered) == 0) stop("No valid matrisome pairs were found in the expression data.")
  n_pairs <- nrow(pair_db_filtered)
  cat("Found", n_pairs, "valid matrisome pairs to score.\n")

  # Progress update: starting calculation
  if (!is.null(update_progress)) {
    update_progress(value = 0.2, text = "Preparing expression data...")
  }

  # 2. Prepare expression matrix with only genes represented in the pairs
  pair_genes <- unique(c(pair_db_filtered$Gene1, pair_db_filtered$Gene2))
  expr_pairs <- expr[pair_genes, , drop = FALSE]

  # Ensure sparse matrix format (dgCMatrix) for efficient operations
  if (!inherits(expr_pairs, "dgCMatrix")) {
    expr_pairs <- as(expr_pairs, "dgCMatrix")
  }

  # Apply sqrt transformation directly on sparse matrix's data slot
  expr_pairs@x <- sqrt(pmax(expr_pairs@x, 0))

  # Sparsify adjacency matrix for efficient multiplication
  adj_matrix <- as(adj_matrix, "dgCMatrix")

  # 3. Compute neighbour-weighted expression for all genes at once
  if (!is.null(update_progress)) {
    update_progress(value = 0.3, text = "Computing neighbour expression...")
  }
  neighbor_expr <- t(adj_matrix %*% t(expr_pairs))

  # Clean up original expression matrix to free memory
  rm(expr)
  gc(verbose = FALSE)

  # 4. Process interactions in chunks to control peak memory usage
  #    This prevents OOM on memory-constrained hosted deployments
  n_chunks <- ceiling(n_pairs / chunk_size)
  coexpression_list <- vector("list", n_chunks)

  for (i in seq_len(n_chunks)) {
    # Calculate chunk boundaries
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, n_pairs)
    chunk_rows <- start_idx:end_idx

    # Update progress
    if (!is.null(update_progress)) {
      progress_val <- 0.3 + (0.6 * i / n_chunks)
      update_progress(value = progress_val, text = sprintf("Processing matrisome pairs (%d/%d)...", i, n_chunks))
    }

    # Get indices for this chunk
    gene_1_idx <- match(pair_db_filtered$Gene1[chunk_rows], rownames(expr_pairs))
    gene_2_idx <- match(pair_db_filtered$Gene2[chunk_rows], rownames(expr_pairs))

    # Compute co-expression scores for this chunk
    chunk_coexpression <- neighbor_expr[gene_1_idx, , drop = FALSE] * expr_pairs[gene_2_idx, , drop = FALSE]
    rownames(chunk_coexpression) <- paste(
      pair_db_filtered$Gene1[chunk_rows],
      pair_db_filtered$Gene2[chunk_rows],
      sep = "-"
    )

    coexpression_list[[i]] <- chunk_coexpression

    # Force garbage collection between chunks
    if (i < n_chunks) gc(verbose = FALSE)
  }

  # 5. Combine chunks into final activity matrix
  if (!is.null(update_progress)) {
    update_progress(value = 0.95, text = "Finalising results...")
  }

  coexpression <- do.call(rbind, coexpression_list)
  rm(coexpression_list, neighbor_expr, expr_pairs)
  gc(verbose = FALSE)

  # 6. Store results in Seurat object
  seurat_obj[["MATRISOMEPAIR"]] <- CreateAssayObject(counts = coexpression)
  seurat_obj <- SetAssayData(
    seurat_obj,
    assay = "MATRISOMEPAIR",
    layer = "data",
    new.data = coexpression
  )

  return(seurat_obj)
}

# --------------------------------------------------------------------------- #
# SpatialExperiment to Seurat Conversion
# Adapted from VisiumStitched package
# --------------------------------------------------------------------------- #

#' Convert a SpatialExperiment object to a Seurat object
#'
#' @param spe A SpatialExperiment object (single or multi-sample)
#' @param spatial_cols Named character vector mapping spatial columns.
#'   Defaults work for standard Visium. Set elements to NA to use fallback values.
#' @param symbol_col Column name in rowData containing gene symbols. Auto-detected if NULL.
#' @param verbose Logical, whether to print progress messages
#' @return A Seurat object with spatial image
#'
spe_to_seurat <- function(
    spe,
    spatial_cols = c(
      "tissue" = "in_tissue",
      "row" = "array_row",
      "col" = "array_col",
      "imagerow" = "pxl_row_in_fullres",
      "imagecol" = "pxl_col_in_fullres"
    ),
    symbol_col = NULL,
    verbose = TRUE
) {

  # ============ SANITY CHECKS ============

  # 1. Check class
  if (!inherits(spe, "SpatialExperiment")) {
    stop("Input must be a SpatialExperiment object")
  }

  # 2. Multi-sample handling: auto-subset first sample
  sample_ids <- unique(spe$sample_id)
  if (length(sample_ids) > 1) {
    if (verbose) message(sprintf(
      "Multi-sample SPE detected (%d samples). Using first: '%s'",
      length(sample_ids), sample_ids[1]
    ))
    spe <- spe[, spe$sample_id == sample_ids[1]]
  }
  sample_id <- unique(spe$sample_id)

  # 3. Check counts assay exists
  if (!"counts" %in% SummarizedExperiment::assayNames(spe)) {
    stop("No 'counts' assay found. Available: ", paste(SummarizedExperiment::assayNames(spe), collapse = ", "))
  }

  # 4. Check required spatial_cols elements
  required_elements <- c("imagerow", "imagecol")
  missing_elements <- setdiff(required_elements, names(spatial_cols))
  if (length(missing_elements) > 0) {
    stop("spatial_cols must contain: ", paste(missing_elements, collapse = ", "))
  }

  # 5. Combine colData and spatialCoords for column lookup
  col_info <- cbind(SummarizedExperiment::colData(spe), SpatialExperiment::spatialCoords(spe))

  # 6. Check imagerow/imagecol exist (required)
  for (coord in c("imagerow", "imagecol")) {
    col_name <- spatial_cols[coord]
    if (is.na(col_name) || !col_name %in% colnames(col_info)) {
      stop(sprintf(
        "'%s' column '%s' not found. Available: %s",
        coord, col_name, paste(colnames(col_info), collapse = ", ")
      ))
    }
  }

  # 7. Check imgData exists with lowres image
  img_data <- SpatialExperiment::imgData(spe)
  if (nrow(img_data) == 0) {
    stop("No image data found in imgData(spe)")
  }
  if (!"lowres" %in% img_data$image_id) {
    stop("No 'lowres' image found. Available: ", paste(img_data$image_id, collapse = ", "))
  }

  # 8. Check for duplicate barcodes
  if (any(duplicated(colnames(spe)))) {
    warning("Duplicate cell barcodes found - may cause issues")
  }

  # 9. Auto-detect or validate symbol column
  rd_cols <- colnames(SummarizedExperiment::rowData(spe))
  common_symbol_cols <- c("symbol", "gene_name", "Symbol", "gene_short_name", "SYMBOL")

  if (is.null(symbol_col)) {
    detected <- intersect(common_symbol_cols, rd_cols)
    if (length(detected) > 0) {
      symbol_col <- detected[1]
      if (verbose) message(sprintf("Auto-detected symbol column: '%s'", symbol_col))
    }
  } else if (!symbol_col %in% rd_cols) {
    warning(sprintf("Specified symbol_col '%s' not found in rowData. Available: %s",
                    symbol_col, paste(rd_cols, collapse = ", ")))
    symbol_col <- NULL
  }

  has_symbols <- !is.null(symbol_col) && !all(is.na(SummarizedExperiment::rowData(spe)[[symbol_col]]))
  if (!has_symbols && verbose) {
    message("No symbol column found - using rownames as gene names")
    message(sprintf("rowData columns: %s", paste(rd_cols, collapse = ", ")))
  }

  # 10. Check for altExps
  if (length(SingleCellExperiment::altExpNames(spe)) > 0 && verbose) {
    message(sprintf("Note: Ignoring %d altExp(s): %s",
                    length(SingleCellExperiment::altExpNames(spe)),
                    paste(SingleCellExperiment::altExpNames(spe), collapse = ", ")))
  }

  # 11. QC stats
  if (verbose) {
    n_genes <- nrow(spe)
    n_spots <- ncol(spe)
    n_zero_genes <- sum(Matrix::rowSums(SummarizedExperiment::assay(spe, "counts")) == 0)
    n_zero_spots <- sum(Matrix::colSums(SummarizedExperiment::assay(spe, "counts")) == 0)
    message(sprintf("Input: %d genes, %d spots", n_genes, n_spots))
    if (n_zero_genes > 0) message(sprintf("Note: %d genes with zero counts", n_zero_genes))
    if (n_zero_spots > 0) warning(sprintf("%d spots with zero counts", n_zero_spots))
  }

  # ============ CONVERSION ============

  SPOT_DIAMETER <- 55e-6

  if (verbose) message("Converting to Seurat object...")

  # Remove altExps to prevent Seurat from using them
  spe_clean <- spe
  if (length(SingleCellExperiment::altExpNames(spe_clean)) > 0) {
    for (ae in SingleCellExperiment::altExpNames(spe_clean)) {
      SingleCellExperiment::altExp(spe_clean, ae) <- NULL
    }
  }

  seur <- Seurat::as.Seurat(spe_clean, counts = "counts", data = NULL)

  # Get assay name dynamically
  assay_name <- Seurat::DefaultAssay(seur)
  if (verbose) message(sprintf("Active assay: %s", assay_name))

  # Update gene names to symbols if available
  if (has_symbols) {
    if (verbose) message(sprintf("Mapping gene symbols from '%s' column...", symbol_col))
    ensembl_ids <- rownames(seur)
    symbols <- SummarizedExperiment::rowData(spe)[[symbol_col]]
    new_names <- ifelse(is.na(symbols) | symbols == "", ensembl_ids, symbols)
    new_names <- make.unique(new_names)

    rownames(seur@assays[[assay_name]]@counts) <- new_names
    rownames(seur@assays[[assay_name]]@data) <- new_names
  }

  if (verbose) message(sprintf("Adding coordinates and image for sample %s...", sample_id))

  # Helper to get column with fallback
  get_col <- function(key, fallback) {
    col_name <- spatial_cols[key]
    if (!is.na(col_name) && col_name %in% colnames(col_info)) {
      return(col_info[[col_name]])
    }
    return(fallback)
  }

  # Build coordinates - use actual values when available, fallback otherwise
  coords <- data.frame(
    tissue = as.integer(get_col("tissue", rep(1L, ncol(spe)))),
    row = get_col("row", seq_len(ncol(spe))),
    col = get_col("col", seq_len(ncol(spe))),
    imagerow = col_info[[spatial_cols["imagerow"]]],
    imagecol = col_info[[spatial_cols["imagecol"]]],
    row.names = colnames(spe)
  )

  # Convert image to array
  this_img <- array(
    t(grDevices::col2rgb(SpatialExperiment::imgRaster(spe))),
    dim = c(dim(SpatialExperiment::imgRaster(spe)), 3)
  ) / 256

  # Get scale factor
  sf <- img_data$scaleFactor[img_data$image_id == "lowres"]

  # Sanitize sample_id for use as key (alphanumeric only)
  safe_key <- paste0(gsub("[^[:alnum:]]", "", sample_id), "_")

  # Create VisiumV1 object
  seur@images[[sample_id]] <- new(
    Class = "VisiumV1",
    image = this_img,
    scale.factors = Seurat::scalefactors(
      spot = NA,
      fiducial = NA,
      hires = NA,
      lowres = sf
    ),
    coordinates = coords,
    spot.radius = SPOT_DIAMETER / sf,
    assay = assay_name,
    key = safe_key
  )

  if (verbose) message("Done!")
  return(seur)
}

#' k-nearest-neighbour search over spot coordinates.
#'
#' Uses RANN::nn2 (a Seurat dependency, so it should be present in every deployment unless they screw up badly). Falls
#' back to an exact chunked search for small objects if RANN is unavailable,
#' and refuses the O(n^2) fallback on large ones rather than stalling a worker forever.
#'
#' @param coords Numeric matrix, n x 2.
#' @param k Number of neighbours to return (including self).
#' @return list(idx, dist), each n x k.
.matrispace_knn <- function(coords, k) {
  n <- nrow(coords)
  k <- as.integer(min(k, n))

  if (requireNamespace("RANN", quietly = TRUE)) {
    nn <- RANN::nn2(data = coords, query = coords, k = k)
    return(list(idx = nn$nn.idx, dist = nn$nn.dists))
  }

  if (n > 20000L) {
    stop("Leakage correction needs the RANN package for objects with more than ",
         "20,000 spots. Install RANN (it ships with Seurat) and reload.")
  }

  idx  <- matrix(0L,  n, k)
  dst  <- matrix(Inf, n, k)
  block <- max(1L, floor(2e7 / n))
  for (start in seq(1L, n, by = block)) {
    rows <- start:min(n, start + block - 1L)
    d2 <- outer(coords[rows, 1L], coords[, 1L], "-")^2 +
          outer(coords[rows, 2L], coords[, 2L], "-")^2
    for (ii in seq_along(rows)) {
      o <- order(d2[ii, ])[seq_len(k)]
      idx[rows[ii], ] <- o
      dst[rows[ii], ] <- sqrt(d2[ii, o])
    }
  }
  list(idx = idx, dist = dst)
}

#' Build a row-stochastic neighbour-weight matrix over spots.
#'
#' Row i holds 1/deg(i) at the columns of i's retained neighbours and zero on
#' the diagonal, so W %*% x is the local neighbour mean of x. Spots with no
#' retained neighbour get an all-zero row and are therefore left uncorrected.
#'
#' @param coords Data frame or matrix of spot coordinates; the first two
#'   columns are used.
#' @param k Neighbours per spot (6 matches the Visium hexagonal ring).
#' @param dist_factor Multiple of the median nearest-neighbour distance beyond
#'   which a neighbour is discarded. Inf disables the cap.
#' @return A sparse n x n dgCMatrix, or NULL if a weight matrix cannot be built.
matrispace_leakage_weights <- function(coords, k = 6L, dist_factor = 1.6) {
  coords <- as.matrix(coords[, 1:2, drop = FALSE])
  storage.mode(coords) <- "double"
  if (anyNA(coords)) return(NULL)

  n <- nrow(coords)
  if (n < 8L) return(NULL)
  k <- as.integer(min(k, n - 1L))
  if (k < 1L) return(NULL)

  nn   <- .matrispace_knn(coords, k = k + 1L)
  idx  <- nn$idx[,  -1L, drop = FALSE]   # column 1 is the spot itself
  dst  <- nn$dist[, -1L, drop = FALSE]

  d1 <- dst[, 1L]
  d1 <- d1[is.finite(d1) & d1 > 0]
  cutoff <- if (length(d1) > 0L && is.finite(dist_factor)) {
    stats::median(d1) * dist_factor
  } else Inf

  keep <- is.finite(dst) & dst <= cutoff & idx >= 1L & idx <= n
  if (!any(keep)) return(NULL)

  i <- row(idx)[keep]
  j <- as.integer(idx[keep])
  deg <- tabulate(i, nbins = n)
  deg[deg == 0L] <- 1L

  Matrix::sparseMatrix(i = i, j = j, x = 1 / deg[i], dims = c(n, n))
}

#' Apply the first-order leakage correction to one sparse expression matrix.
#'
#' @param mat Sparse genes x spots matrix.
#' @param Wt Transposed spot weight matrix from matrispace_leakage_weights().
#' @param alpha Leakage coefficient; 0 disables the correction.
#' @param log_scale TRUE for a log1p-scaled layer ("data"), FALSE for counts.
#' @param max_block_nnz Approximate ceiling on the non-zeros held in one block.
#' @return Corrected sparse matrix with the original dimnames.
matrispace_correct_matrix <- function(mat, Wt, alpha,
                                      log_scale = FALSE,
                                      max_block_nnz = 3e7) {
  if (is.null(mat) || !is.finite(alpha) || alpha <= 0) return(mat)
  if (!inherits(mat, "sparseMatrix")) return(mat)

  mat <- methods::as(mat, "CsparseMatrix")
  ng  <- nrow(mat)
  nnz <- length(mat@x)
  if (ng == 0L || nnz == 0L) return(mat)

  rows_per_block <- max(1L, min(ng, as.integer(floor(ng * max_block_nnz /
                                                     max(nnz * 7, 1)))))
  starts <- seq(1L, ng, by = rows_per_block)
  out <- vector("list", length(starts))

  for (b in seq_along(starts)) {
    rows <- starts[b]:min(ng, starts[b] + rows_per_block - 1L)
    sub  <- methods::as(mat[rows, , drop = FALSE], "CsparseMatrix")
    if (log_scale) sub@x <- expm1(sub@x)

    res <- methods::as(sub - alpha * (sub %*% Wt), "CsparseMatrix")
    res@x[res@x < 0] <- 0
    res <- Matrix::drop0(res)
    if (log_scale) res@x <- log1p(res@x)

    out[[b]] <- res
  }

  res <- if (length(out) == 1L) out[[1L]] else do.call(rbind, out)
  dimnames(res) <- dimnames(mat)
  res
}

#' Splice corrected rows back into a sparse matrix without changing its pattern.
#'
#' Rewriting the @x slot in place is exact and costs one vectorised pass, instead of the sparse re-assembly that `mat[rows, ] <- ...` would trigger.
#'
#' @param mat Original sparse genes x spots matrix (dgCMatrix).
#' @param ridx Integer row indices of `mat` that were corrected.
#' @param csub Corrected submatrix, length(ridx) x ncol(mat), rows in the order
#'   of `ridx`.
#' @param max_dense_cells Ceiling on the block held at once.
#' @return A new dgCMatrix; `mat` is left untouched.
.matrispace_splice_rows <- function(mat, ridx, csub, max_dense_cells = 2e7) {
  if (!inherits(mat, "dgCMatrix")) {
    mat[ridx, ] <- csub
    return(mat)
  }
  nnz <- length(mat@x)
  if (nnz == 0L || length(ridx) == 0L) return(mat)

  pos <- integer(nrow(mat))
  pos[ridx] <- seq_along(ridx)

  ent_row <- pos[mat@i + 1L]
  ent <- which(ent_row > 0L)
  if (length(ent) == 0L) return(mat)

  ent_row <- ent_row[ent]
  ent_col <- rep.int(seq_len(ncol(mat)), diff(mat@p))[ent]

  newx <- mat@x
  ncell <- as.numeric(length(ridx)) * ncol(mat)
  if (ncell <= max_dense_cells) {
    cd <- as.matrix(csub)
    newx[ent] <- cd[cbind(ent_row, ent_col)]
  } else {
    block <- max(1L, as.integer(floor(max_dense_cells / max(length(ridx), 1L))))
    for (start in seq(1L, ncol(mat), by = block)) {
      stop_c <- min(ncol(mat), start + block - 1L)
      sel <- ent_col >= start & ent_col <= stop_c
      if (!any(sel)) next
      cd <- as.matrix(csub[, start:stop_c, drop = FALSE])
      newx[ent[sel]] <- cd[cbind(ent_row[sel], ent_col[sel] - start + 1L)]
    }
  }

  result <- methods::new("dgCMatrix", i = mat@i, p = mat@p, x = newx,
                         Dim = mat@Dim, Dimnames = mat@Dimnames)
  Matrix::drop0(result)
}

#' SpotClean-style leakage correction of a spatial Seurat object.
#'
#' Pass `genes = NULL` to correct every row (much slower, not endorsed).
#'
#' Metadata columns (pre-computed UCell signature scores, ScType ECM domain
#' annotation) are never recomputed here, but they carry over leakage across which is supposed to be stochastic
#' and should therefore cancel out in normalization vs housekeeping genes.
#'
#' @param obj Seurat object with spatial coordinates.
#' @param genes Character vector of features to correct, or NULL for all.
#' @param alpha Leakage coefficient. 0 returns the object untouched.
#' @param k Neighbors per spot.
#' @param dist_factor Distance cap, as a multiple of the median NN distance.
#' @param assays Assays to correct; NULL means every assay in the object.
#' @param verbose progress messages.
#' @return list(object, applied, n_spots, n_corrected, n_genes, layers, message)
matrispace_spotclean <- function(obj, genes = NULL, alpha = 0.10, k = 6L,
                                 dist_factor = 1.6, assays = NULL,
                                 verbose = FALSE) {
  fail <- function(msg) list(object = obj, applied = FALSE, n_spots = NA_integer_,
                             n_corrected = 0L, n_genes = 0L,
                             layers = character(0), message = msg)

  if (!inherits(obj, "Seurat")) return(fail("Not a Seurat object."))
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha > 1) {
    return(fail("Leakage coefficient must be greater than zero and at most one."))
  }
  if (length(k) != 1L || !is.finite(k) || k < 1) {
    return(fail("Neighbour count must be a positive integer."))
  }
  k <- as.integer(k)
  if (length(dist_factor) != 1L || is.na(dist_factor) || dist_factor <= 0) {
    return(fail("Distance factor must be positive."))
  }

  scoped <- !is.null(genes)
  if (scoped) {
    genes <- unique(as.character(genes))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (length(genes) == 0L) return(fail("No features selected for correction."))
  }

  cells <- colnames(obj)
  n <- length(cells)

  coords <- tryCatch(align_spatial_coordinates(obj, cells = cells),
                     error = function(e) NULL)
  if (is.null(coords) || nrow(coords) < 8L) {
    return(fail("Spatial coordinates unavailable; leakage correction skipped."))
  }
  coords <- as.data.frame(coords)

  xy <- if (all(c("imagecol", "imagerow") %in% colnames(coords))) {
    coords[, c("imagecol", "imagerow"), drop = FALSE]
  } else if (all(c("x", "y") %in% colnames(coords))) {
    coords[, c("x", "y"), drop = FALSE]
  } else if (all(c("col", "row") %in% colnames(coords))) {
    coords[, c("col", "row"), drop = FALSE]
  } else NULL
  if (is.null(xy)) return(fail("No usable coordinate columns; correction skipped."))

  xy <- data.frame(x = suppressWarnings(as.numeric(xy[[1L]])),
                   y = suppressWarnings(as.numeric(xy[[2L]])))
  ok <- is.finite(xy$x) & is.finite(xy$y)
  common <- intersect(cells, rownames(coords)[ok])
  if (length(common) < 8L) return(fail("Too few spots with valid coordinates."))

  W_sub <- matrispace_leakage_weights(xy[match(common, rownames(coords)), , drop = FALSE],
                                      k = k, dist_factor = dist_factor)
  if (is.null(W_sub)) return(fail("Could not build a spot neighbourhood graph."))

  pos <- match(common, cells)
  Ws  <- methods::as(W_sub, "TsparseMatrix")
  W   <- Matrix::sparseMatrix(i = pos[Ws@i + 1L], j = pos[Ws@j + 1L], x = Ws@x,
                              dims = c(n, n))
  Wt  <- Matrix::t(W)
  rm(W_sub, Ws, W)

  target <- if (is.null(assays)) SeuratObject::Assays(obj) else
    intersect(assays, SeuratObject::Assays(obj))
  touched <- character(0)
  hit_genes <- character(0)

  for (a in target) {
    for (lyr in c("counts", "data")) {
      m <- tryCatch(SeuratObject::LayerData(obj, assay = a, layer = lyr),
                    error = function(e) NULL)
      if (is.null(m) || length(dim(m)) != 2L) next
      if (nrow(m) == 0L || ncol(m) != n) next
      if (!inherits(m, "sparseMatrix")) next

      # The spatial weights are ordered against colnames(obj). Never apply
      # them positionally to a layer whose spots cannot be aligned exactly.
      layer_cells <- colnames(m)
      if (is.null(layer_cells) || !setequal(layer_cells, cells)) next
      if (!identical(layer_cells, cells)) {
        m <- m[, cells, drop = FALSE]
      }

      log_scale <- identical(lyr, "data")

      ridx <- if (scoped) {
        r <- match(genes, rownames(m))
        r[!is.na(r)]
      } else integer(0)
      if (scoped && length(ridx) == 0L) next
      if (scoped) hit_genes <- union(hit_genes, rownames(m)[ridx])

      corrected <- tryCatch({
        if (!scoped) {
          matrispace_correct_matrix(m, Wt, alpha, log_scale = log_scale)
        } else {
          csub <- matrispace_correct_matrix(m[ridx, , drop = FALSE], Wt, alpha,
                                            log_scale = log_scale)
          .matrispace_splice_rows(m, ridx, csub)
        }
      }, error = function(e) { warning(sprintf(
        "Leakage correction failed on %s/%s (%s); layer left unchanged.",
        a, lyr, conditionMessage(e))); NULL })
      if (is.null(corrected)) next

      ok_set <- tryCatch({
        SeuratObject::LayerData(obj, assay = a, layer = lyr) <- corrected
        TRUE
      }, error = function(e) { warning(sprintf(
        "Could not write corrected %s/%s (%s); layer left unchanged.",
        a, lyr, conditionMessage(e))); FALSE })

      if (isTRUE(ok_set)) {
        touched <- c(touched, paste(a, lyr, sep = "/"))
        if (verbose) message(sprintf("Leakage correction applied to %s/%s.", a, lyr))
      }
      rm(m, corrected)
    }
  }

  if (length(touched) == 0L) {
    return(fail(if (scoped)
      "None of the selected features were found in a correctable layer."
    else "No sparse counts/data layer could be corrected."))
  }

  list(object      = obj,
       applied     = TRUE,
       n_spots     = n,
       n_corrected = length(common),
       n_genes     = if (scoped) length(hit_genes) else nrow(obj),
       layers      = touched,
       message     = sprintf("%s feature%s, \u03b1 = %.2f, k = %d, %d/%d spots in-graph",
                             if (scoped) format(length(hit_genes), big.mark = ",") else "all",
                             if (!scoped || length(hit_genes) != 1L) "s" else "",
                             alpha, k, length(common), n))
}
