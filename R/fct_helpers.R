# helpers.R
# This file contains all helper functions used in the server logic.

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
standardize_gene_symbols <- function(seurat_obj) {
  # Collect all unique features across every assay
  all_features <- unique(unlist(lapply(SeuratObject::Assays(seurat_obj), function(a) {
    rownames(seurat_obj[[a]])
  })))

  # Get HGNC mapping (returns a data frame)
  map_df <- suppressWarnings(
    scCustomize::Updated_HGNC_Symbols(all_features, verbose = FALSE, case_check_as_warn = TRUE)
  )

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
  spot_coords <- GetTissueCoordinates(obj)
  metadata <- obj@meta.data
  data_plot <- cbind(spot_coords, metadata)
  data_plot
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

#' Add feature analysis results and LISA statistics to Seurat object.
addfeat <- function(obj, feat1, sel1, feat2, sel2) {
  m <- obj@meta.data

  # Handle primary feature (feature1)
  if (feat1 == "") {
    m$feature1 <- 0
  } else if (sel1 == "matrisome gene" || sel1 == "any gene") {
    feature1_data <- safe_get_assay_data(obj)
    if (!is.null(feature1_data) && feat1 %in% rownames(feature1_data)) {
      m$feature1 <- feature1_data[feat1, ]
    } else {
      warning("Could not retrieve data for feature1")
      m$feature1 <- 0
    }
  } else if (sel1 == "matrisome signature") {
    if (feat1 %in% colnames(obj@meta.data)) {
      m$feature1 <- obj@meta.data[, feat1]
    } else {
      m$feature1 <- 0
    }
  }

  # Handle secondary feature (feature2)
  if (is.null(sel2) || feat2 == "" || sel2 == "none") {
    m$feature2 <- 0
  } else if (sel2 == "matrisome gene" || sel2 == "any gene") {
    feature2_data <- safe_get_assay_data(obj)
    if (!is.null(feature2_data) && feat2 %in% rownames(feature2_data)) {
      m$feature2 <- feature2_data[feat2, ]
    } else {
      warning("Could not retrieve data for feature2")
      m$feature2 <- 0
    }
  } else if (sel2 == "matrisome signature") {
    if (feat2 %in% colnames(obj@meta.data)) {
      m$feature2 <- obj@meta.data[, feat2]
    } else {
      m$feature2 <- 0
    }
  }

  # Get spatial coordinates using GetTissueCoordinates for consistency
  coords <- GetTissueCoordinates(obj)
  if (!all(c("row", "col") %in% names(coords))) {
    # Handle VisiumV1 vs VisiumV2 coordinate structures
    if ("coordinates" %in% slotNames(obj@images[[1]])) {
      raw_coords <- obj@images[[1]]@coordinates
      coords$row <- raw_coords$row
      coords$col <- raw_coords$col
    } else {
      # VisiumV2: use x/y as fallback
      coords$row <- coords$y
      coords$col <- coords$x
    }
  }
  # Ensure imagerow/imagecol always exist (needed for downstream plotting/export)
  if (!all(c("imagerow", "imagecol") %in% names(coords))) {
    coords$imagerow <- coords$x
    coords$imagecol <- coords$y
  }
  m <- cbind(m, coords)

  # Clean up negative values
  m$feature1[m$feature1 < 0] <- 0
  m$feature2[m$feature2 < 0] <- 0

  # MODIFIED & IMPROVED LISA CALCULATION
  # Only run if a secondary feature is present and has variance
  if (!is.null(sel2) && sel2 != "none" && feat2 != "" && var(m$feature2, na.rm = TRUE) > 0) {
    # Calculate LISA statistics
    Wij <- as.matrix(dist(as.matrix(cbind(m$row, m$col))))
    Wij[Wij == 0] <- 1e-9
    Wij <- 1 / Wij
    diag(Wij) <- 0
    if (sum(Wij) > 0) Wij <- Wij / sum(Wij)
    diag(Wij) <- 0

    # Standardize features (Z-score)
    x <- m$feature1
    y <- m$feature2

    x_scaled <- if(sd(x, na.rm = TRUE) > 0) as.numeric(scale(x)) else rep(0, length(x))
    y_scaled <- if(sd(y, na.rm = TRUE) > 0) as.numeric(scale(y)) else rep(0, length(y))
    x_scaled[is.na(x_scaled)] <- 0
    y_scaled[is.na(y_scaled)] <- 0

    # Calculate LISA clustering
    lisa.clust <- as.character(interaction(x_scaled > 0, Wij %*% y_scaled > 0))
    lisa.clust <- gsub("TRUE", "High", lisa.clust)
    lisa.clust <- gsub("FALSE", "Low", lisa.clust)
    m$LISA <- lisa.clust
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
process_matrisome_expression <- function(gene_list, seurat_obj, name, counts_matrix, matrisome = NULL) {
  # Use the pre-fetched counts matrix passed as argument
  counts <- counts_matrix

  available_genes <- c()

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
  coords <- GetTissueCoordinates(seurat_obj)

  # Calculate metrics
  base_counts <- colSums(counts[available_genes, , drop = FALSE])
  log_scaled <- {
    sums_log <- log1p(base_counts)
    (sums_log - min(sums_log)) / (max(sums_log) - min(sums_log))
  }
  robust <- (base_counts - median(base_counts)) / IQR(base_counts)

  # Create data frame - handle VisiumV1 vs VisiumV2 coordinate structures
  if (all(c("imagerow", "imagecol") %in% colnames(coords))) {
    df <- data.frame(
      spot = names(base_counts),
      base_counts = base_counts,
      log_scaled = log_scaled,
      robust = robust,
      imagerow = coords$imagerow,
      imagecol = coords$imagecol
    )
  } else {
    # VisiumV2: use x/y as fallback
    df <- data.frame(
      spot = names(base_counts),
      base_counts = base_counts,
      log_scaled = log_scaled,
      robust = robust,
      imagerow = coords$x,
      imagecol = coords$y
    )
  }
  return(df)
}


#' Generate spatial plots for matrisome categories.
create_matrisome_plots <- function(seurat_obj, display_name, internal_name, type, annotation_col) {

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

  # --- Data Preparation (with Annotation) ---
  coords <- GetTissueCoordinates(seurat_obj)

  # Fetch scores AND the current annotation column
  plot_data <- FetchData(seurat_obj, vars = c(robust_feature_name, log_scaled_feature_name, annotation_col))
  df <- cbind(coords, plot_data)

  # Define the correct Seurat default palette
  seurat_spectral_palette <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
  seurat_color_ramp <- colorRamp(seurat_spectral_palette)

  # --- Create the interactive spatial distribution plot ---

  robust_values <- df[[robust_feature_name]]
  robust_min <- min(robust_values, na.rm = TRUE)
  robust_max <- max(robust_values, na.rm = TRUE)

  # Handle the edge case where all values are the same (prevents division by zero)
  if (robust_max == robust_min) {
    robust_colors <- rep(rgb(seurat_color_ramp(0), maxColorValue = 255), nrow(df))
  } else {
    robust_normalized <- (robust_values - robust_min) / (robust_max - robust_min)
    robust_colors <- apply(seurat_color_ramp(robust_normalized), 1, function(rgb_vals) {
      if(any(is.na(rgb_vals))) return("#808080")
      rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], maxColorValue = 255)
    })
  }

  robust_ticks <- pretty(c(robust_min, robust_max), n = 5)
  robust_tick_text <- round(robust_ticks, 1)

  interactive_spatial <- plot_ly(df, x = ~imagecol, y = ~imagerow, type = 'scatter', mode = 'markers',
                                 marker = list(
                                   color = robust_colors,
                                   colorbar = list(
                                     title = "",
                                     tickvals = robust_ticks,
                                     ticktext = robust_tick_text,
                                     len = 0.5,
                                     thickness = 15
                                   ),
                                   cmin = robust_min,
                                   cmax = robust_max,
                                   colorscale = seurat_spectral_palette,
                                   showscale = F,
                                   size = 5
                                 ),
                                 text = as.formula(paste0("~paste('Value:', round(`", robust_feature_name, "`, 3), '<br>Annotation:', `", annotation_col, "`)")),
                                 hoverinfo = 'text') %>%
    layout(xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
           yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE,
                        scaleanchor = "x", scaleratio = 1, autorange = "reversed"),
           plot_bgcolor = 'rgba(0,0,0,0)', paper_bgcolor = 'rgba(0,0,0,0)') %>%
    config(responsive = TRUE)

  # --- Create the interactive hotspot plot ---

  hotspot_values <- df[[log_scaled_feature_name]]
  hotspot_min <- min(hotspot_values, na.rm = TRUE)
  hotspot_max <- max(hotspot_values, na.rm = TRUE)

  if (hotspot_max == hotspot_min) {
    hotspot_colors <- rep(rgb(seurat_color_ramp(0), maxColorValue = 255), nrow(df))
  } else {
    hotspot_normalized <- (hotspot_values - hotspot_min) / (hotspot_max - hotspot_min)
    hotspot_colors <- apply(seurat_color_ramp(hotspot_normalized), 1, function(rgb_vals) {
      if(any(is.na(rgb_vals))) return("#808080")
      rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], maxColorValue = 255)
    })
  }

  hotspot_ticks <- pretty(c(hotspot_min, hotspot_max), n = 5)
  hotspot_tick_text <- round(hotspot_ticks, 1)

  interactive_hotspot <- plot_ly(df, x = ~imagecol, y = ~imagerow, type = 'scatter', mode = 'markers',
                                 marker = list(
                                   color = hotspot_colors,
                                   colorbar = list(
                                     title = "",
                                     tickvals = hotspot_ticks,
                                     ticktext = hotspot_tick_text,
                                     len = 0.5,
                                     thickness = 15
                                   ),
                                   cmin = hotspot_min,
                                   cmax = hotspot_max,
                                   colorscale = seurat_spectral_palette,
                                   showscale = F,
                                   size = 5
                                 ),
                                 text = as.formula(paste0("~paste('Value:', round(`", log_scaled_feature_name, "`, 3), '<br>Annotation:', `", annotation_col, "`)")),
                                 hoverinfo = 'text') %>%
    layout(xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
           yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE,
                        scaleanchor = "x", scaleratio = 1, autorange = "reversed"),
           plot_bgcolor = 'rgba(0,0,0,0)', paper_bgcolor = 'rgba(0,0,0,0)') %>%
    config(responsive = TRUE)

  # Return the full list of objects.
  return(list(
    seurat_object = seurat_obj,
    interactive_spatial = interactive_spatial,
    interactive_hotspot = interactive_hotspot,
    robust_feature = robust_feature_name,
    log_scaled_feature = log_scaled_feature_name
  ))
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
  if (!file.exists(extdata_path("ECM_domains_transformed4ScType.xlsx"))) {
    warning("ECM domain database not found. Skipping annotation.")
    return(seurat_obj)
  }
  # Prepare gene sets
  gs_list <- gene_sets_prepare(extdata_path("ECM_domains_transformed4ScType.xlsx"), "ECM")
  gs_list$gs_positive[["Vascular"]] <- NULL
  gs_list$gs_negative[["Vascular"]] <- NULL

  # Calculate scores
  es.max <- sctype_score(scRNAseqData = seurat_obj[["SCT"]]@scale.data, scaled = TRUE,
                         gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

  # Get top score for each spot
  spot_results <- do.call("rbind", lapply(rownames(seurat_obj@meta.data), function(coord) {
    es.max.subset <- es.max[, coord, drop = FALSE]
    es.max.coord <- sort(rowSums(as.matrix(es.max.subset)), decreasing = TRUE)
    head(data.frame(spot = coord, type = names(es.max.coord), scores = es.max.coord), 1)
  }))

  spot_results$type[as.numeric(as.character(spot_results$scores)) <= 0] <- "not.assigned"
  spot_results$type[grepl("Vascular", spot_results$type, ignore.case = TRUE)] <- "not.assigned"

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
    result <- suppressMessages(suppressWarnings(
      scCustomize::Updated_HGNC_Symbols(original_genes, verbose = FALSE, case_check_as_warn = TRUE)
    ))
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

  # SCTransform (only if needed for ScType and no existing SCT data)
  needs_sctype <- !"ecm_domain_annotation" %in% colnames(seurat_obj@meta.data)
  has_sct <- "SCT" %in% SeuratObject::Assays(seurat_obj)
  has_scale_data <- has_sct && nrow(seurat_obj[["SCT"]]@scale.data) > 0

  if (needs_sctype && !has_scale_data) {
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

  # --- Step 3: Check & Add ECM Domain Annotation ---
  if (!"ecm_domain_annotation" %in% colnames(seurat_obj@meta.data)) {
    showNotification("ECM domain annotation not found. Running ScType...", type = "message", duration = 8)
    seurat_obj <- sctype_annotate_ecm(seurat_obj)
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
create_feature_plot <- function(data, feature_col, feature_name, annotation_vector) {
  if (is.null(data) || !feature_col %in% colnames(data)) {
    return(plot_ly(type = 'scatter', mode = 'markers') %>%
             add_annotations(
               text = "No data available",
               showarrow = FALSE,
               font = list(size = 16)
             ))
  }

  colors <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  max_val <- max(data[[feature_col]], na.rm = TRUE)
  data$Annotation <- annotation_vector

  plot_ly(data,
          x = ~imagecol,
          y = ~imagerow,
          type = 'scatter',
          mode = 'markers',
          marker = list(
            size = 5,
            color = data[[feature_col]],
            colorscale = list(
              list(0, colors[1]),
              list(0.5, colors[6]),
              list(1, colors[11])
            ),
            showscale = TRUE,
            cmin = 0,
            cmax = max_val,
            colorbar = list(
              title = list(text = feature_name, font = list(size = 14)),
              thickness = 15, len = 0.5, y = 0.5
            )
          ),
          text = ~paste(
            feature_name, ":", round(.data[[feature_col]], 3),
            "<br>Annotation:", .data$Annotation
          ),
          hoverinfo = 'text') %>%
    layout(
      xaxis = list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
      yaxis = list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE,
                   scaleanchor = "x", scaleratio = 1, autorange = "reversed"),
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
    y_range <- range(data_df$feature1, na.rm = TRUE)
    p <- p %>% layout(yaxis = list(title = 'Expression', range = y_range))
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

  rng <- c(data_df$feature1, data_df$feature2)
  rng <- c(0, max(rng, na.rm = TRUE))

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
    p <- p %>% layout(yaxis = list(title = 'Expression', range = rng))
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

#' Find Enriched Markers in an Activity Assay
#'
#' This function calculates enrichment statistics (avg_log2FC and percentage
#' difference) for features in a given activity assay, comparing each group
#' (cluster) to all other groups.
#'
#' @param seurat_obj A Seurat object containing the activity assay.
#' @param assay_name A string specifying the name of the assay to analyze
#'        (e.g., "LRACTIVITY" or "LIGAND_ACTIVITY").
#' @param group_by A string specifying the metadata column to use for grouping
#'        cells (e.g., "seurat_clusters", "ecm_domain_annotation").
#'
#' @return A tidy dataframe with enrichment statistics for each feature in each group.
#'         Columns: feature, cluster, avg_log2FC, perc_difference, pct.1, pct.2
find_lr_markers <- function(seurat_obj, assay_name, group_by) {

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


#' Coalesce Reciprocal Interaction Pairs into Axes
#'
#' Identifies reciprocal pairs (e.g., A-B and B-A) in statistics and mean score
#' matrices and combines them into a single "Communication Axis" (e.g., A-B_AXIS).
#' Non-reciprocal pairs are kept as is.
#'
#' @param stats_df A long-format dataframe of statistics from `find_lr_markers`.
#'        Must contain 'feature' and 'cluster' columns.
#' @param means_matrix A wide-format matrix of mean scores, with features as
#'        rows and clusters as columns.
#'
#' @return A list containing the coalesced `stats` dataframe and `means` matrix.
#'
coalesce_reciprocal_pairs <- function(stats_df, means_matrix) {

  cat("Coalescing reciprocal interaction pairs into communication axes...\n")

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
  # Sort L-R pair components alphabetically to identify reciprocal pairs
  # Example: "TGFB1-TGFBR1" and "TGFBR1-TGFB1" both become "TGFB1-TGFBR1"
  original_features <- rownames(means_matrix)
  canonical_names <- sapply(original_features, function(name) {
    paste(sort(strsplit(name, "-")[[1]]), collapse = "-")
  })

  # STEP 3: Identify reciprocal pairs
  # Count occurrences - if canonical name appears twice, it's reciprocal
  name_counts <- table(canonical_names)
  reciprocal_ids <- names(name_counts[name_counts > 1])

  # STEP 4: Create new feature names
  # Add "_AXIS" suffix to reciprocal pairs to indicate bidirectional communication
  # Non-reciprocal pairs keep original names (unidirectional signaling)
  new_feature_names <- ifelse(
    canonical_names %in% reciprocal_ids,
    paste0(canonical_names, "_AXIS"),  # Bidirectional axis
    original_features                   # Unidirectional interaction
  )

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

#' Calculate Spatial Ligand-Receptor Activity Scores (Single-Core with Progress)
#'
#' This function calculates interaction scores on a single core, designed to be
#' called from within a Shiny observer with progress feedback (via shinybusy).
#' It avoids all parallel backend complexity to ensure stability in any environment.
#
#' @param seurat_obj A Seurat object with expression data.
#' @param lr_db A ligand-receptor database dataframe with columns 'Ligand' and 'Receptor'.
#' @param adj_matrix A sparse adjacency matrix defining spatial neighbors.
#' @param assay The assay to use for expression data (default: "SCT").
#' @param layer The layer to use within the assay (default: "data").
#' @param update_progress Optional progress update function (e.g., from shinybusy).
#'        Called at key steps with signature: update_progress(value, text)
#' @param chunk_size Number of interactions to process per batch (default: 5000).
#'        Lower values reduce peak memory usage but may be slightly slower.
#'
#' @return A Seurat object with the new "LRACTIVITY" assay added.
#'
compute_spatial_lr_scores_single_core <- function(seurat_obj, lr_db, adj_matrix, assay = "SCT", layer = "data", update_progress = NULL, chunk_size = 5000) {

  # 1. Prepare data
  expr <- LayerData(seurat_obj, assay = assay, layer = layer)

  # Ensure adjacency matrix matches expression matrix spot order
  if (!all(colnames(expr) == colnames(adj_matrix))) {
    adj_matrix <- adj_matrix[colnames(expr), colnames(expr)]
  }

  # Filter L-R database to only include genes present in this dataset
  genes_in_data <- rownames(expr)
  lr_db_filtered <- lr_db %>%
    filter(Ligand %in% genes_in_data & Receptor %in% genes_in_data)

  if (nrow(lr_db_filtered) == 0) stop("No valid ligand-receptor pairs found in expression data.")
  n_interactions <- nrow(lr_db_filtered)
  cat("Found", n_interactions, "valid interaction rows to score.\n")

  # Progress update: starting calculation
  if (!is.null(update_progress)) {
    update_progress(value = 0.2, text = "Preparing expression data...")
  }

  # 2. Prepare expression matrix with only LR genes
  lr_genes <- unique(c(lr_db_filtered$Ligand, lr_db_filtered$Receptor))
  expr_lr <- expr[lr_genes, , drop = FALSE]

  # Ensure sparse matrix format (dgCMatrix) for efficient operations
  if (!inherits(expr_lr, "dgCMatrix")) {
    expr_lr <- as(expr_lr, "dgCMatrix")
  }

  # Apply sqrt transformation directly on sparse matrix's data slot
  expr_lr@x <- sqrt(pmax(expr_lr@x, 0))

  # Sparsify adjacency matrix for efficient multiplication
  adj_matrix <- as(adj_matrix, "dgCMatrix")

  # 3. Compute neighbour-weighted expression for all genes at once
  if (!is.null(update_progress)) {
    update_progress(value = 0.3, text = "Computing neighbour expression...")
  }
  neighbor_expr <- t(adj_matrix %*% t(expr_lr))

  # Clean up original expression matrix to free memory
  rm(expr)
  gc(verbose = FALSE)

  # 4. Process interactions in chunks to control peak memory usage
  #    This prevents OOM on memory-constrained servers (e.g., shinyapps.io)
  n_chunks <- ceiling(n_interactions / chunk_size)
  activity_list <- vector("list", n_chunks)

  for (i in seq_len(n_chunks)) {
    # Calculate chunk boundaries
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, n_interactions)
    chunk_rows <- start_idx:end_idx

    # Update progress
    if (!is.null(update_progress)) {
      progress_val <- 0.3 + (0.6 * i / n_chunks)
      update_progress(value = progress_val, text = sprintf("Processing interactions (%d/%d)...", i, n_chunks))
    }

    # Get indices for this chunk
    lig_idx <- match(lr_db_filtered$Ligand[chunk_rows], rownames(expr_lr))
    rec_idx <- match(lr_db_filtered$Receptor[chunk_rows], rownames(expr_lr))

    # Compute activity for this chunk
    chunk_activity <- neighbor_expr[lig_idx, , drop = FALSE] * expr_lr[rec_idx, , drop = FALSE]
    rownames(chunk_activity) <- paste(lr_db_filtered$Ligand[chunk_rows], lr_db_filtered$Receptor[chunk_rows], sep = "-")

    activity_list[[i]] <- chunk_activity

    # Force garbage collection between chunks
    if (i < n_chunks) gc(verbose = FALSE)
  }

  # 5. Combine chunks into final activity matrix
  if (!is.null(update_progress)) {
    update_progress(value = 0.95, text = "Finalising results...")
  }

  activity <- do.call(rbind, activity_list)
  rm(activity_list, neighbor_expr, expr_lr)
  gc(verbose = FALSE)

  # 6. Store results in Seurat object
  seurat_obj[["LRACTIVITY"]] <- CreateAssayObject(counts = activity)
  seurat_obj <- SetAssayData(seurat_obj, assay = "LRACTIVITY", layer = "data", new.data = activity)

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
