#!/usr/bin/env Rscript

# Technical benchmark requested in response to Reviewer 1, Major Comment 1.
#
# Primary comparison:
#   1. Recompute ECM niche scores exactly as STapp_beta does for this object:
#      - use SCT when present;
#      - use the assay's data layer;
#      - update/filter signature symbols;
#      - call UCell::AddModuleScore_UCell(..., name = "_UCell").
#   2. Compare UCell with simple sums/means from the identical SCT expression
#      matrix and with literal raw UMI sums from the Spatial counts layer.
#
# Optional stress test (--downsample):
#   Apply binomial thinning to the raw Spatial count matrix and compare UCell,
#   raw UMI sums, and log-normalized means on the identical thinned matrices.
#   This is an algorithm-level capture-depth test, not biological validation.

args <- commandArgs(trailingOnly = TRUE)
full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
if (length(file_arg) != 1L) stop("Could not determine the script location.")
script_path <- normalizePath(
  sub("^--file=", "", file_arg),
  winslash = "/",
  mustWork = TRUE
)
script_dir <- dirname(script_path)

arg_value <- function(name, default) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) == 0L) return(default)
  sub(paste0("^", prefix), "", hit[[1]])
}

input_path <- normalizePath(
  arg_value(
    "input",
    file.path(script_dir, "data", "Breast_cancer_Ductal_Carcinoma_In_Situ_Invasive.RDS")
  ),
  winslash = "/",
  mustWork = TRUE
)
signature_path <- normalizePath(
  arg_value(
    "signatures",
    file.path(script_dir, "data", "ecm_ucell_signatures.rds")
  ),
  winslash = "/",
  mustWork = TRUE
)
hgnc_path <- normalizePath(
  arg_value(
    "hgnc",
    file.path(script_dir, "data", "hgnc_long_data.rds")
  ),
  winslash = "/",
  mustWork = TRUE
)
output_dir <- normalizePath(
  arg_value("output", file.path(script_dir, "results")),
  winslash = "/",
  mustWork = FALSE
)

run_downsampling <- "--downsample" %in% args
run_strength_tests <- "--strength-tests" %in% args
downsample_reps <- as.integer(arg_value("reps", "20"))
downsample_rates <- as.numeric(strsplit(arg_value("rates", "0.75,0.50,0.25"), ",")[[1]])
strength_reps <- as.integer(arg_value("strength-reps", "20"))
random_seed <- as.integer(arg_value("seed", "264998"))

if (!is.finite(downsample_reps) || downsample_reps < 1L) {
  stop("--reps must be a positive integer.")
}
if (any(!is.finite(downsample_rates)) || any(downsample_rates <= 0) ||
    any(downsample_rates >= 1)) {
  stop("Every --rates value must be strictly between 0 and 1.")
}
if (!is.finite(strength_reps) || strength_reps < 1L) {
  stop("--strength-reps must be a positive integer.")
}

required_packages <- c(
  "BiocParallel", "ggplot2", "Matrix", "Seurat", "SeuratObject", "UCell"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(random_seed)

message("Reading Seurat object...")
obj <- readRDS(input_path)
if (!inherits(obj, "Seurat")) stop("Input is not a Seurat object.")

signature_list <- readRDS(signature_path)
hgnc_long_data <- readRDS(hgnc_path)
if (!is.list(signature_list) || is.null(names(signature_list))) {
  stop("Signature file must contain a named list of gene vectors.")
}

# Exact copies of the relevant STapp_beta helper logic.
update_hgnc_symbols <- function(genes, hgnc_long_data) {
  approved <- unique(hgnc_long_data$symbol)
  rename_tbl <- hgnc_long_data[
    !hgnc_long_data$prev_symbol %in% approved,
    ,
    drop = FALSE
  ]
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

translate_gene_signatures <- function(seurat_obj, signatures, hgnc_data) {
  object_genes <- rownames(seurat_obj)
  translated <- list()
  for (signature_name in names(signatures)) {
    corrected <- update_hgnc_symbols(
      signatures[[signature_name]],
      hgnc_data
    )$Output_Features
    valid <- unique(corrected[corrected %in% object_genes])
    if (length(valid) > 0L) translated[[signature_name]] <- valid
  }
  translated
}

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || length(unique(x[keep])) < 2L ||
      length(unique(y[keep])) < 2L) {
    return(NA_real_)
  }
  unname(stats::cor(x[keep], y[keep], method = "spearman"))
}

top_indices <- function(x, fraction = 0.10) {
  keep <- which(is.finite(x))
  n_top <- max(1L, ceiling(length(keep) * fraction))
  keep[order(x[keep], decreasing = TRUE)[seq_len(n_top)]]
}

jaccard <- function(a, b) {
  union_n <- length(union(a, b))
  if (union_n == 0L) return(NA_real_)
  length(intersect(a, b)) / union_n
}

z_score <- function(x) {
  spread <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(spread) || spread == 0) return(rep(0, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / spread)
}

write_gzip_csv <- function(x, path) {
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  utils::write.csv(x, connection, row.names = FALSE)
}

assays <- SeuratObject::Assays(obj)
score_assay <- if ("SCT" %in% assays) "SCT" else SeuratObject::DefaultAssay(obj)
SeuratObject::DefaultAssay(obj) <- score_assay

translated_signatures <- translate_gene_signatures(
  obj,
  signature_list,
  hgnc_long_data
)
if (length(translated_signatures) == 0L) {
  stop("None of the signature genes are present in the Seurat object.")
}

coverage_rows <- lapply(names(signature_list), function(signature_name) {
  original <- unique(signature_list[[signature_name]])
  mapping <- update_hgnc_symbols(original, hgnc_long_data)
  corrected <- unique(mapping$Output_Features)
  detected <- intersect(corrected, rownames(obj))
  data.frame(
    signature = signature_name,
    requested_genes = length(original),
    corrected_unique_genes = length(corrected),
    detected_genes = length(detected),
    coverage_fraction = length(detected) / length(corrected),
    missing_genes = paste(setdiff(corrected, detected), collapse = ";"),
    stringsAsFactors = FALSE
  )
})
coverage <- do.call(rbind, coverage_rows)
utils::write.csv(
  coverage,
  file.path(output_dir, "signature_coverage.csv"),
  row.names = FALSE
)

score_columns <- paste0(names(translated_signatures), "_UCell")
precomputed_scores <- obj@meta.data[, intersect(score_columns, colnames(obj@meta.data)), drop = FALSE]

message("Recomputing UCell scores with the STapp_beta assay/layer behavior...")
scoring_obj <- obj
scoring_obj@meta.data <- scoring_obj@meta.data[
  ,
  setdiff(colnames(scoring_obj@meta.data), score_columns),
  drop = FALSE
]
scored_obj <- UCell::AddModuleScore_UCell(
  scoring_obj,
  features = translated_signatures,
  maxRank = 1500,
  chunk.size = 1000,
  BPPARAM = BiocParallel::SerialParam(),
  assay = score_assay,
  slot = "data",
  name = "_UCell"
)
ucell_scores <- scored_obj@meta.data[, score_columns, drop = FALSE]
ucell_scores <- ucell_scores[colnames(obj), , drop = FALSE]

score_expression <- SeuratObject::LayerData(
  obj,
  assay = score_assay,
  layer = "data"
)
score_expression <- score_expression[, colnames(obj), drop = FALSE]

raw_assay <- if ("Spatial" %in% assays) "Spatial" else score_assay
raw_counts <- SeuratObject::LayerData(
  obj,
  assay = raw_assay,
  layer = "counts"
)
raw_counts <- raw_counts[, colnames(obj), drop = FALSE]

library_size <- Matrix::colSums(raw_counts)
detected_features <- Matrix::colSums(raw_counts > 0)

spot_rows <- vector("list", length(translated_signatures))
comparison_rows <- vector("list", length(translated_signatures))

for (index in seq_along(translated_signatures)) {
  signature_name <- names(translated_signatures)[[index]]
  genes <- translated_signatures[[index]]
  ucell_column <- paste0(signature_name, "_UCell")

  ucell <- as.numeric(ucell_scores[[ucell_column]])
  same_input <- score_expression[genes, , drop = FALSE]
  simple_sum <- as.numeric(Matrix::colSums(same_input))
  simple_mean <- as.numeric(Matrix::colMeans(same_input))
  raw_genes <- intersect(genes, rownames(raw_counts))
  raw_sum <- as.numeric(Matrix::colSums(raw_counts[raw_genes, , drop = FALSE]))

  existing_correlation <- NA_real_
  existing_max_abs_difference <- NA_real_
  if (ucell_column %in% colnames(precomputed_scores)) {
    existing <- as.numeric(precomputed_scores[colnames(obj), ucell_column])
    existing_correlation <- safe_spearman(ucell, existing)
    existing_max_abs_difference <- max(abs(ucell - existing), na.rm = TRUE)
  }

  comparison_rows[[index]] <- data.frame(
    signature = signature_name,
    spots = ncol(obj),
    detected_signature_genes = length(genes),
    spearman_ucell_vs_sct_sum = safe_spearman(ucell, simple_sum),
    spearman_ucell_vs_sct_mean = safe_spearman(ucell, simple_mean),
    spearman_ucell_vs_raw_umi_sum = safe_spearman(ucell, raw_sum),
    spearman_ucell_vs_library_size = safe_spearman(ucell, library_size),
    spearman_sct_sum_vs_library_size = safe_spearman(simple_sum, library_size),
    spearman_raw_umi_sum_vs_library_size = safe_spearman(raw_sum, library_size),
    spearman_ucell_vs_detected_features = safe_spearman(ucell, detected_features),
    spearman_sct_sum_vs_detected_features = safe_spearman(simple_sum, detected_features),
    spearman_raw_umi_sum_vs_detected_features = safe_spearman(raw_sum, detected_features),
    top10_jaccard_ucell_vs_sct_sum = jaccard(
      top_indices(ucell),
      top_indices(simple_sum)
    ),
    top10_jaccard_ucell_vs_raw_umi_sum = jaccard(
      top_indices(ucell),
      top_indices(raw_sum)
    ),
    precomputed_ucell_spearman = existing_correlation,
    precomputed_ucell_max_abs_difference = existing_max_abs_difference,
    stringsAsFactors = FALSE
  )

  spot_rows[[index]] <- data.frame(
    spot_id = colnames(obj),
    signature = signature_name,
    library_size = as.numeric(library_size),
    detected_features = as.numeric(detected_features),
    ucell_stapp = ucell,
    sct_sum = simple_sum,
    sct_mean = simple_mean,
    raw_umi_sum = raw_sum,
    stringsAsFactors = FALSE
  )
}

method_comparison <- do.call(rbind, comparison_rows)
spot_scores <- do.call(rbind, spot_rows)
utils::write.csv(
  method_comparison,
  file.path(output_dir, "method_comparison.csv"),
  row.names = FALSE
)
write_gzip_csv(
  spot_scores,
  file.path(output_dir, "spot_scores.csv.gz")
)

# Scatterplots compare the exact STapp UCell score with both simple baselines.
scatter_rows <- lapply(split(spot_scores, spot_scores$signature), function(d) {
  rbind(
    data.frame(
      signature = d$signature,
      comparator = "SCT expression sum",
      comparator_score = d$sct_sum,
      ucell_score = d$ucell_stapp
    ),
    data.frame(
      signature = d$signature,
      comparator = "Raw UMI sum",
      comparator_score = d$raw_umi_sum,
      ucell_score = d$ucell_stapp
    )
  )
})
scatter_data <- do.call(rbind, scatter_rows)
scatter_plot <- ggplot2::ggplot(
  scatter_data,
  ggplot2::aes(x = comparator_score, y = ucell_score)
) +
  ggplot2::geom_point(alpha = 0.25, size = 0.6, color = "#28536B") +
  ggplot2::facet_grid(signature ~ comparator, scales = "free") +
  ggplot2::labs(
    x = "Simple aggregate-expression score",
    y = "STapp UCell score",
    title = "UCell versus simple ECM-signature aggregation"
  ) +
  ggplot2::theme_bw(base_size = 9)
ggplot2::ggsave(
  file.path(output_dir, "score_comparison.pdf"),
  scatter_plot,
  width = 8.2,
  height = 5.8,
  units = "in"
)

# Spatial maps use standardized values only for visualization; all numerical
# comparisons above use the untransformed scores.
coordinates <- tryCatch(
  Seurat::GetTissueCoordinates(obj),
  error = function(e) NULL
)
if (is.null(coordinates) && length(obj@images) > 0L) {
  # Seurat 4.3.3 cannot dispatch GetTissueCoordinates() directly on some
  # VisiumV2 image objects. Their centroid boundary retains the spot table.
  image_attributes <- attributes(obj@images[[1]])
  centroids <- tryCatch(
    image_attributes$boundaries$centroids,
    error = function(e) NULL
  )
  if (!is.null(centroids)) {
    coordinates <- tryCatch(
      Seurat::GetTissueCoordinates(centroids),
      error = function(e) NULL
    )
    if (!is.null(coordinates) && "cell" %in% colnames(coordinates)) {
      rownames(coordinates) <- coordinates$cell
    }
  }
}
if (!is.null(coordinates) && nrow(coordinates) > 0L) {
  if (all(c("imagecol", "imagerow") %in% colnames(coordinates))) {
    coordinate_table <- data.frame(
      spot_id = rownames(coordinates),
      x = coordinates$imagecol,
      y = -coordinates$imagerow,
      stringsAsFactors = FALSE
    )
  } else {
    numeric_columns <- which(vapply(coordinates, is.numeric, logical(1)))
    if (length(numeric_columns) >= 2L) {
      coordinate_table <- data.frame(
        spot_id = rownames(coordinates),
        x = coordinates[[numeric_columns[[1]]]],
        y = -coordinates[[numeric_columns[[2]]]],
        stringsAsFactors = FALSE
      )
    } else {
      coordinate_table <- NULL
    }
  }

  if (!is.null(coordinate_table)) {
    map_rows <- lapply(split(spot_scores, spot_scores$signature), function(d) {
      rbind(
        data.frame(
          spot_id = d$spot_id,
          signature = d$signature,
          method = "UCell (STapp)",
          score_z = z_score(d$ucell_stapp)
        ),
        data.frame(
          spot_id = d$spot_id,
          signature = d$signature,
          method = "SCT expression sum",
          score_z = z_score(d$sct_sum)
        ),
        data.frame(
          spot_id = d$spot_id,
          signature = d$signature,
          method = "Raw UMI sum",
          score_z = z_score(d$raw_umi_sum)
        )
      )
    })
    map_data <- merge(
      do.call(rbind, map_rows),
      coordinate_table,
      by = "spot_id",
      all.x = TRUE,
      sort = FALSE
    )
    map_plot <- ggplot2::ggplot(
      map_data,
      ggplot2::aes(x = x, y = y, color = score_z)
    ) +
      ggplot2::geom_point(size = 0.55) +
      ggplot2::facet_grid(signature ~ method) +
      ggplot2::coord_equal() +
      ggplot2::scale_color_gradient2(
        low = "#2166AC",
        mid = "#F7F7F7",
        high = "#B2182B",
        midpoint = 0,
        name = "Score z"
      ) +
      ggplot2::labs(title = "Spatial comparison of ECM-signature scores") +
      ggplot2::theme_void(base_size = 9) +
      ggplot2::theme(
        strip.text = ggplot2::element_text(color = "black"),
        plot.title = ggplot2::element_text(hjust = 0.5)
      )
    ggplot2::ggsave(
      file.path(output_dir, "spatial_score_comparison.pdf"),
      map_plot,
      width = 9.2,
      height = 6.0,
      units = "in"
    )
  }
}

thin_counts <- function(count_matrix, retention) {
  thinned <- methods::as(count_matrix, "dgCMatrix")
  thinned@x <- as.numeric(stats::rbinom(
    length(thinned@x),
    size = as.integer(round(thinned@x)),
    prob = retention
  ))
  Matrix::drop0(thinned)
}

score_raw_matrix <- function(count_matrix, signatures) {
  ucell <- as.data.frame(
    UCell::ScoreSignatures_UCell(
      count_matrix,
      features = signatures,
      maxRank = 1500,
      chunk.size = 1000,
      BPPARAM = BiocParallel::SerialParam()
    ),
    check.names = FALSE
  )
  library <- Matrix::colSums(count_matrix)
  safe_library <- pmax(library, 1)
  normalized <- count_matrix %*% Matrix::Diagonal(x = 10000 / safe_library)
  normalized <- log1p(normalized)

  output <- list()
  for (signature_name in names(signatures)) {
    genes <- intersect(signatures[[signature_name]], rownames(count_matrix))
    output[[signature_name]] <- list(
      ucell = as.numeric(ucell[[paste0(signature_name, "_UCell")]]),
      raw_sum = as.numeric(Matrix::colSums(count_matrix[genes, , drop = FALSE])),
      lognorm_mean = as.numeric(Matrix::colMeans(normalized[genes, , drop = FALSE]))
    )
  }
  output
}

if (run_downsampling) {
  message(
    "Running optional downsampling: ",
    downsample_reps,
    " repetitions at retention rates ",
    paste(downsample_rates, collapse = ", ")
  )
  full_depth <- score_raw_matrix(raw_counts, translated_signatures)
  simulation_rows <- list()
  row_index <- 1L

  for (retention in downsample_rates) {
    for (replicate_id in seq_len(downsample_reps)) {
      thinned <- thin_counts(raw_counts, retention)
      reduced <- score_raw_matrix(thinned, translated_signatures)
      for (signature_name in names(translated_signatures)) {
        for (method_name in c("ucell", "raw_sum", "lognorm_mean")) {
          original <- full_depth[[signature_name]][[method_name]]
          perturbed <- reduced[[signature_name]][[method_name]]
          simulation_rows[[row_index]] <- data.frame(
            signature = signature_name,
            method = method_name,
            retention = retention,
            replicate = replicate_id,
            spearman_vs_full_depth = safe_spearman(original, perturbed),
            top10_jaccard_vs_full_depth = jaccard(
              top_indices(original),
              top_indices(perturbed)
            ),
            stringsAsFactors = FALSE
          )
          row_index <- row_index + 1L
        }
      }
    }
  }

  simulation_results <- do.call(rbind, simulation_rows)
  utils::write.csv(
    simulation_results,
    file.path(output_dir, "downsampling_replicates.csv"),
    row.names = FALSE
  )

  groups <- split(
    simulation_results,
    interaction(
      simulation_results$signature,
      simulation_results$method,
      simulation_results$retention,
      drop = TRUE
    )
  )
  downsampling_summary <- do.call(rbind, lapply(groups, function(d) {
    data.frame(
      signature = d$signature[[1]],
      method = d$method[[1]],
      retention = d$retention[[1]],
      repetitions = nrow(d),
      spearman_median = stats::median(d$spearman_vs_full_depth, na.rm = TRUE),
      spearman_q025 = stats::quantile(
        d$spearman_vs_full_depth,
        0.025,
        na.rm = TRUE,
        names = FALSE
      ),
      spearman_q975 = stats::quantile(
        d$spearman_vs_full_depth,
        0.975,
        na.rm = TRUE,
        names = FALSE
      ),
      top10_jaccard_median = stats::median(
        d$top10_jaccard_vs_full_depth,
        na.rm = TRUE
      ),
      top10_jaccard_q025 = stats::quantile(
        d$top10_jaccard_vs_full_depth,
        0.025,
        na.rm = TRUE,
        names = FALSE
      ),
      top10_jaccard_q975 = stats::quantile(
        d$top10_jaccard_vs_full_depth,
        0.975,
        na.rm = TRUE,
        names = FALSE
      ),
      stringsAsFactors = FALSE
    )
  }))
  rownames(downsampling_summary) <- NULL
  utils::write.csv(
    downsampling_summary,
    file.path(output_dir, "downsampling_summary.csv"),
    row.names = FALSE
  )

  downsampling_plot <- ggplot2::ggplot(
    downsampling_summary,
    ggplot2::aes(
      x = retention,
      y = spearman_median,
      color = method,
      group = method
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = spearman_q025, ymax = spearman_q975, fill = method),
      alpha = 0.15,
      color = NA
    ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(~signature) +
    ggplot2::scale_x_continuous(
      breaks = sort(unique(downsampling_summary$retention)),
      labels = function(x) paste0(round(100 * x), "%")
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x = "UMIs retained",
      y = "Spearman correlation with full-depth score",
      color = "Method",
      fill = "Method",
      title = "Score stability under controlled count downsampling"
    ) +
    ggplot2::theme_bw(base_size = 9)
  ggplot2::ggsave(
    file.path(output_dir, "downsampling_stability.pdf"),
    downsampling_plot,
    width = 7.4,
    height = 4.5,
    units = "in"
  )
}

if (run_strength_tests) {
  message("Running focused UCell strength tests...")

  rank_fraction <- function(x) {
    rank(x, ties.method = "average", na.last = "keep") / sum(is.finite(x))
  }

  stability_row <- function(reference, perturbed) {
    data.frame(
      spearman = safe_spearman(reference, perturbed),
      top10_jaccard = jaccard(top_indices(reference), top_indices(perturbed)),
      median_abs_percentile_rank_error = stats::median(
        abs(rank_fraction(reference) - rank_fraction(perturbed)),
        na.rm = TRUE
      ),
      max_abs_score_change = max(abs(reference - perturbed), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  # 1. maxRank sensitivity. The tested values span the STapp default and the
  # approximate feature-detection depth of this Visium section.
  maxrank_values <- c(1500L, 3000L, 6000L, 9000L)
  maxrank_scores <- list()
  for (maxrank in maxrank_values) {
    message("  maxRank = ", maxrank)
    maxrank_scores[[as.character(maxrank)]] <- as.data.frame(
      UCell::ScoreSignatures_UCell(
        score_expression,
        features = translated_signatures,
        maxRank = maxrank,
        chunk.size = 1000,
        BPPARAM = BiocParallel::SerialParam()
      ),
      check.names = FALSE
    )
  }

  maxrank_spot_rows <- list()
  maxrank_comparison_rows <- list()
  row_index <- 1L
  comparison_index <- 1L
  for (maxrank in maxrank_values) {
    for (signature_name in names(translated_signatures)) {
      column <- paste0(signature_name, "_UCell")
      values <- as.numeric(maxrank_scores[[as.character(maxrank)]][[column]])
      maxrank_spot_rows[[row_index]] <- data.frame(
        spot_id = colnames(obj),
        signature = signature_name,
        maxRank = maxrank,
        score = values,
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L

      reference <- as.numeric(maxrank_scores[["1500"]][[column]])
      metrics <- stability_row(reference, values)
      maxrank_comparison_rows[[comparison_index]] <- cbind(
        data.frame(
          signature = signature_name,
          reference_maxRank = 1500L,
          tested_maxRank = maxrank,
          stringsAsFactors = FALSE
        ),
        metrics
      )
      comparison_index <- comparison_index + 1L
    }
  }
  maxrank_spot_scores <- do.call(rbind, maxrank_spot_rows)
  maxrank_comparison <- do.call(rbind, maxrank_comparison_rows)
  write_gzip_csv(
    maxrank_spot_scores,
    file.path(output_dir, "strength_maxrank_spot_scores.csv.gz")
  )
  utils::write.csv(
    maxrank_comparison,
    file.path(output_dir, "strength_maxrank_comparison.csv"),
    row.names = FALSE
  )
  detected_quantiles <- stats::quantile(
    detected_features,
    probs = c(0, 0.05, 0.25, 0.50, 0.75, 0.95, 1),
    names = FALSE
  )
  utils::write.csv(
    data.frame(
      quantile = c(0, 0.05, 0.25, 0.50, 0.75, 0.95, 1),
      detected_features = as.numeric(detected_quantiles)
    ),
    file.path(output_dir, "strength_detected_feature_quantiles.csv"),
    row.names = FALSE
  )

  # 2. Dominant-gene sensitivity. Quantify how much of each aggregate is
  # supplied by its largest genes, amplify the three largest contributors,
  # and perform a complete leave-one-signature-gene-out analysis.
  dominance_rows <- list()
  amplification_rows <- list()
  leave_one_out_rows <- list()
  dominance_index <- 1L
  amplification_index <- 1L
  loo_index <- 1L

  for (signature_name in names(translated_signatures)) {
    message("  dominant-gene test: ", signature_name)
    genes <- translated_signatures[[signature_name]]
    sct_matrix <- score_expression[genes, , drop = FALSE]
    raw_genes <- intersect(genes, rownames(raw_counts))
    raw_matrix <- raw_counts[raw_genes, , drop = FALSE]

    summarize_concentration <- function(x, representation) {
      dense <- as.matrix(x)
      totals <- colSums(dense)
      ordered <- apply(dense, 2, sort, decreasing = TRUE)
      top1 <- ifelse(totals > 0, ordered[1, ] / totals, NA_real_)
      top3 <- ifelse(
        totals > 0,
        colSums(ordered[seq_len(min(3L, nrow(ordered))), , drop = FALSE]) / totals,
        NA_real_
      )
      data.frame(
        signature = signature_name,
        representation = representation,
        spots_with_positive_total = sum(totals > 0),
        median_top1_fraction = stats::median(top1, na.rm = TRUE),
        q90_top1_fraction = stats::quantile(top1, 0.90, na.rm = TRUE, names = FALSE),
        median_top3_fraction = stats::median(top3, na.rm = TRUE),
        q90_top3_fraction = stats::quantile(top3, 0.90, na.rm = TRUE, names = FALSE),
        stringsAsFactors = FALSE
      )
    }
    dominance_rows[[dominance_index]] <- summarize_concentration(
      sct_matrix,
      "SCT_data"
    )
    dominance_index <- dominance_index + 1L
    dominance_rows[[dominance_index]] <- summarize_concentration(
      raw_matrix,
      "raw_UMI"
    )
    dominance_index <- dominance_index + 1L

    contribution_order <- genes[order(
      as.numeric(Matrix::rowSums(sct_matrix)),
      decreasing = TRUE
    )]
    amplified_genes <- head(contribution_order, 3L)
    baseline_ucell <- as.numeric(
      maxrank_scores[["1500"]][[paste0(signature_name, "_UCell")]]
    )
    baseline_sct_mean <- as.numeric(Matrix::colMeans(sct_matrix))
    baseline_raw_sum <- as.numeric(Matrix::colSums(raw_matrix))

    for (gene in amplified_genes) {
      for (factor_value in c(2, 5, 10)) {
        changed_sct <- score_expression
        changed_sct[gene, ] <- changed_sct[gene, ] * factor_value
        changed_ucell_table <- as.data.frame(
          UCell::ScoreSignatures_UCell(
            changed_sct,
            features = setNames(list(genes), signature_name),
            maxRank = 1500,
            chunk.size = 1000,
            BPPARAM = BiocParallel::SerialParam()
          ),
          check.names = FALSE
        )
        changed_ucell <- as.numeric(
          changed_ucell_table[[paste0(signature_name, "_UCell")]]
        )
        changed_sct_mean <- as.numeric(Matrix::colMeans(
          changed_sct[genes, , drop = FALSE]
        ))
        changed_raw <- raw_counts
        changed_raw[gene, ] <- changed_raw[gene, ] * factor_value
        changed_raw_sum <- as.numeric(Matrix::colSums(
          changed_raw[raw_genes, , drop = FALSE]
        ))

        method_values <- list(
          ucell = cbind(baseline_ucell, changed_ucell),
          sct_mean = cbind(baseline_sct_mean, changed_sct_mean),
          raw_sum = cbind(baseline_raw_sum, changed_raw_sum)
        )
        for (method_name in names(method_values)) {
          metrics <- stability_row(
            method_values[[method_name]][, 1],
            method_values[[method_name]][, 2]
          )
          amplification_rows[[amplification_index]] <- cbind(
            data.frame(
              signature = signature_name,
              amplified_gene = gene,
              amplification_factor = factor_value,
              method = method_name,
              stringsAsFactors = FALSE
            ),
            metrics
          )
          amplification_index <- amplification_index + 1L
        }
        rm(changed_sct, changed_raw)
        invisible(gc())
      }
    }

    loo_signatures <- lapply(genes, function(gene) setdiff(genes, gene))
    names(loo_signatures) <- paste0("minus_", genes)
    loo_ucell <- as.data.frame(
      UCell::ScoreSignatures_UCell(
        score_expression,
        features = loo_signatures,
        maxRank = 1500,
        chunk.size = 1000,
        BPPARAM = BiocParallel::SerialParam()
      ),
      check.names = FALSE
    )
    for (gene in genes) {
      remaining <- setdiff(genes, gene)
      raw_remaining <- intersect(remaining, rownames(raw_counts))
      perturbed_values <- list(
        ucell = as.numeric(loo_ucell[[paste0("minus_", gene, "_UCell")]]),
        sct_mean = as.numeric(Matrix::colMeans(
          score_expression[remaining, , drop = FALSE]
        )),
        raw_sum = as.numeric(Matrix::colSums(
          raw_counts[raw_remaining, , drop = FALSE]
        ))
      )
      reference_values <- list(
        ucell = baseline_ucell,
        sct_mean = baseline_sct_mean,
        raw_sum = baseline_raw_sum
      )
      for (method_name in names(reference_values)) {
        metrics <- stability_row(
          reference_values[[method_name]],
          perturbed_values[[method_name]]
        )
        leave_one_out_rows[[loo_index]] <- cbind(
          data.frame(
            signature = signature_name,
            removed_gene = gene,
            method = method_name,
            stringsAsFactors = FALSE
          ),
          metrics
        )
        loo_index <- loo_index + 1L
      }
    }
  }
  utils::write.csv(
    do.call(rbind, dominance_rows),
    file.path(output_dir, "strength_gene_concentration.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    do.call(rbind, amplification_rows),
    file.path(output_dir, "strength_dominant_gene_amplification.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    do.call(rbind, leave_one_out_rows),
    file.path(output_dir, "strength_leave_one_gene_out.csv"),
    row.names = FALSE
  )

  # 3. Dataset-composition dependence. Score the same 500 evaluation spots
  # alone, in three distinct biological contexts, and in the full section.
  # Only the other columns supplied to the scoring procedure change.
  if (!"cell_type" %in% colnames(obj@meta.data)) {
    stop("The composition test requires the sample's cell_type annotation.")
  }
  set.seed(random_seed)
  evaluation_candidates <- rownames(obj@meta.data)[
    as.character(obj@meta.data$cell_type) == "Unassigned"
  ]
  evaluation_spots <- sample(
    evaluation_candidates,
    min(500L, length(evaluation_candidates)),
    replace = FALSE
  )
  cell_type <- as.character(obj@meta.data$cell_type)
  names(cell_type) <- rownames(obj@meta.data)
  context_groups <- list(
    evaluation_only = character(0),
    evaluation_plus_malignant = names(cell_type)[cell_type %in% c(
      "DCIS1", "DCIS2", "Invasive"
    )],
    evaluation_plus_stromal = names(cell_type)[cell_type %in% c(
      "Adipocytes", "Myoepithelial", "Stromal"
    )],
    evaluation_plus_immune = names(cell_type)[cell_type == "Immune"],
    full_section = colnames(obj)
  )
  composition_score_rows <- list()
  composition_index <- 1L
  for (context_name in names(context_groups)) {
    message("  composition context: ", context_name)
    context_cells <- if (context_name == "full_section") {
      colnames(obj)
    } else {
      union(evaluation_spots, context_groups[[context_name]])
    }
    context_cells <- intersect(colnames(obj), context_cells)
    context_expression <- score_expression[, context_cells, drop = FALSE]
    context_raw <- raw_counts[, context_cells, drop = FALSE]
    context_ucell <- as.data.frame(
      UCell::ScoreSignatures_UCell(
        context_expression,
        features = translated_signatures,
        maxRank = 1500,
        chunk.size = 1000,
        BPPARAM = BiocParallel::SerialParam()
      ),
      check.names = FALSE
    )
    rownames(context_ucell) <- context_cells

    context_obj <- subset(obj, cells = context_cells)
    SeuratObject::DefaultAssay(context_obj) <- score_assay
    context_obj <- Seurat::AddModuleScore(
      context_obj,
      features = translated_signatures,
      assay = score_assay,
      name = "CompositionAMS",
      seed = random_seed,
      nbin = 16,
      slot = "data"
    )
    for (signature_index in seq_along(translated_signatures)) {
      signature_name <- names(translated_signatures)[[signature_index]]
      genes <- translated_signatures[[signature_name]]
      raw_genes <- intersect(genes, rownames(context_raw))
      method_scores <- list(
        ucell = context_ucell[evaluation_spots, paste0(signature_name, "_UCell")],
        addmodulescore = context_obj@meta.data[
          evaluation_spots,
          paste0("CompositionAMS", signature_index)
        ],
        sct_mean = Matrix::colMeans(
          context_expression[genes, evaluation_spots, drop = FALSE]
        ),
        raw_sum = Matrix::colSums(
          context_raw[raw_genes, evaluation_spots, drop = FALSE]
        )
      )
      for (method_name in names(method_scores)) {
        composition_score_rows[[composition_index]] <- data.frame(
          spot_id = evaluation_spots,
          signature = signature_name,
          context = context_name,
          context_spots = length(context_cells),
          method = method_name,
          score = as.numeric(method_scores[[method_name]]),
          stringsAsFactors = FALSE
        )
        composition_index <- composition_index + 1L
      }
    }
    rm(context_obj, context_expression, context_raw, context_ucell)
    invisible(gc())
  }
  composition_scores <- do.call(rbind, composition_score_rows)
  composition_reference <- composition_scores[
    composition_scores$context == "full_section",
    ,
    drop = FALSE
  ]
  composition_summary_rows <- list()
  composition_summary_index <- 1L
  for (signature_name in unique(composition_scores$signature)) {
    for (method_name in unique(composition_scores$method)) {
      reference <- composition_reference[
        composition_reference$signature == signature_name &
          composition_reference$method == method_name,
        ,
        drop = FALSE
      ]
      reference <- reference[match(evaluation_spots, reference$spot_id), ]
      for (context_name in setdiff(unique(composition_scores$context), "full_section")) {
        tested <- composition_scores[
          composition_scores$signature == signature_name &
            composition_scores$method == method_name &
            composition_scores$context == context_name,
          ,
          drop = FALSE
        ]
        tested <- tested[match(evaluation_spots, tested$spot_id), ]
        metrics <- stability_row(reference$score, tested$score)
        composition_summary_rows[[composition_summary_index]] <- cbind(
          data.frame(
            signature = signature_name,
            method = method_name,
            reference_context = "full_section",
            tested_context = context_name,
            evaluation_spots = length(evaluation_spots),
            tested_context_spots = tested$context_spots[[1]],
            median_score_shift = stats::median(
              tested$score - reference$score,
              na.rm = TRUE
            ),
            stringsAsFactors = FALSE
          ),
          metrics
        )
        composition_summary_index <- composition_summary_index + 1L
      }
    }
  }
  write_gzip_csv(
    composition_scores,
    file.path(output_dir, "strength_composition_spot_scores.csv.gz")
  )
  utils::write.csv(
    do.call(rbind, composition_summary_rows),
    file.path(output_dir, "strength_composition_summary.csv"),
    row.names = FALSE
  )

  # 4. Heterogeneous capture efficiency. Each replicate independently assigns
  # 25%, 50%, 75%, or 100% retention to individual spots, then thins counts.
  thin_counts_by_spot <- function(count_matrix, retention_by_spot) {
    thinned <- methods::as(count_matrix, "dgCMatrix")
    entry_columns <- rep.int(
      seq_len(ncol(thinned)),
      diff(thinned@p)
    )
    thinned@x <- as.numeric(stats::rbinom(
      length(thinned@x),
      size = as.integer(round(thinned@x)),
      prob = retention_by_spot[entry_columns]
    ))
    Matrix::drop0(thinned)
  }
  message("  heterogeneous capture test: ", strength_reps, " repetitions")
  retention_levels <- c(0.25, 0.50, 0.75, 1.00)
  full_depth <- score_raw_matrix(raw_counts, translated_signatures)
  heterogeneous_rows <- list()
  heterogeneous_group_rows <- list()
  hetero_index <- 1L
  hetero_group_index <- 1L
  for (replicate_id in seq_len(strength_reps)) {
    set.seed(random_seed + replicate_id)
    retention <- sample(
      rep(retention_levels, length.out = ncol(raw_counts)),
      ncol(raw_counts),
      replace = FALSE
    )
    thinned <- thin_counts_by_spot(raw_counts, retention)
    reduced <- score_raw_matrix(thinned, translated_signatures)
    for (signature_name in names(translated_signatures)) {
      for (method_name in c("ucell", "raw_sum", "lognorm_mean")) {
        reference <- full_depth[[signature_name]][[method_name]]
        perturbed <- reduced[[signature_name]][[method_name]]
        metrics <- stability_row(reference, perturbed)
        heterogeneous_rows[[hetero_index]] <- cbind(
          data.frame(
            signature = signature_name,
            method = method_name,
            replicate = replicate_id,
            score_change_vs_retention_spearman = safe_spearman(
              perturbed - reference,
              retention
            ),
            stringsAsFactors = FALSE
          ),
          metrics
        )
        hetero_index <- hetero_index + 1L

        reference_rank <- rank_fraction(reference)
        perturbed_rank <- rank_fraction(perturbed)
        for (retention_value in retention_levels) {
          keep <- retention == retention_value
          heterogeneous_group_rows[[hetero_group_index]] <- data.frame(
            signature = signature_name,
            method = method_name,
            replicate = replicate_id,
            retention = retention_value,
            spots = sum(keep),
            within_group_spearman = safe_spearman(
              reference[keep],
              perturbed[keep]
            ),
            median_abs_percentile_rank_error = stats::median(
              abs(reference_rank[keep] - perturbed_rank[keep]),
              na.rm = TRUE
            ),
            median_percentile_rank_shift = stats::median(
              perturbed_rank[keep] - reference_rank[keep],
              na.rm = TRUE
            ),
            stringsAsFactors = FALSE
          )
          hetero_group_index <- hetero_group_index + 1L
        }
      }
    }
    rm(thinned, reduced)
    invisible(gc())
  }
  heterogeneous_replicates <- do.call(rbind, heterogeneous_rows)
  heterogeneous_by_depth <- do.call(rbind, heterogeneous_group_rows)
  utils::write.csv(
    heterogeneous_replicates,
    file.path(output_dir, "strength_heterogeneous_capture_replicates.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    heterogeneous_by_depth,
    file.path(output_dir, "strength_heterogeneous_capture_by_depth.csv"),
    row.names = FALSE
  )

  summarize_numeric_groups <- function(data, group_columns, value_columns) {
    groups <- split(
      data,
      interaction(data[, group_columns, drop = FALSE], drop = TRUE)
    )
    do.call(rbind, lapply(groups, function(d) {
      identifiers <- d[1, group_columns, drop = FALSE]
      summaries <- list()
      for (column in value_columns) {
        summaries[[paste0(column, "_median")]] <- stats::median(
          d[[column]],
          na.rm = TRUE
        )
        summaries[[paste0(column, "_q025")]] <- stats::quantile(
          d[[column]], 0.025, na.rm = TRUE, names = FALSE
        )
        summaries[[paste0(column, "_q975")]] <- stats::quantile(
          d[[column]], 0.975, na.rm = TRUE, names = FALSE
        )
      }
      cbind(identifiers, as.data.frame(summaries, check.names = FALSE))
    }))
  }
  heterogeneous_summary <- summarize_numeric_groups(
    heterogeneous_replicates,
    c("signature", "method"),
    c(
      "spearman",
      "top10_jaccard",
      "median_abs_percentile_rank_error",
      "score_change_vs_retention_spearman"
    )
  )
  heterogeneous_depth_summary <- summarize_numeric_groups(
    heterogeneous_by_depth,
    c("signature", "method", "retention"),
    c(
      "within_group_spearman",
      "median_abs_percentile_rank_error",
      "median_percentile_rank_shift"
    )
  )
  rownames(heterogeneous_summary) <- NULL
  rownames(heterogeneous_depth_summary) <- NULL
  utils::write.csv(
    heterogeneous_summary,
    file.path(output_dir, "strength_heterogeneous_capture_summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    heterogeneous_depth_summary,
    file.path(output_dir, "strength_heterogeneous_capture_depth_summary.csv"),
    row.names = FALSE
  )
}

run_metadata <- data.frame(
  field = c(
    "timestamp",
    "input_path",
    "signature_path",
    "hgnc_path",
    "spots",
    "genes",
    "score_assay",
    "score_layer",
    "raw_assay",
    "ucell_version",
    "seurat_version",
    "random_seed",
    "downsampling_run",
    "downsampling_repetitions",
    "downsampling_retention_rates",
    "strength_tests_run",
    "strength_test_repetitions"
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    input_path,
    signature_path,
    hgnc_path,
    ncol(obj),
    nrow(obj),
    score_assay,
    "data",
    raw_assay,
    as.character(utils::packageVersion("UCell")),
    as.character(utils::packageVersion("Seurat")),
    random_seed,
    run_downsampling,
    if (run_downsampling) downsample_reps else NA,
    if (run_downsampling) paste(downsample_rates, collapse = ";") else NA,
    run_strength_tests,
    if (run_strength_tests) strength_reps else NA
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(
  run_metadata,
  file.path(output_dir, "run_metadata.csv"),
  row.names = FALSE
)

message("Benchmark complete. Results written to: ", output_dir)
