#!/usr/bin/env Rscript

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- normalizePath(
  sub("^--file=", "", file_arg[[1]]),
  winslash = "/",
  mustWork = TRUE
)
script_dir <- dirname(script_path)
results_dir <- file.path(script_dir, "results")

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 is required.")
}

method_labels <- c(
  ucell = "UCell",
  sct_mean = "Summed normalized expression",
  lognorm_mean = "Log-normalized mean",
  raw_sum = "Raw UMI sum",
  addmodulescore = "AddModuleScore"
)

method_colors <- c(
  "UCell" = "#0072B2",
  "Summed normalized expression" = "#009E73",
  "Log-normalized mean" = "#009E73",
  "Raw UMI sum" = "#D55E00",
  "AddModuleScore" = "#CC79A7"
)

base_theme <- ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom",
    strip.background = ggplot2::element_rect(fill = "grey94"),
    plot.title.position = "plot"
  )

save_figure <- function(plot, stem, width, height) {
  ggplot2::ggsave(
    file.path(results_dir, paste0(stem, ".pdf")),
    plot,
    width = width,
    height = height,
    units = "in"
  )
  ggplot2::ggsave(
    file.path(results_dir, paste0(stem, ".png")),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300
  )
}

# Figure 1: maxRank sensitivity relative to the STapp value of 1500.
maxrank <- utils::read.csv(
  file.path(results_dir, "strength_maxrank_comparison.csv"),
  stringsAsFactors = FALSE
)
maxrank_long <- rbind(
  data.frame(
    signature = maxrank$signature,
    maxRank = maxrank$tested_maxRank,
    metric = "Spearman correlation",
    value = maxrank$spearman
  ),
  data.frame(
    signature = maxrank$signature,
    maxRank = maxrank$tested_maxRank,
    metric = "Top-10% hotspot Jaccard",
    value = maxrank$top10_jaccard
  )
)
maxrank_long$signature <- factor(
  maxrank_long$signature,
  levels = c("Interstitial", "Basement")
)
p_maxrank <- ggplot2::ggplot(
  maxrank_long,
  ggplot2::aes(x = maxRank, y = value, color = signature, shape = signature)
) +
  ggplot2::geom_hline(yintercept = 1, linewidth = 0.35, color = "grey70") +
  ggplot2::geom_vline(
    xintercept = 5712,
    linewidth = 0.45,
    linetype = "dashed",
    color = "grey45"
  ) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.3) +
  ggplot2::facet_wrap(~metric, ncol = 1) +
  ggplot2::scale_x_continuous(breaks = c(1500, 3000, 6000, 9000)) +
  ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  ggplot2::labs(
    x = "UCell maxRank",
    y = "Agreement with maxRank = 1500",
    color = "Signature",
    shape = "Signature",
    title = "ECM UCell scores are sensitive to maxRank",
    subtitle = "Dashed line: median of 5,712 detected genes per spot"
  ) +
  base_theme
save_figure(p_maxrank, "figure_strength_maxrank", 6.6, 6.2)

# Figure 2: stability after amplifying the three largest SCT contributors.
amplification <- utils::read.csv(
  file.path(results_dir, "strength_dominant_gene_amplification.csv"),
  stringsAsFactors = FALSE
)
amplification$method <- unname(method_labels[amplification$method])
amplification_long <- rbind(
  data.frame(
    signature = amplification$signature,
    method = amplification$method,
    factor = amplification$amplification_factor,
    gene = amplification$amplified_gene,
    metric = "Spearman correlation",
    value = amplification$spearman
  ),
  data.frame(
    signature = amplification$signature,
    method = amplification$method,
    factor = amplification$amplification_factor,
    gene = amplification$amplified_gene,
    metric = "Top-10% hotspot Jaccard",
    value = amplification$top10_jaccard
  )
)
amplification_groups <- split(
  amplification_long,
  interaction(
    amplification_long$signature,
    amplification_long$method,
    amplification_long$factor,
    amplification_long$metric,
    drop = TRUE
  )
)
amplification_summary <- do.call(rbind, lapply(amplification_groups, function(d) {
  data.frame(
    signature = d$signature[[1]],
    method = d$method[[1]],
    factor = d$factor[[1]],
    metric = d$metric[[1]],
    median = stats::median(d$value),
    minimum = min(d$value),
    maximum = max(d$value)
  )
}))
rownames(amplification_summary) <- NULL
p_amplification <- ggplot2::ggplot(
  amplification_summary,
  ggplot2::aes(x = factor, y = median, color = method, group = method)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = minimum, ymax = maximum, fill = method),
    alpha = 0.12,
    color = NA
  ) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2) +
  ggplot2::facet_grid(metric ~ signature) +
  ggplot2::scale_x_continuous(breaks = c(2, 5, 10)) +
  ggplot2::scale_y_continuous(limits = c(0.4, 1), breaks = seq(0.4, 1, 0.1)) +
  ggplot2::scale_color_manual(values = method_colors) +
  ggplot2::scale_fill_manual(values = method_colors) +
  ggplot2::labs(
    x = "Expression multiplier applied to a dominant gene",
    y = "Agreement with unperturbed score",
    color = "Method",
    fill = "Method",
    title = "UCell resists domination by highly amplified ECM genes",
    subtitle = "Lines show the median and bands the range across the three largest SCT contributors"
  ) +
  base_theme
save_figure(p_amplification, "figure_strength_dominant_gene", 8.2, 6.0)

# Figure 3: stability of fixed evaluation spots across dataset compositions.
composition <- utils::read.csv(
  file.path(results_dir, "strength_composition_summary.csv"),
  stringsAsFactors = FALSE
)
composition$method <- unname(method_labels[composition$method])
context_labels <- c(
  evaluation_only = "Evaluation only",
  evaluation_plus_malignant = "+ malignant",
  evaluation_plus_stromal = "+ stromal",
  evaluation_plus_immune = "+ immune"
)
composition$context <- factor(
  unname(context_labels[composition$tested_context]),
  levels = unname(context_labels)
)
composition_long <- rbind(
  data.frame(
    signature = composition$signature,
    method = composition$method,
    context = composition$context,
    metric = "Spearman correlation",
    value = composition$spearman
  ),
  data.frame(
    signature = composition$signature,
    method = composition$method,
    context = composition$context,
    metric = "Top-10% hotspot Jaccard",
    value = composition$top10_jaccard
  )
)
p_composition <- ggplot2::ggplot(
  composition_long,
  ggplot2::aes(x = context, y = value, color = method, group = method, shape = method)
) +
  ggplot2::geom_hline(yintercept = 1, linewidth = 0.35, color = "grey70") +
  ggplot2::geom_line(linewidth = 0.75, position = ggplot2::position_dodge(width = 0.12)) +
  ggplot2::geom_point(size = 2.1, position = ggplot2::position_dodge(width = 0.12)) +
  ggplot2::facet_grid(metric ~ signature) +
  ggplot2::scale_y_continuous(limits = c(0.80, 1.002), breaks = seq(0.8, 1, 0.05)) +
  ggplot2::scale_color_manual(values = method_colors) +
  ggplot2::labs(
    x = "Dataset supplied while rescoring the same 500 spots",
    y = "Agreement with full-section score",
    color = "Method",
    shape = "Method",
    title = "UCell and per-spot aggregates are composition-independent"
  ) +
  base_theme +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
save_figure(p_composition, "figure_strength_composition", 8.4, 6.0)

# Figure 4: spot-rank error under heterogeneous capture efficiency.
heterogeneous <- utils::read.csv(
  file.path(results_dir, "strength_heterogeneous_capture_depth_summary.csv"),
  stringsAsFactors = FALSE
)
heterogeneous$method <- unname(method_labels[heterogeneous$method])
heterogeneous_long <- rbind(
  data.frame(
    signature = heterogeneous$signature,
    method = heterogeneous$method,
    retention = heterogeneous$retention,
    metric = "Absolute percentile-rank error",
    median = heterogeneous$median_abs_percentile_rank_error_median,
    q025 = heterogeneous$median_abs_percentile_rank_error_q025,
    q975 = heterogeneous$median_abs_percentile_rank_error_q975
  ),
  data.frame(
    signature = heterogeneous$signature,
    method = heterogeneous$method,
    retention = heterogeneous$retention,
    metric = "Percentile-rank shift",
    median = heterogeneous$median_percentile_rank_shift_median,
    q025 = heterogeneous$median_percentile_rank_shift_q025,
    q975 = heterogeneous$median_percentile_rank_shift_q975
  )
)
p_heterogeneous <- ggplot2::ggplot(
  heterogeneous_long,
  ggplot2::aes(x = retention, y = median, color = method, group = method)
) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.35, color = "grey65") +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = q025, ymax = q975, fill = method),
    alpha = 0.12,
    color = NA
  ) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2) +
  ggplot2::facet_grid(metric ~ signature, scales = "free_y") +
  ggplot2::scale_x_continuous(
    breaks = c(0.25, 0.50, 0.75, 1),
    labels = c("25%", "50%", "75%", "100%")
  ) +
  ggplot2::scale_color_manual(values = method_colors) +
  ggplot2::scale_fill_manual(values = method_colors) +
  ggplot2::labs(
    x = "UMIs retained independently for each spot",
    y = "Change from full-depth spot ranking",
    color = "Method",
    fill = "Method",
    title = "Normalized aggregation is most stable under heterogeneous capture",
    subtitle = "Lines show medians and bands the 2.5th to 97.5th percentiles across 20 simulations"
  ) +
  base_theme
save_figure(p_heterogeneous, "figure_strength_heterogeneous_capture", 8.2, 6.0)

# Reviewer figure 1: direct rank comparison with the requested summed-expression
# baseline. Normalized means and sums have identical ranks because every spot
# uses the same fixed signature genes.
spot_scores <- utils::read.csv(
  gzfile(file.path(results_dir, "spot_scores.csv.gz")),
  stringsAsFactors = FALSE
)
method_comparison <- utils::read.csv(
  file.path(results_dir, "method_comparison.csv"),
  stringsAsFactors = FALSE
)
percentile_rank <- function(x) {
  100 * (rank(x, ties.method = "average") - 1) / (length(x) - 1)
}
reviewer_scatter_rows <- lapply(split(spot_scores, spot_scores$signature), function(d) {
  rbind(
    data.frame(
      signature = d$signature,
      comparator = "Summed normalized expression",
      comparator_percentile = percentile_rank(d$sct_sum),
      ucell_percentile = percentile_rank(d$ucell_stapp)
    ),
    data.frame(
      signature = d$signature,
      comparator = "Raw UMI sum",
      comparator_percentile = percentile_rank(d$raw_umi_sum),
      ucell_percentile = percentile_rank(d$ucell_stapp)
    )
  )
})
reviewer_scatter <- do.call(rbind, reviewer_scatter_rows)
reviewer_scatter$comparator <- factor(
  reviewer_scatter$comparator,
  levels = c("Summed normalized expression", "Raw UMI sum")
)
reviewer_annotations <- do.call(rbind, lapply(seq_len(nrow(method_comparison)), function(i) {
  d <- method_comparison[i, ]
  rbind(
    data.frame(
      signature = d$signature,
      comparator = "Summed normalized expression",
      label = sprintf(
        "Spearman rho = %.3f\nTop-10%% Jaccard = %.3f",
        d$spearman_ucell_vs_sct_sum,
        d$top10_jaccard_ucell_vs_sct_sum
      )
    ),
    data.frame(
      signature = d$signature,
      comparator = "Raw UMI sum",
      label = sprintf(
        "Spearman rho = %.3f\nTop-10%% Jaccard = %.3f",
        d$spearman_ucell_vs_raw_umi_sum,
        d$top10_jaccard_ucell_vs_raw_umi_sum
      )
    )
  )
}))
reviewer_annotations$comparator <- factor(
  reviewer_annotations$comparator,
  levels = levels(reviewer_scatter$comparator)
)
p_reviewer_scatter <- ggplot2::ggplot(
  reviewer_scatter,
  ggplot2::aes(x = comparator_percentile, y = ucell_percentile)
) +
  ggplot2::geom_abline(
    intercept = 0,
    slope = 1,
    color = "grey55",
    linewidth = 0.55
  ) +
  ggplot2::geom_point(alpha = 0.20, size = 0.55, color = "#28536B") +
  ggplot2::geom_label(
    data = reviewer_annotations,
    ggplot2::aes(x = 3, y = 97, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3,
    label.size = 0.25,
    fill = "white"
  ) +
  ggplot2::facet_grid(signature ~ comparator) +
  ggplot2::coord_equal(xlim = c(0, 100), ylim = c(0, 100)) +
  ggplot2::scale_x_continuous(breaks = seq(0, 100, 25)) +
  ggplot2::scale_y_continuous(breaks = seq(0, 100, 25)) +
  ggplot2::labs(
    x = "Aggregate-expression percentile among spots",
    y = "UCell percentile among spots",
    title = "UCell agrees with normalized ECM summation but not raw summation",
    subtitle = "Each point is one of 4,992 spots; the diagonal indicates identical ranking"
  ) +
  base_theme +
  ggplot2::theme(legend.position = "none")
save_figure(p_reviewer_scatter, "figure_reviewer_ucell_vs_sum", 8.4, 6.4)

# Reviewer figure 2: spatial comparison on the same percentile scale.
obj <- readRDS(file.path(
  script_dir,
  "data",
  "Breast_cancer_Ductal_Carcinoma_In_Situ_Invasive.RDS"
))
coordinates <- tryCatch(
  Seurat::GetTissueCoordinates(obj),
  error = function(e) NULL
)
if (is.null(coordinates) && length(obj@images) > 0L) {
  centroids <- tryCatch(
    attributes(obj@images[[1]])$boundaries$centroids,
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
if (is.null(coordinates)) stop("Could not recover spatial coordinates.")
if (all(c("imagecol", "imagerow") %in% colnames(coordinates))) {
  coordinate_table <- data.frame(
    spot_id = rownames(coordinates),
    x = coordinates$imagecol,
    y = -coordinates$imagerow
  )
} else {
  numeric_columns <- which(vapply(coordinates, is.numeric, logical(1)))
  coordinate_table <- data.frame(
    spot_id = rownames(coordinates),
    x = coordinates[[numeric_columns[[1]]]],
    y = -coordinates[[numeric_columns[[2]]]]
  )
}
spatial_rows <- lapply(split(spot_scores, spot_scores$signature), function(d) {
  rbind(
    data.frame(
      spot_id = d$spot_id,
      signature = d$signature,
      method = "UCell",
      percentile = percentile_rank(d$ucell_stapp)
    ),
    data.frame(
      spot_id = d$spot_id,
      signature = d$signature,
      method = "Summed normalized expression",
      percentile = percentile_rank(d$sct_sum)
    ),
    data.frame(
      spot_id = d$spot_id,
      signature = d$signature,
      method = "Raw UMI sum",
      percentile = percentile_rank(d$raw_umi_sum)
    )
  )
})
spatial_data <- merge(
  do.call(rbind, spatial_rows),
  coordinate_table,
  by = "spot_id",
  all.x = TRUE,
  sort = FALSE
)
spatial_data$method <- factor(
  spatial_data$method,
  levels = c("UCell", "Summed normalized expression", "Raw UMI sum")
)
p_reviewer_spatial <- ggplot2::ggplot(
  spatial_data,
  ggplot2::aes(x = x, y = y, color = percentile)
) +
  ggplot2::geom_point(size = 0.60) +
  ggplot2::facet_grid(signature ~ method) +
  ggplot2::coord_equal() +
  ggplot2::scale_color_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 50,
    limits = c(0, 100),
    name = "Spot\npercentile"
  ) +
  ggplot2::labs(
    title = "Spatial comparison of UCell and summed ECM expression",
    subtitle = "Scores are displayed as within-method percentiles for direct visual comparison"
  ) +
  ggplot2::theme_void(base_size = 10) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(color = "black"),
    strip.background = ggplot2::element_rect(fill = "grey94", color = "grey50"),
    plot.title = ggplot2::element_text(hjust = 0),
    plot.subtitle = ggplot2::element_text(hjust = 0),
    legend.position = "right"
  )
save_figure(p_reviewer_spatial, "figure_reviewer_spatial_percentiles", 9.4, 6.2)

message("Strength-test figures written to: ", results_dir)
