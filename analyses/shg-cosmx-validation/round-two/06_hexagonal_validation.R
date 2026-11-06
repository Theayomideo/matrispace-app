options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
analysis_root <- if (length(args) >= 1L) args[[1L]] else "orthogonal_validation"
derived_dir <- file.path(analysis_root, "derived")
results_dir <- file.path(analysis_root, "results")
figure_dir <- file.path(analysis_root, "figures", "individual_svgs")
qa_dir <- file.path(figure_dir, "qa")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

required <- c("data.table", "ggplot2", "hexbin", "svglite", "cowplot", "ragg")
missing <- required[!vapply(required, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

library(data.table)
library(ggplot2)
library(hexbin)

cells <- fread(file.path(derived_dir, "section04_cell_collagen_scores.csv.gz"))
image <- readRDS(file.path(derived_dir, "section03_shg_collagen_image.rds"))
storage.mode(image) <- "double"
image <- image / max(image, na.rm = TRUE)

field_width_um <- 4000
pixel_size_um <- field_width_um / ncol(image)
xbins <- 37L
min_cells <- 10L
block_size_um <- 400
n_boot <- 1000L
n_perm <- 999L
seed <- 15240431L

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || sd(x[keep]) == 0 || sd(y[keep]) == 0) return(NA_real_)
  suppressWarnings(cor(x[keep], y[keep], method = "spearman"))
}

jaccard_top <- function(x, y, fraction = 0.20) {
  keep <- is.finite(x) & is.finite(y)
  x_cut <- quantile(x[keep], 1 - fraction, names = FALSE, type = 8)
  y_cut <- quantile(y[keep], 1 - fraction, names = FALSE, type = 8)
  xi <- which(keep & x >= x_cut)
  yi <- which(keep & y >= y_cut)
  length(intersect(xi, yi)) / length(union(xi, yi))
}

partial_rank_correlation <- function(x, y, covariates) {
  keep <- is.finite(x) & is.finite(y) &
    apply(covariates, 1L, function(z) all(is.finite(z)))
  z <- as.data.frame(covariates[keep, , drop = FALSE])
  z[] <- lapply(z, rank, ties.method = "average")
  rx <- resid(lm(rank(x[keep], ties.method = "average") ~ ., data = z))
  ry <- resid(lm(rank(y[keep], ties.method = "average") ~ ., data = z))
  cor(rx, ry)
}

message("Assigning all CosMx cells to area-matched hexagons...")
hb_cells <- hexbin(
  cells$x_um, cells$y_um,
  xbins = xbins,
  xbnds = c(0, field_width_um),
  ybnds = c(0, field_width_um),
  shape = 1,
  IDs = TRUE
)

cell_summary <- data.table(
  hex_id = hb_cells@cID,
  collagens_ucell = cells$collagens_ucell,
  fibrillar_collagens_ucell = cells$fibrillar_collagens_ucell,
  nCount_Spatial = cells$nCount_Spatial,
  nFeature_Spatial = cells$nFeature_Spatial
)[, .(
  n_cells = .N,
  collagens_ucell = mean(collagens_ucell, na.rm = TRUE),
  fibrillar_collagens_ucell = mean(fibrillar_collagens_ucell, na.rm = TRUE),
  nCount_Spatial = mean(nCount_Spatial, na.rm = TRUE),
  nFeature_Spatial = mean(nFeature_Spatial, na.rm = TRUE)
), by = hex_id]

centers_xy <- hcell2xy(hb_cells)
centers <- data.table(
  hex_id = hb_cells@cell,
  x_um = centers_xy$x,
  y_um = centers_xy$y
)

message("Assigning the full SHG raster to the same hexagons...")
x_coords <- (seq_len(ncol(image)) - 0.5) * pixel_size_um
y_coords <- (seq_len(nrow(image)) - 0.5) * pixel_size_um
shg_points <- data.table(
  x_um = rep(x_coords, each = nrow(image)),
  y_um = rep(y_coords, times = ncol(image)),
  intensity = as.vector(image)
)

hb_shg <- hexbin(
  shg_points$x_um, shg_points$y_um,
  xbins = xbins,
  xbnds = c(0, field_width_um),
  ybnds = c(0, field_width_um),
  shape = 1,
  IDs = TRUE
)

shg_summary <- data.table(
  hex_id = hb_shg@cID,
  intensity = shg_points$intensity
)[, .(
  shg_pixels = .N,
  shg_mean = mean(intensity, na.rm = TRUE),
  shg_median = median(intensity, na.rm = TRUE)
), by = hex_id]
rm(shg_points)
invisible(gc())

hexes <- merge(cell_summary, shg_summary, by = "hex_id")
hexes <- merge(hexes, centers, by = "hex_id")
hexes <- hexes[n_cells >= min_cells]

horizontal_width_um <- field_width_um / xbins
hex_side_um <- horizontal_width_um / sqrt(3)
hex_height_um <- 2 * hex_side_um
row_spacing_um <- 1.5 * hex_side_um
hex_area_um2 <- sqrt(3) / 2 * horizontal_width_um^2
equivalent_square_side_um <- sqrt(hex_area_um2)

# Offset-grid and axial coordinates for topology-preserving toroidal shifts.
hexes[, hex_row := as.integer(round(y_um / row_spacing_um))]
hexes[, hex_col := as.integer(round(
  (x_um - (hex_row %% 2L) * horizontal_width_um / 2) / horizontal_width_um
))]
hexes[, axial_q := hex_col - ((hex_row - (hex_row %% 2L)) %/% 2L)]
hexes[, axial_q_mod := axial_q %% hb_cells@dimen[[2L]]]

hexes[, cosmx_percentile := (frank(collagens_ucell, ties.method = "average") - 0.5) / .N]
hexes[, shg_percentile := (frank(shg_mean, ties.method = "average") - 0.5) / .N]
hexes[, joint_collagen_index :=
  2 * cosmx_percentile * shg_percentile /
  pmax(cosmx_percentile + shg_percentile, 1e-8)]

observed <- safe_spearman(hexes$collagens_ucell, hexes$shg_mean)

block_bootstrap_ci <- function(d) {
  z <- copy(d)
  z[, block_x := floor(x_um / block_size_um)]
  z[, block_y := floor(y_um / block_size_um)]
  z[, block_id := paste(block_x, block_y, sep = "_")]
  block_ids <- unique(z$block_id)
  rows_by_block <- split(seq_len(nrow(z)), z$block_id)
  set.seed(seed + 101L)
  boot <- vapply(seq_len(n_boot), function(i) {
    sampled <- sample(block_ids, length(block_ids), replace = TRUE)
    index <- unlist(rows_by_block[sampled], use.names = FALSE)
    safe_spearman(z$collagens_ucell[index], z$shg_mean[index])
  }, numeric(1L))
  list(
    ci = quantile(boot, c(0.025, 0.975), na.rm = TRUE, names = FALSE, type = 8),
    values = boot
  )
}

hex_circular_shift <- function(d) {
  ny <- hb_cells@dimen[[1L]]
  nx <- hb_cells@dimen[[2L]]
  score_matrix <- matrix(NA_real_, nrow = ny, ncol = nx)
  shg_matrix <- matrix(NA_real_, nrow = ny, ncol = nx)
  index <- cbind(d$hex_row + 1L, d$axial_q_mod + 1L)
  score_matrix[index] <- d$collagens_ucell
  shg_matrix[index] <- d$shg_mean

  matrix_observed <- safe_spearman(as.vector(score_matrix), as.vector(shg_matrix))
  if (!isTRUE(all.equal(observed, matrix_observed, tolerance = 1e-12))) {
    stop("Hex coordinate conversion changed the observed correlation")
  }

  set.seed(seed + 202L)
  shifts <- data.table(
    dr = sample.int(ny, n_perm, replace = TRUE) - 1L,
    dq = sample.int(nx, n_perm, replace = TRUE) - 1L
  )
  zero <- shifts$dr == 0L & shifts$dq == 0L
  while (any(zero)) {
    shifts$dr[zero] <- sample.int(ny, sum(zero), replace = TRUE) - 1L
    shifts$dq[zero] <- sample.int(nx, sum(zero), replace = TRUE) - 1L
    zero <- shifts$dr == 0L & shifts$dq == 0L
  }

  null <- vapply(seq_len(n_perm), function(i) {
    row_order <- ((seq_len(ny) - 1L + shifts$dr[[i]]) %% ny) + 1L
    col_order <- ((seq_len(nx) - 1L + shifts$dq[[i]]) %% nx) + 1L
    shifted <- shg_matrix[row_order, col_order, drop = FALSE]
    safe_spearman(as.vector(score_matrix), as.vector(shifted))
  }, numeric(1L))

  p <- (1 + sum(abs(null) >= abs(matrix_observed), na.rm = TRUE)) /
    (1 + sum(is.finite(null)))
  list(null = null, p = p, shifts = shifts)
}

bootstrap <- block_bootstrap_ci(hexes)
shift <- hex_circular_shift(hexes)

inference <- data.table(
  spatial_unit = "area-matched regular hexagon",
  xbins = xbins,
  horizontal_width_um = horizontal_width_um,
  vertical_height_um = hex_height_um,
  hex_area_um2 = hex_area_um2,
  equivalent_square_side_um = equivalent_square_side_um,
  minimum_cells = min_cells,
  retained_hexagons = nrow(hexes),
  spearman = observed,
  block_bootstrap_ci_low = bootstrap$ci[[1L]],
  block_bootstrap_ci_high = bootstrap$ci[[2L]],
  circular_shift_p = shift$p,
  top20_jaccard = jaccard_top(hexes$collagens_ucell, hexes$shg_mean),
  partial_spearman_adjusted_cell_density_and_nFeature = partial_rank_correlation(
    hexes$collagens_ucell,
    hexes$shg_mean,
    as.matrix(hexes[, .(n_cells, nFeature_Spatial)])
  )
)

fwrite(hexes, file.path(results_dir, "primary_hex_data_area_matched.csv.gz"))
fwrite(inference, file.path(results_dir, "primary_hex_inference.csv"))
fwrite(
  data.table(null_spearman = shift$null, dr = shift$shifts$dr, dq = shift$shifts$dq),
  file.path(results_dir, "hex_spatial_shift_null.csv.gz")
)
fwrite(
  data.table(bootstrap_spearman = bootstrap$values),
  file.path(results_dir, "hex_block_bootstrap.csv.gz")
)

metadata <- data.table(
  field = c(
    "timestamp", "field_width_um", "image_pixels", "pixel_size_um",
    "shg_normalization", "hex_xbins", "hex_horizontal_width_um",
    "hex_vertical_height_um", "hex_area_um2", "equivalent_square_side_um",
    "minimum_cells", "spatial_null", "spatial_null_permutations",
    "bootstrap", "bootstrap_repetitions", "bootstrap_block_size_um"
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), field_width_um,
    paste(dim(image), collapse = "x"), pixel_size_um,
    "divide by image maximum", xbins, horizontal_width_um,
    hex_height_um, hex_area_um2, equivalent_square_side_um,
    min_cells, "toroidal circular shifts on axial hex coordinates", n_perm,
    "spatial block bootstrap", n_boot, block_size_um
  )
)
fwrite(metadata, file.path(results_dir, "hex_run_metadata.csv"))

# Plot-only geometry.
angles <- pi / 6 + (0:5) * pi / 3
offsets <- data.table(
  vertex = 1:6,
  dx = hex_side_um * cos(angles),
  dy = hex_side_um * sin(angles)
)
polygons <- hexes[, .(
  vertex = offsets$vertex,
  x_vertex = x_um + offsets$dx,
  y_vertex = y_um + offsets$dy,
  joint_collagen_index = joint_collagen_index
), by = hex_id]

collagen_palette <- c(
  "#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4", "#E6F598",
  "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"
)

map_theme <- theme_void(base_family = "sans") +
  theme(
    plot.margin = margin(0, 0, 0, 0),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.position = "none"
  )

p_hex <- ggplot(polygons, aes(x_vertex, y_vertex, group = hex_id,
                              fill = joint_collagen_index)) +
  geom_polygon(colour = NA) +
  scale_fill_gradientn(colours = collagen_palette, limits = c(0, 1), na.value = "white") +
  annotate(
    "segment", x = 180, xend = 680, y = 3760, yend = 3760,
    colour = "black", linewidth = 1.35, lineend = "butt"
  ) +
  scale_y_reverse(limits = c(field_width_um, 0), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0, field_width_um), expand = c(0, 0)) +
  coord_fixed(clip = "on") +
  map_theme

p_scatter <- ggplot(hexes, aes(collagens_ucell, shg_mean)) +
  geom_point(alpha = 0.32, size = 0.85, colour = "#355C7D") +
  geom_smooth(
    method = "lm", formula = y ~ x, se = TRUE,
    colour = "#B53636", fill = "#E9A1A1", linewidth = 0.8
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 10.5, base_family = "sans") +
  theme(
    plot.margin = margin(5, 7, 5, 5),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8.5, colour = "black"),
    axis.line = element_line(linewidth = 0.45),
    axis.ticks = element_line(linewidth = 0.4),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA)
  )

save_svg <- function(plot, filename, width, height) {
  svglite::svglite(file.path(figure_dir, filename), width = width, height = height,
                   bg = "transparent")
  print(plot)
  dev.off()
}

save_png <- function(plot, filename, width, height) {
  ragg::agg_png(file.path(qa_dir, filename), width = width, height = height,
                units = "in", res = 180, background = "white")
  print(plot)
  dev.off()
}

save_svg(p_hex, "03_joint_collagen_index_hex_map.svg", 4.6, 4.6)
save_svg(p_scatter, "04_hex_concordance_scatter.svg", 5.2, 4.3)
save_png(p_hex, "qa_03_joint_collagen_index_hex_map.png", 4.6, 4.6)
save_png(p_scatter, "qa_04_hex_concordance_scatter.png", 5.2, 4.3)

legend_plot <- ggplot(data.frame(x = c(0, 1), y = c(0, 1)), aes(x, y, colour = y)) +
  geom_point() +
  scale_colour_gradientn(
    colours = collagen_palette, limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = scales::label_number(accuracy = 0.01),
    name = NULL,
    guide = guide_colourbar(
      barheight = grid::unit(42, "mm"),
      barwidth = grid::unit(6.5, "mm"),
      ticks = TRUE, frame.colour = "black", frame.linewidth = 0.45
    )
  ) +
  theme_void(base_family = "sans") +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 9, colour = "black"),
    legend.key.height = grid::unit(42, "mm"),
    legend.key.width = grid::unit(6.5, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(0, 0, 0, 0),
    plot.background = element_rect(fill = "transparent", colour = NA)
  )
legend <- cowplot::get_legend(legend_plot)

svglite::svglite(
  file.path(figure_dir, "legend_03_joint_collagen_index_hex.svg"),
  width = 0.90, height = 1.90, bg = "transparent"
)
grid::grid.newpage()
grid::grid.draw(legend)
dev.off()

ragg::agg_png(
  file.path(qa_dir, "qa_legend_03_joint_collagen_index_hex.png"),
  width = 0.90, height = 1.90, units = "in", res = 220, background = "white"
)
grid::grid.newpage()
grid::grid.draw(legend)
dev.off()

cat("Area-matched hexagonal validation complete\n")
print(inference)
