args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: Rscript run_leakage_ranking.R input.rds [output_dir]")
input <- normalizePath(args[1], mustWork = TRUE)
out <- if (length(args) > 1L) args[2] else "results"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
app <- asNamespace("matrispace.app")
ref <- app$initialize_reference_data()
obj <- app$standardize_gene_symbols(readRDS(input))
source <- app$default_expression_source(obj)
counts <- app$safe_get_assay_data(obj, assay = source$assay, slot = "counts")
corrected <- app$matrispace_spotclean(obj, genes = ref$matrisome$gene,
                                    alpha = 0.10, k = 6L, assays = source$assay)
if (!corrected$applied) stop(corrected$message)
adjusted <- app$safe_get_assay_data(corrected$object, assay = source$assay, slot = "counts")
stopifnot(identical(dimnames(counts), dimnames(adjusted)))
results <- do.call(rbind, lapply(names(ref$matrisome_feature_signatures)[1:6], function(category) {
  genes <- intersect(ref$matrisome_feature_signatures[[category]], rownames(counts))
  before <- Matrix::colSums(counts[genes, , drop = FALSE])
  after <- Matrix::colSums(adjusted[genes, , drop = FALSE])
  data.frame(category = category, genes = length(genes), spots = length(before),
             spearman = cor(before, after, method = "spearman"))
}))
write.csv(results, file.path(out, "canonical_leakage_ranking.csv"), row.names = FALSE)
write.csv(data.frame(input = basename(input), md5 = unname(tools::md5sum(input)),
                     assay = source$assay, layer = "counts", alpha = 0.10, k = 6L,
                     detected_canonical = sum(rownames(obj) %in% ref$mobj),
                     app_version = as.character(utils::packageVersion("matrispace.app"))),
          file.path(out, "run_metadata.csv"), row.names = FALSE)
writeLines(trimws(capture.output(sessionInfo()), which = "right"), file.path(out, "sessionInfo.txt"))
print(results)
