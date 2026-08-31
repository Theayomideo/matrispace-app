# Suppress R CMD check notes for dplyr/ggplot2 NSE variables
utils::globalVariables(c(
  ".axis_gene_1", ".axis_gene_2", "Database_Score", "Interaction_Type", "Source",
  "LISA_dynamic", "Ligand", "Receptor", "annotation", "annotation_group",
  "avg_log2FC", "clust", "cluster", "color", "domain", "feature",
  "feature1", "feature2", "imagecol", "imagerow", "mean_score",
  "perc_difference", "proportion", "scalefactors", "score", "signature_name"
))
