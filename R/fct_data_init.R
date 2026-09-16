#' Initialize all reference data for MatriSpace
#'
#' Loads matrisome reference data, standardizes gene symbols,
#' builds feature signatures, and defines ECM quotes.
#' Called once at server startup.
#'
#' @return A named list containing all reference data objects:
#'   matrisome, mobj, recs, matrisome_pairs, ecm_ucell_signatures,
#'   matrisome_feature_signatures, ecm_quotes
#' @noRd
initialize_reference_data <- function() {

  # Load curated ECM niche signatures (not derivable from matrisome_v2)
  ecm_ucell_signatures <- readRDS(extdata_path("ecm_ucell_signatures.rds"))

  # Load the MatriComDB pair database for spatial co-expression analysis
  matrisome_pairs <- readRDS(extdata_path("matrisome_pairs.rds"))

  # MatriComDB contains multiple annotated interaction classes, including
  # structural heteromeric interactions and homomeric assemblies. MatriSpace
  # scores pairwise spatial co-expression rather than directionality, so
  # exclude homomeric rows and collapse reciprocal rows to one unique pair.
  matrisome_pairs <- matrisome_pairs %>%
    filter(
      !is.na(Gene1), !is.na(Gene2),
      nzchar(Gene1), nzchar(Gene2),
      Gene1 != Gene2
    ) %>%
    mutate(
      .pair_gene_1 = pmin(Gene1, Gene2),
      .pair_gene_2 = pmax(Gene1, Gene2)
    ) %>%
    group_by(.pair_gene_1, .pair_gene_2) %>%
    summarise(
      Gene1 = dplyr::first(.pair_gene_1),
      Gene2 = dplyr::first(.pair_gene_2),
      Interaction_Type = paste(sort(unique(Interaction_Type)), collapse = "; "),
      Source = paste(sort(unique(Source)), collapse = "; "),
      Database_Score = max(Database_Score),
      .groups = "drop"
    ) %>%
    select(Gene1, Gene2, Interaction_Type, Source, Database_Score)

  # Matrisome reference (single source of truth for all matrisome data)
  matrisome <- readRDS(extdata_path("matrisome_v2.rds"))
  recs <- unique(readRDS(extdata_path("receptors.RDS"))$to)
  mobj <- unique(matrisome$gene)

  # Build matrisome_feature_signatures from matrisome data frame
  # Main categories: filter by notes column
  matrisome_feature_signatures <- list(
    ecm_glycoproteins = matrisome$gene[matrisome$notes == "ECM Glycoproteins"],
    collagens = matrisome$gene[matrisome$notes == "Collagens"],
    proteoglycans = matrisome$gene[matrisome$notes == "Proteoglycans"],
    ecm_regulators = matrisome$gene[matrisome$notes == "ECM Regulators"],
    secreted_factors = matrisome$gene[matrisome$notes == "Secreted Factors"],
    `ecm-affiliated_proteins` = matrisome$gene[matrisome$notes == "ECM-affiliated Proteins"]
  )

  # Subcategories + families: filter by ecm_subcategory column (semicolon-separated)
  .subcategory_map <- c(
    "basement_membrane" = "Basement Membrane",
    "hemostasis" = "Hemostasis",
    "elastic_fibers" = "Elastic fibers",
    "growth_factor-binding" = "Growth Factor-binding",
    "laminin_-_basement_membrane" = "Laminins",
    "matricellular" = "Matricellular proteins",
    "syndecan" = "Syndecan",
    "glypican" = "Glypican",
    "annexin" = "Annexin",
    "cathepsin" = "Cathepsin",
    "ccn_family" = "CCN Family",
    "cystatin" = "Cystatin",
    "facit" = "FACIT",
    "fibulin" = "Fibulin",
    "galectin" = "Galectin",
    "mucin" = "Mucin",
    "plexin" = "Plexin",
    "semaphorin" = "Semaphorin",
    "perivascular" = "Peri-vascular ECM"
  )

  for (.sig_name in names(.subcategory_map)) {
    matrisome_feature_signatures[[.sig_name]] <- matrisome$gene[
      grepl(.subcategory_map[[.sig_name]], matrisome$ecm_subcategory, fixed = TRUE)
    ]
  }

  # Remove NAs (genes without ecm_subcategory)
  matrisome_feature_signatures <- lapply(matrisome_feature_signatures, function(g) g[!is.na(g)])

  # Static data for UI elements
  ecm_quotes <- c(
    "More than an inert scaffold, the ECM constitutes a dynamic repository of biological information, actively directing cell fate.",
    "The ECM is a social network, storing and presenting signals essential to instruct cell phenotypes.",
    "The physical and mechanical properties of the ECM, can be converted into potent and instructive biochemical signals.",
    "Analyzing the ECM of microenvironments is critical to understand cellular context.",
    "Each tissue possesses a unique matrisome\u2014a specific repertoire of ECM proteins that dictates cellular identity and function.",
    "Dysregulation of the ECM is not only a consequence of disease, but a central driver of pathologies like fibrosis and cancer.",
    "Cells physically interact with the ECM via cell-surface receptors that can translate extracellular topography into intracellular commands.",
    "The ECM is a long-term storage of the history of a tissue. Its composition is the sum of events that have occurred during tissue development, injury, and repair.",
    "In the tumor context, the ECM acts first as a barrier to be breached by cells for effective dissemination, but then as a permissive highway.",
    "The molecular diversity of the ECM, from the collagen backbone to proteoglycan cushions, enables the emergence of complex cellular assemblies and function.",
    "Once viewed as the 'glue holding cells together', we now recognize the ECM as a master regulator of tissue homeostasis.",
    "MatriSpace: Unraveling the spatial code of the matrisome, one transcript at a time."
  )

  # Return all data objects as a named list
  list(
    matrisome = matrisome,
    mobj = mobj,
    recs = recs,
    matrisome_pairs = matrisome_pairs,
    ecm_ucell_signatures = ecm_ucell_signatures,
    matrisome_feature_signatures = matrisome_feature_signatures,
    ecm_quotes = ecm_quotes
  )
}
