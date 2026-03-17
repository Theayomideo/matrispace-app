#' Initialize all reference data for MatriSpace
#'
#' Loads matrisome reference data, standardizes gene symbols,
#' builds feature signatures, and defines ECM quotes.
#' Called once at server startup.
#'
#' @return A named list containing all reference data objects:
#'   matrisome, mobj, recs, lr_db, ecm_ucell_signatures,
#'   matrisome_feature_signatures, ecm_quotes
#' @noRd
initialize_reference_data <- function() {

  # Load curated ECM niche signatures (not derivable from matrisome_v2)
  ecm_ucell_signatures <- readRDS(extdata_path("ecm_ucell_signatures.rds"))

  # Load LR database for spatial co-expression analysis
  lr_db <- readRDS(extdata_path("ultimate_ecm_interactions_DEDUPLICATED.rds"))

  # Matrisome reference (single source of truth for all matrisome data)
  matrisome <- readRDS(extdata_path("matrisome_v2.rds"))
  recs <- unique(readRDS(extdata_path("receptors.RDS"))$to)

  # Standardize gene symbols in reference data to current HGNC nomenclature
  # Helper: standardize a vector that may contain duplicates
  .standardize_genes <- function(genes) {
    uniq <- unique(genes)
    result <- scCustomize::Updated_HGNC_Symbols(uniq, verbose = FALSE, case_check_as_warn = TRUE)
    mapping <- setNames(result$Output_Features, result$input_features)
    unname(mapping[genes])
  }

  suppressMessages(suppressWarnings({
    # Standardize matrisome gene symbols
    matrisome$gene <- .standardize_genes(matrisome$gene)

    # Derive mobj from matrisome (single source of truth)
    mobj <- unique(matrisome$gene)

    # Standardize other reference data
    recs <- unique(.standardize_genes(recs))

    # LR database columns - then drop duplicate interaction pairs
    lr_db$Ligand <- .standardize_genes(lr_db$Ligand)
    lr_db$Receptor <- .standardize_genes(lr_db$Receptor)
    lr_db <- lr_db[!duplicated(lr_db[, c("Ligand", "Receptor")]), ]

    # ECM niche signatures - unique within each signature
    ecm_ucell_signatures <- lapply(ecm_ucell_signatures, function(g) unique(.standardize_genes(g)))
  }))

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
    lr_db = lr_db,
    ecm_ucell_signatures = ecm_ucell_signatures,
    matrisome_feature_signatures = matrisome_feature_signatures,
    ecm_quotes = ecm_quotes
  )
}
