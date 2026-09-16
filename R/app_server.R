# App Server

app_server <- function(input, output, session) {

  # Central reactiveValues object for analysis results
  rv <- reactiveValues(
    analysis_results = NULL,
    matrisome_results = NULL,
    raw_uploaded_object = NULL,
    active_seurat_object = NULL,
    active_sample_name = NULL,
    loaded_collection = NULL,
    data_is_loaded = FALSE,
    feature_analysis_done = FALSE,
    matrisome_analysis_done = FALSE,
    matrisome_pair_results = NULL,
    is_visiumv2_origin = FALSE,
    custom_gene_sets = list(),
    primary_contrast_state = "none",
    secondary_contrast_state = "none"
  )


  # ---------------------------------------------------------------------------
  # INITIAL SETUP & OBSERVERS
  # ---------------------------------------------------------------------------
  # Helper function: Convert to sentence case (capitalize first letter only)
  to_sentence_case <- function(text) {
    if (is.na(text) || text == "") return(text)
    paste0(toupper(substring(text, 1, 1)), substring(text, 2))
  }

  # Load shinyjs
  useShinyjs()

  # Initialize reference data (replaces global.R loading)
  ref_data <- initialize_reference_data()
  matrisome <- ref_data$matrisome
  mobj <- ref_data$mobj
  recs <- ref_data$recs
  matrisome_pairs <- ref_data$matrisome_pairs
  ecm_ucell_signatures <- ref_data$ecm_ucell_signatures
  matrisome_feature_signatures <- ref_data$matrisome_feature_signatures
  ecm_quotes <- ref_data$ecm_quotes

  # --- Guided Tour Definition ---

  guide <- Cicerone$
    new()$
    step(
      el = "app-branding",
      title = "Welcome to MatriSpace",
      description = "An interactive platform for analysing spatial transcriptomics data. This guide will introduce the key steps of the workflow to link matrisome gene expression levels to spatial contexts."
    )$
    step(
      el = "seurat_file-label",
      title = "1. Upload Your Data",
      description = "Upload a pre-processed Seurat or SpatialExperiment object (.rds file) to begin your session. MatriSpace will automatically detect and process your data."
    )$
    step(
      el = "load_upload",
      title = "Load Your Data",
      description = "After selecting your file, click here to load the data. MatriSpace will check if preprocessing is needed and guide you through the process."
    )$
    step(
      el = "run_matrisome_profile",
      title = "Workflow 1: Matrisome Profiling",
      description = "Choose this discovery-oriented workflow to gain a holistic view of ECM gene expression. This analysis calculates scores for pre-defined matrisome profiles, revealing overall matrisome gene and gene set expression patterns and identifying spatial hotspots of matrisome gene/gene set expression."
    )$
    step(
      el = "primary-feature-selection",
      title = "Workflow 2: Feature Analysis",
      description = "Choose this hypothesis-driven workflow to investigate the expression of specific genes or gene sets. Select a primary feature and an optional secondary feature to perform co-expression, co-localisation, and spatial relationship analyses."
    )$
    step(
      el = "run_button",
      title = "Run Feature Analysis",
      description = "After selecting your feature(s) in Workflow 2, click this button. The application will generate a suite of interactive visualisations for spatially-resolved quantitative analyses."
    )$
    step(
      el = "export_dropdown_btn",
      title = "Export Results",
      description = "Once your analysis is complete, use this button to download all data tables and high-resolution plots (PNG and PDF) packaged in a single ZIP archive for your records or publications."
    )$
    step(
      el = "about",
      title = "Tour Complete",
      description = "You have now completed this guided tour. For further information on data sources and citations, please consult the 'About' section. You can restart this tour at any time from the sidebar."
    )

  # --- Tour Trigger ---
  observeEvent(input$guided_tour, {
    # Open all accordions for tour, then reset to initial state when tour ends
    shinyjs::runjs("
      window.matrispaceTourActive = true;

      // Reset UI to initial state (only Load Data open)
      window.resetToInitialState = function() {
        document.querySelectorAll('#main_accordion .accordion-button').forEach(function(btn) {
          if(btn.getAttribute('aria-expanded') === 'true') {
            btn.click();
          }
        });
        setTimeout(function() {
          const loadDataBtn = Array.from(document.querySelectorAll('#main_accordion .accordion-button'))
            .find(el => el.textContent.trim().includes('Load Data'));
          if(loadDataBtn && loadDataBtn.getAttribute('aria-expanded') === 'false') {
            loadDataBtn.click();
          }
        }, 50);
        window.matrispaceTourActive = false;
      };

      // MutationObserver to detect when tour overlay is removed
      if(window.matrispaceTourObserver) {
        window.matrispaceTourObserver.disconnect();
      }
      window.matrispaceTourObserver = new MutationObserver(function(mutations) {
        for(const mutation of mutations) {
          for(const node of mutation.removedNodes) {
            if(node.nodeType === 1 && node.id === 'driver-page-overlay') {
              if(window.matrispaceTourActive) {
                window.resetToInitialState();
                window.matrispaceTourObserver.disconnect();
              }
              return;
            }
          }
        }
      });
      window.matrispaceTourObserver.observe(document.body, { childList: true, subtree: true });

      // Open all accordion panels for the tour
      document.querySelectorAll('#main_accordion .accordion-button').forEach(function(btn) {
        if(btn.getAttribute('aria-expanded') === 'false') {
          btn.click();
        }
      });
    ")

    # Start the tour
    guide$init()$start()
  })

  # Accordion completion system using bslib functions
  success_icon <- bsicons::bs_icon("check-circle-fill", class="text-success", title = "Completed")
  original_load_icon <- bsicons::bs_icon("database-add")
  original_profile_matrisome_icon <- bsicons::bs_icon("bezier2")
  original_feature_icon <- bsicons::bs_icon("ui-checks-grid")

  # Observer to enable/disable main division point size slider
  observeEvent(input$toggle_main_divisions, {
    # toggleState is from shinyjs
    toggleState(id = "main_pt_size", condition = !input$toggle_main_divisions)
  })

  # Observer to enable/disable subdivision point size slider
  observeEvent(input$toggle_subdivisions, {
    toggleState(id = "sub_pt_size", condition = !input$toggle_subdivisions)
  })

  # Track completion state
  completed_sections <- reactiveValues(
    load_data = FALSE,
    profile_matrisome = FALSE,
    feature_selection = FALSE
  )

  # Centralized observer for successful data loading
  observeEvent(rv$data_is_loaded, {
    req(rv$data_is_loaded == TRUE)
    sendSweetAlert(session, title = "Success!", text = "Data loaded.", type = "success")

    # 2. Update the accordion to show completion and move to the next step.
    accordion_panel_update("main_accordion", "Load Data", icon = success_icon)
    completed_sections$load_data <- TRUE
    accordion_panel_set("main_accordion", "Profile Matrisome")
  }, ignoreInit = TRUE) # ignoreInit prevents this from firing on app startup

  # NOTE: rv$data_is_loaded is now set directly in the data loading observers
  # (No need for a separate observer since the universal container approach handles this)

  # Run Analysis button handler
  observeEvent(input$run_button, {
    # Update Select Feature panel with success icon
    accordion_panel_update("main_accordion", "Select Feature", icon = success_icon)
    completed_sections$feature_selection <- TRUE
    # Auto-open Select Feature panel to show results
    accordion_panel_set("main_accordion", "Select Feature")
    # ALSO navigate to the main content panel
    nav_select("main_content", "feature_selection")
    # Navigate to Feature Analysis tab within Feature Selection (where feature plots are)
    nav_select("feature_tabs", "Feature Analysis")
  })

  # Combined handler for green tick removal and navigation
  observeEvent(input$active_panel, {

    if (!is.null(input$active_panel)) {
      # FIRST: Handle green tick removal when reopening completed sections
      if (input$active_panel == "load_data" && completed_sections$load_data) {
        accordion_panel_update("main_accordion", "Load Data", icon = original_load_icon)
        completed_sections$load_data <- FALSE
      } else if (input$active_panel == "profile_matrisome" && completed_sections$profile_matrisome) {
        accordion_panel_update("main_accordion", "Profile Matrisome", icon = original_profile_matrisome_icon)
        completed_sections$profile_matrisome <- FALSE
      } else if (input$active_panel == "feature_selection" && completed_sections$feature_selection) {
        accordion_panel_update("main_accordion", "Select Feature", icon = original_feature_icon)
        completed_sections$feature_selection <- FALSE
      }

      # SECOND: Handle main content navigation
      panel_map <- list(
        "load_data" = "load_data",
        "profile_matrisome" = "profile_matrisome",
        "feature_selection" = "feature_selection",
        "spatial_analysis" = "spatial_analysis",
        "signature_analysis" = "signature_analysis"
        # Note: "about" is handled by a separate observeEvent(input$about,...)
      )

      if (input$active_panel %in% names(panel_map)) {
        nav_select("main_content", panel_map[[input$active_panel]])
      }
    }
  })

  # Also ensure the "About" button navigation is present
  observeEvent(input$about, {
    nav_select(id = "main_content", selected = "about")
  })

  shinyjs::runjs("$('#export_dropdown_btn').prop('disabled', true);") # Disable export button initially
  shinyjs::disable("run_matrisome_profile") # Disable matrisome profile button initially
  shinyjs::disable("sel1")
  shinyjs::disable("sel2")
  shinyjs::disable("gene1")
  shinyjs::disable("gene2")
  shinyjs::disable("run_button")
  shinyjs::disable("viopt")
  shinyjs::disable("minexp")
  shinyjs::disable("kerd")
  shinyjs::disable("cs1")
  shinyjs::disable("cs2")
  shinyjs::disable("extraviz")
  shinyjs::disable("cellTypeSelector")
  shinyjs::disable("export")
  shinyjs::disable("export2")
  shinyjs::disable("scol")
  shinyjs::disable("fsig")
  shinyjs::disable("simsrc")
  shinyjs::disable("corsel")
  shinyjs::disable("trg")
  shinyjs::disable("extinfo")
  shinyjs::disable("run_matrisome_pair_analysis")

  # ---------------------------------------------------------------------------
  # DATA LOADING & PROCESSING
  # ---------------------------------------------------------------------------

  # Observer for loading UPLOADED data
  observeEvent(input$load_upload, {
      req(input$seurat_file)

      # 1. Hard Reset of the entire app state
      session$sendCustomMessage("clearViz", TRUE)
      rv$data_is_loaded <- FALSE
      rv$active_seurat_object <- NULL
      rv$active_sample_name <- NULL
      rv$loaded_collection <- NULL
      rv$analysis_results <- NULL
      rv$feature_analysis_done <- FALSE
      rv$matrisome_analysis_done <- FALSE
      rv$matrisome_pair_results <- NULL
      nav_select("main_content", "load_data")
      shinyjs::runjs("$('#export_dropdown_btn').prop('disabled', true).addClass('initially-disabled');")
      updateRadioButtons(session, "sel1", selected = character(0))
      updateRadioButtons(session, "sel2", selected = character(0))
      shinyjs::runjs("$('.feature-card').removeClass('selected');")

      # Reject uploads >2 GB on disk: a deserialized Seurat is ~2-4x its file
      # size and the worker silently OOMs during readRDS with no log traceback.
      upload_size_bytes <- tryCatch(file.size(input$seurat_file$datapath),
                                    error = function(e) NA_real_)
      max_upload_bytes <- 2 * 1024^3
      if (!is.na(upload_size_bytes) && upload_size_bytes > max_upload_bytes) {
          showNotification(
            sprintf(paste(
              "Uploaded file is %.2f GB; this exceeds the safe load size",
              "(%.2f GB) for the current server tier. Please subset the",
              "object locally or contact the maintainers for a larger tier."),
              upload_size_bytes / 1024^3,
              max_upload_bytes  / 1024^3
            ),
            type = "error", duration = 15
          )
          return()
      }

      # Read the uploaded file
      show_modal_spinner(spin = "self-building-square", color = "blue", text = "Loading uploaded data...")
      # Reclaim memory before allocating the new object; otherwise two large
      # Seurat objects briefly coexist and can OOM the worker. Spinner is up
      # so the user sees feedback during the 1-3s gc pause.
      gc(verbose = FALSE)
      obj <- tryCatch(readRDS(input$seurat_file$datapath), error = function(e) {
          showNotification(paste("Error reading .rds file:", e$message), type = "error"); NULL
      })
      remove_modal_spinner()
      gc(verbose = FALSE)

      if (is.null(obj)) {
          showNotification("Failed to read .rds file.", type = "error")
          return()
      }

      # Handle SpatialExperiment -> convert to Seurat
      if (inherits(obj, "SpatialExperiment")) {
          show_modal_spinner(spin = "self-building-square", color = "blue",
                             text = "Converting SpatialExperiment to Seurat...")
          obj <- tryCatch({
              spe_to_seurat(obj, verbose = FALSE)
          }, error = function(e) {
              showNotification(paste("SPE conversion failed:", e$message), type = "error")
              NULL
          })
          remove_modal_spinner()

          if (is.null(obj)) return()

          showNotification("SpatialExperiment converted successfully!", type = "message")
      }

      if (!inherits(obj, "Seurat")) {
          showNotification("File must be a Seurat or SpatialExperiment object.", type = "error")
          return()
      }

      # Handle VisiumV2 -> convert to VisiumV1 for compatibility with Seurat 5.0.x
      # Access @images directly to avoid Images() method which fails on VisiumV2
      img_names_check <- names(obj@images)
      has_visiumv2 <- tryCatch({
        length(img_names_check) > 0 && any(sapply(img_names_check, function(nm) {
          inherits(obj@images[[nm]], "VisiumV2")
        }))
      }, error = function(e) FALSE)

      rv$is_visiumv2_origin <- has_visiumv2  # Track for point size adjustment

      if (has_visiumv2) {
          show_modal_spinner(spin = "self-building-square", color = "blue",
                             text = "Converting VisiumV2 to VisiumV1 format...")
          obj <- tryCatch({
              convert_visiumv2_to_v1(obj, verbose = TRUE)
          }, error = function(e) {
              showNotification(paste("VisiumV2 conversion failed:", e$message), type = "error")
              obj  # Return original object on failure
          })
          remove_modal_spinner()
          showNotification("VisiumV2 converted to VisiumV1 format.", type = "message")
      }

      # Standardize gene symbols to current HGNC nomenclature
      obj <- standardize_gene_symbols(obj)

      # 3. Check if preprocessing is needed
      domain_ok <- tryCatch("ecm_domain_annotation" %in% colnames(obj@meta.data), error = function(e) FALSE)
      sct_ok <- tryCatch({
          has_sct <- "SCT" %in% SeuratObject::Assays(obj)
          has_scale_data <- has_sct && nrow(LayerData(obj, assay = "SCT", layer = "scale.data")) > 0
          has_sct_model <- has_sct && length(obj[["SCT"]]@SCTModel.list) > 0
          has_sct && (has_scale_data || has_sct_model)
      }, error = function(e) FALSE)
      features_ok <- tryCatch("collagens" %in% colnames(obj@meta.data), error = function(e) FALSE)
      ucell_ok <- tryCatch("Interstitial_UCell" %in% colnames(obj@meta.data), error = function(e) FALSE)

      if (!all(sct_ok, features_ok, domain_ok, ucell_ok)) {
          # Object needs processing, so ask the user via modal
          showModal(modalDialog(
              title = "Additional Processing Required",
              "Your Seurat object is missing key components for full analysis (e.g., SCT data, matrisome scores, ECM annotations).",
              br(), br(),
              "MatriSpace can automatically run these steps. This may take several minutes.",
              br(), br(), strong("Do you want to proceed?"),
              footer = tagList(
                  actionButton("cancel_processing", "Cancel"),
                  actionButton("start_processing", "Proceed", class = "btn-primary")
              )
          ))
          rv$raw_uploaded_object <- obj
      } else {
          # Object is ready, populate the UNIVERSAL source of truth
          obj <- sctype_annotate_ecm(obj)
          if (!"ecm_domain_annotation" %in% colnames(obj@meta.data)) {
            showNotification(
              "ECM niche annotations are unavailable because this dataset does not contain the required marker genes.",
              type = "warning", duration = 10
            )
          }
          obj@meta.data[is.na(obj@meta.data)] <- "not.assigned"
          obj <- slim_seurat_for_app(obj)
          gc(verbose = FALSE)
          rv$active_seurat_object <- obj
          rv$active_sample_name <- tools::file_path_sans_ext(input$seurat_file$name)
          rv$loaded_collection <- "upload"  # Track as uploaded data
          rv$data_is_loaded <- TRUE # This triggers the success observer
      }
  })

  # Observers for the preprocessing modal buttons
  observeEvent(input$start_processing, {
      removeModal()
      req(rv$raw_uploaded_object)

      show_modal_spinner(
        spin = "self-building-square", color = "blue",
        text = tagList(h4("Processing Object..."), p("This may take several minutes."))
      )

      processed_obj <- tryCatch({
        prepare_uploaded_object(
          rv$raw_uploaded_object,
          mat_feat_sigs = matrisome_feature_signatures,
          ecm_ucell_sigs = ecm_ucell_signatures
        )
      }, error = function(e) {
        showNotification(paste("Error during preprocessing:", e$message), type = "error")
        NULL
      })

      # On success, populate the UNIVERSAL source of truth
      if (!is.null(processed_obj)) {
        # Keep the spinner up through slim+gc — 1-3s of blocking work otherwise
        # leaves the user staring at a blank screen between modal close and
        # the first reactive render.
        processed_obj <- slim_seurat_for_app(processed_obj)
        gc(verbose = FALSE)

        rv$active_seurat_object <- processed_obj
        rv$active_sample_name <- tools::file_path_sans_ext(input$seurat_file$name)
        rv$loaded_collection <- "upload"  # Track as uploaded data
        rv$data_is_loaded <- TRUE
      }
      remove_modal_spinner()
  })

  observeEvent(input$cancel_processing, {
      removeModal()
      raw_obj <- rv$raw_uploaded_object
      raw_obj@meta.data[is.na(raw_obj@meta.data)] <- "not.assigned"

      raw_obj <- slim_seurat_for_app(raw_obj)
      gc(verbose = FALSE)

      # Populate the UNIVERSAL source of truth with the raw object
      rv$active_seurat_object <- raw_obj
      rv$active_sample_name <- tools::file_path_sans_ext(input$seurat_file$name)
      rv$loaded_collection <- "upload"  # Track as uploaded data
      rv$data_is_loaded <- TRUE
      showNotification("Processing cancelled. Some analysis tabs may be unavailable.", type = "warning", duration = 8)
  })

  # --- Spot leakage (bleed) correction ---------------------------------------
  # The switch lives in the Load Data panel and is OFF by default.
  #
  # The correction is GENE-SCOPED: it is applied only to the features the user
  # actually selected, at the moment an analysis is launched. This is exact,
  # not an approximation of a whole-matrix correction; the neighbour mean of a
  # gene involves that gene only, so correcting row g in isolation gives the
  # same numbers as correcting the full matrix and reading row g back out.
  # It is simply orders of magnitude cheaper: the whole-matrix path costs ~14 s
  # per layer on a 20k x 5k object, the scoped path a fraction of a second.
  #
  # rv$active_seurat_object and d0() therefore always hold the PRISTINE object;
  # correction happens inside leak_apply() at the three analysis entry points.

  # Sliders are debounced so that dragging alpha does not thrash the observers
  # that clear cached results.
  leak_alpha <- reactive({
    a <- input$leak_alpha
    if (is.null(a) || !is.finite(a)) 0.10 else as.numeric(a)
  })
  leak_alpha_d <- debounce(leak_alpha, 900)

  leak_k <- reactive({
    kk <- input$leak_k
    if (is.null(kk) || !is.finite(kk)) 6L else as.integer(kk)
  })
  leak_k_d <- debounce(leak_k, 900)

  leak_status_val <- reactiveVal(NULL)

  #' Correct `genes` in `obj` if the switch is on; otherwise return `obj`.
  #' Called from observers only, so writing to leak_status_val here is safe.
  leak_apply <- function(obj, genes, context = "analysis") {
    if (!isTRUE(input$leak_correct)) return(obj)

    genes <- unique(as.character(genes))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (length(genes) == 0L) {
      leak_status_val(list(ok = FALSE,
                           msg = sprintf("No correctable features in the %s.", context)))
      return(obj)
    }

    t0 <- Sys.time()
    res <- tryCatch(
      matrispace_spotclean(obj, genes = genes,
                           alpha = leak_alpha(), k = leak_k()),
      error = function(e) list(applied = FALSE, message = conditionMessage(e))
    )

    if (!isTRUE(res$applied)) {
      leak_status_val(list(ok = FALSE, msg = res$message))
      showNotification(paste("Leakage correction not applied:", res$message),
                       type = "warning", duration = 10)
      return(obj)
    }

    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    leak_status_val(list(
      ok  = TRUE,
      msg = paste0(res$message, " (", context, "): ",
                   formatC(elapsed, format = "f", digits = 2), " s")))
    res$object
  }

  output$leak_status <- renderUI({
    st <- leak_status_val()
    if (is.null(st)) {
      return(tags$span(style = "color:#7880a0;",
                       "Applied to the selected features when an analysis is run."))
    }
    tags$span(style = sprintf("color:%s;", if (isTRUE(st$ok)) "#1e7a45" else "#c85050"),
              if (isTRUE(st$ok)) paste("Applied:", st$msg) else st$msg)
  })
  outputOptions(output, "leak_status", suspendWhenHidden = FALSE)

  # Any change to the correction settings invalidates every cached result.
  observeEvent(list(input$leak_correct, leak_alpha_d(), leak_k_d()), {
    leak_status_val(NULL)
    if (!isTRUE(rv$data_is_loaded)) return()
    rv$analysis_results        <- NULL
    rv$matrisome_results       <- NULL
    rv$matrisome_pair_results  <- NULL
    rv$feature_analysis_done   <- FALSE
    rv$matrisome_analysis_done <- FALSE
    showNotification(
      sprintf("Leakage correction %s. Previous results cleared; please re-run the analyses.",
              if (isTRUE(input$leak_correct)) "updated" else "switched off"),
      type = "message", duration = 7)
  }, ignoreInit = TRUE)

  # Data accessors (always the uncorrected object)
  d0 <- reactive({
      req(rv$active_seurat_object)
  })

  output$ecm_annotation_status <- renderUI({
    req(rv$data_is_loaded)
    if ("ecm_domain_annotation" %in% colnames(d0()@meta.data)) return(NULL)
    div(
      class = "alert alert-warning",
      "ECM niche annotations are unavailable because this dataset does not contain the required marker genes."
    )
  })

  nm <- reactive({
      req(rv$active_sample_name)
  })

  # VisiumV2 point size adjustment using multiplier (slider remains responsive)
  effective_pt_size <- function(base_size) {
    if (isTRUE(rv$is_visiumv2_origin)) {
      base_size * 3  # Scale up proportionally for VisiumV2
    } else {
      base_size
    }
  }

  # Enable load button for file uploads
  observeEvent(input$seurat_file, {
    req(input$seurat_file)
    if (grepl("\\.rds$", input$seurat_file$name, ignore.case = TRUE)) {
      shinyjs::enable("load_upload")
    } else {
      shinyjs::disable("load_upload")
      showNotification("Please upload an .rds file", type = "error")
    }
  })

  # Master controller for FEATURE analysis triggered by the run button
  observeEvent(input$run_button, {
    if (length(input$sel1) == 0 || !nzchar(input$sel1)) {
      showNotification("Please choose a primary feature type before running.",
                       type = "warning")
      return()
    }

    sel1_val <- input$sel1
    sel2_val <- if (length(input$sel2) == 0 || !nzchar(input$sel2)) "none" else input$sel2
    custom1 <- identical(sel1_val, "matrisome signature") && identical(input$list_source1, "custom")
    custom2 <- identical(sel2_val, "matrisome signature") && identical(input$list_source2, "custom")
    allowed_genes <- available_matrisome_genes(d0(), mobj)
    custom_sets <- list()

    if (custom1) {
      genes <- unique(input$custom_genes1)
      if (length(setdiff(genes, allowed_genes)) > 0) {
        showNotification("The primary custom list contains genes outside the detected matrisome.", type = "error")
        return()
      }
      if (length(genes) < 2) {
        showNotification("Select at least two genes for the primary custom matrisome list.", type = "warning")
        return()
      }
      custom_sets$customlist1 <- genes
    }

    if (custom2) {
      genes <- unique(input$custom_genes2)
      if (length(setdiff(genes, allowed_genes)) > 0) {
        showNotification("The secondary custom list contains genes outside the detected matrisome.", type = "error")
        return()
      }
      if (length(genes) < 2) {
        showNotification("Select at least two genes for the secondary custom matrisome list.", type = "warning")
        return()
      }
      custom_sets$customlist2 <- genes
    }

    # Curated matrisome signatures are read from pre-computed metadata columns,
    # which the leakage correction does not touch. Say so rather than silently
    # mixing corrected and uncorrected features in the same analysis.
    if (isTRUE(input$leak_correct) &&
        ((identical(sel1_val, "matrisome signature") && identical(input$list_source1, "curated")) ||
         (identical(sel2_val, "matrisome signature") && identical(input$list_source2, "curated")))) {
      showNotification(
        paste("Leakage correction is ON, but curated matrisome signature scores come from",
              "pre-computed metadata and are NOT corrected. To score a gene list against the",
              "corrected matrix, use the custom matrisome list option instead."),
        type = "warning", duration = 14)
    }

    gene1_val <- if (custom1) "Custom matrisome list" else if (length(input$gene1) == 0) "" else input$gene1
    gene2_val <- if (custom2) "Custom matrisome list" else if (length(input$gene2) == 0) "" else input$gene2

    if (!nzchar(gene1_val)) {
      showNotification("Please select a primary feature before running.", type = "warning")
      return()
    }

    # Show the main loading modal with a random ECM quote.
    random_quote <- sample(ecm_quotes, 1)
    show_modal_spinner(
      spin = "self-building-square",
      color = "blue",
      text = tagList(
        h4("Running Feature Analysis..."),
        p("This may take a moment... or two?"),
        hr(),
        tags$i(random_quote) # Display the quote in italics
      )
    )
    on.exit(remove_modal_spinner(), add = TRUE)

    rv$analysis_results <- NULL
    rv$custom_gene_sets <- list()
    rv$primary_contrast_state <- "none"
    rv$secondary_contrast_state <- "none"
    gc(verbose = FALSE)

    # Gene-scoped leakage correction: only the features actually selected are
    # corrected, and only when the switch is on. Curated signatures read
    # pre-computed metadata, so there is nothing to correct for them (warned above).
    leak_genes <- character(0)
    if (sel1_val %in% c("matrisome gene", "any gene")) leak_genes <- c(leak_genes, gene1_val)
    if (sel2_val %in% c("matrisome gene", "any gene")) leak_genes <- c(leak_genes, gene2_val)
    if (length(custom_sets) > 0) leak_genes <- c(leak_genes, unlist(custom_sets, use.names = FALSE))
    obj_run <- if (length(leak_genes) > 0) {
      leak_apply(d0(), leak_genes, context = "feature analysis")
    } else d0()

    run_result <- tryCatch({
      score_result <- if (length(custom_sets) > 0) {
        score_custom_gene_sets(obj_run, custom_sets)
      } else {
        NULL
      }

      feature1_values <- if ("customlist1" %in% names(custom_sets)) {
        values <- score_result$scores[["customlist1"]]
        names(values) <- rownames(score_result$scores)
        values
      } else NULL

      feature2_values <- if ("customlist2" %in% names(custom_sets)) {
        values <- score_result$scores[["customlist2"]]
        names(values) <- rownames(score_result$scores)
        values
      } else NULL

      list(
        data = addfeat(
          obj_run, gene1_val, sel1_val, gene2_val, sel2_val,
          feature1_values = feature1_values,
          feature2_values = feature2_values
        ),
        scoring = score_result
      )
    }, error = function(e) {
        showNotification(
          paste("Feature analysis failed:", conditionMessage(e)),
          type = "error", duration = 10
        )
        NULL
      }
    )
    if (is.null(run_result)) return()
    results_df <- run_result$data

    # Clean annotation column: handle NAs properly
    original_annot_col_name <- ann()

    if (original_annot_col_name %in% colnames(results_df)) {
      # Convert to character to allow modification
      results_df[[original_annot_col_name]] <- as.character(results_df[[original_annot_col_name]])
      # Replace NAs with "not.assigned"
      results_df[[original_annot_col_name]][is.na(results_df[[original_annot_col_name]])] <- "not.assigned"
    }

    # Store results in our reactiveVal container
    rv$analysis_results <- list(
      data  = results_df,
      gene1 = gene1_val,
      gene2 = gene2_val,
      sel2  = sel2_val,
      custom_gene_sets = custom_sets,
      custom_score_assay = if (is.null(run_result$scoring)) NULL else run_result$scoring$assay,
      custom_score_layer = if (is.null(run_result$scoring)) NULL else run_result$scoring$layer,
      custom_score_diagnostics = if (is.null(run_result$scoring)) list() else run_result$scoring$diagnostics,
      custom_contrast_scores = if (is.null(run_result$scoring)) list() else run_result$scoring$contrast_scores
    )
    rv$custom_gene_sets <- custom_sets
    if (isTRUE(rv$analysis_results$custom_score_diagnostics$customlist1$flagged)) {
      rv$primary_contrast_state <- "warning"
    }
    if (isTRUE(rv$analysis_results$custom_score_diagnostics$customlist2$flagged)) {
      rv$secondary_contrast_state <- "warning"
    }

    rv$feature_analysis_done <- TRUE

    print_celltype_colors(d0(), original_annot_col_name)
  })

  # Separate observer for MATRISOME PROFILE analysis
  observeEvent(input$run_matrisome_profile, {
    # Select a random quote from the collection
    random_quote <- sample(ecm_quotes, 1)

    # Show a loading modal
    show_modal_spinner(
      spin = "self-building-square", color = "blue",
      text = tagList(
        h4("Profiling the Extracellular Matrix..."),
        p("This may take a moment... or two?"),
        hr(),
        tags$i(random_quote)
      )
    )
    # Always remove the spinner and roll back matrisome_analysis_done on
    # failure so downstream `req(rv$matrisome_analysis_done)` reactives
    # don't fire against half-built data.
    matrisome_succeeded <- FALSE
    on.exit({
      remove_modal_spinner()
      if (!matrisome_succeeded) {
        rv$matrisome_analysis_done <- FALSE
      }
    }, add = TRUE)

    err <- tryCatch({

      # Every matrisome score this observer builds is derived from matrisome
      # genes only, so scope the correction to the matrisome genes detected in
      # the object rather than the whole transcriptome.
      seurat_obj_with_scores <- leak_apply(
        d0(),
        available_matrisome_genes(d0(), mobj),
        context = "matrisome profile"
      )
      raw_counts <- safe_get_assay_data(seurat_obj_with_scores, slot = "counts")
      if (is.null(raw_counts)) {
        showNotification("Could not access counts data for matrisome scoring.",
                         type = "error", duration = 10)
        return(invisible())
      }

      main_divisions_map <- c("ECM Glycoproteins" = "glycoprotein", "Collagens" = "collagen", "Proteoglycans" = "proteoglycan", "ECM-affiliated Proteins" = "affliated_protein", "ECM Regulators" = "regulator", "Secreted Factors" = "secreted_factor")

      # Map for Matrisome subcategories card (functional categories)
      ecm_subcategories_map <- c(
        "Perivascular" = "perivascular",
        "Hemostasis" = "hemostasis",
        "Elastic fibers" = "elastic_fibers",
        "Growth-factor binding" = "gf_binding"
      )

      # Map for Matrisome gene families card
      ecm_gene_families_map <- c(
        "Laminins" = "laminin",
        "Matricellular proteins" = "matricellular",
        "Syndecans" = "syndecan",
        "Glypicans" = "glypican"
      )

      # Map new display names to matrisome data lookup keys (preserves compatibility with matrisome.rds)
      data_lookup_map <- c(
        "Perivascular" = "Peri-vascular ECM",
        "Growth-factor binding" = "Growth Factor-binding",
        "Laminins" = "Laminin",
        "Syndecans" = "Syndecan",
        "Glypicans" = "Glypican"
      )

      # Accumulate all 28 scores in a plain data.frame and merge into
      # @meta.data in one AddMetaData call: direct `[[<-` writes would
      # copy-on-write the whole meta.data slot once per column.
      n_spots <- ncol(seurat_obj_with_scores)
      score_frame <- data.frame(row.names = colnames(seurat_obj_with_scores))

      # Use the detailed progress bar for this long-running task
      withProgress(message = 'Calculating Matrisome Scores', value = 0, {
        # Combine all sub-maps for progress calculation
        all_sub_maps <- c(ecm_subcategories_map, ecm_gene_families_map)
        total_items <- length(main_divisions_map) + length(all_sub_maps)

        process_one <- function(display_name, internal_name, gene_list) {
          df <- tryCatch(
            process_matrisome_expression(gene_list, seurat_obj_with_scores,
                                         internal_name, counts_matrix = raw_counts, matrisome = matrisome),
            error = function(e) {
              warning(sprintf("Matrisome scoring failed for %s: %s",
                              display_name, conditionMessage(e)))
              NULL
            }
          )
          if (!is.null(df)) {
            score_frame[[paste0(internal_name, "_robust_score")]]    <<- df$robust
            score_frame[[paste0(internal_name, "_log_scaled_score")]] <<- df$log_scaled
          }
          incProgress(1/total_items, detail = paste("Processing", display_name))
        }

        # Process main divisions
        for (display_name in names(main_divisions_map)) {
          internal_name <- main_divisions_map[[display_name]]
          gene_list <- matrisome$gene[matrisome$notes == display_name]
          process_one(display_name, internal_name, gene_list)
        }

        # Process all subcategories and gene families
        for (display_name in names(all_sub_maps)) {
          internal_name <- all_sub_maps[[display_name]]
          # Use grepl on ecm_subcategory column to handle semicolon-separated tags
          lookup_name <- if (display_name %in% names(data_lookup_map)) data_lookup_map[[display_name]] else display_name
          gene_list <- matrisome$gene[grepl(lookup_name, matrisome$ecm_subcategory, fixed = TRUE)]
          process_one(display_name, internal_name, gene_list)
        }
      })

      # Drop the counts reference before plot construction so it doesn't
      # outlive the loop above.
      raw_counts <- NULL
      gc(verbose = FALSE)

      current_annotation_vector <- cleaned_annotation()
      score_frame$Annotation <- align_to_cells(current_annotation_vector, rownames(score_frame))
      seurat_obj_with_scores <- SeuratObject::AddMetaData(
        seurat_obj_with_scores, metadata = score_frame
      )
      rm(score_frame); gc(verbose = FALSE)

      autocorr_results <- tryCatch(calculate_autocorrelation(seurat_obj_with_scores),
                                   error = function(e) {
                                     warning(sprintf("Autocorrelation failed: %s", conditionMessage(e)))
                                     NULL
                                   })

      # Compute ECM-niche hotspot scores before plot data so the plot helper
      # runs on a fully-prepared object (no post-hoc re-attaching).
      ecm_sig_cols <- c("Interstitial_UCell", "Basement_UCell")
      hotspot_frame <- data.frame(row.names = colnames(seurat_obj_with_scores))
      for (col_name in ecm_sig_cols) {
        if (col_name %in% colnames(seurat_obj_with_scores@meta.data)) {
          raw_scores <- seurat_obj_with_scores@meta.data[[col_name]]
          log_scores <- log1p(raw_scores)
          ls_min <- min(log_scores, na.rm = TRUE)
          ls_max <- max(log_scores, na.rm = TRUE)
          hotspot_scores <- if (isTRUE(ls_max > ls_min))
            (log_scores - ls_min) / (ls_max - ls_min)
          else
            rep(0, length(log_scores))
          new_col_name <- sub("_UCell", "_Hotspot_Score", col_name)
          hotspot_frame[[new_col_name]] <- hotspot_scores
        }
      }
      if (ncol(hotspot_frame) > 0) {
        seurat_obj_with_scores <- SeuratObject::AddMetaData(
          seurat_obj_with_scores, metadata = hotspot_frame
        )
      }
      rm(hotspot_frame); gc(verbose = FALSE)

      # Build lightweight plot-data per category (no plotly objects yet);
      # actual plotly construction is deferred to renderPlotly() so only the
      # visible category materializes. Cache coords once for all 14 calls.
      coords_cache <- align_spatial_coordinates(seurat_obj_with_scores)

      plots_list <- list()
      all_maps <- c(main_divisions_map, ecm_subcategories_map, ecm_gene_families_map)

      for (display_name in names(all_maps)) {
        internal_name <- all_maps[[display_name]]
        robust_col <- paste0(internal_name, "_robust_score")

        if (robust_col %in% colnames(seurat_obj_with_scores@meta.data)) {
          type <- if (display_name %in% names(main_divisions_map)) "category" else "subcategory"
          plots_list[[internal_name]] <- tryCatch(
            compute_matrisome_plot_data(
              seurat_obj_with_scores, display_name, internal_name,
              type = type, annotation_col = "Annotation",
              coords_cache = coords_cache
            ),
            error = function(e) {
              warning(sprintf("Plot data computation failed for %s: %s",
                              display_name, conditionMessage(e)))
              NULL
            }
          )
        } else {
          plots_list[[internal_name]] <- NULL
        }
      }
      rm(coords_cache); gc(verbose = FALSE)

      # Attach a Seurat reference into each plot entry for backward compat
      # with consumers that read rv$matrisome_results$<cat>$seurat_object.
      # Cheap: all entries share the same underlying SEXP.
      for (nm in names(plots_list)) {
        if (!is.null(plots_list[[nm]])) {
          plots_list[[nm]]$seurat_object <- seurat_obj_with_scores
        }
      }

      rv$matrisome_results <- c(
        plots_list,
        list(
          seurat_object  = seurat_obj_with_scores,   # canonical top-level reference
          autocorrelation = autocorr_results
        )
      )

      # Update UI and state
      rv$matrisome_analysis_done <- TRUE
      matrisome_succeeded <- TRUE
      accordion_panel_update("main_accordion", "Profile Matrisome", icon = success_icon)
      completed_sections$profile_matrisome <- TRUE
      nav_select("main_content", "profile_matrisome")

      NULL
    },
    error = function(e) {
      showNotification(
        paste("Matrisome profiling failed:", conditionMessage(e)),
        type = "error", duration = 12
      )
      conditionMessage(e)
    })

    invisible(err)
  })

  # --- [ D3.js INTERACTIVE ANNOTATION PLOT LOGIC - MULTI-INSTANCE ] ---

  observe({
    req(rv$feature_analysis_done, input$run_button, ann())
    seurat_obj <- d0()
    annotation_col <- ann()
    colors_list <- as.list(color_map())
    showNotification("Updating main interactive view...", type = "message", duration = 2)

    tryCatch({
      spatial_image <- seurat_obj@images[[1]]
      req(inherits(spatial_image, "VisiumV1"))
      img_raster <- spatial_image@image
      scale_factor_lowres <- spatial_image@scale.factors$lowres

      temp_image_path <- tempfile(fileext = ".png")
      png(temp_image_path, width = ncol(img_raster), height = nrow(img_raster))
      par(mar = c(0, 0, 0, 0)); plot(as.raster(img_raster)); dev.off()
      image_uri <- knitr::image_uri(temp_image_path)
      unlink(temp_image_path)

      coords <- spatial_image@coordinates; coords$barcode <- rownames(coords)
      meta <- seurat_obj@meta.data; meta$barcode <- rownames(meta)
      if (!"in_tissue" %in% colnames(coords)) {
        coords$in_tissue <- if ("tissue" %in% colnames(coords)) coords$tissue else if ("in_tissue" %in% colnames(meta)) meta[match(coords$barcode, meta$barcode), "in_tissue"] else 1
      }
      meta$in_tissue <- NULL
      spot_data <- merge(coords, meta, by = "barcode")

      spot_data[[annotation_col]] <- as.character(spot_data[[annotation_col]])
      spot_data[[annotation_col]][is.na(spot_data[[annotation_col]])] <- "not.assigned"

      payload <- list(
        target_id = "main",
        spots = spot_data, imageUrl = image_uri, clusterColors = colors_list,
        imageDims = list(width = ncol(img_raster), height = nrow(img_raster)),
        lowresSF = scale_factor_lowres, annotationName = annotation_col
      )
      session$sendCustomMessage(type = "renderVizData", message = jsonlite::toJSON(payload, auto_unbox = TRUE, na = "string"))

    }, error = function(e) {
      showNotification(paste("Error updating main interactive plot:", e$message), type = "error", duration = 8)
    })
  })

  observeEvent(input$run_matrisome_profile, {
    req(rv$matrisome_analysis_done)
    seurat_obj <- d0()

    annotation_col <- "ecm_domain_annotation"
    if (!annotation_col %in% colnames(seurat_obj@meta.data)) {
      showNotification("ECM data for interactive plot not found.", type = "warning", duration = 5)
      return()
    }

    colors_list <- as.list(create_ecm_color_map(unique(seurat_obj@meta.data[[annotation_col]])))

    tryCatch({
      spatial_image <- seurat_obj@images[[1]]
      img_raster <- spatial_image@image
      scale_factor_lowres <- spatial_image@scale.factors$lowres

      temp_image_path <- tempfile(fileext = ".png")
      png(temp_image_path, width = ncol(img_raster), height = nrow(img_raster))
      par(mar = c(0, 0, 0, 0)); plot(as.raster(img_raster)); dev.off()
      image_uri <- knitr::image_uri(temp_image_path)
      unlink(temp_image_path)

      coords <- spatial_image@coordinates; coords$barcode <- rownames(coords)
      meta <- seurat_obj@meta.data; meta$barcode <- rownames(meta)
      if (!"in_tissue" %in% colnames(coords)) {
        coords$in_tissue <- if ("tissue" %in% colnames(coords)) coords$tissue else if ("in_tissue" %in% colnames(meta)) meta[match(coords$barcode, meta$barcode), "in_tissue"] else 1
      }
      meta$in_tissue <- NULL
      spot_data <- merge(coords, meta, by = "barcode")

      spot_data[[annotation_col]] <- as.character(spot_data[[annotation_col]])
      spot_data[[annotation_col]][is.na(spot_data[[annotation_col]])] <- "not.assigned"

      payload <- list(
        target_id = "ecm",
        spots = spot_data, imageUrl = image_uri, clusterColors = colors_list,
        imageDims = list(width = ncol(img_raster), height = nrow(img_raster)),
        lowresSF = scale_factor_lowres, annotationName = annotation_col
      )
      session$sendCustomMessage(type = "renderVizData", message = jsonlite::toJSON(payload, auto_unbox = TRUE, na = "string"))

    }, error = function(e) {
      showNotification(paste("Error rendering ECM interactive plot:", e$message), type = "error", duration = 8)
    })
  }, ignoreInit = TRUE)

  observeEvent(input$selectedCellTypes_main, {
    req(rv$feature_analysis_done)
    session$sendCustomMessage("updateVisibleAnnotations", list(target_id = "main", annotations = input$selectedCellTypes_main))
  }, ignoreNULL = FALSE)

  observeEvent(input$spot_opacity_main, {
    session$sendCustomMessage("updateSpotStyle", list(target_id = "main", style = "opacity", value = input$spot_opacity_main))
  }, ignoreInit = TRUE)

  observeEvent(input$spot_size_main, {
    session$sendCustomMessage("updateSpotStyle", list(target_id = "main", style = "size", value = input$spot_size_main))
  }, ignoreInit = TRUE)

  observeEvent(input$spot_borders_main, {
    session$sendCustomMessage("updateSpotStyle", list(target_id = "main", style = "borders", value = input$spot_borders_main))
  }, ignoreInit = TRUE)

  observeEvent(input$spot_opacity_ecm, {
    session$sendCustomMessage("updateSpotStyle", list(target_id = "ecm", style = "opacity", value = input$spot_opacity_ecm))
  }, ignoreInit = TRUE)

  observeEvent(input$spot_size_ecm, {
    session$sendCustomMessage("updateSpotStyle", list(target_id = "ecm", style = "size", value = input$spot_size_ecm))
  }, ignoreInit = TRUE)

  observeEvent(input$spot_borders_ecm, {
    session$sendCustomMessage("updateSpotStyle", list(target_id = "ecm", style = "borders", value = input$spot_borders_ecm))
  }, ignoreInit = TRUE)

  output$legend_plot_ecm <- renderPlot({
    req(rv$matrisome_analysis_done)
    seurat_obj <- d0()

    annotation_col <- "ecm_domain_annotation"
    if (!annotation_col %in% colnames(seurat_obj@meta.data)) return()

    tryCatch({
      domain_levels <- sort(unique(seurat_obj@meta.data[[annotation_col]]))
      color_map <- create_ecm_color_map(domain_levels)

      legend_df <- data.frame(domain = factor(domain_levels, levels = domain_levels))

      p <- ggplot(legend_df, aes(x = 1, y = domain, color = domain)) +
        geom_point(alpha = 0) +
        scale_color_manual(values = color_map, name = "ECM niches") +
        guides(color = guide_legend(override.aes = list(shape = 15, size = 10, alpha = 1))) +
        theme_void() +
        theme(
          legend.position = "left",
          legend.title = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 14)
        )

      print(p)

    }, error = function(e) {
      return(NULL)
    })
  }, bg = "transparent")

  # --- [ END D3 LOGIC ] ---

  # Note: Clearing of plots when uploading files is now handled in the main load_upload observer above
  # The KS killswitch has been replaced by rv$feature_analysis_done and rv$matrisome_analysis_done flags

  # Interactive cell type/annotation type selector
  #----------------------------------------------------------------------------------
  cellTypesReactive <- reactive({
    req(cleaned_annotation()) # Depend on the clean vector
    counts <- table(cleaned_annotation())
    counts[order(names(counts))]
  })

  output$cellTypeSelector_main <- renderUI({
    req(cellTypesReactive(), color_map())
    cell_type_counts <- cellTypesReactive()
    new_cell_types <- names(cell_type_counts)
    stable_colors <- color_map()
    current_selection <- isolate(input$selectedCellTypes_main)

    new_selection <- if (is.null(current_selection)) {
      new_cell_types
    } else {
      preserved <- intersect(current_selection, new_cell_types)
      if (length(preserved) == 0) new_cell_types else preserved
    }

    choices_html <- sapply(new_cell_types, function(type) {
      color <- stable_colors[type]
      count <- cell_type_counts[type]
      sprintf("<div><span class='color-swatch' style='background-color:%s;'></span> %s (%s)</div>", color, type, format(count, big.mark=","))
    })

    selectizeInput(
      "selectedCellTypes_main", NULL,
      choices = setNames(new_cell_types, choices_html),
      selected = new_selection, multiple = TRUE,
      options = list(
        plugins = list('remove_button'), create = FALSE, placeholder = 'Click to select annotations',
        render = I("{ item: function(item, escape) { return item.label; }, option: function(item, escape) { return item.label; } }")
      )
    )
  })

  selectedCellTypes <- reactive({ input$selectedCellTypes_main })

  observeEvent(input$selectAllAnnotations_main, {
    all_types <- names(cellTypesReactive())
    updateSelectizeInput(session, "selectedCellTypes_main", selected = all_types)
  })

  observeEvent(input$deselectAllAnnotations_main, {
    updateSelectizeInput(session, "selectedCellTypes_main", selected = character(0))
  })
  #----------------------------------------------------------------------------------

  # Fetch metadata columns from uploaded objects
  available_metadata <- reactive({
    req(d0())
    meta_df <- d0()@meta.data
    MAX_CATEGORIES <- 100

    valid_cols <- sapply(colnames(meta_df), function(col) {
      is_categorical <- is.character(meta_df[[col]]) || is.factor(meta_df[[col]])
      if (!is_categorical) return(FALSE)
      num_unique <- length(unique(na.omit(meta_df[[col]])))
      return(num_unique > 1 && num_unique < MAX_CATEGORIES)
    })

    valid_names <- names(which(valid_cols))
    exclude_cols <- c("orig.ident", "seurat_clusters", grep("_LISA$", valid_names, value = TRUE))
    return(valid_names[!valid_names %in% exclude_cols])
  })

  # Reactive for selected annotation column
  ann <- reactive({
    selected <- input$annot
    req(length(selected) == 1L, nzchar(selected), selected %in% colnames(d0()@meta.data))
    selected
  })

  # Clean annotation vector (handles NAs)
  cleaned_annotation <- reactive({
    req(d0(), ann())
    annot_vector <- as.character(d0()@meta.data[[ann()]])
    annot_vector[is.na(annot_vector)] <- "not.assigned"
    return(annot_vector)
  })

  # Stable color map for annotations
  color_map <- reactive({
    req(d0(), ann())

    # Get all possible unique annotation levels for the selected column
    all_levels <- unique(as.character(d0()@meta.data[[ann()]]))

    # Use the new helper function to generate the map
    create_custom_color_map(all_levels)
  })

  # Track data source (always upload in offline mode)
  loaded_data_source <- reactiveVal("upload")

  # Metadata selector UI
  output$metadata_selector <- renderUI({
    req(d0())
    selectInput(
      "annot",
      "Choose annotation type:",
      choices = available_metadata(),
      selected = available_metadata()[1]
    )
  })


  # Value box reactives
  sample_name_val <- reactive({
    req(d0(), nm())
    if(is.null(nm())) return("No data")
    name <- tools::toTitleCase(gsub("_", " ", nm()))
    if(nchar(name) > 25 && grepl(" - ", name)) {
      HTML(sub(" - ", "<br>", name))
    } else {
      name
    }
  })

  spot_info_val <- reactive({
    req(d0())
    tryCatch({
      n_spots <- nrow(d0()@meta.data)
      n_features <- safe_get_nfeatures(d0(), assay = DefaultAssay(d0()))
      paste0(formatC(n_spots, format="d", big.mark=","), " / ",
             formatC(n_features, format="d", big.mark=","))
    }, error = function(e) { "No data" })
  })

  matrisome_count_val <- reactive({
    req(d0())
    tryCatch({
      genes <- safe_get_rownames(d0(), assay = DefaultAssay(d0()))
      formatC(sum(genes %in% mobj), format="d", big.mark=",")
    }, error = function(e) { "0" })
  })


  # --- Conditionally render the value boxes ---
  output$value_boxes_ui <- renderUI({

    if (isTRUE(rv$data_is_loaded)) {
      layout_column_wrap(
        width = 1/3,
        flippable_value_box(
          front_ui = bslib::tooltip(
            value_box(
              title = tags$div(
                class = "d-flex justify-content-between align-items-center w-100",
                "Tissue", # Main title
                bsicons::bs_icon("arrow-repeat", class = "opacity-50") # Flip icon
              ),
              value = sample_name_val(),
              showcase = bsicons::bs_icon("binoculars"),
              theme = "primary",
              full_screen = FALSE
            ),
            sample_name_val(), # The tooltip shows the full name on hover
            placement = "bottom"
          ),
          back_ui = uiOutput("tissue_info_back")
        ),
        value_box(
          title = "Spots / Genes",
          value = spot_info_val(),
          showcase = bsicons::bs_icon("bounding-box-circles"),
          theme = "secondary",
          full_screen = FALSE
        ),
        value_box(
          title = "Matrisome genes",
          value = matrisome_count_val(),
          showcase = bsicons::bs_icon("bezier"),
          theme = "info",
          class = "info-bg-force text-white-force",
          full_screen = FALSE
        )
      )

    } else {
      # --- NO DATA: Render a static, blank, and equal-height version ---
      # Using HTML non-breaking space ensures the boxes maintain their height.
      layout_column_wrap(
        width = 1/3,
        div(class = "flip-card disabled",
            value_box(title = "Tissue", value = HTML("&nbsp;"), showcase = bsicons::bs_icon("binoculars"), theme = "primary")
        ),
        div(class = "flip-card disabled",
            value_box(title = "Spots / Genes", value = HTML("&nbsp;"), showcase = bsicons::bs_icon("bounding-box-circles"), theme = "secondary")
        ),
        div(class = "flip-card disabled",
            value_box(
              title = "Matrisome genes",
              value = HTML("&nbsp;"),
              showcase = bsicons::bs_icon("bezier"),
              theme = "info",
              class = "info-bg-force text-white-force"
            )
        )
      )
    }
  })

  # Render the back of the flippable tissue card
  output$tissue_info_back <- renderUI({
    req(nm(), d0())
    n_spots <- nrow(d0()@meta.data)
    n_features <- safe_get_nfeatures(d0(), assay = DefaultAssay(d0()))
    tagList(
      tags$div(tags$strong("Sample:"), nm()),
      tags$div(tags$strong("Spots:"), formatC(n_spots, format="d", big.mark=",")),
      tags$div(tags$strong("Genes:"), formatC(n_features, format="d", big.mark=","))
    )
  })


  # --- Reactive expression that generates the tissue plot ONCE ---
  # This is the "recipe" for the plot. It's cached for performance.
  generate_tissue_plot <- reactive({
    req(d0())

    # Generate the base plot with invisible spots and return it
    SpatialDimPlot(d0(), label = FALSE, alpha = 0.0, crop = FALSE) + NoLegend()

  }) %>% bindCache(nm()) # Cache the plot based on the sample name

  # --- Output for the "Load Data" view ---
  output$tissue_image_plot <- renderPlot({
    req(isTRUE(rv$data_is_loaded))
    withProgress(message = 'Generating tissue image...', value = 0, {
      print(generate_tissue_plot())
    })
  })

  cluster.plot <- function(){
    req(d0())
    req(nm())

    # Find a suitable cluster column
    meta_cols <- colnames(d0()@meta.data)
    group_col <- NULL

    if ("seurat_clusters" %in% meta_cols) {
      group_col <- "seurat_clusters"
    } else {
      # Look for any column with "cluster" in the name (case-insensitive)
      cluster_cols <- meta_cols[grepl("cluster", meta_cols, ignore.case = TRUE)]
      if (length(cluster_cols) > 0) {
        group_col <- cluster_cols[1]
      }
    }

    # If no cluster column found, return a placeholder plot
    if (is.null(group_col)) {
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5, label = "No cluster column found",
                        size = 5, hjust = 0.5, color = "grey50") +
               theme_void())
    }

    p0 <- SpatialDimPlot(d0(), group.by = group_col, pt.size.factor = effective_pt_size(1.6)) +
      theme(legend.position = "right",
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 8))
    return(p0)
  }

  output$cluster_plot <- renderPlot({
    req(d0())
    req(nm())
    withProgress(message = 'Loading tissue cluster', value = 0, {
      req(cluster.plot)
      print(cluster.plot())
    })
  }) %>% bindCache(nm())

  # Generate the tissue link button (no manifest in offline mode)
  generate_tissue_link <- reactive({
    req(d0(), nm())
    tags$button(
      "Learn more about this tissue sample",
      class = "btn btn-outline-primary btn-sm",
      style = "margin: 0; display: inline-block;",
      disabled = TRUE
    )
  })

  # --- Output for the "Load Data" view ---
  output$tissue_link <- renderUI({
    generate_tissue_link()
  })

  # Second output for the Feature Selection view
  output$tissue_link_features <- renderUI({
    generate_tissue_link()
  })

  # After loading the sample, enable the card selectors
  observeEvent(d0(), {
    req(d0())
    shinyjs::enable("annot")
    shinyjs::enable("sel1")
    shinyjs::enable("sel2")
    shinyjs::enable("run_button")
    shinyjs::enable("run_matrisome_profile")
    if ("ecm_domain_annotation" %in% colnames(d0()@meta.data)) {
      shinyjs::enable("run_matrisome_pair_analysis")
    } else {
      shinyjs::disable("run_matrisome_pair_analysis")
    }

    # Reset selections when new data is loaded
    updateRadioButtons(session, "sel1", selected = character(0))
    updateRadioButtons(session, "sel2", selected = character(0))
    rv$custom_gene_sets <- list()
  })

  # Conditionally render the PRIMARY feature dropdown
  output$ui1 <- renderUI({
    # Check if data is loaded and a selection has been made
    if (!is.null(input$sel1) && !is.null(d0())) {
      if (input$sel1 == "matrisome gene") {
        nmv <- safe_get_rownames(d0(), assay = DefaultAssay(d0()))
        nmv <- nmv[nmv %in% mobj]
        nmv <- sort(nmv)
        selectizeInput("gene1", "Select a primary feature:",
                       choices = nmv,
                       options = list(create = FALSE, placeholder = 'Type to search...'))
      } else if (input$sel1 == "matrisome signature") {
        tagList(
          radioButtons(
            "list_source1", "Gene list source:",
            choices = c("Curated matrisome list" = "curated", "Custom matrisome list" = "custom"),
            selected = "curated", inline = TRUE
          ),
          conditionalPanel(
            "input.list_source1 == 'curated'",
            selectizeInput("gene1", "Select Matrisome gene list:",
                           choices = list(
                             `Matrisome categories` = c("ECM Glycoproteins", "Collagens", "Proteoglycans",
                                                         "ECM-affiliated Proteins", "ECM Regulators", "Secreted Factors"),
                             `Matrisome subcategories` = c("Perivascular", "Hemostasis", "Elastic fibers",
                                                      "Growth-factor binding"),
                             `Matrisome gene families` = c("Laminins", "Matricellular proteins", "Syndecans", "Glypicans"),
                             `Others` = c("Annexins", "Cathepsins", "CCNs", "Cystatins", "FACITs", "Fibulins",
                                          "Galectins", "Mucins", "Plexins", "Semaphorins")
                           ),
                           options = list(create = FALSE, placeholder = 'Select a gene list...'))
          ),
          conditionalPanel(
            "input.list_source1 == 'custom'",
            selectizeInput(
              "custom_genes1", "Select matrisome genes:", choices = NULL,
              multiple = TRUE,
              options = list(create = FALSE, placeholder = 'Type to search detected matrisome genes...',
                             plugins = list('remove_button'))
            ),
            div(class = "small text-muted", textOutput("custom_count1", inline = TRUE))
          )
        )
      } else {
        # Show empty dropdown when no valid selection
        selectizeInput("gene1", "Select a primary feature:",
                       choices = "",
                       options = list(create = FALSE, placeholder = 'Select a feature type above...'))
      }
    } else {
      # Show empty dropdown when nothing selected or no data
      selectizeInput("gene1", "Select a primary feature:",
                     choices = "",
                     options = list(create = FALSE, placeholder = 'Select a feature type above...'))
    }
  })

  # Conditionally render the SECONDARY feature dropdown
  output$ui2 <- renderUI({
    req(input$sel2, d0())

    if (input$sel2 == "none") {
      return(NULL)
    }

    # For matrisome gene and signature, the client-side approach is fine as lists are small
    if (input$sel2 == "matrisome gene") {
      nmv <- safe_get_rownames(d0(), assay = DefaultAssay(d0()))
      nmv <- nmv[nmv %in% mobj]
      nmv <- sort(nmv)
      return(selectizeInput("gene2", "Select a secondary feature:",
                     choices = nmv,
                     options = list(create = FALSE, placeholder = 'Type to search...')))
    }

    if (input$sel2 == "matrisome signature") {
      return(tagList(
        radioButtons(
          "list_source2", "Gene list source:",
          choices = c("Curated matrisome list" = "curated", "Custom matrisome list" = "custom"),
          selected = "curated", inline = TRUE
        ),
        conditionalPanel(
          "input.list_source2 == 'curated'",
          selectizeInput("gene2", "Select Matrisome gene list:",
                         choices = list(
                           `Matrisome categories` = c("ECM Glycoproteins", "Collagens", "Proteoglycans",
                                                       "ECM-affiliated Proteins", "ECM Regulators", "Secreted Factors"),
                           `Matrisome subcategories` = c("Perivascular", "Hemostasis", "Elastic fibers",
                                                    "Growth-factor binding"),
                           `Matrisome gene families` = c("Laminins", "Matricellular proteins", "Syndecans", "Glypicans"),
                           `Others` = c("Annexins", "Cathepsins", "CCNs", "Cystatins", "FACITs", "Fibulins",
                                        "Galectins", "Mucins", "Plexins", "Semaphorins")
                         ),
                         options = list(create = FALSE, placeholder = 'Select a gene list...'))
        ),
        conditionalPanel(
          "input.list_source2 == 'custom'",
          selectizeInput(
            "custom_genes2", "Select matrisome genes:", choices = NULL,
            multiple = TRUE,
            options = list(create = FALSE, placeholder = 'Type to search detected matrisome genes...',
                           plugins = list('remove_button'))
          ),
          div(class = "small text-muted", textOutput("custom_count2", inline = TRUE))
        )
      ))
    }

    # "Any gene" option: server-side selectizeInput with NULL initial choices
    if (input$sel2 == "any gene") {
      return(selectizeInput("gene2", "Select a secondary feature:",
                     choices = NULL, # Choices are now managed on the server
                     options = list(create = FALSE,
                                    placeholder = 'Type to search any gene...')))
    }
  })

  custom_matrisome_choices <- reactive({
    req(d0())
    available_matrisome_genes(d0(), mobj)
  })

  observe({
    req(identical(input$sel1, "matrisome signature"), identical(input$list_source1, "custom"))
    choices <- custom_matrisome_choices()
    updateSelectizeInput(session, "custom_genes1",
                         choices = choices,
                         options = list(placeholder = if (length(choices) > 0) {
                           "Type to search detected matrisome genes..."
                         } else "No matrisome genes detected in the active assay"),
                         server = TRUE)
  })

  observe({
    req(identical(input$sel2, "matrisome signature"), identical(input$list_source2, "custom"))
    choices <- custom_matrisome_choices()
    updateSelectizeInput(session, "custom_genes2",
                         choices = choices,
                         options = list(placeholder = if (length(choices) > 0) {
                           "Type to search detected matrisome genes..."
                         } else "No matrisome genes detected in the active assay"),
                         server = TRUE)
  })

  output$custom_count1 <- renderText({
    paste(length(input$custom_genes1), "genes selected")
  })

  output$custom_count2 <- renderText({
    paste(length(input$custom_genes2), "genes selected")
  })

  # Observer for server-side "any gene" selectizeInput
  observe({
    req(input$sel2 == "any gene", d0()) # Only run when "any gene" is active

    all_genes <- safe_get_rownames(d0(), assay = DefaultAssay(d0()))

    updateSelectizeInput(
      session,
      "gene2",
      choices = sort(all_genes),
      server = TRUE
    )
  })


  feature_contrast_details <- function(role) {
    results <- rv$analysis_results
    if (is.null(results)) return(NULL)
    id <- if (role == "primary") "customlist1" else "customlist2"
    diagnostic <- results$custom_score_diagnostics[[id]]
    if (is.null(diagnostic) || !isTRUE(diagnostic$flagged)) return(NULL)
    list(id = id, diagnostic = diagnostic)
  }

  any_feature_contrast_warning <- reactive({
    !is.null(feature_contrast_details("primary")) ||
      !is.null(feature_contrast_details("secondary"))
  })

  contrast_badge <- function(role) {
    details <- feature_contrast_details(role)
    if (is.null(details)) return(NULL)
    state <- rv[[paste0(role, "_contrast_state")]]
    if (identical(state, "contrast")) {
      span("Enhanced contrast", class = "feature-contrast-badge active")
    } else {
      span("Limited UCell contrast", class = "feature-contrast-badge warning")
    }
  }

  contrast_notice <- function(role) {
    if (!any_feature_contrast_warning()) return(NULL)
    details <- feature_contrast_details(role)
    if (is.null(details)) {
      return(div(class = "feature-contrast-slot feature-contrast-placeholder"))
    }

    state <- rv[[paste0(role, "_contrast_state")]]
    diagnostic <- details$diagnostic
    prefix <- if (role == "primary") "primary" else "secondary"
    if (identical(state, "contrast")) {
      return(div(
        class = "feature-contrast-slot",
        div(
          class = "feature-contrast-notice active",
          div(
            class = "feature-contrast-copy",
            strong("Relative contrast shown"),
            span("Colours show within-tissue percentiles; hover and statistics retain the original UCell score.")
          ),
          actionButton(paste0("return_", prefix, "_standard"), "Standard scale",
                       class = "btn btn-sm btn-outline-secondary")
        )
      ))
    }

    detail_text <- paste0(
      "Scores are tightly clustered: ", round(100 * diagnostic$upper_fraction),
      "% of spots score at least ", diagnostic$upper_cutoff,
      " (IQR ", format(round(diagnostic$iqr, 3), nsmall = 3), ")."
    )

    if (identical(state, "standard")) {
      return(div(
        class = "feature-contrast-slot",
        div(
          class = "feature-contrast-notice retained",
          div(
            class = "feature-contrast-copy",
            strong("Showing standard UCell scale"),
            span(detail_text)
          ),
          actionButton(paste0("show_", prefix, "_contrast"), "Show relative contrast",
                       class = "btn btn-sm btn-outline-secondary")
        )
      ))
    }

    div(
      class = "feature-contrast-slot",
      div(
        class = "feature-contrast-notice warning",
        div(
          class = "feature-contrast-copy",
          strong("UCell colours have limited contrast"),
          span(paste(detail_text, "Values and statistics will not change."))
        ),
        div(
          class = "feature-contrast-actions",
          actionButton(paste0("keep_", prefix, "_ucell"), "Keep standard",
                       class = "btn btn-sm btn-outline-secondary"),
          actionButton(paste0("show_", prefix, "_contrast"), "Show relative contrast",
                       class = "btn btn-sm btn-primary")
        )
      )
    )
  }

  output$primary_contrast_badge <- renderUI(contrast_badge("primary"))
  output$secondary_contrast_badge <- renderUI(contrast_badge("secondary"))
  output$primary_contrast_notice <- renderUI(contrast_notice("primary"))
  output$secondary_contrast_notice <- renderUI(contrast_notice("secondary"))

  observeEvent(input$keep_primary_ucell, rv$primary_contrast_state <- "standard")
  observeEvent(input$keep_secondary_ucell, rv$secondary_contrast_state <- "standard")
  observeEvent(input$show_primary_contrast, rv$primary_contrast_state <- "contrast")
  observeEvent(input$show_secondary_contrast, rv$secondary_contrast_state <- "contrast")
  observeEvent(input$return_primary_standard, rv$primary_contrast_state <- "standard")
  observeEvent(input$return_secondary_standard, rv$secondary_contrast_state <- "standard")

  # Plot outputs
  output$primary_feature_plot <- renderPlotly({
    req(rv$feature_analysis_done) # Analysis complete check
    analysis_results <- rv$analysis_results
    req(analysis_results)

    # Get full, unfiltered data
    full_plot_data <- analysis_results$data
    selected_types <- selectedCellTypes()
    annotation_vector <- cleaned_annotation()
    annotation_vector <- annotation_vector[rownames(full_plot_data)]

    # Create a logical index to filter both data and annotation consistently
    row_indices_to_keep <- annotation_vector %in% selected_types
    filtered_plot_data <- full_plot_data[row_indices_to_keep, ]
    filtered_annotation_vector <- annotation_vector[row_indices_to_keep]
    show_contrast <- identical(rv$primary_contrast_state, "contrast")
    colour_values <- if (show_contrast) {
      analysis_results$custom_contrast_scores$customlist1
    } else NULL

    # Pass the matched filtered data and filtered vector
    create_feature_plot(
      filtered_plot_data,
      "feature1",
      analysis_results$gene1,
      filtered_annotation_vector,
      colour_values = colour_values,
      relative_contrast = show_contrast
    )
  })

  output$secondary_feature_plot <- renderPlotly({
    req(rv$feature_analysis_done) # Analysis complete check
    analysis_results <- rv$analysis_results
    req(analysis_results)

    if (analysis_results$sel2 != "none" && !is.null(analysis_results$gene2) && analysis_results$gene2 != "") {

      # Get full, unfiltered data
      full_plot_data <- analysis_results$data
      selected_types <- selectedCellTypes()
      annotation_vector <- cleaned_annotation()
      annotation_vector <- annotation_vector[rownames(full_plot_data)]

      # Create a logical index to filter both data and annotation consistently
      row_indices_to_keep <- annotation_vector %in% selected_types
      filtered_plot_data <- full_plot_data[row_indices_to_keep, ]
      filtered_annotation_vector <- annotation_vector[row_indices_to_keep]
      show_contrast <- identical(rv$secondary_contrast_state, "contrast")
      colour_values <- if (show_contrast) {
        analysis_results$custom_contrast_scores$customlist2
      } else NULL

      # Pass the matched filtered data and filtered vector
      create_feature_plot(
        filtered_plot_data,
        "feature2",
        analysis_results$gene2,
        filtered_annotation_vector,
        colour_values = colour_values,
        relative_contrast = show_contrast
      )
    } else {
      # Return an informative empty plot if no secondary feature was part of the analysis
      plot_ly() %>%
        add_annotations(
          text = "No secondary feature selected",
          showarrow = FALSE,
          font = list(size = 16)
        )
    }
  })

  # Feature autocorrelation statistics (Moran's I)
  feature_autocorrelation_stats <- eventReactive(rv$analysis_results, {

    # Ensure analysis results exist before proceeding
    if (is.null(rv$analysis_results) || is.null(rv$analysis_results$data)) return(list(primary = NULL, secondary = NULL))

    results <- list(primary = NULL, secondary = NULL)
    data <- rv$analysis_results$data

    if (nrow(data) < 3) return(results)

    # Use pixel coords if array coords are collinear (SPE conversion edge case)
    coords <- data[, c("row", "col")]

    # Check if row and col are identical sequential integers (collinear fallback case)
    if (all(coords$row == coords$col)) {
      # Fall back to pixel coordinates which are always valid
      coords <- data[, c("imagerow", "imagecol")]
    }

    weight_matrix <- MERINGUE::getSpatialNeighbors(coords)

    if (var(data$feature1, na.rm = TRUE) != 0) {
      tryCatch({
        primary_feature_vector <- setNames(data$feature1, rownames(data))
        primary_test_result <- MERINGUE::moranTest(primary_feature_vector, weight_matrix)
        results$primary <- list(morans_I = primary_test_result["observed"], p_value = primary_test_result["p.value"])

      }, error = function(e) {})
    }

    if (!is.null(rv$analysis_results$sel2) && rv$analysis_results$sel2 != "none" && !is.null(data$feature2) && var(data$feature2, na.rm = TRUE) != 0) {
      tryCatch({
        secondary_feature_vector <- setNames(data$feature2, rownames(data))
        secondary_test_result <- MERINGUE::moranTest(secondary_feature_vector, weight_matrix)
        results$secondary <- list(morans_I = secondary_test_result["observed"], p_value = secondary_test_result["p.value"])

      }, error = function(e) {})
    }

    return(results)
  })

  # ECM Niche Autocorrelation reactive
  ecm_autocorrelation_stats <- eventReactive(rv$matrisome_results, {
    req(rv$matrisome_results)
    seurat_obj <- rv$matrisome_results$glycoprotein$seurat_object # Get object with scores
    req(seurat_obj)

    stats_list <- list()
    ecm_sig_cols <- c("Interstitial_UCell", "Basement_UCell")

    # Check if there's enough data
    coords <- GetTissueCoordinates(seurat_obj)
    if (nrow(coords) < 3) return(NULL)

    # Handle VisiumV2 coordinate column names (x/y instead of imagecol/imagerow)
    if(!all(c("imagerow","imagecol")%in%colnames(coords))){
      coords$imagerow <- coords$x
      coords$imagecol <- coords$y
    }
    coords <- coords[,colnames(coords)%in%c("imagerow","imagecol")]

    weight_matrix <- MERINGUE::getSpatialNeighbors(coords)

    for (col_name in ecm_sig_cols) {
      prefix <- tolower(sub("_UCell", "", col_name))
      if (col_name %in% colnames(seurat_obj@meta.data)) {
        scores <- seurat_obj@meta.data[[col_name]]
        if (var(scores, na.rm = TRUE) > 0) {
          tryCatch({
            test_result <- MERINGUE::moranTest(setNames(scores, rownames(seurat_obj@meta.data)), weight_matrix)
            stats_list[[prefix]] <- list(morans_I = test_result["observed"], p_value = test_result["p.value"])
          }, error = function(e) {
            stats_list[[prefix]] <- NULL
          })
        }
      }
    }
    return(stats_list)
  })

  # RENDER UI for the PRIMARY feature card footer
  output$primary_feature_autocorrelation_footer <- renderUI({
    req(rv$analysis_results) # Wait for data
    stats <- feature_autocorrelation_stats()

    if (!is.null(stats$primary)) {
      morans_i <- stats$primary$morans_I
      p_value <- stats$primary$p_value

      if (p_value >= 0.05) {
        badge_text <- "RANDOM"; badge_class <- "badge-random"
      } else if (morans_i > 0) {
        badge_text <- "CLUSTERED"; badge_class <- "badge-clustered"
      } else {
        badge_text <- "DISPERSED"; badge_class <- "badge-dispersed"
      }

      tags$div(
        # Use justify-content-between to push items to both ends
        class = "footer-content d-flex justify-content-between align-items-center",
        # Left-aligned label
        tags$span(
          class = "footer-label",
          bsicons::bs_icon("bar-chart-line-fill"),
          "Spatial Autocorrelation"
        ),
        # Right-aligned stats
        div(
          class = "footer-stats d-flex align-items-center",
          style = "gap: 0.75rem;",
          tags$span(paste0("Moran's I: ", format(round(morans_i, 2), nsmall = 2))),
          tags$span(class = "p-value", paste0("p = ", format.pval(p_value, digits = 2, eps = 0.001))),
          tags$span(class = paste("badge", badge_class), badge_text)
        )
      )
    } else {
      tags$div(class = "footer-content", "Spatial autocorrelation not available")
    }
  })

  # RENDER UI for the SECONDARY feature card footer
  output$secondary_feature_autocorrelation_footer <- renderUI({
    req(rv$analysis_results) # Wait for data
    stats <- feature_autocorrelation_stats()

    if (!is.null(stats$secondary)) {
      morans_i <- stats$secondary$morans_I
      p_value <- stats$secondary$p_value

      if (p_value >= 0.05) {
        badge_text <- "RANDOM"; badge_class <- "badge-random"
      } else if (morans_i > 0) {
        badge_text <- "CLUSTERED"; badge_class <- "badge-clustered"
      } else {
        badge_text <- "DISPERSED"; badge_class <- "badge-dispersed"
      }

      tags$div(
        class = "footer-content d-flex justify-content-between align-items-center",
        tags$span(
          class = "footer-label",
          bsicons::bs_icon("bar-chart-line-fill"),
          "Spatial Autocorrelation"
        ),
        div(
          class = "footer-stats d-flex align-items-center",
          style = "gap: 0.75rem;",
          tags$span(paste0("Moran's I: ", format(round(morans_i, 2), nsmall = 2))),
          tags$span(class = "p-value", paste0("p = ", format.pval(p_value, digits = 2, eps = 0.001))),
          tags$span(class = paste("badge", badge_class), badge_text)
        )
      )
    } else {
      NULL
    }
  })

  # p3a function is now in helpers.R

  # p3b function is now in helpers.R

  output$primary_expression_plot <- renderPlotly({
    req(rv$feature_analysis_done) # Analysis complete check
    analysis_results <- rv$analysis_results
    req(analysis_results)
    # Pass both the data and the correct gene name to the function
    p3a(
      data_df = analysis_results$data,
      feature_name = analysis_results$gene1,
      stable_color_map = color_map(),
      selected_types = selectedCellTypes(),
      cleaned_annot_vector = cleaned_annotation(),
      use_violin = input$violin_plot,
      independent_y = input$independent_y
    )
  })

  output$secondary_expression_plot <- renderPlotly({
    req(rv$feature_analysis_done) # Analysis complete check
    analysis_results <- rv$analysis_results
    req(analysis_results)

    # Check based on analysis results (not live input values) for consistency
    if (analysis_results$sel2 != "none" && !is.null(analysis_results$gene2) && analysis_results$gene2 != "") {
      # Pass the data and the gene name used in the analysis
      p3b(
        data_df = analysis_results$data,
        feature_name = analysis_results$gene2,
        stable_color_map = color_map(),
        selected_types = selectedCellTypes(),
        cleaned_annot_vector = cleaned_annotation(),
        use_violin = input$violin_plot,
        independent_y = input$independent_y
      )
    } else {
      # If no secondary feature was run, return an empty plot
      return(NULL)
    }
  })

  # Co-expression Plot
  output$coexpression_plot <- renderPlot({
    req(rv$feature_analysis_done) # Analysis complete check
    analysis_results <- rv$analysis_results
    req(analysis_results)
    req(input$coex_palette, input$coex_alpha) # Ensure the new inputs are available

    # Check if a secondary feature was selected for the analysis
    if (analysis_results$sel2 == "none" || is.null(analysis_results$gene2) || analysis_results$gene2 == "") {
      return(
        ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = "A secondary feature is required", size = 5) +
          theme_void()
      )
    }

    # Create a temporary Seurat object and add feature values to metadata
    temp_seurat_obj <- d0()
    cell_order <- intersect(Cells(temp_seurat_obj), rownames(analysis_results$data))
    req(length(cell_order) > 0)
    temp_seurat_obj <- temp_seurat_obj[, cell_order]
    temp_seurat_obj$feature_blend_1 <- analysis_results$data[cell_order, "feature1"]
    temp_seurat_obj$feature_blend_2 <- analysis_results$data[cell_order, "feature2"]

    colors_to_use <- coex_blend_colors(input$coex_palette)

    sfp_args <- list(
      pt.size.factor = effective_pt_size(input$coex_pt_size),
      alpha = input$coex_alpha
    )

    # 5. Generate the full plot object using the temporary object and the new feature names
    full_plot_object <- SpatialFeaturePlotBlend(
      object = temp_seurat_obj,
      features = c("feature_blend_1", "feature_blend_2"),
      # Pass original names for the legend labels
      feature_1_alt_name = analysis_results$gene1,
      feature_2_alt_name = analysis_results$gene2,
      combine = TRUE,
      bottom_left = colors_to_use$bottom_left,
      bottom_right = colors_to_use$bottom_right,
      top_left = colors_to_use$top_left,
      top_right = colors_to_use$top_right,
      sfp_extra_arguments = sfp_args
    )

    # 6. Extract, modify, and combine the desired plots
    inner_plots <- full_plot_object$patches$plots[[1]]$patches$plots
    blended_plot <- inner_plots[[3]] + theme(plot.title = element_blank()) # Remove title
    legend_plot <- inner_plots[[4]]

    final_plot <- wrap_plots(blended_plot, legend_plot, nrow = 1, widths = c(0.75, 0.25))

    # 7. Return the final plot
    final_plot
  })

  crosscorrelation_stats <- eventReactive(rv$analysis_results, {

    # Ensure analysis results exist and a secondary feature exists before proceeding
    if (is.null(rv$analysis_results) || is.null(rv$analysis_results$sel2) || rv$analysis_results$sel2 == "none") return(NULL)

    results <- list(pearson = NULL, spatial = NULL)
    data <- rv$analysis_results$data

    # Check for sufficient data and variance in both features
    if (nrow(data) < 3 || var(data$feature1, na.rm = TRUE) == 0 || var(data$feature2, na.rm = TRUE) == 0) {
      return(results)
    }

    # --- Pearson Correlation ---
    tryCatch({
      pearson_test <- cor.test(data$feature1, data$feature2)
      results$pearson <- list(cor = pearson_test$estimate, pval = pearson_test$p.value)
    }, error = function(e) {})

    # --- Spatial Cross-Correlation ---
    tryCatch({
      coords <- data[, c("row", "col")]
      # Check for collinear coordinates (VisiumV2 fallback case)
      if (all(coords$row == coords$col)) {
        coords <- data[, c("imagerow", "imagecol")]
      }
      weight_matrix <- MERINGUE::getSpatialNeighbors(coords)

      # Use named vectors to ensure correct alignment
      primary_vec <- setNames(data$feature1, rownames(data))
      secondary_vec <- setNames(data$feature2, rownames(data))

      spatial_test <- spatial_cross_cor_test(primary_vec, secondary_vec, weight_matrix, n = 999)
      results$spatial <- list(cor = spatial_test$cor, pval = spatial_test$pval)

    }, error = function(e) {})

    return(results)
  })

  # RENDER HTML for the co-expression card footer
  output$coexpression_crosscor_footer <- renderUI({
    req(rv$analysis_results)
    stats <- crosscorrelation_stats()

    if (is.null(stats) || (is.null(stats$pearson) && is.null(stats$spatial))) {
      return(tags$div(class = "footer-content", "Cross-correlation not available"))
    }

    tagList(
      if (!is.null(stats$pearson)) {
        tags$div(
          class = "footer-content",
          tags$span(class = "footer-label", "Pearson Correlation"),
          tags$div(
            class = "footer-stats",
            tags$span(paste0("r = ", format(round(stats$pearson$cor, 2), nsmall = 2))),
            tags$span(class = "p-value", paste0("p = ", format.pval(stats$pearson$pval, digits = 2, eps = 0.001)))
          )
        )
      },
      if (!is.null(stats$spatial)) {
        tags$div(
          class = "footer-content",
          tags$span(class = "footer-label", "Spatial Cross-Correlation"),
          tags$div(
            class = "footer-stats",
            tags$span(format(round(stats$spatial$cor, 2), nsmall = 2)),
            if (!is.null(stats$spatial$pval) && !is.na(stats$spatial$pval)) {
              tags$span(class = "p-value", paste0("p = ", format.pval(stats$spatial$pval, digits = 2, eps = 1e-3)))
            }
          )
        )
      }
    )
  })

  output$lisa_clustering_plot <- renderPlot({
    req(rv$feature_analysis_done, rv$analysis_results)
    analysis_data <- rv$analysis_results$data

    plot_lisa_clustering(
      k = analysis_data,
      selectedCellTypes = selectedCellTypes(),
      ann_col = ann(),
      f1_name = rv$analysis_results$gene1,
      f2_name = rv$analysis_results$gene2
    )
  })

  # --- [ MATRISOME PLOT RENDERING LOGIC ] ---
  # This section programmatically creates all static/interactive renderers
  # for the Matrisome tab based on the conditionalPanel UI.

  # Helper function to generate the server-side renderers for a single plot
  # It uses local() to correctly handle variables within the loop.
  create_plot_renderers <- function(matrisome_group, plot_name, pt_size_input_id) {
    local({
      current_matrisome_group <- matrisome_group

      # Helper function: Creates a standardized ggplot for unavailable data
      render_not_available_plot <- function(message = "Data not available for this category") {
        ggplot() +
          annotate("text", x = 1, y = 1, label = message, size = 5, color = "grey50") +
          theme_void()
      }

      # Creates a standardized plotly object for when data is unavailable.
      render_not_available_plotly <- function(message = "Data not available for this category") {
        plot_ly() %>%
          layout(
            xaxis = list(visible = FALSE),
            yaxis = list(visible = FALSE),
            annotations = list(
              text = message, showarrow = FALSE, font = list(size = 14, color = "grey50"),
              xref = "paper", yref = "paper", x = 0.5, y = 0.5
            )
          )
      }
      # Master plot renderer with validation and sanitization
      render_static_matrisome_plot <- function(feature_type) {
        renderPlot({
          # Check analysis flag before tryCatch to stop silently on app launch
          req(rv$matrisome_analysis_done)

          tryCatch({
            # The other reqs can stay inside for post-run data validation.
            req(rv$matrisome_results, input[[pt_size_input_id]])

            # Validate data existence
            plot_data <- rv$matrisome_results[[current_matrisome_group]]
            if (is.null(plot_data) || is.null(plot_data$seurat_object)) {
              return(render_not_available_plot())
            }

            full_seurat_obj <- plot_data$seurat_object
            feature_col <- if (feature_type == "distribution") plot_data$robust_feature else plot_data$log_scaled_feature

            if (!feature_col %in% colnames(full_seurat_obj@meta.data)) {
              return(render_not_available_plot())
            }

            # Check data sufficiency
            feature_values <- full_seurat_obj@meta.data[[feature_col]]
            valid_indices <- !is.na(feature_values)
            non_na_count <- sum(valid_indices)
            total_count <- length(feature_values)

            # A plot is only meaningful if it has a reasonable number of data points.
            if (non_na_count < 20 || non_na_count < (total_count * 0.01)) {
              return(render_not_available_plot("Insufficient data for visualization"))
            }

            # Subset to spots with valid data only
            valid_spot_names <- colnames(full_seurat_obj)[valid_indices]
            seurat_obj_for_plotting <- full_seurat_obj[, valid_spot_names]

            # Determine the correct legend title based on the plot type
            legend_title <- if (feature_type == "distribution") {
              "Expression Z-Score"
            } else {
              "Expression Intensity (Log)"
            }

            SpatialFeaturePlot(seurat_obj_for_plotting, features = feature_col, pt.size.factor = effective_pt_size(input[[pt_size_input_id]])) +
              theme(legend.position = "right") +
              labs(fill = legend_title)

          }, error = function(e) {
            # Final Safety Net: Catches any other unexpected errors.
            render_not_available_plot("An error occurred during plotting")
          })
        }, height = 400)
      }

      # Distribution Plot Rendering
      static_id <- paste0(plot_name, "_dist_static")
      plotly_id <- paste0(plot_name, "_dist_plotly")

      # Call the new robust rendering function
      output[[static_id]] <- render_static_matrisome_plot("distribution")

      # On-demand plotly: paired with suspendWhenHidden below so only the
      # visible category materializes a plotly object.
      output[[plotly_id]] <- renderPlotly({
        req(rv$matrisome_analysis_done, rv$matrisome_results)
        plot_data <- rv$matrisome_results[[current_matrisome_group]]
        if (is.null(plot_data) || is.null(plot_data$df) || is.null(plot_data$robust_colors)) {
          render_not_available_plotly()
        } else {
          build_matrisome_plotly(plot_data, kind = "spatial")
        }
      })

      # --- Hotspot Plot Rendering ---
      static_id_hotspot <- paste0(plot_name, "_hotspot_static")
      plotly_id_hotspot <- paste0(plot_name, "_hotspot_plotly")

      # Call the new robust rendering function
      output[[static_id_hotspot]] <- render_static_matrisome_plot("hotspot")

      output[[plotly_id_hotspot]] <- renderPlotly({
        req(rv$matrisome_analysis_done, rv$matrisome_results)
        plot_data <- rv$matrisome_results[[current_matrisome_group]]
        if (is.null(plot_data) || is.null(plot_data$df) || is.null(plot_data$hotspot_colors)) {
          render_not_available_plotly()
        } else {
          build_matrisome_plotly(plot_data, kind = "hotspot")
        }
      })

      outputOptions(output, static_id, suspendWhenHidden = TRUE)
      outputOptions(output, plotly_id, suspendWhenHidden = TRUE)
      outputOptions(output, static_id_hotspot, suspendWhenHidden = TRUE)
      outputOptions(output, plotly_id_hotspot, suspendWhenHidden = TRUE)
    })
  }

  # Loop to generate Main Division Plots
  main_divisions <- c(
    "glycoprotein" = "glycoprotein", "collagen" = "collagen", "proteoglycan" = "proteoglycan",
    "regulator" = "regulator", "secreted_factor" = "secreted_factor", "affliated_protein" = "affliated_protein"
  )
  for (group in names(main_divisions)) {
    create_plot_renderers(group, main_divisions[group], pt_size_input_id = "main_pt_size")
  }

  # --- Loop to generate Matrisome Subcategories Plots ---
  ecm_subcategories_plots <- c(
    "perivascular" = "perivascular",
    "hemostasis" = "hemostasis",
    "elastic_fibers" = "elastic_fibers",
    "gf_binding" = "gf_binding"
  )
  for (group in names(ecm_subcategories_plots)) {
    create_plot_renderers(group, ecm_subcategories_plots[group], pt_size_input_id = "subcategories_pt_size")
  }

  # --- Loop to generate Matrisome Gene Families Plots ---
  ecm_gene_families_plots <- c(
    "laminin" = "laminin",
    "matricellular" = "matricellular",
    "syndecan" = "syndecan",
    "glypican" = "glypican"
  )
  for (group in names(ecm_gene_families_plots)) {
    create_plot_renderers(group, ecm_gene_families_plots[group], pt_size_input_id = "families_pt_size")
  }

  # Helper function to generate the server-side renderers for autocorrelation footers
  create_footer_renderer <- function(plot_name) {
    local({
      current_plot_name <- plot_name
      output_id <- paste0(current_plot_name, "_autocorrelation_footer")

      output[[output_id]] <- renderUI({
        # Ensure analysis is complete
        req(rv$matrisome_analysis_done)

        # Check 1: Does the matrisome data exist for this plot at all?
        if (is.null(rv$matrisome_results[[current_plot_name]])) {
          return(tags$div(class = "footer-content", "Autocorrelation not applicable (no data)."))
        }

        autocorr_data <- rv$matrisome_results$autocorrelation

        # Check 2: Was autocorrelation data calculated?
        if (is.null(autocorr_data)) {
          return(tags$div(class = "footer-content", "Autocorrelation calculation failed."))
        }

        # Get the correct category name for lookup
        category_id <- paste0(current_plot_name, "_robust_score")
        result <- autocorr_data[autocorr_data$category == category_id, , drop = FALSE]

        # Check 3: Is there a result for this specific category?
        if (nrow(result) == 0) {
          return(tags$div(class = "footer-content", "Autocorrelation not available for this category (uniform expression)."))
        }

        # Extract values
        morans_i <- result$morans_I[1]
        p_value <- result$p_value[1]

        # Check for NA p-value before conditional statements
        if (is.na(p_value)) {
          badge_text <- "N/A"
          badge_class <- "badge-random"
        } else if (p_value >= 0.05) {
          badge_text <- "RANDOM"
          badge_class <- "badge-random"
        } else if (morans_i > 0) {
          badge_text <- "CLUSTERED"
          badge_class <- "badge-clustered"
        } else {
          badge_text <- "DISPERSED"
          badge_class <- "badge-dispersed"
        }

        # Build the final UI
        tags$div(
          class = "footer-content d-flex justify-content-between align-items-center",
          tags$span(
            class = "footer-label",
            bsicons::bs_icon("bar-chart-line-fill"),
            "Spatial Autocorrelation"
          ),
          tags$div(
            class = "footer-stats",
            tags$span(paste0("Moran's I: ", format(round(morans_i, 2), nsmall = 2))),
            tags$span(
              class = "p-value",
              paste0("p = ", format.pval(p_value, digits = 2, eps = 0.001))
            ),
            tags$span(class = paste("badge", badge_class), badge_text)
          )
        )
      })
    })
  }

  # --- Loop to generate Main Division Footers ---
  for (group in names(main_divisions)) {
    create_footer_renderer(main_divisions[group])
  }

  # --- Loop to generate Matrisome Subcategories Footers ---
  for (group in names(ecm_subcategories_plots)) {
    create_footer_renderer(ecm_subcategories_plots[group])
  }

  # --- Loop to generate Matrisome Gene Families Footers ---
  for (group in names(ecm_gene_families_plots)) {
    create_footer_renderer(ecm_gene_families_plots[group])
  }

  # --- [ ECM NICHES RENDERING LOGIC ] ---

  # Helper functions to create "not available" plots
  render_not_available_plot <- function(message = "Data not available for this signature") {
    ggplot() +
      annotate("text", x = 1, y = 1, label = message, size = 5, color = "grey50") +
      theme_void()
  }

  render_not_available_plotly <- function(message = "Data not available for this signature") {
    plot_ly() %>%
      layout(
        xaxis = list(visible = FALSE),
        yaxis = list(visible = FALSE),
        annotations = list(
          text = message, showarrow = FALSE, font = list(size = 14, color = "grey50"),
          xref = "paper", yref = "paper", x = 0.5, y = 0.5
        )
      )
  }

  # Programmatic Rendering for Signature Cards
  ecm_signatures <- list(
    Interstitial = "interstitial",
    Basement = "basement"
  )

  for (sig_name in names(ecm_signatures)) {
    local({
      # Necessary to ensure loop variables are correctly captured
      current_sig_name <- sig_name
      prefix <- ecm_signatures[[current_sig_name]]

      # Define column names
      dist_col <- paste0(current_sig_name, "_UCell")
      hotspot_col <- paste0(current_sig_name, "_Hotspot_Score")

      # Helper to get the final seurat object
      get_final_obj <- reactive({
        req(rv$matrisome_results)
        # Safely find the seurat object from any available result
        for (item in rv$matrisome_results) {
          if (!is.null(item) && !is.null(item$seurat_object)) return(item$seurat_object)
        }
        return(NULL)
      })

      # Static Distribution Plot
      output[[paste0(prefix, "_dist_static")]] <- renderPlot({
        req(rv$matrisome_analysis_done, get_final_obj())
        obj <- get_final_obj()
        if (!dist_col %in% colnames(obj@meta.data)) {
          return(render_not_available_plot())
        }
        SpatialFeaturePlot(obj, features = dist_col, pt.size.factor = effective_pt_size(input$ecm_pt_size)) +
          theme(legend.position = "right") + labs(fill = "Niche Score")
      })

      # Interactive Distribution Plot
      output[[paste0(prefix, "_dist_plotly")]] <- renderPlotly({
        req(rv$matrisome_analysis_done, get_final_obj())
        obj <- get_final_obj()
        if (!dist_col %in% colnames(obj@meta.data)) {
          return(render_not_available_plotly())
        }

        # Get the required annotation column name
        annotation_col <- ann()
        # Fetch BOTH the signature score and the annotation data
        plot_data <- FetchData(obj, vars = c(dist_col, annotation_col))
        # Combine with coordinates to create the final data frame
        coords <- align_spatial_coordinates(obj, cells = rownames(plot_data))
        plot_data <- plot_data[rownames(coords), , drop = FALSE]
        df <- cbind(coords, plot_data)
        create_feature_plot(df, dist_col, current_sig_name, cleaned_annotation())
      })

      # Static Hotspot Plot
      output[[paste0(prefix, "_hotspot_static")]] <- renderPlot({
        req(rv$matrisome_analysis_done, get_final_obj())
        obj <- get_final_obj()
        if (!hotspot_col %in% colnames(obj@meta.data)) {
          return(render_not_available_plot())
        }
        SpatialFeaturePlot(obj, features = hotspot_col, pt.size.factor = effective_pt_size(input$ecm_pt_size)) +
          theme(legend.position = "right") + labs(fill = "Niche Score Intensity (Log)")
      })

      # Interactive Hotspot Plot
      output[[paste0(prefix, "_hotspot_plotly")]] <- renderPlotly({
        req(rv$matrisome_analysis_done, get_final_obj())
        obj <- get_final_obj()
        if (!hotspot_col %in% colnames(obj@meta.data)) {
          return(render_not_available_plotly())
        }
        plot_data <- FetchData(obj, vars = hotspot_col)
        coords <- align_spatial_coordinates(obj, cells = rownames(plot_data))
        plot_data <- plot_data[rownames(coords), , drop = FALSE]
        df <- cbind(coords, plot_data)
        create_feature_plot(df, hotspot_col, paste(current_sig_name, "Hotspot"), cleaned_annotation())
      })

      # Autocorrelation Footer (using renderText as requested)
      output[[paste0(prefix, "_autocorrelation_footer")]] <- renderText({
        req(ecm_autocorrelation_stats())
        stats <- ecm_autocorrelation_stats()[[prefix]]

        if (is.null(stats)) {
          return(as.character(tags$div(class = "footer-content", "Autocorrelation not available.")))
        }

        morans_i <- stats$morans_I
        p_value <- stats$p_value

        if (is.na(p_value) || p_value >= 0.05) {
          badge_text <- "RANDOM"; badge_class <- "badge-random"
        } else if (morans_i > 0) {
          badge_text <- "CLUSTERED"; badge_class <- "badge-clustered"
        } else {
          badge_text <- "DISPERSED"; badge_class <- "badge-dispersed"
        }

        # Build HTML using tags and convert to string for renderText
        html_content <- tags$div(
          class = "footer-content d-flex justify-content-between align-items-center",
          tags$span(class = "footer-label", bsicons::bs_icon("bar-chart-line-fill"), "Spatial Autocorrelation"),
          tags$div(
            class = "footer-stats d-flex align-items-center",
            style = "gap: 0.75rem;",
            tags$span(paste0("Moran's I: ", format(round(morans_i, 2), nsmall = 2))),
            tags$span(class = "p-value", paste0("p = ", format.pval(p_value, digits = 2, eps = 0.001))),
            tags$span(class = paste("badge", badge_class), badge_text)
          )
        )
        as.character(html_content)
      })
    })
  }

  # Prepare data for the ECM signature distribution plot
  signature_distribution_data <- eventReactive(rv$matrisome_results, {
    req(rv$matrisome_results, ann())

    # Safely get the final seurat object with all scores
    seurat_obj <- NULL
    for (item in rv$matrisome_results) {
      if (!is.null(item) && !is.null(item$seurat_object)) {
        seurat_obj <- item$seurat_object
        break
      }
    }
    req(seurat_obj)

    # Define needed columns and check if they exist
    annot_col <- ann()
    sig_cols <- c("Interstitial_UCell", "Basement_UCell")

    # Filter to only include signature columns that actually exist
    available_sig_cols <- sig_cols[sig_cols %in% colnames(seurat_obj@meta.data)]
    req(annot_col %in% colnames(seurat_obj@meta.data), length(available_sig_cols) > 0)

    # Process the data using dplyr and tidyr
    plot_df <- seurat_obj@meta.data %>%
      # 1. Select only the necessary columns (annotation + available signatures)
      select(all_of(c(annot_col, available_sig_cols))) %>%
      # 2. Reshape from wide to long format
      tidyr::pivot_longer(
        cols = all_of(available_sig_cols),
        names_to = "signature_name",
        values_to = "score"
      ) %>%
      # Clean up the signature names for the legend (use full names)
      mutate(signature_name = case_when(
        signature_name == "Interstitial_UCell" ~ "Interstitial ECM",
        signature_name == "Basement_UCell" ~ "Basement membrane",
        TRUE ~ sub("_UCell", "", signature_name)
      ))

    # Rename the dynamic annotation column to a static name
    names(plot_df)[names(plot_df) == annot_col] <- "annotation_group"

    plot_df %>%
      # 3. Group by annotation and signature
      group_by(annotation_group, signature_name) %>%
      # 4. Calculate the average score for each group
      summarise(mean_score = mean(score, na.rm = TRUE), .groups = 'drop') %>%
      # Sort within each group
      group_by(annotation_group) %>%
      arrange(desc(mean_score), .by_group = TRUE) %>%
      ungroup() %>%
      mutate(signature_name = forcats::fct_inorder(signature_name))
  })

  # Render the ECM signature distribution bar plot
  output$ecm_signature_distribution_plot <- renderPlot({
    # Only render after matrisome analysis is complete
    req(rv$matrisome_analysis_done)

    plot_data <- signature_distribution_data()
    req(plot_data)

    ggplot(plot_data, aes(x = annotation_group, y = mean_score, fill = signature_name)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
      scale_fill_manual(
        values = c(
          "Basement membrane" = "#6a3d9aff",
          "Interstitial ECM" = "#2b9e2bff"
        )
      ) +
      labs(
        x = NULL,
        y = "Average Niche Score",
        fill = NULL
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_blank(),
        legend.position = "top"
      )
  })

  # ========================================================================== #
  #       MATRISOME-PAIR ANALYSIS: CLUSTER SELECTOR DROPDOWN                  #
  # ========================================================================== #
  # Dynamically populates the niche selection dropdown for the volcano plot.
  # Only shows ECM domain annotations (excludes "not.assigned" spots).
  # ========================================================================== #

  output$matrisome_pair_cluster_selector_ui <- renderUI({
    req(rv$data_is_loaded)

    # Validate that ECM annotation column exists (added during preprocessing or manually)
    validate(need("ecm_domain_annotation" %in% colnames(d0()@meta.data),
                  "ECM Domain Annotation not found in data."))

    # Get unique ECM niche annotations
    cluster_choices <- unique(as.character(d0()@meta.data[["ecm_domain_annotation"]]))

    # Filter out unassigned spots (not meaningful for enrichment analysis)
    # Handle both dot and underscore variants
    cluster_choices <- cluster_choices[!grepl("not.assigned|not_assigned", cluster_choices, ignore.case = TRUE)]

    # Create dropdown (sorted alphabetically for easy navigation)
    selectInput(
      "matrisome_pair_cluster_select",
      "Focus volcano plot on:",
      choices = sort(cluster_choices),
      selected = cluster_choices[1]  # Default to first alphabetically
    )
  })

  # ========================================================================== #
  #                SPATIAL MATRISOME-PAIR CO-EXPRESSION                       #
  # ========================================================================== #
  # Calculates enrichment of MatriComDB pair co-expression within ECM niches.
  # Results are stored for the heatmap and volcano plot.
  # ========================================================================== #

  observeEvent(input$run_matrisome_pair_analysis, {
    req(d0(), "ecm_domain_annotation" %in% colnames(d0()@meta.data), matrisome_pairs)

    # Disable the button to prevent re-clicks
    shinyjs::disable("run_matrisome_pair_analysis")

    # Show progress modal (shinybusy)
    show_modal_progress_line(text = "Starting matrisome pair analysis...")

    # Use tryCatch to ensure the modal is always removed, even on error
    tryCatch({

      update_modal_progress(value = 0.1, text = "Preparing data: building spatial neighbour graph...")

      # Scope the correction to MatriComDB pair genes that are present in the
      # assay used for this analysis.
      pair_assay <- if ("SCT" %in% SeuratObject::Assays(d0())) "SCT" else DefaultAssay(d0())
      seurat_obj_local <- leak_apply(
        d0(),
        intersect(
          unique(c(matrisome_pairs$Gene1, matrisome_pairs$Gene2)),
          safe_get_rownames(d0(), assay = pair_assay)
        ),
        context = "matrisome pair analysis"
      )

      # ECM domain annotation column name (added by preprocessing or ScType)
      grouping_var <- "ecm_domain_annotation"

      # Standardize annotation names: replace spaces and hyphens with underscores
      # This ensures consistent matching across different data sources
      seurat_obj_local[[grouping_var]] <- gsub(" |\\-", "_", seurat_obj_local[[grouping_var, drop = TRUE]])

      # Build spatial adjacency matrix from tissue coordinates
      pair_coords <- GetTissueCoordinates(seurat_obj_local)
      # Handle VisiumV2 coordinate column names (x/y instead of imagecol/imagerow)
      if(!all(c("imagerow","imagecol")%in%colnames(pair_coords))){
        pair_coords$imagerow <- pair_coords$x
        pair_coords$imagecol <- pair_coords$y
      }
      pair_coords <- pair_coords[,colnames(pair_coords)%in%c("imagerow","imagecol")]
      adj_matrix <- getSpatialNeighbors(pair_coords)

      seurat_with_scores <- compute_spatial_matrisome_pair_scores_single_core(
        seurat_obj = seurat_obj_local,
        pair_db    = matrisome_pairs,
        adj_matrix = adj_matrix,
        assay      = pair_assay,
        update_progress = update_modal_progress  # Pass shinybusy function for progress updates
      )

      rm(seurat_obj_local, adj_matrix, pair_coords)
      gc(verbose = FALSE)

      update_modal_progress(value = 0.85, text = "Calculating enrichment: comparing niches...")
      interaction_stats <- find_matrisome_pair_enrichment(
        seurat_with_scores,
        "MATRISOMEPAIR",
        grouping_var
      )

      update_modal_progress(value = 0.95, text = "Finalising results...")
      mean_interaction_scores <- AverageExpression(
        seurat_with_scores,
        assays = "MATRISOMEPAIR",
        group.by = grouping_var
      )[[1]]

      interaction_results <- list(
        stats = interaction_stats,
        means = mean_interaction_scores
      )

      rm(seurat_with_scores, interaction_stats, mean_interaction_scores)
      gc(verbose = FALSE)

      # Store results for plotting (triggers heatmap and volcano plot rendering)
      rv$matrisome_pair_results <- interaction_results

    }, error = function(e) {
      # Include context to help debug failures
      error_msg <- paste0(
        "Matrisome pair analysis failed at sample '", nm(), "'\n",
        "Error: ", e$message, "\n",
        "Please reload the sample and try the analysis again."
      )
      showNotification(error_msg, type = "error", duration = 15)
      rv$matrisome_pair_results <- NULL
    }, finally = {
      remove_modal_progress()
      shinyjs::enable("run_matrisome_pair_analysis")
    })
  }, ignoreInit = TRUE)

  # ====================================================================== #
  #               MATRISOME-PAIR HEATMAP: TOP PAIRS BY NICHE             #
  # ====================================================================== #
  #
  # Shows the top 10 most enriched matrisome pairs per niche.
  # Helps identify shared and niche-specific pair co-expression patterns.
  #
  output$matrisome_pair_heatmap_plot <- renderPlot({
    req(rv$matrisome_pair_results)

    # Extract enrichment statistics and mean expression values
    stats_data <- rv$matrisome_pair_results$stats
    means_data <- rv$matrisome_pair_results$means

    # Filter by significance: >5% spatial enrichment, >0.5 log2FC, top 10 per niche
    top_axes <- stats_data %>%
      filter(!grepl("not.assigned|not_assigned", cluster, ignore.case = TRUE)) %>%  # Remove unassigned spots
      filter(perc_difference > 0.05 & avg_log2FC > 0.5) %>%  # Significance filters
      group_by(cluster) %>%
      slice_max(order_by = avg_log2FC, n = 10)  # Top 10 per niche

    # Stop and show message if no significant results are found
    validate(need(nrow(top_axes) > 0, "No significantly enriched matrisome pairs found with current filters."))

    # --- MATRIX PREPARATION ---
    # Build matrix: rows = matrisome pairs, columns = ECM niches
    # Values = mean co-expression scores across all spots in each niche
    matrix_to_plot <- means_data[unique(top_axes$feature), , drop = FALSE]
    # Remove "not assigned" column (handles both dot and underscore variants)
    cols_to_keep <- !grepl("not.assigned|not_assigned", colnames(matrix_to_plot), ignore.case = TRUE)
    matrix_to_plot_filtered <- matrix_to_plot[, cols_to_keep, drop = FALSE]

    # Clean up column names for display (replace underscores with spaces)
    clean_colnames <- gsub("_", " ", colnames(matrix_to_plot_filtered))

    # --- PHEATMAP PARAMETERS ---
    pheatmap(
      matrix_to_plot_filtered,
      labels_col = clean_colnames,
      scale = "row",           # Z-score normalization per pair
                               # Allows comparison of relative enrichment across niches
                               # without being biased by absolute expression levels
      main = "",
      fontsize_row = 8,        # Small font to fit many pair names
      fontsize_col = 10,       # Slightly larger for niche names (fewer columns)
      border_color = "grey60",
      cluster_rows = TRUE,     # Groups pairs with similar enrichment patterns
      cluster_cols = TRUE,     # Groups niches with similar pair profiles
      angle_col = 45,          # Diagonal labels prevent overlap
      cellwidth = 30,          # Wide enough to read niche names
      cellheight = 10          # Compact rows to fit all pairs
    )
  })

  # ====================================================================== #
  #        MATRISOME-PAIR VOLCANO PLOT: NICHE-SPECIFIC ENRICHMENT        #
  # ====================================================================== #
  #
  # Detailed view of all matrisome pairs for the selected niche.
  # X-axis: Spatial enrichment (% difference in co-expressing spots)
  # Y-axis: Expression enrichment (log2 fold change in co-expression score)
  #
  output$matrisome_pair_volcano_plot <- renderPlot({
    req(rv$matrisome_pair_results, input$matrisome_pair_cluster_select)

    # Standardize niche name to match analysis results (spaces/hyphens → underscores)
    selected_cluster_clean <- gsub(" |\\-", "_", input$matrisome_pair_cluster_select)
    display_niche <- trimws(gsub("\\s+", " ", gsub("_+", " ", input$matrisome_pair_cluster_select)))

    # Filter statistics to selected niche
    stats_subset <- rv$matrisome_pair_results$stats %>% filter(cluster == selected_cluster_clean)
    validate(need(
      nrow(stats_subset) > 0,
      paste("No data available for niche:", input$matrisome_pair_cluster_select)
    ))

    # --- IDENTIFY TOP HITS FOR LABELING ---
    # Apply same significance filters as heatmap (perc_difference > 0.05, avg_log2FC > 0.5)
    # Label top 20 to avoid overcrowding (ggrepel handles overlap automatically)
    top_hits_to_label <- stats_subset %>%
      filter(avg_log2FC > 0.5 & perc_difference > 0.05) %>%
      slice_max(order_by = avg_log2FC, n = 20)

    # Symmetric y-axis centered at zero
    y_extreme <- max(abs(quantile(stats_subset$avg_log2FC, c(0.01, 0.99), na.rm = TRUE)))
    y_limit <- min(y_extreme * 1.1, 10)
    y_limits <- c(-y_limit, y_limit)

    # --- VOLCANO PLOT CONSTRUCTION ---
    ggplot(stats_subset, aes(x = perc_difference, y = avg_log2FC)) +
      geom_point(color = "grey", alpha = 0.5, size = 0.6) +
      geom_point(data = top_hits_to_label, color = "red", alpha = 0.8, size = 1.2) +
      geom_label_repel(
        data = top_hits_to_label, aes(label = feature),
        size = 2.5, max.overlaps = Inf, box.padding = 0.6,
        force = 3, force_pull = 0.5, min.segment.length = 0,
        fill = "white", label.size = 0.1
      ) +
      # Reference lines marking significance thresholds
      geom_hline(yintercept = 0.5, linetype = "dashed") +  # avg_log2FC threshold
      geom_vline(xintercept = 0.05, linetype = "dashed") +  # perc_difference threshold
      theme_bw(base_size = 14) +
      labs(
        title = paste("Matrisome Pair Enrichment in the", display_niche, "Niche"),
        x = "Percentage Point Difference (In-Niche vs. Out-of-Niche)",
        y = "Enrichment Score (Avg. Log2 Fold Change)"
      ) +
      coord_cartesian(ylim = y_limits)
  })

  # --- [ EXPORT BUTTON STATE CONTROL ] ---
  # Enable the export dropdown if either analysis has been completed
  observe({
    if (isTRUE(rv$feature_analysis_done) || isTRUE(rv$matrisome_analysis_done)) {
      shinyjs::runjs("$('#export_dropdown_btn').prop('disabled', false).removeClass('initially-disabled');")
    }
  })

  # --- [ EXPORT ALL (ZIP) ] ---
  matrisome_coexpression_stats_filename <- "matrisome_coexpression_stats.csv"
  matrisome_coexpression_means_filename <- "matrisome_coexpression_means.csv"
  matrisome_coexpression_heatmap_name <- "30_matrisome_coexpression_heatmap"
  matrisome_coexpression_volcano_prefix <- "31_matrisome_coexpression_volcano_"

  output$export_all <- downloadHandler(
    filename = function() {
      paste0("MatriSpace_Export_", gsub(" ", "_", nm()), "_", Sys.Date(), ".zip")
    },
    content = function(file) {
      # --- MASTER ERROR HANDLER ---
      tryCatch({

        withProgress(message = 'Preparing Export Package', value = 0, {

          # 1. Setup temporary directory structure
          incProgress(0.1, detail = "Setting up directories...")
          temp_dir <- file.path(tempdir(), "matrispace_export")
          if (dir.exists(temp_dir)) unlink(temp_dir, recursive = TRUE)
          dir.create(temp_dir)

          dir_plots_png <- file.path(temp_dir, "Plots", "PNG_for_Viewing")
          dir_plots_pdf <- file.path(temp_dir, "Plots", "PDF_for_Publication")
          dir_data <- file.path(temp_dir, "Data")
          dir.create(dir_plots_png, recursive = TRUE)
          dir.create(dir_plots_pdf, recursive = TRUE)
          dir.create(dir_data)

          # --- ROBUSTLY GET THE FINAL SEURAT OBJECT ---
          # Use matrisome results object if available (has all scores), otherwise use d0()
          seurat_obj_for_export <- NULL
          if (isTRUE(rv$matrisome_analysis_done) && !is.null(rv$matrisome_results$seurat_object)) {
            seurat_obj_for_export <- rv$matrisome_results$seurat_object
          } else {
            seurat_obj_for_export <- d0()
          }
          if (is.null(seurat_obj_for_export)) {
              stop("Could not find a valid Seurat object to export.")
          }

          # --- 2. EXPORT DATA ---
          incProgress(0.1, detail = "Exporting data tables...")

          # FEATURE DATA EXPORT - Only if feature analysis was done
          if (!is.null(rv$analysis_results)) {
            tryCatch({
              df_main <- rv$analysis_results$data[, c("imagerow", "imagecol", ann(), "feature1", "feature2")]
              colnames(df_main) <- c("ImageRow", "ImageCol", "Annotation", "Primary_Feature", "Secondary_Feature")
              write.csv(df_main, file.path(dir_data, "main_feature_data.csv"))
            }, error = function(e){ warning("Failed to export main_feature_data.csv") })

            if (length(rv$analysis_results$custom_gene_sets) > 0) {
              tryCatch({
                assay_name <- if (is.null(rv$analysis_results$custom_score_assay)) "" else rv$analysis_results$custom_score_assay
                layer_name <- if (is.null(rv$analysis_results$custom_score_layer)) "" else rv$analysis_results$custom_score_layer
                custom_rows <- lapply(names(rv$analysis_results$custom_gene_sets), function(id) {
                  display_mode <- if (
                    (id == "customlist1" && identical(rv$primary_contrast_state, "contrast")) ||
                    (id == "customlist2" && identical(rv$secondary_contrast_state, "contrast"))
                  ) "relative_contrast" else "standard_ucell"
                  data.frame(
                    internal_id = id,
                    feature_role = if (id == "customlist1") "primary" else "secondary",
                    gene = rv$analysis_results$custom_gene_sets[[id]],
                    assay = assay_name,
                    layer = layer_name,
                    scoring_method = "UCell",
                    display_mode = display_mode,
                    stringsAsFactors = FALSE
                  )
                })
                write.csv(do.call(rbind, custom_rows), file.path(dir_data, "custom_gene_sets.csv"), row.names = FALSE)
              }, error = function(e){ warning("Failed to export custom_gene_sets.csv") })
            }

            tryCatch({
              stats_list <- list()
              f_stats <- feature_autocorrelation_stats()
              if(!is.null(f_stats$primary)) stats_list[['Primary_Morans_I']] <- f_stats$primary$morans_I
              if(!is.null(f_stats$primary)) stats_list[['Primary_Morans_p_value']] <- f_stats$primary$p_value
              if(!is.null(f_stats$secondary)) stats_list[['Secondary_Morans_I']] <- f_stats$secondary$morans_I
              if(!is.null(f_stats$secondary)) stats_list[['Secondary_Morans_p_value']] <- f_stats$secondary$p_value

              c_stats <- crosscorrelation_stats()
              if(!is.null(c_stats$pearson)) stats_list[['Pearson_Correlation_r']] <- c_stats$pearson$cor
              if(!is.null(c_stats$pearson)) stats_list[['Pearson_Correlation_p_value']] <- c_stats$pearson$pval
              if(!is.null(c_stats$spatial)) stats_list[['Spatial_Cross_Correlation']] <- c_stats$spatial$cor

              df_stats <- data.frame(Metric = names(stats_list), Value = unlist(stats_list))
              write.csv(df_stats, file.path(dir_data, "spatial_statistics_summary.csv"), row.names = FALSE)
            }, error = function(e){ warning("Failed to export spatial_statistics_summary.csv") })
          }

          # MATRISOME DATA EXPORT - Only if matrisome analysis was done
          if (isTRUE(rv$matrisome_analysis_done)) {
            tryCatch({
              score_cols <- grep("_score|_UCell", colnames(seurat_obj_for_export@meta.data), value = TRUE)
              df_scores <- seurat_obj_for_export@meta.data[, score_cols, drop = FALSE]
              write.csv(df_scores, file.path(dir_data, "matrisome_scores.csv"))
            }, error = function(e){ warning("Failed to export matrisome_scores.csv") })
          }

          # Matrisome-pair co-expression data export
          if (!is.null(rv$matrisome_pair_results)) {
            tryCatch({
              write.csv(rv$matrisome_pair_results$stats, file.path(dir_data, matrisome_coexpression_stats_filename), row.names = FALSE)
              write.csv(rv$matrisome_pair_results$means, file.path(dir_data, matrisome_coexpression_means_filename))
            }, error = function(e){ warning("Failed to export matrisome-pair co-expression data") })
          }

          # FULL METADATA EXPORT - Always export
          tryCatch({
            write.csv(seurat_obj_for_export@meta.data, file.path(dir_data, "spot_metadata.csv"))
          }, error = function(e){ warning("Failed to export spot_metadata.csv") })

          # --- 3. EXPORT PLOTS ---
          save_plot <- function(plot_obj, name, width = 7, height = 6) {
            ggsave(file.path(dir_plots_png, paste0(name, ".png")), plot = plot_obj, width = width, height = height, dpi = 300)
            ggsave(file.path(dir_plots_pdf, paste0(name, ".pdf")), plot = plot_obj, width = width, height = height)
          }

          # Base Tissue Image with Scale Bar
          incProgress(0.1, detail = "Exporting Base Tissue Image...")
          tryCatch({
            p_tissue_base <- SpatialDimPlot(seurat_obj_for_export, label = FALSE, alpha = 0.0, crop = FALSE) + NoLegend()
            p_tissue_with_bar <- add_scale_bar(p_tissue_base, seurat_obj_for_export)
            save_plot(p_tissue_with_bar, "00_tissue_image_with_scale_bar")
          }, error = function(e){ warning("Failed to export Base Tissue Image") })

          # Annotation Plot
          incProgress(0.1, detail = "Exporting Annotation plot...")
          tryCatch({
            # Build annotation data from Seurat object
            coords <- GetTissueCoordinates(seurat_obj_for_export)
            annot_col <- ann()
            coords$annotation <- seurat_obj_for_export@meta.data[[annot_col]]

            p_annot <- ggplot(coords, aes(x = imagecol, y = imagerow, colour = annotation)) +
              geom_point(size = 1.5) + scale_color_manual(values = color_map()) +
              coord_fixed() + theme_void() + scale_y_reverse() + labs(colour = "Annotation")

            p_annot_with_bar <- add_scale_bar(p_annot, seurat_obj_for_export)
            save_plot(p_annot_with_bar, "01_annotation_plot_with_scale_bar")

          }, error = function(e){ warning("Failed to export Annotation Plot") })

          # Feature Plots - Only if feature analysis was done
          if (!is.null(rv$analysis_results)) {
            incProgress(0.1, detail = "Exporting Feature plots...")
            tryCatch({
              plot_data <- rv$analysis_results$data
              pal <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
              primary_relative <- identical(rv$primary_contrast_state, "contrast") &&
                !is.null(rv$analysis_results$custom_contrast_scores$customlist1)
              plot_data$.feature1_display <- if (primary_relative) {
                rv$analysis_results$custom_contrast_scores$customlist1[rownames(plot_data)]
              } else plot_data$feature1
              p_feat1 <- ggplot(plot_data, aes(x = imagecol, y = imagerow, colour = .feature1_display)) +
                geom_point(size = 1.5) + scale_colour_gradientn(colours = pal, limits = if (primary_relative) c(0, 1) else NULL) +
                coord_fixed() + theme_void() + scale_y_reverse() +
                labs(colour = if (primary_relative) "UCell relative contrast" else rv$analysis_results$gene1)
              save_plot(p_feat1, "02_primary_feature")

              if (rv$analysis_results$sel2 != "none") {
                secondary_relative <- identical(rv$secondary_contrast_state, "contrast") &&
                  !is.null(rv$analysis_results$custom_contrast_scores$customlist2)
                plot_data$.feature2_display <- if (secondary_relative) {
                  rv$analysis_results$custom_contrast_scores$customlist2[rownames(plot_data)]
                } else plot_data$feature2
                p_feat2 <- ggplot(plot_data, aes(x = imagecol, y = imagerow, colour = .feature2_display)) +
                  geom_point(size = 1.5) + scale_colour_gradientn(colours = pal, limits = if (secondary_relative) c(0, 1) else NULL) +
                  coord_fixed() + theme_void() + scale_y_reverse() +
                  labs(colour = if (secondary_relative) "UCell relative contrast" else rv$analysis_results$gene2)
                save_plot(p_feat2, "03_secondary_feature")
              }
            }, error = function(e){ warning("Failed to export Feature Plots") })

            # Co-expression Plot
            incProgress(0.1, detail = "Exporting Co-expression plot...")
            tryCatch({
              if (rv$analysis_results$sel2 != "none") {
                temp_seurat_obj <- d0()
                cell_order <- intersect(Cells(temp_seurat_obj), rownames(rv$analysis_results$data))
                if (length(cell_order) > 0) {
                  temp_seurat_obj <- temp_seurat_obj[, cell_order]
                  temp_seurat_obj$feature_blend_1 <- rv$analysis_results$data[cell_order, "feature1"]
                  temp_seurat_obj$feature_blend_2 <- rv$analysis_results$data[cell_order, "feature2"]
                  colors_to_use <- coex_blend_colors(input$coex_palette)
                  full_plot_object <- SpatialFeaturePlotBlend(
                    object = temp_seurat_obj,
                    features = c("feature_blend_1", "feature_blend_2"),
                    feature_1_alt_name = rv$analysis_results$gene1,
                    feature_2_alt_name = rv$analysis_results$gene2,
                    combine = TRUE,
                    bottom_left = colors_to_use$bottom_left,
                    bottom_right = colors_to_use$bottom_right,
                    top_left = colors_to_use$top_left,
                    top_right = colors_to_use$top_right,
                    sfp_extra_arguments = list(pt.size.factor = effective_pt_size(input$coex_pt_size), alpha = input$coex_alpha)
                  )
                  inner_plots <- full_plot_object$patches$plots[[1]]$patches$plots
                  blended_plot <- inner_plots[[3]] + theme(plot.title = element_blank())
                  legend_plot <- inner_plots[[4]]
                  final_coex_plot <- wrap_plots(blended_plot, legend_plot, nrow = 1, widths = c(0.75, 0.25))
                  save_plot(final_coex_plot, "04_co_expression", width = 10, height = 7)
                }
              }
            }, error = function(e){ warning("Failed to export Co-expression Plot") })

            tryCatch({
              if (rv$analysis_results$sel2 != "none") {
                p_lisa <- plot_lisa_clustering(
                  k = rv$analysis_results$data,
                  selectedCellTypes = unique(rv$analysis_results$data[[ann()]]),
                  ann_col = ann(),
                  f1_name = rv$analysis_results$gene1,
                  f2_name = rv$analysis_results$gene2
                )
                save_plot(p_lisa, "05_lisa_clustering", width = 12, height = 6)
              }
            }, error = function(e){ warning("Failed to export LISA Plot") })
          }

          # Matrisome & ECM Plots - Only if matrisome analysis was done
          if (isTRUE(rv$matrisome_analysis_done)) {
            incProgress(0.3, detail = "Exporting Matrisome & ECM plots...")
            tryCatch({
              main_divisions <- c("glycoprotein", "collagen", "proteoglycan", "regulator", "secreted_factor", "affliated_protein")
              for (group in main_divisions) {
                if (!is.null(rv$matrisome_results[[group]])) {
                  p_dist <- SpatialFeaturePlot(seurat_obj_for_export, features = paste0(group, "_robust_score"), pt.size.factor = effective_pt_size(input$main_pt_size)) + theme(legend.position = "right")
                  p_hot <- SpatialFeaturePlot(seurat_obj_for_export, features = paste0(group, "_log_scaled_score"), pt.size.factor = effective_pt_size(input$main_pt_size)) + theme(legend.position = "right")
                  save_plot(p_dist, paste0("10_matrisome_", group, "_distribution"))
                  save_plot(p_hot, paste0("11_matrisome_", group, "_hotspot"))
                }
              }

              if ("ecm_domain_annotation" %in% colnames(seurat_obj_for_export@meta.data)) {
                p_ecm_domain <- SpatialDimPlot(seurat_obj_for_export, group.by = "ecm_domain_annotation", pt.size.factor = effective_pt_size(input$ecm_domain_pt_size))
                save_plot(p_ecm_domain, "20_ecm_domain_annotation")
              }

              ecm_sigs <- c("Interstitial", "Basement")
              for(sig in ecm_sigs) {
                p_sig <- SpatialFeaturePlot(seurat_obj_for_export, features = paste0(sig, "_UCell"), pt.size.factor = effective_pt_size(input$ecm_pt_size))
                save_plot(p_sig, paste0("21_ecm_niche_score_", tolower(sig)))
              }

              ecm_subcategories <- c("perivascular", "hemostasis", "elastic_fibers", "gf_binding")
              for (subcat in ecm_subcategories) {
                if (paste0(subcat, "_robust_score") %in% colnames(seurat_obj_for_export@meta.data)) {
                  p_dist <- SpatialFeaturePlot(seurat_obj_for_export,
                                               features = paste0(subcat, "_robust_score"),
                                               pt.size.factor = effective_pt_size(input$subcategories_pt_size)) +
                    theme(legend.position = "right")
                  p_hot <- SpatialFeaturePlot(seurat_obj_for_export,
                                              features = paste0(subcat, "_log_scaled_score"),
                                              pt.size.factor = effective_pt_size(input$subcategories_pt_size)) +
                    theme(legend.position = "right")
                  save_plot(p_dist, paste0("12_matrisome_subcategory_", subcat, "_distribution"))
                  save_plot(p_hot, paste0("13_matrisome_subcategory_", subcat, "_hotspot"))
                }
              }

              ecm_gene_families <- c("laminin", "matricellular", "syndecan", "glypican")
              for (family in ecm_gene_families) {
                if (paste0(family, "_robust_score") %in% colnames(seurat_obj_for_export@meta.data)) {
                  p_dist <- SpatialFeaturePlot(seurat_obj_for_export,
                                               features = paste0(family, "_robust_score"),
                                               pt.size.factor = effective_pt_size(input$families_pt_size)) +
                    theme(legend.position = "right")
                  p_hot <- SpatialFeaturePlot(seurat_obj_for_export,
                                              features = paste0(family, "_log_scaled_score"),
                                              pt.size.factor = effective_pt_size(input$families_pt_size)) +
                    theme(legend.position = "right")
                  save_plot(p_dist, paste0("14_matrisome_gene_family_", family, "_distribution"))
                  save_plot(p_hot, paste0("15_matrisome_gene_family_", family, "_hotspot"))
                }
              }

              tryCatch({
                sig_data <- signature_distribution_data()
                if (!is.null(sig_data) && nrow(sig_data) > 0) {
                  p_niche_dist <- ggplot(sig_data, aes(x = annotation_group, y = mean_score, fill = signature_name)) +
                    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
                    scale_fill_manual(values = c("Basement membrane" = "#6a3d9aff",
                                                 "Interstitial ECM" = "#2b9e2bff")) +
                    labs(x = NULL, y = "Average Niche Score", fill = NULL) +
                    theme_minimal(base_size = 14) +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                          legend.position = "top")
                  save_plot(p_niche_dist, "22_niche_distribution", width = 10, height = 6)
                }
              }, error = function(e) warning("Failed to export niche distribution plot"))
            }, error = function(e){ warning("Failed to export Matrisome/ECM Plots") })
          }

          # Matrisome-pair co-expression plots
          if (!is.null(rv$matrisome_pair_results)) {
            incProgress(0.1, detail = "Exporting matrisome-pair co-expression plots...")
            tryCatch({
              stats_data <- rv$matrisome_pair_results$stats
              means_data <- rv$matrisome_pair_results$means

              # Matrisome-pair heatmap
              top_axes <- stats_data %>%
                filter(!grepl("not.assigned|not_assigned", cluster, ignore.case = TRUE)) %>%
                filter(perc_difference > 0.05 & avg_log2FC > 0.5) %>%
                group_by(cluster) %>%
                slice_max(order_by = avg_log2FC, n = 10)

              if (nrow(top_axes) > 0) {
                matrix_to_plot <- means_data[unique(top_axes$feature), , drop = FALSE]
                cols_to_keep <- !grepl("not.assigned|not_assigned", colnames(matrix_to_plot), ignore.case = TRUE)
                matrix_to_plot_filtered <- matrix_to_plot[, cols_to_keep, drop = FALSE]
                clean_colnames <- gsub("_", " ", colnames(matrix_to_plot_filtered))

                # Save heatmap
                png(file.path(dir_plots_png, paste0(matrisome_coexpression_heatmap_name, ".png")), width = 10, height = 8, units = "in", res = 300)
                pheatmap(matrix_to_plot_filtered, labels_col = clean_colnames, scale = "row",
                         fontsize_row = 8, fontsize_col = 10, border_color = "grey60",
                         cluster_rows = TRUE, cluster_cols = TRUE, angle_col = 45,
                         cellwidth = 30, cellheight = 10)
                dev.off()

                pdf(file.path(dir_plots_pdf, paste0(matrisome_coexpression_heatmap_name, ".pdf")), width = 10, height = 8)
                pheatmap(matrix_to_plot_filtered, labels_col = clean_colnames, scale = "row",
                         fontsize_row = 8, fontsize_col = 10, border_color = "grey60",
                         cluster_rows = TRUE, cluster_cols = TRUE, angle_col = 45,
                         cellwidth = 30, cellheight = 10)
                dev.off()
              }

              # Matrisome-pair volcano plots for each niche
              niches <- unique(stats_data$cluster)
              niches <- niches[!grepl("not.assigned|not_assigned", niches, ignore.case = TRUE)]

              for (niche in niches) {
                stats_subset <- stats_data %>% filter(cluster == niche)
                if (nrow(stats_subset) == 0) next

                top_hits <- stats_subset %>%
                  filter(avg_log2FC > 0.5 & perc_difference > 0.05) %>%
                  slice_max(order_by = avg_log2FC, n = 20)

                niche_clean <- trimws(gsub("\\s+", " ", gsub("_+", " ", niche)))

                y_extreme <- max(abs(quantile(stats_subset$avg_log2FC, c(0.01, 0.99), na.rm = TRUE)))
                y_limit <- min(y_extreme * 1.1, 10)
                y_limits <- c(-y_limit, y_limit)

                p_volcano <- ggplot(stats_subset, aes(x = perc_difference, y = avg_log2FC)) +
                  geom_point(color = "grey", alpha = 0.5, size = 0.6) +
                  geom_point(data = top_hits, color = "red", alpha = 0.8, size = 1.2) +
                  geom_label_repel(data = top_hits, aes(label = feature),
                                  size = 2.5, max.overlaps = Inf, box.padding = 0.6,
                                  force = 3, force_pull = 0.5, min.segment.length = 0,
                                  fill = "white", label.size = 0.1) +
                  geom_hline(yintercept = 0.5, linetype = "dashed") +
                  geom_vline(xintercept = 0.05, linetype = "dashed") +
                  theme_bw(base_size = 14) +
                  labs(title = paste("Matrisome Pair Enrichment in the", niche_clean, "Niche"),
                       x = "% Difference (In vs Out)", y = "Avg Log2 FC") +
                  coord_cartesian(ylim = y_limits)

                save_plot(p_volcano, paste0(matrisome_coexpression_volcano_prefix, gsub(" ", "_", niche)), width = 8, height = 6)
              }
            }, error = function(e){ warning("Failed to export matrisome-pair co-expression plots") })
          }

          # --- EXPORT SEURAT OBJECT ---
          incProgress(0.1, detail = "Exporting Seurat object...")
          tryCatch({
            seurat_filename <- paste0(gsub(" ", "_", nm()), "_MatriSpace.rds")
            saveRDS(seurat_obj_for_export, file.path(dir_data, seurat_filename))
          }, error = function(e){ warning("Failed to export Seurat object") })

          # --- 4. CREATE README ---
          incProgress(0.1, detail = "Finalizing package...")
          seurat_filename <- paste0(gsub(" ", "_", nm()), "_MatriSpace.rds")
          seurat_export_path <- file.path(dir_data, seurat_filename)
          data_readme_lines <- c(
            "- Data/: Contains the raw data tables generated by the analysis."
          )
          if (file.exists(file.path(dir_data, "spot_metadata.csv"))) {
            data_readme_lines <- c(
              data_readme_lines,
              "  - spot_metadata.csv: Complete metadata for each spot (coordinates, annotations, scores)."
            )
          }
          if (!is.null(rv$analysis_results) && file.exists(file.path(dir_data, "main_feature_data.csv"))) {
            data_readme_lines <- c(
              data_readme_lines,
              "  - main_feature_data.csv: Expression values for the selected features for each spot."
            )
          }
          if (file.exists(file.path(dir_data, "custom_gene_sets.csv"))) {
            data_readme_lines <- c(
              data_readme_lines,
              "  - custom_gene_sets.csv: Matrisome genes selected for custom primary or secondary feature scores."
            )
          }
          if (isTRUE(rv$matrisome_analysis_done) && file.exists(file.path(dir_data, "matrisome_scores.csv"))) {
            data_readme_lines <- c(
              data_readme_lines,
              "  - matrisome_scores.csv: All calculated matrisome gene list scores for each spot."
            )
          }
          if (!is.null(rv$analysis_results) && file.exists(file.path(dir_data, "spatial_statistics_summary.csv"))) {
            data_readme_lines <- c(
              data_readme_lines,
              "  - spatial_statistics_summary.csv: A summary of all calculated correlation statistics."
            )
          }
          if (!is.null(rv$matrisome_pair_results) && file.exists(file.path(dir_data, matrisome_coexpression_stats_filename))) {
            data_readme_lines <- c(
              data_readme_lines,
              paste0("  - ", matrisome_coexpression_stats_filename, ": Matrisome-pair co-expression enrichment statistics per niche.")
            )
          }
          if (!is.null(rv$matrisome_pair_results) && file.exists(file.path(dir_data, matrisome_coexpression_means_filename))) {
            data_readme_lines <- c(
              data_readme_lines,
              paste0("  - ", matrisome_coexpression_means_filename, ": Mean co-expression scores for each matrisome pair per niche.")
            )
          }
          if (file.exists(seurat_export_path)) {
            data_readme_lines <- c(
              data_readme_lines,
              paste0("  - ", seurat_filename, ": Complete Seurat object with all analysis results.")
            )
          }
          readme_text <- c(
            "MatriSpace Analysis Export",
            "==========================",
            paste("Sample Name:", nm()),
            paste("Export Date:", Sys.Date()),
            "",
            "This ZIP file contains the following directories:",
            "- Plots/: Contains all generated plots in two formats.",
            "  - PNG_for_Viewing/: High-resolution PNGs for presentations and quick viewing.",
            "  - PDF_for_Publication/: Vector-based PDFs for publication and editing.",
            data_readme_lines
          )
          writeLines(readme_text, file.path(temp_dir, "README.txt"))

          # --- 5. ZIP AND SERVE (ROBUST METHOD) ---
          # This pattern avoids Windows file locking issues by changing the
          # working directory before creating the archive.

          # Save the current working directory
          old_wd <- getwd()
          # Change to the temporary directory
          setwd(temp_dir)
          # Get a list of all files with relative paths (e.g., "Data/file.csv")
          files_to_zip <- list.files(".", recursive = TRUE)
          # Create the zip archive. 'file' is the full path to the destination
          zip::zip(zipfile = file, files = files_to_zip)
          # Restore the original working directory
          setwd(old_wd)
        })
      }, error = function(e) {
        # This code runs if any part of the tryCatch fails
        showNotification(
          paste("Export failed:", e$message),
          type = "error",
          duration = 10
        )
      })
    },
    contentType = "application/zip"
  )

  # --- [ EXPORT METADATA (CSV) ] ---
  output$export_metadata <- downloadHandler(
    filename = function() {
      paste0("MatriSpace_Metadata_", gsub(" ", "_", nm()), "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      tryCatch({
        # Use matrisome results object if available (has all scores)
        seurat_obj <- NULL
        if (isTRUE(rv$matrisome_analysis_done) && !is.null(rv$matrisome_results$seurat_object)) {
          seurat_obj <- rv$matrisome_results$seurat_object
        } else {
          seurat_obj <- d0()
        }
        req(seurat_obj)
        write.csv(seurat_obj@meta.data, file)
      }, error = function(e) {
        showNotification(paste("Export failed:", e$message), type = "error", duration = 10)
      })
    },
    contentType = "text/csv"
  )

  # --- [ EXPORT MATRISOME SCORES (CSV) ] ---
  output$export_scores <- downloadHandler(
    filename = function() {
      paste0("MatriSpace_Scores_", gsub(" ", "_", nm()), "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      tryCatch({
        req(rv$matrisome_analysis_done, rv$matrisome_results$seurat_object)
        seurat_obj <- rv$matrisome_results$seurat_object
        score_cols <- grep("_score|_UCell", colnames(seurat_obj@meta.data), value = TRUE)
        if (length(score_cols) == 0) {
          showNotification("No matrisome scores found. Run Matrisome Profile first.", type = "warning")
          return(NULL)
        }
        df_scores <- seurat_obj@meta.data[, score_cols, drop = FALSE]
        write.csv(df_scores, file)
      }, error = function(e) {
        showNotification(paste("Export failed:", e$message), type = "error", duration = 10)
      })
    },
    contentType = "text/csv"
  )

  # --- [ EXPORT PLOTS ONLY (ZIP) ] ---
  output$export_plots <- downloadHandler(
    filename = function() {
      paste0("MatriSpace_Plots_", gsub(" ", "_", nm()), "_", Sys.Date(), ".zip")
    },
    content = function(file) {
      tryCatch({
        withProgress(message = 'Exporting Plots', value = 0, {
          # Use matrisome results object if available (has all scores), otherwise use d0()
          seurat_obj <- NULL
          if (isTRUE(rv$matrisome_analysis_done) && !is.null(rv$matrisome_results$seurat_object)) {
            seurat_obj <- rv$matrisome_results$seurat_object
          } else {
            seurat_obj <- d0()
          }
          req(seurat_obj)

          # Setup temp directory
          temp_dir <- file.path(tempdir(), "matrispace_plots_export")
          if (dir.exists(temp_dir)) unlink(temp_dir, recursive = TRUE)
          dir.create(temp_dir)

          dir_png <- file.path(temp_dir, "PNG")
          dir_pdf <- file.path(temp_dir, "PDF")
          dir.create(dir_png)
          dir.create(dir_pdf)

          save_plot <- function(plot_obj, name, width = 7, height = 6) {
            ggsave(file.path(dir_png, paste0(name, ".png")), plot = plot_obj, width = width, height = height, dpi = 300)
            ggsave(file.path(dir_pdf, paste0(name, ".pdf")), plot = plot_obj, width = width, height = height)
          }

          # Base tissue image
          incProgress(0.2, detail = "Tissue image...")
          tryCatch({
            p_tissue <- SpatialDimPlot(seurat_obj, label = FALSE, alpha = 0.0, crop = FALSE) + NoLegend()
            p_tissue <- add_scale_bar(p_tissue, seurat_obj)
            save_plot(p_tissue, "00_tissue_image")
          }, error = function(e) warning("Failed: tissue image"))

          # Annotation plot
          incProgress(0.2, detail = "Annotation plot...")
          tryCatch({
            # Build annotation data from Seurat object
            coords <- GetTissueCoordinates(seurat_obj)
            annot_col <- ann()
            coords$annotation <- seurat_obj@meta.data[[annot_col]]

            p_annot <- ggplot(coords, aes(x = imagecol, y = imagerow, colour = annotation)) +
              geom_point(size = 1.5) + scale_color_manual(values = color_map()) +
              coord_fixed() + theme_void() + scale_y_reverse() + labs(colour = "Annotation")
            p_annot <- add_scale_bar(p_annot, seurat_obj)
            save_plot(p_annot, "01_annotation_plot")
          }, error = function(e) warning("Failed: annotation plot"))

          # Feature plots (if feature analysis done)
          if (!is.null(rv$analysis_results)) {
            incProgress(0.2, detail = "Feature plots...")
            tryCatch({
              plot_data <- rv$analysis_results$data
              pal <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
              primary_relative <- identical(rv$primary_contrast_state, "contrast") &&
                !is.null(rv$analysis_results$custom_contrast_scores$customlist1)
              plot_data$.feature1_display <- if (primary_relative) {
                rv$analysis_results$custom_contrast_scores$customlist1[rownames(plot_data)]
              } else plot_data$feature1
              p_feat1 <- ggplot(plot_data, aes(x = imagecol, y = imagerow, colour = .feature1_display)) +
                geom_point(size = 1.5) + scale_colour_gradientn(colours = pal, limits = if (primary_relative) c(0, 1) else NULL) +
                coord_fixed() + theme_void() + scale_y_reverse() +
                labs(colour = if (primary_relative) "UCell relative contrast" else rv$analysis_results$gene1)
              save_plot(p_feat1, "02_primary_feature")

              if (rv$analysis_results$sel2 != "none") {
                secondary_relative <- identical(rv$secondary_contrast_state, "contrast") &&
                  !is.null(rv$analysis_results$custom_contrast_scores$customlist2)
                plot_data$.feature2_display <- if (secondary_relative) {
                  rv$analysis_results$custom_contrast_scores$customlist2[rownames(plot_data)]
                } else plot_data$feature2
                p_feat2 <- ggplot(plot_data, aes(x = imagecol, y = imagerow, colour = .feature2_display)) +
                  geom_point(size = 1.5) + scale_colour_gradientn(colours = pal, limits = if (secondary_relative) c(0, 1) else NULL) +
                  coord_fixed() + theme_void() + scale_y_reverse() +
                  labs(colour = if (secondary_relative) "UCell relative contrast" else rv$analysis_results$gene2)
                save_plot(p_feat2, "03_secondary_feature")

                temp_seurat_obj <- seurat_obj
                cell_order <- intersect(Cells(temp_seurat_obj), rownames(rv$analysis_results$data))
                if (length(cell_order) > 0) {
                  temp_seurat_obj <- temp_seurat_obj[, cell_order]
                  temp_seurat_obj$feature_blend_1 <- rv$analysis_results$data[cell_order, "feature1"]
                  temp_seurat_obj$feature_blend_2 <- rv$analysis_results$data[cell_order, "feature2"]
                  colors_to_use <- coex_blend_colors(input$coex_palette)
                  full_plot_object <- SpatialFeaturePlotBlend(
                    object = temp_seurat_obj,
                    features = c("feature_blend_1", "feature_blend_2"),
                    feature_1_alt_name = rv$analysis_results$gene1,
                    feature_2_alt_name = rv$analysis_results$gene2,
                    combine = TRUE,
                    bottom_left = colors_to_use$bottom_left,
                    bottom_right = colors_to_use$bottom_right,
                    top_left = colors_to_use$top_left,
                    top_right = colors_to_use$top_right,
                    sfp_extra_arguments = list(pt.size.factor = effective_pt_size(input$coex_pt_size), alpha = input$coex_alpha)
                  )
                  inner_plots <- full_plot_object$patches$plots[[1]]$patches$plots
                  blended_plot <- inner_plots[[3]] + theme(plot.title = element_blank())
                  legend_plot <- inner_plots[[4]]
                  final_coex_plot <- wrap_plots(blended_plot, legend_plot, nrow = 1, widths = c(0.75, 0.25))
                  save_plot(final_coex_plot, "04_co_expression", width = 10, height = 7)
                }

                tryCatch({
                  p_lisa <- plot_lisa_clustering(
                    k = rv$analysis_results$data,
                    selectedCellTypes = unique(rv$analysis_results$data[[ann()]]),
                    ann_col = ann(),
                    f1_name = rv$analysis_results$gene1,
                    f2_name = rv$analysis_results$gene2
                  )
                  save_plot(p_lisa, "05_lisa_clustering", width = 12, height = 6)
                }, error = function(e) warning("Failed to export LISA plot"))
              }
            }, error = function(e) warning("Failed: feature plots"))
          }

          # Matrisome plots (if matrisome analysis done)
          if (isTRUE(rv$matrisome_analysis_done)) {
            incProgress(0.3, detail = "Matrisome plots...")
            tryCatch({
              main_divisions <- c("glycoprotein", "collagen", "proteoglycan", "regulator", "secreted_factor", "affliated_protein")
              for (group in main_divisions) {
                if (!is.null(rv$matrisome_results[[group]])) {
                  p_dist <- SpatialFeaturePlot(seurat_obj, features = paste0(group, "_robust_score"), pt.size.factor = effective_pt_size(input$main_pt_size)) + theme(legend.position = "right")
                  p_hot <- SpatialFeaturePlot(seurat_obj, features = paste0(group, "_log_scaled_score"), pt.size.factor = effective_pt_size(input$main_pt_size)) + theme(legend.position = "right")
                  save_plot(p_dist, paste0("10_matrisome_", group, "_distribution"))
                  save_plot(p_hot, paste0("11_matrisome_", group, "_hotspot"))
                }
              }

              if ("ecm_domain_annotation" %in% colnames(seurat_obj@meta.data)) {
                p_ecm <- SpatialDimPlot(seurat_obj, group.by = "ecm_domain_annotation", pt.size.factor = effective_pt_size(input$ecm_domain_pt_size))
                save_plot(p_ecm, "20_ecm_domain_annotation")
              }

              for (sig in c("Interstitial", "Basement")) {
                p_sig <- SpatialFeaturePlot(seurat_obj, features = paste0(sig, "_UCell"), pt.size.factor = effective_pt_size(input$ecm_pt_size))
                save_plot(p_sig, paste0("21_ecm_niche_score_", tolower(sig)))
              }

              ecm_subcategories <- c("perivascular", "hemostasis", "elastic_fibers", "gf_binding")
              for (subcat in ecm_subcategories) {
                if (paste0(subcat, "_robust_score") %in% colnames(seurat_obj@meta.data)) {
                  p_dist <- SpatialFeaturePlot(seurat_obj,
                                               features = paste0(subcat, "_robust_score"),
                                               pt.size.factor = effective_pt_size(input$subcategories_pt_size)) +
                    theme(legend.position = "right")
                  p_hot <- SpatialFeaturePlot(seurat_obj,
                                              features = paste0(subcat, "_log_scaled_score"),
                                              pt.size.factor = effective_pt_size(input$subcategories_pt_size)) +
                    theme(legend.position = "right")
                  save_plot(p_dist, paste0("12_matrisome_subcategory_", subcat, "_distribution"))
                  save_plot(p_hot, paste0("13_matrisome_subcategory_", subcat, "_hotspot"))
                }
              }

              ecm_gene_families <- c("laminin", "matricellular", "syndecan", "glypican")
              for (family in ecm_gene_families) {
                if (paste0(family, "_robust_score") %in% colnames(seurat_obj@meta.data)) {
                  p_dist <- SpatialFeaturePlot(seurat_obj,
                                               features = paste0(family, "_robust_score"),
                                               pt.size.factor = effective_pt_size(input$families_pt_size)) +
                    theme(legend.position = "right")
                  p_hot <- SpatialFeaturePlot(seurat_obj,
                                              features = paste0(family, "_log_scaled_score"),
                                              pt.size.factor = effective_pt_size(input$families_pt_size)) +
                    theme(legend.position = "right")
                  save_plot(p_dist, paste0("14_matrisome_gene_family_", family, "_distribution"))
                  save_plot(p_hot, paste0("15_matrisome_gene_family_", family, "_hotspot"))
                }
              }

              tryCatch({
                sig_data <- signature_distribution_data()
                if (!is.null(sig_data) && nrow(sig_data) > 0) {
                  p_niche_dist <- ggplot(sig_data, aes(x = annotation_group, y = mean_score, fill = signature_name)) +
                    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
                    scale_fill_manual(values = c("Basement membrane" = "#6a3d9aff",
                                                 "Interstitial ECM" = "#2b9e2bff")) +
                    labs(x = NULL, y = "Average Niche Score", fill = NULL) +
                    theme_minimal(base_size = 14) +
                    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
                          legend.position = "top")
                  save_plot(p_niche_dist, "22_niche_distribution", width = 10, height = 6)
                }
              }, error = function(e) warning("Failed to export niche distribution plot"))
            }, error = function(e) warning("Failed: matrisome plots"))
          }

          # Matrisome-pair co-expression plots
          if (!is.null(rv$matrisome_pair_results)) {
            incProgress(0.1, detail = "Matrisome-pair co-expression plots...")
            tryCatch({
              stats_data <- rv$matrisome_pair_results$stats
              means_data <- rv$matrisome_pair_results$means

              # Matrisome-pair heatmap
              top_axes <- stats_data %>%
                filter(!grepl("not.assigned|not_assigned", cluster, ignore.case = TRUE)) %>%
                filter(perc_difference > 0.05 & avg_log2FC > 0.5) %>%
                group_by(cluster) %>%
                slice_max(order_by = avg_log2FC, n = 10)

              if (nrow(top_axes) > 0) {
                matrix_to_plot <- means_data[unique(top_axes$feature), , drop = FALSE]
                cols_to_keep <- !grepl("not.assigned|not_assigned", colnames(matrix_to_plot), ignore.case = TRUE)
                matrix_to_plot_filtered <- matrix_to_plot[, cols_to_keep, drop = FALSE]
                clean_colnames <- gsub("_", " ", colnames(matrix_to_plot_filtered))

                png(file.path(dir_png, paste0(matrisome_coexpression_heatmap_name, ".png")), width = 10, height = 8, units = "in", res = 300)
                pheatmap(matrix_to_plot_filtered, labels_col = clean_colnames, scale = "row",
                         fontsize_row = 8, fontsize_col = 10, border_color = "grey60",
                         cluster_rows = TRUE, cluster_cols = TRUE, angle_col = 45,
                         cellwidth = 30, cellheight = 10)
                dev.off()

                pdf(file.path(dir_pdf, paste0(matrisome_coexpression_heatmap_name, ".pdf")), width = 10, height = 8)
                pheatmap(matrix_to_plot_filtered, labels_col = clean_colnames, scale = "row",
                         fontsize_row = 8, fontsize_col = 10, border_color = "grey60",
                         cluster_rows = TRUE, cluster_cols = TRUE, angle_col = 45,
                         cellwidth = 30, cellheight = 10)
                dev.off()
              }

              # Matrisome-pair volcano plots for each niche
              niches <- unique(stats_data$cluster)
              niches <- niches[!grepl("not.assigned|not_assigned", niches, ignore.case = TRUE)]

              for (niche in niches) {
                stats_subset <- stats_data %>% filter(cluster == niche)
                if (nrow(stats_subset) == 0) next

                top_hits <- stats_subset %>%
                  filter(avg_log2FC > 0.5 & perc_difference > 0.05) %>%
                  slice_max(order_by = avg_log2FC, n = 20)

                niche_clean <- trimws(gsub("\\s+", " ", gsub("_+", " ", niche)))

                y_extreme <- max(abs(quantile(stats_subset$avg_log2FC, c(0.01, 0.99), na.rm = TRUE)))
                y_limit <- min(y_extreme * 1.1, 10)
                y_limits <- c(-y_limit, y_limit)

                p_volcano <- ggplot(stats_subset, aes(x = perc_difference, y = avg_log2FC)) +
                  geom_point(color = "grey", alpha = 0.5, size = 0.6) +
                  geom_point(data = top_hits, color = "red", alpha = 0.8, size = 1.2) +
                  geom_label_repel(data = top_hits, aes(label = feature),
                                  size = 2.5, max.overlaps = Inf, box.padding = 0.6,
                                  force = 3, force_pull = 0.5, min.segment.length = 0,
                                  fill = "white", label.size = 0.1) +
                  geom_hline(yintercept = 0.5, linetype = "dashed") +
                  geom_vline(xintercept = 0.05, linetype = "dashed") +
                  theme_bw(base_size = 14) +
                  labs(title = paste("Matrisome Pair Enrichment in the", niche_clean, "Niche"),
                       x = "% Difference (In vs Out)", y = "Avg Log2 FC") +
                  coord_cartesian(ylim = y_limits)

                save_plot(p_volcano, paste0(matrisome_coexpression_volcano_prefix, gsub(" ", "_", niche)), width = 8, height = 6)
              }
            }, error = function(e) warning("Failed: matrisome-pair plots"))
          }

          # Create ZIP
          incProgress(0.1, detail = "Creating archive...")
          old_wd <- getwd()
          setwd(temp_dir)
          files_to_zip <- list.files(".", recursive = TRUE)
          zip::zip(zipfile = file, files = files_to_zip)
          setwd(old_wd)
        })
      }, error = function(e) {
        showNotification(paste("Export failed:", e$message), type = "error", duration = 10)
      })
    },
    contentType = "application/zip"
  )

}
