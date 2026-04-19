# app_ui.R
# User Interface definition for MatriSpace Shiny application

#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {

# User Interface ----

page_fillable(
  theme = bs_theme(version = 5),
  gap = 0,

  useShinyjs(),
  use_cicerone(),

  # Add native bslib busy indicators
  useBusyIndicators(),
  busyIndicatorOptions(),

  golem_add_external_resources(),


  # Full-width header spanning both sidebar and main content
  div(
    class = "app-header app-scale-target border-bottom",
    div(
      class = "d-flex justify-content-between align-items-center py-3 px-4",
      # Left - App name with enhanced styling
      div(
        id = "app-branding",
        class = "d-flex align-items-center gap-2",
        tags$img(
          src = "www/matrispace_logo.png",
          class = "matrispace-logo",
          alt = "MatriSpace Logo",
          style = "height: 50px; margin-right: 0.25rem;"
        ),
        h1("MatriSpace",
           class = "app-title m-0 me-2",
           style = "font-size: 2.2rem; font-weight: 700;"
        )
      ),
      # Right - About button
      div(
        class = "d-flex align-items-center gap-4",
        actionButton("about", "About",
                     class = "btn btn-outline-primary btn-lg",
                     style = "font-weight: 600; padding: 0.5rem 1.5rem; border-width: 2px;"
        )
      )
    )
  ),

  htmltools::tagAppendAttributes(
  layout_sidebar(
    gap = 0,
    sidebar = sidebar(
      width = "400px",

      # Developed by section
      div(
        class = "developed-by-section mb-2",
        h5("Developed by", style = "margin-bottom: 0.75rem; font-weight: 600; font-size: 1rem; color: var(--color-text-primary);"),
        div(
          class = "d-flex justify-content-center align-items-center gap-3 mb-3",
          tags$a(
            href = "https://www.oulu.fi/en/research-groups/izzi-group",
            target = "_blank",
            tags$img(src = "www/izzilab_logo.png", alt = "Izzi Lab", style = "height: 60px; cursor: pointer;")
          ),
          tags$a(
            href = "https://naba.lab.uic.edu/",
            target = "_blank",
            tags$img(src = "www/naba_lab_logo.png", alt = "Naba Lab for ECM Research", style = "height: 60px; cursor: pointer;")
          )
        ),
        p(
          style = "font-size: 0.85rem; font-style: italic; margin-bottom: 0.75rem;",
          "If you use MatriSpace in your publication, please cite our ",
          tags$span("manuscript", style = "color: var(--color-primary-accent); font-weight: 500;"),
          " (coming soon)"
        ),
        p(
          class = "d-flex align-items-center justify-content-start",
          style = "font-size: 0.9rem; margin-bottom: 0.75rem;",
          "Follow us at ",
          tags$svg(
            style = "width: 18px; height: 18px; margin: 0 0.25rem; vertical-align: text-bottom;",
            viewBox = "0 0 568 501",
            fill = "#1185fe",
            xmlns = "http://www.w3.org/2000/svg",
            tags$path(d = "M123.121 33.6637C188.241 82.5526 258.281 181.681 284 234.873C309.719 181.681 379.759 82.5526 444.879 33.6637C491.866 1.61183 568 -28.9064 568 57.9464C568 75.2916 558.055 203.659 552.222 224.501C531.947 296.954 458.067 315.434 392.347 304.249C507.222 323.8 536.444 388.56 473.333 453.32C353.473 576.312 301.061 422.461 287.631 383.039C285.169 375.812 284.017 372.431 284 375.306C283.983 372.431 282.831 375.812 280.369 383.039C266.939 422.461 214.527 576.312 94.6667 453.32C31.5556 388.56 60.7778 323.8 175.653 304.249C109.933 315.434 36.0535 296.954 15.7778 224.501C9.94525 203.659 0 75.2916 0 57.9464C0 -28.9064 76.1345 1.61183 123.121 33.6637Z")
          ),
          tags$a(
            href = "https://bsky.app/profile/matrisome.bsky.social",
            target = "_blank",
            "The Matrisome Project",
            style = "color: var(--color-primary-accent); margin-left: 0.25rem;"
          )
        ),
        div(
          class = "d-flex gap-2",
          tags$a(
            href = "mailto:matrisomeproject@gmail.com",
            class = "btn btn-outline-primary btn-sm sidebar-action-btn",
            style = "font-size: 0.8rem; padding: 0.2rem 0.5rem; text-decoration: none; white-space: nowrap;",
            bsicons::bs_icon("envelope"), " Contact us"
          ),
          actionButton("guided_tour",
                      tags$span(bsicons::bs_icon("info-circle"), " Take a guided tour of MatriSpace"),
                      class = "btn btn-outline-primary btn-sm sidebar-action-btn",
                      style = "font-size: 0.8rem; padding: 0.2rem 0.5rem; flex: 1; white-space: nowrap;")
        )
      ),

      accordion(
        id = "main_accordion",
        # Data Loading Panel
        accordion_panel(
          "Load Data",
          icon = bsicons::bs_icon("database-add"),

          # File upload UI (no conditionalPanel wrapper)
          fileInput("seurat_file", "Upload Seurat or SpatialExperiment (.rds)",
                    accept = ".rds",
                    buttonLabel = "Browse...",
                    placeholder = "No file selected"),
          div(
            style = "margin-top: 10px;",
            actionButton("load_upload", "Load Object",
                         icon = icon("upload"),
                         class = "btn-primary w-100")
          )

        ),

        # Profile Matrisome Panel
        accordion_panel(
          "Profile Matrisome",
          icon = bsicons::bs_icon("bezier2"),

          # This was moved from Feature Selection
          uiOutput("metadata_selector"),

          # Add a new run button for this section
          div(
            class = "mt-4 d-grid gap-2",
            actionButton(
              "run_matrisome_profile",
              "Run Matrisome Profile",
              icon = icon("play-circle"),
              class = "btn-primary"
            )
          )
        ),

        # Feature Selection Panel
        accordion_panel(
          "Select Feature",
          icon = bsicons::bs_icon("ui-checks-grid"),

          # The metadata_selector is now gone from here.

          # Primary Feature Selection with cards
          div(
            id = "primary-feature-selection",
            class = "mb-3",
            p(class = "mb-2", "Primary Feature Type"),
            div(
              style = "display: flex; gap: 8px;",
              # Hidden radio buttons (NO DEFAULT SELECTION)
              div(
                style = "display: none;",
                radioButtons(
                  "sel1",
                  NULL,
                  choices = c("matrisome gene", "matrisome signature"),
                  selected = character(0) # IMPORTANT: No default
                )
              ),

              # Matrisome Gene Card
              div(
                id = "card_gene",
                class = "feature-card",
                div(
                  class = "text-center p-2",
                  bsicons::bs_icon("diagram-3", size = "1.5rem"),
                  div(class = "small mt-1", "Matrisome gene")
                )
              ),

              # Matrisome Signature Card
              div(
                id = "card_signature",
                class = "feature-card",
                div(
                  class = "text-center p-2",
                  bsicons::bs_icon("grid-3x3", size = "1.5rem"),
                  div(class = "small mt-1", "Matrisome gene list")
                )
              )
            )
          ),

          # The UI for the dropdown will be rendered here by the server
          uiOutput("ui1"),

          # Secondary Feature Selection with cards
          div(
            class = "mb-3",
            p(class = "mb-2", "Secondary Feature Type"),
            div(
              style = "display: flex; flex-wrap: wrap; gap: 8px;",
              # Hidden radio buttons (NO DEFAULT SELECTION)
              div(
                style = "display: none;",
                radioButtons(
                  "sel2",
                  NULL,
                  choices = c("none", "matrisome gene", "any gene", "matrisome signature"),
                  selected = character(0) # IMPORTANT: No default
                )
              ),

              # Card definitions
              div(id = "card_none", class = "feature-card secondary-card", div(class = "text-center p-2", bsicons::bs_icon("x-circle", size = "1.5rem"), div(class = "small mt-1", "None"))),
              div(id = "card_gene2", class = "feature-card secondary-card", div(class = "text-center p-2", bsicons::bs_icon("diagram-3", size = "1.5rem"), div(class = "small mt-1", "Matrisome gene"))),
              div(id = "card_any", class = "feature-card secondary-card", div(class = "text-center p-2", bsicons::bs_icon("search", size = "1.5rem"), div(class = "small mt-1", "Any gene"))),
              div(id = "card_signature2", class = "feature-card secondary-card", div(class = "text-center p-2", bsicons::bs_icon("grid-3x3", size = "1.5rem"), div(class = "small mt-1", "Matrisome gene list")))
            )
          ),

          # The UI for the secondary dropdown will be rendered here by the server
          uiOutput("ui2"),
          # Run Feature Analysis button
          div(
            class = "mt-4 d-grid gap-2",
            actionButton(
              "run_button",
              "Run Feature Analysis",
              icon = icon("chart-line"),
              class = "btn-primary"
            )
          )
        )
      ),

      # Export dropdown - outside accordion, always visible
      div(
        class = "dropdown d-grid",
        tags$button(
          id = "export_dropdown_btn",
          class = "btn btn-success dropdown-toggle initially-disabled",
          type = "button",
          `data-bs-toggle` = "dropdown",
          `aria-expanded` = "false",
          disabled = "disabled",
          icon("download"), " Export Results"
        ),
        tags$ul(
          class = "dropdown-menu w-100",
          `aria-labelledby` = "export_dropdown_btn",
          tags$li(downloadLink("export_all", class = "dropdown-item", icon("file-zipper"), " Export All (ZIP)")),
          tags$li(tags$hr(class = "dropdown-divider")),
          tags$li(downloadLink("export_metadata", class = "dropdown-item", icon("table"), " Spot Metadata (CSV)")),
          tags$li(downloadLink("export_scores", class = "dropdown-item", icon("chart-bar"), " Matrisome Scores (CSV)")),
          tags$li(tags$hr(class = "dropdown-divider")),
          tags$li(downloadLink("export_plots", class = "dropdown-item", icon("images"), " Plots Only (ZIP)"))
        )
      )
    ),

    # Main content area with different pages
    navset_hidden(
      id = "main_content",

      # Data Loading View
      nav_panel(
        "load_data", #Data Overview

        div(class = "mb-4"), # Spacer required for value box rendering

        # This UI is now rendered dynamically by the server
        uiOutput("value_boxes_ui"),


        # Image and Spatial Plots
        layout_column_wrap(
          width = 1/2,

          # Tissue Image Card
          card(
            class = "info-accent-card",
            card_header(
              "Tissue image",
              uiOutput("tissue_link"),
              class = "d-flex align-items-center justify-content-between"
            ),
            card_body(
              plotOutput("tissue_image_plot")
            )
          ),

          # Spatial Plot Card
          card(
            class = "info-accent-card",
            full_screen = TRUE,
            card_header(
              "Cluster map",
              class = "d-flex align-items-center justify-content-between",

            ),
            card_body(
              plotOutput("cluster_plot")
            )
          )
        )
      ),

      # Profile Matrisome View - NEW MAIN PAGE
      nav_panel(
        "profile_matrisome",
        navset_card_tab(
          id = "matrisome_tabs",
          nav_panel(
            "Matrisome profile",
            # --- Main Divisions Card ---
            card(
              class = "section-card",
              card_header(
                class = "d-flex justify-content-between align-items-center",
                span(
                  "Matrisome categories",
                  bslib::tooltip(
                    bs_icon("info-circle", class = "ms-1"),
                    HTML(paste(
                      "Visualises a matrisome gene set using two distinct statistical views. Each view helps answer a different biological question.<br><br>",
                      "&bull; <strong>Spatial Distribution (Relative Activity):</strong> This view uses a <em>Robust Z-score</em> to show how much a spot's expression level deviates from the tissue's typical (median) level. Use this to find regions with statistically up- or down-regulated expression level compared to the tissue's baseline.<br><br>",
                      "&bull; <strong>Hotspots (Absolute Activity):</strong> This view uses a <em>log-transformed score</em> to effectively show the full expression range. Use this to identify spots with the highest expression level without letting extreme values dominate the colour scale."
                    )),
                    options = list(customClass = "matrispace-tooltip")
                  )
                ),
                # Controls wrapper
                div(
                  class = "d-flex align-items-center",
                  style = "gap: 1.5rem;",
                  div(
                    style = "display: flex; align-items: center; height: 38px;",
                    input_switch(
                      id = "toggle_main_divisions",
                      label = "Interactive Plots",
                      value = FALSE
                    )
                  ),
                  div(
                    class = "d-flex align-items-center",
                    style = "gap: 0.5rem; height: 38px;",
                    tags$label(
                      "Spot Size:",
                      `for` = "main_pt_size",
                      style = "margin-bottom: 0; line-height: 38px;"
                    ),
                    div(
                      style = "width: 120px;",
                      sliderInput(
                        "main_pt_size",
                        label = NULL,
                        min = 1, max = 8, value = 2, step = 1, ticks = FALSE
                      )
                    )
                  )
                )
              ),
              card_body(
                layout_column_wrap(
                  width = 1/2,
                  create_matrisome_ui("Glycoproteins", "glycoprotein"),
                  create_matrisome_ui("Collagens", "collagen"),
                  create_matrisome_ui("Proteoglycans", "proteoglycan"),
                  create_matrisome_ui("ECM regulators", "regulator"),
                  create_matrisome_ui("Secreted factors", "secreted_factor"),
                  create_matrisome_ui("ECM-affiliated", "affliated_protein")
                )
              )
            ),

            # --- Matrisome Subcategories Card (Functional Categories) ---
            card(
              class = "section-card",
              card_header(
                class = "d-flex justify-content-between align-items-center",
                span(
                  "Matrisome subcategories",
                  bslib::tooltip(
                    bs_icon("info-circle", class = "ms-1"),
                    HTML(paste(
                      "Visualises a matrisome gene set using two distinct statistical views. Each view helps answer a different biological question.<br><br>",
                      "&bull; <strong>Spatial Distribution (Relative Activity):</strong> This view uses a <em>Robust Z-score</em> to show how much a spot's expression level deviates from the tissue's typical (median) level. Use this to find regions that are statistically up- or down-regulated compared to the tissue's baseline.<br><br>",
                      "&bull; <strong>Hotspots (Absolute Activity):</strong> This view uses a <em>log-transformed score</em> to effectively show the full expression range. Use this to identify spots with the highest expression level without letting extreme values dominate the colour scale."
                    )),
                    options = list(customClass = "matrispace-tooltip")
                  )
                ),
                # INDEPENDENT CONTROLS for subcategories
                div(
                  class = "d-flex align-items-center",
                  style = "gap: 1.5rem;",
                  div(
                    style = "display: flex; align-items: center; height: 38px;",
                    input_switch(
                      id = "toggle_subcategories",
                      label = "Interactive Plots",
                      value = FALSE
                    )
                  ),
                  div(
                    class = "d-flex align-items-center",
                    style = "gap: 0.5rem; height: 38px;",
                    tags$label(
                      "Spot Size:",
                      `for` = "subcategories_pt_size",
                      style = "margin-bottom: 0; line-height: 38px;"
                    ),
                    div(
                      style = "width: 120px;",
                      sliderInput(
                        "subcategories_pt_size",
                        label = NULL,
                        min = 1, max = 8, value = 2, step = 1, ticks = FALSE
                      )
                    )
                  )
                )
              ),
              card_body(
                layout_column_wrap(
                  width = 1/2,
                  create_subdivision_ui("Perivascular", "perivascular", "subcategories"),
                  create_subdivision_ui("Hemostasis", "hemostasis", "subcategories"),
                  create_subdivision_ui("Elastic fibers", "elastic_fibers", "subcategories"),
                  create_subdivision_ui("Growth-factor binding", "gf_binding", "subcategories")
                )
              )
            ),

            # --- Matrisome Gene Families Card ---
            card(
              class = "section-card",
              card_header(
                class = "d-flex justify-content-between align-items-center",
                span(
                  "Matrisome gene families",
                  bslib::tooltip(
                    bs_icon("info-circle", class = "ms-1"),
                    HTML(paste(
                      "Visualises a matrisome gene set using two distinct statistical views. Each view helps answer a different biological question.<br><br>",
                      "&bull; <strong>Spatial Distribution (Relative Activity):</strong> This view uses a <em>Robust Z-score</em> to show how much a spot's expression level deviates from the tissue's typical (median) level. Use this to find regions that are statistically up- or down-regulated compared to the tissue's baseline.<br><br>",
                      "&bull; <strong>Hotspots (Absolute Activity):</strong> This view uses a <em>log-transformed score</em> to effectively show the full expression range. Use this to identify spots with the highest expression level without letting extreme values dominate the colour scale."
                    )),
                    options = list(customClass = "matrispace-tooltip")
                  )
                ),
                # INDEPENDENT CONTROLS for gene families
                div(
                  class = "d-flex align-items-center",
                  style = "gap: 1.5rem;",
                  div(
                    style = "display: flex; align-items: center; height: 38px;",
                    input_switch(
                      id = "toggle_families",
                      label = "Interactive Plots",
                      value = FALSE
                    )
                  ),
                  div(
                    class = "d-flex align-items-center",
                    style = "gap: 0.5rem; height: 38px;",
                    tags$label(
                      "Spot Size:",
                      `for` = "families_pt_size",
                      style = "margin-bottom: 0; line-height: 38px;"
                    ),
                    div(
                      style = "width: 120px;",
                      sliderInput(
                        "families_pt_size",
                        label = NULL,
                        min = 1, max = 8, value = 2, step = 1, ticks = FALSE
                      )
                    )
                  )
                )
              ),
              card_body(
                layout_column_wrap(
                  width = 1/2,
                  create_subdivision_ui("Laminins", "laminin", "families"),
                  create_subdivision_ui("Matricellular proteins", "matricellular", "families"),
                  create_subdivision_ui("Syndecans", "syndecan", "families"),
                  create_subdivision_ui("Glypicans", "glypican", "families")
                )
              )
            )
          ), # end of Matrisome Profile nav_panel

          nav_panel(
            "ECM niches",
            interactive_spatial_viewer_ui(
              id_prefix = "ecm",
              title = "ECM niche annotations",
              show_sidebar = FALSE,
              tooltip_text = paste(
                "An interactive map classifying each spot into a dominant ECM niche based on marker gene expression.<br><br>",
                "This provides a functional, high-level annotation of the tissue's microenvironment, revealing the underlying ECM architecture."
              )
            ),
            card(
              class = "section-card",
              card_header(
                class = "d-flex justify-content-between align-items-center",
                span(
                  "Niche score",
                  bslib::tooltip(
                    bs_icon("info-circle", class = "ms-1"),
                    HTML(paste(
                      "A spatial plot visualising the activity of a specific ECM niche (e.g., Interstitial). The score is calculated from the collective expression of each gene from the gene set in each spot.<br><br>",
                      "Unlike the discrete 'ECM niche annotations' map, this continuous score can reveal gradients of signature expression levels and areas where different ECM gene expression programmes may overlap."
                    )),
                    options = list(customClass = "matrispace-tooltip")
                  )
                ),
                div(
                  class = "d-flex align-items-center",
                  style = "gap: 1.5rem;",
                  div(
                    style = "display: flex; align-items: center; height: 38px;",
                    input_switch(
                      id = "toggle_ecm_interactive",
                      label = "Interactive Plots",
                      value = FALSE
                    )
                  ),
                  div(
                    class = "d-flex align-items-center",
                    style = "gap: 0.5rem; height: 38px;",
                    tags$label(
                      "Spot Size:",
                      `for` = "ecm_pt_size",
                      style = "margin-bottom: 0; line-height: 38px;"
                    ),
                    div(
                      style = "width: 120px;",
                      sliderInput(
                        "ecm_pt_size",
                        label = NULL,
                        min = 1, max = 8, value = 3, step = 1, ticks = FALSE
                      )
                    )
                  )
                )
              ),
              card_body(
                layout_column_wrap(
                  width = 1/2,
                  create_ecm_sig_ui("Interstitial ECM", "interstitial"),
                  create_ecm_sig_ui("Basement membrane", "basement")
                )
              )
            ),

            # --- ECM Signature Distribution Plot ---
            card(
              class = "section-card",
              full_screen = TRUE,
              card_header(
                span(
                  "Niche distribution",
                  bslib::tooltip(
                    bs_icon("info-circle", class = "ms-1"),
                    HTML(paste(
                      "This plot summarizes the average expression level of each ECM niche across all annotation groups.<br><br>",
                      "This feature can be used to quantitatively identify which cell types or tissue regions are the primary contributors to the interstitial ECM or the basement membrane programs."
                    )),
                    options = list(customClass = "matrispace-tooltip")
                  )
                )
              ),
              card_body(
                plotOutput("ecm_signature_distribution_plot", height = "500px")
              )
            ),
            card(
              class = "section-card",
              full_screen = TRUE,
              card_header(
                span(
                  "Spatial LR co-expression analysis",
                  bslib::tooltip(
                    bsicons::bs_icon("info-circle", class = "ms-1"),
                    HTML(paste(
                      "This analysis scores potential signaling pathways by measuring the co-expression of an ECM ligand in one spot and its corresponding receptor in neighboring spots.<br><br>",
                      "&bull; <b>Heatmap:</b> Displays the top enriched LR 'axes' (ligand-receptor pairs) across all annotated niches, providing a high-level overview of the most active communication patterns.<br>",
                      "&bull; <b>Volcano Plot:</b> Provides a detailed view of all LR pairs for a selected niche, highlighting which pairs are the most significantly enriched in that specific environment."
                    )),
                    options = list(customClass = "matrispace-tooltip")
                  )
                )
              ),
              card_body(
                layout_sidebar(
                  sidebar = sidebar(
                    width = "250px",
                    h5("Plot options"),

                    # Dynamic dropdown: populated by server.R output$lr_cluster_selector_ui
                    # Filters volcano plot to show enrichment for selected ECM niche
                    uiOutput("lr_cluster_selector_ui"),

                    # Warning alert: informs user of computation time
                    div(
                      class = "alert alert-warning text-center p-2 mt-3 lr-warning-alert",
                      tags$strong(bsicons::bs_icon("hourglass-split"), " Computation Time"),
                      p(
                        class = "mt-1 mb-0",
                        "This analysis is intensive and may take several minutes. The app will be blocked while running."
                      )
                    ),

                    # Analysis trigger button
                    # Disabled during analysis via shinyjs (see server.R observeEvent)
                    actionButton(
                      "run_lr_analysis",
                      "Calculate LR Scores",
                      icon = icon("cogs"),
                      class = "btn-primary w-100 mt-3"
                    )
                  ),

                  # Main plot area: side-by-side heatmap and volcano plot
                  layout_column_wrap(
                    width = 1/2,  # Equal width columns

                    # Left: Heatmap showing top enriched L-R axes across ALL niches
                    card(
                      card_header("Top enriched LR axes"),
                      card_body(
                        shinycssloaders::withSpinner(  # Loading indicator during rendering
                          plotOutput("lr_heatmap_plot", height = "600px")
                        )
                      )
                    ),

                    # Right: Volcano plot for SELECTED niche (from dropdown)
                    card(
                      card_header("Interaction volcano plot"),
                      card_body(
                        shinycssloaders::withSpinner(
                          plotOutput("lr_volcano_plot", height = "600px")
                        )
                      )
                    )
                  )
                )
              )
            )
          ) # end of ECM Niches nav_panel
        )
      ),

      # Feature Selection View - with its own tabs
      nav_panel(
        "feature_selection",
        navset_card_tab(
          id = "feature_tabs",
          # Only the "Feature Analysis" panel remains
          nav_panel(
            "Feature analysis",
            interactive_spatial_viewer_ui(
              id_prefix = "main",
              title = "Tissue overview",
              show_sidebar = TRUE,
              tooltip_text = paste(
                "An interactive, overlaid view of the H&E tissue image and its corresponding annotations.<br><br>",
                "Use the Annotation Selector to filter spots. Hover over the image to activate a magnifier, and hover over spots for details."
              )
            ),

            card(
              full_screen = TRUE,
              card_header(
                class = "d-flex justify-content-between align-items-center",
                span(
                  "Feature expression",
                  bslib::tooltip(
                    bs_icon("info-circle", class = "ms-1"),
                    HTML(paste(
                      "A quantitative comparison of selected feature's expression levels across the different annotated groups.<br><br>",
                      "Each violin or box visualises the distribution of expression values within a group, helping you identify which groups are the main source of the feature. The ANOVA p-value tests for statistically significant differences between these groups."
                    )),
                    options = list(customClass = "matrispace-tooltip")
                  )
                ),
                # Main controls container for pills
                div(
                  class = "d-flex align-items-center",
                  style = "gap: 0.75rem;", # Standard gap between pills

                  # Pill 1: Violin Plot toggle
                  div(
                    input_switch("violin_plot", "Violin Plot", FALSE)
                  ),

                  # Pill 2: Independent Y-Axis toggle
                  div(
                    input_switch("independent_y", "Independent Y-Axis", FALSE)
                  )
                )
              ),
              card_body(
                height = "500px",
                plotlyOutput("primary_expression_plot"),
                plotlyOutput("secondary_expression_plot")
              )
            ),

            # Feature Plots
            layout_column_wrap(
              width = 1/2,

              # Feature Plot 1
              card(
                full_screen = FALSE,
                card_header(
                  span(
                    "Primary feature",
                    bslib::tooltip(
                      bs_icon("info-circle", class = "ms-1"),
                      HTML(paste(
                        "A spatial plot that pinpoints exactly where your selected feature is expressed in the tissue.<br><br>",
                        "Spots are coloured from blue (low expression) to red (high expression), allowing the visualization of potential correlation in feature expression patterns across different regions. Hover over spots to see their precise value."
                      )),
                      options = list(customClass = "matrispace-tooltip")
                    )
                  ),
                  class = "d-flex align-items-center justify-content-between",
                ),
                card_body(
                  style = "min-height: 500px;",
                  plotlyOutput("primary_feature_plot", height = "100%")
                ),
                card_footer(
                  # Use htmlOutput for better performance
                  htmlOutput("primary_feature_autocorrelation_footer")
                )
              ),

              # Feature Plot 2
              card(
                full_screen = FALSE,
                card_header(
                  span(
                    "Secondary feature",
                    bslib::tooltip(
                      bs_icon("info-circle", class = "ms-1"),
                      HTML(paste(
                        "A spatial plot that pinpoints exactly where your selected feature is active in the tissue.<br><br>",
                        "Spots are coloured from blue (low expression) to red (high expression), allowing the visualization of potential correlation in feature expression patterns across different regions. Hover over spots to see their precise value."
                      )),
                      options = list(customClass = "matrispace-tooltip")
                    )
                  ),
                  class = "d-flex align-items-center justify-content-between",
                ),
                card_body(
                  style = "min-height: 500px;",
                  plotlyOutput("secondary_feature_plot", height = "100%")
                ),
                card_footer(
                  # Use htmlOutput for better performance
                  htmlOutput("secondary_feature_autocorrelation_footer")
                )
              )
            ),

            # Co-expression Plot
            card(
              full_screen = TRUE,
              card_header(
                class = "d-flex justify-content-between align-items-center",
                span(
                  "Co-expression",
                  bslib::tooltip(
                    bs_icon("info-circle", class = "ms-1"),
                    HTML(paste(
                      "A blended visualisation for identifying areas where two features are co-expressed in the same spot, which can suggest functional association.<br><br>",
                      "The colour of each spot is a mix of the two feature colours, as shown in the legend. Bright, mixed colours indicate high expression of both features, highlighting potential 'hotspots' of co-expression."
                    )),
                    options = list(customClass = "matrispace-tooltip")
                  )
                ),
                # Plot settings popover
                bslib::popover(
                  title = "Plot Settings",
                  # The button that triggers the popover
                  trigger = actionButton(
                    "coex_settings_btn",
                    label = NULL,
                    icon = icon("gear"),
                    class = "btn-sm btn-light"
                  ),
                  # Controls inside the popover
                  radioButtons(
                    inputId = "coex_palette",
                    label = "Color Palette:",
                    choices = c("Classic", "Vibrant", "Microscopy"),
                    selected = "Vibrant",
                    inline = TRUE
                  ),
                  sliderInput(
                    inputId = "coex_alpha",
                    label = "Opacity:",
                    ticks = FALSE,
                    min = 0.1, max = 1.0, value = 0.7, step = 0.1
                  ),
                  sliderInput(
                    inputId = "coex_pt_size",
                    label = "Point Size:",
                    ticks = FALSE,
                    min = 1, max = 8, value = 3, step = 0.5
                  )
                )
              ),
              card_body(
                style = "min-height: 500px;",
                # Use plotOutput for ggplot objects
                plotOutput("coexpression_plot", height = "500px")
              ),
              card_footer(
                htmlOutput("coexpression_crosscor_footer")
              )
            ),

            card(
              full_screen = TRUE,
              card_header(
                span(
                  "Spatial association",
                  bslib::tooltip(
                    bs_icon("info-circle", class = "ms-1"),
                    HTML(paste(
                      "Analyses the relationship between a feature's expression at a specific spot and the expression of another feature in its immediate neighborhood.<br><br>",
                      "<b>Spatial Plot:</b> Maps the four types of association across the tissue.<br>",
                      "<b>Composition Plot:</b> Shows the proportion of cell annotations that make up each association type, sorted from most to least common."
                    )),
                    options = list(customClass = "matrispace-tooltip")
                  )
                )
              ),
              card_body(
                style = "min-height: 500px;",
                plotOutput("lisa_clustering_plot", height = "500px")
              )
            )
          )
        )
      ),

      nav_panel(
        "about",
        card(
          card_header("About MatriSpace"),
          card_body(
            class = "about-content",

            h4("Funding support"),
            tags$ul(
              tags$li(strong("A. Naba:"), " NIH/NHGRI Human Biomolecular Atlas Program (HuBMAP) - Grant U01HG012680"),
              tags$li(strong("A. Naba:"), " NIH/NCI Innovative Molecular Analysis Technologies Program - Grant R21CA261642"),
              tags$li(strong("V. Izzi:"), " Cancer Foundation Finland (2023-2024)"),
              tags$li(strong("V. Izzi:"), " DigiHealth-project, a strategic profiling project at the University of Oulu, and the Infotech Institute of the University of Oulu"),
              tags$li(strong("V. Izzi:"), " The European Union CARES project")
            ),

            h4("Data resources"),
            p(strong("MatriSpace is powered by the following databases:"),
            tags$ul(
              tags$li("The matrisome gene lists are available via The Matrisome Project: ",
                      tags$a(href = "https://matrisome.org", target = "_blank", "https://matrisome.org")),
              tags$li("Experimental ECM gene sets are available via the Molecular Signatures Database (MSigDB): ",
                      tags$a(href = "https://www.gsea-msigdb.org/gsea/msigdb/index.jsp", target = "_blank", "https://www.gsea-msigdb.org/gsea/msigdb/index.jsp")),
              tags$li("Matrisome gene subcategories and families are sourced from Gene Ontology: ",
                      tags$a(href = "https://geneontology.org/", target = "_blank", "https://geneontology.org/"))
            )),

            p(strong("Open-access (OA) data included in the online version of MatriSpace are sourced from:")),
            tags$ul(
              tags$li(tags$a(href = "https://www.10xgenomics.com/", target = "_blank", "10X Genomics website")),
              tags$li(tags$a(href = "https://humantumoratlas.org/", target = "_blank", "HTAN (Human Tumor Atlas Network)")),
              tags$li(tags$a(href = "https://www.ncbi.nlm.nih.gov/geo/", target = "_blank", "GEO (Gene Expression Omnibus, NCBI)"))
            ),

            h4("Source code"),
            p("The MatriSpace Shiny app is available at ",
              tags$a(href = "https://github.com/izzilab/matrispace-app", target = "_blank", "https://github.com/izzilab/matrispace-app")),
            p("The MatriSpace R package is available at ",
              tags$a(href = "https://github.com/izzilab/matrispace", target = "_blank", "https://github.com/izzilab/matrispace")),

            h4("Development & Technical help"),
            p("MatriSpace is developed and maintained by the Naba Lab for ECM Research at the University of Illinois Chicago and the Izzi Lab at Oulu University.
              You can reach us via the \"Contact Us\" button on the home page of the app")
          )
        )
      )
    )
  ),
  class = "app-scale-target"
  )
)

}
