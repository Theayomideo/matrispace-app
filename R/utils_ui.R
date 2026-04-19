#' Create Interactive Spatial Viewer UI
#'
#' Reusable D3.js spatial viewer with zoom/magnifier, annotation filtering,
#' and spot controls. Multi-instance support via id_prefix.
#'
#' @param id_prefix Prefix for input/output IDs
#' @param title Card title
#' @param show_sidebar Whether to show sidebar controls
#' @param tooltip_text HTML tooltip text
#' @return A bslib card UI element
#' @noRd
interactive_spatial_viewer_ui <- function(id_prefix, title, show_sidebar = TRUE, tooltip_text) {
  ns_id <- function(name) paste0(name, "_", id_prefix)

  viz_content <- div(
    class = "d3-viz-instance-container",
    `data-prefix` = id_prefix,
    div(
      class = "viz-and-zoom-layout",
      div(id = ns_id("d3_viz_container"), div(id = ns_id("centerViz"), tags$svg(id = ns_id("d3_svg")))),
      div(id = ns_id("zoom-panel"), div(class = "zoom_container hidden",
        tags$span(id = ns_id("verticalLine"), class = "verticalLine"),
        tags$span(id = ns_id("horizontalLine"), class = "horizontalLine"),
        tags$img(id = ns_id("tissue_original"), class = "original", src = "")
      ))
    )
  )

  card(
    class = "section-card",
    full_screen = TRUE,
    card_header(
      class = "d-flex align-items-center justify-content-between",
      span(
        title,
        bslib::tooltip(
          bs_icon("info-circle", class = "ms-1"),
          HTML(tooltip_text),
          options = list(customClass = "matrispace-tooltip")
        )
      ),
      bslib::popover(
        title = "Visual Settings",
        trigger = actionButton(
          ns_id("viz_settings_btn"),
          label = NULL, icon = icon("gear"), class = "btn-sm btn-light"
        ),
        sliderInput(ns_id("spot_opacity"), "Spot Opacity:", min = 0, max = 1, value = 0.8, step = 0.1, ticks = FALSE),
        sliderInput(ns_id("spot_size"), "Spot Size:", min = 1, max = 8, value = 3.5, step = 0.5, ticks = FALSE),
        input_switch(ns_id("spot_borders"), "Show Borders", value = TRUE)
      )
    ),

    if (show_sidebar) {
      layout_sidebar(
        sidebar = sidebar(
          width = "275px",
          card_title("Annotation selector:"),
          div(
            class = "d-flex justify-content-start gap-2 mb-3",
            actionButton(ns_id("selectAllAnnotations"), "Select All", class = "btn-outline-primary btn-sm"),
            actionButton(ns_id("deselectAllAnnotations"), "Deselect All", class = "btn-outline-secondary btn-sm")
          ),
          uiOutput(ns_id("cellTypeSelector"))
        ),
        viz_content
      )
    } else {
      div(
        class = "d-flex align-items-start",
        style = "gap: 1rem; padding: 1rem;",
        div(
          style = "flex: 0 0 200px;",
          plotOutput(ns_id("legend_plot"), height = "250px")
        ),
        div(
          style = "flex: 1 1 auto;",
          viz_content
        )
      )
    }
  )
}

#' Create Matrisome Category UI Card
#'
#' @param title Display title
#' @param name_prefix Internal name prefix for plot IDs
#' @return A navset_card_tab UI element
#' @noRd
create_matrisome_ui <- function(title, name_prefix) {
  navset_card_tab(
    full_screen = FALSE, title = title,
    nav_panel(
      tagList(
        "Spatial distribution",
        bslib::tooltip(
          bs_icon("info-circle", class = "ms-1"),
          "Allows the visualisation of gene expression level relative to the tissue's baseline using a Z-score, highlighting statistically significant regional differences. Ideal for identifying areas of relative up- or down-regulation.",
          options = list(customClass = "matrispace-tooltip")
        )
      ),
      conditionalPanel("input.toggle_main_divisions == false", plotOutput(paste0(name_prefix, "_dist_static"))),
      conditionalPanel("input.toggle_main_divisions == true", plotlyOutput(paste0(name_prefix, "_dist_plotly")))
    ),
    nav_panel(
      tagList(
        "Hotspots",
        bslib::tooltip(
          bs_icon("info-circle", class = "ms-1"),
          "Highlights spots with the highest absolute expression using a log-scaled score, allowing the visualisation of the full signal range. Ideal for pinpointing primary sources of a gene set.",
          options = list(customClass = "matrispace-tooltip")
        )
      ),
      conditionalPanel("input.toggle_main_divisions == false", plotOutput(paste0(name_prefix, "_hotspot_static"))),
      conditionalPanel("input.toggle_main_divisions == true", plotlyOutput(paste0(name_prefix, "_hotspot_plotly")))
    ),
    footer = card_footer(
      uiOutput(paste0(name_prefix, "_autocorrelation_footer"))
    )
  )
}

#' Create Subdivision/Family UI Card
#'
#' @param title Display title
#' @param name_prefix Internal name prefix for plot IDs
#' @param toggle_prefix Toggle input prefix ("subcategories" or "families")
#' @return A navset_card_tab UI element
#' @noRd
create_subdivision_ui <- function(title, name_prefix, toggle_prefix) {
  navset_card_tab(
    full_screen = FALSE, title = title,
    nav_panel(
      tagList(
        "Spatial distribution",
        bslib::tooltip(
          bs_icon("info-circle", class = "ms-1"),
          "Allows the visualisation of gene expression level relative to the tissue's baseline using a Z-score, highlighting statistically significant regional differences. Ideal for identifying areas of relative up- or down-regulation.",
          options = list(customClass = "matrispace-tooltip")
        )
      ),
      conditionalPanel(paste0("input.toggle_", toggle_prefix, " == false"), plotOutput(paste0(name_prefix, "_dist_static"))),
      conditionalPanel(paste0("input.toggle_", toggle_prefix, " == true"), plotlyOutput(paste0(name_prefix, "_dist_plotly")))
    ),
    nav_panel(
      tagList(
        "Hotspots",
        bslib::tooltip(
          bs_icon("info-circle", class = "ms-1"),
          "Highlights spots with the highest absolute expression using a log-scaled score, allowing the visualisation of the full signal range. Ideal for pinpointing primary sources of a gene set.",
          options = list(customClass = "matrispace-tooltip")
        )
      ),
      conditionalPanel(paste0("input.toggle_", toggle_prefix, " == false"), plotOutput(paste0(name_prefix, "_hotspot_static"))),
      conditionalPanel(paste0("input.toggle_", toggle_prefix, " == true"), plotlyOutput(paste0(name_prefix, "_hotspot_plotly")))
    ),
    footer = card_footer(
      uiOutput(paste0(name_prefix, "_autocorrelation_footer"))
    )
  )
}

#' Create ECM Niche Signature UI Card
#'
#' @param title Display title
#' @param name_prefix Internal name prefix for plot IDs
#' @return A navset_card_tab UI element
#' @noRd
create_ecm_sig_ui <- function(title, name_prefix) {
  navset_card_tab(
    full_screen = FALSE, title = title,
    nav_panel(
      tagList(
        "Spatial distribution",
        bslib::tooltip(
          bs_icon("info-circle", class = "ms-1"),
          "Allows the visualisation of gene expression level relative to the tissue's baseline using a Z-score, highlighting statistically significant regional differences. Ideal for identifying areas of relative up- or down-regulation.",
          options = list(customClass = "matrispace-tooltip")
        )
      ),
      conditionalPanel("input.toggle_ecm_interactive == false", plotOutput(paste0(name_prefix, "_dist_static"))),
      conditionalPanel("input.toggle_ecm_interactive == true", plotlyOutput(paste0(name_prefix, "_dist_plotly")))
    ),
    nav_panel(
      tagList(
        "Hotspots",
        bslib::tooltip(
          bs_icon("info-circle", class = "ms-1"),
          "Highlights spots with the highest absolute expression using a log-scaled score, allowing the visualisation of the full signal range. Ideal for pinpointing primary sources of a gene set.",
          options = list(customClass = "matrispace-tooltip")
        )
      ),
      conditionalPanel("input.toggle_ecm_interactive == false", plotOutput(paste0(name_prefix, "_hotspot_static"))),
      conditionalPanel("input.toggle_ecm_interactive == true", plotlyOutput(paste0(name_prefix, "_hotspot_plotly")))
    ),
    footer = card_footer(
      htmlOutput(paste0(name_prefix, "_autocorrelation_footer"))
    )
  )
}

#' Create Flippable Value Box
#'
#' Interactive metadata cards with front/back via CSS 3D transform.
#'
#' @param front_ui Front-side UI content
#' @param back_ui Back-side UI content
#' @return A div with flip-card structure
#' @noRd
flippable_value_box <- function(front_ui, back_ui) {
  tags$div(
    class = "flip-card",
    tags$div(
      class = "flip-card-inner",
      tags$div(class = "flip-card-front", front_ui),
      tags$div(class = "flip-card-back", back_ui)
    )
  )
}
