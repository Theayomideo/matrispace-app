#' Turn an R list into an HTML list
#'
#' @param list An R list
#' @param class a CSS class for the list
#' @return an HTML list
#' @noRd
list_to_li <- function(list, class = NULL) {
  if (is.null(class)) {
    tagList(lapply(list, tags$li))
  } else {
    res <- lapply(list, tags$li)
    tagList(lapply(res, function(x) tagAppendAttributes(x, class = class)))
  }
}

#' @noRd
jq_hide <- function(id) {
  tags$script(sprintf("$(function(){$('#%s').hide();})", id))
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path("www", app_sys("app/www"))

  tags$head(
    # Bootstrap Icons CDN
    tags$link(rel = "stylesheet", href = "https://cdn.jsdelivr.net/npm/bootstrap-icons@1.11.3/font/bootstrap-icons.min.css"),

    # Google Fonts - Inter font family
    tags$link(rel = "preconnect", href = "https://fonts.googleapis.com"),
    tags$link(rel = "preconnect", href = "https://fonts.gstatic.com", crossorigin = NA),
    tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap"),

    # App stylesheet
    tags$link(rel = "stylesheet", type = "text/css", href = "www/styles.css"),

    # D3.js library
    tags$script(src = "www/d3.v4.js"),

    # Custom JavaScript
    tags$script(src = "www/custom.js")
  )
}
