#' Invoke a golem JavaScript handler
#'
#' @param name The name of the handler
#' @param ... other params passed to session$sendCustomMessage
#' @noRd
invoke_js <- function(name, ...) {
  session <- shiny::getDefaultReactiveDomain()
  session$sendCustomMessage(type = name, message = list(...))
}
