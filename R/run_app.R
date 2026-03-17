#' Run the MatriSpace Application
#'
#' @param host Host address (default "127.0.0.1" for local, "0.0.0.0" for Docker)
#' @param port Port number (default NULL lets Shiny choose)
#' @param ... Additional arguments passed to golem_opts.
#'   See `?golem::get_golem_options` for more details.
#' @export
run_app <- function(host = "127.0.0.1", port = NULL, ...) {
  # Set 15GB upload limit
  options(shiny.maxRequestSize = 15000 * 1024^2)

  with_golem_options(
    app = shinyApp(
      ui = app_ui,
      server = app_server,
      options = list(host = host, port = port)
    ),
    golem_opts = list(...)
  )
}
