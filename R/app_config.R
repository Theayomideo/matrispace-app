#' Access files in the current app
#'
#' @param ... character vectors, specifying subdirectory and file(s)
#'   within your package. The default, none, returns the root of the app.
#' @noRd
app_sys <- function(...) {
  system.file(..., package = "matrispace.app")
}

#' Access files in inst/extdata
#'
#' Convenience wrapper for accessing bundled reference data files.
#'
#' @param ... character vectors, specifying file path within inst/extdata
#' @return Full system path to the file
#' @export
extdata_path <- function(...) {
  app_sys("extdata", ...)
}

#' Read App Config
#'
#' @param value Value to retrieve from the config file.
#' @param config R_CONFIG_ACTIVE value.
#' @param use_parent Logical, scan the parent directory for config file.
#' @noRd
get_golem_config <- function(
    value,
    config = Sys.getenv("R_CONFIG_ACTIVE", "default"),
    use_parent = TRUE
) {
  config::get(
    value = value,
    config = config,
    file = app_sys("golem-config.yml"),
    use_parent = use_parent
  )
}
