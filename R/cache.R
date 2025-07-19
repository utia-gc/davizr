# cache.R
# Utilities related to caching results


#' Wrapper around [Hmisc::runifChanged()].
#'
#' @description Wraps [Hmisc::runifChanged()]. Introduces convenience whereby
#'   the directory for `file` is created if it doesn't already exist.
#'
#' @param file Path to file to store cached results. `dirname(file)` is created
#'   if it does not already exist.
#' @param fun The usually slow function to run.
#' @param ... Additional arguments and objects passed to
#'   [Hmisc::runifChanged()].
#'
#' @returns Result of `fun`.
#' @export
run_if_changed <- function(file, fun, ...) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)

  Hmisc::runifChanged(fun = fun, ..., file = file)
}
