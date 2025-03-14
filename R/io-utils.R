#' Warn that a file was not overwritten
#'
#' @param path A character path to an existing file or directory
#'
#' @return (Invisible) `NULL`
warn_no_overwrite <- function(path) {
  msg <- glue::glue("file already exists at `{path}`, and `overwrite = FALSE`, so `{path}` not overwritten")

  rlang::warn(message = msg, class = "warning_no_overwrite")

  return(invisible())
}


#' Create directory to path
#'
#' The goal is to mimic the behavior of `mkdir -p` command. That is, create the
#' full directory path recursively without any errors or warnings if the
#' directory already exists.
#'
#' @param path A character path to a directory
#'
#' @return A logical indicating if the directory creation succeeded
#' @export
mkdir_p <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
}
