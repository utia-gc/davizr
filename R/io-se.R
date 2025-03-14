#' Write a `SummarizedExperiment` object to file as an RDS.
#'
#' @param se A `SummarizedExperiment`
#' @param path A character path to write the `SummarizedExperiment` file to
#' @param overwrite A logical controlling whether a file should be overwritten
#'   if `path` already exists
#'
#' @return (Invisible) a logical indicating whether the `se` was successfully
#'   written to `path`
#' @export
#'
#' @examples
#' # don't run
#' if (FALSE) {
#'   # create SummarizedExperiment object
#'   counts <- random_matrix()
#'   col_data <- random_col_data()
#'   se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(counts = counts),
#'     colData = col_data
#'   )
#'
#'   # write the SummarizedExperiment
#'   write_se(se, path = "se.rds")
#' }
write_se <- function(se, path, overwrite = FALSE) {
  if (file.exists(path) & !overwrite) {
    warn_no_overwrite(path)
    return(invisible(FALSE))
  }

  saveRDS(object = se, file = path)

  return(invisible(TRUE))
}


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
