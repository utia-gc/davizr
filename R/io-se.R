#' Write a `SummarizedExperiment` object to file as an RDS.
#'
#' This is essentially a thin wrapper around `saveRDS()` with additional logic
#' for things like directory creation and control of overwriting files.
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
  # if file already exists at path and overwrite turned off print a warning and return
  if (file.exists(path) & !overwrite) {
    warn_no_overwrite(path)
    return(invisible(FALSE))
  }

  # create the base directory if it doesn't already exist
  mkdir_p(dirname(path))

  # write the SummarizedExperiment as an RDS file
  saveRDS(object = se, file = path)

  return(invisible(TRUE))
}
