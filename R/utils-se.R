#' Get samples data (column data) from a `SummarizedExperiment`
#'
#' @inheritParams write_se
#' @param rownames_var A character vector for the name of the variable to store
#'   row names as. If `NULL`, row names stay in the row names and no new
#'   variable in the data table is created.
#'
#' @return A data.frame of samples data.
#' @export
get_samples_data <- function(se, rownames_var = "sample") {
  samples_data <- se |>
    SummarizedExperiment::colData() |>
    as.data.frame()

  if (!is.null(rownames_var)) {
    samples_data <- tibble::rownames_to_column(samples_data, var = rownames_var)
  }

  return(samples_data)
}


#' Throw an error for a variable that does not exist in the sample data of a
#' `SummarizedExperiment`
#'
#' @param var The sample data variable that does not exist in the
#'   `SummarizedExperiment`
#' @param se The `SummarizedExperiment` that does not contain the sample data
#'   variable
#'
#' @return An `rlang_error` of custom class `error_no_sample_variable`
abort_no_sample_variable <- function(var, se) {
  msg <- c(
    "Can't find variable {.var {var}} in sample data of {.var {se}}."
  )

  cli::cli_abort(
    message = msg,
    class = "error_no_sample_variable"
  )
}
