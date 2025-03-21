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


#' Set a modeling variable in the metadata of a `SummarizedExperiment`
#'
#' @details The character vectors in `values` must all be names of variables
#' available in the sample data. Otherwise, an `error_no_sample_variable` is
#' thrown.
#'
#' @inheritParams write_se
#' @param variable A character vector for which variable to set
#' @param values A character vector of which variables to set the modeling
#'   variable as
#'
#' @return The `SummarizedExperiment` with the modeling variable set to values
#'   in the metadata slot
#' @export
set_modeling_variable <- function(se, variable = c("interest", "nuisance", "screening"), values) {

  variable <- rlang::arg_match(variable)

  # throw error if attempt to supply a value that isn't a column in sample data
  sample_data_variables <- colnames(get_samples_data(se))
  if (!all(values %in% sample_data_variables)) {
    invalid_sample_variables <- values[!(values %in% sample_data_variables)]
    abort_no_sample_variable(var = invalid_sample_variables, se = "se")
  }

  S4Vectors::metadata(se)[["modeling"]][["variables"]][[variable]] <- values

  return(se)
}


#' Get a modeling variable in the metadata of a `SummarizedExperiment`
#'
#' @inheritParams set_modeling_variable
#'
#' @return A character vector of modeling variables.
#' @export
get_modeling_variable <- function(se, variable = c("interest", "nuisance", "screening")) {

  variable <- rlang::arg_match(variable)

  values <- S4Vectors::metadata(se)[["modeling"]][["variables"]][[variable]]

  return(values)
}


#' Get all modeling variables from a `SummarizedExperiment`
#'
#' @inheritParams set_modeling_variable
#'
#' @return A character vector of all modeling variables.
#' @export
get_modeling_all_variables <- function(se) {

  values <- Reduce(
    union,
    c(
      get_modeling_variable(se, variable = "interest"),
      get_modeling_variable(se, variable = "nuisance"),
      get_modeling_variable(se, variable = "screening")
    )
  )

  return(values)
}


#' Get explanatory modeling variables from a `SummarizedExperiment`
#'
#' @description Explanatory variables are the union of the interest and nuisance
#' variables.
#'
#' @inheritParams set_modeling_variable
#'
#' @return A character vector of explanatory variables.
#' @export
get_modeling_explanatory_variables <- function(se) {

  values <- union(
    get_modeling_variable(se, variable = "interest"),
    get_modeling_variable(se, variable = "nuisance")
  )

  return(values)
}


#' Throw an error when `SummarizedExperiment` has no modeling variables
#'
#' @param se The `SummarizedExperiment` that does not modeling variables
#'   variable
#'
#' @return An `rlang_error` of custom class `error_no_modeling_variables`
#' @export
abort_no_modeling_variables <- function(se) {
  msg <- c(
    "Must have modeling variables in {.var {se}}."
  )

  cli::cli_abort(
    message = msg,
    class = "error_no_modeling_variables"
  )
}
