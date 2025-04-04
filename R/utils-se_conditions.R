# utils-se_conditions.R
# Signal conditions associated with SummarizedExperiment objects


#' Throw an error for a variable that does not exist in the sample data of a
#' `SummarizedExperiment`
#'
#' @param var The sample data variable that does not exist in the
#'   `SummarizedExperiment`
#' @param se The `SummarizedExperiment` that does not contain the sample data
#'   variable
#'
#' @return An `rlang_error` of custom class `error_no_sample_variable`
#' @noRd
abort_no_sample_variable <- function(var, se) {
  msg <- c(
    "Can't find variable {.var {var}} in sample data of {.var {se}}."
  )

  cli::cli_abort(
    message = msg,
    class = "error_no_sample_variable"
  )
}


#' Throw an error when `SummarizedExperiment` has no modeling variables
#'
#' @param se The `SummarizedExperiment` that does not modeling variables
#'   variable
#'
#' @return An `rlang_error` of custom class `error_no_modeling_variables`
#' @noRd
abort_no_modeling_variables <- function(se) {
  msg <- c(
    "Must have modeling variables in {.var {se}}."
  )

  cli::cli_abort(
    message = msg,
    class = "error_no_modeling_variables"
  )
}


#' Throw an error when `SummarizedExperiment` does not contain a requested assay
#'
#' @param assay A character assay name that is not in the assay names of the
#'   `SummarizedExperiment`
#' @param se The `SummarizedExperiment` that does not contain the requested
#'   `assay`
#'
#' @return An `rlang_error` of custom class `error_no_assay`
#' @noRd
abort_no_assay <- function(assay, se) {
  msg <- c(
    "Can't find assay {.arg {assay}} in assay names of {.arg {se}}.",
    "i" = "Did you check for valid assays with {.code SummarizedExperiment::assayNames(se)}?"
  )

  cli::cli_abort(
    message = msg,
    class = "error_no_assay"
  )
}


#' Throw an error when `SummarizedExperiment` does not contain a requested model
#' formula name
#'
#' @param name A character formula name that is not in the modeling formulas
#'   metadata of the `SummarizedExperiment`.
#' @param se The `SummarizedExperiment` that does not contain the requested
#'   formula `name`.
#'
#' @return An `rlang_error` of custom class `error_no_formula`.
#' @noRd
abort_no_formula <- function(name, se) {
  msg <- c(
    "Can't find formula with name {.arg {name}} in formulas slot of {.arg {se}} metadata.",
    "i" = "Did you check for valid formula names with {.code S4Vectors::metadata({se})[[\"metadata\"]][[\"formulas\"]]}."
  )

  cli::cli_abort(
    message = msg,
    class = "error_no_formula"
  )
}
