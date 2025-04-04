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


#' Set aesthetics metadata in a `SummarizedExperiment`
#'
#' @inheritParams write_se
#' @param variable A character vector for the sample data variable to set
#'   aesthetic values for.
#' @param values A named vector mapping values in the sample data `variable` to
#'   values of an aesthetic.
#' @param aesthetic A vector indicating the aesthetic to set.
#'
#' @return The `SummarizedExperiment` with the aesthetics set in the metadata
#'   slot.
#' @export
set_aesthetics <- function(se, variable, values, aesthetic = c("color", "shape")) {

  # throw error if values isn't a named vector
  if (is.null(names(values))) {
    cli::cli_abort(
      message = "{.arg values} must be a named vector.",
      class = "error_unnamed_values"
    )
  }

  variable_data <- se[[variable]]

  # throw error if variable not in sample data
  if (is.null(variable_data)) {
    abort_no_sample_variable(var = variable, se = "se")
  }
  # throw error if value has a name that's not in the variable data
  if (!all(names(values) %in% variable_data)) {
    cli::cli_abort(
      message = "Every name in {.arg values} must be present in se[[variable]].",
      class = "error_invalid_variable_data"
    )
  }

  S4Vectors::metadata(se)[["aesthetics"]][[aesthetic]][[variable]] <- values

  return(se)
}


#' Set a formula in the metadata of a `SummarizedExperiment`
#'
#' @inheritParams write_se
#' @param name A character name for the formula.
#' @param formula A formula
#'
#' @return The `SummarizedExperiment` with the formula added to the metadata
#'   slot.
#' @export
set_formula <- function(se, name, formula) {

  # throw error if a variable in the formula isn't a column in sample data
  formula_variables <- all.vars(formula)
  sample_data_variables <- colnames(get_samples_data(se))
  if (!all(formula_variables %in% sample_data_variables)) {
    invalid_sample_variables <- setdiff(formula_variables, sample_data_variables)
    abort_no_sample_variable(var = invalid_sample_variables, se = "se")
  }

  # add the formula to the metadata slot for modeling formulas
  S4Vectors::metadata(se)[["modeling"]][["formulas"]][[name]] <- formula

  return(se)
}


#' Get a formula from the metadata of a `SummarizedExperiment` by its name
#'
#' @inheritParams set_formula
#'
#' @return The requested formula
#' @export
get_formula <- function(se, name) {

  formula <- S4Vectors::metadata(se)[["modeling"]][["formulas"]][[name]]

  return(formula)
}

