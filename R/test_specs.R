# test_specs.R
# Specify tests for differential analysis


#' Specify tests for differential analysis
#'
#' @param path A character path to YAML of test specs
#' @inheritParams new_TestSpecs
#'
#' @returns A `TestSpecs` object
#' @export
TestSpecs <- function(path, alpha, lfc_threshold) {
  # check that file at path exists
  if (!file.exists(path)) abort_file_does_not_exist(path)

  # read tests data from YAML
  tests_data <- read_tests_yaml(path)
  # read contrasts data from YAML
  contrasts_data <- read_contrasts_yaml(path)

  # construct TestSpecs object
  test_specs <- new_TestSpecs(
    tests = tests_data,
    contrasts = contrasts_data,
    alpha = alpha,
    lfc_threshold = lfc_threshold
  )

  # validate
  validate_TestSpecs(test_specs)

  return(test_specs)
}


#' Low level `TestSpecs` constructor
#'
#' @param tests A list of tests data
#' @param contrasts A list of contrasts data
#' @param alpha A numeric alpha threshold for differential features
#' @param lfc_threshold A numeric log2 fold change threshold for differential features
#'
#' @returns A `TestSpecs` object
new_TestSpecs <- function(tests, contrasts, alpha, lfc_threshold) {
  test_specs <- structure(
    list(
      tests = tests,
      contrasts = contrasts,
      thresholds = list(
        alpha = alpha,
        lfc = lfc_threshold
      )
    ),
    class = c("TestSpecs", "list")
  )

  return(test_specs)
}


#' Validate a `TestSpecs` object
#'
#' @param test_specs A `TestSpecs` object
#'
#' @returns (Invisibly) a `TestSpecs` object
validate_TestSpecs <- function(test_specs) {
  # validate that any contrast names listed in the `tests` field is present in the `contrasts` field
  # build vector of contrast names
  contrast_names <- names(get_contrasts(test_specs))
 # check that all contrast names in tests are also in contrast names
  contrast_names_not_in_contrasts <- c()
  for (test in test_specs[["tests"]]) {
    contrast_names_not_in_contrasts <- append(contrast_names_not_in_contrasts, setdiff(test[["contrast_names"]], contrast_names))
  }

  if (length(contrast_names_not_in_contrasts) > 0) {
    msg <- c(
      "{.cls TestSpecs} has contrast names in {.field tests} that are not names in {.field contrasts}.",
      "i" = "All contrast names listed in {.field tests} must be named contrasts in {.field contrasts}.",
      "i" = "Problematic contrast names: {.str {contrast_names_not_in_contrasts}}."
    )
    cli::cli_abort(
      message = msg,
      class = "error_contrast_name_not_in_contrasts"
    )
  }

  return(invisible(test_specs))
}


#' Get all contrasts from a `TestSpecs` object
#'
#' @param test_specs A `TestSpecs` object
#'
#' @returns A named character vector. Names are contrast names, values are contrasts.
#' @export
get_contrasts <- function(test_specs) {
  # build vector of contrast names
  contrasts <- c()
  for (contrast in test_specs[["contrasts"]]) {
    named_contrast <- stats::setNames(contrast[["contrast"]], contrast[["name"]])
    contrasts <- append(contrasts, named_contrast)
  }

  return(contrasts)
}


#' Return alpha value of a `TestSpecs` object
#'
#' @inheritParams get_contrasts
#'
#' @returns A numeric alpha threshold
#' @export
get_alpha <- function(test_specs) {
  alpha <- test_specs[["thresholds"]][["alpha"]]

  return(alpha)
}


#' Return log2 fold change threshold of a `TestSpecs` object
#'
#' @inheritParams get_contrasts
#'
#' @returns A numeric log2 fold change threshold
#' @export
get_lfc_threshold <- function(test_specs) {
  lfc_threshold <- test_specs[["thresholds"]][["lfc"]]

  return(lfc_threshold)
}


read_tests_yaml <- function(path) {
  # read tests data from YAML file
  tests_yaml <- yaml::read_yaml(path)[["tests"]]

  # check contrasts key exists in test specs YAML
  if (is.null(tests_yaml)) abort_malformed_test_specs_yaml(path)

  #  preallocate the list size based on the number of tests in the YAML file
  tests_data <- vector("list", length = length(tests_yaml))
  # populate tests data with a list of test data for each test
  tests_names <- names(tests_yaml)
  for (i in 1:length(tests_yaml)) {
    tests_data[[i]] <- list(
      name = tests_names[[i]],
      contrast_names = tests_yaml[[i]][["contrasts"]],
      description = tests_yaml[[i]][["description"]]
    )
  }

  return(tests_data)
}


read_contrasts_yaml <- function(path) {
  # read contrasts data from YAML
  contrasts_yaml <- yaml::read_yaml(path)[["contrasts"]]

  # check contrasts key exists in test specs YAML
  if (is.null(contrasts_yaml)) abort_malformed_test_specs_yaml(path)

  # preallocate the list size based on the number of contrasts in the YAML file
  contrasts_data <- vector("list", length = length(contrasts_yaml))
  # populate contrasts data with a list of contrast data for each contrast
  contrasts_names <- names(contrasts_yaml)
  for (i in 1:length(contrasts_yaml)) {
    contrasts_data[[i]] <- list(
      name = contrasts_names[[i]],
      contrast = contrasts_yaml[[i]][["contrast"]],
      description = contrasts_yaml[[i]][["description"]]
    )
  }

  return(contrasts_data)
}


abort_file_does_not_exist <- function(path) {
  cli::cli_abort(
    message = c("File {.path {path}} does not exist"),
    class = "error_file_does_not_exist"
  )
}


abort_malformed_test_specs_yaml <- function(path) {
  msg <- c(
    "Test specs YAML file {.path {path}} is malformed.",
    "i" = "Test specs YAML file must contain top level keys: {.field tests} and {.field contrasts}."
  )

  cli::cli_abort(
    message = msg,
    class = "error_malformed_test_specs_yaml"
  )
}
