# test_specs.R
# Specify tests for differential analysis


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

  return(test_specs)
}


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


get_alpha <- function(test_specs) {
  alpha <- test_specs[["thresholds"]][["alpha"]]

  return(alpha)
}


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
