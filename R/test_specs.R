# test_specs.R
# Specify tests for differential analysis


TestSpecs <- function(path, alpha, lfc_threshold) {
  # read tests data from YAML
  tests_yaml <- yaml::read_yaml(path)[["tests"]]
  # preallocate the list size based on the number of tests in the YAML file
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

  # read contrasts data from YAML
  contrasts_yaml <- yaml::read_yaml(path)[["contrasts"]]
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


get_alpha <- function(x, ...) {
  UseMethod("get_alpha")
}


get_alpha.TestSpecs <- function(test_specs) {
  alpha <- test_specs[["thresholds"]][["alpha"]]

  return(alpha)
}


get_lfc_threshold <- function(x, ...) {
  UseMethod("get_lfc_threshold")
}


get_lfc_threshold.TestSpecs <- function(test_specs) {
  lfc_threshold <- test_specs[["thresholds"]][["lfc"]]

  return(lfc_threshold)
}
