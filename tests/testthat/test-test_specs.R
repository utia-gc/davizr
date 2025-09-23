# test-test_specs.R
# Test functionality related to test specifications.

test_that("TestSpecs() with a valid test specs YAML file returns a TestSpecs object with working accessors", {
  # create test specs object from path to YAML
  test_specs_yaml <- test_path("data", "test-specs", "test-specs-valid.yml")
  test_specs <- TestSpecs(
    test_specs_yaml,
    alpha = 0.05,
    lfc_threshold = 0
  )

  # test that TestSpecs class is returned
  expect_s3_class(test_specs, "TestSpecs")

  # test that TestSpecs has expected names
  expect_setequal(
    names(test_specs),
    c("tests", "contrasts", "thresholds")
  )

  # test that tests field is a list of test records formatted as expected
  expected_tests <- list(
    list(
      name = "Treated vs Untreated",
      contrasts = c("treated_vs_untreated"),
      description = "Difference between treated and untreated."
    )
  )
  expect_equal(test_specs[["tests"]], expected_tests)

  # test that contrasts field is a list of contrast records formatted as expected
  expected_contrasts <- list(
    list(
      name = "treated_vs_untreated",
      contrast = "conditiontreated - conditionuntreated",
      description = "Difference between treated and untreated."
    )
  )
  expect_equal(test_specs[["contrasts"]], expected_contrasts)

  # test accessors
  expect_equal(
    get_alpha(test_specs),
    0.05
  )
  expect_equal(
    get_lfc_threshold(test_specs),
    0
  )
})


test_that("TestSpecs() throws error when test specs YAML file does not exist", {
  # create test specs object from path to YAML
  test_specs_yaml <- test_path("data", "test-specs", "test-specs-nonexistent_file.yml")
  expect_error(
    test_specs <- TestSpecs(
      test_specs_yaml,
      alpha = 0.05,
      lfc_threshold = 0
    ),
    class = "error_file_does_not_exist"
  )
})


test_that("TestSpecs() throws error when 'tests' key does not exist in test specs YAML file", {
  # create test specs object from path to YAML
  test_specs_yaml <- test_path("data", "test-specs", "test-specs-invalid-no_tests.yml")
  expect_error(
    test_specs <- TestSpecs(
      test_specs_yaml,
      alpha = 0.05,
      lfc_threshold = 0
    ),
    class = "error_malformed_test_specs_yaml"
  )
})


test_that("TestSpecs() throws error when 'contrasts' key does not exist in test specs YAML file", {
  # create test specs object from path to YAML
  test_specs_yaml <- test_path("data", "test-specs", "test-specs-invalid-no_contrasts.yml")
  expect_error(
    test_specs <- TestSpecs(
      test_specs_yaml,
      alpha = 0.05,
      lfc_threshold = 0
    ),
    class = "error_malformed_test_specs_yaml"
  )
})
