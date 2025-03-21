test_that("construct_se() makes the correct SummarizedExperiment when counts colnames and samples rownames agree", {
  # create expected objects to test against
  se_expected <- random_se()

  # test that the SummarizedExperiment object is created as expected
  expect_no_error(
    se <- construct_se(counts = random_matrix(), samples = random_col_data())
  )
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(se, se_expected)
  # library size variable is correct
  expect_in("library_size", colnames(SummarizedExperiment::colData(se)))
  expect_true(is.numeric(se[["library_size"]]))
})


test_that("construct_se() makes the correct SummarizedExperiment when counts colnames and samples rownames in different order", {
  # create expected objects to test against
  se_expected <- random_se()

  # construct counts and samples data with column and row order shuffled, respectively
  counts <- random_matrix()[, withr::with_seed(0, sample(ncol(random_matrix())))]
  col_data <- random_col_data()[withr::with_seed(1, sample(nrow(random_col_data()))), ]

  # test that the SummarizedExperiment object is created as expected
  expect_no_error(
    se <- construct_se(counts = counts, samples = col_data)
  )
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(se, se_expected)
  # library size variable is correct
  expect_in("library_size", colnames(SummarizedExperiment::colData(se)))
  expect_true(is.numeric(se[["library_size"]]))
})


test_that("construct_se() makes the correct SummarizedExperiment when counts has a sample not in samples", {
  # create expected objects to test against
  se_expected <- random_se()

  # construct counts and samples data where counts has a sample not found in samples data
  counts <- cbind(random_matrix(), sample9 = c(0))
  col_data <- random_col_data()

  # test that the SummarizedExperiment object is created as expected
  expect_no_error(
    se <- construct_se(counts = counts, samples = col_data)
  )
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(se, se_expected)
  # library size variable is correct
  expect_in("library_size", colnames(SummarizedExperiment::colData(se)))
  expect_true(is.numeric(se[["library_size"]]))
})


test_that("construct_se() makes the correct SummarizedExperiment when samples has a sample not in counts", {
  # create expected objects to test against
  se_expected <- random_se()

  # construct counts and samples data where counts has a sample not found in samples data
  counts <- random_matrix()
  col_data <- rbind(random_col_data(), sample9 = c("treat", 1))

  # test that the SummarizedExperiment object is created as expected
  expect_no_error(
    se <- construct_se(counts = counts, samples = col_data)
  )
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(se, se_expected)
  # library size variable is correct
  expect_in("library_size", colnames(SummarizedExperiment::colData(se)))
  expect_true(is.numeric(se[["library_size"]]))
})


test_that("construct_se() makes the correct SummarizedExperiment with correct library size variable", {
  # create expected objects to test against
  se_expected <- random_se()
  se_expected[["library_size_raw"]] <- se_expected[["library_size"]]
  se_expected[["library_size"]] <- NULL

  # test that the SummarizedExperiment object is created as expected
  expect_no_error(
    se <- construct_se(counts = random_matrix(), samples = random_col_data(), library_size_var = "library_size_raw")
  )
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(se, se_expected)
  # custom library size variable is correct
  expect_in("library_size_raw", colnames(SummarizedExperiment::colData(se)))
  expect_true(is.numeric(se[["library_size_raw"]]))
})


test_that("construct_se() makes the correct SummarizedExperiment without library size variable when `library_size_var = NULL`", {
  # create expected objects to test against
  se_expected <- random_se()
  se_expected[["library_size"]] <- NULL

  # test that the SummarizedExperiment object is created as expected
  expect_no_error(
    se <- construct_se(counts = random_matrix(), samples = random_col_data(), library_size_var = NULL)
  )
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(se, se_expected)
  # custom library size variable is correct
  expect_false("library_size" %in% colnames(SummarizedExperiment::colData(se)))
})


test_that("construct_se() makes the correct SummarizedExperiment with modeling variables set", {
  # create expected objects to test against
  se_expected <- random_se()
  S4Vectors::metadata(se_expected)[["modeling"]][["variables"]][["interest"]] <- c("condition")
  S4Vectors::metadata(se_expected)[["modeling"]][["variables"]][["nuisance"]] <- c("batch")
  S4Vectors::metadata(se_expected)[["modeling"]][["variables"]][["screening"]] <- c("library_size")

  # test that the SummarizedExperiment object is created as expected
  expect_no_error(
    se <- construct_se(
      counts = random_matrix(),
      samples = random_col_data(),
      interest_variables = c("condition"),
      nuisance_variables = c("batch"),
      screening_variables = c("library_size")
    )
  )
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(se, se_expected)
  # library size variable is correct
  expect_in("library_size", colnames(SummarizedExperiment::colData(se)))
  expect_true(is.numeric(se[["library_size"]]))
  # modeling variables are correct
  expect_identical(
    S4Vectors::metadata(se)[["modeling"]][["variables"]][["interest"]],
    c("condition")
  )
  expect_identical(
    S4Vectors::metadata(se)[["modeling"]][["variables"]][["nuisance"]],
    c("batch")
  )
  expect_identical(
    S4Vectors::metadata(se)[["modeling"]][["variables"]][["screening"]],
    c("library_size")
  )
})


test_that("construct_se() throws error when attempting to set interest variables to variable that doesn't exist in sample data", {
  expect_error(
    construct_se(
      counts = random_matrix(),
      samples = random_col_data(),
      interest_variables = c("disease")
    ),
    class = "error_no_sample_variable"
  )
})


test_that("construct_se() throws error when attempting to set nuisance variables to variable that doesn't exist in sample data", {
  expect_error(
    construct_se(
      counts = random_matrix(),
      samples = random_col_data(),
      nuisance_variables = c("environment")
    ),
    class = "error_no_sample_variable"
  )
})


test_that("construct_se() throws error when attempting to set screening variables to variable that doesn't exist in sample data", {
  expect_error(
    construct_se(
      counts = random_matrix(),
      samples = random_col_data(),
      screening_variables = c("weight", "height")
    ),
    class = "error_no_sample_variable"
  )
})


test_that("construct_se() throws error when counts and sample data have no samples in common", {
  counts <- random_matrix()
  colnames(counts) <- paste0("samp", 1:8)
  col_data <- random_col_data()
  rownames(col_data) <- c(paste0("ctrl", 1:4), paste0("treat", 1:4))

  expect_error(
    construct_se(
      counts = counts, samples = col_data
    ),
    class = "error_no_common_samples"
  )
})
