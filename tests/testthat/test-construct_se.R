test_that("construct_se() makes the correct SummarizedExperiment when counts colnames and samples rownames agree", {
  # test that the SummarizedExperiment object is created as expected
  expect_no_error(
    se <- construct_se(counts = random_matrix(), samples = random_col_data())
  )
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(se, random_se())
})


test_that("construct_se() makes the correct SummarizedExperiment when counts colnames and samples rownames in different order", {
  # construct counts and samples data with column and row order shuffled, respectively
  counts <- random_matrix()[, withr::with_seed(0, sample(ncol(random_matrix())))]
  col_data <- random_col_data()[withr::with_seed(1, sample(nrow(random_col_data()))), ]

  # test that the SummarizedExperiment object is created as expected
  expect_no_error(
    se <- construct_se(counts = counts, samples = col_data)
  )
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(se, random_se())
})


test_that("construct_se() makes the correct SummarizedExperiment when counts has a sample not in samples", {
  # construct counts and samples data where counts has a sample not found in samples data
  counts <- cbind(random_matrix(), sample9 = c(0))
  col_data <- random_col_data()

  # test that the SummarizedExperiment object is created as expected
  expect_no_error(
    se <- construct_se(counts = counts, samples = col_data)
  )
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(se, random_se())
})


test_that("construct_se() makes the correct SummarizedExperiment when samples has a sample not in counts", {
  # construct counts and samples data where counts has a sample not found in samples data
  counts <- random_matrix()
  col_data <- rbind(random_col_data(), sample9 = c("treat", 1))

  # test that the SummarizedExperiment object is created as expected
  expect_no_error(
    se <- construct_se(counts = counts, samples = col_data)
  )
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(se, random_se())
})


test_that("construct_se() errors when counts and colData have no samples in common", {
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
