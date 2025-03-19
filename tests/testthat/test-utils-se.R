test_that("get_samples_data() returns a data.frame of samples data with sample names as a column called 'sample' by default", {
  # create objects for running tests
  se <- random_se()
  sample_names <- paste0("sample", seq_along(1:8))

  # expect function works without signalling anything
  expect_no_condition(
    samples_data <- get_samples_data(se)
  )
  # returns data.frame object
  expect_s3_class(samples_data, "data.frame")
  # 'sample' variable in data.frame
  expect_in("sample", colnames(samples_data))
  expect_setequal(sample_names, samples_data[["sample"]])
})


test_that("get_samples_data() returns a data.frame of samples data with sample names as a column called 'custom_sample_variable' when `rownames_var = 'custom_sample_variable'", {
  # create objects for running tests
  se <- random_se()
  sample_names <- paste0("sample", seq_along(1:8))

  # expect function works without signalling anything
  expect_no_condition(
    samples_data <- get_samples_data(se, rownames_var = "custom_sample_variable")
  )
  # returns data.frame object
  expect_s3_class(samples_data, "data.frame")
  # 'sample' variable in data.frame
  expect_in("custom_sample_variable", colnames(samples_data))
  expect_setequal(sample_names, samples_data[["custom_sample_variable"]])
})


test_that("get_samples_data() returns a data.frame of samples data with sample names as rownames when `rownames_var = NULL`", {
  # create objects for running tests
  se <- random_se()
  sample_names <- paste0("sample", seq_along(1:8))

  # expect function works without signalling anything
  expect_no_condition(
    samples_data <- get_samples_data(se, rownames_var = NULL)
  )
  # returns data.frame object
  expect_s3_class(samples_data, "data.frame")
  # 'sample' variable not in data.frame
  expect_setequal(c("condition", "batch", "library_size"), colnames(samples_data))
  expect_setequal(sample_names, rownames(samples_data))
})
