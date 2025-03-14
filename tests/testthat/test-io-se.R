test_that("write_se() writes a SummarizedExperiment object to a file", {
  # create SummarizedExperiment object
  counts <- random_matrix()
  col_data <- random_col_data()
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = col_data
  )

  # set the file path to write the SE directly into a temp directory
  dir <- withr::local_tempdir()
  se_path <- file.path(dir, "se.rds")

  # test that write_se returns TRUE for success and created a file at the correct path
  expect_false(file.exists(se_path))
  expect_no_warning(
    written <- write_se(se, path = se_path)
  )
  expect_true(written)
  expect_true(file.exists(se_path))
  expect_s4_class(readRDS(se_path), "SummarizedExperiment")
})


test_that("write_se() throws a warning when overwrite set to false and file already exists at path", {
  # create SummarizedExperiment object
  counts <- random_matrix()
  col_data <- random_col_data()
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = col_data
  )

  # set the file path to write the SE directly into a temp directory
  dir <- withr::local_tempdir()
  se_path <- file.path(dir, "se.rds")

  # create a file at SE path
  file.create(se_path)

  # test that write_se throws a warning when file exists at se_path and overwrite is false
  expect_true(file.exists(se_path))
  expect_warning(written <- write_se(se, path = se_path), class = "warning_no_overwrite")
  expect_false(written)
  expect_true(file.exists(se_path))
  expect_error(readRDS(se_path), "error reading from connection")
})


test_that("write_se() overwrites existing file when overwrite set to true and file already exists at path", {
  # create SummarizedExperiment object
  counts <- random_matrix()
  col_data <- random_col_data()
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = col_data
  )

  # set the file path to write the SE directly into a temp directory
  dir <- withr::local_tempdir()
  se_path <- file.path(dir, "se.rds")

  # create a file at SE path
  file.create(se_path)

  # test that write_se throws a warning when file exists at se_path and overwrite is false
  expect_true(file.exists(se_path))
  expect_no_warning(written <- write_se(se, path = se_path, overwrite = TRUE))
  expect_true(written)
  expect_true(file.exists(se_path))
  expect_s4_class(readRDS(se_path), "SummarizedExperiment")
})
