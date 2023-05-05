test_that("construct_matrix_dds() makes the correct DESeqDataSet when counts colnames and colData rownames agree", {
  counts <- random_matrix()
  col_data <- random_col_data()
  design <- "~ batch + condition"

  expect_equal(
    construct_matrix_dds(
      countData = counts, colData = col_data, design = design
    ),
    random_dds(),
    ignore_function_env = TRUE,
    ignore_formula_env = TRUE
  )
})

test_that("construct_matrix_dds() makes the correct DESeqDataSet when counts colnames and colData rownames in different order", {
  counts <- random_matrix()[, withr::with_seed(0, sample(ncol(random_matrix())))]
  col_data <- random_col_data()[withr::with_seed(1, sample(nrow(random_col_data()))), ]
  design <- "~ batch + condition"

  expect_equal(
    construct_matrix_dds(
      countData = counts, colData = col_data, design = design
    ),
    random_dds(),
    ignore_function_env = TRUE,
    ignore_formula_env = TRUE
  )
})

test_that("construct_matrix_dds() makes the correct DESeqDataSet when counts has a sample not in colData", {
  counts <- cbind(random_matrix(), sample9 = c(0))
  col_data <- random_col_data()
  design <- "~ batch + condition"

  expect_equal(
    construct_matrix_dds(
      countData = counts, colData = col_data, design = design
    ),
    random_dds(),
    ignore_function_env = TRUE,
    ignore_formula_env = TRUE
  )
})

test_that("construct_matrix_dds() makes the correct DESeqDataSet when colData has a sample not in counts", {
  counts <- random_matrix()
  col_data <- rbind(random_col_data(), sample9 = c("treat", 1))
  design <- "~ batch + condition"

  expect_equal(
    construct_matrix_dds(
      countData = counts, colData = col_data, design = design
    ),
    random_dds(),
    ignore_function_env = TRUE,
    ignore_formula_env = TRUE
  )
})

test_that("construct_matrix_dds() errors when counts and colData have no samples in common", {
  counts <- random_matrix()
  colnames(counts) <- paste0("samp", 1:8)
  col_data <- random_col_data()
  rownames(col_data) <- c(paste0("ctrl", 1:4), paste0("treat", 1:4))
  design <- "~ batch + condition"

  expect_error(
    construct_matrix_dds(
      countData = counts, colData = col_data, design = design
    ),
    class = "no_common_samples"
  )
})
