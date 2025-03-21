test_that("plot_modeling_variables_pairs() returns a scatter plot matrix of sample data for all modeling variables", {
  # construct SummarizedExperiment to test
  se <- random_se()
  se <- set_modeling_variable(se, variable = "interest", values = c("condition"))
  se <- set_modeling_variable(se, variable = "nuisance", values = c("batch"))
  se <- set_modeling_variable(se, variable = "screening", values = c("library_size"))

  expect_no_error(
    modeling_variables_pairs_plot <- plot_modeling_variables_pairs(se)
  )
  # returns ggmatrix object
  expect_s3_class(modeling_variables_pairs_plot, "ggmatrix")
  # test plot is the same
  vdiffr::expect_doppelganger(
    "modeling variables pairs plot",
    modeling_variables_pairs_plot
  )
})


test_that("plot_modeling_variables_pairs() throws error when there are no modeling variables", {
  # construct SummarizedExperiment to test
  se <- random_se()

  expect_error(
    modeling_variables_pairs_plot <- plot_modeling_variables_pairs(se),
    class = "error_no_modeling_variables"
  )
})


test_that("plot_canonical_correlation() returns a canonical correlation plot for explanatory modeling variables", {
  # construct SummarizedExperiment to test
  se <- random_se()
  se <- set_modeling_variable(se, variable = "interest", values = c("condition"))
  se <- set_modeling_variable(se, variable = "nuisance", values = c("batch"))
  se <- set_modeling_variable(se, variable = "screening", values = c("library_size"))

  expect_no_error(
    explanatory_cc_plot <- plot_canonical_correlation(se, "explanatory")
  )
  # returns Heatmap object
  expect_s4_class(explanatory_cc_plot, "Heatmap")
  # plot is reproducible
  vdiffr::expect_doppelganger(
    "explanatory canonical correlation plot",
    explanatory_cc_plot
  )
})


test_that("plot_canonical_correlation() throws error when there are no explanatory modeling variables", {
  # construct SummarizedExperiment to test
  se <- random_se()

  expect_error(
    explanatory_cc_plot <- plot_canonical_correlation(se, "explanatory"),
    class = "error_no_modeling_variables"
  )
})


test_that("plot_canonical_correlation() returns a canonical correlation plot for all modeling variables", {
  # construct SummarizedExperiment to test
  se <- random_se()
  se <- set_modeling_variable(se, variable = "interest", values = c("condition"))
  se <- set_modeling_variable(se, variable = "nuisance", values = c("batch"))
  se <- set_modeling_variable(se, variable = "screening", values = c("library_size"))

  expect_no_error(
    all_cc_plot <- plot_canonical_correlation(se, "all")
  )
  # returns Heatmap object
  expect_s4_class(all_cc_plot, "Heatmap")
  # plot is reproducible
  vdiffr::expect_doppelganger(
    "all canonical correlation plot",
    all_cc_plot
  )
})


test_that("plot_canonical_correlation() throws error when there are no modeling variables", {
  # construct SummarizedExperiment to test
  se <- random_se()

  expect_error(
    all_cc_plot <- plot_canonical_correlation(se, "all"),
    class = "error_no_modeling_variables"
  )
})


test_that("tbl_sample_data() returns datatable of sample data from SummarizedExperiment by default", {
  expect_no_error(
    sample_data <- tbl_sample_data(random_se())
  )
  # returns datatables object
  expect_s3_class(sample_data, "datatables")
  expect_s3_class(sample_data, "htmlwidget")
  # underlying data has expected shape and names
  expect_equal(
    dim(sample_data$x$data),
    c(8L, 4L)
  )
  expect_in(
    c("sample", "condition", "batch", "library_size"),
    colnames(sample_data$x$data)
  )
})


test_that("tbl_sample_data() returns datatable of sample data from SummarizedExperiment when `interactive = TRUE`", {
  expect_no_error(
    sample_data <- tbl_sample_data(random_se(), interactive = TRUE)
  )
  # returns datatables object
  expect_s3_class(sample_data, "datatables")
  expect_s3_class(sample_data, "htmlwidget")
  # underlying data has expected shape and names
  expect_equal(
    dim(sample_data$x$data),
    c(8L, 4L)
  )
  expect_in(
    c("sample", "condition", "batch", "library_size"),
    colnames(sample_data$x$data)
  )
})


test_that("tbl_sample_data() returns datatable of sample data from SummarizedExperiment when `interactive = FALSE`", {
  expect_no_error(
    sample_data <- tbl_sample_data(random_se(), interactive = FALSE)
  )
  # returns data.frame object
  expect_s3_class(sample_data, "data.frame")
  # underlying data has expected shape and names
  expect_equal(
    dim(sample_data),
    c(8L, 4L)
  )
  expect_in(
    c("sample", "condition", "batch", "library_size"),
    colnames(sample_data)
  )
})


test_that("tbl_interest_variables_freq() returns frequency table of variables of interest from SummarizedExperiment as datatable by default", {
  # construct expected frequency table to test against
  expected_freq_tbl <- data.frame(
    condition = factor(c("control", "treat")),
    Freq = c(4L, 4L)
  )

  # construct SummarizedExperiment to test
  se <- random_se()
  se <- set_modeling_variable(se, variable = "interest", values = c("condition"))

  expect_no_error(
    actual_freq_tbl <- tbl_interest_variables_freq(se)
  )
  # returns datatable object
  expect_s3_class(actual_freq_tbl, "datatables")
  expect_s3_class(actual_freq_tbl, "htmlwidget")
  # underlying data has expected shape and names
  expect_equal(
    dim(actual_freq_tbl$x$data),
    c(2L, 2L)
  )
  expect_setequal(
    c("condition", "Freq"),
    colnames(actual_freq_tbl$x$data)
  )
  # underlying data matches expected
  expect_identical(
    actual_freq_tbl$x$data,
    expected_freq_tbl
  )
})


test_that("tbl_interest_variables_freq() returns frequency table of multiple variables of interest from SummarizedExperiment as datatable by default", {
  # construct expected frequency table to test against
  expected_freq_tbl <- data.frame(
    condition = factor(c("control", "treat", "control", "treat")),
    batch = factor(c(1, 1, 2, 2)),
    Freq = c(2L, 2L, 2L, 2L)
  )

  # construct SummarizedExperiment to test
  se <- random_se()
  se <- set_modeling_variable(se, variable = "interest", values = c("condition", "batch"))

  expect_no_error(
    actual_freq_tbl <- tbl_interest_variables_freq(se)
  )
  # returns datatable object
  expect_s3_class(actual_freq_tbl, "datatables")
  expect_s3_class(actual_freq_tbl, "htmlwidget")
  # underlying data has expected shape and names
  expect_equal(
    dim(actual_freq_tbl$x$data),
    c(4L, 3L)
  )
  expect_setequal(
    c("condition", "batch", "Freq"),
    colnames(actual_freq_tbl$x$data)
  )
  # underlying data matches expected
  expect_identical(
    actual_freq_tbl$x$data,
    expected_freq_tbl
  )
})


test_that("tbl_interest_variables_freq() returns frequency table of variables of interest from SummarizedExperiment as datatable when `table_format = 'datatable'", {
  # construct expected frequency table to test against
  expected_freq_tbl <- data.frame(
    condition = factor(c("control", "treat")),
    Freq = c(4L, 4L)
  )

  # construct SummarizedExperiment to test
  se <- random_se()
  se <- set_modeling_variable(se, variable = "interest", values = c("condition"))

  expect_no_error(
    actual_freq_tbl <- tbl_interest_variables_freq(se, table_format = "datatable")
  )
  # returns datatable object
  expect_s3_class(actual_freq_tbl, "datatables")
  expect_s3_class(actual_freq_tbl, "htmlwidget")
  # underlying data has expected shape and names
  expect_equal(
    dim(actual_freq_tbl$x$data),
    c(2L, 2L)
  )
  expect_setequal(
    c("condition", "Freq"),
    colnames(actual_freq_tbl$x$data)
  )
  # underlying data matches expected
  expect_identical(
    actual_freq_tbl$x$data,
    expected_freq_tbl
  )
})


test_that("tbl_interest_variables_freq() returns frequency table of variables of interest from SummarizedExperiment as data.frame when `table_format = 'data.frame'", {
  # construct expected frequency table to test against
  expected_freq_tbl <- data.frame(
    condition = factor(c("control", "treat")),
    Freq = c(4L, 4L)
  )

  # construct SummarizedExperiment to test
  se <- random_se()
  se <- set_modeling_variable(se, variable = "interest", values = c("condition"))

  expect_no_error(
    actual_freq_tbl <- tbl_interest_variables_freq(se, table_format = "data.frame")
  )
  # returns datatable object
  expect_s3_class(actual_freq_tbl, "data.frame")
  # underlying data has expected shape and names
  expect_equal(
    dim(actual_freq_tbl),
    c(2L, 2L)
  )
  expect_setequal(
    c("condition", "Freq"),
    colnames(actual_freq_tbl)
  )
  # underlying data matches expected
  expect_identical(
    actual_freq_tbl,
    expected_freq_tbl
  )
})
