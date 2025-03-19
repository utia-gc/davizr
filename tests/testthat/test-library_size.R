test_that("tbl_library_size() returns datatable of read counts by default", {
  expect_no_error(
    library_sizes <- tbl_library_sizes(random_matrix())
  )
  # returns datatables object
  expect_s3_class(library_sizes, "datatables")
  expect_s3_class(library_sizes, "htmlwidget")
  # underlying data has expected shape and names
  expect_equal(
    dim(library_sizes$x$data),
    c(8L, 2L)
  )
  expect_in(
    c("Sample", "Library size"),
    colnames(library_sizes$x$data)
  )
})


test_that("tbl_library_size() returns datatable of read counts when interactive argument is TRUE", {
  expect_no_error(
    library_sizes <- tbl_library_sizes(random_matrix(), interactive = TRUE)
  )
  # returns datatables object
  expect_s3_class(library_sizes, "datatables")
  expect_s3_class(library_sizes, "htmlwidget")
  # underlying data has expected shape and names
  expect_equal(
    dim(library_sizes$x$data),
    c(8L, 2L)
  )
  expect_in(
    c("Sample", "Library size"),
    colnames(library_sizes$x$data)
  )
})


test_that("tbl_library_size() returns data.frame of read counts when interactive argument is FALSE", {
  expect_no_error(
    library_sizes <- tbl_library_sizes(random_matrix(), interactive = FALSE)
  )
  # returns data.frame object
  expect_s3_class(library_sizes, "data.frame")
  # underlying data has expected shape and names
  expect_equal(
    dim(library_sizes),
    c(8L, 2L)
  )
  expect_in(
    c("Sample", "Library size"),
    colnames(library_sizes)
  )
})


test_that("plot_logcpm_density() returns interactive plotly plot of read counts density by default", {
  logcpm_density_plot <- plot_logcpm_density(random_se())

  # returns plotly object
  expect_s3_class(logcpm_density_plot, "plotly")
  expect_s3_class(logcpm_density_plot, "htmlwidget")
})


test_that("plot_logcpm_density() returns interactive plotly plot of read counts density when interactive argument is TRUE", {
  logcpm_density_plot <- plot_logcpm_density(random_se(), interactive = TRUE)

  # returns plotly object
  expect_s3_class(logcpm_density_plot, "plotly")
  expect_s3_class(logcpm_density_plot, "htmlwidget")
})


test_that("plot_logcpm_density() returns static ggplot plot of read counts density when interactive argument is FALSE", {
  logcpm_density_plot <- plot_logcpm_density(random_se(), interactive = FALSE)

  # returns ggplot object
  expect_s3_class(logcpm_density_plot, "gg")
  expect_s3_class(logcpm_density_plot, "ggplot")
  # test plot is the same
  vdiffr::expect_doppelganger(
    "logcpm density plot noninteractive",
    logcpm_density_plot
  )
})


test_that("plot_library_sizes() throws an error when the library size variable isn't in SummarizedExperiment", {
  # create objects for testing
  se <- random_se()

  expect_error(
    plot_library_sizes(se, library_size_var = "raw_library_sizes"),
    class = "error_no_sample_variable"
  )
})


test_that("plot_library_sizes() returns interactive plotly plot of library sizes by default", {
  # create objects for testing
  se <- random_se()

  expect_no_error(
    library_sizes_plot <- plot_library_sizes(se, library_size_var = "library_size")
  )
  # returns plotly object
  expect_s3_class(library_sizes_plot, "plotly")
  expect_s3_class(library_sizes_plot, "htmlwidget")
})


test_that("plot_library_sizes() returns interactive plotly plot of library sizes by default with custom library size variable", {
  # create objects for testing
  se <- random_se()
  se <- add_library_size(se, library_size_var = "custom_library_size")

  expect_no_error(
    library_sizes_plot <- plot_library_sizes(se, library_size_var = "custom_library_size")
  )
  # returns plotly object
  expect_s3_class(library_sizes_plot, "plotly")
  expect_s3_class(library_sizes_plot, "htmlwidget")
})


test_that("plot_library_sizes() returns interactive plotly plot of library sizes when `interactive = TRUE`", {
  # create objects for testing
  se <- random_se()

  expect_no_error(
    library_sizes_plot <- plot_library_sizes(se, library_size_var = "library_size", interactive = TRUE)
  )
  # returns plotly object
  expect_s3_class(library_sizes_plot, "plotly")
  expect_s3_class(library_sizes_plot, "htmlwidget")
})


test_that("plot_library_sizes() returns static ggplot plot of library sizes when `interactive = FALSE`", {
  # create objects for testing
  se <- random_se()

  expect_no_error(
    library_sizes_plot <- plot_library_sizes(se, library_size_var = "library_size", interactive = FALSE)
  )
  # returns ggplot object
  expect_s3_class(library_sizes_plot, "gg")
  expect_s3_class(library_sizes_plot, "ggplot")
  # test plot is the same
  vdiffr::expect_doppelganger(
    "library sizes plot noninteractive",
    library_sizes_plot
  )
})
