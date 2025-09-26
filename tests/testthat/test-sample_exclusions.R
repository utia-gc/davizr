# test-sample_exclusions.R
# Test behaviors related to excluding samples from an analysis

test_that("plot_sample_exclusions_heatmap plots a heatmap of excluded samples", {
  # construct inputs
  se <- random_se()
  # build list of samples to exclude
  sample_exclusions <- list(
    "Low depth" = c("sample1", "sample3"),
    "High rRNA" = c("sample7", "sample3")
  )

  # expect plotting is successful
  expect_no_error(
    sample_exclusions_heatmap <- plot_sample_exclusions_heatmap(se, sample_exclusions)
  )
  # expect returns plotly object
  expect_s3_class(sample_exclusions_heatmap, "plotly")
  expect_s3_class(sample_exclusions_heatmap, "htmlwidget")
})


test_that("drop_excluded_samples() drops excluded samples from a SummarizedExperiment", {
   # construct inputs
  se <- random_se()
  # build list of samples to exclude
  sample_exclusions <- list(
    "Low depth" = c("sample1", "sample3"),
    "High rRNA" = c("sample7", "sample3")
  )

  # expect sample dropping is successful
  expect_no_error(
    se <- drop_excluded_samples(se, sample_exclusions)
  )
  # expect filtered SummarizedExperiment is returned
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(ncol(se), 5L)
  expect_setequal(
    colnames(se),
    c("sample2", "sample4", "sample5", "sample6", "sample8")
  )
})


test_that("drop_excluded_samples() throws error when attempt to exclude sample that isn't in SummarizedExperiment", {
   # construct inputs
  se <- random_se()
  # build list of samples to exclude
  sample_exclusions <- list(
    # "sampleA" is not a column in se
    "Low depth" = c("sampleA", "sample3"),
    "High rRNA" = c("sample7", "sample3")
  )

  # expect correct error thrown
  expect_error(
    drop_excluded_samples(se, sample_exclusions),
    class = "error_no_sample"
  )
})
