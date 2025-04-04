test_that("perform_variance_partition_analysis() throws 'error_no_assay' when specified assay not present in SummarizedExperiment", {
  # create objects for running tests
  se <- random_se() |>
    set_formula(name = "vp_explanatory", formula = ~ (1 | condition) + (1 | batch))

  # expect throws error
  expect_error(
    perform_variance_partition_analysis(se, assay = "logcpm", formula = "vp_explanatory"),
    class = "error_no_assay"
  )
})


test_that("perform_variance_partition_analysis() throws 'error_no_formula' when specified formula not present in SummarizedExperiment", {
  # create objects for running tests
  se <- random_se()
  SummarizedExperiment::assay(se, "vsd") <- DESeq2::varianceStabilizingTransformation(SummarizedExperiment::assay(se, "counts"))

  # expect throws error
  expect_error(
    perform_variance_partition_analysis(se, assay = "vsd", formula = "vp_explanatory"),
    class = "error_no_formula"
  )
})


test_that("perform_variance_partition_analysis() performs a variance partitioning analysis", {
  # create objects for running tests
  se <- random_se()
  SummarizedExperiment::assay(se, "vsd") <- DESeq2::varianceStabilizingTransformation(SummarizedExperiment::assay(se, "counts"))
  se <- set_formula(se, name = "vp_explanatory", formula = ~ (1 | condition) + (1 | batch))

  # expect run successfully
  expect_no_error(
    vp <- perform_variance_partition_analysis(se, assay = "vsd", formula = "vp_explanatory")
  )
  # expect returns varPartResults object
  expect_s4_class(vp, "varPartResults")
  # expect condition and batch variables are in results
  expect_contains(colnames(vp), c("condition", "batch", "Residuals"))
  # expect that bar plots and variance partitioning plots can be made from the results
  expect_no_error(
    vp_bar_plot <- variancePartition::plotPercentBars(vp[1:5, ])
  )
  vdiffr::expect_doppelganger(
    "vp bar plot",
    vp_bar_plot
  )
  expect_no_error(
    vp_plot <- variancePartition::plotVarPart(vp)
  )
  vdiffr::expect_doppelganger(
    "vp violin plot",
    vp_plot
  )
})
