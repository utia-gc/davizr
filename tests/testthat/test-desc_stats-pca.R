test_that("perform_pca() throws 'error_no_assay' when specified assay not present in SummarizedExperiment", {
  # create objects for running tests
  se <- random_se()

  # expect throws error
  expect_error(
    perform_pca(se, assay = "logcpm"),
    class = "error_no_assay"
  )
})


test_that("perform_pca() performs PCA", {
  # create objects for running tests
  se <- random_se()
  SummarizedExperiment::assay(se, "vld") <- suppressMessages(DESeq2::varianceStabilizingTransformation(SummarizedExperiment::assay(se, "counts")))

  # expect runs successfully
  expect_no_error(
    pca <- perform_pca(se, assay = "vld")
  )
  # expect returns a PCA object
  expect_s3_class(pca, "PCA")
  expect_named(
    pca,
    c("scores", "variance_explained", "prcomp", "metadata")
  )
  # expect PCA is reproducible
  expect_snapshot(pca)
})


test_that("plot_scores() returns a plotly plot of PC scores for PC1 and PC2 by default", {
  # create objects for running tests
  se <- random_se()
  SummarizedExperiment::assay(se, "vld") <- suppressMessages(DESeq2::varianceStabilizingTransformation(SummarizedExperiment::assay(se, "counts")))
  pca <- perform_pca(se, assay = "vld")

  # expect runs successfully
  expect_no_error(
    scores_plot <- plot_scores(pca)
  )
  # returns plotly object
  expect_s3_class(scores_plot, "plotly")
  expect_s3_class(scores_plot, "htmlwidget")
})


test_that("plot_scores() returns a plotly plot of PC scores for PC1 and PC2 by default", {
  # create objects for running tests
  se <- random_se()
  SummarizedExperiment::assay(se, "vld") <- suppressMessages(DESeq2::varianceStabilizingTransformation(SummarizedExperiment::assay(se, "counts")))
  pca <- perform_pca(se, assay = "vld")

  # expect runs successfully
  expect_no_error(
    scores_plot <- plot_scores(pca, x = "PC2", y = "PC3")
  )
  # returns plotly object
  expect_s3_class(scores_plot, "plotly")
  expect_s3_class(scores_plot, "htmlwidget")
})


test_that("plot_scores() returns a plotly plot of PC scores for PC1 and PC2 by default when `plot_type = 'plotly'`", {
  # create objects for running tests
  se <- random_se()
  SummarizedExperiment::assay(se, "vld") <- suppressMessages(DESeq2::varianceStabilizingTransformation(SummarizedExperiment::assay(se, "counts")))
  pca <- perform_pca(se, assay = "vld")

  # expect runs successfully
  expect_no_error(
    scores_plot <- plot_scores(pca, plot_type = "plotly")
  )
  # returns plotly object
  expect_s3_class(scores_plot, "plotly")
  expect_s3_class(scores_plot, "htmlwidget")
})


test_that("plot_scores() returns a ggplot plot of PC scores for PC1 and PC2  when `plot_type = 'ggplot'`", {
  # create objects for running tests
  se <- random_se()
  SummarizedExperiment::assay(se, "vld") <- suppressMessages(DESeq2::varianceStabilizingTransformation(SummarizedExperiment::assay(se, "counts")))
  pca <- perform_pca(se, assay = "vld")

  # expect runs successfully
  expect_no_error(
    scores_plot <- plot_scores(pca, plot_type = "ggplot")
  )
  # expect returns ggplot object
  expect_s3_class(scores_plot, "gg")
  expect_s3_class(scores_plot, "ggplot")
  # test plot is the same
  vdiffr::expect_doppelganger(
    "pca scores ggplot",
    scores_plot
  )
})


test_that("plot_scores() returns a ggplot plot of PC scores for PC1 and PC2 with no margins plot when `plot_type = 'ggplot', margin_plot = 'none'`", {
  # create objects for running tests
  se <- random_se()
  SummarizedExperiment::assay(se, "vld") <- suppressMessages(DESeq2::varianceStabilizingTransformation(SummarizedExperiment::assay(se, "counts")))
  pca <- perform_pca(se, assay = "vld")

  # expect runs successfully
  expect_no_error(
    scores_plot <- plot_scores(pca, margin_plot = "none", plot_type = "ggplot")
  )
  # expect returns ggplot object
  expect_s3_class(scores_plot, "gg")
  expect_s3_class(scores_plot, "ggplot")
  # test plot is the same
  vdiffr::expect_doppelganger(
    "pca scores no margin ggplot",
    scores_plot
  )
})


test_that("plot_scores() returns a plotly plot with color and shape aesthetics set from the `SummarizedExperiment` metadata", {
  # create objects for running tests
  se <- random_se()
  se <- set_aesthetics(se, variable = "condition", values = c("control" = "pink", "treat" = "black"), aesthetic = "color")
  se <- set_aesthetics(se, variable = "batch", values = c("1" = 12, "2" = 8), aesthetic = "shape")
  SummarizedExperiment::assay(se, "vld") <- suppressMessages(DESeq2::varianceStabilizingTransformation(SummarizedExperiment::assay(se, "counts")))
  pca <- perform_pca(se, assay = "vld")

  # expect runs successfully
  expect_no_error(
    scores_plot <- plot_scores(pca, color = "condition", shape = "batch")
  )
  # expect returns plotly object
  expect_s3_class(scores_plot, "plotly")
  expect_s3_class(scores_plot, "htmlwidget")
})


test_that("plot_scores() returns a ggplot plot with color and shape aesthetics set from the `SummarizedExperiment` metadata", {
  # create objects for running tests
  se <- random_se()
  se <- set_aesthetics(se, variable = "condition", values = c("control" = "pink", "treat" = "black"), aesthetic = "color")
  se <- set_aesthetics(se, variable = "batch", values = c("1" = 12, "2" = 8), aesthetic = "shape")
  SummarizedExperiment::assay(se, "vld") <- suppressMessages(DESeq2::varianceStabilizingTransformation(SummarizedExperiment::assay(se, "counts")))
  pca <- perform_pca(se, assay = "vld")

  # expect runs successfully
  expect_no_error(
    scores_plot <- plot_scores(pca, plot_type = "ggplot", color = "condition", shape = "batch")
  )
  # expect returns ggplot object
  expect_s3_class(scores_plot, "gg")
  expect_s3_class(scores_plot, "ggplot")
  # test plot is the same
  vdiffr::expect_doppelganger(
    "pca scores color shape ggplot",
    scores_plot
  )
})


test_that("plot_scree() returns a plotly scree plot of cumulative and individual proportion variance explained for each PC by default", {
  # create objects for running tests
  se <- random_se()
  SummarizedExperiment::assay(se, "vld") <- suppressMessages(DESeq2::varianceStabilizingTransformation(SummarizedExperiment::assay(se, "counts")))
  pca <- perform_pca(se, assay = "vld")

  # expect runs successfully
  expect_no_error(
    scree_plot <- plot_scree(pca)
  )
  # returns plotly object
  expect_s3_class(scree_plot, "plotly")
  expect_s3_class(scree_plot, "htmlwidget")
})


test_that("plot_scree() returns a plotly scree plot of cumulative and individual proportion variance explained for each PC by default when `plot_type = 'plotly'`", {
  # create objects for running tests
  se <- random_se()
  SummarizedExperiment::assay(se, "vld") <- suppressMessages(DESeq2::varianceStabilizingTransformation(SummarizedExperiment::assay(se, "counts")))
  pca <- perform_pca(se, assay = "vld")

  # expect runs successfully
  expect_no_error(
    scree_plot <- plot_scree(pca, plot_type = "plotly")
  )
  # returns plotly object
  expect_s3_class(scree_plot, "plotly")
  expect_s3_class(scree_plot, "htmlwidget")
})


test_that("plot_scree() returns a ggplot2 scree plot of cumulative and individual proportion variance explained for each PC by default when `plot_type = 'ggplot'`", {
  # create objects for running tests
  se <- random_se()
  SummarizedExperiment::assay(se, "vld") <- suppressMessages(DESeq2::varianceStabilizingTransformation(SummarizedExperiment::assay(se, "counts")))
  pca <- perform_pca(se, assay = "vld")

  # expect runs successfully
  expect_no_error(
    scree_plot <- plot_scree(pca, plot_type = "ggplot")
  )
  # expect returns ggplot object
  expect_s3_class(scree_plot, "gg")
  expect_s3_class(scree_plot, "ggplot")
  # test plot is the same
  vdiffr::expect_doppelganger(
    "pca scree ggplot",
    scree_plot
  )
})
