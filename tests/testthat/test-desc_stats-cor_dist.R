test_that("plot_cor_heatmap() throws `error_no_assay` when assay not in SummarizedExperiment", {
  # create objects for running tests
  se <- random_se()

  # expect runs successfully
  expect_error(
    cor_heatmap_plot <- plot_cor_heatmap(se, assay = "vsd", method = "pearson"),
    class = "error_no_assay"
  )
})


test_that("plot_cor_heatmap() returns a ComplexHeatmap for Spearman correlation with top and left annotations by default", {
  # create objects for running tests
  se <- random_se()
  se <- set_aesthetics(se, variable = "condition", values = c("control" = "pink", "treat" = "black"), aesthetic = "color")
  se <- set_aesthetics(se, variable = "batch", values = c("1" = "blue", "2" = "green"), aesthetic = "color")

  # expect runs successfully
  expect_no_error(
    cor_heatmap_plot <- plot_cor_heatmap(se, assay = "counts", method = "spearman")
  )
  # expect returns Heatmap object
  expect_s4_class(cor_heatmap_plot, "Heatmap")
  # test plot is the same
  vdiffr::expect_doppelganger(
    "cor heatmap spearman default",
    cor_heatmap_plot
  )
})


test_that("plot_cor_heatmap() returns a ComplexHeatmap for Spearman correlation with top and left annotations when `annotations = 'both'`", {
  # create objects for running tests
  se <- random_se()
  se <- set_aesthetics(se, variable = "condition", values = c("control" = "pink", "treat" = "black"), aesthetic = "color")
  se <- set_aesthetics(se, variable = "batch", values = c("1" = "blue", "2" = "green"), aesthetic = "color")

  # expect runs successfully
  expect_no_error(
    cor_heatmap_plot <- plot_cor_heatmap(se, assay = "counts", method = "spearman", annotations = "both")
  )
  # expect returns Heatmap object
  expect_s4_class(cor_heatmap_plot, "Heatmap")
  # test plot is the same
  vdiffr::expect_doppelganger(
    "cor heatmap spearman annotations both",
    cor_heatmap_plot
  )
})


test_that("plot_cor_heatmap() returns a ComplexHeatmap for Spearman correlation with top annotations when `annotations = 'top'`", {
  # create objects for running tests
  se <- random_se()
  se <- set_aesthetics(se, variable = "condition", values = c("control" = "pink", "treat" = "black"), aesthetic = "color")
  se <- set_aesthetics(se, variable = "batch", values = c("1" = "blue", "2" = "green"), aesthetic = "color")

  # expect runs successfully
  expect_no_error(
    cor_heatmap_plot <- plot_cor_heatmap(se, assay = "counts", method = "spearman", annotations = "top")
  )
  # expect returns Heatmap object
  expect_s4_class(cor_heatmap_plot, "Heatmap")
  # test plot is the same
  vdiffr::expect_doppelganger(
    "cor heatmap spearman annotations top",
    cor_heatmap_plot
  )
})


test_that("plot_cor_heatmap() returns a ComplexHeatmap for Spearman correlation with left annotations when `annotations = 'left'`", {
  # create objects for running tests
  se <- random_se()
  se <- set_aesthetics(se, variable = "condition", values = c("control" = "pink", "treat" = "black"), aesthetic = "color")
  se <- set_aesthetics(se, variable = "batch", values = c("1" = "blue", "2" = "green"), aesthetic = "color")

  # expect runs successfully
  expect_no_error(
    cor_heatmap_plot <- plot_cor_heatmap(se, assay = "counts", method = "spearman", annotations = "left")
  )
  # expect returns Heatmap object
  expect_s4_class(cor_heatmap_plot, "Heatmap")
  # test plot is the same
  vdiffr::expect_doppelganger(
    "cor heatmap spearman annotations left",
    cor_heatmap_plot
  )
})


test_that("plot_cor_heatmap() returns a ComplexHeatmap for Spearman correlation with no annotations when `annotations = 'none'`", {
  # create objects for running tests
  se <- random_se()
  se <- set_aesthetics(se, variable = "condition", values = c("control" = "pink", "treat" = "black"), aesthetic = "color")
  se <- set_aesthetics(se, variable = "batch", values = c("1" = "blue", "2" = "green"), aesthetic = "color")

  # expect runs successfully
  expect_no_error(
    cor_heatmap_plot <- plot_cor_heatmap(se, assay = "counts", method = "spearman", annotations = "none")
  )
  # expect returns Heatmap object
  expect_s4_class(cor_heatmap_plot, "Heatmap")
  # test plot is the same
  vdiffr::expect_doppelganger(
    "cor heatmap spearman annotations none",
    cor_heatmap_plot
  )
})


test_that("plot_cor_heatmap() returns a ComplexHeatmap for Spearman correlation with diagonal retained when `diagonal = 'as_is'`", {
  # create objects for running tests
  se <- random_se()
  se <- set_aesthetics(se, variable = "condition", values = c("control" = "pink", "treat" = "black"), aesthetic = "color")
  se <- set_aesthetics(se, variable = "batch", values = c("1" = "blue", "2" = "green"), aesthetic = "color")

  # expect runs successfully
  expect_no_error(
    cor_heatmap_plot <- plot_cor_heatmap(se, assay = "counts", method = "spearman", diagonal = "as_is")
  )
  # expect returns Heatmap object
  expect_s4_class(cor_heatmap_plot, "Heatmap")
  # test plot is the same
  vdiffr::expect_doppelganger(
    "cor heatmap spearman diagonal as_is",
    cor_heatmap_plot
  )
})


test_that("plot_cor_heatmap() returns a ComplexHeatmap for pearson correlation", {
  # create objects for running tests
  se <- random_se()
  se <- set_aesthetics(se, variable = "condition", values = c("control" = "pink", "treat" = "black"), aesthetic = "color")
  se <- set_aesthetics(se, variable = "batch", values = c("1" = "blue", "2" = "green"), aesthetic = "color")
  SummarizedExperiment::assay(se, "vsd") <- suppressMessages(DESeq2::varianceStabilizingTransformation(SummarizedExperiment::assay(se, "counts")))

  # expect runs successfully
  expect_no_error(
    cor_heatmap_plot <- plot_cor_heatmap(se, assay = "vsd", method = "pearson")
  )
  # expect returns Heatmap object
  expect_s4_class(cor_heatmap_plot, "Heatmap")
  # test plot is the same
  vdiffr::expect_doppelganger(
    "cor heatmap pearson default",
    cor_heatmap_plot
  )
})


test_that("plot_cor_heatmap() returns a ComplexHeatmap for Spearman correlation when SummarizedExperiment has no aesthetics", {
  # create objects for running tests
  se <- random_se()
  SummarizedExperiment::assay(se, "vsd") <- suppressMessages(DESeq2::varianceStabilizingTransformation(SummarizedExperiment::assay(se, "counts")))

  # expect runs successfully
  expect_no_error(
    cor_heatmap_plot <- plot_cor_heatmap(se, assay = "vsd", method = "spearman")
  )
  # expect returns Heatmap object
  expect_s4_class(cor_heatmap_plot, "Heatmap")
  # test plot is the same
  vdiffr::expect_doppelganger(
    "cor heatmap spearman no aesthetics",
    cor_heatmap_plot
  )
})


test_that("plot_dist_heatmap() throws `error_no_assay` when assay not in SummarizedExperiment", {
  # create objects for running tests
  se <- random_se()

  # expect runs successfully
  expect_error(
    dist_heatmap_plot <- plot_dist_heatmap(se, assay = "vsd", method = "euclidean"),
    class = "error_no_assay"
  )
})


test_that("plot_dist_heatmap() returns a ComplexHeatmap for euclidean distance with top and left annotations by default", {
  # create objects for running tests
  se <- random_se()
  se <- set_aesthetics(se, variable = "condition", values = c("control" = "pink", "treat" = "black"), aesthetic = "color")
  se <- set_aesthetics(se, variable = "batch", values = c("1" = "blue", "2" = "green"), aesthetic = "color")

  # expect runs successfully
  expect_no_error(
    dist_heatmap_plot <- plot_dist_heatmap(se, assay = "counts", method = "euclidean")
  )
  # expect returns Heatmap object
  expect_s4_class(dist_heatmap_plot, "Heatmap")
  # test plot is the same
  vdiffr::expect_doppelganger(
    "dist heatmap euclidean default",
    dist_heatmap_plot
  )
})


test_that("plot_dist_heatmap() returns a ComplexHeatmap for euclidean distance with top and left annotations when `annotations = 'both'`", {
  # create objects for running tests
  se <- random_se()
  se <- set_aesthetics(se, variable = "condition", values = c("control" = "pink", "treat" = "black"), aesthetic = "color")
  se <- set_aesthetics(se, variable = "batch", values = c("1" = "blue", "2" = "green"), aesthetic = "color")

  # expect runs successfully
  expect_no_error(
    dist_heatmap_plot <- plot_dist_heatmap(se, assay = "counts", method = "euclidean", annotations = "both")
  )
  # expect returns Heatmap object
  expect_s4_class(dist_heatmap_plot, "Heatmap")
  # test plot is the same
  vdiffr::expect_doppelganger(
    "dist heatmap euclidean annotations both",
    dist_heatmap_plot
  )
})


test_that("plot_dist_heatmap() returns a ComplexHeatmap for euclidean distance with top annotations when `annotations = 'top'`", {
  # create objects for running tests
  se <- random_se()
  se <- set_aesthetics(se, variable = "condition", values = c("control" = "pink", "treat" = "black"), aesthetic = "color")
  se <- set_aesthetics(se, variable = "batch", values = c("1" = "blue", "2" = "green"), aesthetic = "color")

  # expect runs successfully
  expect_no_error(
    dist_heatmap_plot <- plot_dist_heatmap(se, assay = "counts", method = "euclidean", annotations = "top")
  )
  # expect returns Heatmap object
  expect_s4_class(dist_heatmap_plot, "Heatmap")
  # test plot is the same
  vdiffr::expect_doppelganger(
    "dist heatmap euclidean annotations top",
    dist_heatmap_plot
  )
})


test_that("plot_dist_heatmap() returns a ComplexHeatmap for euclidean distance with left annotations when `annotations = 'left'`", {
  # create objects for running tests
  se <- random_se()
  se <- set_aesthetics(se, variable = "condition", values = c("control" = "pink", "treat" = "black"), aesthetic = "color")
  se <- set_aesthetics(se, variable = "batch", values = c("1" = "blue", "2" = "green"), aesthetic = "color")

  # expect runs successfully
  expect_no_error(
    dist_heatmap_plot <- plot_dist_heatmap(se, assay = "counts", method = "euclidean", annotations = "left")
  )
  # expect returns Heatmap object
  expect_s4_class(dist_heatmap_plot, "Heatmap")
  # test plot is the same
  vdiffr::expect_doppelganger(
    "dist heatmap euclidean annotations left",
    dist_heatmap_plot
  )
})


test_that("plot_dist_heatmap() returns a ComplexHeatmap for euclidean distance with no annotations when `annotations = 'none'`", {
  # create objects for running tests
  se <- random_se()
  se <- set_aesthetics(se, variable = "condition", values = c("control" = "pink", "treat" = "black"), aesthetic = "color")
  se <- set_aesthetics(se, variable = "batch", values = c("1" = "blue", "2" = "green"), aesthetic = "color")

  # expect runs successfully
  expect_no_error(
    dist_heatmap_plot <- plot_dist_heatmap(se, assay = "counts", method = "euclidean", annotations = "none")
  )
  # expect returns Heatmap object
  expect_s4_class(dist_heatmap_plot, "Heatmap")
  # test plot is the same
  vdiffr::expect_doppelganger(
    "dist heatmap euclidean annotations none",
    dist_heatmap_plot
  )
})


test_that("plot_dist_heatmap() returns a ComplexHeatmap for euclidean distance with diagonal retained when `diagonal = 'as_is'`", {
  # create objects for running tests
  se <- random_se()
  se <- set_aesthetics(se, variable = "condition", values = c("control" = "pink", "treat" = "black"), aesthetic = "color")
  se <- set_aesthetics(se, variable = "batch", values = c("1" = "blue", "2" = "green"), aesthetic = "color")

  # expect runs successfully
  expect_no_error(
    dist_heatmap_plot <- plot_dist_heatmap(se, assay = "counts", method = "euclidean", diagonal = "as_is")
  )
  # expect returns Heatmap object
  expect_s4_class(dist_heatmap_plot, "Heatmap")
  # test plot is the same
  vdiffr::expect_doppelganger(
    "dist heatmap euclidean diagonal as_is",
    dist_heatmap_plot
  )
})


test_that("plot_dist_heatmap() returns a ComplexHeatmap for manhattan distance", {
  # create objects for running tests
  se <- random_se()
  se <- set_aesthetics(se, variable = "condition", values = c("control" = "pink", "treat" = "black"), aesthetic = "color")
  se <- set_aesthetics(se, variable = "batch", values = c("1" = "blue", "2" = "green"), aesthetic = "color")
  SummarizedExperiment::assay(se, "vsd") <- suppressMessages(DESeq2::varianceStabilizingTransformation(SummarizedExperiment::assay(se, "counts")))

  # expect runs successfully
  expect_no_error(
    dist_heatmap_plot <- plot_dist_heatmap(se, assay = "vsd", method = "manhattan")
  )
  # expect returns Heatmap object
  expect_s4_class(dist_heatmap_plot, "Heatmap")
  # test plot is the same
  vdiffr::expect_doppelganger(
    "dist heatmap manhattan default",
    dist_heatmap_plot
  )
})


test_that("plot_dist_heatmap() returns a ComplexHeatmap for euclidean distance when SummarizedExperiment has no aesthetics", {
  # create objects for running tests
  se <- random_se()
  SummarizedExperiment::assay(se, "vsd") <- suppressMessages(DESeq2::varianceStabilizingTransformation(SummarizedExperiment::assay(se, "counts")))

  # expect runs successfully
  expect_no_error(
    dist_heatmap_plot <- plot_dist_heatmap(se, assay = "vsd", method = "euclidean")
  )
  # expect returns Heatmap object
  expect_s4_class(dist_heatmap_plot, "Heatmap")
  # test plot is the same
  vdiffr::expect_doppelganger(
    "dist heatmap euclidean no aesthetics",
    dist_heatmap_plot
  )
})

