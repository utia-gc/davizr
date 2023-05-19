test_that("perform_pca() returns DFrame with prcomp object from vst-transformed DESeqTransform", {
  vsd <- DESeq2::varianceStabilizingTransformation(example_dds())

  expect_equal(
    perform_pca(vsd),
    vsd_pca_data()
  )
})

test_that("perform_pca() returns DFrame with prcomp object from rlog-transformed DESeqTransform", {
  rld <- DESeq2::rlog(example_dds())

  expect_equal(
    perform_pca(rld),
    rld_pca_data()
  )
})

test_that("perform_pca() returns DFrame with prcomp object from 100 genes", {
  vsd <- DESeq2::varianceStabilizingTransformation(example_dds())

  expect_equal(
    perform_pca(vsd, ntop = 100),
    vsd_pca_data_100()
  )
})

test_that("get_prcomp() accesses the prcomp object in a PCA DFrame", {
  pca <- S4Vectors::metadata(vsd_pca_data_100())$prcomp

  expect_equal(
    get_prcomp(vsd_pca_data_100()),
    pca
  )
})

test_that("get_variance_explained() accesses the variance_explained df in a PCA DFrame", {
  vsd <- vsd_pca_data_100()
  variance_explained <- S4Vectors::metadata(vsd)$variance_explained

  expect_equal(
    get_variance_explained(vsd),
    variance_explained
  )
})

test_that("plot_scree() makes a scree plot from a PCA DFrame", {
  pca <- rld_pca_data()
  p <- plot_scree(pca)

  vdiffr::expect_doppelganger(
    "rld scree plot",
    p
  )
})

test_that("plot_biplot() plots first two principal components by default", {
  pca <- rld_pca_data()
  p <- plot_biplot(pca)

  vdiffr::expect_doppelganger(
    "rld biplot",
    p
  )
})

test_that("plot_biplot() plots supplied principal components", {
  pca <- rld_pca_data()
  p <- plot_biplot(pca, x = PC2, y = PC4)

  vdiffr::expect_doppelganger(
    "rld biplot PC4 vs PC2",
    p
  )
})

test_that("plot_biplot() applies aesthetic mappings to points", {
  pca <- rld_pca_data()
  p <- plot_biplot(pca, point_mapping = ggplot2::aes(color = sizeFactor, shape = condition))

  vdiffr::expect_doppelganger(
    "rld biplot colored sizeFactor shaped condition",
    p
  )
})
