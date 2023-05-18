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
