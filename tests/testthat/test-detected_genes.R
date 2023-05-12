test_that("count_experiment_detected_genes() counts detected genes when no `sizeFactors`", {
  dds <- random_dds()

  expect_equal(
    count_experiment_detected_genes(dds),
    100L
  )
})

test_that("count_experiment_detected_genes() counts detected genes when `sizeFactors` available", {
  dds <- DESeq2::estimateSizeFactors(random_dds())

  expect_equal(
    count_experiment_detected_genes(dds),
    100L
  )
})

test_that("count_experiment_detected_genes() counts detected genes with high threshold", {
  dds <- DESeq2::estimateSizeFactors(random_dds())

  expect_equal(
    count_experiment_detected_genes(dds, threshold = 25, covariates = "condition"),
    95L
  )
})

test_that("count_library_detected_genes() counts detected genes in each library when no `sizeFactors`", {
  dds <- random_dds()

  expect_equal(
    count_library_detected_genes(dds),
    library_detected_genes()
  )
})

test_that("count_library_detected_genes() counts detected genes in each library when `sizeFactors` available", {
  dds <- DESeq2::estimateSizeFactors(random_dds())

  expect_equal(
    count_library_detected_genes(dds),
    library_detected_genes()
  )
})

test_that("count_library_detected_genes() counts detected genes in each library with low threshold", {
  dds <- DESeq2::estimateSizeFactors(random_dds())

  expect_equal(
    count_library_detected_genes(dds, threshold = 1),
    library_detected_genes_low_thresh()
  )
})

test_that("plot_detected_genes() produces plot when detected_genes column is present", {
  dds <- random_dds()
  dds$detected_genes <- library_detected_genes()
  p <- plot_detected_genes(dds)

  vdiffr::expect_doppelganger(
    "detected genes", p
  )
})

test_that("plot_detected_genes() produces plot when detected_genes column is absent", {
  p <- plot_detected_genes(random_dds())

  vdiffr::expect_doppelganger(
    "detected genes", p
  )
})
