test_that("count_experiment_detected_genes() counts the correct number of detected genes", {
  expect_equal(count_experiment_detected_genes(DESeq2::estimateSizeFactors(random_dds())), 100L)
})

test_that("count_experiment_detected_genes() counts the correct number of detected genes with high threshold", {
  expect_equal(
    count_experiment_detected_genes(DESeq2::estimateSizeFactors(random_dds()), threshold = 25, covariates = "condition"),
    95L
  )
})

test_that("count_library_detected_genes() counts the correct number of detected genes in each sample", {
  expect_equal(
    count_library_detected_genes(DESeq2::estimateSizeFactors(random_dds())),
    library_detected_genes()
  )
})

test_that("count_library_detected_genes() counts the correct number of detected genes in each sample with low threshold", {
  expect_equal(
    count_library_detected_genes(DESeq2::estimateSizeFactors(random_dds()), threshold = 1),
    library_detected_genes_low_thresh()
  )
})
