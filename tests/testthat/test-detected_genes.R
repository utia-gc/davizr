test_that("count_experiment_detected_genes() counts the correct number of detected genes", {
  expect_equal(count_experiment_detected_genes(DESeq2::estimateSizeFactors(random_dds())), 100L)
})

test_that("count_experiment_detected_genes() counts the correct number of detected genes with high threshold", {
  expect_equal(
    count_experiment_detected_genes(DESeq2::estimateSizeFactors(random_dds()), threshold = 25, covariates = "condition"),
    95L
  )
})
