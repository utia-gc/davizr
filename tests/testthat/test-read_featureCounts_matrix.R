test_that("read_featureCounts_matrix() returns counts matrix from vector of file paths", {
  featureCounts_paths <- file.path(
    test_path("data"),
    list.files(test_path("data"), pattern = "^sample_[12]_counts\\.txt")
  )
  expect_equal(
    read_featureCounts_matrix(featureCounts_paths),
    featureCounts_matrix()
  )
})
