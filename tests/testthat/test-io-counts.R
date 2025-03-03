test_that("read_counts_file() returns counts matrix from of file path", {
  counts_path <- file.path(
    test_path("data"),
    "rnaseq_nf-test_counts.csv"
  )
  counts <- matrix(
    data = c(0, 0, 1, 0, 0, 0, 0, 4, 3, 0),
    nrow = 5,
    ncol = 2,
    dimnames = list(
      c("YAL069W", "YAL068W-A", "YAL068C", "YAL067W-A", "YAL067C"),
      c("srr1066657", "srr6924569")
    )
  )

  expect_equal(
    read_counts_file(counts_path),
    counts
  )
})
