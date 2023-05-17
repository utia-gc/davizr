test_that("extract_results_summary_data() gets summary data from results", {
  results <- extract_shrunken_results(example_dds(), name = "condition_B_vs_A", alpha = 0.05)

  expect_equal(
    extract_results_summary_data(results),
    results_summary_data()
  )
})
