# test-diff_dream.R
# Test dream pipeline

test_that("run_dream() transforms w/ dream-voom, fits LMM w/ dream, and computes moderated stats w/ eBayes", {
  # construct SummarizedExperiment w/ mixed model formula
  n_genes <- 250L
  se <- vp_de_se()[1:n_genes, ]
  se <- set_formula(se, "design", as.formula("~ 0 + disease + (1 | individual)"))
  # make contrasts matrix
  contrasts_matrix <- variancePartition::makeContrastsDream(
    formula = get_formula(se, "design"),
    data = SummarizedExperiment::colData(se),
    contrasts = c("disease1 - disease0")
  )

  # test that dream workflow is successful
  expect_no_error(
    dream_wf_result <- run_dream(se, contrasts_matrix = contrasts_matrix)
  )

  # test that results list has expected names
  expect_named(dream_wf_result, c("voom_obj", "fit_mm", "se"))

  # test that voom object is the expected EList
  expect_s4_class(dream_wf_result[["voom_obj"]], "EList")
  expect_identical(ncol(dream_wf_result[["voom_obj"]]), 24L)
  expect_lte(nrow(dream_wf_result[["voom_obj"]]), n_genes)

  # test that SummarizedExperiment object is the expected SummarizedExperiment
  expect_s4_class(dream_wf_result[["se"]], "SummarizedExperiment")
  expect_identical(ncol(dream_wf_result[["se"]]), 24L)
  expect_lte(nrow(dream_wf_result[["se"]]), n_genes)
  # test that SummarizedExperiment has same dimensions as EList
  expect_identical(dim(dream_wf_result[["se"]]), dim(dream_wf_result[["voom_obj"]]))
  # test that SummarizedExperiment has dream-voom expression matrix
  expect_contains(SummarizedExperiment::assayNames(dream_wf_result[["se"]]), c("voom_dream"))

  # test that mixed model fit is expected output from dream + eBayes
  expect_s4_class(dream_wf_result[["fit_mm"]], "MArrayLM2")
  expect_lte(nrow(dream_wf_result[["se"]]), n_genes)
})
