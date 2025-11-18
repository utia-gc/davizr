# test-diff_dream.R
# Test dream pipeline

test_that("run_dream() throws `error_no_formula` when design formula not set in `SummarizedExperiment`", {
  # construct SummarizedExperiment w/ no design formula
  n_genes <- 250L
  se <- vp_de_se()[1:n_genes, ]
  # make contrasts matrix
  contrasts_matrix <- variancePartition::makeContrastsDream(
    formula = as.formula("~ 0 + disease + (1 | individual)"),
    data = SummarizedExperiment::colData(se),
    contrasts = c(
      "disease1_vs_disease0" = "disease1 - disease0"
    )
  )

  # test that dream workflow throws error when `SummarizedExperiment` has no design formula
  expect_error(
    dream_wf_result <- run_dream(se, contrasts_matrix = contrasts_matrix),
    class = "error_no_formula"
  )
})


test_that("run_dream() transforms w/ dream-voom, fits LMM w/ dream, and computes moderated stats w/ eBayes for mixed model formula", {
  # construct SummarizedExperiment w/ mixed model formula
  n_genes <- 250L
  se <- vp_de_se()[1:n_genes, ]
  se <- set_formula(se, "design", as.formula("~ 0 + disease + (1 | individual)"))
  # make contrasts matrix
  contrasts_matrix <- variancePartition::makeContrastsDream(
    formula = get_formula(se, "design"),
    data = SummarizedExperiment::colData(se),
    contrasts = c(
      "disease1_vs_disease0" = "disease1 - disease0"
    )
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

  # test that differential expression analysis results can be extracted with `topTable()`
  expect_no_condition(
    disease_de_result <- variancePartition::topTable(
      dream_wf_result$fit_mm,
      coef = "disease1_vs_disease0",
      number = 1,
      sort.by = "B"
    )
  )
  # test that a t-test was run
  expect_contains(colnames(disease_de_result), "t")
  # test that TDRKH is most significant gene
  expect_equal(rownames(disease_de_result[1, ]), "ENST00000525790.1 gene=TDRKH")
})


test_that("run_dream() transforms w/ dream-voom, fits LMM w/ dream, and computes moderated stats w/ eBayes for more complex mixed model formula", {
  # construct SummarizedExperiment w/ mixed model formula
  n_genes <- 250L
  se <- vp_de_se()[1:n_genes, ]
  se <- set_formula(se, "design", as.formula("~ 0 + disease_subtype + sex + (1 | individual)"))
  # make contrasts matrix
  contrasts_matrix <- variancePartition::makeContrastsDream(
    formula = get_formula(se, "design"),
    data = SummarizedExperiment::colData(se),
    contrasts = c(
      "subtype1_vs_subtype0" = "disease_subtype1 - disease_subtype0",
      "subtype2_vs_subtype0" = "disease_subtype2 - disease_subtype0",
      "diseaseavg_vs_healthy" = "(disease_subtype2 + disease_subtype1)/2 - disease_subtype0"
    )
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

  # test that differential expression analysis results can be extracted with `topTable()`
  # test for disease average vs healthy
  expect_no_condition(
    disease_de_result <- variancePartition::topTable(
      dream_wf_result$fit_mm,
      coef = "diseaseavg_vs_healthy",
      number = 1,
      sort.by = "B"
    )
  )
  # test that a t-test was run
  expect_contains(colnames(disease_de_result), "t")
  # test that TDRKH is most significant gene
  expect_equal(rownames(disease_de_result[1, ]), "ENST00000525790.1 gene=TDRKH")

  # test for any disease vs healthy differences with an F test
  expect_no_condition(
    disease_de_result <- variancePartition::topTable(
      dream_wf_result$fit_mm,
      coef = c("subtype1_vs_subtype0", "subtype2_vs_subtype0"),
      number = 1,
      sort.by = "B"
    )
  )
  # test that a t-test was run
  expect_contains(colnames(disease_de_result), c("F", "F.std"))
  # test that TDRKH is most significant gene
  expect_equal(rownames(disease_de_result[1, ]), "ENST00000525790.1 gene=TDRKH")
})
