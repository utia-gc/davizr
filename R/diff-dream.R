run_dream <- function(se, contrasts_matrix) {
  # get data for input to dream-voom
  formula = get_formula(se, "design")
  samples_data = se |>
    SummarizedExperiment::colData() |>
    as.data.frame()

  # convert SummarizedExperiment to DGEList and compute scaling factors
  dge <- se |>
    edgeR::SE2DGEList() |>
    edgeR::calcNormFactors()

  # transform counts to log2CPM with observation level weights
  voom_obj <- variancePartition::voomWithDreamWeights(
    counts = dge,
    formula = formula,
    data = samples_data,
    plot = FALSE,
    save.plot = TRUE
  )

  # fit per gene mixed models with dream method followed by computing moderated
  # statistics
  fit_mm <- variancePartition::dream(
    voom_obj,
    formula = get_formula(se, "design"),
    data = voom_obj[["targets"]],
    L = contrasts_matrix,
    hideErrorsInBackend = TRUE
  ) |>
    variancePartition::eBayes()

  # add dream-voom transformed expression values to SummarizedExperiment
  # voomWithDreamWeights drops genes where model does not converge, so these
  # genes must also be dropped from the SummarizedExperiment to keep dimensions
  # the same
  se <- se[rownames(voom_obj[["E"]]), , drop = FALSE]
  SummarizedExperiment::assay(se, "voom_dream") <- voom_obj[["E"]]

  # make list of objects to output
  dream_wf_result <- list(
    "voom_obj" = voom_obj,
    "fit_mm" = fit_mm,
    "se" = se
  )

  return(dream_wf_result)
}
