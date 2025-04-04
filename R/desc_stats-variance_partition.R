#' Performance variance partitioning analysis
#'
#' @inheritParams perform_pca
#' @param formula_name A character name for the formula in formula models slot
#'   of `SummarizedExperiment` metadata.
#' @param ... Additional parameters passed to
#'   [variancePartition::fitExtractVarPartModel]. The `exprObj`, `formula`, and
#'   `data` arguments are already set from the `SummarizedExperiment`.
#'
#' @return A `varPartResults` from [variancePartition::fitExtractVarPartModel].
#' @export
perform_variance_partition_analysis <- function(se, assay, formula_name, ...) {
  # fail if assay does not exist in SummarizedExperiment
  if (!assay %in% SummarizedExperiment::assayNames(se)) {
    abort_no_assay(assay, "se")
  }
  # fail if formula does not exist in SummarizedExperiment
  if (!formula_name %in% names(S4Vectors::metadata(se)[["modeling"]][["formulas"]])) {
    abort_no_formula(formula_name, "se")
  }

  # make a list of variance partitioning arguments
  vp_args <- list(
    exprObj = SummarizedExperiment::assay(se, "vsd"),
    formula = get_formula(se, formula_name),
    data = get_samples_data(se, rownames_var = NULL)
  )
  additional_vp_args <- list(...)

  # run the variance partitioning analysis
  variance_fractions <- do.call(
    variancePartition::fitExtractVarPartModel,
    c(vp_args, additional_vp_args)
  )
  variance_fractions <- variancePartition::sortCols(variance_fractions)

  return(variance_fractions)
}
