#' Count detected genes for each library
#'
#' @param dds A `DESeqDataSet`.
#' @param threshold Threshold for minimum normalized counts.
#'
#' @return A named vector of counts of detected genes for each sample.
#' @export
count_library_detected_genes <- function(dds, threshold = 5) {
  base::colSums(
    DESeq2::counts(dds, normalized = TRUE) >= threshold
  )
}

#' Count detected genes for a whole experiment
#'
#' @description Count the number of genes with normalized counts above a specified threshold in the smallest condition set.
#'
#' @param dds A `DESeqDataSet`.
#' @param threshold Threshold for minimum normalized counts.
#' @inheritParams count_samples_each_condition
#'
#' @return Count of experiment-wide detected genes.
#' @export
count_experiment_detected_genes <- function(dds, threshold = 5, covariates = NULL) {
  counts <- DESeq2::counts(dds, normalized = TRUE)

  sum(
    base::rowSums(counts >= threshold) >= min(count_samples_each_condition(dds, covariates = covariates)$count)
  )
}

#' Count the number of samples in each condition
#'
#' @param dds A `DESeqDataSet`.
#' @param covariates A vector of covariates found in `dds` to count the number of samples within.
#'   If `NULL`, the covariates from the `dds` design formula will be used.
#'
#' @return A data frame with a count of the number of samples in each condition.
#' @export
count_samples_each_condition <- function(dds, covariates = NULL) {
  if (is.null(covariates)) {
    covariates <- get_design_covariates(dds)
  }

  dds %>%
    SummarizedExperiment::colData() %>%
    as.data.frame() %>%
    dplyr::group_by(dplyr::pick({{ covariates }})) %>%
    dplyr::summarize(count = dplyr::n(), .groups = 'keep')
}

#' Get covariates from design formulas
#'
#' @param dds A `DESeqDataSet`.
#'
#' @return A character vector of covariates found in the design formula.
#' @export
get_design_covariates <- function(dds) {
  as.character(DESeq2::design(dds))[2] %>%
    stringr::str_split_1(pattern = " \\+ ")
}
