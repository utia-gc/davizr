#' Plot detected gene counts for each library
#'
#' @description Generate a bar plot with the count of detected genes for each library.
#'
#'   If detected gene counts are not present in `dds` (i.e. no `detected_genes` column in `dds` column data),
#'   detected genes are counted by [count_library_detected_genes()].
#'
#' @param dds A `DESeqDataSet`.
#' @inheritParams count_library_detected_genes
#'
#' @return A ggplot object bar plot of counts of detected genes for each library.
#' @export
plot_detected_genes <- function(dds, threshold = 5) {
  if (!has_detected_genes(dds)) {
    dds$detected_genes <- count_library_detected_genes(dds, threshold = threshold)
  }

  dds %>%
    SummarizedExperiment::colData() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "sample_name") %>%
    ggplot2::ggplot() +
    ggplot2::geom_bar(ggplot2::aes(x = sample_name, y = detected_genes), stat = "identity")
}

#' Count detected genes for each library
#'
#' @param dds A `DESeqDataSet`.
#' @param threshold Threshold for minimum normalized counts.
#'
#' @return A named vector of counts of detected genes for each sample.
#' @export
count_library_detected_genes <- function(dds, threshold = 5) {
  if (!has_sizeFactor(dds)) {
    dds <- DESeq2::estimateSizeFactors(dds)
  }

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
  if (!has_sizeFactor(dds)) {
    dds <- DESeq2::estimateSizeFactors(dds)
  }
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

#' Check if a `DESeqDataSet` has a `sizeFactor`
#'
#' @param dds A `DESeqDataSet`.
#'
#' @return boolean for whether `dds` has a `sizeFactor` column.
#' @noRd
has_sizeFactor <- function(dds) {
  "sizeFactor" %in% colnames(SummarizedExperiment::colData(dds))
}

#' Check if a `DESeqDataSet` has detected genes counts
#'
#' @param dds
#'
#' @return boolean for whether `dds` has a `detected_genes` column.
#' @noRd
has_detected_genes <- function(dds) {
  "detected_genes" %in% colnames(SummarizedExperiment::colData(dds))
}
