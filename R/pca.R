#' Perform PCA and return data for transformed data
#'
#' @inheritParams DESeq2::plotPCA
#'
#' @return A `DFrame` with a `prcomp` object in the metadata.
#' @export
perform_pca <- function(object, ntop = 500) {
  gene_indexes <- which_pca_genes(object, ntop)
  assay_matrix <- prepare_pca_assay_matrix(object, gene_indexes)

  pca <- stats::prcomp(assay_matrix)

  x <- pca$x %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "sample_name")

  col_data <- object %>%
    SummarizedExperiment::colData() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "sample_name")

  pca_df <- dplyr::full_join(
    col_data,
    x,
    by = "sample_name"
  ) %>%
    tibble::column_to_rownames(var = "sample_name") %>%
    S4Vectors::DataFrame()
  S4Vectors::metadata(pca_df)$prcomp <- pca

  pca_df
}

#' Access prcomp from a PCA DFrame
#'
#' @description An accessor method for getting the `prcomp` object from a PCA DFrame
#'   such as that produced by [perform_pca()]
#'
#' @param pca_dframe A `Dframe` of PCA data
#'
#' @return A `prcomp` object
#' @export
get_prcomp <- function(pca_dframe) {
  S4Vectors::metadata(pca_dframe)$prcomp
}

#' Prepare input matrix for PCA
#'
#' @inheritParams perform_pca
#' @param indexes Vector of gene indexes to subset assay matrix
#'
#' @return Matrix of assay data transposed with samples as rows and genes as columns.
#' @noRd
prepare_pca_assay_matrix <- function(object, indexes) {
  t(SummarizedExperiment::assay(object)[indexes, ])
}

#' Get indexes of which genes to select for PCA
#'
#' @inheritParams perform_pca
#'
#' @return A numeric vector of indices
#' @noRd
which_pca_genes <- function(object, ntop) {
  row_vars <- MatrixGenerics::rowVars(SummarizedExperiment::assay(object))

  order(row_vars, decreasing = TRUE)[seq_len(min(ntop, length(row_vars)))]
}
