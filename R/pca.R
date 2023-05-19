#' Make a scree plot for a PCA DFrame
#'
#' @param pca_df A `Dframe` of PCA data
#'
#' @return A ggplot object scree plot
#' @export
plot_scree <- function(pca_df) {
  pca_df %>%
    get_variance_explained() %>%
    ggplot2::ggplot(ggplot2::aes(x = principal_component)) +
    ggplot2::geom_col(ggplot2::aes(y = proportion_var)) +
    ggplot2::geom_line(ggplot2::aes(y = cumulative_proportion_var)) +
    ggplot2::geom_point(ggplot2::aes(y = cumulative_proportion_var)) +
    ggplot2::scale_x_continuous(
      name = "Principal component",
      breaks = scales::breaks_pretty()
    ) +
    ggplot2::scale_y_continuous(
      name = "Proportion variance",
      sec.axis = ggplot2::sec_axis(trans = ~.*1, name = "Cumulative proportion variance")
    )
}

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
  S4Vectors::metadata(pca_df)$variance_explained <- generate_variance_explained_summary(pca)

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

#' Access variance explained from a PCA DFrame
#'
#' @description An accessor method for getting the variance explained data from
#'   a PCA DFrame such as that produced by [perform_pca()]
#'
#' @inheritParams get_prcomp
#'
#' @return A data frame with the proportion of variance explained and cumulative
#'   proportion of variance explained for each principal component.
#' @export
get_variance_explained <- function(pca_dframe) {
  S4Vectors::metadata(pca_dframe)$variance_explained
}

#' Generate a summary of variance explained by PCA
#'
#' @param prcomp A `prcomp` object
#'
#' @return A data frame with the proportion of variance explained and cumulative
#'   proportion of variance explained for each principal component.
#' @noRd
generate_variance_explained_summary <- function(prcomp) {
  data.frame(
    principal_component = seq_len(length(prcomp$sdev)),
    proportion_var = prcomp$sdev^2 / sum(prcomp$sdev^2),
    cumulative_proportion_var = cumsum(prcomp$sdev^2 / sum(prcomp$sdev^2))
  )
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
