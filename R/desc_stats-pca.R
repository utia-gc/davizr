#' Perform principal components analysis (PCA) on an assay from
#' `SummarizedExperiment`
#'
#' @description Perform PCA for an assay from `SummarizedExperiment`, and return
#'   the results alongside sample data from the `SummarizedExperiment`.
#'
#' @details `perform_pca()` is essentially a convenience around running
#'   [stats::prcomp] on a particular assay matrix from a `SummarizedExperiment`,
#'   joining the PC scores to the sample data, and also returning some data
#'   about the proportion variance explained.
#'
#'   Note that there is no additional transformation of the expression data
#'   aside from subsetting the expression matrix to include only the `ntop` most
#'   variable genes and coercing the matrix to the expected shape for
#'   [stats::prcomp]. If you wish to center or scale data, for example, you must
#'   do so prior to running `perform_pca()` and store the transformed matrix in
#'   a new element in the `assay` slot of your `SummarizedExperiment`.
#'
#' @inheritParams write_se
#' @param assay A character vector for the assay from `SummarizedExperiment` to
#'   use as input to the PCA.
#' @param ntop An integer for the number of most variable genes to use as input
#'   to the PCA.
#'
#' @return A list of type `PCA`. Contains PC scores joined to sample data, a
#'   data.frame of proportion of variance explained by each PC, the `prcomp`
#'   object returned by [stats::prcomp], and the metadata slot from the
#'   `SummarizedExperiment`.
#' @export
perform_pca <- function(se, assay, ntop = 500) {
  # fail if assay does not exist in SummarizedExperiment
  if (!assay %in% SummarizedExperiment::assayNames(se)) {
    abort_no_assay(assay, "se")
  }

  # extract the matrix of expression values that will be used for PCA
  expr_matrix <- SummarizedExperiment::assay(se, assay)

  # subset the expression matrix and transpose it for PCA
  row_vars <- MatrixGenerics::rowVars(expr_matrix)
  highly_variable_gene_indexes <- order(row_vars, decreasing = TRUE)[seq_len(min(ntop, length(row_vars)))]
  expr_matrix <- t(expr_matrix[highly_variable_gene_indexes, ])

  # run the pca
  prcomp <- stats::prcomp(expr_matrix)

  # prepare scores matrix for joining with sample data
  scores <- as.data.frame(prcomp$x)
  scores[["sample"]] <- rownames(scores)

  # prepare sample data for joining with pc scores
  samples <- davizr::get_samples_data(se)

  # join scores and sample data
  sample_scores <- dplyr::full_join(
    samples,
    scores,
    by = dplyr::join_by("sample" == "sample")
  )
  rownames(sample_scores) <- sample_scores[["sample"]]

  # get the variance explained
  variance_explained <- explain_variance(prcomp)

  # construct the PCA object with all relevant info
  pca <- list(
    scores = sample_scores,
    variance_explained = variance_explained,
    prcomp = prcomp,
    metadata = S4Vectors::metadata(se)
  )
  pca <- structure(pca, class = "PCA")

  return(pca)
}


#' Compute the proportion of variance explained from a principal components
#' analysis
#'
#' @param prcomp A `prcomp` object, i.e. the object returned by [stats::prcomp].
#'
#' @return A data.frame of the proportion of variance explained and cumulative
#'   proportion of variance explained for each principal_component.
#' @noRd
explain_variance <- function(prcomp) {
  variance_explained <- data.frame(
    principal_component = paste0("PC", seq_len(length(prcomp$sdev))),
    proportion_var = prcomp$sdev^2 / sum(prcomp$sdev^2),
    cumulative_proportion_var = cumsum(prcomp$sdev^2 / sum(prcomp$sdev^2))
  )

  return(variance_explained)
}


#' Plot a 2D scatter plot of PCA scores
#'
#' @inheritParams plot_scree
#' @param x A character vector of the principal component to plot on the x-axis.
#' @param y A character vector of the principal component to plot on the y-axis.
#' @param margin_plot A character vector indicating the plot type to make in the
#'   margins.
#' @param label A character vector indicating the sample data variable to use
#'   for labeling samples. Labels don't show up on ggplots, but will be visible
#'   on the plotly tooltip.
#' @param color A character vector indicating the sample data variable to use
#'   for coloring the plot points.
#' @param shape A character vector indicating the sample data variable to use
#'   for shape of the plot points.
#'
#' @return A plot of PC scores.
#' @export
plot_scores <- function(
    pca,
    x = "PC1",
    y = "PC2",
    plot_type = c("plotly", "ggplot"),
    margin_plot = c("rug", "none"),
    label = "sample",
    color = NULL,
    shape = NULL
) {

  plot_type <- rlang::arg_match(plot_type)
  margin_plot <- rlang::arg_match(margin_plot)

  # make basic scores plot
  scores_plot <- ggplot2::ggplot(
    pca[["scores"]],
    ggplot2::aes(x = .data[[x]], y = .data[[y]])
  ) +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = make_percent_var_label(pca, x),
      y = make_percent_var_label(pca, y)
    )

  # add aesthetics that are set in options
  if (!is.null(label)) {
    scores_plot <- scores_plot +
      ggplot2::aes(label = .data[[label]])
  }
  if (!is.null(color)) {
    scores_plot <- scores_plot + ggplot2::aes(color = .data[[color]])

    color_values <- pca[["metadata"]][["aesthetics"]][["color"]][[color]]
    if (!is.null(color_values)) {
      scores_plot <- scores_plot + ggplot2::scale_color_manual(values = color_values)
    }
  }
  if (!is.null(shape)) {
    scores_plot <- scores_plot + ggplot2::aes(shape = .data[[shape]])

    shape_values <- pca[["metadata"]][["aesthetics"]][["shape"]][[shape]]
    if (!is.null(shape_values)) {
      scores_plot <- scores_plot + ggplot2::scale_shape_manual(values = shape_values)
    }
  }


  # add margin plot
  scores_plot <- switch(margin_plot,
    rug = scores_plot + ggplot2::geom_rug(),
    none = scores_plot
  )

  # coerce the plot to the requested format
  scores_plot <- switch(plot_type,
    plotly = plotly::ggplotly(scores_plot),
    ggplot = scores_plot
  )

  return(scores_plot)
}


#' Make an axis label to display the percent variance explained by a principal
#' component
#'
#' @inheritParams plot_scree
#' @param pc A character vector indicating the principal component to construct
#'   a label for.
#'
#' @return A character vector of the principal component and the percent
#'   variation it explains.
#' @noRd
make_percent_var_label <- function(pca, pc) {

  # get the percent variance explained for the specified principal component from the PCA object
  percent_variance <- pca[["variance_explained"]][pca[["variance_explained"]] == pc, "proportion_var", drop = TRUE] * 100

  # construct the label for the percent variance explained
  percent_variance_label <- stringr::str_glue("{pc}: {round(percent_variance, 2)}% variation")

  return(percent_variance_label)
}


#' Plot a PCA scree plot
#'
#' @param pca A `PCA` object returned by [perform_pca].
#' @param plot_type A character vector indicating the plot type to return.
#'
#' @return A scree plot.
#' @export
plot_scree <- function(pca, plot_type = c("plotly", "ggplot")) {

  plot_type <- rlang::arg_match(plot_type)

  # prepare data for plotting
  var_explained <- pca[["variance_explained"]]
  var_explained[["pc"]] <- as.integer(gsub("^PC", "", var_explained[["principal_component"]]))

  # make the scree plot
  scree_plot <- ggplot2::ggplot(
    var_explained,
    ggplot2::aes(x = .data[["pc"]])
  ) +
    ggplot2::geom_col(ggplot2::aes(y = proportion_var)) +
    ggplot2::geom_line(ggplot2::aes(y = cumulative_proportion_var)) +
    ggplot2::geom_point(ggplot2::aes(y = cumulative_proportion_var)) +
    ggplot2::scale_x_continuous(
      name = "Principal compoment",
      breaks = scales::breaks_pretty()
    ) +
    ggplot2::scale_y_continuous(
      name = "Proportion variance",
      sec.axis = ggplot2::sec_axis(trans = ~.*1, name = "Cumulative proportion variance")
    )

  # coerce the plot to the requested format
  scree_plot <- switch(plot_type,
    plotly = plotly::ggplotly(scree_plot),
    ggplot = scree_plot
  )

  return(scree_plot)
}
