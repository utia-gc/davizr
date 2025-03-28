#' Compute and plot pairwise sample--sample correlations for a
#' `SummarizedExperiment` `assay`
#'
#' @details If there are aesthetics set -- e.g. by [set_aesthetics] -- then
#'   heatmap annotations are automatically created from those aesthetics.
#'   Annotation position can be controlled by `annotations`.
#'
#' @inheritParams write_se
#' @param assay A character vector indicating the assay to compute correlation
#'   for.
#' @param method A character vector indicating the correlation method. Passed to
#'   [stats::cor]. Also sets a name for the heatmap unless a `name` argument is
#'   supplied.
#' @param ... Other arguments passed to [ComplexHeatmap::Heatmap].
#' @param annotations A character vector indicating which column/row annotations
#'   to include.
#' @param diagonal A character vector indicating what should happen to the
#'   matrix diagonal. `"remove"` plots the diagonal in white, effectively
#'   removing it from the plot. `"as_is"` plots the diagonal as is.
#'
#' @return A Heatmap plot of pairwise sample--sample correlations.
#' @export
plot_cor_heatmap <- function(
    se,
    assay,
    method = c("pearson", "kendall", "spearman"),
    ...,
    annotations = c("both", "top", "left", "none"),
    diagonal = c("remove", "as_is")
) {

  method <- rlang::arg_match(method)
  annotations <- rlang::arg_match(annotations)
  diagonal <- rlang::arg_match(diagonal)

  # compute correlations
  correlation_matrix <- compute_cor(se, assay, method)

  # set remove diagonal function
  rm_diag <- switch(
    diagonal,
    remove = remove_diagonal,
    as_is = NULL
  )

  # set heatmap name
  # if the name keyword variable is given as an argument, use that
  # otherwise, generate one given the correlation method
  heatmap_name <- ifelse(
    exists("name"),
    name,
    switch(
      method,
      pearson = "Pearson\nr",
      kendall = "Kendall\ntau",
      spearman = "Spearman\nrho"
    )
  )

  # top and left annotations
  top_anno <- switch(
    annotations,
    both = construct_heatmap_annotation(se, which = "column", legends_names = TRUE),
    top = construct_heatmap_annotation(se, which = "column", legends_names = TRUE),
    left = NULL,
    none = NULL
  )
  left_anno <- switch(
    annotations,
    both = construct_heatmap_annotation(se, which = "row", legends_names = FALSE),
    top = NULL,
    left = construct_heatmap_annotation(se, which = "row", legends_names = TRUE),
    none = NULL
  )

  # construct color mapping function
  color_map <- circlize::colorRamp2(
    c(min(correlation_matrix), mean(c(min(correlation_matrix), max(correlation_matrix))), max(correlation_matrix)),
    rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")[c(1, 4, 7)])
  )

  # plot correlation heatmap
  correlation_heatmap <- ComplexHeatmap::Heatmap(
    correlation_matrix,
    col = color_map,
    name = heatmap_name,
    cell_fun = rm_diag,
    top_annotation = top_anno,
    left_annotation = left_anno,
    ...
  )

  return(correlation_heatmap)
}


#' Compute and plot pairwise sample--sample distances for a
#' `SummarizedExperiment` `assay`
#'
#' @inheritParams plot_cor_heatmap
#'
#' @param method A character vector passed to [stats::dist]. Also sets a name for the heatmap unless a `name` argument is
#'   supplied.
#'
#' @return A Heatmap plot of pairwise sample--sample correlations.
#' @export
plot_dist_heatmap <- function(
    se,
    assay,
    method = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
    ...,
    annotations = c("both", "top", "left", "none"),
    diagonal = c("remove", "as_is")
) {

  method <- rlang::arg_match(method)
  annotations <- rlang::arg_match(annotations)
  diagonal <- rlang::arg_match(diagonal)

  # compute ditsances
  distance_matrix <- compute_dist(se, assay, method)

  # set remove diagonal function
  rm_diag <- switch(
    diagonal,
    remove = remove_diagonal,
    as_is = NULL
  )

  # set heatmap name
  # if the name keyword variable is given as an argument, use that
  # otherwise, generate one given the correlation method
  heatmap_name <- ifelse(exists("name"), name, method)

  # top and left annotations
  top_anno <- switch(
    annotations,
    both = construct_heatmap_annotation(se, which = "column", legends_names = TRUE),
    top = construct_heatmap_annotation(se, which = "column", legends_names = TRUE),
    left = NULL,
    none = NULL
  )
  left_anno <- switch(
    annotations,
    both = construct_heatmap_annotation(se, which = "row", legends_names = FALSE),
    top = NULL,
    left = construct_heatmap_annotation(se, which = "row", legends_names = TRUE),
    none = NULL
  )

  # construct color mapping function
  color_map <- circlize::colorRamp2(
    c(
      0,
      mean(c(min(distance_matrix), max(distance_matrix))),
      max(distance_matrix)
    ),
    rev(RColorBrewer::brewer.pal(n = 7, name = "Blues")[c(2, 4, 7)])
  )

  # plot correlation heatmap
  distance_heatmap <- ComplexHeatmap::Heatmap(
    distance_matrix,
    col = color_map,
    name = heatmap_name,
    cell_fun = rm_diag,
    top_annotation = top_anno,
    left_annotation = left_anno,
    ...
  )

  return(distance_heatmap)
}


#' Compute pairwise sample--sample correlations for a `SummarizedExperiment`
#' `assay`
#'
#' @inheritParams plot_cor_heatmap
#'
#' @return A matrix of pairwise sample correlations.
compute_cor <- function(se, assay, method) {

  # fail if assay does not exist in SummarizedExperiment
  if (!assay %in% SummarizedExperiment::assayNames(se)) {
    abort_no_assay(assay = assay, se = "se")
  }

  # extract the matrix of expression values that will be used for correlation analysis
  expr_matrix <- SummarizedExperiment::assay(se, assay)
  # compute correlations
  correlation_matrix <- stats::cor(expr_matrix, method = method)

  return(correlation_matrix)
}


#' Compute pairwise sample--sample distances for a `SummarizedExperiment`
#' `assay`
#'
#' @inheritParams plot_dist_heatmap
#'
#' @return A matrix of pairwise sample--sample distances.
compute_dist <- function(se, assay, method = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) {

  # fail if assay does not exist in SummarizedExperiment
  if (!assay %in% SummarizedExperiment::assayNames(se)) {
    abort_no_assay(assay = assay, se = "se")
  }

  # extract the matrix of expression values that will be used for distance analysis
  expr_matrix <- SummarizedExperiment::assay(se, assay)
  # compute correlations
  distance_matrix <- stats::dist(t(expr_matrix), method = method) |>
    as.matrix()

  return(distance_matrix)
}


#' Construct heatmap annotations from `SummarizedExperiment` aesthetics
#'
#' @description Creates heatmap annotation with
#' [ComplexHeatmap::HeatmapAnnotation] from aesthetics in a
#' `SummarizedExperiment`, e.g. those created by [set_aesthetics].
#'
#'
#' @inheritParams write_se
#' @inheritParams ComplexHeatmap::HeatmapAnnotation
#' @param legends_names A logical controlling whether or not to show legend and
#'   annotation name
#'
#' @return A `HeatmapAnnotation` if aesthetics are set in `se`, otherwise
#'   `NULL`.
#' @noRd
construct_heatmap_annotation <- function(se, which = c("column", "row"), legends_names = TRUE) {

  which <- rlang::arg_match(which)

  color_map <- S4Vectors::metadata(se)[["aesthetics"]][["color"]]

  if (is.null(color_map)) {
    return(NULL)
  }

  annotations_names <- purrr::set_names(names(color_map))
  annotations <- purrr::map(annotations_names, ~ se[[.x]])

  heatmap_annotation_args <- c(
    annotations,
    list(
      col = color_map,
      which = which,
      show_legend = legends_names,
      show_annotation_name = legends_names,
      gp = grid::gpar(col = "white")
    )
  )

  heatmap_annotation <- do.call(
    ComplexHeatmap::HeatmapAnnotation,
    heatmap_annotation_args
  )

  return(heatmap_annotation)
}


#' Remove diagonal from [ComplexHeatmap::Heatmap]
#' @export
#' @noRd
remove_diagonal <- function(j, i, x, y, w, h, fill) {
  if(i == j) {
    grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = "white", col = "white"))
  }
}
