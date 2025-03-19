#' Add a library size variable to a `SummarizedExperiment`
#'
#' Library size is computed as the column sums of the counts matrix.
#'
#' @inheritParams construct_se
#' @param library_size_var A character vector name for the library size
#'   variable.
#'
#' @return The `SummarizedExperiment` object passed as input with an additional
#'   variable in the sample data for library size.
#' @export
add_library_size <- function(se, library_size_var = "library_size") {
  se[[library_size_var]] <-  colSums(SummarizedExperiment::assay(se, "counts"))

  return(se)
}


#' Make table of library sizes from a matrix of counts
#'
#' @inheritParams construct_se
#' @param interactive A logical controlling whether to return an interactive
#'   data table with [DT::datatable()] or a `data.frame` from [tibble::tbl_df].
#'
#' @return A table of library sizes (matrix column sums) for each sample ordered
#'   by increasing library size. Returns an interactive data table with
#'   [DT::datatable()] or a `data.frame` from [tibble::tbl_df].
#' @export
tbl_library_sizes <- function(counts, interactive = TRUE) {
  library_sizes <- counts |>
    colSums() |>
    tibble::enframe(name = "sample", value = "library_size") |>
    dplyr::arrange(library_size) |>
    dplyr::mutate(
      library_size = prettyNum(library_size, big.mark = ",", scientific = FALSE)
    )
  colnames(library_sizes) <- c("Sample", "Library size")

  if (interactive == TRUE) {
    library_sizes <- DT::datatable(library_sizes, rownames = FALSE)
  }

  return(library_sizes)
}


#' Plot the library size for each sample
#'
#' @inheritParams construct_se
#' @param library_size_var A character of the variable in the
#'   `SummarizedExperiment` sample data to plot. This variable must exist in the
#'   `SummarizedExperiment` sample data.
#' @param interactive A logical controlling whether to return an interactive
#'   plotly plot or a static ggplot.
#'
#' @return A dot plot of library size for each sample
#' @export
plot_library_sizes <- function(se, library_size_var, interactive = TRUE) {
  # check that SummarizedExperiment contains the library sizes variable
  if (is.null(se[[library_size_var]])) {
    abort_no_sample_variable(var = library_size_var, se = "se")
  }
  # get samples data
  samples_data <- get_samples_data(se)

  # plot library sizes
  library_sizes_plot <- ggplot2::ggplot(
    samples_data,
    ggplot2::aes(x = .data[["sample"]], y = .data[[library_size_var]], label = .data[["sample"]])
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6)) +
    ggplot2::labs(y = "Library size (Millions)")

  if (interactive == TRUE) {
    library_sizes_plot <- plotly::ggplotly(library_sizes_plot)
  }

  return(library_sizes_plot)
}


#' Plot the log CPM density of each sample
#'
#' @importFrom rlang .data
#'
#' @inheritParams write_se
#' @inheritParams plot_library_sizes
#'
#' @return A density plot of the log CPM each sample
#' @export
plot_logcpm_density <- function(se, interactive = TRUE) {
  # compute log cpm and make long dataframe
  logcpm_long <- se |>
    edgeR::cpm(log = TRUE) |>
    as.data.frame() |>
    tidyr::pivot_longer(
      cols = tidyselect::everything(),
      names_to = "sample",
      values_to = "logcpm"
    )

  # plot density of logcpm for each library
  logcpm_density_plot <- logcpm_long |>
    ggplot2::ggplot(ggplot2::aes(x = .data[["logcpm"]], color = .data[["sample"]])) +
    ggplot2::geom_density()

  if (interactive == TRUE) {
    logcpm_density_plot <- plotly::ggplotly(logcpm_density_plot)
  }

  return(logcpm_density_plot)
}
