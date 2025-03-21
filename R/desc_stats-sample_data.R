#' Display table of sample data
#'
#' @inheritParams tbl_library_sizes
#'
#' @return A table of sample data. Returns an interactive data table with
#'   [DT::datatable()] or a `data.frame` from [tibble::tbl_df].
#' @export
tbl_sample_data <- function(se, interactive = TRUE) {
  sample_data <- get_samples_data(se)

  if (interactive == TRUE) {
    sample_data <- DT::datatable(sample_data, rownames = FALSE)
  }

  return(sample_data)
}


#' Display a frequency table of number of samples per combination of interest
#' variables in a `SummarizedExperiment`
#'
#' @inheritParams write_se
#' @param table_format A character vector of the frequency table format to
#'   return
#'
#' @return A frequency table of number of samples per combination of interest
#'   variables
#' @export
tbl_interest_variables_freq <- function(se, table_format = c("datatable", "data.frame")) {

  table_format <- rlang::arg_match(table_format)

  # get sample data for interest variables
  interest_variables <- get_modeling_variable(se, "interest")
  interest_variables_sample_data <- get_samples_data(se)[, interest_variables, drop = FALSE]

  # construct frequency table of interest variables sample data
  freq_table <- as.data.frame(table(interest_variables_sample_data))

  # coerce the frequency table to the requested format
  freq_table <- switch(table_format,
    datatable = DT::datatable(freq_table, rownames = FALSE),
    data.frame = as.data.frame(freq_table)
  )

  return(freq_table)
}


#' Plot pairs scatter plot matrix for all modeling variables of
#' `SummarizedExperiment`
#'
#' @inheritParams write_se
#' @param seed A numeric seed to control the position of the point jittering.
#'   Used for reproducible results.
#'
#' @return A scatter plot matrix comparing pairs of all sample variables set in
#'   model variables.
#' @export
plot_modeling_variables_pairs <- function(se, seed = 2025) {

  # get sample data for interest variables
  all_modeling_variables <- get_modeling_all_variables(se)
  sample_data <- get_samples_data(se)

  # throw error if no modeling variables
  if (is.null(all_modeling_variables)) {
    abort_no_modeling_variables("se")
  }

  # make pairs plot
  modeling_variables_pairs_plot <- GGally::ggpairs(
    sample_data,
    columns = all_modeling_variables,
    upper = list(
      combo = GGally::wrap("dot", position = ggplot2::position_jitter(seed = seed))
    ),
    lower = list(
      combo = "facetdensitystrip"
    ),
    diag = list(
      continuous = GGally::wrap("densityDiag", alpha = 0.5)
    )
  )

  return(modeling_variables_pairs_plot)
}


#' Plot canonical correlation heatmap for sample data of specific variable
#' classes from `SummarizedExperiment`
#'
#' @inheritParams write_se
#' @param modeling_variable_class A character for the class of modeling
#'   variables to plot canonical correlation for.
#'
#' @return A [ComplexHeatmap::Heatmap] of canonical correlation values for the
#'   sample data variables for a specified class of modeling variables.
#' @export
plot_canonical_correlation <- function(se, modeling_variable_class = c("all", "explanatory")) {

  modeling_variable_class <- rlang::arg_match(modeling_variable_class)

  # get sample data for interest variables
  modeling_variables <- switch(modeling_variable_class,
    all = get_modeling_all_variables(se),
    explanatory = get_modeling_explanatory_variables(se)
  )
  sample_data <- get_samples_data(se)

  # throw error if no modeling variables
  if (is.null(modeling_variables)) {
    abort_no_modeling_variables("se")
  }

  # compute canonical correlation
  can_cor_formula <- stats::reformulate(modeling_variables)
  can_cor_matrix <- variancePartition::canCorPairs(
    formula = can_cor_formula,
    data = sample_data
  )
  # plot canoncical correlation heatmap
  can_cor_heatmap <- ComplexHeatmap::Heatmap(
    matrix = can_cor_matrix,
    col = circlize::colorRamp2(seq(0, 1, by = 0.1), scales::viridis_pal(option = "viridis")(11)),
    # col = scales::viridis_pal(option = "viridis"),
    name = "Canonical\ncor.",
    clustering_distance_columns = function(m) stats::as.dist(1 - m),
    clustering_distance_row = function(m) stats::as.dist(1 - m)
  )

  return(can_cor_heatmap)
}
