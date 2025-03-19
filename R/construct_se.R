#' Construct a `SummarizedExperiment` from a counts matrix and sample data data
#' frame
#'
#' @description `construct_se()` is essentially a wrapper around
#'   [SummarizedExperiment::SummarizedExperiment()] that offers additional
#'   validation of sample order in the counts matrix and samples data as
#'   required to construct the `SummarizedExperiment` object.
#'
#'   `construct_se()` subsets the counts matrix and column data by the sorted
#'   intersection of counts matrix column names and sample data row names, then
#'   constructs a `SummarizedExperiment` object from the subset of the counts
#'   matrix and sample data.
#'
#'   `construct_se()` by default adds a variable to the sample data for library
#'   size, that is, the column sums of the counts matrix. Controlled by
#'   the `library_size_var` argument.
#'
#' @param counts An integer matrix with genes as rows and samples as columns.
#'   Row names are gene names. Column names are sample names.
#' @param samples A data frame of samples data. Row names are sample names.
#' @param library_size_var A character vector name for the library size
#'   variable. If `NULL`, library size is not computed.
#' @param ... Additional arguments passed to
#'   [SummarizedExperiment::SummarizedExperiment()].
#'
#' @return A `SummarizedExperiment` object. `counts` assay contains counts data.
#' @export
construct_se <- function(counts, samples, library_size_var = "library_size", ...) {
  # get ordered samples in common
  common_samples <- intersect(colnames(counts), rownames(samples))
  samples_order <- stringr::str_sort(common_samples, numeric = TRUE)

  # throw error if there are no samples in common
  if (length(samples_order) == 0) {
    abort_no_common_samples()
  }

  # subset counts matrix and samples data by the ordered samples and construct SummarizedExperiment
  counts_matrix <- counts[, samples_order]
  samples_data <- samples[samples_order, ]

  # construct the SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts_matrix),
    colData = samples_data,
    ...
  )
  # compute library sizes variable
  if (!is.null(library_size_var)) {
    se <- add_library_size(se, library_size_var = library_size_var)
  }

  return(se)
}


#' Abort when `countData` and `colData` have no samples in common.
#'
#' @return A condition of class `no_common_samples`
#' @noRd
abort_no_common_samples <- function() {
  cli::cli_abort(
    c("{.var counts} and {.var samples} must have samples in common.",
      "x" = "There are 0 shared sample names between {.var counts} and {.var samples}.",
      "i" = "Did you verify {.var colnames(counts)} and {.var rownames(samples)} have elements in common?"),
    "error_no_common_samples"
  )
}
