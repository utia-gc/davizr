#' Construct a DESeqDataSet from a counts matrix and column data
#'
#' @description `construct_matrix_dds()` is essentially a wrapper around [DESeq2::DESeqDataSetFromMatrix()]
#'   that offers additional validation of sample order in the counts matrix and column data as required to
#'   construct the `DESeqDataSet` object.
#'
#'   `construct_matrix_dds()` subsets the counts matrix and column data by the sorted intersection of counts matrix
#'   column names and column data row names, then constructs a `DESeqDataSet` object from the subset of the counts
#'   matrix and column data.
#'
#' @inheritParams DESeq2::DESeqDataSetFromMatrix
#'
#' @return A `DESeqDataSet` object.
#' @export
construct_matrix_dds <- function(countData, colData, design = ~ 1, ...) {
  samples_order <- stringr::str_sort(BiocGenerics::intersect(colnames(countData), rownames(colData)), numeric = TRUE)

  if (length(samples_order) == 0) {
    abort_no_common_samples()
  }

  counts <- countData[, samples_order]
  col_data <- colData[samples_order, ]

  DESeq2::DESeqDataSetFromMatrix(
    countData = counts, colData = col_data, design = stats::formula(design), ...
  )
}

#' Abort when `countData` and `colData` have no samples in common.
#'
#' @return A condition of class `no_common_samples`
#' @noRd
abort_no_common_samples <- function() {
  cli::cli_abort(
    c("{.var countData} and {.var colData} must have samples in common.",
      "x" = "There are 0 shared sample names between {.var countData} and {.var colData}.",
      "i" = "Did you verify {.var colnames(countData)} and {.var rownames(colData)} have elements in common?"),
    "no_common_samples"
  )
}
