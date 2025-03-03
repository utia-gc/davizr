#' Read counts matrix from file
#'
#' @param path A string path to a counts file in CSV format such as the one
#'   created by [utia-gc/rnaseq](https://github.com/utia-gc/rnaseq).
#'
#' @return An integer matrix with genes as rows and samples as columns. Row
#'   names are gene names. Column names are sample names.
#' @export
read_counts_file <- function(path) {
  counts <- data.table::fread(path, header = TRUE, data.table = FALSE) |>
    janitor::clean_names() |>
    tibble::column_to_rownames(var = "geneid") |>
    dplyr::relocate(sort(tidyselect::everything())) |>
    as.matrix()

  return(counts)
}
