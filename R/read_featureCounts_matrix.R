#' Read counts matrix from featureCounts file(s)
#'
#' @param paths Paths to featureCounts files
#' @param pattern Pattern to look for and remove from sample names.
#'   All instances of this pattern will be removed from sample names.
#'
#'   The default interpretation is a regular expression, as described in stringi::about_search_regex.
#'
#' @return A matrix of counts from featureCounts.
#'   Columns correspond to samples with sample names as column names, and
#'   rows correspond to features with feature names (content of `Geneid` column) as row names
#' @export
read_featureCounts_matrix <- function(paths, pattern = "\\.bam") {
  read_join_featureCounts(paths) %>%
    extract_counts_from_featureCounts() %>%
    clean_matrix_colnames(pattern = pattern)
}

#' Return counts matrix from a featureCounts data frame
#'
#' @param featureCounts A featureCounts data frame.
#'
#' @return A matrix of counts from featureCounts.
#'   Columns correspond to samples with sample names as column names, and
#'   rows correspond to features with feature names (content of `Geneid` column) as row names.
#' @noRd
extract_counts_from_featureCounts <- function(featureCounts) {
  featureCounts %>%
    dplyr::select(
      -c('Chr', 'Start', 'End', 'Strand', 'Length')
    ) %>%
    tibble::column_to_rownames(var = 'Geneid') %>%
    as.matrix()
}

#' Read and join featureCounts counts file(s)
#'
#' @param paths Path(s) to featureCounts counts file(s).
#'
#' @return A data frame with columns corresponding to the output of featureCounts -
#'   `Geneid, Chr, Start, End, Strand, Length` - and counts for all samples.
#' @noRd
#'
#' @examples
#' read_join_featureCounts("data/counts/sample_1.txt")
#' read_join_featureCounts(c("data/counts/sample_1.txt", "data/counts/sample_2.txt"))
read_join_featureCounts <- function(paths) {
  join_featureCounts(
    lapply(paths, read_featureCounts)
  )
}

#' Read a featureCounts counts file
#'
#' @param path Path to featureCounts counts file.
#'
#' @return A data frame with columns corresponding to the output of featureCounts -
#'   `Geneid, Chr, Start, End, Strand, Length, <sample_counts>`
#' @noRd
#'
#' @examples
#' read_feaureCounts("data/counts/sample.txt")
read_featureCounts <- function(path) {
  utils::read.table(
    path,
    header = TRUE,
    skip = 1
  )
}

#' Join featureCounts counts data frames into a single data frame with all sample info
#'
#' @param featureCounts A list of data frames with columns corresponding to the output of featureCounts.
#'   Such data frames can be produced by [read_featureCounts()].
#'
#' @return A data frame with all features and counts from all samples provided.
#' @noRd
join_featureCounts <- function(featureCounts) {
  purrr::reduce(
    featureCounts,
    dplyr::full_join,
    by = c('Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length')
  )
}
