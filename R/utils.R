#' Clean the column names of a matrix
#'
#' @description A thin wrapper around `stringr::str_replace_all()` to clean matrix column names.
#'
#' @param m A matrix.
#' @inheritParams stringr::str_replace_all
#'
#' @return Matrix with cleaned column names.
#' @export
clean_matrix_colnames <- function(m, pattern, replacement = "") {
  colnames(m) <- stringr::str_replace_all(
    colnames(m), pattern = pattern, replacement = replacement
  )

  m
}
