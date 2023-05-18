#' Plot results summary data
#'
#' @param res A `DESeqResults` object.
#'
#' @return A `ggplot` object.
#' @export
plot_results_summary <- function(res) {
  res %>%
    extract_results_summary_data() %>%
    as.data.frame() %>%
    dplyr::mutate(
      percent_label = percent_not_all_zeros %>%
        round(2) %>%
        paste0("%"),
      degs = dplyr::case_when(
        category %in% c("up", "down") ~ TRUE,
        category %in% c("outliers", "low counts", "all zeros") ~ FALSE
      )
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = category, y = count)) +
    ggplot2::geom_col() +
    ggplot2::geom_text(ggplot2::aes(label = count), vjust = -0.2) +
    ggplot2::geom_text(ggplot2::aes(label = percent_label), vjust = 1.5, color = "white") +
    ggplot2::facet_grid(~degs, scales = "free_x") +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank()
    )
}

#' Extract summary data for `DESeqResults`
#'
#' @description This is essentially a rewrite of `DESeq2::summary()` that returns
#'   a useful `S4Vectors::DF` object instead of just printing summary info.
#'
#' @param object A `DESeqResults` object.
#'
#' @return A `DataFrame` that closely mimics the output of `DESeq2::summary()`.
#' @export
extract_results_summary_data <- function(object) {
  metadata <- extract_results_summary_metadata(object)

  test_column <- set_summary_test_column(object)

  all_zero <- count_all_zeros(object)
  not_all_zero <- nrow(object) - all_zero

  up <- sum(
    object[[test_column$column]] < metadata$alpha & object$log2FoldChange > metadata$lfc_threshold,
    na.rm = TRUE
  )
  down <- sum(
    object[[test_column$column]] < metadata$alpha & object$log2FoldChange < metadata$lfc_threshold,
    na.rm = TRUE
  )
  outlier <- sum(object$baseMean > 0 & is.na(object$pvalue))
  filtered_low <- sum(!is.na(object$pvalue) & is.na(object$padj))

  df <- data.frame(
    category = c("up", "down", "outliers", "low counts", "all zeros"),
    count = c(up, down, outlier, filtered_low, all_zero)
  ) %>%
    dplyr::mutate(
      percent_total = count / nrow(object) * 100,
      percent_not_all_zeros = count / not_all_zero * 100,
      percent_not_all_zeros = replace(percent_not_all_zeros, category == "all zeros", NA)
    )

  results_summary_data <- S4Vectors::DataFrame(
    df
  )
  S4Vectors::metadata(results_summary_data) <- metadata

  return(results_summary_data)
}

#' Extract summary metadata for `DESeqResults`
#'
#' @inheritParams extract_results_summary_data
#'
#' @return A list of metadata.
#' @noRd
extract_results_summary_metadata <- function(object) {
  metadata <- list(
    alpha = set_summary_alpha(object),
    lfc_threshold = set_summary_lfc_threshold(object),
    filter_threshold = set_summary_filter_threshold(object)
  )
}

#' Set the summary test column
#'
#' @param object
#'
#' @return A list of summary metadata.
#' @noRd
set_summary_test_column <- function(object) {
  if(has_sval(object)) {
    list("column" = "svalue", "name" = "s-value")
  } else {
    list("column" = "padj", "name" = "adjusted p-value")
  }
}

#' Count the number of genes with 0 counts from `DESeqResults`
#'
#' @inheritParams extract_results_summary_data
#'
#' @return A count of the number of genes with 0 counts.
#' @noRd
count_all_zeros <- function(object) {
  sum(object$baseMean == 0)
}

#' Set the summary alpha level from `DESeqResults`
#'
#' @inheritParams extract_results_summary_data
#'
#' @return alpha level.
#' @noRd
set_summary_alpha <- function(object) {
  if (has_sval(object)) {
    alpha <- 0.005
  } else {
    if (is.numeric(S4Vectors::metadata(object)$alpha)) {
      alpha <- S4Vectors::metadata(object)$alpha
    } else {
      alpha <- 0.1
    }
  }

  return(alpha)
}

#' Set the summary log 2 fold-change threshold from `DESeqResults`
#'
#' @inheritParams extract_results_summary_data
#'
#' @return log2 fold change threshold.
#' @noRd
set_summary_lfc_threshold <- function(object) {
  ifelse(
    is.numeric(S4Vectors::metadata(object)$lfcThreshold),
    S4Vectors::metadata(object)$lfcThreshold,
    0
  )
}

#' Set the summary filtering threshold threshold from `DESeqResults`
#'
#' @inheritParams extract_results_summary_data
#'
#' @return Filtering threshold.
#' @noRd
set_summary_filter_threshold <- function(object) {
  if (!has_sval(object)) {
    if (is.numeric(S4Vectors::metadata(object)$filterThreshold)) {
      filter_threshold <- round(S4Vectors::metadata(object)$filterThreshold[[1]], 2)
    } else {
      filter_threshold <- 0
    }
  }

  return(filter_threshold)
}

#' Check if a `DESeqResults` object has an `sval` column
#'
#' @inheritParams extract_results_summary_data
#'
#' @return boolean for whether `object` has an `sval` column.
#' @noRd
has_sval <- function(object) {
  "svalue" %in% names(object)
}

#' Wrap together [DESeq2::results()] and [DESeq2::lfcShrink()] for convenience
#'
#' @description A common workflow in `DESeq2` involves extracting a `DESeqResults`
#'   object with [DESeq2::results()] followed by shrinkage of log2 fold changes
#'   with [DESeq2::lfcShrink()]. These are kept as separate functions within
#'   `DESeq2` so that additional parameters can be passed to the functions for
#'   shrinkage estimation. However, I do not often do this, so I wrap these two
#'   functions together for convenience here.
#'
#' @param dds A `DESeqDataSet` object
#' @param name The name of the coefficient/comparison to extract results for
#'   (see [DESeq2::resultsNames()]). This is passed to the `name` argument of
#'   [DESeq2::results()] and the `coef` argument of [DESeq2::lfcShrink()].
#' @inheritParams DESeq2::results
#' @inheritParams DESeq2::lfcShrink
#'
#' @return A `DESeqResults` object with shrunken log2 fold changes.
#' @export
extract_shrunken_results <- function(
    dds,
    name,
    alpha = 0.05,
    lfcThreshold = 0,
    type = "apeglm",
    quiet = TRUE) {
  res <- DESeq2::results(
    dds,
    name = name,
    alpha = alpha,
    lfcThreshold = lfcThreshold
  )

  DESeq2::lfcShrink(
    dds,
    coef = name,
    res = res,
    type = type,
    lfcThreshold = lfcThreshold,
    quiet = quiet
  )
}
