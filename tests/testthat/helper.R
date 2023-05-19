example_dds <- function(...) {
  withr::with_seed(
    2023,
    DESeq2::makeExampleDESeqDataSet(betaSD = 0.75, ...) %>%
      DESeq2::DESeq(quiet = TRUE)
  )
}

counts_matrix <- function() {
  matrix(
    rep(0, 16), ncol = 2,
    dimnames = list(
      c("YAL069W", "YAL068W-A", "YAL068C", "YAL067W-A",
        "YAL067C", "YAL066W", "YAL065C", "YAL064W-B"),
      c("sample_1_sorted", "sample_2_sorted")
    )
  )
}

random_matrix <- function(n = 100, m = 8) {
  withr::with_seed(
    0,
    matrix(
      sample.int(100, size = m * n, replace = TRUE),
      nrow = n, ncol = m,
      dimnames = list(paste0("gene", 1:n), paste0("sample", 1:m))
    )
  )
}

random_col_data <- function(m = 8) {
  data.frame(
    row.names = paste0("sample", 1:m),
    condition = factor(c(rep("control", m/2), rep("treat", m/2))),
    batch = factor(rep(c(1, 2), m/2))
  )
}

random_dds <- function() {
  DESeq2::DESeqDataSetFromMatrix(
    countData = random_matrix(),
    colData = random_col_data(),
    design = formula("~ batch + condition")
  )
}

library_detected_genes <- function() {
  x <- c(97, 95, 93, 95, 91, 89, 95, 95)
  names(x) <- paste0("sample", 1:8)

  return(x)
}

library_detected_genes_low_thresh <- function() {
  x <- c(98, 98, 99, 99, 98, 100, 98, 99)
  names(x) <- paste0("sample", 1:8)

  return(x)
}

results_summary_data <- function() {
  metadata <- list(
    alpha = 0.05,
    lfc_threshold = 0,
    filter_threshold = 3.37
  )
  df <- S4Vectors::DataFrame(
    category = c("up", "down", "outliers", "low counts", "all zeros"),
    count = c(66, 62, 2, 154, 0),
    percent_total = c(6.6, 6.2, 0.2, 15.4, 0),
    percent_not_all_zeros = c(6.6, 6.2, 0.2, 15.4, NA)
  )
  S4Vectors::metadata(df) <- metadata

  df
}

vsd_pca_data <- function() {
  object <- DESeq2::varianceStabilizingTransformation(example_dds())

  ntop <- 500

  row_vars <- MatrixGenerics::rowVars(SummarizedExperiment::assay(object))
  select <- order(row_vars, decreasing = TRUE)[seq_len(min(ntop, length(row_vars)))]

  pca <- stats::prcomp(t(SummarizedExperiment::assay(object)[select, ]))

  x <- pca$x %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "sample_name")

  col_data <- object %>%
    SummarizedExperiment::colData() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "sample_name")

  pca_df <- dplyr::full_join(
    col_data,
    x,
    by = "sample_name"
  ) %>%
    tibble::column_to_rownames(var = "sample_name") %>%
    S4Vectors::DataFrame()
  S4Vectors::metadata(pca_df)$prcomp <- pca
  S4Vectors::metadata(pca_df)$variance_explained <- data.frame(
    principal_component = seq_len(length(pca$sdev)),
    proportion_var = pca$sdev^2 / sum(pca$sdev^2),
    cumulative_proportion_var = cumsum(pca$sdev^2 / sum(pca$sdev^2))
  )

  return(pca_df)
}

vsd_pca_data_100 <- function() {
  object <- DESeq2::varianceStabilizingTransformation(example_dds())

  ntop <- 100

  row_vars <- MatrixGenerics::rowVars(SummarizedExperiment::assay(object))
  select <- order(row_vars, decreasing = TRUE)[seq_len(min(ntop, length(row_vars)))]

  pca <- stats::prcomp(t(SummarizedExperiment::assay(object)[select, ]))

  x <- pca$x %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "sample_name")

  col_data <- object %>%
    SummarizedExperiment::colData() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "sample_name")

  pca_df <- dplyr::full_join(
    col_data,
    x,
    by = "sample_name"
  ) %>%
    tibble::column_to_rownames(var = "sample_name") %>%
    S4Vectors::DataFrame()
  S4Vectors::metadata(pca_df)$prcomp <- pca
  S4Vectors::metadata(pca_df)$variance_explained <- data.frame(
    principal_component = seq_len(length(pca$sdev)),
    proportion_var = pca$sdev^2 / sum(pca$sdev^2),
    cumulative_proportion_var = cumsum(pca$sdev^2 / sum(pca$sdev^2))
  )

  return(pca_df)
}

rld_pca_data <- function() {
  object <- DESeq2::rlog(example_dds())

  ntop <- 500

  row_vars <- MatrixGenerics::rowVars(SummarizedExperiment::assay(object))
  select <- order(row_vars, decreasing = TRUE)[seq_len(min(ntop, length(row_vars)))]

  pca <- stats::prcomp(t(SummarizedExperiment::assay(object)[select, ]))

  x <- pca$x %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "sample_name")

  col_data <- object %>%
    SummarizedExperiment::colData() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "sample_name")

  pca_df <- dplyr::full_join(
    col_data,
    x,
    by = "sample_name"
  ) %>%
    tibble::column_to_rownames(var = "sample_name") %>%
    S4Vectors::DataFrame()
  S4Vectors::metadata(pca_df)$prcomp <- pca
  S4Vectors::metadata(pca_df)$variance_explained <- data.frame(
    principal_component = seq_len(length(pca$sdev)),
    proportion_var = pca$sdev^2 / sum(pca$sdev^2),
    cumulative_proportion_var = cumsum(pca$sdev^2 / sum(pca$sdev^2))
  )

  return(pca_df)
}
