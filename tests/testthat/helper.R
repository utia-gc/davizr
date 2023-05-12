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
