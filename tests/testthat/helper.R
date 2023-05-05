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
