featureCounts_matrix <- function() {
  matrix(
    rep(0, 16), ncol = 2,
    dimnames = list(
      c("YAL069W", "YAL068W-A", "YAL068C", "YAL067W-A",
        "YAL067C", "YAL066W", "YAL065C", "YAL064W-B"),
      c("sample_1_sorted", "sample_2_sorted")
    )
  )
}
