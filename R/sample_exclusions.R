# sample_exclusions.R
# Exclude samples from an analysis

#' Drop excluded samples from a `SummarizedExperiment`
#'
#' @inheritParams plot_sample_exclusions_heatmap
#'
#' @returns A `SummarizedExperiment` object without the columns that overlap
#'   sample names found in `sample_exclusions`
#' @export
drop_excluded_samples <- function(se, sample_exclusions) {
  # get vector of sample names to exclude by pulling names out of sample exclusions
  samples_exclude <- sample_exclusions |>
    unlist() |>
    unname() |>
    unique()

  # throw error if there is a sample to exclude which isn't a sample in the SummarizedExperiment
  nonexistent_samples_exclude <- setdiff(samples_exclude, colnames(se))
  if (length(nonexistent_samples_exclude) > 0L) abort_no_sample(nonexistent_samples_exclude, "se")

  # get the samples to include
  samples_include <- setdiff(colnames(se), samples_exclude)
  # filter SummarizedExperiment to only the samples to include
  se <- se[, samples_include]

  return(se)
}


#' Plot a heatmap of sample exclusions and the reasons for exclusions
#'
#' @param se A `SummarizedExperiment` object
#' @param sample_exclusions A `list` of key:value pairs mapping reason for
#'   exclusion to samples excluded.
#'
#' @returns A `plotly` object heatmap of samples to be excluded and reasons for
#'   exclusion
#' @export
plot_sample_exclusions_heatmap <- function(se, sample_exclusions) {
  # construct matrix where each column is a reason for exclusion and each row is a sample
  # element[i, j] indicates whether or not sample i is excluded for reason j
  binary_exclusions_all <- lapply(sample_exclusions, \(sample_exclusion) colnames(se) %in% sample_exclusion) |>
    list2DF() |>
    as.matrix()
  rownames(binary_exclusions_all) <- colnames(se)

  # filter exclusions matrix to only those samples that are excluded for at least one reason
  binary_exclusions <- binary_exclusions_all[rowSums(binary_exclusions_all) > 0, , drop = FALSE]
  # convert to numeric matrix for plotting
  binary_exclusions[] <- apply(binary_exclusions, 2, as.numeric)

  # plot the binary heatmap for excluded samples
  exclusions_heatmap <- plotly::plot_ly(
    x = colnames(binary_exclusions),
    y = rownames(binary_exclusions),
    xgap = 2,
    ygap = 2,
    z = binary_exclusions,
    type = "heatmap",
    colorscale = data.frame(
      z = c(0, 0.5, 0.5, 1),
      col = c("#f7f7f7", "#f7f7f7", "#CF142B", "#CF142B")
    ),
    showscale = TRUE,
    hovertemplate = paste(
      "<b>Sample:</b> %{y}<br>",
      "<b>Reason:</b> %{x}<br>",
      "<b>Excluded:</b> %{z}<br>",
      "<extra></extra>"
    ),
    customdata = binary_exclusions,
    colorbar = list(
      title = "Excluded",
      tickmode = "array",
      tickvals = c(0, 1),
      ticktext = c("FALSE", "TRUE")
    )
  ) |>
    plotly::layout(
      xaxis = list(title = "Reason"),
      yaxis = list(title = "Sample")
    )

  return(exclusions_heatmap)
}

