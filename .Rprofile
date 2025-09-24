if (file.exists("renv/activate.R")) {
  source("renv/activate.R")

  # set up Bioconductor repos
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    renv::install("BiocManager")
  }
  options(repos = BiocManager::repositories())
}

# attach devtools in all interactive sessions
if (interactive()) {
  suppressMessages(require(devtools))
}
