% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/library_size.R
\name{add_library_size}
\alias{add_library_size}
\title{Add a library size variable to a \code{SummarizedExperiment}}
\usage{
add_library_size(se, library_size_var = "library_size")
}
\arguments{
\item{library_size_var}{A character vector name for the library size
variable.}
}
\value{
The \code{SummarizedExperiment} object passed as input with an additional
variable in the sample data for library size.
}
\description{
Library size is computed as the column sums of the counts matrix.
}
