% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-se.R
\name{abort_no_modeling_variables}
\alias{abort_no_modeling_variables}
\title{Throw an error when \code{SummarizedExperiment} has no modeling variables}
\usage{
abort_no_modeling_variables(se)
}
\arguments{
\item{se}{The \code{SummarizedExperiment} that does not modeling variables
variable}
}
\value{
An \code{rlang_error} of custom class \code{error_no_modeling_variables}
}
\description{
Throw an error when \code{SummarizedExperiment} has no modeling variables
}
