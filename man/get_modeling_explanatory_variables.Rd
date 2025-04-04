% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-se.R
\name{get_modeling_explanatory_variables}
\alias{get_modeling_explanatory_variables}
\title{Get explanatory modeling variables from a \code{SummarizedExperiment}}
\usage{
get_modeling_explanatory_variables(se)
}
\arguments{
\item{se}{A \code{SummarizedExperiment}}
}
\value{
A character vector of explanatory variables.
}
\description{
Explanatory variables are the union of the interest and nuisance
variables.
}
