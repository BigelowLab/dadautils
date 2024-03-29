% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/primer.R
\name{primer_counts}
\alias{primer_counts}
\title{Compute primer counts}
\usage{
primer_counts(FWD, REV, fn_FWD, fn_REV, form = c("matrix", "table")[2])
}
\arguments{
\item{FWD}{character, forward primer}

\item{REV}{character, reverse primer}

\item{fn_FWD}{XString forward subject to match}

\item{fn_REV}{XString, reverse subject to match}

\item{form}{character, desired output format - either 'matrix' with row names or "table" (tibble)}
}
\value{
table or matrix of counts
}
\description{
Compute primer counts
}
