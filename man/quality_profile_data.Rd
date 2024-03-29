% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quality_profile.R
\name{quality_profile_data}
\alias{quality_profile_data}
\title{Compute datasets for quality profile}
\usage{
quality_profile_data(
  fl = c(system.file("extdata", "sam1F.fastq.gz", package = "dada2"),
    system.file("extdata", "sam1R.fastq.gz", package = "dada2")),
  n = 5e+05,
  stats = fastq_stats(fl, n = n),
  aggregate = FALSE
)
}
\arguments{
\item{fl}{(Required). \code{character}.
File path(s) to fastq or fastq.gz file(s).}

\item{n}{(Optional). Default 500,000.
The number of records to sample from the fastq file.}

\item{stats}{list of two elements (forward and reverse) with output
of \code{\link{fastq_stats}}}

\item{aggregate}{(Optional). Default \code{FALSE}.
If TRUE, compute an aggregate quality profile for all fastq files provided.}
}
\value{
as list with the following elements
\itemize{
  \item{files character vector of input filenames}
  \item{opts list of input arguments}
  \item{statdf} table of statistical results
  \item{anndf} table of annotation info
  \item{plotdf} table of stats transformed for plotting
  \item{plotdf.summary} table of plot data for plotting (NULL if aggregate is \code{FALSE})
  \item{statdf.summary} table of summary stats for plotting  (NULL if aggregate is \code{FALSE})
  \item{quality_stats} per file quality stats as per \code{\link{fastq_stats}}
}
}
\description{
This a thin wrapper around \code{\link[dada2]{plotQualityProfile}} but produces a list
that includes tables of ancillary statistics.
}
\examples{
\dontrun{
  files <- c(system.file("extdata", "sam1F.fastq.gz", package="dada2"),
             system.file("extdata", "sam1R.fastq.gz", package="dada2"))
  result <- quality_profile_data(files)
  print(result)
 }
}
\seealso{
\code{\link[dada2]{plotQualityProfile}}
}
