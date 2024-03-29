% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quality_profile.R
\name{quality_profile_pairs}
\alias{quality_profile_pairs}
\title{Compute pairwise (forward-reverse) quality profiles}
\usage{
quality_profile_pairs(
  filelist = list(forward = c(system.file("extdata", "sam1F.fastq.gz", package =
    "dada2"), system.file("extdata", "sam2F.fastq.gz", package = "dada2")), reverse =
    c(system.file("extdata", "sam1R.fastq.gz", package = "dada2"), system.file("extdata",
    "sam2R.fastq.gz", package = "dada2"))),
  n = 5e+05,
  aggregate = FALSE,
  amplicon_length = 400,
  min_overlap = 20,
  plot_filename = file.path(dirname(dirname(filelist$forward[1])),
    "quality_profiles.pdf"),
  overlap_filename = file.path(dirname(dirname(filelist$forward[1])), "overlap.csv"),
  ...
)
}
\arguments{
\item{filelist}{list of forward and reverse fastq files}

\item{n}{(Optional). Default 500,000.
The number of records to sample from the fastq file.}

\item{aggregate}{(Optional). Default FALSE.
If TRUE, compute an aggregate quality profile for all fastq files provided.}

\item{amplicon_length}{numeric, the expected amplicon length, used to compute expected overlap}

\item{min_overlap}{numeric, issue warnings for any overlap less than this value}

\item{plot_filename}{character or NULL, if provided save the quality profile plots as PDF
NULL to skip writing the file}

\item{overlap_filename}{character or NULL, if provided save the overlap table to this file as CSV
NULL to skip writing the file}

\item{...}{arguments for \code{\link{quality_profile_cutoff}}}
}
\value{
complex list of forward, reverse and merged quality stats
}
\description{
This a thin wrapper around \code{\link[dada2]{plotQualityProfile}} but produces a list
that includes tables of ancillary statistics in addition to a printable plot object.
}
\seealso{
\code{\link[dada2]{plotQualityProfile}}
}
