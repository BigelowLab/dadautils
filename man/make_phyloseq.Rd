% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phyloseq.R
\name{make_phyloseq}
\alias{make_phyloseq}
\title{A simple wrapper around \code{\link[phyloseq]{phyloseq}}}
\usage{
make_phyloseq(
  otu_table = NULL,
  sample_data = NULL,
  tax_table = NULL,
  merge_refseq = TRUE,
  refseq_fun = Biostrings::DNAStringSet
)
}
\arguments{
\item{otu_table}{table of otu sequence counts}

\item{sample_data}{table of ancillary sample info}

\item{tax_table}{table of taxonomy}

\item{merge_refseq}{logical, if TRUE that merge a refseq object}

\item{refseq_fun}{the name of the function to call, by default \code{\link[Biostrings]{DNAStringSet}}.
Must accept the output of \code{\link[phyloseq]{taxa_names}}}
}
\value{
a \code{\link[phyloseq]{phyloseq}} object
}
\description{
A simple wrapper around \code{\link[phyloseq]{phyloseq}}
}
