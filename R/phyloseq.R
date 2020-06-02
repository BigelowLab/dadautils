#' A simple wrapper around \code{\link[phyloseq]{phyloseq}}
#'
#' @export
#' @param otu_table table of otu sequence counts
#' @param sample_data table of ancillary sample info 
#' @param tax_table table of taxonomy
#' @param merge_refseq logical, if TRUE that merge a refseq object
#' @param refseq_fun the name of the function to call, by default \code{\link[Biostrings]{DNAStringSet}}.
#'   Must accept the output of \code{\link[phyloseq]{taxa_names}}
#' @return a \code{\link[phyloseq]{phyloseq}} object
make_phyloseq <- function(otu_table = NULL, sample_data = NULL, tax_table = NULL,
  merge_refseq = TRUE, refseq_fun = Biostrings::DNAStringSet){


  order_magnitude <- function(x) {
    pmax(1,ceiling(log10(x)))
  }

  x <- phyloseq::phyloseq(
    phyloseq::otu_table(otu_table, taxa_are_rows=FALSE, errorIfNULL = FALSE), 
    phyloseq::sample_data(sample_data, errorIfNULL = FALSE))
    
  if (!is.null(tax_table)) phyloseq::tax_table(x) <- tax_table
    
  if (merge_refseq){
    ref <- refseq_fun(phyloseq::taxa_names(x))
    names(ref) <- phyloseq::taxa_names(x)
    x <- phyloseq::merge_phyloseq(x, ref)
    ix <- seq_len(phyloseq::ntaxa(x))
    n <- order_magnitude(length(ix))
    phyloseq::taxa_names(x) <- sprintf(paste0("ASV%0.",n, "i"), ix)
  }
    
  x

}