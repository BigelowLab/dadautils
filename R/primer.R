#' Count the number of matches of the primer in the subject
#' 
#' @export
#' @param primer character, pattern string
#' @param fn XString subject to match
#' @return numeric, count of hits
primer_hits <- function(primer, fn) {
  nhits <- Biostrings::vcountPattern(primer, 
                                     ShortRead::sread(ShortRead::readFastq(fn)), 
                                     fixed = FALSE)
  return(sum(nhits > 0))
}


#' Compute primer counts
#'
#' @export
#' @param FWD character, forward primer
#' @param REV character, reverse primer
#' @param fn_FWD XString forward subject to match
#' @param fn_REV XString, reverse subject to match
#' @param form character, desired output format - either 'matrix' with row names or "table" (tibble)
#' @return table or matrix of counts
primer_counts <- function(FWD, REV, fn_FWD, fn_REV,
  form = c("matrix", "table")[2]){
  
  x <- rbind(
    FWD.ForwardReads = sapply(FWD, primer_hits, fn = fn_FWD[[1]]),
    FWD.ReverseReads = sapply(FWD, primer_hits, fn = fn_REV[[1]]),
    REV.ForwardReads = sapply(REV, primer_hits, fn = fn_FWD[[1]]),
    REV.ReverseReads = sapply(REV, primer_hits, fn = fn_REV[[1]]))
  
  if (tolower(form[1]) == 'table'){
    x <- dplyr::as_tibble(x, rownames = "name")
  }
  
  x
}


#' Given a primer, compute all of the possible orientations
#' 
#' @export
#' @param primer input sequence, see \code{\link[Biostrings]{DNAString}}
#' @return name character vector
all_orients <- function(primer = "TTGAAAA-CTC-N")  {
  dna <- Biostrings::DNAString(primer)
  orients <- c(Forward = dna, 
             Complement =  Biostrings::complement(dna), 
             Reverse = IRanges::reverse(dna), 
             RevComp = Biostrings::reverseComplement(dna)) 
  return(sapply(orients, toString))
} 

#' Remove primers
#'
#' @export
#' @param filelist list of forward and reverse fastq files
#' @param output_path character, the output path
#' @param compress logical, see \code{\link[dada2]{removePrimers}}
#' @param ... other arguments for \code{\link[dada2]{removePrimers}}
#' @param save_results logical, save CSV if TRUE
#' @return two element list of tibbles see \code{\link[dada2]{removePrimers}}
remove_primers <- function(filelist,
                           output_path = file.path(dirname(filelist$forward[1]),'no_primers'),
                           save_results = TRUE,
                           compress = TRUE,
                           ...){
  
  fprim <- file.path(output_path, basename(filelist$forward))
  
  norev <- length(filelist$reverse) == 0
  
  rprim <- if(norev){
    NULL
  } else {
    file.path(output_path, basename(filelist$reverse))
  }
  
  if (compress){
    fprim <- charlier::add_extension(fprim, ext = ".gz", no_dup = TRUE)
    rprim <- if(!norev) charlier::add_extension(rprim, ext = ".gz", no_dup = TRUE)
  }
  
  x <- list(
    forward = dada2::removePrimers(filelist$forward,
                                   fprim,
                                   compress = compress,
                                   ...) %>% 
              dplyr::as_tibble(rownames = "name"),
    reverse = if (norev) {
              NULL
              } else {
              dada2::removePrimers(filelist$reverse,
                           rprim,
                           compress = compress,
                           ...) %>% 
              dplyr::as_tibble(rownames = "name")
              })
    
  
  if (save_results) {
    readr::write_csv(x$forward, file.path(output_path, "primer-results_forward.csv"))
    if (!norev) readr::write_csv(x$reverse, file.path(output_path, "primer-results_reverse.csv"))
  }
  x
}
