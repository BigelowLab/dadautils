library(rlang)
library(dplyr)
library(ggplot2)
library(charlier)
#library(dadautils)
devtools::load_all()
library(ShortRead)

PATH <- "/mnt/storage/data/edna/dada/projects/ben_foo/weeds-workshop"
files <- list_filepairs(file.path(PATH, "cutadapt"))


    


#' Build a table of expected error info for paired files
#'
#' @export
#' @param filelist a list of \code{forward} and \code{reverse} filename pairs
#' @return a table of 
#' \itemize{
# \item{name_f, name_r filenames}
# \item{data_f, data_r list column of fastq data as \code{ShortReadQ}}
# \item{score_f, score_r list column of vectors returned by \code{ShortRead::alphabetScore}}
# \item{ee_f, ee_r list column of expected error}  
#'}
list_maxEE <- function(filelist = example_filepairs()){
  fqlist <- sapply(filelist,
      function(x){
        lapply(x, ShortRead::readFastq)
      }, simplify = FALSE) 
  scorelist <- sapply(fq,
    function(x){
      sapply(x, ShortRead::alphabetScore)
    }, simplify = FALSE)
  
}