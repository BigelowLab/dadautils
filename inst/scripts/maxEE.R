library(rlang)
library(dplyr)
library(ggplot2)
library(charlier)
#library(dadautils)
devtools::load_all("/mnt/storage/data/edna/packages/dadautils")

library(ShortRead)
library(Rsubread)

#PATH <- "/mnt/storage/data/edna/dada/projects/ben_foo/weeds-workshop"
#files <- list_filepairs(file.path(PATH, "cutadapt"))

PATH <- "/mnt/storage/data/edna/dada/projects/robin_foo/cyanos/cutadapt/lowquality"
files <- list_filepairs(PATH)
ee = paired_quality_scores(files)

pp = paired_ee_per_read(ee)

tt =  paired_ee_threshold()
  

  
quality_scores <- function(filelist = example_filepairs(), nreads = -1, ...) {
  sapply(filelist,
    function(files){
      lapply(files, Rsubread::qualityScores, nreads = nreads, ...)
    }, simplify = FALSE)
}


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