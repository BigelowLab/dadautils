#' Filter and trim
#' 
#' @export
#' @param filelist list of forward and reverse fastq files
#' @param output_path character, the output path
#' @param compress logical, see \code{\link[dada2]{filterAndTrim}}
#' @param multithread numeric, the number of cores to use. Defaults to \code{\link{count_cores}}
#' @param ... other arguments for \code{\link[dada2]{filterAndTrim}}
#' @param save_results logical, save CSV if TRUE
#' @return integer matrix as tibble see \code{\link[dada2]{filterAndTrim}}
filter_and_trim <- function(filelist, 
                            output_path = file.path(dirname(filelist$forward[1]),'filterAndTrim'),
                            save_results = FALSE,
                            compress = TRUE,
                            multithread = count_cores(),
                            ...){
  
  ffilt <- file.path(output_path, basename(filelist$forward))
  rfilt <- file.path(output_path, basename(filelist$reverse))
  
  if (compress){
    ffilt <- add_extension(ffilt, ext = ".gz", no_dup = TRUE)
    rfilt <- add_extension(rfilt, ext = ".gz", no_dup = TRUE)
  }
  
  x <- dada2::filterAndTrim(filelist$forward, 
                            ffilt, 
                            rev = filelist$reverse, 
                            filt.rev = rfilt,
                            compress = compress,
                            multithread = multithread, 
                            ...)

  x <- dplyr::as_tibble(x, rownames = "name")
  if (save_results) {
    x <- readr::write_csv(x, file.path(output_path, "filter_and_trim-results.csv"))
  }   
  x
}


#' Run dada2::learnErrors on a set of fastq files
#'
#' @export
#' @param filelist list of forward and reverse fastq files
#' @param multithread numeric, the number of cores to use. Defaults to \code{\link{count_cores}}
#' @param ... arguments for \code{\link[dada2]{learnErrors}}
#' @param output_path character, the output path
#' @param save_output logical, if TRUE save the output to the specified output_path
#' @param save_graphics logical, if TRUE try to capture quality plots form the resulting cut_files
learn_errors <- function(filelist, 
  multithread = count_cores(),
  output_path = dirname(filelist$forward[1]),
  save_output = FALSE, 
  save_graphics = FALSE,
  ...){

  errs <- list(
      forward =  dada2::learnErrors(filelist$forward, multithread = multithread, ...),
      reverse =  dada2::learnErrors(filelist$reverse, multithread = multithread, ...)
    )

  if (save_output){
    saveRDS(errs, file = file.path(output_path, "learn_errors.rds"))
  }
  
  if (save_graphics){
    pforward <- dada2::plotErrors(errs$forward, nominalQ=TRUE) + 
      ggplot2::ggtitle("Forward")
    preverse <- dada2::plotErrors(errs$reverse, nominalQ=TRUE) +
      ggplot2::ggtitle("Reverse")
    ofile <- file.path(output_path, "learn_errors.pdf")
    grDevices::pdf(ofile, height = 7.5, width = 10.5)    
    try(
      print(pforward + preverse)
    )
    grDevices::dev.off()
  }

  errs
}


#' Run dada
#'
#' @export
#' @param filelist list of forward and reverse fastq files
#' @param errs list of forward and reverse outputs of learnErrors
#' @param multithread numeric, the number of cores to use. Defaults to \code{\link{count_cores}}
#' @param ... arguments for \code{\link[dada2]{dada}}
#' @return list with elements for forward and reverse as returned by \code{\link[dada2]{dada}}
run_dada <- function(filelist, errs, 
  multithread = count_cores(),
  ...){
  

  filelist <- lapply(filelist, dada2::derepFastq)

  x <- list(
    forward = dada2::dada(filelist$forward, errs$forward, multithread = multithread, ...),
    reverse = dada2::dada(filelist$reverse, errs$reverse, multithread = multithread, ...)
  )
  x
}

#' Merge pairs ala dada2::mergePairs
#'
#' @export
#' @param filelist list of forward and reverse fastq files
#' @param dada_r  list of dada2::dada results
#' @param ... arguments for \code{\link[dada2]{mergePairs}}
#' @return as returned by \code{\link[dada2]{mergePairs}}
merge_pairs <- function(filelist, dada_r, ...){
  filelist <- lapply(filelist, dada2::derepFastq)
  x <- dada2::mergePairs(
    dada_r$forward, 
    filelist$forward,
    dada_r$reverse, 
    filelist$reverse, 
    ...) 

  x
}


#' Count uniques
#'
#' @export
#' @param x object from which uniques-vector can be extracted
#' @param ... further arguments for  \code{\link[dada2]{getUniques}}
#' @return integer number of uniques
count_uniques <- function(x, ...){
  sum(dada2::getUniques(x, ...)) 
}



#' Plot quality profiles for one or more FASTQs
#'
#' @export
#' @param x character vector of fastq filenames 
#' @param ofile character, the name of the PDF file to generate of NA
#' @param ... further arguments for \code{\link[dada2]{plotQualityProfile}}
plot_qualityProfile <- function(x, 
  ofile = c(NA, "qualityProfile.pdf")[1],
  ...){
  if (!is.na(ofile)) grDevices::pdf(ofile)
  ok <- dada2::plotQualityProfile(x, ...)
  if (!is.na(ofile)) grDevices::dev.off()
  ok
}



#' Convert a table (like a tibble or data.frame) to a matrix
#'
#' Useful for reconstructing a \code{dada2} matrix from an eDNA table
#' 
#' @export
#' @param x table object like a tibble or data.frame
#' @param rowname the column (variable) to use as source of rownames (before optional transpose)
#' @param transpose logical, if TRUE transpose just before returning.  
#' @return a matrix with column and row names
table_as_matrix <- function(x, 
  rowname = 1,
  transpose = TRUE){
  
  cnames <- colnames(x)
  rnames <- x[[rowname]]
  m <- x %>%
    dplyr::select(-tidyselect::all_of(rowname)) %>%
    as.matrix()
  rownames(m) <- rnames
  if (transpose) m <- t(m)
  return(m)
}



#' Convert an eDNA sequence counts table to a \code{dada2} friendly matrix
#'
#' @export
#' @param x filename or table, if a filename we try to read it as CSV
#' @return \code{dada2} friendly matrix
seqtab_to_matrix <- function(x){
  if (is.character(x) && file.exists(x[1])){ 
    x <- suppressMessages(readr::read_csv(x))
  } 
  table_as_matrix(x)
}

#' Convert an eDNA taxonomy table to a \code{dada2} friendly matrix
#'
#' @export
#' @param x filename or table, if a filename we try to read it as CSV
#' @param ... further arguments for \code{table_to_matrix}
#' @return \code{dada2} friendly matrix
taxtable_to_matrix <- function(x, ...){
   if (is.character(x) && file.exists(x[1])){ 
    x <- suppressMessages(readr::read_csv(x))
  } 
  keep <- c("ASV", "Kingdom", "Supergroup", "Division", "Phylum", "Class", "Order", 
            "Family", "Genus", "Species")
  x %>%
   dplyr::select(tidyselect::any_of(keep)) %>%
   table_as_matrix(transpose = FALSE)
}
