#' Filter and trim
#'
#' @export
#' @param filelist list of forward and reverse fastq files
#' @param output_path character, the output path
#' @param compress logical, see \code{\link[dada2]{filterAndTrim}}
#' @param multithread numeric, the number of cores to use. Defaults to \code{\link{count_cores}}
#' @param truncLen character or numeric.  If character and "auto" then auto-compute trunLen, default "auto"
#' @param cutoff_params list, in the event that trunLen is "auto" then pass these params to 
#'    \code{\link{quality_profile_cutoff}}
#' @param ... other arguments for \code{\link[dada2]{filterAndTrim}}
#' @param save_results logical, save CSV if TRUE
#' @param verbose logical, if TRUE output informational messages
#' @return integer matrix as tibble see \code{\link[dada2]{filterAndTrim}}
filter_and_trim <- function(filelist,
                            output_path = file.path(dirname(filelist$forward[1]),'filterAndTrim'),
                            save_results = FALSE,
                            compress = TRUE,
                            multithread = charlier::count_cores(),
                            truncLen = "auto",
                            cutoff_params = list(score = 30, model = "Mean ~ poly(Cycle, 2)", quantile_min = 0.99),
                            verbose = FALSE, 
                            ...){

  if(FALSE){
    output_path = file.path(dirname(filelist$forward[1]),'filterAndTrim')
    save_results = FALSE
    compress = TRUE
    multithread = charlier::count_cores()
    truncLen = "auto"
    cutoff_params = list(score = 30, model = "Mean ~ poly(Cycle, 2)", quantile_min = 0.99)
    verbose = FALSE
  }

  ffilt <- file.path(output_path, basename(filelist$forward))
  
  norev <- length(filelist$reverse) == 0
  
  rfilt <- if(norev){
      NULL
    } else {
      file.path(output_path, basename(filelist$reverse))
    }
  
  if (compress){
    ffilt <- charlier::add_extension(ffilt, ext = ".gz", no_dup = TRUE)
    rfilt <- if(!norev) charlier::add_extension(rfilt, ext = ".gz", no_dup = TRUE)
  }

  
  if (inherits(truncLen, 'character') && tolower(truncLen[1]) == "auto"){
    # if here, compute cutoffs 'truncLen' per file
    # apply filter and trim per file pair
    
    cutoffs <- sapply(filelist,
      function(files){
        if (length(files) > 0){
          r <- quality_profile_data(files) %>%
              quality_profile_cutoff(params = cutoff_params) %>%
              `[[`("cutoff")
        } else {
          r <- NULL
        }
        r
      }, simplify = FALSE)
    
      x <- lapply( seq_along(filelist[[1]]),
        function(i){
          trunLen <- c(0, 0)
          if (norev){
            trunc_len <- cutoffs$forward$Cycle[i]
          } else {
            trunc_len <- c(cutoffs$forward$Cycle[i], cutoffs$reverse$Cycle[i])
          }
          
          if (verbose){
            cat("filterAndTrim:\n")
            cat(sprintf("  forward file: %s", basename(filelist$forward[i])), "\n")
            if (norev){
              cat("  reverse file: none\n")
            } else {
              cat(sprintf("  reverse file: %s", basename(filelist$reverse[i])), "\n")
            }
            cat(sprintf("  compress: %s", compress), "\n")
            cat(sprintf("  multithread: %s", multithread), "\n")
            cat(sprintf("  truncLen: %s", paste(trunc_len, collapse = " ")), "\n")
            dots <- list(...)
            for (n in names(dots)){
              cat(sprintf("  %s: %s", n, paste(dots[[n]], collapse = " ")), "\n")
            }
          }
          
          dada2::filterAndTrim(filelist$forward[i],
                               ffilt[i],
                               rev = if (norev) { NULL } else {filelist$reverse[i]},
                               filt.rev = rfilt[i],
                               compress = compress,
                               multithread = multithread,
                               truncLen = trunc_len,
                               ...) %>%
              dplyr::as_tibble(rownames = "name")
        }) %>%
        dplyr::bind_rows()
    
  } else {
    
    if (verbose){
      cat("filterAndTrim:\n")
      cat(sprintf("  first forward file: %s", basename(filelist$forward[1])), "\n")
      if (norev){
        cat("  reverse file: none\n")
      } else {
        cat(sprintf("  first reverse file: %s", basename(filelist$reverse[1])), "\n")
      }
      cat(sprintf("  compress: %s", compress), "\n")
      cat(sprintf("  multithread: %s", multithread), "\n")
      cat(sprintf("  truncLen: %s", paste(truncLen, collapse = " ")), "\n")
      dots <- list(...)
      for (n in names(dots)){
        cat(sprintf("  %s: %s", n, paste(dots[[n]], collapse = " ")), "\n")
      }
    }
    
    x <- dada2::filterAndTrim(filelist$forward,
                              ffilt,
                              rev = if (norev) { NULL } else { filelist$reverse},
                              filt.rev = rfilt,
                              compress = compress,
                              multithread = multithread,
                              truncLen = truncLen,
                              ...) %>%
    dplyr::as_tibble(x, rownames = "name")
  }
  
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
#' @return list with elements for forward and reverse as returned by \code{\link[dada2]{learnErrors}}. 
#'    The reverse element may be NULL.
learn_errors <- function(filelist,
  multithread = count_cores(),
  output_path = dirname(filelist$forward[1]),
  save_output = FALSE,
  save_graphics = FALSE,
  ...){

  norev <- length(filelist$reverse) == 0
  
  errs <- list(
      forward =  dada2::learnErrors(filelist$forward, multithread = multithread, ...),
      reverse =  if (norev){
                   NULL
                 } else {
                   dada2::learnErrors(filelist$reverse, multithread = multithread, ...)
                 }
    )

  if (save_output){
    saveRDS(errs, file = file.path(output_path, "learn_errors.rds"))
  }

  if (save_graphics){
    pforward <- dada2::plotErrors(errs$forward, nominalQ=TRUE) +
      ggplot2::ggtitle("Forward")
    if (!norev) preverse <- dada2::plotErrors(errs$reverse, nominalQ=TRUE) +
                   ggplot2::ggtitle("Reverse")
    ofile <- file.path(output_path, "learn_errors.pdf")
    grDevices::pdf(ofile, height = 7.5, width = 10.5)
    try(
      if (!norev) {print(pforward + preverse)} else {print(pforward)}
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
#' @return list with elements for forward and reverse as returned by \code{\link[dada2]{dada}}. 
#'    The reverse element may be NULL
run_dada <- function(filelist, errs,
  multithread = count_cores(),
  ...){

  norev <- length(filelist$reverse) == 0
  if (norev){
    filelist <- dada2::derepFastq(filelist$forward)
    x <- list(
      forward = dada2::dada(filelist, errs$forward, multithread = multithread, ...),
      reverse = NULL
    )
  } else {
    filelist <- lapply(filelist, dada2::derepFastq)
    x <- list(
      forward = dada2::dada(filelist$forward, errs$forward, multithread = multithread, ...),
      reverse = dada2::dada(filelist$reverse, errs$reverse, multithread = multithread, ...)
    )
  }
  x
}

#' Merge pairs ala dada2::mergePairs
#'
#' @export
#' @param filelist list of forward and reverse fastq files
#' @param dada_r  list of dada2::dada results
#' @param save_output logical, if TRUE save the output as RDS
#' @param ... arguments for \code{\link[dada2]{mergePairs}}
#' @return as returned by \code{\link[dada2]{mergePairs}}
merge_pairs <- function(filelist, dada_r,
  save_output = TRUE,
  ...){
  ff <- lapply(filelist, dada2::derepFastq)
  x <- dada2::mergePairs(
    dada_r$forward,
    ff$forward,
    dada_r$reverse,
    ff$reverse,
    ...)
    if (save_output){
      output_path <- dirname(filelist$forward[1])
      saveRDS(x, file = file.path(output_path, "mergers.rds"))
    }
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


#' Plot a suite of quality profiles. Similar to \code{\link{plot_qualityProfile}}
#' but works with paired forward/reverse files.
#'
#' @export
#' @param x a list with file pairings as character vectors
#' @param n integer, the number of profiles (forward and reverse) to plot
#' @param filename character, the name of the output PDF file or NULL to skip file handling
plot_qualityProfiles <- function(x,
  n = 2,
  filename = "quality_profiles.pdf"){

  norev <- length(x$reverse) == 0
  ix <- seq_len(n)
  if (!charlier::is.nullna(filename[1])) { grDevices::pdf(filename) }
  print(dadautils::plot_qualityProfile(x$forward[ix]), filename = NA)
  if (!norev) print(dadautils::plot_qualityProfile(x$reverse[ix]), filename = NA)
  if (!charlier::is.nullna(filename[1]) ){ grDevices::dev.off() }
}


#' Plot quality profiles for one or more FASTQs
#'
#' @export
#' @param x character vector of fastq filenames
#' @param filename character, the name of the PDF file to generate or NA
#' @param ... further arguments for \code{\link[dada2]{plotQualityProfile}}
plot_qualityProfile <- function(x,
  filename = c(NA, "qualityProfile.pdf")[1],
  ...){
  if (!charlier::is.nullna(filename[1])) grDevices::pdf(filename[1])
  ok <- dada2::plotQualityProfile(x, ...)
  if (!charlier::is.nullna(filename[1])) grDevices::dev.off()
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

#' Generate a FASTA object from a vector of sequences
#'
#' @export
#' @param x character vector or a sequence count matrix with named columns
#' @param prefix character, by default 'ASV'
#' @param sep character, by default "_"
#' @param pad logical, if TRUE zero pad the ASV identifier count
#' @param file character or NA.  If not NA then save to the file specified
#' @return \code{Biostrings::BStringSet}
asv_fasta <- function(x, prefix = "ASV", sep = "_", pad = TRUE,
  file = NA){
  if (inherits(x, 'matrix')) x <- colnames(x)
  nx <- length(x)
  id <- if (pad){
    n <- nchar(length(x))
    pattern <- paste0("%s%s%0.",n,"i")
    sprintf(pattern, prefix[1], sep[1], seq_len(nx))
  } else {
    sprintf("%s%s%i", prefix[1], sep[1], seq_len(nx))
  }
  names(x) <- id
  r <- Biostrings::BStringSet(x, use.names = TRUE)
  if (!is.na(file)){
    Biostrings::writeXStringSet(r, file[1], append=FALSE,
                    compress=FALSE, compression_level=NA, format="fasta")
  }
  r
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


#' Count the charcater lengths of sequences
#'
#' @export
#' @param x matrix table of sequences as per \code{\link[dada2]{removeBimeraDenovo}}
#' @param ofile character or NA, if not NA save result to a CSV file
#' @return character counts
sequence_lengths <- function(x,
  ofile = c(NA, "sequence_lengths.csv")[1]){

  n <- nchar(dada2::getSequences(x))
  if (!is.na(ofile)){
    readr::write_csv(dplyr::as_tibble(n), ofile)
  }
  n
}
