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
    ffilt <- charlier::add_extension(ffilt, ext = ".gz", no_dup = TRUE)
    rfilt <- charlier::add_extension(rfilt, ext = ".gz", no_dup = TRUE)
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
#' @param ofile character, the name of the output PDF file
plot_qualityProfiles <- function(x,
  n = 2,
  ofile = "quality_profiles.pdf"){

  ix <- seq_len(n)
  pdf(ofile)
  print(dadautils::plot_qualityProfile(x$forward[ix]), filename = NA)
  print(dadautils::plot_qualityProfile(x$reverse[ix]), filename = NA)
  dev.off()
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
  if (!is.na(filename[1])) grDevices::pdf(filename[1])
  ok <- dada2::plotQualityProfile(x, ...)
  if (!is.na(filename[1])) grDevices::dev.off()
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
    writeXStringSet(r, file[1], append=FALSE,
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


#' Compute quality profiles for one or more FastQ files.
#'
#' This a thin wrapper around \code{\link[dada2]{plotQualityProfile}} but produces a list
#' that includes tables of ancillary statistics in addition to a printable plot object.
#'
#' @export
#' @seealso \code{\link[dada2]{plotQualityProfile}}
#' @param fl (Required). \code{character}.
#'        File path(s) to fastq or fastq.gz file(s).
#' @param n (Optional). Default 500,000.
#'        The number of records to sample from the fastq file.
#' @param aggregate (Optional). Default FALSE.
#'        If TRUE, compute an aggregate quality profile for all fastq files provided.
#' @return as list with the following elements
#' \itemize{
#'   \item{files character vector of input filenames}
#'   \item{opts list of inut arguments}
#'   \item{stat_data} table of statistical results
#'   \item{ann_data} table of annotation info
#'   \item{plot_data} table of stats transformed for plotting
#'   \item{plot_summary} table of summary stats for plotting
#'   \item{plot} ggplot2 object suitable for printing
#' }
#' @examples
#' \dontrun{
#'   result <- quality_profile(system.file("extdata", "sam1F.fastq.gz", package="dada2"))
#'   print(result$plot)
#'  }
quality_profile <- function(fl, n = 500000, aggregate = FALSE){
  statdf <- data.frame(
      Cycle=integer(0),
      Mean=numeric(0),
      Q25=numeric(0),
      Q50=numeric(0),
      Q75=numeric(0),
      Cum=numeric(0),
      file=character(0))
  anndf <- data.frame(
      minScore=numeric(0),
      label=character(0),
      rclabel=character(0),
      rc=numeric(0),
      file=character(0))

  FIRST <- TRUE
  for(f in fl[!is.na(fl)]) {
    srqa <- ShortRead::qa(f, n = n)
    df <- srqa[["perCycle"]]$quality
    rc <- sum(srqa[["readCounts"]]$read) # Handle aggregate form from qa of a directory
    if (rc >= n) {
      rclabel <- paste("Reads >= ", n)
    } else {
      rclabel <- paste("Reads: ", rc)
    }
    # Calculate summary statistics at each position
    means <- rowsum(df$Score*df$Count, df$Cycle)/rowsum(df$Count, df$Cycle)
    get_quant <- function(xx, yy, q) { xx[which(cumsum(yy)/sum(yy) >=q)][[1]] }
    q25s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.25), simplify=TRUE)
    q50s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.5), simplify=TRUE)
    q75s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.75), simplify=TRUE)
    cums <- by(df, df$Cycle, function(foo) sum(foo$Count), simplify=TRUE)
    if(!all(sapply(list(names(q25s),
                        names(q50s),
                        names(q75s),
                        names(cums)),
                        identical,
                        rownames(means)))) {
      stop("Calculated quantiles/means weren't compatible.")
    }
    if(FIRST) {
      plotdf <- cbind(df, file=basename(f))
      FIRST <- FALSE
    } else {
      plotdf <- rbind(plotdf, cbind(df, file=basename(f)))
    }

    statdf <- rbind(statdf,
      data.frame(
        Cycle=as.integer(rownames(means)),
        Mean=means,
        Q25=as.vector(q25s),
        Q50=as.vector(q50s),
        Q75=as.vector(q75s),
        Cum=10*as.vector(cums)/rc,
        file=basename(f)))
    anndf <- rbind(anndf,
      data.frame(
        minScore=min(df$Score),
        label=basename(f),
        rclabel=rclabel,
        rc=rc,
        file=basename(f)))
  }
  anndf$minScore <- min(anndf$minScore)
  statdf.summary <- plotdf.summary <- NULL

  if (aggregate) {
    	plotdf.summary <- aggregate(.data$Count ~ .data$Cycle + .data$Score, data = plotdf, sum)
    	plotdf.summary$label <- paste(nrow(anndf), "files (aggregated)")
    	means <- rowsum(plotdf.summary$Score*plotdf.summary$Count, plotdf.summary$Cycle)/
    	  rowsum(plotdf.summary$Count, plotdf.summary$Cycle)
      q25s <- by(plotdf.summary, plotdf.summary$Cycle,
        function(foo) get_quant(foo$Score, foo$Count, 0.25), simplify=TRUE)
      q50s <- by(plotdf.summary, plotdf.summary$Cycle,
        function(foo) get_quant(foo$Score, foo$Count, 0.5), simplify=TRUE)
      q75s <- by(plotdf.summary, plotdf.summary$Cycle,
        function(foo) get_quant(foo$Score, foo$Count, 0.75), simplify=TRUE)
      cums <- by(plotdf.summary, plotdf.summary$Cycle,
        function(foo) sum(foo$Count), simplify=TRUE)
      statdf.summary <- data.frame(
        Cycle=as.integer(rownames(means)),
        Mean=means,
        Q25=as.vector(q25s),
        Q50=as.vector(q50s),
        Q75=as.vector(q75s),
        Cum=10*as.vector(cums)/sum(pmin(anndf$rc, n)))
      p <- ggplot2::ggplot(data=plotdf.summary, ggplot2::aes(x=.data$Cycle, y=.data$Score)) +
        ggplot2::geom_tile(ggplot2::aes(fill=.data$Count)) +
  		  ggplot2::scale_fill_gradient(low="#F5F5F5", high="black") +
  		  ggplot2::geom_line(data=statdf.summary, ggplot2::aes(y=.data$Mean), color="#66C2A5") +
  		  ggplot2::geom_line(data=statdf.summary, ggplot2::aes(y=.data$Q25),
  		    color="#FC8D62", size=0.25, linetype="dashed") +
  		  ggplot2::geom_line(data=statdf.summary, ggplot2::aes(y=.data$Q50),
  		    color="#FC8D62", size=0.25) +
  		  ggplot2::geom_line(data=statdf.summary, ggplot2::aes(y=.data$Q75),
  		    color="#FC8D62", size=0.25, linetype="dashed") +
  		  ggplot2::ylab("Quality Score") +
  		  ggplot2::xlab("Cycle") +
  		  ggplot2::annotate("text", x=0, y=0,
  		    label=sprintf("Total reads: %d", sum(anndf$rc)), color="red", hjust=0) +
  		  ggplot2::theme_bw() +
  		  ggplot2::theme(panel.grid=ggplot2::element_blank()) +
  		  ggplot2::guides(fill=FALSE) +
        ggplot2::facet_wrap(~label)
      if(length(unique(statdf$Cum))>1) {
        p <- p +
          ggplot2::geom_line(data=statdf.summary, ggplot2::aes(y=.data$Cum),
            color="red", size=0.25, linetype="solid") +
          ggplot2::scale_y_continuous(limits = c(0,NA),
          sec.axis=ggplot2::sec_axis(~.*10, breaks=c(0,100), labels=c("0%", "100%"))) +
          ggplot2::theme(axis.text.y.right = ggplot2::element_text(color = "red"),
            axis.title.y.right = ggplot2::element_text(color = "red"))
      } else {
        p <- p + ggplot2::ylim(c(0,NA))
      }
    } else {
    	p <- ggplot2::ggplot(data=plotdf, ggplot2::aes(x=.data$Cycle, y=.data$Score)) +
    	  ggplot2::geom_tile(ggplot2::aes(fill=.data$Count)) +
  		  ggplot2::scale_fill_gradient(low="#F5F5F5", high="black") +
  		  ggplot2::geom_line(data=statdf, ggplot2::aes(y=.data$Mean), color="#66C2A5") +
  		  ggplot2::geom_line(data=statdf, ggplot2::aes(y=.data$Q25), color="#FC8D62", size=0.25,
  		    linetype="dashed") +
  		  ggplot2::geom_line(data=statdf, ggplot2::aes(y=.data$Q50), color="#FC8D62", size=0.25) +
        ggplot2::geom_line(data=statdf, ggplot2::aes(y=.data$Q75), color="#FC8D62", size=0.25,
          linetype="dashed") +
        ggplot2::ylab("Quality Score") +
        ggplot2::xlab("Cycle") +
  		  ggplot2::theme_bw() +
  		  ggplot2::theme(panel.grid=ggplot2::element_blank()) +
  		  ggplot2::guides(fill=FALSE) +
  		  ggplot2::geom_text(data=anndf, ggplot2::aes(x=0, label=rclabel, y=0),
  		    color="red", hjust=0) +
        ggplot2::facet_wrap(~file)
      if(length(unique(statdf$Cum))>1) {
        p <- p +
          ggplot2::geom_line(data=statdf, ggplot2::aes(y=.data$Cum),
            color="red", size=0.25, linetype="solid") +
          ggplot2::scale_y_continuous(limits = c(0,NA),
            sec.axis=ggplot2::sec_axis(~.*10, breaks=c(0,100), labels=c("0%", "100%"))) +
          ggplot2::theme(axis.text.y.right = ggplot2::element_text(color = "red"),
            axis.title.y.right = ggplot2::element_text(color = "red"))
      } else {
        p <- p + ggplot2::ylim(c(0,NA))
      }
    }

    list(
      files         = fl,
      opts          = list(n = n, aggregate = aggregate),
      stat_data     = dplyr::as_tibble(statdf),
      ann_data      = dplyr::as_tibble(anndf),
      plot_data     = dplyr::as_tibble(plotdf),
      plot_summary  = dplyr::as_tibble(plotdf.summary),
      plot          = p
    )
  }
