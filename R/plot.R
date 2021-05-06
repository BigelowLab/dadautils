#' Compute datasets for quality profile
#'
#' This a thin wrapper around \code{\link[dada2]{plotQualityProfile}} but produces a list
#' that includes tables of ancillary statistics.
#'
#' @export
#' @seealso \code{\link[dada2]{plotQualityProfile}}
#' @param fl (Required). \code{character}.
#'        File path(s) to fastq or fastq.gz file(s).
#' @param n (Optional). Default 500,000.
#'        The number of records to sample from the fastq file.
#' @param aggregate (Optional). Default \code{FALSE}.
#'        If TRUE, compute an aggregate quality profile for all fastq files provided.
#' @return as list with the following elements
#' \itemize{
#'   \item{files character vector of input filenames}
#'   \item{opts list of input arguments}
#'   \item{statdf} table of statistical results
#'   \item{anndf} table of annotation info
#'   \item{plotdf} table of stats transformed for plotting
#'   \item{plotdf.summary} table of plot data for plotting (NULL if aggregate is \code{FALSE})
#'   \item{statdf.summary} table of summary stats for plotting  (NULL if aggregate is \code{FALSE})
#' }
#' @examples
#' \dontrun{
#'   files <- c(system.file("extdata", "sam1F.fastq.gz", package="dada2"),
#'              system.file("extdata", "sam1R.fastq.gz", package="dada2"))
#'   result <- quality_profile_data(files)
#'   print(result)
#'  }
quality_profile_data <- function(
  fl = c(system.file("extdata", "sam1F.fastq.gz", package="dada2"),
         system.file("extdata", "sam1R.fastq.gz", package="dada2")), 
         n = 500000, aggregate = FALSE){
  
  if (FALSE){
    fl <- c(system.file("extdata", "sam1F.fastq.gz", package="dada2"),
            system.file("extdata", "sam1R.fastq.gz", package="dada2"))
    n <- 500000
    aggregate <- FALSE
  }
  
  get_quant <- function(xx, yy, q) { xx[which(cumsum(yy)/sum(yy) >= q)][[1]] }
  
  # compute per file
  xx <- lapply(fl,
    function(f){
      srqa <- ShortRead::qa(f, n = n)
      df <- srqa[["perCycle"]]$quality %>%
        dplyr::as_tibble()
      rc <- sum(srqa[["readCounts"]]$read) # Handle aggregate form from qa of a directory
      if (rc >= n) {
        rclabel <- paste("Reads >= ", n)
      } else {
        rclabel <- paste("Reads: ", rc)
      }
      # Calculate summary statistics at each position
      means <- rowsum(df$Score*df$Count, df$Cycle)/rowsum(df$Count, df$Cycle)    
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
      plotdf <- df %>%
        dplyr::mutate(file = basename(f))
      statdf <- dplyr::tibble(
          Cycle=as.integer(rownames(means)),
          Mean=means[,1],
          Q25=as.vector(q25s),
          Q50=as.vector(q50s),
          Q75=as.vector(q75s),
          Cum=10*as.vector(cums)/rc,
          file=basename(f))
    
      anndf <- dplyr::tibble(
          minScore=min(df$Score),
          label=basename(f),
          rclabel=rclabel,
          rc=rc,
          file=basename(f))
      list(plotdf = plotdf, statdf = statdf, anndf = anndf)
    }) # lapply through fl 

    plotdf <- lapply(xx, "[[", "plotdf") %>%
      dplyr::bind_rows()
    statdf <- lapply(xx, "[[", "statdf") %>%
      dplyr::bind_rows()
    anndf <- lapply(xx, "[[", "anndf") %>%
      dplyr::bind_rows() %>%
      dplyr::mutate(minScore = min(.data$minScore))
      
    if (aggregate && length(fl) > 1){
    	plotdf.summary <- aggregate(Count ~ Cycle + Score, data = plotdf, sum) %>%
         dplyr::as_tibble()
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
      statdf.summary <- dplyr::tibble(
        Cycle=as.integer(rownames(means)),
        Mean=means[,1],
        Q25=as.vector(q25s),
        Q50=as.vector(q50s),
        Q75=as.vector(q75s),
        Cum=10*as.vector(cums)/sum(pmin(anndf$rc, n)))
    } else {
      plotdf.summary <- statdf.summary <- NULL
    }
    list(files = fl,
         opts = list(n = n, aggregate = aggregate),
         plotdf = plotdf, 
         statdf = statdf, 
         anndf = anndf, 
         plotdf.summary = plotdf.summary,
         statdf.summary = statdf.summary)
}


#' Generates a plottable object from \code{quality_profile_data} object
#' 
#' @export
#' @param x list as generated by \code{quality_profile_data}
#' @return the input list with a "plot" element added.  
#' @examples
#' \dontrun{
#'   files <- c(system.file("extdata", "sam1F.fastq.gz", package="dada2"),
#'              system.file("extdata", "sam1R.fastq.gz", package="dada2"))
#'   result <- quality_profile_data(files)
#'   result_and_plot <- draw_quality_profile(result)
#'   print(result$plot)
#'  }
draw_quality_profile <- function(x = quality_profile_data()){
  
  if (x$opts$aggregate) {
    	
      p <- ggplot2::ggplot(data=x$plotdf.summary, ggplot2::aes(x=.data$Cycle, y=.data$Score)) +
        ggplot2::geom_tile(ggplot2::aes(fill=.data$Count)) +
  		  ggplot2::scale_fill_gradient(low="#F5F5F5", high="black") +
  		  ggplot2::geom_line(data=x$statdf.summary, ggplot2::aes(y=.data$Mean), color="#66C2A5") +
  		  ggplot2::geom_line(data=x$statdf.summary, ggplot2::aes(y=.data$Q25),
  		    color="#FC8D62", size=0.25, linetype="dashed") +
  		  ggplot2::geom_line(data=x$statdf.summary, ggplot2::aes(y=.data$Q50),
  		    color="#FC8D62", size=0.25) +
  		  ggplot2::geom_line(data=x$statdf.summary, ggplot2::aes(y=.data$Q75),
  		    color="#FC8D62", size=0.25, linetype="dashed") +
  		  ggplot2::ylab("Quality Score") +
  		  ggplot2::xlab("Cycle") +
  		  ggplot2::annotate("text", x=0, y=0,
  		    label=sprintf("Total reads: %d", sum(x$anndf$rc)), color="red", hjust=0) +
  		  ggplot2::theme_bw() +
  		  ggplot2::theme(panel.grid=ggplot2::element_blank()) +
  		  ggplot2::guides(fill=FALSE) +
        ggplot2::facet_wrap(~label)
      if(length(unique(x$statdf$Cum))>1) {
        p <- p +
          ggplot2::geom_line(data=x$statdf.summary, ggplot2::aes(y=.data$Cum),
            color="red", size=0.25, linetype="solid") +
          ggplot2::scale_y_continuous(limits = c(0,NA),
          sec.axis=ggplot2::sec_axis(~.*10, breaks=c(0,100), labels=c("0%", "100%"))) +
          ggplot2::theme(axis.text.y.right = ggplot2::element_text(color = "red"),
            axis.title.y.right = ggplot2::element_text(color = "red"))
      } else {
        p <- p + ggplot2::ylim(c(0,NA))
      }
    } else {
    	p <- ggplot2::ggplot(data=x$plotdf, ggplot2::aes(x=.data$Cycle, y=.data$Score)) +
    	  ggplot2::geom_tile(ggplot2::aes(fill=.data$Count)) +
  		  ggplot2::scale_fill_gradient(low="#F5F5F5", high="black") +
  		  ggplot2::geom_line(data=x$statdf, ggplot2::aes(y=.data$Mean), color="#66C2A5") +
  		  ggplot2::geom_line(data=x$statdf, ggplot2::aes(y=.data$Q25), color="#FC8D62", size=0.25,
  		    linetype="dashed") +
  		  ggplot2::geom_line(data=x$statdf, ggplot2::aes(y=.data$Q50), color="#FC8D62", size=0.25) +
        ggplot2::geom_line(data=x$statdf, ggplot2::aes(y=.data$Q75), color="#FC8D62", size=0.25,
          linetype="dashed") +
        ggplot2::ylab("Quality Score") +
        ggplot2::xlab("Cycle") +
  		  ggplot2::theme_bw() +
  		  ggplot2::theme(panel.grid=ggplot2::element_blank()) +
  		  ggplot2::guides(fill=FALSE) +
  		  ggplot2::geom_text(data=x$anndf, ggplot2::aes(x=0, label=.data$rclabel, y=0),
  		    color="red", hjust=0) +
        ggplot2::facet_wrap(~file)
      if(length(unique(x$statdf$Cum))>1) {
        p <- p +
          ggplot2::geom_line(data=x$statdf, ggplot2::aes(y=.data$Cum),
            color="red", size=0.25, linetype="solid") +
          ggplot2::scale_y_continuous(limits = c(0,NA),
            sec.axis=ggplot2::sec_axis(~.*10, breaks=c(0,100), labels=c("0%", "100%"))) +
          ggplot2::theme(axis.text.y.right = ggplot2::element_text(color = "red"),
            axis.title.y.right = ggplot2::element_text(color = "red"))
      } else {
        p <- p + ggplot2::ylim(c(0,NA))
      }
    }
    
    x[['plot']] <- p
    invisible(x)
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
quality_profile <- function(
  fl = c(system.file("extdata", "sam1F.fastq.gz", package="dada2"),
         system.file("extdata", "sam1R.fastq.gz", package="dada2")), 
  n = 500000, 
  aggregate = FALSE){
  
  if (FALSE){
    fl <- c(system.file("extdata", "sam1F.fastq.gz", package="dada2"),
            system.file("extdata", "sam1R.fastq.gz", package="dada2"))
    n <- 500000
    aggregate <- FALSE
  }
 
 quality_profile_data(fl, n = n, aggregate = aggregate) %>%
   draw_quality_profile()
}



#quality_profile_old <- function(fl, n = 500000, aggregate = FALSE){
#  
#  if (FALSE){
#    fl <- c(system.file("extdata", "sam1F.fastq.gz", package="dada2"),
#            system.file("extdata", "sam1R.fastq.gz", package="dada2"))
#    n <- 500000
#    aggregate <- FALSE
#  }
#  statdf <- data.frame(
#      Cycle=integer(0),
#      Mean=numeric(0),
#      Q25=numeric(0),
#      Q50=numeric(0),
#      Q75=numeric(0),
#      Cum=numeric(0),
#      file=character(0))
#  anndf <- data.frame(
#      minScore=numeric(0),
#      label=character(0),
#      rclabel=character(0),
#      rc=numeric(0),
#      file=character(0))
#
#  FIRST <- TRUE
#  
#  for(f in fl[!is.na(fl)]) {
#    srqa <- ShortRead::qa(f, n = n)
#    df <- srqa[["perCycle"]]$quality
#    rc <- sum(srqa[["readCounts"]]$read) # Handle aggregate form from qa of a directory
#    if (rc >= n) {
#      rclabel <- paste("Reads >= ", n)
#    } else {
#      rclabel <- paste("Reads: ", rc)
#    }
#    # Calculate summary statistics at each position
#    means <- rowsum(df$Score*df$Count, df$Cycle)/rowsum(df$Count, df$Cycle)
#    get_quant <- function(xx, yy, q) { xx[which(cumsum(yy)/sum(yy) >= q)][[1]] }
#    q25s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.25), simplify=TRUE)
#    q50s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.5), simplify=TRUE)
#    q75s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.75), simplify=TRUE)
#    cums <- by(df, df$Cycle, function(foo) sum(foo$Count), simplify=TRUE)
#    if(!all(sapply(list(names(q25s),
#                        names(q50s),
#                        names(q75s),
#                        names(cums)),
#                        identical,
#                        rownames(means)))) {
#      stop("Calculated quantiles/means weren't compatible.")
#    }
#    if(FIRST) {
#      plotdf <- cbind(df, file=basename(f))
#      FIRST <- FALSE
#    } else {
#      plotdf <- rbind(plotdf, cbind(df, file=basename(f)))
#    }
#
#    statdf <- rbind(statdf,
#      data.frame(
#        Cycle=as.integer(rownames(means)),
#        Mean=means,
#        Q25=as.vector(q25s),
#        Q50=as.vector(q50s),
#        Q75=as.vector(q75s),
#        Cum=10*as.vector(cums)/rc,
#        file=basename(f)))
#    anndf <- rbind(anndf,
#      data.frame(
#        minScore=min(df$Score),
#        label=basename(f),
#        rclabel=rclabel,
#        rc=rc,
#        file=basename(f)))
#  } # loop through f (input fl files)
#  anndf$minScore <- min(anndf$minScore)
#  statdf.summary <- plotdf.summary <- NULL
#  
#  if (aggregate){
#  	plotdf.summary <- aggregate(Count ~ Cycle + Score, data = plotdf, sum) %>%
#      dplyr::as_tibble()
#    plotdf.summary$label <- paste(nrow(anndf), "files (aggregated)")
#    means <- rowsum(plotdf.summary$Score*plotdf.summary$Count, plotdf.summary$Cycle)/
#      rowsum(plotdf.summary$Count, plotdf.summary$Cycle)
#    q25s <- by(plotdf.summary, plotdf.summary$Cycle,
#      function(foo) get_quant(foo$Score, foo$Count, 0.25), simplify=TRUE)
#    q50s <- by(plotdf.summary, plotdf.summary$Cycle,
#      function(foo) get_quant(foo$Score, foo$Count, 0.5), simplify=TRUE)
#    q75s <- by(plotdf.summary, plotdf.summary$Cycle,
#      function(foo) get_quant(foo$Score, foo$Count, 0.75), simplify=TRUE)
#    cums <- by(plotdf.summary, plotdf.summary$Cycle,
#      function(foo) sum(foo$Count), simplify=TRUE)
#    statdf.summary <- data.frame(
#      Cycle=as.integer(rownames(means)),
#      Mean=means,
#      Q25=as.vector(q25s),
#      Q50=as.vector(q50s),
#      Q75=as.vector(q75s),
#      Cum=10*as.vector(cums)/sum(pmin(anndf$rc, n)))
#        
#  }
#
#
#  if (aggregate) {
#    	plotdf.summary <- aggregate(Count ~ Cycle + Score, data = plotdf, sum)
#    	plotdf.summary$label <- paste(nrow(anndf), "files (aggregated)")
#    	means <- rowsum(plotdf.summary$Score*plotdf.summary$Count, plotdf.summary$Cycle)/
#    	  rowsum(plotdf.summary$Count, plotdf.summary$Cycle)
#      q25s <- by(plotdf.summary, plotdf.summary$Cycle,
#        function(foo) get_quant(foo$Score, foo$Count, 0.25), simplify=TRUE)
#      q50s <- by(plotdf.summary, plotdf.summary$Cycle,
#        function(foo) get_quant(foo$Score, foo$Count, 0.5), simplify=TRUE)
#      q75s <- by(plotdf.summary, plotdf.summary$Cycle,
#        function(foo) get_quant(foo$Score, foo$Count, 0.75), simplify=TRUE)
#      cums <- by(plotdf.summary, plotdf.summary$Cycle,
#        function(foo) sum(foo$Count), simplify=TRUE)
#      statdf.summary <- data.frame(
#        Cycle=as.integer(rownames(means)),
#        Mean=means,
#        Q25=as.vector(q25s),
#        Q50=as.vector(q50s),
#        Q75=as.vector(q75s),
#        Cum=10*as.vector(cums)/sum(pmin(anndf$rc, n)))
#      p <- ggplot2::ggplot(data=plotdf.summary, ggplot2::aes(x=.data$Cycle, y=.data$Score)) +
#        ggplot2::geom_tile(ggplot2::aes(fill=.data$Count)) +
#  		  ggplot2::scale_fill_gradient(low="#F5F5F5", high="black") +
#  		  ggplot2::geom_line(data=statdf.summary, ggplot2::aes(y=.data$Mean), color="#66C2A5") +
#  		  ggplot2::geom_line(data=statdf.summary, ggplot2::aes(y=.data$Q25),
#  		    color="#FC8D62", size=0.25, linetype="dashed") +
#  		  ggplot2::geom_line(data=statdf.summary, ggplot2::aes(y=.data$Q50),
#  		    color="#FC8D62", size=0.25) +
#  		  ggplot2::geom_line(data=statdf.summary, ggplot2::aes(y=.data$Q75),
#  		    color="#FC8D62", size=0.25, linetype="dashed") +
#  		  ggplot2::ylab("Quality Score") +
#  		  ggplot2::xlab("Cycle") +
#  		  ggplot2::annotate("text", x=0, y=0,
#  		    label=sprintf("Total reads: %d", sum(anndf$rc)), color="red", hjust=0) +
#  		  ggplot2::theme_bw() +
#  		  ggplot2::theme(panel.grid=ggplot2::element_blank()) +
#  		  ggplot2::guides(fill=FALSE) +
#        ggplot2::facet_wrap(~label)
#      if(length(unique(statdf$Cum))>1) {
#        p <- p +
#          ggplot2::geom_line(data=statdf.summary, ggplot2::aes(y=.data$Cum),
#            color="red", size=0.25, linetype="solid") +
#          ggplot2::scale_y_continuous(limits = c(0,NA),
#          sec.axis=ggplot2::sec_axis(~.*10, breaks=c(0,100), labels=c("0%", "100%"))) +
#          ggplot2::theme(axis.text.y.right = ggplot2::element_text(color = "red"),
#            axis.title.y.right = ggplot2::element_text(color = "red"))
#      } else {
#        p <- p + ggplot2::ylim(c(0,NA))
#      }
#    } else {
#    	p <- ggplot2::ggplot(data=plotdf, ggplot2::aes(x=.data$Cycle, y=.data$Score)) +
#    	  ggplot2::geom_tile(ggplot2::aes(fill=.data$Count)) +
#  		  ggplot2::scale_fill_gradient(low="#F5F5F5", high="black") +
#  		  ggplot2::geom_line(data=statdf, ggplot2::aes(y=.data$Mean), color="#66C2A5") +
#  		  ggplot2::geom_line(data=statdf, ggplot2::aes(y=.data$Q25), color="#FC8D62", size=0.25,
#  		    linetype="dashed") +
#  		  ggplot2::geom_line(data=statdf, ggplot2::aes(y=.data$Q50), color="#FC8D62", size=0.25) +
#        ggplot2::geom_line(data=statdf, ggplot2::aes(y=.data$Q75), color="#FC8D62", size=0.25,
#          linetype="dashed") +
#        ggplot2::ylab("Quality Score") +
#        ggplot2::xlab("Cycle") +
#  		  ggplot2::theme_bw() +
#  		  ggplot2::theme(panel.grid=ggplot2::element_blank()) +
#  		  ggplot2::guides(fill=FALSE) +
#  		  ggplot2::geom_text(data=anndf, ggplot2::aes(x=0, label=rclabel, y=0),
#  		    color="red", hjust=0) +
#        ggplot2::facet_wrap(~file)
#      if(length(unique(statdf$Cum))>1) {
#        p <- p +
#          ggplot2::geom_line(data=statdf, ggplot2::aes(y=.data$Cum),
#            color="red", size=0.25, linetype="solid") +
#          ggplot2::scale_y_continuous(limits = c(0,NA),
#            sec.axis=ggplot2::sec_axis(~.*10, breaks=c(0,100), labels=c("0%", "100%"))) +
#          ggplot2::theme(axis.text.y.right = ggplot2::element_text(color = "red"),
#            axis.title.y.right = ggplot2::element_text(color = "red"))
#      } else {
#        p <- p + ggplot2::ylim(c(0,NA))
#      }
#    }
#
#    list(
#      files         = fl,
#      opts          = list(n = n, aggregate = aggregate),
#      stat_data     = dplyr::as_tibble(statdf),
#      ann_data      = dplyr::as_tibble(anndf),
#      plot_data     = dplyr::as_tibble(plotdf),
#      plot_summary  = dplyr::as_tibble(plotdf.summary),
#      plot          = p
#    )
#  }
#