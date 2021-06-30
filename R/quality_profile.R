#' Compute cutoff levels \code{trunLen} for \code{\link[dada2]{filterAndTrim}}
#'
#' @export
#' @param x list of quality score data as per \code{quality_profile_data}
#' @param method character the threshold method to use - just "ruler" method now
#' @param params list (or NULL) of parameters needed to implement the specified method 
#' @param form character either "full" or "reduced".  If reduce return just the threshold row.
#' @return the input list with an added cutoff tibble of 
#' \itemize{
#'  \item{Cycle the predicted Cycle value - the computed cutoff}
#'  \item{Score the fitted Score value}
#'  \item{overlap the expected overlap at the cutoff (forward cutoff len + )}
#'  \item{model the model as a list column (for example an lm model)}
#'  \item{file the name of the file}
#'  }
#' @examples
#' \dontrun{
#'   files <- c(system.file("extdata", "sam1F.fastq.gz", package="dada2"),
#'              system.file("extdata", "sam1R.fastq.gz", package="dada2"))
#'   result <- quality_profile_data(files) %>%
#'    quality_profile_cutoff()
#'  }
quality_profile_cutoff <- function(x = quality_profile_data(),
  method = "ruler",
  params = list(score = 30, model = "Mean ~ poly(Cycle, 2)"),
  form = c("full", "reduced")[2]){
   
    r <- NULL
    if (tolower(method[1]) == "ruler"){
      
      #' function to compute model and do cutoffs using Pete/Robin horizontal ruler
      #' x1 tibble 
      #' key tibble - ignored
      #' threshold - where the ruler is placed in Quality Score values
      #' model - model as formula or character to be cast as formula
      #' form character (if 'reduced' then retrive just the place where the cutoff occurs)
      fit_ruler <- function(x1, key, 
                            threshold = 30, 
                            model = stats::as.formula("Mean ~ poly(Cycle, 2)"),
                            form = "reduced"){
        if (!inherits(model, "formula")) model <- as.formula(model)
        f <- lm(model, data = x1)
        p <- predict(f) %>% 
          dplyr::as_tibble(rownames = "Cycle") %>%
          dplyr::mutate(Cycle = as.numeric(.data$Cycle)) %>%
          dplyr::rename(Score = .data$value) 
        ix <- which(findInterval(p$Score, threshold) > 0)
        ix <- ix[length(ix)]
        p <- p %>%  
          dplyr::mutate(model = list(f), 
                        file = x1$file)
        if (tolower(form[1]) == "reduced"){
          p <- p %>%
            dplyr::slice(ix)
        }
        p
      }
      # group statdf by file
      # fit Mean ~ Cycle
      # evaluate at params$Score
      r <- x$statdf %>%
        dplyr::group_by(.data$file) %>%
        dplyr::group_map(fit_ruler, .keep = TRUE, 
                         threshold = params$score,
                         model = params$model,
                         form = form) %>%
        dplyr::bind_rows()
    } else {
      stop("method not known:", method[1])
    }
    x[['cutoff']] <- r
    x
 }


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
  
  # retrieve quantile (xx value) where yy quantile is at or above q
  get_quant <- function(xx, yy, q) { xx[which(cumsum(yy)/sum(yy) >= q)][[1]] }
  
  # retrieve list of Cycle, Q25, Q50, Q75, Cum, Mean per Cycle-group data.frame
  # this replaces the get_quant and repeated by() statements which are awesome and fast
  # but also require repeated computation of the same vectors (sumsum, y, sum)
  # so this reduces that load.  There is a tradeoff as we construct a (named) list each for 
  # each Cycle, but the ease of programmatic structures seems like a nice trade
  get_quants <- function(x, key, qs = c(0.25, 0.50, 0.75)) { 
    cumsumx <- cumsum(x$Count)
    sumx <- cumsumx[length(cumsumx)]
    y <- cumsumx/sumx
    ix <- sapply(qs, function(qs) which(y >= qs)[1])
    #cum <- sum(x$Count)
    meen <- sum(x$Score * x$Count)/sumx
    as.list(c(key$Cycle, x$Score[ix], sumx, meen)) %>%
      setNames(c("Cycle", sprintf("Q%i", qs*100), "Cum", "Mean"))
  }
  
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
      
      # Calculate summary statistics at each Cycle
      if (FALSE){
        statdf <- df %>%
          dplyr::group_by(.data$Cycle) %>%
          dplyr::group_map(get_quants, .keep = TRUE) %>%
          dplyr::bind_rows() %>%
          dplyr::mutate(file = basename(f))
        if (!all(unique(df$Cycle)) %in% qs$Cycle)
          stop("One or more Cycles failed to generate quantiles, means and cums")
        #plotdf <- df %>%
        #  dplyr::mutate(file = basename(f))
      } else {
        means <- rowsum(df$Score*df$Count, df$Cycle)/rowsum(df$Count, df$Cycle)  
        cums <- by(df, df$Cycle, function(foo) sum(foo$Count), simplify=TRUE)  
        q25s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.25), simplify=TRUE)
        q50s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.5), simplify=TRUE)
        q75s <- by(df, df$Cycle, function(foo) get_quant(foo$Score, foo$Count, 0.75), simplify=TRUE)
        if(!all(sapply(list(names(q25s),
                            names(q50s),
                            names(q75s),
                            names(cums)),
                            identical,
                            rownames(means)))) {
            stop("Calculated quantiles/means weren't compatible.")
          }
        #plotdf <- cbind(df, file=basename(f))
        statdf <- dplyr::tibble(
            Cycle=as.integer(rownames(means)),
            Mean=means[,1],
            Q25=as.vector(q25s),
            Q50=as.vector(q50s),
            Q75=as.vector(q75s),
            Cum=10*as.vector(cums)/rc,
            file=basename(f))
      } # development block
      
      plotdf <- df %>%
        dplyr::mutate(file = basename(f))

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
#'   result_and_plot <- quality_profile_drawing(result)
#'   print(result$plot)
#'  }
quality_profile_drawing <- function(x = quality_profile_data()){
  
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
        ggplot2::geom_smooth(data=x$statdf, ggplot2::aes(y=.data$Mean), 
                             method = lm, formula = y ~ poly(x, 2), se = TRUE) + 
  		  ggplot2::geom_line(data=x$statdf, ggplot2::aes(y=.data$Q25), color="#FC8D62", size=0.25,
  		    linetype="dashed") +
  		  ggplot2::geom_line(data=x$statdf, ggplot2::aes(y=.data$Q50), color="#66C2A5", size=0.25) +
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
#' @param x (Required). \code{character}.
#'        File path(s) to fastq or fastq.gz file(s).
#' @param n (Optional). Default 500,000.
#'        The number of records to sample from the fastq file.
#' @param aggregate (Optional). Default FALSE.
#'        If TRUE, compute an aggregate quality profile for all fastq files provided.
#' @param ... arguments for \code{\link{quality_profile_cutoff}}
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
  x = c(system.file("extdata", "sam1F.fastq.gz", package="dada2"),
         system.file("extdata", "sam1R.fastq.gz", package="dada2")), 
  n = 500000, 
  aggregate = FALSE, 
  ...){
  
  if (FALSE){
    x <- c(system.file("extdata", "sam1F.fastq.gz", package="dada2"),
            system.file("extdata", "sam1R.fastq.gz", package="dada2"))
    n <- 500000
    aggregate <- FALSE
  }
 
 quality_profile_data(x, n = n, aggregate = aggregate) %>%
   quality_profile_cutoff(...) %>%
   quality_profile_drawing()
}

#' Compute pairwise (forward-reverse) quality profiles
#'
#' This a thin wrapper around \code{\link[dada2]{plotQualityProfile}} but produces a list
#' that includes tables of ancillary statistics in addition to a printable plot object.
#'
#' @export
#' @seealso \code{\link[dada2]{plotQualityProfile}}
#' @param filelist list of forward and reverse fastq files
#' @param n (Optional). Default 500,000.
#'        The number of records to sample from the fastq file.
#' @param aggregate (Optional). Default FALSE.
#'        If TRUE, compute an aggregate quality profile for all fastq files provided.
#' @param amplicon_length numeric, the expected amplicon length, used to compute expected overlap
#' @param ... arguments for \code{\link{quality_profile_cutoff}}
quality_profile_pairs <- function(
  filelist = list(
    forward = c(system.file("extdata", "sam1F.fastq.gz", package="dada2"),
                system.file("extdata", "sam2F.fastq.gz", package="dada2")),
    reverse = c(system.file("extdata", "sam1R.fastq.gz", package="dada2"),
                system.file("extdata", "sam2R.fastq.gz", package="dada2"))),
  n = 500000, 
  aggregate = FALSE, 
  amplicon_length = 400,
  ...){

  if (FALSE){
    filelist = list(
      forward = c(system.file("extdata", "sam1F.fastq.gz", package="dada2"),
                  system.file("extdata", "sam2F.fastq.gz", package="dada2")),
      reverse = c(system.file("extdata", "sam1R.fastq.gz", package="dada2"),
                  system.file("extdata", "sam2R.fastq.gz", package="dada2")))
    n = 500000
    aggregate = FALSE
    amplicon_length = 400
    xx <- lapply(filelist, quality_profile, n = n, aggregate = aggregate)
  }
    xx <- lapply(filelist, quality_profile, n = n, aggregate = aggregate, ...)
    
    overlap <- function(f = 250, r = 189, a = 400){f - (a-r)}
    
    xx[['overlap']] <- lapply(names(xx), 
      function(n){
        x <- xx[[n]]$cutoff
        setNames(x, paste0(substring(n, 1, 1), names(x)))
        }) %>%
        dplyr::bind_cols() %>%
        dplyr::mutate(overlap = overlap(f = .data$fCycle, r = .data$rCycle, a = amplicon_length))
    invisible(xx)
}
