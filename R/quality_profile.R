#' Retrieve statistics from a FastQ file
#'
#' @export
#' @param filename character, one or more filenames
#' @param n (Optional). Default 500,000.
#'        The number of records to sample from the fastq file.
#' @return tibble 
#'   \describe{
#'      \item{filename}{character, the name of the file}
#'      \item{nread numeric}{the number of reads in the file}
#'      \item{ncycle}{list column of numeric vectors, the number of cycles per read}
#'      \item{quantile_cycles}{list columns of numeric vector, quantiles of ncycles from 0.01 to 0.99 by 0.01}
#'      \item{per_cycle_quality}{list columns of tibbles, quality per cycle metrics}
#'    }
fastq_stats <- function(filename = unname(unlist(example_filepairs())) , n = 500000){
  
  if (length(filename) > 1) {
    r <- parallel::mclapply(filename, fastq_stats, n = n) |>
      dplyr::bind_rows()
    return(r)
  }
  
  stopifnot(file.exists(filename[1]))
  x <- try(ShortRead::readFastq(filename[1]))
  if (inherits(x, "try-error")){
    print(x)
    return(NULL)
  }
  
  w <- x |>
      ShortRead::width()
      
  srqa <- ShortRead::qa(list(x) |> setNames(basename(filename[1])), n = n)
    
  dplyr::tibble(filename = filename[1], 
       nread = length(x),
       ncycles = list(w),
       quantile_cycles = list(quantile(w, probs = seq(from = 0.01, to = 0.99, by = 0.01))),
       per_cycle_quality = list(srqa[["perCycle"]]$quality |> dplyr::as_tibble()))
}



#' Compute expected overlap of read pairs
#' 
#' @export
#' @param f numeric, length of forward read
#' @param r numeric, length of reverse read
#' @param a numeric, expected amplicon length
#' @return numeric, estimated overlap of reads
reads_overlap <- function(f = 250, r = 189, a = 400){f - (a-r)}

#' Compute cutoff levels \code{trunLen} for \code{\link[dada2]{filterAndTrim}}
#'
#' @export
#' @param x list of quality score data as per \code{quality_profile_data}
#' @param method character the threshold method to use - just "ruler" method now
#' @param params list (or NULL) of parameters needed to implement the specified method 
#' \itemize{
#'   \item{score numeric,  the threshold score to use with "ruler" method}
#'   \item{model character, the model}
#'   \item{quantile_min numeric, quantile minimum used as a back stop if ruler method fails, set to NA to skip.  Only
#'         used if form is "reduced".}
#' }
#' @param form character either "full" or "reduced".  If reduce return just the threshold row.
#' @param verbose logical, if TRUE then output messages for debugging purposes
#' @return the input list with an added cutoff tibble of 
#' \itemize{
#'  \item{Cycle the predicted Cycle value - the computed cutoff}
#'  \item{Score the fitted Score value}
#'  \item{status "a" for autothreshold, "f" for fixed threshold, "p" for percentile threshold}
#'  \item{model the model as a list column (for example an lm model)}
#'  \item{file the name of the file}
#'  }
#' @examples
#' \dontrun{
#'   files <- c(system.file("extdata", "sam1F.fastq.gz", package="dada2"),
#'              system.file("extdata", "sam1R.fastq.gz", package="dada2"))
#'   result <- quality_profile_data(files) |>
#'    quality_profile_cutoff()
#'  }
quality_profile_cutoff <- function(x = quality_profile_data(),
  method = "ruler",
  params = list(score = 30, model = "Mean ~ poly(Cycle, 2)", quantile_min = 0.9),
  form = c("full", "reduced")[2],
  verbose = FALSE){
  
    # used for devel
    if (FALSE){
      x = quality_profile_data()
      method = "ruler"
      params = list(score = 30, model = "Mean ~ poly(Cycle, 2)", quantile_min = 0.9)
      form = c("full", "reduced")[2]
    }
  
    r <- NULL
    
    if (tolower(method[1]) == "ruler"){
      
      # function to compute model and do cutoffs using Pete/Robin horizontal ruler
      # x1 tibble 
      # key tibble - ignored
      # threshold - where the ruler is placed in Quality Score values
      # model - model as formula or character to be cast as formula
      # form character (if 'reduced' then retrive just the place where the cutoff occurs)
      # qstat as per \code{fastq_stats}
      # min_fraction_above_threshold numeric, minimum fraction of quality scores above the threshold
      # cutoff_adjustment numeric, adjust the cutoff by this amount after it is computed
      fit_ruler <- function(x1, key, 
                            threshold = 30, 
                            quantile_min = 0.9, 
                            model = stats::as.formula("Mean ~ poly(Cycle, 2)"),
                            form = "reduced",
                            qstat = NULL,
                            min_fraction_above_threshold = 0.7,
                            cutoff_adjustment = 0){
                              
                              
        if (FALSE){
          # for devel
          threshold = 30 
          quantile_min = 0.9 
          model = stats::as.formula("Mean ~ poly(Cycle, 2)")
          form = "reduced"
        }
        
        if (verbose) {
          message("fit_ruler: implementing ruler method")
          message(sprintf(" filename: %s", basename(x1$file[1])))
          message(sprintf(" threshold: %0.1f", threshold))
          message(sprintf(" quantile_min: %0.2f", quantile_min))
          message(sprintf(" min_fraction_above_threshold: %0.2f", min_fraction_above_threshold))
          #message(sprintf(" cutoff_adjustment: %0.2f", cutoff_adjustment))
        }
        
        
        if (!inherits(model, "formula")) model <- as.formula(model)
          
        if (!is.null(qstat)) qstat <- qstat |>
                              dplyr::filter(.data$filename == x1$file[1])
            
        f <- lm(model, data = x1)
        p <- predict(f) |> 
          dplyr::as_tibble(rownames = "Cycle") |>
          dplyr::mutate(Cycle = as.numeric(.data$Cycle)) |>
          dplyr::rename(Score = .data$value) 
           
        p <- p |>  
          dplyr::mutate(model = list(f), 
                        file = x1$file,
                        status = "a")
        
        iz <- findInterval(p$Score, threshold[1]) > 0
        
        if ( (sum(iz)/length(iz)) < min_fraction_above_threshold){
          ix <- 1
          p$status <- "a_fail"
        } else {
          ix <- which(iz)
          ix <- ix[length(ix)]
        }
                  
        if (tolower(form[1]) == "reduced"){
          if (!charlier::is.nullna(quantile_min[1])){
            # compute the quantile
            # transform trunLen to an index (into Cycle dimension)
            # compare to ruler method
            # select more permissive trunLen of the two

            q_min <- qstat$quantile_cycles[[1]][sprintf("%i%%", quantile_min*100)]
            iy <- which(p$Cycle <= q_min)
            if (length(iy) > 0){
              iy <- iy[length(iy)]
            } else {
              message("  fit_ruler: all of the Cycles are below q_min")
              iy <- Inf
            }
            
            if (is.infinite(iy)){
              
              p <- p |>
                dplyr::mutate(status = sprintf("p_%0.2f_fail", quantile_min[1]))
              
            } else if (iy < ix){
              
              ix <- iy      # replace the cutoff - update status
              p <- p |>
                dplyr::mutate(status = sprintf("p_%0.2f", quantile_min[1]))
                
            }
          } 
          # at this point ix is either autoselected or user thresholded
          # and status is one of "a", "p_n.nn" or "p_n.nn_fail"
          p <- p |>
            dplyr::slice(ix)
          if (verbose) message(sprintf(" status: %s", p$status[1]))
        }
        if(verbose){
          message(sprintf("  nrows returned: %i", nrow(p)))
          message(sprintf("  Cycle returned: %0.0f", p$Cycle[1]))
          message(sprintf("  Score returned: %0.0f", p$Score[1]))
        }
        
        p
      } # fit ruler
      
      
      # group statdf by file
      # fit Mean ~ Cycle
      # evaluate at params$Score
      # compare to Nth quantile if user provides a value
      # use the lesser of the two
      
      # make full filenames for subsequent ShortRead::qa() 
      # the number if cycles per file may vary, so we need to use a per-file LUT 
      # to cionnect full file path spec for each row in statdf
      full_file_names <- x$files
      names(full_file_names) <- basename(full_file_names)
      
      r <- x$statdf |>
        dplyr::mutate(file = full_file_names[.data$file]) |> 
        dplyr::group_by(.data$file) |>
        dplyr::group_map(fit_ruler, 
                         .keep = TRUE, 
                         threshold = params$score,
                         model = params$model,
                         quantile = params$quantile_min,
                         form = form,
                         qstat = x$quality_stats) |>
        dplyr::bind_rows() |>
        dplyr::mutate(file = basename(.data$file))
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
#' @param stats list of two elements (forward and reverse) with output
#'   of \code{\link{fastq_stats}}
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
#'   \item{quality_stats} per file quality stats as per \code{\link{fastq_stats}}
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
  n = 500000, 
  stats = fastq_stats(fl, n = n), 
  aggregate = FALSE){
  
  if (FALSE){
    fl <- c(system.file("extdata", "sam1F.fastq.gz", package="dada2"),
            system.file("extdata", "sam1R.fastq.gz", package="dada2"))
    n <- 500000
    aggregate <- FALSE
  }
  
  if (!inherits(fl, "character")){
    if (is.list(fl) && all(names(fl) %in% c("forward", "reverse"))) stop("input is a list, did you mean to call quality_profile_pairs()?")
    stop("Input must be vector of file names")
  } else{
    if (!all(file.exists(fl))) stop("one or more files in input not found") 
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
    as.list(c(key$Cycle, x$Score[ix], sumx, meen)) |>
      setNames(c("Cycle", sprintf("Q%i", qs*100), "Cum", "Mean"))
  }
  
  
  
  # compute per file
  xx <- lapply(fl,
    function(f){
      
      s <- stats |>
        dplyr::filter(.data$filename == f)
      df <- (s |>
        dplyr::pull(.data$per_cycle_quality))[[1]]
      rc <- s$nread 
      
      #srqa <- ShortRead::qa(f, n = n)
      #df <- srqa[["perCycle"]]$quality |>
      #  dplyr::as_tibble()
      #rc <- sum(srqa[["readCounts"]]$read) # Handle aggregate form from qa of a directory
      
      if (rc >= n) {
        rclabel <- paste("Reads >= ", n)
      } else {
        rclabel <- paste("Reads: ", rc)
      }
      
      # Calculate summary statistics at each Cycle
      if (FALSE){
        statdf <- df |>
          dplyr::group_by(.data$Cycle) |>
          dplyr::group_map(get_quants, .keep = TRUE) |>
          dplyr::bind_rows() |>
          dplyr::mutate(file = basename(f))
        if (!all(unique(df$Cycle)) %in% qs$Cycle)
          stop("One or more Cycles failed to generate quantiles, means and cums")
        #plotdf <- df |>
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
            Cycle  = as.integer(rownames(means)),
            Mean   = means[,1],
            Q25    = as.vector(q25s),
            Q50    = as.vector(q50s),
            Q75    = as.vector(q75s),
            Cum    = 10*as.vector(cums)/rc,
            file   = basename(f)) 
      } # development block
      
      plotdf <- df |>
        dplyr::mutate(file = basename(f))
      
      # annotations
      anndf <- dplyr::tibble(
          minScore  = min(df$Score),
          label     = basename(f),
          rclabel   = rclabel,
          rc        = rc,
          file      = basename(f))
      
      list(plotdf = plotdf, statdf = statdf, anndf = anndf)
    }) # lapply through fl 

    plotdf <- lapply(xx, "[[", "plotdf") |>
      dplyr::bind_rows()
    statdf <- lapply(xx, "[[", "statdf") |>
      dplyr::bind_rows()
    anndf <- lapply(xx, "[[", "anndf") |>
      dplyr::bind_rows() |>
      dplyr::mutate(minScore = min(.data$minScore))
      
    if (aggregate && length(fl) > 1){
    	plotdf.summary <- aggregate(Count ~ Cycle + Score, data = plotdf, sum) |>
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
         statdf.summary = statdf.summary,
         quality_stats = stats)
}


#' Generates a plottable object from \code{quality_profile_data} object
#'
#' This is not pair savvy.  Plots are faceted by files in the order the come in \code{x$files}
#' 
#' @export
#' @param x list as generated by \code{quality_profile_data}
#' @return a ggplot object
#' @examples
#' \dontrun{
#'   files <- c(system.file("extdata", "sam1F.fastq.gz", package="dada2"),
#'              system.file("extdata", "sam1R.fastq.gz", package="dada2"))
#'   result <- quality_profile_data(files)
#'   result <- quality_profile_drawing(result)
#'   print(result)
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
  		  ggplot2::guides(fill="none") +
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
  		  ggplot2::guides(fill="none") +
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
    
    invisible(p)
}

#' Generates a plottable objects from \code{quality_profile_data} object by filepair (forward/reverse)
#'
#' This is pair savvy.  Plots are faceted by files in the order the come in each element of \code{x$files}
#' (In theory they will be paired horizontally [forward, reverse])
#' 
#' @export
#' @param X named two element list of profile - "forward" and "reverse"
#' @return list ggplot objects  
quality_profile_pairs_drawing <- function(X){
  
  if (!all(names(X) %in% c("forward", "reverse"))){
    stop("input list have two elements:'forward' and 'reverse'")
  }
  npairs <- length(X$forward$files)
  norev <- length(X$reverse$files) == 0
  if (!norev) if(length(X$reverse$files) != npairs) stop("'forward' and 'reverse' must be same number of files")
  if (norev) return(quality_profile_drawing(X[[1]]) ) 
    
    
    
  pp <- lapply(seq_len(npairs),
   function(ipair){
     # slice plotdf, statdf and anndf for each file-pair
     # merge
     # call quality_profile_drawing
     ffile <- basename(X$forward$files[ipair])
     rfile  <-  basename(X$reverse$files[ipair])

     # build a temporary list that binds the f/r pairs together
     x <- list(files = c(ffile , rfile ),
               opts = X$forward$opts,  #presumably forward and reverse are the same
               statdf = X$forward$statdf |> 
                         dplyr::filter(.data$file == ffile) |>
                         dplyr::bind_rows(X$reverse$statdf |> 
                           dplyr::filter(.data$file == rfile)),
               anndf = X$forward$anndf |> 
                         dplyr::filter(.data$file == ffile) |>
                         dplyr::bind_rows(X$reverse$anndf |> 
                           dplyr::filter(.data == rfile)),
               plotdf = X$forward$plotdf |> 
                         dplyr::filter(.data$file == ffile) |>
                         dplyr::bind_rows(X$reverse$plotdf |> 
                           dplyr::filter(.data$file == rfile)))
     
     quality_profile_drawing(x)  
   })
  
   invisible(pp)
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
 
  x <- quality_profile_data(x, n = n, aggregate = aggregate) |>
    quality_profile_cutoff(...) 
  x$plot <- quality_profile_drawing(x)
  invisible(x)
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
#' @param min_overlap numeric, issue warnings for any overlap less than this value 
#' @param plot_filename character or NULL, if provided save the quality profile plots as PDF
#'   NULL to skip writing the file
#' @param overlap_filename character or NULL, if provided save the overlap table to this file as CSV
#'   NULL to skip writing the file
#' @param ... arguments for \code{\link{quality_profile_cutoff}}
#' @return complex list of forward, reverse and merged quality stats 
quality_profile_pairs <- function(
  filelist = list(
    forward = c(system.file("extdata", "sam1F.fastq.gz", package="dada2"),
                system.file("extdata", "sam2F.fastq.gz", package="dada2")),
    reverse = c(system.file("extdata", "sam1R.fastq.gz", package="dada2"),
                system.file("extdata", "sam2R.fastq.gz", package="dada2"))),
  n = 500000, 
  aggregate = FALSE, 
  amplicon_length = 400,
  min_overlap = 20, 
  plot_filename = file.path(dirname(dirname(filelist$forward[1])), "quality_profiles.pdf"),
  overlap_filename = file.path(dirname(dirname(filelist$forward[1])), "overlap.csv"),
  ...){

    if (FALSE){
      n = 500000
      aggregate = FALSE
      amplicon_length = 400
      min_overlap = 20
    }



    xx <- lapply(filelist,
      function(filenames, n = n, aggregate = aggregate, ... ){
        if (length(filenames) > 0){
          x <- quality_profile_data(filenames, n = n, aggregate = aggregate ) |>
            quality_profile_cutoff(...)
        } else {
          x <- NULL
        }
        return(x)
        }, n = n, aggregate = aggregate, ...)
    xx$pairs_plot <- quality_profile_pairs_drawing(xx)
    
    
    # write the plots
    if (!charlier::is.nullna(plot_filename[1])){
      grDevices::pdf(plot_filename[1])
      dummy <- lapply(xx$pairs_plot, print)
      grDevices::dev.off()
    }
    
    
    xx[['overlap']] <- lapply(c("forward", "reverse"), 
      function(n){
        x <- xx[[n]]$cutoff
        setNames(x, paste0(substring(n, 1, 1), names(x)))
        }) |>
        dplyr::bind_cols() |>
        dplyr::mutate(overlap = reads_overlap(f = .data$fCycle, r = .data$rCycle, a = amplicon_length))
     
    ix <- xx[['overlap']]$overlap < min_overlap[1]
    if (any(ix)){
      warning(sprintf("one or more overlaps fall below minimum of %i", min_overlap))
    }
        
    # write a modified overlap file if provided a filename
    if (!charlier::is.nullna(overlap_filename[1])){
      dummy <- xx$overlap |>
        dplyr::select(-.data$fmodel, -.data$rmodel) |>
        # dplyr::mutate(amplicon_length = amplicon_length) |>
        readr::write_csv(overlap_filename[1])
    }  

    invisible(xx)
}
