#' Count the length of filepairs
#'
#' @export
#' @param x a list with file pairings as character vectors
#' @return see \code{\link[auntie]{count_filepairs}}
count_filepairs <- function(x){
  .Deprecated("count_filepairs", package = "auntie", msg = "Please use auntie::count_filepairs()")
  auntie::count_filepairs(x)
}

#' Verify that a list of file pairs is suitable for processing.
#'
#' @export
#' @param x a list with file pairings as character vectors
#' @param ... arguments for \code{\link[auntie]{verify_filepairs}}
#' @return see \code{\link[auntie]{verify_filepairs}}
verify_filepairs <- function(x, ...){
  
  .Deprecated("verify_filepairs", package = "auntie", msg = "Please use auntie::verify_filepairs()")
  auntie::verify_filepairs(x, ...)
}

#' Given a filepair list, make a vector of alternating filenames
#'
#' @export
#' @param x a filepairs list (two element list for forard and reverse filenames)
#' @return see \code{\link[auntie]{interleave_filepairs}}
interleave_filepairs <- function(x){
  .Deprecated("interleave_filepairs", package = "auntie", msg = "Please use auntie::interleave_filepairs()")
  auntie::interleave_filepairs(x)
}

#' Given a even numbered vector of filenames, deinterleaf into a filepairs list
#'
#' @export
#' @param x a vector of items
#' @param ... arguments for \code{\link[auntie]{deinterleave_filepairs}}
#' @return see \code{\link[auntie]{deinterleave_filepairs}}
deinterleave_filepairs <- function(x, ...){
  .Deprecated("deinterleave_filepairs", package = "auntie", msg = "Please use auntie::deinterleave_filepairs()")
  auntie::deinterleave_filepairs(x, ...)
}

#' Compute and test the size of filepairs
#' 
#' @export
#' @param x named list of filepairs, typically 'forward' and 'reverse'
#' @param ... arguments for \code{\link[auntie]{size_filepairs}}
#' @return see \code{\link[auntie]{size_filepairs}}
size_filepairs <- function(x, ...){
  .Deprecated("size_filepairs", package = "auntie", msg = "Please use auntie::size_filepairs()")
  auntie::size_filepairs(x, ...)  
}


#' List files and separate into forward and reverse file lists based upon filenaming patterns.
#'
#' If reverse files are not present (ala \code{PacBio}) then that element is set to a zero-length
#' character element.
#'
#' @export
#' @param path character, the input path
#' @param ... arguments for \code{\link[auntie]{list_filepairs}}
#' @return see \code{\link[auntie]{list_filepairs}}
list_filepairs <- function(path, ...){

  .Deprecated("list_filepairs", package = "auntie", msg = "Please use auntie::list_filepairs()")
  auntie::list_filepairs(path, ...)  
}

#' Retrieve an example filepairs listing using dada2 package example data
#'
#' @export
#' @return named list of sorted forward and reverse filenames
example_filepairs  <- function(){
  list(
    forward = c(system.file("extdata", "sam1F.fastq.gz", package="dada2"),
                system.file("extdata", "sam2F.fastq.gz", package="dada2")),
    reverse = c(system.file("extdata", "sam1R.fastq.gz", package="dada2"),
                system.file("extdata", "sam2R.fastq.gz", package="dada2")))
}

#' Read fastq files in a filepairs listing
#'
#' @export
#' @param filelist named list of sorted forward and reverse filenames (reverse may be empty)
#' @return named list of \code{\link[ShortRead]{ShortReadQ}} objects for foreward and reverse
read_fastq_paired <- function(filelist = example_filepairs()){
  sapply(filelist,
      function(x){
        if (length(x) > 0){
          r <- parallel::mclapply(x, ShortRead::readFastq) %>%
            setNames(basename(x))
        } else {
            r <- NULL
        }
        r
      }, simplify = FALSE)
}

#' Compute quality scores using \code{\link[Rsubread]{qualityScores}}
#' 
#' @export
#' @param filelist a paired list forward/reverse filenames (one or both may be empty)
#' @param nreads numeric, see \code{\link[Rsubread]{qualityScores}}
#' @param cleanup logical, if TRUE tidy up the detritus left be Rsubread
#' @return a paired list of quality score matrices or NULL if that component of inputs is empty
paired_quality_scores <- function(filelist = example_filepairs(), nreads = -1,
  cleanup = TRUE) {
  
  rr <- sapply(filelist,
    function(files, nreads = -1){
      if (length(files) > 0){
        r <- parallel::mclapply(files, 
          function(file, nreads = -1) {
            x <- Rsubread::qualityScores(file, nreads = nreads)
          }, nreads = nreads) %>%
          stats::setNames(basename(files))
      } else {
        r <- NULL
      }
      r
    }, nreads = nreads, simplify = FALSE)
    
    if (cleanup) ok <- purge_files()
  rr
}

#' Compute EE per read for \code{\link[Rsubread]{qualityScores}} (one or both may be empty)
#'
#' @export
#' @param x paired list of quality scores as per \code{quality_scores}
#' @param fun the function to use for summary (by default sum)
#' @param na.rm logical, if TRUE then remove NAs before computing \code{fun}
#' @return paired list of numeric vectors with quality per read or NULL where inputs are empty
paired_ee_per_read <- function(x = paired_quality_scores(),
  fun = sum, na.rm = TRUE){
  
  sapply(x,
    function(scores_dir, fun = sum, na.rm = TRUE){
      if (length(scores_dir) > 0){
        r <- lapply(scores_dir,
          function(s, fun = sum, na.rm = TRUE){
            ee <- 10^(-1*s/10)
            apply(ee, 1, fun, na.rm = na.rm)
          }, fun = fun, na.rm = na.rm)
      } else {
        r <- NULL
      }
      r
    }, fun = fun, na.rm = na.rm, simplify = FALSE)
}

#' Threshold paired expected error computations
#' 
#' @export
#' @param x paried list of expected errors - see \code{paired_ee_per_read} (one or both may be empty)
#' @param thresholds numeric vector of one or more thresholds
#' @param op the comparative operator to use for the thrrsholding process (backquoted)
#' @param form character, either 'list' (two elements, table for forward and reverse each) or "table" for a single table
#' @param sample_names character or NULL prepends output with provided names (unless NULL)
#' @param filename character of NULL if not NULL then write a sigle CSV to file, ignored unless \code{form} is 'table'
#' @return a two element list of tibbles with threshold counts of errors-per-read above the threshold, or 
#'   a single table of the two bound together by row or NULL where the input is empty
#' @examples
#' \dontrun{
#'  example_filepairs() %>%
#'    paired_quality_scores() %>%
#'    paired_ee_per_read() %>%
#'    paired_ee_threshold(sample_names  = c("foo", "bar"), filename = "~/my_ee_stuff.csv")
#'  }
paired_ee_threshold <- function(x = paired_ee_per_read(),
  thresholds = c(1,2,3,4,5),
  op = `<=`,
  form = c("list", "table")[2],
  sample_names = NULL,
  filename = NULL){
    
  xx <- sapply(x,
    function(xx, op = `<=`){
      if (length(xx) > 0){
        r <- sapply(xx,
          function(e, op = `<=`, thresholds = c(1,2,3,4,5)){
            sapply(thresholds,
              function(thresh, op = `<=`) sum(op(e,thresh))/length(e) , simplify = FALSE) %>%
              setNames(paste("t", thresholds, sep = "_")) %>%
              dplyr::as_tibble()
          }, op = op, thresholds = thresholds, simplify = FALSE) %>%
          dplyr::bind_rows() %>%
          dplyr::mutate(file = names(xx)) %>%
          dplyr::relocate(file, .before = 1)
      } else {
        r <- NULL
      }
      r
    }, simplify = FALSE)
    
  if (!is.null(sample_names)){
    xx <- sapply(xx,
      function(x, sample_names = "unknown"){
        dplyr::mutate(x, sample = sample_names) %>%
        dplyr::relocate(sample, .before = 1)
      }, 
      simplify = FALSE,
      sample_names = sample_names)
  }
  
  if (tolower(form[1]) == "table"){
    xx <- sapply(names(xx),
      function(name){
        dplyr::mutate(xx[[name]], direction = name) %>%
        dplyr::relocate(.data$direction, .before = 1)
      }, simplify = FALSE) %>%
      dplyr::bind_rows()
      
      if (!is.null(filename[1])) {
        xx <- readr::write_csv(xx, filename[1])
      }
  }
  xx
}

