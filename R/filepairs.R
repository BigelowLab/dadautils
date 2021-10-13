#' Given a filepair list, make a vector of alternating filenames
#'
#' @export
#' @param x a filepairs list (two element list for forard and reverse filenames)
#' @return a vector of interleaved filenames
interleave_filepairs <- function(x){
  nf <- length(x[[1]])
  nr <- length(x[[2]])
  if (nf != nr) stop("input must two element listing of files of equal length")
  n <- nf + nr
  fi <- seq(from = 1, to = n, by = 2)
  ri <- seq(from = 2, to = n, by = 2)
  r <- rep("", n)
  r[fi] <- x[[1]]
  r[ri] <- x[[2]]
  r
}

#' Given a even numbered vector of filenames, deinterleaf into a filepairs list
#'
#' @export
#' @param x a vector of items
#' @param elements character,the names of the elements
#' @return return a filepairs list
deinterleave_filepairs <- function(x, elements = c("forward", "reverse")){
  n <- length(x)
  if (as.logical(n %% 2)) stop("input must have even number of elements")
  fi <- seq(from = 1, to = n, by = 2)
  ri <- seq(from = 2, to = n, by = 2)
  list(x[fi], x[ri]) %>%
    stats::setNames(elements)
}



#' Verify that a list of file pairs is suitable for processing.
#'
#' This will throw an error if there is a mismatch between the number of
#' elements in the list elements. Otherwise it simply passes the data through.
#'
#' @export
#' @param x a list with file pairings as character vectors
#' @param elements character, the names of the file pair elements to test
#' @param require_reverse logical, if TRUE then insist that reverse files must be present (ala Illumina).
#'   If FALSE (the default) then allow reverse files to be absent (ala PacBio)
#' @param min_size numeric, files with fewer bytes than this number are considered empty.  This is added to 
#'   detect output files that are written but have no content (well, have no meaningful content).
#'   Set to \code{NA} or \code{NULL} to skip this assessment, set to \code{1} to test for at least one byte.
#' @return the input list
verify_filepairs <- function(x, 
  elements = c("forward", "reverse"),
  require_reverse = FALSE,
  min_size = NA){

  if (!all((elements %in% names(x)))) stop("input is missing one or more required elements")

  ll <- lengths(x)

  if (length(ll[2]) == 0){
    if (require_reverse) stop("reverse elements are required")
  } else {
    if (!all(ll %in% ll[1])) stop("elements of input must be the same length")
  }
  
  if (!charlier::is.nullna(min_size[1])){
    ok <- lapply(x, 
      function(x){
        fi <- file.info(x, extra_cols = FALSE)
        if (any(fi$isdir)) stop("one or more elements is a directory - all elements should be files")
        if (any(fi$size < min_size[1])) stop("one or more elements is smaller than min_size:", min_size[1])
      })
  }
  x
}

#' List files and separate into forward and reverse file lists based upon filenaming patterns.
#'
#' If reverse files are not present (ala \code{PacBio}) then that element is set to a zero-length
#' character element.
#'
#' @export
#' @param path character, the input path
#' @param pattern_forward file pattern, the default is to match "^.*R1_001"
#' @param pattern_reverse file pattern, the default is to match either "^.*R2_001"
#' @param patterns_exclude one or more file patterns to exclude, default is "^.*\\.cutadapt_output\\.txt$"
#'        Unlike \code{pattern_forward} and \code{pattern_reverse} this can have multiple elements.  Set to 
#'        NULL of NA to skip.
#' @param glob logical, if \code{TRUE} the input patterns are considered file globs like "*_R1_001.fastq"
#'        and will be converted to regex patterns using \code{\link[utils]{glob2rx}}.
#'        If glob is \code{FALSE} then the the patterns are passed to \code{\link[base]{list.files}}
#'        as provided by the user.
#' @param verify logical, if TRUE test that the filepairs are the same length
#' @param ... further arguments for \code{\link[base]{list.files}}
#' @return named list of sorted forward and reverse filenames (reverse possibly zero-length character)
list_filepairs <- function(path,
                       pattern_forward = "^.*R1_001",
                       pattern_reverse = "^.*R2_001",
                       patterns_exclude = "^.*\\.cutadapt_output\\.txt$",
                       glob = FALSE,
                       verify = TRUE,
                       ...){

  if (glob){
    pattern_forward = utils::glob2rx(pattern_forward)
    pattern_reverse = utils::glob2rx(pattern_reverse)
  }

  x <- list(
    forward = sort(list.files(path, pattern = pattern_forward, full.names = TRUE, ...)),
    reverse = sort(list.files(path, pattern = pattern_reverse, full.names = TRUE, ...)) 
  )
  
  if (!charlier::is.nullna(patterns_exclude[1])){
    if (glob) patterns_exclude = utils::glob2rx(patterns_exclude)
    x <- sapply(x,
      function(x){
        ix <- charlier::mgrepl(patterns_exclude, x, fixed = glob)
        x[!ix]
      }, simplify = FALSE)
  }
  

  if (verify) x <- verify_filepairs(x)

  x
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

