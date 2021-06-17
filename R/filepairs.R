
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
#' @return the input list
verify_filepairs <- function(x, elements = c("forward", "reverse"),
  require_reverse = FALSE){

  if (!all((elements %in% names(x)))) stop("input is missing one or more required elements")

  ll <- lengths(x)

  if (length(ll[2]) == 0){
    if (require_reverse) stop("reverse elements are required")
  } else {
    if (!all(ll %in% ll[1])) stop("elements of input must be the same length")
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
#' @param pattern_forward file pattern, the default is to match either "_R1_001.fastq" or "_R1_001.fastq.gz"
#' @param pattern_reverse file pattern, the default is to match either "_R2_001.fastq" or "_R2_001.fastq.gz"
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
          r <- lapply(x, ShortRead::readFastq) %>%
            setNames(basename(x))
        } else {
            r <- NULL
        }
        r
      }, simplify = FALSE)
}

#' Compute scores for \code{\link[ShortRead]{ShortReadQ}} objects listing
#'
#' @seealso \code{\link[ShortRead]{alphabetScore}}
#'
#' @export
#' @param x named list of sorted forward and reverse \code{\link[ShortRead]{ShortReadQ}} objects (one or both maybe empty)
#' @param fun the function to score by  either \code{ShortRead::alphabetScore} or \code{ShortRead::alphabetByCycle}
#' @return named list of quality score vectors as returned by \code{\link[ShortRead]{alphabetScore}} OR
#'  named list of quality score arrays as returned by \code{\link[ShortRead]{alphabetByCycle}} OR
#'  named list of empty elements where input obects may be empty
score_paired <- function(x = read_fastq_paired(),
  fun = ShortRead::alphabetScore){
  sapply(x,
    function(x){
      if (length(x) > 0){
        r <- sapply(x, fun, simplify = FALSE)
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
#' @param ... other arguments for \code{\link[Rsubread]{qualityScores}}
#' @return a paired list of quality score matrices or NULL if that component of inputs is empty
paired_quality_scores <- function(filelist = example_filepairs(), nreads = -1, ...) {
  sapply(filelist,
    function(files){
      if (length(files) > 0){
        r <- apply(files, 
          function(file) {
            suppressMessages(Rsubread::qualityScores(file, nreads = nreads, ...))
          }) %>%
          stats::setNames(basename(files))
      } else {
        r <- NULL
      }
      r
    }, simplify = FALSE)
}

#' Compute EE per read for \code{\link[Rsubread]{qualityScores}} (one or both may be empty)
#'
#' @export
#' @param x paired list of quality scores as per \code{quality_scores}
#' @param fun the function to use for summary (by default sum)
#' @param na.rm logical, if TRUE then remove NAs before computing \code{fun}
#' @return paried list of numeric vectors with quality per read or NULL where inputs are empty
paired_ee_per_read <- function(x = paired_quality_scores(),
  fun = sum, na.rm = TRUE){
  
  sapply(x,
    function(scores_dir){
      if (length(scores_dir) > 0){
        r <- lapply(scores_dir,
          function(s){
            ee <- 10^(-s/10)
            apply(ee, 1, fun, na.rm = na.rm)
          })
      } else {
        r <- NULL
      }
      r
    }
    ,simplify = FALSE)
}

#' Threshold paired expected error computations
#' 
#' @export
#' @param x paried list of expected errors - see \code{paired_ee_per_read} (one or both may be empty)
#' @param thresholds numeric vector of one or more thresholds
#' @param op the comparative operator to use for the thrrsholding process (backquoted)
#' @return a two element list of tibbles with threshold counts of errors-per-read above the threshold
#'   or NULL where the input is empty
paired_ee_threshold <- function(x = paired_ee_per_read(),
  thresholds = c(1,2,3,4,5),
  op = `<=`){
  # sapply(ee, function(ee) sapply(ee, function(e) sum(e <= 2)), simplify = F)
  sapply(x,
    function(xx){
      if (length(xx) > 0){
        r <- sapply(xx,
          function(e){
            sapply(thresholds,
              function(thresh) sum(op(e,thresh)) , simplify = FALSE) %>%
              setNames(paste("t", thresholds, sep = "_")) %>%
              dplyr::as_tibble()
          }, simplify = FALSE) %>%
          dplyr::bind_rows() %>%
          dplyr::mutate(file = names(xx)) %>%
          dplyr::relocate(file, .before = 1)
      } else {
        r <- NULL
      }
      r
    }, simplify = FALSE)
}




# #' Compute expected errors from a vector - presumably \code{\link[ShortRead]{alphabetScore}} vector
# #'
# #' @seealso \code{\link[dada2]{filterAndTrim}}
# #' @export
# #' @param x numeric vector such as from \code{\link[ShortRead]{alphabetScore}}
# #' @return numeric expected error 
# expected_error <- function(x){
#   sum(10^(-x/10))
# }
# 
# #' Compute expected errors for \code{\link[ShortRead]{alphabetScore}} vectors
# #'
# #' @seealso \code{\link[dada2]{filterAndTrim}}
# #'
# #' @export
# #' @param x named list of quality score vectors as returned by \code{\link[ShortRead]{alphabetScore}}
# #' @return named list expected errors
# ee_paired <- function(x = score_paired()){
#   sapply(x,
#     function(x){
#       sapply(x, expected_error)
#     }, simplify = FALSE)
# }    
# 
# #' Prepare a summary of expected errors
# #'
# #' @export
# #' @param filelist named list of sorted forward and reverse filenames
# #' @return named list of expected error vectors (one for forward, one for reverse)
# #' @examples 
# #' \dontrun{
# #'  example_filepairs() %>%
# #'    expected_error_paired()
# #'
# #'  # $forward
# #'  # sam1F.fastq.gz sam2F.fastq.gz 
# #'  #              0              0 
# #'  # 
# #'  # $reverse
# #'  # sam1R.fastq.gz sam2R.fastq.gz 
# #'                 0              0                
# #' }
# expected_error_paired <- function(filelist = example_filepairs()){
#   filelist %>%
#     read_fastq_paired() %>%
#     score_paired() %>%
#     ee_paired()
# }