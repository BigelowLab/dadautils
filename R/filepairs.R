#' Verify that a list of file pairs is suitable for processing.
#'
#' @export
#' @param x a list with file pairings as character vectors
#' @param min_size numeric, files with fewer bytes than this number are considered empty. Set to NA or NULL to ignore file size.
#' @param elements character, the names of the file pair elements to test
#' @param require_reverse logical, if TRUE then insist that reverse files must be present (ala Illumina). If FALSE (the default) then allow reverse files to be absent (ala PacBio)
#' @return the input list with possibly all elements removed if all files fail the size test
verify_filepairs <- function(x, 
  min_size = NULL, 
  elements = c("forward", "reverse"), 
  require_reverse = FALSE){
  
  if (!all((elements %in% names(x)))) stop("input is missing one or more required elements")

  ll <- lengths(x)

  if (length(ll[2]) == 0){
    if (require_reverse) stop("reverse elements are required")
  } else {
    if (!all(ll %in% ll[1])) stop("elements of input must be the same length")
  }
  
  # test each file pair against min_size
  # for any pair with at least one is small, drop it
  # and provide warning
  if (!charlier::is.nullna(min_size[1])){
    
    y <- size_filepairs(x, min_size = min_size[1], verify = FALSE)
    small_ix <- is.na(y$size_pass) | !y$size_pass
    if (all(small_ix)){
      warning("All files are smaller than min_size: ", min_size[1])
      # empty the list
      x <- sapply(x, function(x) {x <- character(0); x}, simplify = FALSE)
    } else if (any(small_ix)){
      x <- y %>% 
        dplyr::filter(!small_ix) %>%
        dplyr::select(dplyr::all_of(elements)) %>%
        as.list()
      y <- y %>% 
        dplyr::filter(small_ix) %>%
        `[[`(1)
      msg <- sprintf("%i small filepairs encountered, dropping: %s", length(y), paste(basename(y), collapse = ", "))
      warning(msg)
    }
  }
  x
}


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





#' Compute and test the size of filepairs
#' 
#' @export
#' @param x named list of filepairs, typically 'forward' and 'reverse'
#' @param min_size numeric or NA, the mininum size of a file required
#' @param verify logical, if TRUE pass the input to \code{verify_filepairs}
#' @param ... other arguments passed to \code{verify_filepairs} except \code{min_size}
#' @return tibble of 
#' \itemize{
#'   \item forward, charcater vector of foreward filenames
#'   \item reverse, charcater vector of reverse filenames
#'   \item forward_size, numeric file sizes in bytes
#'   \item reverse_size, numeric file sizes in bytes
#'   \item size_pass logical, TRUE/FALSE if the size passes and NA when min_size is NA
#' }
size_filepairs <- function(x, min_size = NA, verify = FALSE, ...){
  # TODO make this work for PACBIO
  if (verify) x <- verify_filepairs(x) 
  nm <- names(x)
  sz <- sapply(x, function(x) file.info(x)$size, simplify = FALSE)
  nmsz <- paste0(nm, "_size")
  size_pass <- sz[[1]] > min_size[1] | sz[[2]] > min_size[1]
  dplyr::as_tibble(x) %>%
    dplyr::mutate(
      !!nmsz[1] := sz[[1]],
      !!nmsz[2] := sz[[2]],
      size_pass = size_pass)
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

