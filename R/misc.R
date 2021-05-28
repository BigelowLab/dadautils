#' Retrieve an integer vector of known logger levels for \code{futile.logger}
#'
#' @export
#' @return a named integer vector
logger_levels <- function(){
  
  .Deprecated("logger_levels", package = "charlier",
    msg = "Please use charlier::logger_levels")
  c(FATAL = 1L, ERROR = 2L, WARN = 4L, INFO = 6L, DEBUG = 8L, TRACE = 9L)
}

#' A tidy wrapper around \code{\link[base]{saveRDS}}
#'
#' @export
#' @param object the R object to be saved
#' @param file character or file connection
#' @param ... other arguments for \code{\link[base]{saveRDS}}
#' @return the original object invisibly
#' @examples
#' \dontrun{
#'   x <- my_big_process() %>%
#'     write_RDS("my_big_data.rds") %>%
#'     dplyr::filter(foo <= 7) %>%
#'     readr::write_csv("my_smaller_data.csv")
#' }
write_RDS <- function(object, file = "", ...){

  .Deprecated("write_RDS", package = "charlier",
    msg = "Please use charlier::write_RDS")
  saveRDS(object, file = file, ...)
  invisible(object)
}

#' A wrapper around \code{\link[base]{readRDS}}
#'
#' @export
#' @param file the filename to read
#' @param ... other arguments for \code{\link[base]{readRDS}}
#' @return R object
read_RDS <- function(file, ...){
  .Deprecated("read_RDS", package = "charlier",
    msg = "Please use charlier::read_RDS")
  readRDS(file, ...)
}



#' Retrieve the PBS JIOBID if available
#'
#' @export
#' @param no_pbs_text character, the text to return if no jobid found
#' @return character jobid or missing_text value
get_pbs_jobid <- function(no_pbs_text = ""){
  
  .Deprecated("get_pbs_jobid", package = "charlier",
    msg = "Please use charlier::get_pbs_jobid")
    
  PBS_JOBID <- Sys.getenv("PBS_JOBID")
  if (nchar(PBS_JOBID) == 0) PBS_JOBID <- no_pbs_text
  PBS_JOBID
}


#' Provide an R session audit
#'
#' @export
#' @param filename character, the name of the file to dump to or "" to print to console
#' @param pbs_jobid character, the OPBS jobid if known
#' @return NULL invisibly
audit <- function(filename = "", pbs_jobid = get_pbs_jobid()){
  
  .Deprecated("audit", package = "charlier",
    msg = "Please use charlier::audit")
  
  cat("Audit date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S", usetz = TRUE), "\n",
    file = filename)
  cat("System PID:", Sys.getpid(), "\n", file = filename, append = TRUE)
  cat("PBS_JOBID:", pbs_jobid, "\n", file = filename, append = TRUE)
  cat("Cores:", count_cores(), "\n", file = filename, append = TRUE)
  cat("R version:", R.version.string, "\n",
    file = filename, append = TRUE)
  cat("libPaths():\n",
    file = filename, append = TRUE)
  for (lp in .libPaths()) cat("    ", lp, "\n",
                              file = filename, append = TRUE)
  x <- installed_packages() %>%
    dplyr::select("Package", "Version",  "LibPath")
  cat("installed.packages():\n",
      file = filename, append = TRUE)
  if (nzchar(filename)){
    conn <- file(filename, open = 'at')
    utils::write.csv(x, file = conn, row.names = FALSE)
    close(conn)
  } else {
    print(x, row.names = FALSE)
  }
  invisible(NULL)
}

#' A wrapper around \code{installed.packages} the returns a table
#'
#' @export
#' @param ... more arguments for \code{installed.packages}
#' @return table
installed_packages <- function(...){
  dplyr::as_tibble(installed.packages(...))
}


#' A wrapper around \code{old.packages} the returns a table
#'
#' @export
#' @param ... more arguments for \code{old.packages}
#' @return table
old_packages <- function(...){
  dplyr::as_tibble(utils::old.packages(...))
}

#' Retrieve a table of packages suitable for updating
#'
#' @export
#' @param lib_pattern character, A glob pattern specification that is used to filter
#'    library locations
#' @return a table of packages identified for update
identify_upgrades <- function(lib_pattern = "^/mnt/modules/bin/dada2"){

  old_packages() %>%
    dplyr::filter(grepl(lib_pattern, .data$LibPath))
}

#' Count the number of CPUs
#'
#' If the number of CPUS has been specified for a PBS session then retrieves the values
#' of environment variable '$NCPUS' otherwise this is wrap of \code{\link[parallel]{detectCores}}
#'
#' @export
#' @return integer count of cores
count_cores <- function(){
  .Deprecated("count_cores", package = "charlier",
    msg = "Please use charlier::count_cores")
  ncpus <- Sys.getenv("NCPUS")
  if (nchar(ncpus) == 0){
   ncpus <- parallel::detectCores()
  } else {
    ncpus <- as.integer(ncpus[1])
  }
  ncpus
}

#' Extract sample names fromn a list of filenames
#'
#' @export
#' @param x character, vector of one or more filenames.  This also could be a
#'   list ala \code{\link{list_filepairs}} in which case it must have a
#'   the specified element name \code{name}
#' @param rule character, the name of the rule to use
#' @param element character, the name of the list element used to find filenames if \code{x} is a list
#' @param custom_fun, a function to do custom extraction, ignored unless \code{rule} is "custom"
#' @return character vector of sample names - one per input
extract_sample_names <- function(x,
  rule = c("before first _","basename","custom")[1],
  element = "forward",
  custom_fun = NULL){

  if (is.list(x)){
    if (!(element[1] %in% names(x))){
      stop("if x is a list it must have a ", element[1],  "element")
    }
    x <- x[[element[1]]]
  }
  sample_names = switch(rule[1],
      "before first _" = sapply(strsplit(basename(x), "_"), `[`, 1),
      "basename" = basename(x),
      "custom" = custom_fun(x),
       stop("rule not known - check with progammer to add a rule")
    )

  sample_names
}



#' Strip the extension(s) off of a filename
#'
#' Note if ext is ".fastq" then ".fastq.gz" and ".fastq.tar.gz" will also be stripped
#'
#' @export
#' @param filename character one or more filenames
#' @param ext character, one or more extension patterns
#' @return filename with extension stripped
strip_extension <- function(
  filename = c("BR2_2016_S216_L001_R2_001.fastq", "foobar.fastq.gz", "fuzzbaz.txt"),
  ext = ".fastq"){

  .Deprecated("strip_extension", package = "charlier",
      msg = "Please use charlier::strip_extension")
      
  ix <- gregexpr(ext, filename, fixed = TRUE)
  sapply(seq_along(ix),
    function(i){
      if (ix[[i]] != -1) {
        s <- substring(filename[i], 1, ix[[i]]-1)
      } else {
        s <- filename[i]
      }
      s
    })
}


#' Add an extension to a filename
#'
#' @export
#' @param filename character, vector of one or more filenames
#' @param ext character, the extension to add
#' @param no_dup logical if TRUE then check first to determine if the filenames already
#'        have the extension, if they do have it do not add it again
#' @return the input filenames, possibly with the spcified extension added
add_extension <- function(
  filename = c("BR2_2016_S216_L001_R2_001.fastq", "foobar.fastq.gz", "fuzzbaz.txt"),
  ext = ".gz",
  no_dup = TRUE){

    .Deprecated("add_extension", package = "charlier",
        msg = "Please use charlier::add_extension")

  if (no_dup){
    pat <- paste0("^.*\\", ext[1], "$")
    ix <- grepl(pat, filename)
    filename[!ix] <- paste0(filename[!ix], ext[1])
  } else {
    filename <- paste0(filename, ext[1])
  }

  filename
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


#' List fastq files and separate into forward and reverse reads
#'
#' @export
#' @param path character, the input path
#' @param pattern_forward file pattern
#' @param pattern_reverse file pattern
#' @return named list of sorted forward and reverse fastq filenames
list_fastq <- function(path,
                       pattern_forward = "^.*_R1_001",
                       pattern_reverse = "^.*_R2_001"){
  .Deprecated('list_filepairs', package = 'dadautils',
    msg = "Please use list_filepairs() which is more versatile")
  list(
    forward = sort(list.files(path, pattern = pattern_forward, full.names = TRUE)),
    reverse = sort(list.files(path, pattern = pattern_reverse, full.names = TRUE)) )
}
