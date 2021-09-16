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


#' Write errors object
#'
#' @export
#' @param x errors object from \code{\link{learn_errors}}
#' @param filename the name of the file to write
#' @return the input errors object
write_errors <- function(x, filename = "learn_errors"){
  charlier::write_RDS(x, filename)
}

#' Read errors object
#'
#' @export
#' @param filename the name of the file to read
#' @return errors object
read_errors <- function(filename = "learn_errors"){
  charlier::read_RDS(filename)
}

#' Write quality_profile_pairs object
#'
#' @export
#' @param x output from \code{\link{quality_profile_pairs}}
#' @param filename the name of the file to write
#' @return the input errors object
write_QPP <- function(x, filename = "qpp"){
  charlier::write_RDS(x, filename)
}

#' Read quality_profile_pairs object
#'
#' @export
#' @param filename the name of the file to read
#' @return \code{\link{quality_profile_pairs}} listing
read_QPP <- function(filename = "qpp"){
  charlier::read_RDS(filename)
}

#' Write mergers object
#'
#' @export
#' @param x output from \code{\link{merge_pairs}}
#' @param filename the name of the file to write
#' @return the input errors object
write_mergers <- function(x, filename = "mergers"){
  charlier::write_RDS(x, filename)
}

#' Read mergers object
#'
#' @export
#' @param filename the name of the file to read
#' @return \code{\link{merge_pairs}} listing
read_mergers <- function(filename = "mergers"){
  charlier::read_RDS(filename)
}

#' Purge files froma folder bearing a name pattern
#'
#' @export
#' @param path character, the path to search
#' @param pattern character, pattern as a regex or if \code{glob = TRUE} then a glob pattern 
#' @param glob logical, if TRUE the pattern is a glob like "foo.*" and not regular expression
#' @return named logical vector where TRUE indicates a file was removed, possibly empty
purge_files <- function(path = ".", pattern = "^\\.Rsubread_qualityScores_score_pid", glob = FALSE){
  
  if (glob) pattern <- utils::glob2rx(pattern)
  ff <- list.files(path, pattern = pattern, full.names = TRUE)
  sapply(ff, file.remove) 
}
