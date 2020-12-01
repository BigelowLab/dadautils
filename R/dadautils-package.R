#' dadautils-package
#'
#' @description Provides added tools for working with the [dada2](https://benjjneb.github.io/dada2/index.html) package in the eDNA project.
#' @docType package
#' @name dadautils-package
#' @aliases dadautils
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#' @importFrom grDevices dev.off pdf 
#' @importFrom utils installed.packages write.csv 
NULL

#' Deprecated functions in dadautils
#' 
#' These functions still work but should be avoided.
#' 
#' \itemize{
#'  \item \code{\link{list_fastq}}: Please use \code{\link{list_filepairs}}
#' }
#' 
#' @name dadautils-deprecated
NULL
