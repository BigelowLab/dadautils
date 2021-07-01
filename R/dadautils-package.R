#' dadautils-package
#'
#' @description Provides added tools for working with the \href{https://benjjneb.github.io/dada2/index.html}{dada2} package in the eDNA project.
#' @docType package
#' @seealso \href{https://github.com/BigelowLab/charlier}{charlier} companion package
#' @name dadautils-package
#' @aliases dadautils
#' @importFrom rlang .data
#' @importFrom dplyr %>%
#' @importFrom grDevices dev.off pdf
#' @importFrom utils installed.packages write.csv
#' @importFrom stats setNames lm predict as.formula
NULL

#' Deprecated functions in dadautils
#'
#' These functions still work but should be avoided.
#'
#' \itemize{
#' \item \code{\link{add_extension}}: Please use \code{\link[charlier]{add_extension}} in \code{charlier} package
#' \item \code{\link{audit}}: Please use \code{\link[charlier]{audit}} in \code{charlier} package
#' \item \code{\link{autopopulate_config}}: Please use \code{\link[charlier]{autopopulate_config}} in \code{charlier} package. 
#' \item \code{link{count_cores}}: Please use \code{\link[charlier]{count_cores}} in \code{charlier} package 
#' \item \code{\link{get_pbs_jobid}}: Please use \code{\link[charlier]{get_pbs_jobid}} in \code{charlier} package
#' \item \code{\link{make_path}}: Please use \code{\link[charlier]{make_path}} in \code{charlier} package
#' \item \code{\link{list_fastq}}: Please use \code{\link{list_filepairs}}
#' \item \code{\link{logger_levels}}: Please use \code{\link[charlier]{logger_levels}} in \code{charlier} package
#' \item \code{\link{read_RDS}}: Please use \code{\link[charlier]{read_RDS}} in \code{charlier} package
#' \item \code{\link{strip_extension}}: Please use \code{\link[charlier]{strip_extension}} in \code{charlier} package
#' \item \code{\link{write_RDS}}: Please use \code{\link[charlier]{write_RDS}} in \code{charlier} package
#' }
#'
#' @name dadautils-deprecated
NULL
