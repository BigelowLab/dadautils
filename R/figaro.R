#' Retrieve default values for figaro arguments
#'
#' @export
#' @return table of 
#' \itemize{
#'  \item{name the verbose name of the argument}
#'  \item{shortname the shortname of the argumement}
#'  \item{value the default value}
#' }
figaro_defaults <- function(){
  dplyr::tibble(
    name = c(
      "python", 
      "app",    
      "inputDirectory",
      "outputDirectory",
      "outputFileName",
      "ampliconLength",
      "forwardPrimerLength",
      "reversePrimerLength",
      "minimumOverlap",
      "subsample",
      "percentile",
      "fileNamingStandard"),
    shortname = c("cmd", "app", "i", "o", "n", "a", "f", "r", "m", "s", "p", "F"),
    value = c(
      "python3",
      "/mnt/storage/data/edna/packages/figaro/figaro.py",
      getwd(),
      getwd(),
      "trimParameters.json",
      "450",
      "20",
      "20",
      "20",
      "-1",
      "83",
      "illumina")
  )
}

#' Run figaro 
#'
#' @export
#' @param x table of arguments - see \code{\link{figaro_defaults}}
#' @param verbose logical, if TRUE output the generated command string
#' @param ... other arguments for \code{\link[base]{system2}}
#' @return hmmm
run_figaro <- function(x, verbose = FALSE, ...){

  cmd <- (x %>%
    dplyr::filter(.data$shortname == 'cmd'))$value
  app <- (x %>%
    dplyr::filter(.data$shortname == 'app'))$value
    
  x <- x %>%
    dplyr::filter(!(.data$shortname %in% c("cmd", "app")))
  args <- c(app, sprintf("--%s %s", x$name, x$value))
  if (verbose){
    cmd_string <- sprintf("%s %s", cmd, paste(args, collapse = " "))
    cat("dadautils::run_figaro()\n")
    cat(cmd_string, "\n")
  }
  ok <- system2(cmd, args = args, ...)
  ok
}
