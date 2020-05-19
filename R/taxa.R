#' Filter by removal of one or more instances of values found in one or more variables.
#'
#' @export
#' @param x a table as produced by \code{\link[dada2]{assignTaxonomy}} or \code{\link[dada2]{addSpecies}}
#' @param vars named list or NULL, with one character vector per variable (column) specifying the instances to remove.  Names of the list must be present in the variable (column) names of the input.  If NULL then return the input as a table.
#' @return a table with records (rows) removed
taxa_remove <- function(x,
  vars = list(Order = "Chloroplast") ){
  
  if (!inherits(x, "data.frame")) x <- dplyr::as_tibble(x)
  
  if (is.null(vars)) return(x)
  
  nms <- names(vars)
  
  ix <- nms %in% colnames(x)
  if (!all(ix)) stop("one or more variables to filter not found in input:", 
    paste(nms[!ix], collapse = " "))
  
  # not a nice iteration - but each loop removes from column !!nm the rows matching vars[[nm]] 
  for (nm in nms) {
    x <- x %>%
      dplyr::filter(!(!!nm %in% vars[[nm]]))
  }
  x
}