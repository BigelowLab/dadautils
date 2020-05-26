#' Filter by removal of one or more instances of values found in one or more variables.
#'
#' @export
#' @param x a table as produced by \code{\link[dada2]{assignTaxonomy}} or \code{\link[dada2]{addSpecies}}
#' @param vars named list or NULL, with one character vector per variable (column) specifying 
#'        the instances to remove.  Names of the list must be present in the variable (column) 
#'        names of the input ar \code{any} in which case the values are searched for in any 
#'        variable.  It is not an error to request filtering of a non-existant column, but 
#'        a warning is generated. If NULL then return the input as a table.
#' @param taxa_levels character, a vector of the column names to search
#' @param index logical, if TRUE return a logical index into the original input where TRUE
#'        indicates rows to be removed.  If FALSE then return a filtered table.
#' @return a table with records (rows) removed or a logical index vector
taxa_remove <- function(x,
  vars = list(any = "Chloroplast"),
  taxa_levels = c("Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species"),
  index = FALSE){
  
  if (!inherits(x, "data.frame")) x <- dplyr::as_tibble(x)
  
  if (is.null(vars)) return(x)

  if (index) {
    INDEX <- rep(TRUE, nrow(x))
    for (nm in names(vars)) {
      if (tolower(nm) == 'any'){
         idx <- sapply(taxa_levels,
            function(vn){
              x[[vn]] %in% vars[[nm]]
            })
          INDEX <- apply(idx, 1, any) & INDEX  
      } else {
        if (nm %in% taxa_levels) {
          INDEX <- INDEX & (x[[nm]] %in% vars[[nm]])
        } else {
          warning(sprintf("variable %s not found among columns of input", nm))
        }
      } 
    }
  } else {
  
    # not a nice iteration - but each loop removes from column !!nm the rows matching vars[[nm]] 
    for (nm in names(vars)) {
      if (nm == 'any'){
        x <- x %>%
          dplyr::filter_all(any_vars(. %in% vars[[nm]]))
      } else {
        if (nm %in% taxa_levels) {
          x <- x %>%
            dplyr::filter(!(!!nm %in% vars[[nm]]))
        } else {
          warning(sprintf("variable %s not found among taxa_levels", nm))
        }
      }
    }
  }
  if (index){
    return(!INDEX)
  } else {
    return(x)
  }
}