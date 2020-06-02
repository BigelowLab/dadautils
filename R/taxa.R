#' A wrapper around assign_taxonomy that optionally populates otherwise truncated 
#' taxonomy levels when running straight up \code{\link[dada2]{assignTaxonomy}}.
#' 
#' @export
#' @param seqs character, see \code{\link[dada2]{assignTaxonomy}}
#' @param taxLevels character, \code{\link[dada2]{assignTaxonomy}}
#' @param populate_truncated logical, if TRUE make sure each level of \code{taxLevel} is 
#' represented in the returned value
#' @param truncated_value character or NA (default), value to assign to truncated levels
#' @param drop_levels character, NA or a vector of levels to drop. Typical to drop 'Species'
#'        if \code{\link[dada2]{addSpecies}} is to be added next.
#' @param ... further arguments for \code{\link[dada2]{assignTaxonomy}}
#' @return character matrix
assign_taxonomy <- function(seqs, 
  taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
  populate_truncated = TRUE,
  truncated_value = NA_character_,
  drop_levels = c(NA_character_, "Species")[1],
  ...){
  
    x <- dada2::assignTaxonomy(seqs, taxLevels = taxLevels, ...)  
    if (populate_truncated){
      ix <- taxLevels %in% colnames(x)
      if (!all(ix)) {
        ix <- taxLevels[!ix]
        y <- matrix(truncated_value, 
          ncol = length(ix), 
          nrow = nrow(x),
          dimnames = list(rownames(x), ix))
        x <- cbind(x,y)
      }
    }  
    if (!is.na(drop_levels[1])){
      ix <- colnames(x) %in% drop_levels
      x <- x[, !ix, drop = FALSE]
    }
    x
  }



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
  
  if (index) {
    x <- !INDEX  
  } else {
    x <- dplyr::slice(x, !INDEX)
  }
  x
}