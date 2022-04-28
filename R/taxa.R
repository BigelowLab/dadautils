#' Verify taxa levels match depoth of taxa in reference database.
#'
#' Verification is by comparing count of taxa levels with the first line in reference
#' database which, in dada2-world, is a delimited list of levels.
#'
#' @export
#' @param taxLevels character, vector of desired taxa levels
#' @param refFasta character, filename of the reference database
#' @param delim character, the delimiter used to split the first contig description line in the refFasta
#' @return logical TRUE if taxLevels and refdb have the same depth of classifications
verify_taxalevels <- function(taxLevels, refFasta, delim = ";"){
  
  stopifnot(inherits(taxLevels, "character"))
  stopifnot(inherits(refFasta, "character"))
  
  
  if (!file.exists(refFasta[1])) stop("refFasta not found:", refFasta[1])
  s <- readLines(refFasta[1], n = 1)
  s <- strsplit(s, delim[1])[[1]]
  
  nt <- length(taxLevels)
  ns <- length(s)
  
  if (nt > ns){
    warning("number of taxLevels is greater than those provided in refFasta")
  }
  if (nt < ns){
    warning("number of taxLevels is less than than those provided in refFasta")
  }
  
  nt==ns
  #return(length(s) == length())
}


#' A wrapper around assign_taxonomy that optionally populates otherwise truncated 
#' taxonomy levels when running straight up \code{\link[dada2]{assignTaxonomy}}.
#' 
#' We compare the number of available levels in the reference database with the number of 
#'   levels in the user specified \code{taxLevels}.  We issue a warning if they are not the same but
#'   we do not throw an error.
#'
#' @export
#' @param seqs character, see \code{\link[dada2]{assignTaxonomy}}
#' @param refFasta filename, a the fully qualified filename to to the reference database
#' @param taxLevels character, \code{\link[dada2]{assignTaxonomy}}
#' @param populate_truncated logical, if TRUE make sure each level of \code{taxLevel} is 
#' represented in the returned value
#' @param truncated_value character or NA (default), value to assign to truncated levels.
#'   If \code{TRUE} then \code{minBoot} is ignored - see \code{\link[dada2]{assignTaxonomy}}
#' @param drop_levels character, NA or a vector of levels to drop. Typical to drop 'Species'
#'        if \code{\link[dada2]{addSpecies}} is to be added next. If \code{TRUE} then \code{minBoot} 
#'        is ignored - see \code{\link[dada2]{assignTaxonomy}}
#' @param save_file logical, if TRUE save the file as CSV. If TRUE then \code{minBoot} is ignored,
#'     see \code{\link[dada2]{assignTaxonomy}}
#' @param filename character, the file to write to if save_file is TRUE 
#' @param ... further arguments for \code{\link[dada2]{assignTaxonomy}}
#' @return a named list with two elements 
#' \itemize{
#'   \item{tax, character matrix}
#'   \item{boot, integer matrix, possibly filled with NA_integer_ if outputBootstaps is FALSE}
#' }
assign_taxonomy <- function(seqs, 
  refFasta = NULL,
  taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
  populate_truncated = TRUE,
  truncated_value = NA_character_,
  drop_levels = c(NA_character_, "Species")[1],
  save_file = FALSE,
  filename = "taxa.csv",
  ...){

    
    if (is.null(refFasta)) stop("refFasta is required")
      
    # here we chack that the levels are the same - we warn but do not throw an error
    ok <- verify_taxalevels(taxLevels, refFasta)
    
    x <- dada2::assignTaxonomy(seqs, refFasta, taxLevels = taxLevels, ...)  
    cat("\n") # because assign taxonomy is verbose without adding a new line
    
    if (!inherits(x, "list")){
      # if the user sets outputBootstraps = FALSE then we simply
      # cast a copy of x (a matrix at this point) to integer
      # this forces each element to NA_integer_, but nicely preserves the
      # matrix dimenions and row/col names
      boot <- x
      suppressWarnings(mode(boot) <- "integer")
      x <- list(tax = x, boot = boot)
    }
    
    if (populate_truncated){
      ix <- taxLevels %in% colnames(x$tax)
      if (!all(ix)) {
        ix <- taxLevels[!ix]
        y <- matrix(truncated_value, 
          ncol = length(ix), 
          nrow = nrow(x$tax),
          dimnames = list(rownames(x$tax), ix))
        x$tax <- cbind(x$tax,y)
      }
      
      if (!is.null(x$boot)){
        ix <- taxLevels %in% colnames(x$boot)
        if (!all(ix)) {
          ix <- taxLevels[!ix]
          y <- matrix(truncated_value, 
            ncol = length(ix), 
            nrow = nrow(x$boot),
            dimnames = list(rownames(x$boot), ix))
          x$boot <- cbind(x$boot,y)
        }
      }
    }  
    if (!is.na(drop_levels[1])){
      ix <- colnames(x$tax) %in% drop_levels
      x$tax <- x$tax[, !ix, drop = FALSE]
      if (!is.null(x$boot)){
        ix <- colnames(x$boot) %in% drop_levels
        x$boot <- x$boot[, !ix, drop = FALSE]
      }
    }
    if (save_file){
      if (!is.null(x$boot)){
        tax <- x$tax ; colnames(tax) <- paste("tax", colnames(tax), sep = ".")
        boot <- x$boot; colnames(boot) <- paste("boot", colnames(boot), sep = ".")
        y <- dplyr::as_tibble(tax, rownames = "seq") %>%
         dplyr::left_join(dplyr::as_tibble(boot, rownames = "seq"), by = "seq") %>%
         readr::write_csv(filename)
      } else {
        readr::write_csv(x$tax %>% dplyr::as_tibble(), filename)
      }
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

  tlevels <- taxa_levels[taxa_levels %in% colnames(x)]
  INDEX <- rep(TRUE, nrow(x))
  
  for (nm in names(vars)) {
    if (tolower(nm) == 'any'){
       idx <- sapply(tlevels,
          function(vn){
            x[[vn]] %in% vars[[nm]]
          })
        INDEX <- apply(idx, 1, any) & INDEX  
    } else {
      if (nm %in% tlevels) {
        INDEX <- INDEX & (x[[nm]] %in% vars[[nm]])
      } else {
        warning(sprintf("variable %s not found among columns of input", nm))
      }
    } 
  }
  
  if (index) {
    x <- !INDEX  
  } else {
    x <- dplyr::filter(x, !INDEX)
  }
  x
}
