#' Retrieve the edna storage path
#'
#' @export
#' @param ... arguments for \code{\link[base]{file.path}}
#' @param root character, the root path
edna_path <- function(..., root = "/mnt/storage/data/edna"){
  file.path(root, ...)
} 

