#' Retrieve the edna storage path
#'
#' @export
#' @param ... arguments for \code{\link[base]{file.path}}
#' @param root character, the root path
edna_path <- function(..., root = "/mnt/storage/data/edna"){
  file.path(root, ...)
} 

#' Given a path - make it if it doesn't exist
#' 
#' @export
#' @param path character, the path to check and/or create
#' @param recursive logical, create paths recursively?
#' @param ... other arguments for \code{\link[base]{dir.create}}
#' @return logical, TRUE if the path exists or is created
make_path <- function(path, recursive = TRUE, ...){
  ok <- dir.exists(path[1])
  if (!ok){
    ok <- dir.create(path, recursive = recursive, ...)
  }
  ok
}