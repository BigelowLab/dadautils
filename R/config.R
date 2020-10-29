#' Autopopulate one or more values within the config with values drawn from the config.
#'
#' The second argument, \code{fields}, identifies the variables that may be found as $PLACEHOLDERS within the config file. For example, if \code{fields = list(io = c("input_directory", "output_directory"))} then we first find those values in the config file, and then search for $IO_INPUT_DIRECTORY and $IO_OUTPUT_DIRECTORY placeholders throughout the config file and replace the specified fields values.
#'
#' @export
#' @param x configuration list
#' @param fields named list, provides the identifiers withion the config to be used as 
autopopulate_config <- function(
  x = read_config("/mnt/storage/data/edna/dada/edna-dada2/config/dada2_example.yml", 
                  autopopulate = FALSE),
  fields = list(global = 
      list(
        data_path = "/mnt/storage/data/edna/dada/projects/16S_Cyano_01",
        reference_path = "/mnt/storage/data/edna/dada/reference")
      )
    ){
    
      # anything to do?
    if (is.null(x$global) || is.null(fields)){ return(x) }
      
    # first get the $PLACEHOLDER=value vector
    v <- unlist(lapply(names(fields$global),
        function(xname){
          if (length(fields$global[[xname]]) == 1){
            keyvals <- x$global[[xname]]
            names(keyvals) <- toupper(xname)
          } else {
            keyvals <- sapply(fields$global[[xname]],
                        function(yname){
                          z <- x$global[[xname]][[yname]]
                          names(z) <- toupper(yname)
                          z
                        }, USE.NAMES = FALSE)
          names(keyvals) <- sprintf("%s_%s", toupper(xname), names(keyvals))
          }
          keyvals
        }))
        
    names(v) <- paste0("$GLOBAL_", names(v))

    tmp_file <- tempfile()
    yaml::write_yaml(x, tmp_file)
    txt <-  readLines(tmp_file)
    unlink(tmp_file)
    
    for (nm in names(v)){
      txt <- stringr::str_replace(txt, stringr::coll(nm), v[[nm]])
    }    

    x <- yaml::read_yaml(text = paste(txt, collapse = "\n"))

    x
}


#' Read the configuration 
#' 
#' @export
#' @param filename configuration filename
#' @param autopopulate logical, if TRUE then autopopulate the configuration
#' @param autodetect_cores logical, if TRUE and multithread is present (and "auto")
#'    then autoset the number of cores.
#' @return list config values or a try-error
read_config <- function(filename,
  autopopulate = TRUE,
  autodetect_cores = TRUE){
  
  cfg <- try(yaml::read_yaml(filename[1]))
  if (inherits(cfg, 'try-error')){
    print(cfg)
    stop("failed to read config file:", filename[1])
  }
  if (autopopulate) cfg <- autopopulate_config(cfg)
  if ("multithread" %in% names(cfg) && autodetect_cores){
      if (cfg$multithread[1] == 'auto') cfg$multithread <- count_cores()
  }
  cfg
}