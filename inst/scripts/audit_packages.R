# create package audit




main <- function(){

  orig_repos <- options("repos")
  options(repos = c(CRAN = "https://cloud.r-project.org" ))
  
  All <- dadautils::installed_packages() %>%
    dplyr::filter(grepl("^/mnt/modules/bin/dada2", LibPath)) %>%
    readr::write_csv("audit_all.csv")
    
  old <- dadautils::old_packages(checkBuilt = FALSE) %>%
    dplyr::filter(grepl("^/mnt/modules/bin/dada2", LibPath)) 
  
  # bioconductor
  bioc <- BiocManager::valid()[['out_of_date']] %>%
    dplyr::as_tibble() %>%
    dplyr::filter(grepl("^/mnt/modules/bin/dada2", LibPath)) %>% 
    readr::write_csv("audit_bioc.csv")
    
  Old <- old %>%
    dplyr::bind_rows(bioc) %>%
    dplyr::distinct()  
  
  # tidyverse
  tidy <- tidyverse::tidyverse_deps() %>%
    readr::write_csv("audit_tidy.csv")
  
  Old <- Old %>%
    #dplyr::filter(!(Package %in% tidy$package)) %>%
    readr::write_csv("audit_old.csv")


  # external
  
  ## cutadapt
  cutadapt_app <- Sys.which("cutadapt") 
  if (cutadapt_app != "") {
    cutadapt_version <- system2(cutadapt_app, args = "--version", stdout = TRUE)
  } else {
    cutadapt_version <- "unknown"
  }


  ## figaro
  #python_app <- Sys.which("python3")
  #if (python_app != "") {
  #  figaro_version <- system2(python_app, args = "/mnt/storage/data/edna/packages/figaro/figaro.py --version", stdout = TRUE)
  #} else {
  #  cutadapt_version <- "unknown"
  #}
  
  
  list(
    tidy = tidy,
    bioc = bioc,
    cutadapt = c(app = cutadapt_app, version = cutadapt_version),
    others = others)

}