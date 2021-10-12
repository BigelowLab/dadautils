#' Run cutadapt
#'
#' @export
#' @param cut_files two element list for forward and reverse cutadapt results to be created
#' @param filt_files two element list of forward and reverse filtered files
#' @param CFG list of configuration
#' @param save_output logical, if TRUE try to capture the output of cutadapt to text files,
#'   one per cutadapt file generated, in the cutadapt output directory
#' @param save_graphics logical, if TRUE try to capture quality plots form the resulting cut_files
#' @return numeric codes, one per cut_file pairing as per \code{system2} where 0 means success
run_cutadapt <- function(
  cut_files,
  filt_files,
  CFG,
  save_output = TRUE,
  save_graphics = FALSE){
  
  
  if (is.numeric(CFG$multithread) && (CFG$multithread > 1)) {
    CFG$cutadapt$more_args <- paste0(CFG$cutadapt$more_args, " --cores=",CFG$multithread )
  }
  
  # here we split - either the user provided PRIMER info in the config (in which case we build the -a and -A args, OR the user drops these in favor of manually providing the these are part of CFG$more_args)
  if (is.null(CFG$primer) || charlier::is.nullna(CFG$primer$FWD) || charlier::is.nullna(CFG$primer$REV)){
    # then just use more_args
    OK <- sapply(seq_along(cut_files$forward),
      function(i){
        if (save_output){
          ofile <- paste0(charlier::strip_extension(cut_files$forward[i]), ".cutadapt_output.txt")
        } else {
          ofile <- ""
        }
        system2(CFG$cutadapt$app, 
          args = c(
            CFG$cutadapt$more_args, 
             "-o", cut_files$forward[i], 
             "-p", cut_files$reverse[i],
             filt_files$forward[i], 
             filt_files$reverse[i]),
           stdout = ofile)
      })
  } else {
    # here we build from CFG$primer
    FWD.RC <- dada2::rc(CFG$primer$FWD)
    REV.RC <- dada2::rc(CFG$primer$REV)
    R1.flags <- paste("-a"," ^", CFG$primer$FWD, "...", REV.RC, sep="")
    R2.flags <- paste("-A", " ^", CFG$primer$REV, "...", FWD.RC, sep="")
    OK <- sapply(seq_along(cut_files$forward),
      function(i){
        if (save_output){
          ofile <- paste0(charlier::strip_extension(cut_files$forward[i]), ".cutadapt_output.txt")
        } else {
          ofile <- ""
        }
        system2(CFG$cutadapt$app, 
          args = c(
            R1.flags, 
            R2.flags, 
            CFG$cutadapt$more_args, 
             "-o", cut_files$forward[i], 
             "-p", cut_files$reverse[i],
             filt_files$forward[i], 
             filt_files$reverse[i]),
           stdout = ofile)
      })
  }
    
  if (all(OK == 0) & save_graphics){
    ix <- seq_len(max(length(cut_files$forward), 2))
    ofile <- paste0(strip_extension(cut_files$forward[1]), ".cutadapt_quality.pdf")
    grDevices::pdf(ofile)
    try(dada2::plotQualityProfile(cut_files$forward[ix]) +  dada2::plotQualityProfile(cut_files$reverse[ix]))
    grDevices::dev.off()
  }
  OK
}
