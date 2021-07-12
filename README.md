# [dadautils](https://github.com/BigelowLab/dadautils)

Provides utility functions for working with [dada2](https://benjjneb.github.io/dada2/index.html) package and the eDNA workflow.

### Requirements

+ [R v4+](https://www.r-project.org/)

+ [rlang](https://CRAN.R-project.org/package=rlang)

+ [dplyr](https://CRAN.R-project.org/package=dplyr)

+ [tidyselect](https://CRAN.R-project.org/package=tidyselect)

+ [stringr](https://CRAN.R-project.org/package=stringr)

+ [readr](https://CRAN.R-project.org/package=readr)

+ [yaml](https://CRAN.R-project.org/package=yaml)

+ [dada2](https://CRAN.R-project.org/package=dada2)

+ [phyloseq](https://bioconductor.org/packages/release/bioc/html/phyloseq.html)

+ [Rsubread](https://bioconductor.org/packages/release/bioc/html/Rsubread.html)

+ [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)

+ [ShortRead](https://bioconductor.org/packages/release/bioc/html/ShortRead.html)

+ [IRanges](https://bioconductor.org/packages/release/bioc/html/IRanges.html)

+ [patchwork](https://CRAN.R-project.org/package=patchwork)

+ [ggplot2](https://CRAN.R-project.org/package=ggplot2)

 + [charlier](https://github.com/BigelowLab/charlier)
	

### Note

All of the requirements are preinstalled in the charlie module `dada2`.  If you are **not** operating within that module then installation of the above (including dependencies) is required.

### Note

This package is [dada2-centric](https://benjjneb.github.io/dada2/index.html) , but you may find the [charlier R package](https://github.com/BigelowLab/charlier) useful for general usage on charlie with in concert with `dadautils`.

### Installation from [github](https://github.com)

You may need a [personal access token](https://github.com/settings/tokens) if you use github.  We set the argument `upgrade = FALSE` so that the installer won't try to install packages in you home library that are already part of the module (but possibly as older versions.)

```
# if you have access to https://github.com/BigelowLab/dadautils
remotes::install_github("BigelowLab/dadautils", upgrade = FALSE)

# or 

# if you have access to /mnt/storage/data/edna/packages/dadautils
devtools::install("/mnt/storage/data/edna/packages/dadautils", upgrade = FALSE)
```

### Establishing a PBS environment and loading modules.

Working within a PBS session is not required, but operationally that is how the eDNA dada2 workflow is deployed (interactive or automated)

```
btupper@cfe1 ~ $ qsub -I -q devel -l walltime=8:00:00 -l ncpus=8,mem=8gb -N ben-edna
qsub: waiting for job 377353.cfe1 to start
qsub: job 377353.cfe1 ready
btupper@c3-6 ~ $
btupper@c3-6 ~ $ module use /mod/bigelow
btupper@c3-6 ~ $ module load dada2
```

### Getting started

Now we can invoke R and load the `dadautils` package from the library.

```
btupper@c3-6 ~ $ R

R version 3.6.2 (2019-12-12) -- "Dark and Stormy Night"
  .
  .
  .
  
library(dadautils)
```

### Getting help

```
?dadautils
```

### Piped workflows with [%>%](https://magrittr.tidyverse.org/)

We try to adopt functions that allow for piped workflows when writing scripts. Using the [magrittr](https://magrittr.tidyverse.org/) pipe operator allows us to pass a data object *implicity* from one function to the next.  That in turn permits the chaining of functions where the output of one function is *implicitly* passed as the first argument to the next function in line.  Think of a bucket brigade where the bucket is implied.

```
new_data <- start_with_this_data %>%
    do_this_to_it() %>%
    then_do_this_to_it() %>%
    and_finally_do_this_one_too()
```

Scripts are not required to use this paradigm, but it is an option (and it makes for tidy code.)


### Miscellaneous utilities


#### Build a path specification

```
edna_path("foo", "bar.csv")
# [1] "/mnt/storage/data/edna/foo/bar.csv"

edna_path("foo", "bar.csv", root = "/home/btupper")
# [1] "/home/btupper/foo/bar.csv"
```

#### File pair listings

Input files are often in pairs (forward and reverse) that are distinguished by naming pattern. This function will list those found in a directory and return a two element list.  These are very convenient and the list is used as a standard input to `dadautils` wrapper functions where appropriate. Some datasets do not have paired files (PacBio for instance), in this case we proceed with a filepair listing approach where one of the filepair elements is empty.

```
fq_files <- list_filepairs("/home/btupper/edna/data/examples/ben_demo_raw/filtN",
  pattern_forward = "*_R1_001.fastq",
  pattern_reverse = "*_R2_001.fastq") 
str(fq_files)
# List of 2
#  $ forward: chr [1:2] 
#     "/home/btupper/edna/data/examples/ben_demo_raw/filtN/BR2_2016_S216_L001_R1_001.fastq"
#     "/home/btupper/edna/data/examples/ben_demo_raw/filtN/BR3_2016_S228_L001_R1_001.fastq"
#  $ reverse: chr [1:2] 
#     "/home/btupper/edna/data/examples/ben_demo_raw/filtN/BR2_2016_S216_L001_R2_001.fastq" 
#     "/home/btupper/edna/data/examples/ben_demo_raw/filtN/BR3_2016_S228_L001_R2_001.fastq"
```

### Wrapper functions

Wrapper functions are intended to make convenient access to other commonly used functions.  Often a wrapper function will set a default value for input arguments, but the user can always override the defaults.  Sometimes wrapper functions will perform ancillary tasks like make and save a graphic, save output files, or the like.

#### Wrappers around functions external to R 

Some functions are external to the R session.  Typically, these are called from the command line with input and output files specified.  We use a wrapper function to hide all of the details form the user.  For example, the `run_cutadapt` function calls [cutadapt](https://cutadapt.readthedocs.io/en/stable/) command.  It optionally will output the `dada2::plotQualityProfile` plots in PDF format.

#### Wrappers for dada2

We have wrapped `dada2::filterAndTrim`,  `dada2::learnErrors`,  `dada2::dada`, and `dada2::mergePairs` to run more simply in a workflow.  We have renamed them following the `tidy_functions_use_underscores` principle and to keep the identities clear. As a convenience, each accepts a paired list of filenames (foreward and reverse) as created above using `list_filepairs`.

See `dadautils::filter_and_trim`, `dadautils::learn_errors`, `dadautils::run_dada` and `dadautils::merge_pairs`.  Each provides conveniences designed to make building workflows simple.