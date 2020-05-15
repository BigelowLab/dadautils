# [dadautils](https://github.com/BigelowLab/dadautils)

Provides utility functions for working with [dada2](https://benjjneb.github.io/dada2/index.html) package and the eDNA workflow.

### Requirements

+ [R v3.6+](https://www.r-project.org/)

+ [rlang](https://CRAN.R-project.org/package=rlang)

+ [dplyr](https://CRAN.R-project.org/package=dplyrg)

+ [readr](https://CRAN.R-project.org/package=readr)

+ [configr](https://CRAN.R-project.org/package=configr)

+ [dada2](https://CRAN.R-project.org/package=dada2)

+ [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)

+ [ShortRead](https://bioconductor.org/packages/release/bioc/html/ShortRead.html)

+ [IRanges](https://bioconductor.org/packages/release/bioc/html/IRanges.html)

+ [patchwork](https://CRAN.R-project.org/package=patchwork)

+ [ggplot2](https://CRAN.R-project.org/package=ggplot2)
	

### Note

All of the requirements are preinstalled in the charlie module `dada2`.  If you are **not** operating within that module then installation of the above (and dependencies) is required.

### Installation from [github](https://github.com)

You may need a [personal access token](https://github.com/settings/tokens).

```
remotes::install_github("BigelowLab/dadautils")
```

### Getting help

```
library(dadautils)
?dadautils
```