library(rlang)
library(dplyr)
library(ggplot2)
library(charlier)
#library(dadautils)
devtools::load_all()
library(ShortRead)


# dashed orange = Q75
# solid orange = Q50
# dashed orange = Q75
# solid green th-ck = mean

PATH <- "/mnt/storage/data/edna/dada/projects/ben_foo/weeds-workshop"
files <- list_filepairs(file.path(PATH, "cutadapt"))

xx <- lapply(files, quality_profile)
x <- xx[[1]]



sp <- SolexaPath(system.file('extdata', package='ShortRead'))
rfq <- readFastq(analysisPath(sp), pattern="s_1_sequence.txt")
sread(rfq)
id(rfq)
quality(rfq)