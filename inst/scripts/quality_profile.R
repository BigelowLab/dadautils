library(rlang)
library(dplyr)
library(ggplot2)
library(charlier)
library(dadautils)
library(ShortRead)


# dashed orange = Q75
# solid orange = Q50
# dashed orange = Q75
# solid green th-ck = mean

PATH <- "/mnt/storage/data/edna/dada/projects/ben_foo/weeds-workshop"
files <- list_filepairs(file.path(PATH, "cutadapt"))

xx <- lapply(files, quality_profile_data)
x <- xx[[1]]
