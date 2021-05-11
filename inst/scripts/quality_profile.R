library(rlang)
library(dplyr)
library(ggplot2)
library(charlier)
library(dadautils)
library(ShortRead)
quality_profile_data
x <- quality_profile_data()

# dashed orange = Q75
# solid orange = Q50
# dashed orange = Q75
# solid green th-ck = mean

n = 500000
f <- system.file("extdata", "sam1F.fastq.gz", package="dada2")
srqa <- ShortRead::qa(f, n = n)

x <- quality_profile_data()
