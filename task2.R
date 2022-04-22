library(GenomicRanges)
source("./src/intervals.R")

filesList = c(
  "data/helas3_ctcf.narrowPeak.gz",
  "data/helas3_jun.narrowPeak.gz",
  "data/hepg2_ctcf.narrowPeak.gz",
  "data/hepg2_jun.narrowPeak.gz"
)

files = lapply(filesList, read.table)
granges = lapply(files, function(x) {
  GRanges(seqnames=x$V1, ranges=IRanges(x$V2, x$V3))
})

pairwiseJaccard(granges)