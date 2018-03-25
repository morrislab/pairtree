library(sciClone)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
files <- args[(1:(length(args)-2))]
resultsfn <- args[(length(args)-1)]
sampsfn <- args[length(args)]

cat(files)
sampnames <- lapply(files, basename)
sampnames <- gsub("\\.dat$", "", sampnames)
inputs <- lapply(files, read.table, header=TRUE)

sc <- sciClone(vafs=inputs, sampleNames=sampnames, useSexChrs=FALSE)
write(toJSON(sampnames), file=sampsfn)
writeClusterTable(sc, resultsfn)
