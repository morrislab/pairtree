library(sciClone)

args <- commandArgs(trailingOnly = TRUE)
files <- args[-length(args)]
outfn <- args[length(args)]

cat(files)
sampnames <- lapply(files, basename)
sampnames <- gsub("\\.dat$", "", sampnames)
samples <- lapply(files, read.table, header=TRUE)

sc <- sciClone(vafs=samples, sampleNames=sampnames)
writeClusterTable(sc, outfn)
