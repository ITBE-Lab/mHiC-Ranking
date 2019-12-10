## This script takes the output of the mHiC pipeline and multiplies all interactions with the core by 0.5.
## This way, we account for the different ploidy in the genome assembly of T. brucei.
## Author: Benedikt G Brink, LMU Munich, 2019

interactions <- read.table(snakemake@input[[1]], sep = "\t")
ids <- grep("core", interactions[,1])
interactions[ids,2] <- floor(interactions[ids,2] * 0.5)
write.table(interactions, file = snakemake@output[[1]], sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
