# script to bin mappability bed file for mHiC ICEing
# Author: Benedikt G Brink, LMU Munich, 2019
library(GenomicRanges)

fai <- read.table(snakemake@input[["fai"]], stringsAsFactors = F)
bedfile <- read.table(snakemake@input[["bed"]])
binsize <- as.numeric(snakemake@params[["resolution"]])

bed_ranges <- GRanges()
ice_frame <- data.frame(chr = vector(), midPoint = vector(), anything = vector(), mappabilityValue = vector())
for (chr in levels(bedfile$V1)) {
  qwe <- bedfile[which(bedfile$V1 == chr),]
  asd <- GRanges(seqnames = qwe$V1, IRanges(start = qwe$V2, end = c(qwe$V2[-1], qwe$V3[length(qwe$V3)])), mappability = qwe$V5)
  chrlength <- fai[which(fai$V1 == chr),2]
  end(asd[length(asd)]) <- chrlength
  bed_ranges <- c(bed_ranges, asd)
}

bins <- slidingWindows(range(bed_ranges), binsize, binsize)
bins <- unlist(bins)
hits <- findOverlaps(bins, bed_ranges)
agg <- aggregate(bed_ranges, hits, score=mean(mappability))
bins$mappability <- agg$score

ice_frame <- data.frame(chr=seqnames(bins), 
                        midPoint=start(bins)+binsize/2, 
                        anything=rep(NA, length(bins)),
                        mappability=bins$mappability)

write.table(ice_frame, file = snakemake@output[[1]], sep = "\t", quote = FALSE, row.names = FALSE)