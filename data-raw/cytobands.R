# From https://www.rdocumentation.org/packages/NetWeaver/versions/0.0.6/topics/ucsc.hg19.cytoband
load("data-raw/ucsc.hg19.cytoband.RData")
load("data-raw/ucsc.hg38.cytoband.RData")

cytobands.hg19 <- ucsc.hg19.cytoband
cytobands.hg38 <- ucsc.hg38.cytoband

colnames(cytobands.hg19) <- colnames(cytobands.hg38) <- c("chrom", "start", "end", "band", "stain")
colnames(cytobands.hg19)
colnames(cytobands.hg38)

usethis::use_data(cytobands.hg19)
usethis::use_data(cytobands.hg38)
