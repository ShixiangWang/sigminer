library(tidyverse)
library(data.table)

## https://github.com/ShixiangWang/sigminer/issues/241

# http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes
# http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/cytoBand.txt.gz
# https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/gap.txt.gz

# Chromosome size ---------------------------------------------------------

download.file(
  "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes",
  "data-raw/mm10.chrom.sizes.txt"
)

chromsize.hg19 %>% str()

chromsize.mm10 = read.table("data-raw/mm10.chrom.sizes.txt", header = FALSE, sep = "\t")
colnames(chromsize.mm10) = c("chrom", "size")
chromsize.mm10 %>% str()

usethis::use_data(chromsize.mm10, overwrite = TRUE)


# Cytobands ---------------------------------------------------------------

download.file(
  "http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/cytoBand.txt.gz",
  "data-raw/mm10.cytoBand.txt.gz"
)

cytobands.hg19 %>% str()

cytobands.mm10 = data.table::fread("data-raw/mm10.cytoBand.txt.gz", data.table = FALSE)
colnames(cytobands.mm10) = c("chrom", "start", "end", "band", "stain")

cytobands.mm10 %>% str()

usethis::use_data(cytobands.mm10, overwrite = TRUE)


# Centromeres -------------------------------------------------------------

download.file(
  "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/gap.txt.gz",
  "data-raw/mm10.gap.txt.gz"
)

centromeres.hg19 %>% str()

centromeres.mm10 = data.table::fread("data-raw/mm10.gap.txt.gz")
centromeres.mm10[V8 == "centromere"]

## The problem is that centromeres do not exist in the mouse. UCSC lists them as 110000-3000000 for every chromosome.
centromeres.mm10 = data.frame(
  chrom = paste0("chr", c(1:19, "X", "Y")),
  left.base = 110000L,
  right.base = 3000000L,
  stringsAsFactors = FALSE
)

centromeres.mm10 %>% str()

usethis::use_data(centromeres.mm10, overwrite = TRUE)


# Transcript --------------------------------------------------------------

download.file(
  "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz",
  "data-raw/mm10.annotation.gtf.gz"
)

library(IRanges)

gtf_mm10 <- data.table::fread("data-raw/mm10.annotation.gtf.gz", skip = 5, sep = "\t", header = FALSE)

head(gtf_mm10)
table(gtf_mm10$V7)
table(gtf_mm10$V1)

extract_col <- function(x, name) {
  library(magrittr)
  stringr::str_extract(x, paste0(name, " ([^;]+);")) %>%
    stringr::str_remove(paste0(name, " ")) %>%
    stringr::str_remove_all("\"") %>%
    stringr::str_remove(";")
}

gtf_mm10[, gene_type := extract_col(V9, "gene_type")]

## Keep only protein coding region
mm10 <- gtf_mm10[V3 == "transcript" & gene_type == "protein_coding", .(V1, V4, V5, V7)]

mm10 <- mm10[, data.table::as.data.table(reduce(IRanges(V4, V5))), by = .(V7, V1)]
colnames(mm10)[1:2] <- c("strand", "chrom")
mm10$width <- NULL

## Save to package
transcript.mm10 <- mm10
usethis::use_data(transcript.mm10, overwrite = TRUE)

## Currently, I don't use gene location data, so don't generate it for now.
