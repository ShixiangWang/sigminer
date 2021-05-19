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

chromsize.mm10 <- read.table("data-raw/mm10.chrom.sizes.txt", header = FALSE, sep = "\t")
colnames(chromsize.mm10) <- c("chrom", "size")
chromsize.mm10 %>% str()

usethis::use_data(chromsize.mm10, overwrite = TRUE)


# Cytobands ---------------------------------------------------------------

download.file(
  "http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/cytoBand.txt.gz",
  "data-raw/mm10.cytoBand.txt.gz"
)

cytobands.hg19 %>% str()

cytobands.mm10 <- data.table::fread("data-raw/mm10.cytoBand.txt.gz", data.table = FALSE)
colnames(cytobands.mm10) <- c("chrom", "start", "end", "band", "stain")

cytobands.mm10 %>% str()

usethis::use_data(cytobands.mm10, overwrite = TRUE)


# Centromeres -------------------------------------------------------------

download.file(
  "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/gap.txt.gz",
  "data-raw/mm10.gap.txt.gz"
)

centromeres.hg19 %>% str()

centromeres.mm10 <- data.table::fread("data-raw/mm10.gap.txt.gz")
centromeres.mm10[V8 == "centromere"]

## The problem is that centromeres do not exist in the mouse. UCSC lists them as 110000-3000000 for every chromosome.
centromeres.mm10 <- data.frame(
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

# Gene --------------------------------------------------------------------

## mm10 gene
gtf_mm10 <- data.table::fread("data-raw/mm10.annotation.gtf.gz", skip = 5, sep = "\t", header = FALSE)

gtf_mm10[, `:=`(
  gene_name = extract_col(V9, "gene_name"),
  gene_id = extract_col(V9, "gene_id"),
  gene_type = extract_col(V9, "gene_type")
)]

gene_mm10 <- gtf_mm10[V3 == "gene", .(V1, V4, V5, V7, gene_name, gene_id, gene_type)]
colnames(gene_mm10)[1:4] <- c("chrom", "start", "end", "strand")

## Save to extdata
saveRDS(gene_mm10, file = "inst/extdata/mouse_mm10_gene_info.rds")



# ### ---------------------------------------------------------------------


# mm9 ---------------------------------------------------------------------


# Chromosome size ---------------------------------------------------------

download.file(
  "http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.chrom.sizes",
  "data-raw/mm9.chrom.sizes"
)


chromsize.mm9 <- read.table("data-raw/mm9.chrom.sizes", header = FALSE, sep = "\t")
colnames(chromsize.mm9) <- c("chrom", "size")
chromsize.mm9 %>% str()

usethis::use_data(chromsize.mm9, overwrite = TRUE)


# Cytobands ---------------------------------------------------------------

download.file(
  "http://hgdownload.cse.ucsc.edu/goldenpath/mm9/database/cytoBand.txt.gz",
  "data-raw/mm9.cytoBand.txt.gz"
)

# cytobands.hg19 %>% str()

cytobands.mm9 <- data.table::fread("data-raw/mm9.cytoBand.txt.gz", data.table = FALSE)
colnames(cytobands.mm9) <- c("chrom", "start", "end", "band", "stain")

cytobands.mm9 %>% str()

usethis::use_data(cytobands.mm9, overwrite = TRUE)


# Centromeres -------------------------------------------------------------
##shell
 # for i in $(seq 1 19) X Y;
 #  do
 #  "https://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/chr${i}_gap.txt.gz";
 #  done


centromeres.mm9 <- data.table::fread("data-raw/mm9.gap")
centromeres.mm9[V8 == "centromere"]

## The problem is that centromeres do not exist in the mouse. UCSC lists them as 0-3000000 for every chromosome.
centromeres.mm9 <- data.frame(
  chrom = paste0("chr", c(1:19, "X", "Y")),
  left.base = 0L,
  right.base = 3000000L,
  stringsAsFactors = FALSE
)

centromeres.mm9 %>% str()

usethis::use_data(centromeres.mm9, overwrite = TRUE)


# Transcript --------------------------------------------------------------

download.file(
  "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M1/gencode.vM1.annotation.gtf.gz",
  "data-raw/mm9.annotation.gtf.gz"
)

library(IRanges)

gtf_mm9 <- data.table::fread("data-raw/mm9.annotation.gtf.gz", skip = 5, sep = "\t", header = FALSE)

head(gtf_mm9)
table(gtf_mm9$V7)
table(gtf_mm9$V1)

extract_col <- function(x, name) {
  library(magrittr)
  stringr::str_extract(x, paste0(name, " ([^;]+);")) %>%
    stringr::str_remove(paste0(name, " ")) %>%
    stringr::str_remove_all("\"") %>%
    stringr::str_remove(";")
}

gtf_mm9[, gene_type := extract_col(V9, "gene_type")]

## Keep only protein coding region
mm9 <- gtf_mm9[V3 == "transcript" & gene_type == "protein_coding", .(V1, V4, V5, V7)]

mm9 <- mm9[, data.table::as.data.table(reduce(IRanges(V4, V5))), by = .(V7, V1)]
colnames(mm9)[1:2] <- c("strand", "chrom")
mm9$width <- NULL

## Save to package
transcript.mm9 <- mm9
usethis::use_data(transcript.mm9, overwrite = TRUE)

# Gene --------------------------------------------------------------------

## mm9 gene
gtf_mm9 <- data.table::fread("data-raw/mm9.annotation.gtf.gz", skip = 5, sep = "\t", header = FALSE)

gtf_mm9[, `:=`(
  gene_name = extract_col(V9, "gene_name"),
  gene_id = extract_col(V9, "gene_id"),
  gene_type = extract_col(V9, "gene_type")
)]

gene_mm9 <- gtf_mm9[V3 == "gene", .(V1, V4, V5, V7, gene_name, gene_id, gene_type)]
colnames(gene_mm9)[1:4] <- c("chrom", "start", "end", "strand")

## Save to extdata
saveRDS(gene_mm9, file = "inst/extdata/mouse_mm9_gene_info.rds")

