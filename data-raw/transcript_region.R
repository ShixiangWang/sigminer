## code to prepare `transcript_region` dataset goes here

download.file(
  "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gtf.gz",
  "data-raw/gencode.v33.hg38.annotation.gtf.gz"
)

download.file(
  "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh37_mapping/gencode.v33lift37.annotation.gtf.gz",
  "data-raw/gencode.v33.hg19.annotation.gtf.gz"
)

library(IRanges)

## hg38
gtf_hg38 <- data.table::fread("data-raw/gencode.v33.hg38.annotation.gtf.gz", skip = 5, sep = "\t", header = FALSE)
table(gtf_hg38$V7)
table(gtf_hg38$V1)

gtf_hg38 <- gtf_hg38[V3 == "transcript", .(V1, V4, V5, V7)]

transcript_hg38 <- gtf_hg38[, data.table::as.data.table(reduce(IRanges(V4, V5), min.gapwidth = 0L)), by = .(V7, V1)]
colnames(transcript_hg38)[1:2] <- c("strand", "chrom")
transcript_hg38$width <- NULL

## hg19
gtf_hg19 <- data.table::fread("data-raw/gencode.v33.hg19.annotation.gtf.gz", skip = 5, sep = "\t", header = FALSE)
table(gtf_hg19$V7)
table(gtf_hg19$V1)

gtf_hg19 <- gtf_hg19[V3 == "transcript", .(V1, V4, V5, V7)]

transcript_hg19 <- gtf_hg19[, data.table::as.data.table(reduce(IRanges(V4, V5), min.gapwidth = 0L)), by = .(V7, V1)]
colnames(transcript_hg19)[1:2] <- c("strand", "chrom")
transcript_hg19$width <- NULL

## Save to package
transcript.hg38 <- transcript_hg38
usethis::use_data(transcript.hg38, overwrite = TRUE)

transcript.hg19 <- transcript_hg19
usethis::use_data(transcript.hg19, overwrite = TRUE)


# Gene location -----------------------------------------------------------

## hg38
gtf_hg38 <- data.table::fread("data-raw/gencode.v33.hg38.annotation.gtf.gz", skip = 5, sep = "\t", header = FALSE)


extract_col <- function(x, name) {
  library(magrittr)
  stringr::str_extract(x, paste0(name, " ([^;]+);")) %>%
    stringr::str_remove(paste0(name, " ")) %>%
    stringr::str_remove_all("\"") %>%
    stringr::str_remove(";")
}

gtf_hg38[, `:=`(
  gene_name = extract_col(V9, "gene_name"),
  gene_id = extract_col(V9, "gene_id"),
  gene_type = extract_col(V9, "gene_type")
)]

gene_hg38 <- gtf_hg38[V3 == "gene", .(V1, V4, V5, V7, gene_name, gene_id, gene_type)]
colnames(gene_hg38)[1:4] <- c("chrom", "start", "end", "strand")

## hg19
gtf_hg19 <- data.table::fread("data-raw/gencode.v33.hg19.annotation.gtf.gz", skip = 5, sep = "\t", header = FALSE)

gtf_hg19[, `:=`(
  gene_name = extract_col(V9, "gene_name"),
  gene_id = extract_col(V9, "gene_id"),
  gene_type = extract_col(V9, "gene_type")
)]

gene_hg19 <- gtf_hg19[V3 == "gene", .(V1, V4, V5, V7, gene_name, gene_id, gene_type)]
colnames(gene_hg19)[1:4] <- c("chrom", "start", "end", "strand")

## Save to extdata
saveRDS(gene_hg19, file = "inst/extdata/human_hg19_gene_info.rds")
saveRDS(gene_hg38, file = "inst/extdata/human_hg38_gene_info.rds")

# fwrite(gene_hg19, file = "inst/extdata/human_hg19_gene_info.tsv.gz")
