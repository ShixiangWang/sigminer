library(tidyverse)
library(data.table)

## https://asia.ensembl.org/Caenorhabditis_elegans/Location/Genome

extract_col <- function(x, name) {
  library(magrittr)
  stringr::str_extract(x, paste0(name, " ([^;]+);")) %>%
    stringr::str_remove(paste0(name, " ")) %>%
    stringr::str_remove_all("\"") %>%
    stringr::str_remove(";")
}

# Gene --------------------------------------------------------------------

gtf_ce11 <- data.table::fread("data-raw/Caenorhabditis_elegans.WBcel235.108.gtf.gz", skip = 5, sep = "\t", header = FALSE)

gtf_ce11[, `:=`(
  gene_name = extract_col(V9, "gene_name"),
  gene_id = extract_col(V9, "gene_id"),
  gene_type = extract_col(V9, "gene_biotype")
)]

gtf_ce11 <- gtf_ce11[V3 == "gene", .(V1, V4, V5, V7, gene_name, gene_id, gene_type)]
colnames(gtf_ce11)[1:4] <- c("chrom", "start", "end", "strand")
gtf_ce11[, chrom := paste0("chr", chrom)]

## Save to extdata
saveRDS(gtf_ce11, file = "inst/extdata/ce11_gene_info.rds")
