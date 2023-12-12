# Source https://github.com/marbl/CHM13
# https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/annotation/
# wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_combined3.gene_annotation_v0.1.gtf

# Chromosome size ---------------------------------------------------------

chromsize.T2T = GenomicRanges::seqinfo(BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0::BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0)
chromsize.T2T
chromsize.T2T = data.frame(chrom = paste0("chr", chromsize.T2T@seqnames), size = chromsize.T2T@seqlengths)

chromsize.T2T %>% str()

usethis::use_data(chromsize.T2T, overwrite = TRUE)


# Cytobands ---------------------------------------------------------------

cytobands.T2T <- data.table::fread("data-raw/chm13v2.0_cytobands_allchrs.bed", data.table = FALSE, header = FALSE)
head(cytobands.T2T)
colnames(cytobands.T2T) <- c("chrom", "start", "end", "band", "stain")

cytobands.T2T %>% str()

usethis::use_data(cytobands.T2T, overwrite = TRUE)


# Centromeres -------------------------------------------------------------

centromeres.hg38

centromeres.T2T <- data.table::fread("data-raw/chm13.draft_v2.0.cen_mask.bed", header = FALSE, data.table = FALSE)
centromeres.T2T
colnames(centromeres.T2T) = c("chrom", "left.base", "right.base")

centromeres.T2T %>% str()

usethis::use_data(centromeres.T2T, overwrite = TRUE)


# Transcript --------------------------------------------------------------

library(IRanges)

gtf_T2T <- data.table::fread("~/../Downloads/chm13v2.0_combined3.gene_annotation_v0.1.gtf", sep = "\t", header = FALSE)

head(gtf_T2T)
table(gtf_T2T$V7)
table(gtf_T2T$V1)

extract_col <- function(x, name) {
  library(magrittr)
  stringr::str_extract(x, paste0(name, " ([^;]+);")) %>%
    stringr::str_remove(paste0(name, " ")) %>%
    stringr::str_remove_all("\"") %>%
    stringr::str_remove(";")
}

#gtf_T2T[, gene_type := extract_col(V9, "gene_type")]

## Keep only protein coding region
gtf_T2T[, gene_id := extract_col(V9, "gene_id")]
gtf_T2T[, gene_name := extract_col(V9, "gene_name")]

hg38_gene = get_genome_annotation("gene", genome_build = "hg38")
hg38_gene
coding_ids = unique(hg38_gene$gene_name[hg38_gene$gene_type == "protein_coding"])

T2T <- gtf_T2T[V3 == "transcript" & gene_name %in% coding_ids, .(V1, V4, V5, V7)]
T2T
T2T <- T2T[, data.table::as.data.table(reduce(IRanges(V4, V5))), by = .(V7, V1)]
colnames(T2T)[1:2] <- c("strand", "chrom")
T2T
T2T$width <- NULL

T2T
## Save to package
transcript.T2T <- T2T
usethis::use_data(transcript.T2T, overwrite = TRUE)

# Gene --------------------------------------------------------------------

## T2T gene # No gene rows here, use merged transcript instead

gene_T2T2 <- gtf_T2T[, .(V1, V4, V5, V7, gene_name, gene_id)]
colnames(gene_T2T2)[1:4] <- c("chrom", "start", "end", "strand")
gene_T2T2
# t2t_gene_id = paste(gene_id, collapse = ",")
gene_T2T2 <- gene_T2T2[!is.na(gene_name), list(start = min(start), end = max(end)), by = .(chrom, strand, gene_name)]
gene_T2T2
gene_T2T2[5, 1:5]

# > gene_T2T2[gene_name == "TP53"]
# chrom strand gene_name     gene_id   start     end
# 1: chr17      -      TP53 XLOC_026526 7565929 7591642
# 2: chr17      -      TP53 XLOC_026527 7589364 7590475

hg38_gene
gene_T2T2_2 = merge(gene_T2T2, unique(hg38_gene[, list(chrom, strand, gene_name, gene_id, gene_type)]), by = c("chrom", "strand", "gene_name"), all.x = TRUE)
gene_T2T2_2

hg38_gene
data.table::setcolorder(gene_T2T2_2, colnames(hg38_gene))
## Save to extdata
saveRDS(gene_T2T2_2, file = "inst/extdata/human_T2T_gene_info.rds")

