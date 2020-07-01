read_vcf <- function(vcfs, samples = NULL, keep_only_pass = TRUE, genome_build = c("hg19", "hg38"), verbose = TRUE) {

  vcfs_name <- vcfs
  message("Reading file(s) ", paste(vcfs, collapse = ", "))
  vcfs <- purrr::map(vcfs, ~data.table::fread(., skip = "#", select = c(1, 2, 4, 5, 7)))

  if (is.null(samples)) {
    names(vcfs) <- basename(vcfs_name)
  } else {
    if (length(samples) != length(vcfs_name)) {
      stop("Unequal files and samples!")
    }
    names(vcfs) <- samples
  }

  vcfs <- data.table::rbindlist(vcfs, use.names = FALSE, idcol = "sample")

  # required.fields = c('Hugo_Symbol', 'Chromosome', 'Start_Position', 'End_Position', 'Reference_Allele', 'Tumor_Seq_Allele2',
  #                     'Variant_Classification', 'Variant_Type', 'Tumor_Sample_Barcode')

  colnames(vcfs) <- c("Tumor_Sample_Barcode", "Chromosome", "Tumor_Sample_Barcode", "Reference_Allele", "Tumor_Seq_Allele2", "filter")

  if (keep_only_pass) {
    vcfs <- vcfs[vcfs$filter == "PASS", ]
  }
  vcfs$filter <- NULL

  ## Annotate gene symbol
  if (genome_build == "hg19") {
    gene_dt <- readRDS(system.file("extdata", "human_hg19_gene_info.rds", package = "sigminer", mustWork = TRUE))
  } else {
    gene_dt <- readRDS(system.file("extdata", "human_hg38_gene_info.rds", package = "sigminer", mustWork = TRUE))
  }

  # structure(list(Tumor_Sample_Barcode = c("201T.cave.annot.vcf",
  #                                         "201T.cave.annot.vcf", "201T.cave.annot.vcf", "201T.cave.annot.vcf",
  #                                         "201T.cave.annot.vcf", "201T.cave.annot.vcf"), Chromosome = c("1",
  #                                                                                                       "1", "1", "1", "1", "1"), Tumor_Sample_Barcode = c(1249680L,
  #                                                                                                                                                          1264136L, 1268403L, 1275688L, 1374508L, 1390967L), Reference_Allele = c("C",
  #                                                                                                                                                                                                                                  "A", "G", "G", "G", "C"), Tumor_Seq_Allele2 = c("T", "T", "A",
  #                                                                                                                                                                                                                                                                                  "A", "A", "A")), row.names = c(NA, -6L), class = c("data.table",
  #                                                                                                                                                                                                                                                                                                                                     "data.frame"), .internal.selfref = <pointer: 0x7f8a3f80cae0>)

}
