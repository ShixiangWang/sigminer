#' Read VCF Files as MAF Object
#'
#' @param vcfs VCF file paths.
#' @param samples sample names for VCF files.
#' @param genome_build genome build version like "hg19".
#' @param keep_only_pass if `TRUE`, keep only 'PASS' mutation for analysis.
#' @param verbose if `TRUE`, print extra info.
#'
#' @return a [MAF].
#' @export
#' @seealso [read_maf], [read_copynumber]
#' @examples
#' vcfs <- list.files(system.file("extdata", package = "sigminer"), "*.vcf", full.names = TRUE)
#' \donttest{
#' maf <- read_vcf(vcfs)
#' maf <- read_vcf(vcfs, keep_only_pass = FALSE)
#' }
#' @testexamples
#' expect_is(maf, "MAF")
read_vcf <- function(vcfs, samples = NULL, genome_build = c("hg19", "hg38", "mm10"), keep_only_pass = TRUE, verbose = TRUE) {
  genome_build <- match.arg(genome_build)
  vcfs_name <- vcfs
  if (verbose) message("Reading file(s): ", paste(vcfs, collapse = ", "))
  vcfs <- purrr::map(vcfs, ~ data.table::fread(., select = c(1, 2, 4, 5, 7)))

  if (is.null(samples)) {
    names(vcfs) <- file_name(vcfs_name, must_chop = ".vcf")
  } else {
    if (length(samples) != length(vcfs_name)) {
      stop("Unequal files and samples!")
    }
    names(vcfs) <- samples
  }

  vcfs <- data.table::rbindlist(vcfs, use.names = FALSE, idcol = "sample")
  colnames(vcfs) <- c("Tumor_Sample_Barcode", "Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2", "filter")
  if (!is.character(vcfs$Chromosome[1])) {
    vcfs$Chromosome <- as.character(vcfs$Chromosome)
  }

  if (keep_only_pass) {
    vcfs <- vcfs[vcfs$filter == "PASS", ]
    if (nrow(vcfs) < 1L) {
      stop("No mutation left after filtering.")
    }
  }
  vcfs$filter <- NULL

  vcfs$Chromosome <- ifelse(startsWith(vcfs$Chromosome, "chr"),
    vcfs$Chromosome,
    paste0("chr", vcfs$Chromosome)
  )
  vcfs$End_Position <- vcfs$Start_Position + pmax(nchar(vcfs$Reference_Allele), nchar(vcfs$Tumor_Seq_Allele2)) - 1L

  if (verbose) message("Annotating Variant Type...")
  vcfs$Variant_Type <- dplyr::case_when(
    nchar(vcfs$Reference_Allele) == 1L & nchar(vcfs$Tumor_Seq_Allele2) == 1L ~ "SNP",
    nchar(vcfs$Reference_Allele) < nchar(vcfs$Tumor_Seq_Allele2) ~ "INS",
    nchar(vcfs$Reference_Allele) > nchar(vcfs$Tumor_Seq_Allele2) ~ "DEL",
    nchar(vcfs$Reference_Allele) == 2L & nchar(vcfs$Tumor_Seq_Allele2) == 2L ~ "DNP",
    nchar(vcfs$Reference_Allele) == 3L & nchar(vcfs$Tumor_Seq_Allele2) == 3L ~ "TNP",
    TRUE ~ "Unknown"
  )

  vcfs$Variant_Classification <- "Unknown"
  vcfs$Hugo_Symbol <- "Unknown"

  # Annotate gene symbol
  gene_file <- switch(
    genome_build,
    mm10 = file.path(
      system.file("extdata", package = "sigminer"),
      "mouse_mm10_gene_info.rds"),
    file.path(
      system.file("extdata", package = "sigminer"),
      paste0("human_", genome_build, "_gene_info.rds")
    ))
  if (!file.exists(gene_file)) query_remote_data(basename(gene_file))
  gene_dt <- readRDS(gene_file)

  if (verbose) message("Annotating mutations to first matched gene based on database ", gene_file, "...")
  dt <- gene_dt[, c("chrom", "start", "end", "gene_name")]
  data.table::setkey(dt, "chrom", "start", "end")

  match_dt <- data.table::foverlaps(vcfs, dt,
    by.x = c("Chromosome", "Start_Position", "End_Position"),
    which = TRUE, nomatch = NULL
  )
  ## Keep only the first gene index
  match_dt <- match_dt[!duplicated(match_dt$xid)]
  vcfs$Hugo_Symbol[match_dt$xid] <- dt$gene_name[match_dt$yid]

  if (verbose) message("Transforming into a MAF object...")
  maftools::read.maf(
    vcfs,
    clinicalData = NULL,
    removeDuplicatedVariants = TRUE,
    useAll = TRUE,
    gisticAllLesionsFile = NULL,
    gisticAmpGenesFile = NULL,
    gisticDelGenesFile = NULL,
    gisticScoresFile = NULL,
    cnLevel = "all",
    cnTable = NULL,
    isTCGA = FALSE,
    vc_nonSyn = "Unknown",
    verbose = verbose
  )
}

#' Read UCSC Xena Variant Format Data as MAF Object
#'
#' @param path a path to variant file.
#'
#' @return a `MAF` object.
#' @export
#'
#' @examples
#' \donttest{
#' if (requireNamespace("UCSCXenaTools")) {
#'   library(UCSCXenaTools)
#'   options(use_hiplot = TRUE)
#'   example_file <- XenaGenerate(subset = XenaDatasets == "mc3/ACC_mc3.txt") %>%
#'     XenaQuery() %>%
#'     XenaDownload()
#'   x <- read_xena_variants(example_file$destfiles)
#'   x@data
#'   y <- sig_tally(x)
#'   y
#' }
#' }
#' @testexamples
#' if (requireNamespace("UCSCXenaTools")) {
#'   expect_is(y, "list")
#' }
read_xena_variants <- function(path) {
  dt <- data.table::fread(path)

  detect_name <- function(x, y, z) {
    if (x %in% z) x else y
  }

  data.table::setnames(
    dt,
    old = c(
      detect_name("Sample_ID", "sample", colnames(dt)),
      "gene",
      detect_name("chrom", "chr", colnames(dt)),
      "start", "end",
      detect_name("ref", "reference", colnames(dt)),
      "alt"
    ),
    new = c(
      "Tumor_Sample_Barcode", "Hugo_Symbol", "Chromosome",
      "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2"
    )
  )

  dt$Variant_Type <- dplyr::case_when(
    nchar(dt$Reference_Allele) == 1L & nchar(dt$Tumor_Seq_Allele2) == 1L ~ "SNP",
    nchar(dt$Reference_Allele) < nchar(dt$Tumor_Seq_Allele2) ~ "INS",
    nchar(dt$Reference_Allele) > nchar(dt$Tumor_Seq_Allele2) ~ "DEL",
    nchar(dt$Reference_Allele) == 2L & nchar(dt$Tumor_Seq_Allele2) == 2L ~ "DNP",
    nchar(dt$Reference_Allele) == 3L & nchar(dt$Tumor_Seq_Allele2) == 3L ~ "TNP",
    TRUE ~ "Unknown"
  )
  dt$Variant_Classification <- "Unknown"
  dt$Hugo_Symbol <- "Unknown"

  maftools::read.maf(
    dt,
    clinicalData = NULL,
    removeDuplicatedVariants = TRUE,
    useAll = TRUE,
    gisticAllLesionsFile = NULL,
    gisticAmpGenesFile = NULL,
    gisticDelGenesFile = NULL,
    gisticScoresFile = NULL,
    cnLevel = "all",
    cnTable = NULL,
    isTCGA = FALSE,
    vc_nonSyn = "Unknown",
    verbose = FALSE
  )
}
