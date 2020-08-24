#' Read VCF files as MAF object
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
#' maf <- read_vcf(vcfs)
#' maf <- read_vcf(vcfs, keep_only_pass = FALSE)
read_vcf <- function(vcfs, samples = NULL, genome_build = c("hg19", "hg38"), keep_only_pass = TRUE, verbose = TRUE) {
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
  if (genome_build == "hg19") {
    annot_file <- system.file("extdata", "human_hg19_gene_info.rds", package = "sigminer", mustWork = TRUE)
    gene_dt <- readRDS(annot_file)
  } else {
    annot_file <- system.file("extdata", "human_hg38_gene_info.rds", package = "sigminer", mustWork = TRUE)
    gene_dt <- readRDS(annot_file)
  }

  if (verbose) message("Annotating mutations to first matched gene based on database ", annot_file, "...")
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
