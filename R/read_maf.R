# Read MAF data -----------------------------------------------------------

#' @title Read MAF Files
#'
#' @description This function is a wrapper of [maftools::read.maf].
#' Useless options in [maftools::read.maf] are dropped here.
#' You can also use [maftools::read.maf] to read the data.
#' All reference alleles and mutation alleles should be recorded in
#' positive strand format.
#'
#' @inheritParams maftools::read.maf
#' @examples
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools", mustWork = TRUE)
#' if (!require("R.utils")) {
#'   message("Please install 'R.utils' package firstly")
#' } else {
#'   laml <- read_maf(maf = laml.maf)
#'   laml
#'
#'   laml_mini <- laml@data[, list(
#'     Tumor_Sample_Barcode, Chromosome,
#'     Start_Position, End_Position,
#'     Reference_Allele, Tumor_Seq_Allele2
#'   )]
#'   laml2 <- read_maf_minimal(laml_mini)
#'   laml2
#' }
#' @testexamples
#' expect_s4_class(laml, "MAF")
#' @export
#' @seealso [read_copynumber] for reading copy number data to [CopyNumber] object.

read_maf <- function(maf, verbose = TRUE) {
  maftools::read.maf(
    maf,
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
    vc_nonSyn = NULL,
    verbose = verbose
  )
}

#' @describeIn read_maf Read Maf data.frame from a minimal maf-like data
#' @param dt A data.frame contains at least the following columns:
#' "Tumor_Sample_Barcode", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2"
#' @export
read_maf_minimal <- function(dt) {
  stopifnot(is.data.frame(dt), c(
    "Tumor_Sample_Barcode", "Chromosome",
    "Start_Position", "End_Position",
    "Reference_Allele", "Tumor_Seq_Allele2"
  ) %in% colnames(dt))

  if (is.null(dt$Variant_Type)) {
    dt$Variant_Type <- dplyr::case_when(
      nchar(dt$Reference_Allele) == 1L & nchar(dt$Tumor_Seq_Allele2) == 1L ~ "SNP",
      nchar(dt$Reference_Allele) < nchar(dt$Tumor_Seq_Allele2) ~ "INS",
      nchar(dt$Reference_Allele) > nchar(dt$Tumor_Seq_Allele2) ~ "DEL",
      nchar(dt$Reference_Allele) == 2L & nchar(dt$Tumor_Seq_Allele2) == 2L ~ "DNP",
      nchar(dt$Reference_Allele) == 3L & nchar(dt$Tumor_Seq_Allele2) == 3L ~ "TNP",
      TRUE ~ "Unknown"
    )
  }

  if (is.null(dt$Hugo_Symbol)) {
    dt$Hugo_Symbol <- "Unknown"
  }

  if (is.null(dt$Variant_Classification)) {
    dt$Variant_Classification <- "Unknown"
    set_vc <- TRUE
  } else {
    set_vc <- FALSE
  }

  message("Please note that the generated MAF object is designed for mutational signature analysis, not recommended for Maftools analysis!")
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
    vc_nonSyn = if (set_vc) "Unknown" else NULL,
    verbose = FALSE
  )
}
