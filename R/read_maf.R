# Read MAF data -----------------------------------------------------------

#' @title Read MAF Files
#'
#' @description This function is a wrapper of [maftools::read.maf].
#' Useless options in [maftools::read.maf] are dropped here.
#' You can also use [maftools::read.maf] to read the data.
#'
#' @inheritParams maftools::read.maf
#' @examples
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools", mustWork = TRUE)
#' if (!require("R.utils")) {
#'   message("Please install 'R.utils' package firstly")
#' } else {
#'   laml <- read_maf(maf = laml.maf)
#'   laml
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
