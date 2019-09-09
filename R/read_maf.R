# Read MAF data -----------------------------------------------------------

#' @inherit maftools::read.maf
#' @examples
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read_maf(maf = laml.maf)
#' @export

read_maf <- function(
  maf, clinicalData = NULL, removeDuplicatedVariants = TRUE,
  useAll = TRUE, gisticAllLesionsFile = NULL,
  gisticAmpGenesFile = NULL, gisticDelGenesFile = NULL,
  gisticScoresFile = NULL, cnLevel = "all", cnTable = NULL,
  isTCGA = FALSE, vc_nonSyn = NULL, verbose = TRUE) {
  maftools::read.maf(
    maf,
    clinicalData = clinicalData,
    removeDuplicatedVariants = removeDuplicatedVariants,
    useAll = useAll,
    gisticAllLesionsFile = gisticAllLesionsFile,
    gisticAmpGenesFile = gisticAmpGenesFile,
    gisticDelGenesFile = gisticDelGenesFile,
    gisticScoresFile = gisticScoresFile,
    cnLevel = cnLevel, cnTable = cnTable,
    isTCGA = isTCGA, vc_nonSyn = vc_nonSyn,
    verbose = verbose
  )
}
