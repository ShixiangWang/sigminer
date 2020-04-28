#' sigminer: Extract, Analyze and Visualize Signatures for Genomic Variations
#'
#' Please go to <https://shixiangwang.github.io/sigminer-doc/> for full vignette.
#'
#' Result visualization for [MAF] is provide by **maftools** package,
#' please read its [vignette](https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html).
#'
#'
#' @docType package
#' @importFrom dplyr progress_estimated
#' @importFrom future plan
#' @importFrom furrr future_pmap
#' @importFrom foreach foreach %dopar%
#' @importFrom graphics hist
#' @importFrom grDevices colorRampPalette
#' @importFrom stats as.formula fisher.test pairwise.table quantile setNames
#' @importFrom utils data download.file head packageDescription
#' @import ggplot2
#' @import cli
#' @name sigminer
NULL
