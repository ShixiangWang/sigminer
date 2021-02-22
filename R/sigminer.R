#' sigminer: Extract, Analyze and Visualize Signatures for Genomic Variations
#'
#' - Author: [Shixiang Wang](https://shixiangwang.github.io/home/) (<w_shixiang@163.com>)
#' - Please go to <https://shixiangwang.github.io/sigminer-doc/> for full vignette.
#' - Please go to <https://shixiangwang.github.io/sigminer/reference/index.html>
#' for organized documentation of functions and datasets.
#' - Result visualization for [MAF] is provide by **maftools** package,
#' please read its [vignette](https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html).
#'
#'
#' @docType package
#' @importFrom dplyr progress_estimated
#' @importFrom future plan
#' @importFrom furrr future_pmap
#' @importFrom graphics hist
#' @importFrom grDevices colorRampPalette pdf dev.off
#' @importFrom stats as.formula t.test fisher.test pairwise.table quantile setNames dist kmeans
#' @importFrom utils data download.file head packageDescription
#' @import NMF
#' @import ggplot2
#' @import cli
#' @name sigminer
NULL

## usethis namespace: start
#' @useDynLib sigminer, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL
