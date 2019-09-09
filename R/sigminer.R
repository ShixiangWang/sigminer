#' sigminer: Capture Genomic Variation Signatures using Non-Negative Matrix Factorization
#'
#' Genomic variations including somatic mutations, copy number variations, etc. are the major
#' force for cancer initialization and development. Due to the specificity of molecular lesions
#' caused by genomic variations, we can generate characteristic variation spectra, called 'signature'.
#' This package helps users to extract, analyze and visualize signatures from genomic variations,
#' thus providing new insight into cancer study.
#'
#' @section read data:
#' Read data into R objects:
#'
#' * [read_copynumber()] - read `data.frame`/files into [CopyNumber] object.
#' * [read_maf()] - read `data.frame`/file into [MAF] object. This is powered by
#' \code{\link[maftools]{read.maf}} function in **maftools** package.
#'
#' @section sig_ functions:
#' Signature analysis for genomic variations:
#'
#'
#' @section draw_ functions:
#' Result visualization for copy number data:
#'
#'
#' Result visualization for signature analysis:
#'
#'
#' Result visualization for [MAF] is provide by **maftools** package,
#' please read its [vignette](https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html).
#'
#'
#' @docType package
#' @importFrom dplyr progress_estimated
#' @import foreach
#' @import flexmix
#' @import doFuture
#' @import future
#' @name sigminer
NULL
