#' sigminer: Extract, Analyze and Visualize Signatures for Genomic Variations
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
#' * [sig_tally()] - generate variation matrix for signature extraction.
#' * [sig_estimate()] - estimate signature number, work with [sig_extract()].
#' * [sig_extract()] - extract signatures throught NMF, a proper signature
#' number can be estimated by [sig_estimate()].
#' * [sig_auto_extract()] - extract signatures throught automatic relevance
#' determination technique.
#'
#' @section get_ functions:
#' Supplementary data clean and analysis functions:
#'
#' * [get_cn_ploidy()]
#' * [get_bayesian_result()]
#' * [get_groups()]
#' * [get_sig_exposure()]
#' * [get_sig_similarity()]
#' * [get_genome_annotation()]
#' * [get_adj_p()]
#' * [get_group_comparison()]
#' * [get_sig_feature_association()]
#' * [get_tidy_association()]
#' * [get_tidy_parameter()]
#' * [handle_hyper_mutation()]
#'
#'
#' @section show_ functions:
#'
#' Common visualization:
#'
#' * [show_sig_number_survey()]
#' * [show_sig_profile()]
#' * [show_sig_exposure()]
#' * [show_group_mapping()]
#' * [show_group_comparison()]
#' * [show_sig_feature_corrplot()]
#'
#'
#' Result visualization for copy number data:
#'
#' * [show_cn_distribution()]
#' * [show_cn_features()]
#' * [show_cn_components()]
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
#' @importFrom doParallel registerDoParallel
#' @importFrom graphics hist
#' @importFrom grDevices colorRampPalette
#' @importFrom stats as.formula fisher.test pairwise.table quantile setNames
#' @importFrom utils data download.file head packageDescription
#' @import flexmix
#' @import ggplot2
#' @import patchwork
#' @name sigminer
NULL
