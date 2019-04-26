#' sigminer: Capture Genomic Variation Signatures using Non-Negative Matrix Factorization
#'
#' sigminer provides functions for identification of
#' copy number signatures (Geoffrey et al. (2018) \doi{10.1038/s41588-018-0179-8}) and
#' mutational signatures (Alexandrov et al. (2013) \doi{10.1038/nature12477})
#' by non-negative matrix factorization, signature analysis and visualization.
#' It can be used to capture signatures of genomic variation, compare genotype or
#' phenotype features of different signatures and thus uncover the relationship
#' between the mechanism of genomic variation and phenotypes in cancer.
#'
#' @section read_ functions:
#' Read data into R objects:
#'
#' * [read_copynumber()] - read `data.frame`/files into [CopyNumber] object.
#' * [read_maf()] - read `data.frame`/file into [MAF] object. This is powered by
#' \code{\link[maftools]{read.maf}} function in **maftools** package.
#' * [read_variation()] - combine a [CopyNumber] object and a [MAF] object as
#' a [GenomicVariation] object. **It is not useful for now**.
#'
#' @section sig_ functions:
#' Signature analysis for genomic variations:
#'
#' * [sig_prepare()] - do preparation step of signature analysis for [CopyNumber]
#' or [MAF] object.
#' * [sig_estimate()] - provides survey plot and consensus map to user
#' for selecting the best signature number.
#' * [sig_extract()] - extract signatures based on specified signature number.
#' * [sig_assign_samples()] - assign samples to signatures. This classify samples
#' into different subgroups based on the dominant signature.
#' * [sig_get_activity()] - obtain signature activity in samples.
#' * [sig_get_correlation()] - obtain correlation matrix between signatures activity.
#' * [sig_get_similarity()] - get similarity between signatures.
#' * [sig_summarize_subtypes()] - get summary of signature subtypes.
#'
#' @section draw_ functions:
#' Result visualization for copy number data:
#'
#' * [draw_cn_distribution()] - plot copy number distribution either by length or chromosome.
#' * [draw_cn_features()] - plot copy number feature distribution.
#' * [draw_cn_components()] - plot mixture fit model components.
#'
#' Result visualization for signature analysis:
#'
#' * [draw_sig_profile()] - plot signature profile.
#' * [draw_sig_activity()] - plot signature activity.
#' * [draw_sig_corrplot()] - plot correlation between signature activities.
#' * [draw_subtypes_comparison()] - plot comparison between signature subtypes.
#'
#' Result visualization for [MAF] is provide by **maftools** package,
#' please read its [vignette](https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html).
#'
#'
#' @docType package
#' @name sigminer
NULL
