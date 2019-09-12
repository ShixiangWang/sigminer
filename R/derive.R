#' Derive Variation Matrix from Mutation Profile
#'
#' Generate a matrix for NMF de-composition and more. This is a generic function,
#' so it can be further extended to other mutation cases.
#'
#' @param object a [CopyNumber] object or [MAF] object.
#' @param ... custom setting for operating object. Detail see S3 method for
#' corresponding class (e.g. `CopyNumber`).
#' @return a `list` contains a `matrix` used for NMF de-composition.
#' @author Shixiang Wang
#' @export
#' @examples
#' \donttest{
#' # Load copy number object
#' load(system.file("extdata", "toy_copynumber.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Prepare copy number signature analysis
#' cn_prepare <- derive(cn)
#' }
derive <- function(object, ...) {
  UseMethod("derive")
}

#' @describeIn derive Derive copy number features, components and component-by-sample matrix
#' @param reference_components default is `FALSE`, calculate mixture components
#' from [CopyNumber] object.
#' @param cores number of compute cores to run this task.
#' You can use [future::availableCores()] function to check how
#' many cores you can use.
#' @param seed seed number.
#' @param min_comp minimal number of components to fit, default is 2.
#' Can also be a vector with length 6, which apply to each feature.
#' @param max_comp maximal number of components to fit, default is 10.
#' Can also be a vector with length 6, which apply to each feature.
#' @param min_prior the minimum relative size of components, default is 0.001.
#' Details about custom setting please refer to **flexmix** package.
#' @param model_selection model selection strategy, default is 'BIC'.
#' Details about custom setting please refer to **flexmix** package.
#' @param threshold default is `0.1`. Sometimes, the result components
#' include adjacent distributions with similar mu
#' (two and more distribution are very close), we use this threshold
#' to obtain a more meaningful fit with less components.
#' @param nrep number of run times for each value of component,
#' keep only the solution with maximum likelihood.
#' @param niter the maximum number of iterations.
#' @param keep_only_matrix if `TRUE`, keep only matrix for signature extraction.
#' @references Macintyre, Geoff, et al. "Copy number signatures and mutational
#' processes in ovarian carcinoma." Nature genetics 50.9 (2018): 1262.
#' @export
derive.CopyNumber <- function(object, reference_components = FALSE,
                              cores = 1, seed = 123456,
                              min_comp = 2, max_comp = 10,
                              min_prior = 0.001,
                              model_selection = "BIC",
                              threshold = 0.1,
                              nrep = 1, niter = 1000,
                              keep_only_matrix = FALSE,
                              ...) {
  stopifnot(is.logical(reference_components) | is.list(reference_components) | is.null(reference_components))

  cn_list <- get_cnlist(object)

  message("=> Step: getting copy number features")

  cn_features <- get_features(
    CN_data = cn_list, cores = cores,
    genome_build = object@genome_build
  )

  message("=> Step: fitting copy number components")
  if (is.null(reference_components) | is.list(reference_components)) {
    cn_components <- reference_components
  } else {
    cn_components <- get_components(
      CN_features = cn_features, seed = seed,
      min_comp = min_comp, max_comp = max_comp,
      min_prior = min_prior,
      model_selection = model_selection,
      threshold = threshold,
      nrep = nrep, niter = niter, cores = cores
    )
  }

  message("=> Step: calculating the sum of posterior probabilities")
  cn_matrix <- get_matrix(cn_features, cn_components, cores = cores)
  if (keep_only_matrix) {
    cn_matrix
  } else {
    list(features = cn_features,
         components = cn_components,
         parameters = get_tidy_parameter(cn_components),
         nmf_matrix = cn_matrix)
  }
}

#' @describeIn derive Derive SBS mutation features, components and component-by-sample matrix
#' @inheritParams maftools::trinucleotideMatrix
#' @references Mayakonda, Anand, et al. "Maftools: efficient and comprehensive analysis of somatic variants in cancer." Genome research 28.11 (2018): 1747-1756.
#' @examples
#' \donttest{
#' # Prepare mutational signature analysis
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read_maf(maf = laml.maf)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' mt_prepare <- derive(
#'   laml,
#'   ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
#'   prefix = "chr", add = TRUE, useSyn = TRUE
#' )
#' }
#' @export
derive.MAF <- function(object, ref_genome = NULL, prefix = NULL,
                       add = TRUE, ignoreChr = NULL, useSyn = TRUE,
                       keep_only_matrix = FALSE,
                       ...) {
  # // TODO: Rewrite this function instead of using maftools
  # // Make result consistent with result from derive.CopyNumber
  res <- maftools::trinucleotideMatrix(
    object,
    ref_genome = ref_genome, prefix = prefix,
    add = add, ignoreChr = ignoreChr, useSyn = useSyn, fn = NULL
  )

  if (keep_only_matrix) {
    res$nmf_matrix
  } else {
    res
  }
}
