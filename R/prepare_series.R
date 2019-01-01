#----------------------------------------------
# Signature analysis prepare functions
#---------------------------------------------

#' Prepare nmf input matrix for copy number signature analysis
#' @inheritParams get_cnlist
#' @inheritParams get_features
#' @inheritParams get_components
#' @inheritParams get_matrix
#' @param reference_components default is `FALSE`, calculate mixture components
#' from [CopyNumber] object.
#' If set it to `NULL`, use pre-compiled components data which come from
#' [CNV signature paper](https://www.nature.com/articles/s41588-018-0179-8).
#' If can also be a `list` contain flexmix object of copy-number features,
#' obtain this from [get_components] function.
#' @return a `matrix` for NMF input
#' @export
#' @family signature analysis prepare function series
prepare_copynumber = function(CopyNumber,
                              reference_components = FALSE,
                              cores = 1, seed = 123456,
                              min_comp = 2, max_comp = 10,
                              min_prior = 0.001,
                              model_selection = "BIC",
                              nrep = 1, niter = 1000, rowIter = 1000) {
  stopifnot(is.logical(reference_components) | is.list(reference_components) | is.null(reference_components))

  cn_list = get_cnlist(CopyNumber)
  cn_features = get_features(CN_data = cn_list, cores = cores,
                             genome_build = CopyNumber@genome_build)

  if (is.null(reference_components) | is.list(reference_components)) {
    cn_components = reference_components
  } else {
    cn_components = get_components(CN_features = cn_features, seed = seed,
                                   min_comp = min_comp, max_comp = max_comp,
                                   min_prior = min_prior,
                                   model_selection = model_selection,
                                   nrep = nrep, niter = niter)
  }

  cn_matrix = get_matrix(cn_features, cn_components, cores = cores, rowIter = rowIter)
  cn_matrix
}


#' @inherit get_context
#' @export
#' @family signature analysis prepare function series
#' @examples
#' \dontrun{
#' laml.tnm <- prepare_maf(maf = laml, ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
#'     prefix = 'chr', add = TRUE, useSyn = TRUE)
#' }
prepare_maf = function(
  maf, ref_genome = NULL, prefix = NULL,
  add = TRUE, ignoreChr = NULL, useSyn = TRUE, fn = NULL
) {
  suppressMessages(get_context(
    maf, ref_genome = ref_genome, prefix = prefix,
    add = add, ignoreChr = ignoreChr, useSyn = useSyn, fn = fn
  ))
}

