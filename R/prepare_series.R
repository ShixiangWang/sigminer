#----------------------------------------------
# Signature analysis prepare functions
#---------------------------------------------


# Prepare nmf matrix for CopyNumber ---------------------------------------

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
#' @author Shixiang Wang
#' @return a `list` contains `matrix` for NMF input, copy number features and components.
#' @export
#' @examples
#' \donttest{
#' # Load copy number object
#' load(system.file("extdata", "toy_copynumber.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' cn_prepare <- prepare_copynumber(cn)
#' }
#' @family signature analysis prepare function series
prepare_copynumber <- function(CopyNumber,
                               reference_components = FALSE,
                               cores = 1, seed = 123456,
                               min_comp = 2, max_comp = 10,
                               min_prior = 0.001,
                               model_selection = "BIC",
                               nrep = 1, niter = 1000, rowIter = 1000) {
  stopifnot(is.logical(reference_components) | is.list(reference_components) | is.null(reference_components))

  cn_list <- get_cnlist(CopyNumber)
  cn_features <- get_features(
    CN_data = cn_list, cores = cores,
    genome_build = CopyNumber@genome_build
  )

  if (is.null(reference_components) | is.list(reference_components)) {
    cn_components <- reference_components
  } else {
    cn_components <- get_components(
      CN_features = cn_features, seed = seed,
      min_comp = min_comp, max_comp = max_comp,
      min_prior = min_prior,
      model_selection = model_selection,
      nrep = nrep, niter = niter
    )
  }

  cn_matrix <- get_matrix(cn_features, cn_components, cores = cores, rowIter = rowIter)
  list(features = cn_features, components = cn_components, nmf_matrix = cn_matrix)
}


# Prepare MAF: Get trinucleotide matrix ------------------------------------

#' Prepare nmf input matrix for mutational signature analysis
#'
#' NMF input matrix here is trinucletiode matrix. This function calls
#' `trinucleotideMatirx` provided by **maftools** to extract 96 mutation motifs.
#' @inherit maftools::trinucleotideMatrix
#' @family signature analysis prepare function series
#' @examples
#' \donttest{
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read_maf(maf = laml.maf)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' laml.tnm <- prepare_maf(
#'   maf = laml, ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
#'   prefix = "chr", add = TRUE, useSyn = TRUE
#' )
#' }
prepare_maf <- function(
                        maf, ref_genome = NULL, prefix = NULL,
                        add = TRUE, ignoreChr = NULL, useSyn = TRUE, fn = NULL) {
  maftools::trinucleotideMatrix(
    maf,
    ref_genome = ref_genome, prefix = prefix,
    add = add, ignoreChr = ignoreChr, useSyn = useSyn, fn = fn
  )
}
