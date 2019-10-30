#' Get Sample Groups from NMF Decomposition Information
#'
#' One of key results from signature analysis is to cluster samples into different
#' groups. This function takes `Signature` object as input
#' and return the membership in each cluster.
#'
#' @param Signature a `Signature` object obtained either from [sig_extract] or [sig_auto_extract].
#' @param method cluster method, could be 'consensus' or 'samples'.
#' @param match_consensus only used when the `method` is 'consensus'.
#' If `TRUE`, the result will match order as shown in consensus map.
#' @return a `data.table` object
#' @import NMF cluster
#' @export
#' @examples
#' \donttest{
#' # Load copy number prepare object
#' load(system.file("extdata", "toy_copynumber_prepare.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Extract copy number signatures
#' library(NMF)
#' sig <- sig_extract(cn_prepare$nmf_matrix, 2, nrun = 10)
#'
#' # Methods 'consensus' and 'samples' are from NMF::predict()
#' get_groups(sig)
#' get_groups(sig, match_consensus = TRUE)
#' get_groups(sig, method = "samples")
#' }
#' @seealso [NMF::predict()]
get_groups <- function(Signature, method = "consensus", match_consensus = FALSE) {
  stopifnot(inherits(Signature, "Signature"))

  if (method == "consensus") {
    if (!"nmf_obj" %in% names(Signature$Raw)) {
      stop("Input Signature object does not contain NMF object, please select other methods")
    }
    nmfObj <- Signature$Raw$nmf_obj
    predict.consensus <- predict(nmfObj, what = "consensus")
    silhouette.consensus <- silhouette(nmfObj, what = "consensus")

    data <- data.frame(
      sample = sampleNames(nmfObj),
      nmf_submethods = predict.consensus,
      sil_width = signif(silhouette.consensus[, "sil_width"], 3)
    )
    colnames(data) <- c("sample", "group", "silhouette_width")
    # If we want to display as we see in consensusmap, we just need to reoder everything.
    # Now re-order data to match consensusmap sample order
    if (match_consensus) {
      sample.order <- attributes(predict.consensus)$iOrd
      data <- data[sample.order, ]
    }
  } else if (method == "samples") {
    if (!"nmf_obj" %in% names(Signature$Raw)) {
      stop("Input Signature object does not contain NMF object, please select other methods")
    }
    nmfObj <- Signature$Raw$nmf_obj
    predict.samples <- predict(nmfObj, what = "samples", prob = T)
    silhouette.samples <- silhouette(nmfObj, what = "samples")
    data <- data.frame(
      sample = names(predict.samples$predict),
      group = predict.samples$predict,
      silhouette_width = signif(silhouette.samples[, "sil_width"], 3),
      prob = signif(predict.samples$prob, 3)
    )
  } else {
    stop(paste("Wrong method:", method, "Possible options are: 'consensus', 'samples' "))
  }

  data <- data.table::as.data.table(data)
  return(data)
}
