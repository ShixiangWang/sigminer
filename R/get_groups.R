#' Get Sample Groups from Signature Decomposition Information
#'
#' One of key results from signature analysis is to cluster samples into different
#' groups. This function takes `Signature` object as input
#' and return the membership in each cluster.
#'
#' @param Signature a `Signature` object obtained either from [sig_extract] or [sig_auto_extract].
#' @param method grouping method, could be one of the following:
#' - 'exposure' - the default method, which assigns a sample into a group whose signature exposure
#' is dominant
#' - 'consensus' - returns the cluster membership based on the hierarchical clustering of the consensus matrix,
#' it can only be used for the result obtained by [sig_extract()] with multiple runs using **NMF** package.
#' - 'samples' - returns the cluster membership based on the contribution of signature to each sample
#' using **NMF** package.
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
#' get_groups(sig, method = "consensus", match_consensus = TRUE)
#' get_groups(sig, method = "samples")
#' }
#' @seealso [NMF::predict()]
get_groups <- function(Signature, method = c("exposure", "consensus", "samples"), match_consensus = FALSE) {
  stopifnot(inherits(Signature, "Signature"))
  method = match.arg(method)

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
      sil_width = signif(silhouette.consensus[, "sil_width"], 3),
      stringsAsFactors = FALSE
    )
    colnames(data) <- c("sample", "group", "silhouette_width")
    # If we want to display as we see in consensusmap, we just need to reoder everything.
    # Now re-order data to match consensusmap sample order
    if (match_consensus) {
      sample.order <- attributes(predict.consensus)$iOrd
      data <- data[sample.order, ]
    }

    data$group = as.character(data$group)
    data = find_enriched_signature(data, Signature)

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
      prob = signif(predict.samples$prob, 3),
      stringsAsFactors = FALSE
    )

    data$group = as.character(data$group)
    data = find_enriched_signature(data, Signature)

  } else if (method == "exposure") {
    expo_df = get_sig_exposure(Signature, type = "relative")
    data = expo_df %>%
      tidyr::gather(key = 'Signature', value = 'Exposure', dplyr::starts_with("Sig")) %>%
      dplyr::group_by(.data$sample) %>%
      dplyr::top_n(1, .data$Exposure) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        enrich_sig = .data$Signature
      ) %>%
      dplyr::rename(group = .data$Signature,
                    weight = .data$Exposure) %>%
      dplyr::mutate(group = sub("Sig", "", .data$group))
  }

  data <- data.table::as.data.table(data)
  return(data)
}
