#' Get Reconstructed Profile Cosine Similarity
#'
#' See [bp_extract_signatures] for examples.
#'
#' @inheritParams sig_extract
#' @param Signature a `Signature` object.
#'
#' @return a `data.table`.
#' @export
get_sig_rec_similarity <- function(Signature, nmf_matrix) {
  stopifnot(inherits(Signature, "Signature"), is.matrix(nmf_matrix))
  raw_catalog <- t(nmf_matrix)
  rec_catalog <- Signature$Signature.norm %*% Signature$Exposure
  cross_samps <- intersect(colnames(raw_catalog), colnames(rec_catalog))
  cross_comps <- intersect(rownames(raw_catalog), rownames(rec_catalog))

  if (length(cross_samps) != ncol(rec_catalog)) {
    message("Filter out some samples because of no raw catalog available.")
  }

  raw_catalog <- raw_catalog[cross_comps, cross_samps, drop = FALSE]
  rec_catalog <- rec_catalog[cross_comps, cross_samps, drop = FALSE]

  sim <- purrr::map2_dbl(
    as.data.frame(raw_catalog),
    as.data.frame(rec_catalog),
    cosineVector
  )

  data.table::data.table(
    sample = cross_samps,
    similarity = sim
  )
}
