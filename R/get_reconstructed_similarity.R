#' Get Reconstructed Profile Cosine Similarity, RSS, etc.
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
  flag_rss <- TRUE
  has_cn <- any(grepl("^CN[^C]", colnames(nmf_matrix)) | startsWith(colnames(nmf_matrix), "copynumber"))

  raw_catalog <- t(nmf_matrix)
  if (has_cn) {
    # CN Method W
    if (!is.null(Signature$Raw$W)) {
      rec_catalog <- Signature$Raw$W %*% Signature$Raw$H
    } else {
      warning("Cannot calculate reconstructed profile without raw W and H for CN 'W'/'M' method.")
      return(invisible(NULL))
    }
  } else {
    rec_catalog <- Signature$Signature.norm %*% Signature$Exposure
  }
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

  if (flag_rss) {
    get_rss <- function(x, y) {
      sum((x - y)^2)
    }
    get_unexplained_variance <- function(x, y) {
      get_rss(x, y) / sum(x^2)
    }

    rss <- purrr::map2_dbl(
      as.data.frame(raw_catalog),
      as.data.frame(rec_catalog),
      get_rss
    )

    unexplained_variance <- purrr::map2_dbl(
      as.data.frame(raw_catalog),
      as.data.frame(rec_catalog),
      get_unexplained_variance
    )
  } else {
    rss <- NULL
    unexplained_variance <- NULL
  }

  data.table::data.table(
    sample = cross_samps,
    similarity = sim,
    rss = rss,
    unexplained_variance = unexplained_variance
  )
}
