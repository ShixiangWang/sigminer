#' Show Signature Consensus Map
#'
#' This function is a wrapper of `NMF::consensusmap()`.
#'
#' @inheritParams NMF::consensusmap
#' @param sig a `Signature` object obtained from [sig_extract].
#' @param lab_row labels for the rows.
#' @param lab_col labels for the columns.
#' @param ... other parameters passing to `NMF::consensusmap()`.
#'
#' @return nothing
#' @export
show_sig_consensusmap <- function(sig,
                                  main = "Consensus matrix",
                                  tracks = c("consensus:", "silhouette:"),
                                  lab_row = NA, lab_col = NA, ...) {
  stopifnot(inherits(sig, "Signature"))

  if (!is.null(sig$Raw$nmf_obj)) {
    NMF::consensusmap(sig$Raw$nmf_obj,
      main = "Consensus matrix",
      tracks = c("consensus:", "silhouette:"),
      labCol = lab_col, labRow = lab_row, ...
    )
  } else {
    stop("Unsupported 'Signature' object, it must be a result from sig_extract().")
  }
}
