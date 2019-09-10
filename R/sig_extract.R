#' Extract variation signatures
#'
#' Do NMF de-composition and then extract signatures.
#'
#' @inheritParams sig_estimate
#' @param n_sig number of signature. Please run [sig_estimate] to select a suitable value.
#' @author Shixiang Wang
#' @references Gaujoux, Renaud, and Cathal Seoighe. "A flexible R package for nonnegative matrix factorization." BMC bioinformatics 11.1 (2010): 367.
#' @references Mayakonda, Anand, et al. "Maftools: efficient and comprehensive analysis of somatic variants in cancer." Genome research 28.11 (2018): 1747-1756.
#' @return a `Signature` object.
#' @import NMF
#' @export
#' @examples
#' \donttest{
#' # Load copy number prepare object
#' load(system.file("extdata", "toy_copynumber_prepare.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Extract copy number signatures
#' library(NMF)
#' res <- sig_extract(cn_prepare$nmf_matrix, 2, nrun = 1)
#' }
sig_extract <- function(nmf_matrix,
                        n_sig,
                        nrun = 10,
                        cores = 1,
                        method = "brunet",
                        pConstant = NULL,
                        seed = 123456) {
  # transpose matrix
  mat <- t(nmf_matrix)

  # To avoid error due to non-conformable arrays
  if (!is.null(pConstant)) {
    if (pConstant < 0 | pConstant == 0) {
      stop("pConstant must be > 0")
    }
    mat <- mat + pConstant
  }

  nmf.res <- NMF::nmf(
    mat,
    n_sig,
    seed = seed,
    nrun = nrun,
    method = method,
    .opt = paste0("p", cores)
  )


  # Signature loading
  W <- NMF::basis(nmf.res)
  # Exposure loading
  H <- NMF::coef(nmf.res)

  K <- ncol(W)
  Signature <- W
  Exposure <- H

  # By scaling the signature loading matrix has all mutation burdens
  # - each signture (column in W) now represents a number of mutations
  # assigned to each signature.
  for (j in seq_len(K)) {
    Signature[, j] <- Signature[, j] * rowSums(H)[j]
    Exposure[j, ] <- Exposure[j, ] * colSums(W)[j]
  }

  Signature.norm <- apply(Signature, 2, function(x) x / sum(x, na.rm = TRUE))
  Exposure.norm <- apply(Exposure, 2, function(x) x / sum(x, na.rm = TRUE))

  sig_names <- paste0("Sig", seq_len(K))
  colnames(W) <- colnames(Signature) <- colnames(Signature.norm) <- sig_names
  rownames(H) <- rownames(Exposure) <- rownames(Exposure.norm) <- sig_names

  res <- list(
    Signature = Signature,
    Signature.norm = Signature.norm,
    Exposure = Exposure,
    Exposure.norm = Exposure.norm,
    K = K,
    raw = list(
      nmfObj = nmf.res,
      W = W,
      H = H
    )
  )
  class(res) <- "Signature"
  attr(res, "nrun") <- nrun
  attr(res, "method") <- method
  attr(res, "pConstant") <- pConstant
  attr(res, "seed") <- seed
  attr(res, "call_method") <- "NMF"

  invisible(res)
}
