#' Extract Signatures through NMF
#'
#' Do NMF de-composition and then extract signatures.
#'
#' @inheritParams sig_estimate
#' @param n_sig number of signature. Please run [sig_estimate] to select a suitable value.
#' @param optimize logical, used for optimizing copy number exposure.
#' NOTE: it can be only used when the `CN` feature is highly variable
#' across signatures and you should also check it by comparing to the result
#' from `optimize=FALSE`.
#' @param ... other arguments passed to [NMF::nmf()].
#' @author Shixiang Wang
#' @references Gaujoux, Renaud, and Cathal Seoighe. "A flexible R package for nonnegative matrix factorization." BMC bioinformatics 11.1 (2010): 367.
#' @references Mayakonda, Anand, et al. "Maftools: efficient and comprehensive analysis of somatic variants in cancer." Genome research 28.11 (2018): 1747-1756.
#' @return a `list` with `Signature` class.
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
#' @seealso [sig_tally] for getting variation matrix,
#' [sig_estimate] for estimating signature number for [sig_extract], [sig_auto_extract] for
#' extracting signatures using automatic relevance determination technique.
sig_extract <- function(nmf_matrix,
                        n_sig,
                        nrun = 10,
                        cores = 1,
                        method = "brunet",
                        optimize = FALSE,
                        pConstant = NULL,
                        seed = 123456, ...) {
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
    .opt = paste0("p", cores),
    ...
  )


  # Signature loading
  W <- NMF::basis(nmf.res)
  # Exposure loading
  H <- NMF::coef(nmf.res)
  # Signature number
  K <- ncol(W)

  has_cn <- grepl("^CN[^C]", rownames(W)) | startsWith(rownames(W), "copynumber")
  scal_res <- helper_scale_nmf_matrix(W, H, K, handle_cn = !optimize & any(has_cn))
  Signature <- scal_res$Signature
  Exposure <- scal_res$Exposure

  if (optimize) {
    ## Handle copy number signature exposure quantification
    ##  https://github.com/ShixiangWang/sigminer/issues/43#issuecomment-574039689
    if (any(has_cn)) {
      mat_cn <- mat[has_cn, ]
      nmf.res_cn <- NMF::nmf(
        mat_cn,
        K,
        seed = seed,
        nrun = nrun,
        method = method,
        .opt = paste0("p", cores),
        ...
      )
      ## Signature loading
      W_cn <- NMF::basis(nmf.res_cn)
      ## Exposure loading
      H_cn <- NMF::coef(nmf.res_cn)

      ## Match the signatures
      ## and assign the exposure
      scal_res_cn <- helper_scale_nmf_matrix(W_cn, H_cn, K, handle_cn = FALSE)

      to_sig <- apply(scal_res_cn$Signature, 2, function(x) x / sum(x, na.rm = TRUE))
      to_ref <- apply(Signature[has_cn,, drop=FALSE], 2, function(x) x / sum(x, na.rm = TRUE))
      #colnames(to_sig) <- colnames(to_ref) <- as.character(seq_len(ncol(to_sig)))
      to_match <- suppressMessages(get_sig_similarity(to_sig, to_ref))
      to_index <- apply(to_match$similarity, 1, which.max) %>% as.integer() %>% order()

      ## Take a check
      if (all(sort(to_index) == seq_len(ncol(to_sig)))) {
        Exposure <- scal_res_cn$Exposure[to_index, , drop = FALSE]
      } else {
        message("=> Copy number signature exposure optimized failedly!")
        message("=> Switch to the old way.")
        scal_res <- helper_scale_nmf_matrix(W, H, K, handle_cn = TRUE)
        Signature <- scal_res$Signature
        Exposure <- scal_res$Exposure
      }
    }
  }

  # Handle hyper mutant samples
  hyper_index <- grepl("_\\[hyper\\]_", colnames(Exposure))
  if (sum(hyper_index) > 0) {
    H.hyper <- Exposure[, hyper_index]
    H.nonhyper <- Exposure[, !hyper_index]
    sample.hyper <- sapply(
      colnames(H.hyper),
      function(x) strsplit(x, "_\\[hyper\\]_")[[1]][[1]]
    )
    unique.hyper <- unique(sample.hyper)
    n.hyper <- length(unique.hyper)
    x.hyper <- array(0, dim = c(nrow(H.hyper), n.hyper))
    for (i in 1:n.hyper) {
      x.hyper[, i] <- rowSums(H.hyper[, sample.hyper %in% unique.hyper[i]])
    }
    colnames(x.hyper) <- unique.hyper
    rownames(x.hyper) <- rownames(Exposure)
    Exposure <- cbind(H.nonhyper, x.hyper)
  }

  Signature.norm <- apply(Signature, 2, function(x) x / sum(x, na.rm = TRUE))
  Exposure.norm <- apply(Exposure, 2, function(x) x / sum(x, na.rm = TRUE))

  # When only one signature
  if (!is.matrix(Exposure.norm)) {
    Exposure.norm <- matrix(Exposure.norm, nrow = 1, dimnames = list(NULL, names(Exposure.norm)))
  }

  if (ncol(Signature) > 1) {
    # Get orders
    sig_orders <- helper_sort_signature(Signature.norm)

    Signature <- Signature[, sig_orders]
    Signature.norm <- Signature.norm[, sig_orders]
    Exposure <- Exposure[sig_orders, ]
    Exposure.norm <- Exposure.norm[sig_orders, ]

    W <- W[, sig_orders]
    H <- H[sig_orders, ]
  }

  sig_names <- paste0("Sig", seq_len(K))
  colnames(W) <- colnames(Signature) <- colnames(Signature.norm) <- sig_names
  rownames(H) <- rownames(Exposure) <- rownames(Exposure.norm) <- sig_names

  res <- list(
    Signature = Signature,
    Signature.norm = Signature.norm,
    Exposure = Exposure,
    Exposure.norm = Exposure.norm,
    K = K,
    Raw = list(
      nmf_obj = nmf.res,
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
