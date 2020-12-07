#' Extract Signatures through NMF
#'
#' Do NMF de-composition and then extract signatures.
#'
#' @inheritParams sig_estimate
#' @inheritParams bp_extract_signatures
#' @param n_sig number of signature. Please run [sig_estimate] to select a suitable value.
#' @param optimize if `TRUE`, then refit the denovo signatures with QP method, see [sig_fit].
#' @param pynmf if `TRUE`, use Python NMF driver [Nimfa](http://nimfa.biolab.si/index.html).
#' The seed currently is not used by this implementation.
#' @param ... other arguments passed to [NMF::nmf()].
#' @author Shixiang Wang
#' @references Gaujoux, Renaud, and Cathal Seoighe. "A flexible R package for nonnegative matrix factorization." BMC bioinformatics 11.1 (2010): 367.
#' @references Mayakonda, Anand, et al. "Maftools: efficient and comprehensive analysis of somatic variants in cancer." Genome research 28.11 (2018): 1747-1756.
#' @return a `list` with `Signature` class.
#' @export
#' @examples
#' \donttest{
#' load(system.file("extdata", "toy_copynumber_tally_M.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Extract copy number signatures
#' res <- sig_extract(cn_tally_M$nmf_matrix, 2, nrun = 1)
#' }
#' @testexamples
#' expect_s3_class(res, "Signature")
#' @seealso [sig_tally] for getting variation matrix,
#' [sig_estimate] for estimating signature number for [sig_extract], [sig_auto_extract] for
#' extracting signatures using automatic relevance determination technique.
sig_extract <- function(nmf_matrix,
                        n_sig,
                        nrun = 10,
                        cores = 1,
                        method = "brunet",
                        optimize = FALSE,
                        pynmf = FALSE,
                        use_conda = TRUE,
                        py_path = "/Users/wsx/anaconda3/bin/python",
                        seed = 123456, ...) {
  eval(parse(text = "suppressMessages(library('NMF'))"))
  # transpose matrix
  mat <- t(nmf_matrix)

  ii <- colSums(mat) < 0.01
  if (any(ii)) {
    message(
      "The follow samples dropped due to null catalogue:\n\t",
      paste0(colnames(mat)[ii], collapse = ", ")
    )
    mat <- mat[, !ii, drop = FALSE]
  }

  if (isFALSE(pynmf)) {
    # To avoid error due to NMF
    mat <- check_nmf_matrix(mat)

    nmf.res <- NMF::nmf(
      mat,
      n_sig,
      seed = seed,
      nrun = nrun,
      method = method,
      .opt = paste0("vp", cores),
      ...
    )

    # Signature loading
    W <- NMF::basis(nmf.res)
    # Exposure loading
    H <- NMF::coef(nmf.res)
    # Signature number
    K <- ncol(W)
  } else {
    message("Calling python as backend...")
    env_install(use_conda, py_path, pkg = "nimfa", pkg_version = "1.4.0")
    reticulate::source_python(system.file("py", "nmf.py", package = "sigminer", mustWork = TRUE))

    result <- MultiNMF(mat, as.integer(n_sig), as.integer(nrun), as.integer(cores))
    W <- result$W
    H <- result$H
    K <- result$K

    rownames(W) <- rownames(mat)
    colnames(H) <- colnames(mat)
  }

  ## has_cn just used for method 'W' and 'M' in copy number signature
  has_cn <- grepl("^CN[^C]", rownames(W)) | startsWith(rownames(W), "copynumber")
  scal_res <- helper_scale_nmf_matrix(W, H, K, handle_cn = any(has_cn))
  Signature <- scal_res$Signature
  Exposure <- scal_res$Exposure

  if (optimize) {
    message("Refit the denovo signatures with QP.")
    ## Optimize signature exposure
    if (any(has_cn)) {
      mat_cn <- mat[has_cn, ]
      W_cn <- Signature[has_cn, ]
      W_cn <- apply(W_cn, 2, function(x) x / sum(x))

      ## Call LCD
      Exposure <- sig_fit(
        catalogue_matrix = mat_cn,
        sig = W_cn,
        method = "QP",
        mode = "copynumber"
      )
    } else {
      ## Call LCD
      Exposure <- sig_fit(
        catalogue_matrix = mat,
        sig = apply(Signature, 2, function(x) x / sum(x)),
        method = "QP"
      )
    }
  }

  # Handle hyper mutant samples
  Exposure <- collapse_hyper_records(Exposure)

  Signature.norm <- apply(Signature, 2, function(x) x / sum(x, na.rm = TRUE))
  Exposure.norm <- apply(Exposure, 2, function(x) x / sum(x, na.rm = TRUE))

  if (optimize) {
    ## Scale the result
    if (any(has_cn)) {
      Signature <- Signature.norm * sum(mat, na.rm = TRUE)
    } else {
      Signature <- Signature.norm * sum(Exposure, na.rm = TRUE)
    }
  }

  # When only one signature
  if (!is.matrix(Exposure.norm)) {
    Exposure.norm <- matrix(Exposure.norm, nrow = 1, dimnames = list(NULL, names(Exposure.norm)))
  }

  if (ncol(Signature) > 1) {
    # Get orders
    sig_orders <- helper_sort_signature(Signature.norm, Exposure)

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
      nmf_obj = if (exists("nmf.res")) nmf.res else NULL,
      W = W,
      H = H
    )
  )
  class(res) <- "Signature"
  attr(res, "nrun") <- nrun
  attr(res, "method") <- method
  attr(res, "seed") <- seed
  attr(res, "call_method") <- "NMF"

  invisible(res)
}
