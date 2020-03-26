#' Get Specified Bayesian NMF Result from Run
#'
#' Sometimes, we may want to use or inspect specified run result from [sig_auto_extract].
#' This function is designed for this purpose.
#'
#' @param run_info a `data.frame` with 1 row and two necessary columns `Run` and `file`.
#'
#' @return a `list`.
#' @author Shixiang Wang
#' @export
#'
#' @examples
#' load(system.file("extdata", "toy_copynumber_tally_W.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#'
#' res <- sig_auto_extract(cn_tally_W$nmf_matrix, result_prefix = "Test_copynumber", nrun = 1)
#'
#' # All run info are stored in res$Raw$summary_run
#' # Obtain result of run 1
#' res_run1 <- get_bayesian_result(res$Raw$summary_run[1, ])
#' @testexamples
#' expect_s3_class(res, "Signature")
#' expect_equal(length(res_run1), 7L)
#'
get_bayesian_result <- function(run_info) {
  necessary_cols <- c("Run", "file")
  stopifnot(all(necessary_cols %in% colnames(run_info)))
  stopifnot(nrow(run_info) == 1)

  result <- readRDS(run_info$file)

  W <- result[[1]]
  H <- result[[2]]
  index <- colSums(W) > 1
  W <- W[, index, drop = FALSE]
  H <- H[index, , drop = FALSE]
  K <- sum(index)

  has_cn <- grepl("^CN[^C]", rownames(W)) | startsWith(rownames(W), "copynumber")
  scal_res <- helper_scale_nmf_matrix(W, H, K, handle_cn = any(has_cn))
  Signature <- scal_res$Signature
  Exposure <- scal_res$Exposure

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

  list(
    Signature = Signature,
    Signature.norm = Signature.norm,
    Exposure = Exposure,
    Exposure.norm = Exposure.norm,
    K = K,
    W = W,
    H = H
  )
}
