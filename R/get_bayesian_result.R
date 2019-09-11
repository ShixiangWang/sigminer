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
#' # Load copy number prepare object
#' load(system.file("extdata", "toy_copynumber_prepare.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' res <- sig_auto_extract(cn_prepare$nmf_matrix, result_prefix = "Test_copynumber", nrun = 1)
#'
#' # All run info are stored in res$Raw$summary_run
#' # Obtain result of run 1
#' res_run1 = get_bayesian_result(res$Raw$summary_run[1, ])
get_bayesian_result <- function(run_info) {
  necessary_cols <- c("Run", "file")
  stopifnot(all(necessary_cols %in% colnames(run_info)))
  stopifnot(nrow(run_info) == 1)

  result <- readRDS(run_info$file)

  W <- result[[1]]
  H <- result[[2]]
  index <- colSums(W) > 1
  W <- W[, index]
  H <- H[index, ]
  K <- sum(index)

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