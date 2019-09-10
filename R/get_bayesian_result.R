get_bayesian_result = function(run_info) {
  necessary_cols = c("Run", "file")
  stopifnot(all(necessary_cols %in% colnames(run_info)))
  stopifnot(nrow(run_info) == 1)

  result = readRDS(run_info$file)

  W = result[[1]]
  H = result[[2]]
  index = colSums(W) > 1
  W = W[, index]
  H = H[index, ]
  K = sum(index)

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
