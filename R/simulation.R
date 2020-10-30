# Simulate catalogue from catalogue/signature pool

# Simulate signatures from signature pool
simulate_signature <- function(x, weights = NULL) {
  stopifnot(is.numeric(x))
  if (is.vector(x)) {
    mat <- matrix(x, nrow = 1, byrow = TRUE)
    colnames(mat) <- names(x)
  } else {
    mat <- x
  }
  mat2 <- t(mat[, , drop = FALSE])
  # Normalize all catalogue/signatures to 1
  mat2 <- apply(mat2, 2, function(x) x / sum(x, na.rm = TRUE))

  if (!is.null(weights)) {
    weights <- rep(1L, ncol(mat2))
  }

  prob <- matrixStats::rowWeightedMeans(mat2, w = weights, na.rm = TRUE)
  names(prob) <- colnames(mat)
  prob
}


simulate_catalogue <- function(x, n, weights = NULL) {
  # x can be a numeric vector or matrix
  # (row represents multiple signatures/samples, col represents multiple catalogs)

  prob <- simulate_signature(x, weights)
  out <- sample(seq_len(length(prob)), size = n, replace = TRUE, prob = prob)
  out <- as.integer(table(out))
  names(out) <- names(prob)
  out
}

