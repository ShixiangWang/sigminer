#' Simulation Analysis
#'
#' - `simulate_signature()` - Simulate signatures from signature pool.
#' - `simulate_catalogue()` - Simulate catalogs from signature/catalog pool.
#' - `simulate_catalogue_matrix()` - Simulate a bootstrapped catalog matrix.
#'
#' @name simulation
#' @param x a numeric vector representing a signature/catalog or matrix with rows representing
#' signatures/samples and columns representing components.
#' @param weights a numeric vector for weights.
#' @param n an integer indicating mutation number to be generated in a catalog.
#' @return a `matrix`.
#' @examples
#' # Generate a catalog
#' set.seed(1234)
#' catalog <- as.integer(table(sample(1:96, 1000, replace = TRUE)))
#' names(catalog) <- paste0("comp", 1:96)
#' # Generate a signature
#' sig <- catalog / sum(catalog)
#'
#' # Simulate catalogs
#' x1 <- simulate_catalogue(catalog, 10) # 10 mutations
#' x1
#' x2 <- simulate_catalogue(catalog, 100) # 100 mutations
#' x2
#' x3 <- simulate_catalogue(catalog, 1000) # 1000 mutations
#' x3
#' # Similar with a signature
#' x4 <- simulate_catalogue(sig, 10) # 10 mutations
#' x4
#'
#' # Load SBS signature
#' load(system.file("extdata", "toy_mutational_signature.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' s <- t(sig2$Signature.norm)
#' # Generate a signature from multiple signatures/catalogs
#' s1 <- simulate_signature(s)
#' s1
#' s2 <- simulate_signature(s, weights = 1:3)
#' s2
#' # Generate a catalog from multiple signatures/catalogs
#' c1 <- simulate_catalogue(s, 100, weights = 1:3)
#' c1
#' @testexamples
#' expect_is(x1, "integer")
#' expect_is(x2, "integer")
#' expect_is(x3, "integer")
#' expect_is(x4, "integer")
#' expect_is(s1, "numeric")
#' expect_is(s2, "numeric")
#' expect_is(c1, "integer")
NULL


#' @rdname simulation
#' @export
simulate_signature <- function(x, weights = NULL) {
  stopifnot(is.numeric(x))
  if (!requireNamespace("matrixStats")) {
    message("'matrixStats' package is not available, installing it.")
    eval(parse(text = "install.packages('matrixStats')"))
  }

  if (is.vector(x)) {
    mat <- matrix(x, nrow = 1, byrow = TRUE)
    colnames(mat) <- names(x)
  } else {
    mat <- x
  }
  mat2 <- t(mat[, , drop = FALSE])
  # Normalize all catalogue/signatures to 1
  mat2 <- apply(mat2, 2, function(x) x / sum(x, na.rm = TRUE))

  if (is.null(weights)) {
    weights <- rep(1L, ncol(mat2))
  }

  prob <- matrixStats::rowWeightedMeans(mat2, w = weights, na.rm = TRUE)
  names(prob) <- colnames(mat)
  prob
}

#' @rdname simulation
#' @export
simulate_catalogue <- function(x, n, weights = NULL) {
  force(n)
  prob <- simulate_signature(x, weights)
  out <- sample(seq_len(length(prob)), size = n, replace = TRUE, prob = prob)
  out <- as.integer(table(factor(out, levels = seq_along(prob))))
  names(out) <- names(prob)
  out
}

#' @rdname simulation
#' @export
simulate_catalogue_matrix <- function(x) {
  stopifnot(is.matrix(x))
  apply(x, 1, function(catalog) {
    simulate_catalogue(catalog, as.integer(sum(catalog)))
  }) %>% t()
}
