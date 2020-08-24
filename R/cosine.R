#' Calculate Cosine Measures
#'
#' @param x a numeric vector or matrix.
#' @param y must be same format as `x`.
#'
#' @return a numeric value or `matrix`.
#' @export
#'
#' @examples
#' x <- c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#' y <- c(0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0)
#' z1 <- cosine(x, y)
#' z1
#' z2 <- cosine(matrix(x), matrix(y))
#' z2
#' @testexamples
#' expect_equal(typeof(z1), typeof(1.0))
#' expect_equal(typeof(z1), typeof(matrix(1.0)))
cosine <- function(x, y) {
  stopifnot(!anyNA(x), !anyNA(y), is.numeric(x), is.numeric(y), is.vector(x) | is.matrix(x))
  if (is.vector(x)) {
    if (!is.vector(y)) {
      stop("y must be same type as x.")
    }
    cosineVector(x, y)
  } else {
    if (!is.matrix(y)) {
      stop("y must be same type as x.")
    }
    cosineMatrix(x, y)
  }
}
