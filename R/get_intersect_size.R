#' Get Overlap Size between Interval x and y
#'
#' @param x.start start position of interval x.
#' @param x.end start position of interval x.
#' @param y.start start position of interval x.
#' @param y.end start position of interval x.
#'
#' @return a numeric vector.
#' @export
#'
#' @examples
#' o1 <- get_intersect_size(1, 5, 3, 20)
#' o1
#' o2 <- get_intersect_size(3, 20, 1, 10)
#' o2
#' o3 <- get_intersect_size(c(1, 2, 1), c(10, 4, 6), c(4, 2, 5), c(10, 3, 22))
#' o3
#'
#' @testexamples
#' expect_equal(o1, 3)
#' expect_equal(o2, 8)
#' expect_equal(o3, c(7, 2, 2))
get_intersect_size <- function(x.start, x.end, y.start, y.end) {
  stopifnot(is.numeric(x.start), is.numeric(x.end), is.numeric(y.start), is.numeric(y.end))

  n <- length(x.start)
  if (n != length(x.end) || n != length(y.start) || n != length(y.end)) {
    stop("Input must be vector with same length.")
  }

  IntersectSize(x.start, x.end, y.start, y.end)
}
