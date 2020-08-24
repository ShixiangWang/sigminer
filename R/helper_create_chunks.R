helper_create_chunks <- function(x, y, n = 2L, ...) {
  n <- as.integer(n)
  v <- round(seq(x, y, length.out = n + 1L))
  x <- data.table::data.table(
    start = v[-length(v)],
    end = v[-1] - 1L,
    ...
  )
  x$end[length(x$end)] <- x$end[length(x$end)] + 1L
  x
}
