# Copy from pryr package
mem_used <- function() {
  show_bytes(sum(gc()[, 1] * c(node_size(), 8)))
}

node_size <- function() {
  bit <- 8L * .Machine$sizeof.pointer
  if (!(bit == 32L || bit == 64L)) {
    stop("Unknown architecture", call. = FALSE)
  }

  if (bit == 32L) 28L else 56L
}

# mem_change <- function(code) {
#   start <- mem_used()
#
#   expr <- substitute(code)
#   eval(expr, parent.frame())
#   rm(code, expr)
#
#   show_bytes(mem_used() - start)
# }

show_bytes <- function(x) {
  structure(x, class = "bytes")
}

print.bytes <- function(x, digits = 3, ...) {
  power <- min(floor(log(abs(x), 1000)), 4)
  if (power < 1) {
    unit <- "B"
  } else {
    unit <- c("kB", "MB", "GB", "TB")[[power]]
    x <- x / (1000^power)
  }

  formatted <- format(signif(x, digits = digits),
    big.mark = ",",
    scientific = FALSE
  )

  cat(formatted, " ", unit, "\n", sep = "")
}
