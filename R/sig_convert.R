#' Convert Signatures between different Genomic Distribution of Components
#'
#' Converts signatures between two representations relative to different sets of mutational opportunities.
#' Currently, only SBS signature is supported.
#'
#' The default opportunity matrix for "human-genome" and "human-exome" comes from COSMIC
#' signature database v2 and v3.
#'
#' @inheritParams sig_fit
#' @param from either one of "human-genome" and "human-exome" or an opportunity matrix
#' (repeated `n` columns with each row represents the total number of mutations for
#' a component, `n` is the number of signature).
#' @param to same as `from`.
#' @references `convert_signatures` function from sigfit package.
#'
#' @return a `matrix`.
#' @export
#' @examples
#' # Load SBS signature
#' load(system.file("extdata", "toy_mutational_signature.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Exome-relative to Genome-relative
#' sig_converted <- sig_convert(sig2,
#'   from = "human-exome",
#'   to = "human-genome"
#' )
#' sig_converted
#'
#' show_sig_profile(sig2, style = "cosmic")
#' show_sig_profile(sig_converted, style = "cosmic")
#' @testexamples
#' expect_is(sig_converted, "matrix")
#'
#' sig_converted2 <- sig_convert(sig2, from = "human-exome", to = "human-exome")
#' expect_true(all.equal(sig2$Signature.norm, sig_converted2))
#'
#' sig_converted2 <- sig_convert(sig2, from = "human-genome", to = "human-genome")
#' expect_true(all.equal(sig2$Signature.norm, sig_converted2))
sig_convert <- function(sig, from = "human-genome", to = "human-exome") {
  stopifnot(!is.null(from), !is.null(to))

  if (inherits(sig, "Signature")) {
    sig_mat <- sig$Signature.norm
  } else if (is.matrix(sig)) {
    sig_mat <- apply(sig, 2, function(x) x / sum(x))
  } else {
    stop("Invalid input for 'sig'", call. = FALSE)
  }

  from <- build_opps_matrix(ncol(sig_mat), nrow(sig_mat), from)
  to <- build_opps_matrix(ncol(sig_mat), nrow(sig_mat), to)

  conv_sigs <- sig_mat %>% t()

  if (any(ncol(from) != ncol(conv_sigs), ncol(to) != ncol(conv_sigs))) {
    stop("Matrix should have same number of components!")
  } else if (!is.null(colnames(conv_sigs))) {
    ## Match the component order
    from <- from[, colnames(conv_sigs)]
    to <- to[, colnames(conv_sigs)]
  }


  for (i in 1:nrow(conv_sigs)) {
    conv_sigs[i, ] <- conv_sigs[i, ] / from[i, ] * to[i, ]
    conv_sigs[i, ] <- conv_sigs[i, ] / sum(conv_sigs[i, ])
  }
  conv_sigs %>% t()
}

build_opps_matrix <- function(nsig, ncat, opps) {
  if (is.character(opps)) {
    strand <- ncat == 192
    opps <- matrix(rep(
      human_trinuc_freqs(opps, strand),
      nsig
    ), nrow = nsig, byrow = TRUE)
    colnames(opps) <- c(
      "A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "C[C>A]A", "C[C>A]C", "C[C>A]G", "C[C>A]T",
      "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T", "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T",
      "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T", "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T",
      "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T", "T[C>G]A", "T[C>G]C", "T[C>G]G", "T[C>G]T",
      "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T",
      "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",
      "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "C[T>A]A", "C[T>A]C", "C[T>A]G", "C[T>A]T",
      "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T", "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T",
      "A[T>C]A", "A[T>C]C", "A[T>C]G", "A[T>C]T", "C[T>C]A", "C[T>C]C", "C[T>C]G", "C[T>C]T",
      "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T", "T[T>C]A", "T[T>C]C", "T[T>C]G", "T[T>C]T",
      "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T", "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T",
      "G[T>G]A", "G[T>G]C", "G[T>G]G", "G[T>G]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T"
    )
  } else {
    opps <- as.matrix(opps) %>% t()
  }

  opps / sum(opps) * nsig
}

human_trinuc_freqs <- function(type = "human-genome", strand = FALSE) {
  if (type == "human-genome") {
    freq <- c(
      1.14e+08, 6.6e+07, 14300000, 91200000, 1.05e+08,
      74600000, 15700000, 1.01e+08, 81700000, 67600000,
      13500000, 79300000, 1.11e+08, 87500000, 12500000,
      1.25e+08, 1.14e+08, 6.6e+07, 14300000, 91200000,
      1.05e+08, 74600000, 15700000, 1.01e+08, 81700000,
      67600000, 13500000, 79300000, 1.11e+08, 87500000,
      12500000, 1.25e+08, 1.14e+08, 6.6e+07, 14300000,
      91200000, 1.05e+08, 74600000, 15700000, 1.01e+08,
      81700000, 67600000, 13500000, 79300000, 1.11e+08,
      87500000, 12500000, 1.25e+08, 1.17e+08, 75700000,
      1.04e+08, 1.41e+08, 73100000, 95500000, 1.15e+08,
      1.13e+08, 64300000, 53600000, 85200000, 82700000,
      1.18e+08, 1.12e+08, 1.07e+08, 2.18e+08, 1.17e+08,
      75700000, 1.04e+08, 1.41e+08, 73100000, 95500000,
      1.15e+08, 1.13e+08, 64300000, 53600000, 85200000,
      82700000, 1.18e+08, 1.12e+08, 1.07e+08, 2.18e+08,
      1.17e+08, 75700000, 1.04e+08, 1.41e+08, 73100000,
      95500000, 1.15e+08, 1.13e+08, 64300000, 53600000,
      85200000, 82700000, 1.18e+08, 1.12e+08, 1.07e+08,
      2.18e+08
    )
  }
  else if (type == "human-exome") {
    freq <- c(
      1940794, 1442408, 514826, 1403756, 2277398,
      2318284, 774498, 2269674, 1740752, 1968596, 631872,
      1734468, 1799540, 1910984, 398440, 2024770, 1940794,
      1442408, 514826, 1403756, 2277398, 2318284, 774498,
      2269674, 1740752, 1968596, 631872, 1734468, 1799540,
      1910984, 398440, 2024770, 1940794, 1442408, 514826,
      1403756, 2277398, 2318284, 774498, 2269674, 1740752,
      1968596, 631872, 1734468, 1799540, 1910984, 398440,
      2024770, 1299256, 1166912, 1555012, 1689928, 978400,
      2119248, 2650754, 1684488, 884052, 1173252, 1993110,
      1251508, 1391660, 1674368, 1559846, 2850934, 1299256,
      1166912, 1555012, 1689928, 978400, 2119248, 2650754,
      1684488, 884052, 1173252, 1993110, 1251508, 1391660,
      1674368, 1559846, 2850934, 1299256, 1166912, 1555012,
      1689928, 978400, 2119248, 2650754, 1684488, 884052,
      1173252, 1993110, 1251508, 1391660, 1674368, 1559846,
      2850934
    )
  }
  else {
    stop("'type' must be either \"human-genome\" or \"human-exome\"")
  }
  if (strand) {
    rep(freq / 2, 2)
  }
  else {
    freq
  }
}
