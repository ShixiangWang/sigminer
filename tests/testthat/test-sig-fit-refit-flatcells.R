## Flat-signature cell family (the documented SBS5/40-style failure mode,
## msuiter docs/research/00-executive-summary.md:30). The P3 acceptance law:
## the two new methods must not regress these cells below the legacy route.
##
## Fixture: hand-built near-flat family over {D uniform, E tilted, F decoy}
## with a 0.2 D + 0.8 E mixture (single column cannot carry the mix -- only
## the pair fits it, cf. /tmp/p3refit/logs/03b-nullcal-probe.log -- so this
## cell punishes any decision layer that drops too much of the flat mass).
##
## A-priori margins fixed before implementation:
##   - capture_new >= capture_legacy - 0.15 (share points; the new methods
##     must stay within 15 points of the legacy route's flat-mass capture)
##   - capture_new >= 0.70 (absolute floor: the planted flat mass must stay
##     majority-owned by the flat family even at the known failure mode)
##   - no non-flat column may eat more than 0.10 of the mix

p3_flat_capture <- function(catalogs, Wsub, method, ...) {
  e <- suppressMessages(suppressWarnings(
    do.call(sig_fit, c(
      list(catalogue_matrix = catalogs, sig = Wsub, method = method),
      list(...)
    ))
  ))
  rel <- p3_normalize(e)
  ## mean share of the mix owned by the flat family columns D and E
  mean(rel["D", ] + rel["E", ])
}

test_that("flat cells: new methods hold flat-mass capture within the margins", {
  W <- p3_refit_sig()
  Wsub <- W[, c("D", "E", "F"), drop = FALSE]
  base::set.seed(40404, kind = "L'Ecuyer-CMRG")
  catalogs <- do.call(cbind, lapply(seq_len(8L), function(i) {
    d <- as.integer(rmultinom(1L, 1000L,
      prob = as.numeric(0.2 * Wsub[, "D"] + 0.8 * Wsub[, "E"])))
    out <- matrix(d, ncol = 1L)
    dimnames(out) <- list(rownames(Wsub), paste0("flat", i))
    out
  }))
  cap_legacy <- p3_flat_capture(catalogs, Wsub, "QP")
  expect_gte(cap_legacy, .85) ## the fixture itself must be well-formed
  for (method in c("likelihood", "lrt")) {
    cap <- p3_flat_capture(catalogs, Wsub, method)
    expect_true(all(is.finite(as.matrix(
      suppressMessages(suppressWarnings(
        sig_fit(catalogs, sig = Wsub, method = method)
      ))
    ))))
    expect_gte(cap, cap_legacy - 0.15, label = paste(method, "vs legacy capture"))
    expect_gte(cap, 0.70, label = paste(method, "absolute flat floor"))
    ## The decoy must not eat the mix.
    e <- suppressMessages(suppressWarnings(
      sig_fit(catalogs, sig = Wsub, method = method)
    ))
    rel <- p3_normalize(e)
    expect_lte(max(rel["F", ]), 0.10, label = paste(method, "decoy share"))
  }
})

test_that("flat cells: likelihood route keeps the calibrated pruning verdict", {
  ## On the deterministic 0.2 D + 0.8 E cell the calibrated defaults must
  ## keep E (removal cost 40.98) and prune D (cost 2.805 < ll_keep 8.0)
  ## and F (cost 0), while all flat mass remains captured on {E} alone.
  W <- p3_refit_sig()
  Wsub <- W[, c("D", "E", "F"), drop = FALSE]
  catalog <- p3_refit_catalog(c(D = .2, E = .8))
  dimnames(catalog)[[1]] <- rownames(Wsub)
  e <- suppressMessages(suppressWarnings(
    sig_fit(catalog, sig = Wsub, method = "likelihood")
  ))
  abs <- as.numeric(e)
  names(abs) <- rownames(e)
  expect_equal(abs[["D"]], 0)
  expect_equal(abs[["F"]], 0)
  expect_gt(abs[["E"]], 990)
})
