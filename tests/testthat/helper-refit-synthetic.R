## ---------------------------------------------------------------------------
## Shared synthetic fixtures for the P3 refit-method test family
## (method = "likelihood" / "lrt", added by the P3 refit lane).
##
## D3 law: every numeric here is hand-built. Nothing in this family derives
## from COSMIC or any other real reference distribution.
## ---------------------------------------------------------------------------

## Six hand-built reference signatures over a 30-component grid:
##   A, B, C  spike-like, clearly distinct supports
##   D        perfectly flat (uniform)
##   E        near-flat (mild monotone tilt, still close to uniform)
##   F        decoy spike on its own component
p3_refit_sig <- function() {
  W <- matrix(0, nrow = 30, ncol = 6)
  dimnames(W) <- list(paste0("c", 1:30), c("A", "B", "C", "D", "E", "F"))
  W["c1", "A"] <- .80; W[c("c2", "c3", "c4", "c5"), "A"] <- .05
  W["c10", "B"] <- .80; W[c("c6", "c7", "c8", "c9"), "B"] <- .05
  W["c20", "C"] <- .75; W[c("c15", "c16", "c17", "c18", "c19"), "C"] <- .05
  W[, "D"] <- 1 / 30
  ## Near-flat shape column: a mild global tilt (the hand-built analogue of
  ## the documented flat-signature family, cf. docs/research/
  ## 00-executive-summary.md:30). Carries real structure so its presence is
  ## statistically detectable, while a uniform+tilt mix remains only
  ## pair-sufficient (measured in /tmp/p3refit/logs/03b-nullcal-probe.log).
  tilt <- 1 + .45 * ((1:30) - 15.5) / 15.5
  W[, "E"] <- tilt / sum(tilt)
  W["c28", "F"] <- .55; W[paste0("c", 21:27), "F"] <- .45 / 7
  W / apply(W, 2, sum)
}

## Deterministic integer catalog: expected reconstruction W %*% e scaled to a
## fixed total and rounded, with the rounding residual parked on the largest
## cell so the column sum is exact. No RNG involved.
p3_refit_catalog <- function(weights, total = 1000L, sample_name = "smp") {
  W <- p3_refit_sig()
  e <- rep(0, ncol(W))
  names(e) <- colnames(W)
  e[names(weights)] <- as.numeric(weights)
  mu <- as.numeric(W %*% e) * total
  d <- round(mu)
  j <- which.max(d)
  d[j] <- d[j] + (total - sum(d))
  out <- matrix(as.integer(d), ncol = 1L)
  dimnames(out) <- list(rownames(W), sample_name)
  out
}

## Draw one null spectrum from a featureless uniform background. This is the
## null-derived fixture shape for the presence-calibration family: it carries
## no signature structure, so any declaration of a reference column by a
## selection-based method is a genuine false positive (the hard FP cases are
## the flat-like columns D/E, cf. the documented flat-signature failure mode
## in docs/research/00-executive-summary.md:30).
p3_refit_null_draw <- function(total = 1000L, sample_name = "null") {
  d <- as.integer(rmultinom(1L, total, prob = rep(1 / 30, 30)))
  out <- matrix(d, ncol = 1L)
  dimnames(out) <- list(rownames(p3_refit_sig()), sample_name)
  out
}

## Count vector for which a signature counts as "declared" in the judgment
## window of the null-calibration family (a-priori, documented in the test).
p3_declare_threshold <- 0.01

## Column-wise share matrix (each sample's exposures over its own total; an
## all-zero sample stays all-zero).
p3_normalize <- function(e) {
  cs <- apply(e, 2L, sum)
  cs[cs == 0] <- 1
  e / rep(cs, each = nrow(e))
}

## Shared method-argument contract for both new refit methods: choices must
## resolve through the house match.arg idiom (case-sensitive exact or
## unique-prefix), and a refusal must enumerate every valid choice.
p3_check_method_argument <- function(catalog, W, exact, partial) {
  expect_silent(suppressMessages(
    sig_fit(catalog, sig = W, method = exact)
  ))
  expect_silent(suppressMessages(
    sig_fit(catalog, sig = W, method = partial)
  ))
  ## Ambiguous prefix ("l" for both new choices) must be refused, naming both.
  err_amb <- tryCatch({
    suppressMessages(sig_fit(catalog, sig = W, method = "l")); "no-error"
  }, error = function(e) conditionMessage(e))
  expect_match(err_amb, "likelihood")
  expect_match(err_amb, "lrt")
  ## Case sensitivity: an unmatched case variant errors, listing all choices.
  err_case <- tryCatch({
    suppressMessages(sig_fit(catalog, sig = W, method = "L")); "no-error"
  }, error = function(e) conditionMessage(e))
  expect_match(err_case, "QP")
  expect_match(err_case, "NNLS")
  expect_match(err_case, "SA")
  expect_match(err_case, "likelihood")
  expect_match(err_case, "lrt")
}

