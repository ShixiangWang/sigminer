## Null-calibration family for the two P3 refit methods (presence test inside
## the judgment window). Reuses the shape of the 198-cell M3b calibration
## protocol (msuiter docs/ROADMAP.md:177; src/rust/src/calibrate.rs:47
## "11 x 3 x 6 grid"): null-derived synthetic fixtures, a-priori documented
## ceilings, measured-not-hand-picked acceptance. This lane runs the cell
## slice that pins its own two methods on hand-built D3 fixtures.
##
## Ceiling law: every bound below was fixed from the fixture-geometry
## calibration in /tmp/p3refit/logs/03b-nullcal-probe.log and
## /tmp/p3refit/logs/03c-stepwise-calibration.log BEFORE any implementation
## ran, from the null false-positive scale of the gates (ll thresholds at 8.0
## and alpha at 0.01 sit at chi-squared tails p ~ 1e-3 and below), never
## from implementation output.
##
## Judgment window: a reference signature counts as "declared" when its
## final relative exposure reaches p3_declare_threshold (1%).

## Run every spectrum of a stack through one method, returning declared flags.
p3_declared_matrix <- function(catalogs, W, method, ...) {
  e <- suppressMessages(suppressWarnings(
    do.call(sig_fit, c(
      list(catalogue_matrix = catalogs, sig = W, method = method),
      list(...)
    ))
  ))
  rel <- p3_normalize(e)
  t(rel >= p3_declare_threshold)
}


test_that("null arm: featureless backgrounds stay below the a-priori FP ceilings", {
  W <- p3_refit_sig()
  base::set.seed(31337, kind = "L'Ecuyer-CMRG")
  nulls <- do.call(cbind, lapply(seq_len(20L), function(i)
    p3_refit_null_draw(total = 1000L, sample_name = paste0("n", i))
  ))
  ## A-priori ceilings (documented above): spike columns A,B,C,F each below
  ## 0.10; the flat-like D,E (the known failure-mode family, exec-summary:30)
  ## each below 0.20; any-declaration per spectrum below 0.30.
  ceilings_spike <- c(A = .10, B = .10, C = .10, F = .10)
  ceilings_flat <- c(D = .20, E = .20)
  for (method in c("likelihood", "lrt")) {
    dec <- p3_declared_matrix(nulls, W, method)
    fp <- colMeans(dec)
    for (k in names(ceilings_spike)) {
      expect_lte(fp[[k]], ceilings_spike[[k]],
        label = paste(method, "spike-column FP", k))
    }
    for (k in names(ceilings_flat)) {
      expect_lte(fp[[k]], ceilings_flat[[k]],
        label = paste(method, "flat-column FP", k))
    }
    expect_lte(mean(rowSums(dec) > 0), .30, label = paste(method, "any-decl on null"))
  }
})

test_that("recall arm: planted signatures are declared above the a-priori floors", {
  W <- p3_refit_sig()
  ## Two-spike mixtures: both A and B must be declared in at least 90% of
  ## spectra; single-spike spectra: A in at least 95% (a-priori floors).
  base::set.seed(71113, kind = "L'Ecuyer-CMRG")
  mix <- do.call(cbind, lapply(seq_len(20L), function(i) {
    d <- as.integer(rmultinom(1L, 1000L,
      prob = as.numeric(0.5 * W[, "A"] + 0.5 * W[, "B"])))
    out <- matrix(d, ncol = 1L)
    dimnames(out) <- list(rownames(W), paste0("mix", i))
    out
  }))
  base::set.seed(71114, kind = "L'Ecuyer-CMRG")
  one <- do.call(cbind, lapply(seq_len(20L), function(i) {
    d <- as.integer(rmultinom(1L, 1000L, prob = as.numeric(W[, "A"])))
    out <- matrix(d, ncol = 1L)
    dimnames(out) <- list(rownames(W), paste0("one", i))
    out
  }))
  for (method in c("likelihood", "lrt")) {
    dec <- p3_declared_matrix(mix, W, method)
    expect_gte(mean(dec[, "A"] & dec[, "B"]), .90,
      label = paste(method, "both-spike recall"))
    expect_gte(mean(p3_declared_matrix(one, W, method)[, "A"]), .95,
      label = paste(method, "single-spike recall"))
  }
})

test_that("all routes speak through the same judgment window on planted spikes", {
  ## A planted spike pair must be declared by every route, legacy included --
  ## the judgment window is method-independent, which is what makes the FP
  ## ceilings above a like-for-like measurement.
  W <- p3_refit_sig()
  catalog <- p3_refit_catalog(c(A = .5, B = .5))
  for (method in c("QP", "likelihood", "lrt")) {
    dec <- p3_declared_matrix(catalog, W, method)
    expect_true(all(dec[1, c("A", "B")]), info = paste(method, "spike declaration"))
  }
})
