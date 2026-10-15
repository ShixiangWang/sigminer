## method = "likelihood" (MuSiCal-branch semantics, P3 refit lane).
##
## Primary-source gate: /tmp/p3refit/logs/00-primary-source-gate.md sec. 1 --
## "Matching and refitting with multinomial likelihood-based sparse NNLS."
## "bidirectional stepwise algorithm"; acceptance thresholds are
## study-optimized parameters ("two likelihood thresholds were optimized",
## "optimal threshold was chosen through a grid search") and are therefore
## exposed as arguments here rather than hard-coded.

## The shared p3_check_method_argument contract helper lives in
## helper-refit-synthetic.R so both refit-method families pin it identically.

test_that("method='likelihood' resolves through the house match.arg idiom", {
  W <- p3_refit_sig()
  catalog <- p3_refit_catalog(c(A = .5, B = .5))
  p3_check_method_argument(catalog, W, exact = "likelihood", partial = "li")
})

test_that("method='likelihood' keeps the two true signatures of a clean mix", {
  W <- p3_refit_sig()
  catalog <- p3_refit_catalog(c(A = .5, B = .5))
  expo <- suppressMessages(suppressWarnings(
    sig_fit(catalog, sig = W, method = "likelihood")
  ))
  expect_true(is.matrix(expo))
  expect_equal(dim(expo), c(6L, 1L))
  expect_equal(rownames(expo), colnames(W))
  rel <- as.numeric(expo) / sum(as.numeric(expo))
  names(rel) <- rownames(expo)
  ## True signatures carry the bulk of the exposure ...
  expect_gt(rel[["A"]], .30)
  expect_gt(rel[["B"]], .30)
  expect_lt(abs(rel[["A"]] - rel[["B"]]), .10)
  ## ... and every irrelevant reference is fully dropped.
  expect_equal(as.numeric(expo[c("C", "D", "E", "F"), ]), rep(0, 4))
})

test_that("method='likelihood' selects the single true spike", {
  W <- p3_refit_sig()
  catalog <- p3_refit_catalog(c(C = 1))
  expo <- suppressMessages(suppressWarnings(
    sig_fit(catalog, sig = W, method = "likelihood")
  ))
  rel <- as.numeric(expo)
  names(rel) <- rownames(expo)
  expect_gt(rel[["C"]], .90)
  ## unname(): the slice inherits the matrix row names; the pin is on values.
  expect_equal(unname(rel[rel != rel[["C"]]]), rep(0, 5))
})

test_that("ll_keep/ll_add are real parameters with the calibrated semantics", {
  W <- p3_refit_sig()
  ## Deterministic flat-family cell (0.2 D + 0.8 E). Fixture refit costs were
  ## measured pre-implementation in /tmp/p3refit/logs/03c-stepwise-calibration
  ## .log: removal-cost of E is 40.98 (>= the calibrated ll_keep 8), removal of
  ## D costs 2.805 and of F 0.0 (both < 8), and the E-alone gain over the
  ## empty reference is 38.17 (>= the calibrated ll_add 8), every other single
  ## column gain is far below it.
  catalog <- p3_refit_catalog(c(D = .2, E = .8))
  kept_names <- function(...) {
    e <- suppressMessages(suppressWarnings(
      sig_fit(catalog, sig = W, method = "likelihood", ...)
    ))
    v <- as.numeric(e)
    names(v) <- rownames(e)
    sort(names(v[v > 0]))
  }
  ## (a) Calibration pin: at the defaults the sweep prunes D and F and keeps
  ## the only member whose removal is significant at the gate: exactly {E}.
  expect_equal(kept_names(), "E")
  ## (b) Never-drop equivalence: with an unfireable keep-gate (ll_keep = -1,
  ## no removal cost is negative) the decision layer keeps the full reference
  ## set, so its refit must be exactly the legacy QP route's answer (the same
  ## QP simplex engine, no selection delta).
  e_full <- suppressMessages(suppressWarnings(
    sig_fit(catalog, sig = W, method = "likelihood", ll_add = 1e9, ll_keep = -1)
  ))
  e_qp <- suppressMessages(suppressWarnings(
    sig_fit(catalog, sig = W, method = "QP")
  ))
  expect_equal(as.numeric(e_full), as.numeric(e_qp), tolerance = 1e-12)
  ## (c) Total-collapse recovery: a keep-gate no member can pass empties the
  ## set; the forward gate then re-admits exactly what clears ll_add against
  ## the empty reference -- per the 03c measurements that is E alone.
  expect_equal(kept_names(ll_keep = 1e9), "E")
})

test_that("ll_add is accepted as a paired threshold argument", {
  W <- p3_refit_sig()
  catalog <- p3_refit_catalog(c(A = .5, B = .5))
  a <- suppressMessages(suppressWarnings(
    sig_fit(catalog, sig = W, method = "likelihood")
  ))
  b <- suppressMessages(suppressWarnings(
    sig_fit(catalog, sig = W, method = "likelihood",
      ll_add = 8.0, ll_keep = 8.0
    )
  ))
  expect_equal(as.numeric(a), as.numeric(b))
})

test_that("method='likelihood' is deterministic and seed-free", {
  W <- p3_refit_sig()
  catalog <- p3_refit_catalog(c(A = .6, C = .4))
  e1 <- suppressMessages(suppressWarnings(
    sig_fit(catalog, sig = W, method = "likelihood")
  ))
  e2 <- suppressMessages(suppressWarnings(
    sig_fit(catalog, sig = W, method = "likelihood")
  ))
  expect_equal(e1, e2)
})

test_that("method='likelihood' handles a zero-count sample like the legacy route", {
  W <- p3_refit_sig()
  zero <- matrix(0L, nrow = 30L, ncol = 1L)
  dimnames(zero) <- list(rownames(W), "zero_smp")
  ez <- suppressMessages(suppressWarnings(
    sig_fit(zero, sig = W, method = "likelihood")
  ))
  expect_equal(as.numeric(ez), rep(0, 6))
})

test_that("method='likelihood' honors return_class and return_error like QP", {
  W <- p3_refit_sig()
  catalog <- p3_refit_catalog(c(B = .5, C = .5))
  res <- suppressMessages(suppressWarnings(
    sig_fit(catalog, sig = W, method = "likelihood",
      return_class = "data.table", return_error = TRUE
    )
  ))
  expect_true(is.list(res))
  expect_true(all(c("expo", "errors", "cosine") %in% names(res)))
  expect_true(data.table::is.data.table(res$expo))
  expect_equal(names(res$expo), c("sample", colnames(W)))
  expect_equal(length(res$errors), 1L)
  expect_equal(length(res$cosine), 1L)
  expect_true(is.finite(res$cosine))
})

test_that("method='likelihood' rejects auto_reduce with an explicit error", {
  W <- p3_refit_sig()
  catalog <- p3_refit_catalog(c(A = .5, B = .5))
  x <- as.numeric(catalog)
  ## House stop idiom: send_stop() raises a bare stop() and routes the human
  ## text through a cli alert, which surfaces as a message condition. Pin the
  ## refusal text at the direct entry point where no purrr wrapper can interfere.
  expect_message(
    expect_error(decompose_LIKELIHOOD(x, 0, "smp", W, auto_reduce = TRUE)),
    "likelihood"
  )
  expect_message(
    expect_error(decompose_LIKELIHOOD(x, 0, "smp", W, auto_reduce = TRUE)),
    "auto_reduce"
  )
  ## Wiring pin: the refusal must also surface through the sig_fit route.
  expect_error(
    suppressMessages(sig_fit(catalog, sig = W,
      method = "likelihood", auto_reduce = TRUE
    ))
  )
})

test_that("method='likelihood' absolute exposures track the planted counts", {
  W <- p3_refit_sig()
  catalog <- p3_refit_catalog(c(A = .6, C = .4), total = 1000L)
  expo <- suppressMessages(suppressWarnings(
    sig_fit(catalog, sig = W, method = "likelihood")
  ))
  abs_expo <- as.numeric(expo)
  names(abs_expo) <- rownames(expo)
  expect_lt(abs(abs_expo[["A"]] - 600), 60)
  expect_lt(abs(abs_expo[["C"]] - 400), 60)
  expect_lt(abs(sum(abs_expo) - 1000), 5)
})
