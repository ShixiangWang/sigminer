## method = "lrt" (PASA/mSigAct-branch semantics, P3 refit lane).
##
## Primary-source gate: /tmp/p3refit/logs/00-primary-source-gate.md sec. 2 --
## verbatim algorithm attestation: the LRT statistic
## "lambda = -2(log P(D|A2) - log P(D|A1))" with A2 a proper subset of A1,
## "which follows a chi-squared distribution with degrees of freedom =
## |A1| - |A2|"; Step 1 "signature presence test" per signature (full minus
## k against full); Step 2 minimal-set semantics "The set is minimal in the
## sense that removing any signature results in a likelihood ratio test
## giving a p value < alpha"; stop condition "The algorithm stops when the
## reconstruction is 'good enough' as assessed by a likelihood ratio test";
## order-independence "the outer loop always selects the signature that
## improves the likelihood the most ... does not depend on the order".
##
## The alpha significance level is a paper-parameter; its default value is
## lane-calibrated (documented), never paper-sourced.

test_that("method='lrt' resolves and shares the argument contract with 'likelihood'", {
  W <- p3_refit_sig()
  catalog <- p3_refit_catalog(c(A = .5, B = .5))
  p3_check_method_argument(catalog, W, exact = "lrt", partial = "lr")
})

test_that("method='lrt' keeps the two true signatures of a clean mix", {
  W <- p3_refit_sig()
  catalog <- p3_refit_catalog(c(A = .5, B = .5))
  expo <- suppressMessages(suppressWarnings(
    sig_fit(catalog, sig = W, method = "lrt")
  ))
  expect_true(is.matrix(expo))
  expect_equal(dim(expo), c(6L, 1L))
  expect_equal(rownames(expo), colnames(W))
  rel <- as.numeric(expo)
  names(rel) <- rownames(expo)
  expect_gt(rel[["A"]], .30)
  expect_gt(rel[["B"]], .30)
  expect_lt(abs(rel[["A"]] - rel[["B"]]), .10)
  expect_equal(as.numeric(expo[c("C", "D", "E", "F"), ]), rep(0, 4))
})

test_that("method='lrt' selects the single true spike", {
  W <- p3_refit_sig()
  catalog <- p3_refit_catalog(c(C = 1))
  expo <- suppressMessages(suppressWarnings(
    sig_fit(catalog, sig = W, method = "lrt")
  ))
  rel <- as.numeric(expo)
  names(rel) <- rownames(expo)
  expect_gt(rel[["C"]], .90)
  ## unname(): the slice inherits the matrix row names; the pin is on values.
  expect_equal(unname(rel[rel != rel[["C"]]]), rep(0, 5))
})

test_that("method='lrt' sparses the flat cell down to the sufficient column", {
  W <- p3_refit_sig()
  ## 0.2 D + 0.8 E: only the pair {D,E} fits as well as E plus its own
  ## column, but removing D costs just lambda = 2.805 -> p = 0.094 >= alpha
  ## (03c calibration log), so a minimal set at alpha = 0.01 keeps E alone.
  catalog <- p3_refit_catalog(c(D = .2, E = .8))
  expo <- suppressMessages(suppressWarnings(
    sig_fit(catalog, sig = W, method = "lrt")
  ))
  abs <- as.numeric(expo)
  names(abs) <- rownames(expo)
  ## The flat-mix mass must stay on the flat-family columns (no-regression
  ## law, docs/research/00-executive-summary.md:30): D + E carry it all.
  expect_gt(abs[["D"]] + abs[["E"]], 990)
  expect_equal(abs[["F"]], 0)
  ## And sparsity holds: the insufficient partner is dropped at alpha = 0.01.
  expect_equal(abs[["D"]], 0)
  expect_gt(abs[["E"]], 990)
})

test_that("alpha = 0 accepts nothing: the presence test clears the board", {
  W <- p3_refit_sig()
  catalog <- p3_refit_catalog(c(A = .5, B = .5))
  expo <- suppressMessages(suppressWarnings(
    sig_fit(catalog, sig = W, method = "lrt", alpha = 0)
  ))
  expect_equal(as.numeric(expo), rep(0, 6))
})

test_that("alpha = 1 keeps everything the full-set refit would keep", {
  W <- p3_refit_sig()
  catalog <- p3_refit_catalog(c(A = .5, B = .5))
  legacy <- suppressMessages(suppressWarnings(
    sig_fit(catalog, sig = W, method = "QP")
  ))
  new <- suppressMessages(suppressWarnings(
    sig_fit(catalog, sig = W, method = "lrt", alpha = 1)
  ))
  ## At alpha = 1 no LRT can ever clear the bar, so the decision layer
  ## degenerates to "keep the full candidate set" and the refit engine is
  ## the same one the legacy QP route drives: identical exposures.
  expect_equal(new, legacy)
})

test_that("method='lrt' is deterministic and seed-free", {
  W <- p3_refit_sig()
  catalog <- p3_refit_catalog(c(A = .6, C = .4))
  e1 <- suppressMessages(suppressWarnings(
    sig_fit(catalog, sig = W, method = "lrt")
  ))
  e2 <- suppressMessages(suppressWarnings(
    sig_fit(catalog, sig = W, method = "lrt")
  ))
  expect_equal(e1, e2)
})

test_that("method='lrt' handles a zero-count sample like the legacy route", {
  W <- p3_refit_sig()
  zero <- matrix(0L, nrow = 30L, ncol = 1L)
  dimnames(zero) <- list(rownames(W), "zero_smp")
  ez <- suppressMessages(suppressWarnings(
    sig_fit(zero, sig = W, method = "lrt")
  ))
  expect_equal(as.numeric(ez), rep(0, 6))
})

test_that("method='lrt' honors return_class and return_error like QP", {
  W <- p3_refit_sig()
  catalog <- p3_refit_catalog(c(B = .5, C = .5))
  res <- suppressMessages(suppressWarnings(
    sig_fit(catalog, sig = W, method = "lrt",
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

test_that("method='lrt' rejects auto_reduce with an explicit error", {
  W <- p3_refit_sig()
  catalog <- p3_refit_catalog(c(A = .5, B = .5))
  x <- as.numeric(catalog)
  ## House stop idiom: send_stop() raises a bare stop() and routes the human
  ## text through a cli alert, which surfaces as a message condition. Pin the
  ## refusal text at the direct entry point where no purrr wrapper can interfere.
  expect_message(
    expect_error(decompose_LRT(x, 0, "smp", W, auto_reduce = TRUE)),
    "lrt"
  )
  expect_message(
    expect_error(decompose_LRT(x, 0, "smp", W, auto_reduce = TRUE)),
    "auto_reduce"
  )
  ## Wiring pin: the refusal must also surface through the sig_fit route.
  expect_error(
    suppressMessages(sig_fit(catalog, sig = W,
      method = "lrt", auto_reduce = TRUE
    ))
  )
})
