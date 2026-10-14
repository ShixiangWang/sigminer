# Tests for the pure K-arbitration core used by the bp de-novo route
# (parity-plan P2: port of the msuiter ms_select_k arbitration -- CV
# argmin -> stability veto -> diagnostics). The expected verdicts are the
# reference outcomes produced by the msuiter reference core itself; the
# provenance of the consumed fixture is recorded beside the fixture file.
#
# The units exercise the namespace-internal core directly (pkgload-loaded
# namespace exposes internals to the test env), so the semantics are
# pinned independently of the stochastic route that feeds it.

# Mirror of the reference frame schema (k, cv + stability columns); the
# core consumes exactly these, the rest of the evidence schema is filled
# by the survey layer.
kcore_ev <- function(k, cv, avg, mn) {
  vr <- .bp_veto_reason(avg, mn)
  data.frame(
    k = as.integer(k),
    cv_test_deviance = cv,
    avg_stability = avg,
    min_cluster_stability = mn,
    passes_veto = nzchar(vr) == FALSE,
    veto_reason = vr,
    is_argmin = rep(FALSE, length(k)),
    selected = rep(FALSE, length(k)),
    stringsAsFactors = FALSE)
}
kcore_set_argmin <- function(ev, k) {
  ev$is_argmin[ev$k == k] <- TRUE
  ev
}

test_that("the layer-2 veto mirrors the documented upstream thresholds", {
  # thresholds: avg >= 0.8 AND min_cluster >= 0.2 AND avg + min >= 1.0
  expect_equal(.bp_veto_reason(0.9, 0.3), "")
  expect_equal(.bp_veto_reason(0.8, 0.2), "")            # boundary inclusive
  expect_equal(.bp_veto_reason(0.8, 0.3), "")
  expect_equal(.bp_veto_reason(0.9, 0.2), "")
  expect_match(.bp_veto_reason(0.75, 0.3),
    "avg_stability 0.75 < 0.8")
  expect_match(.bp_veto_reason(0.9, 0.15),
    "min_cluster_stability 0.15 < 0.2")
  # avg 0.9 passes; min 0.05 fails the min bound AND the sum bound (0.95)
  expect_equal(.bp_veto_reason(0.9, 0.05),
    paste("min_cluster_stability 0.05 < 0.2",
      "avg + min_cluster_stability 0.95 < 1.0", sep = "; "))
  # NA stability fails conservatively with a single reason
  expect_equal(.bp_veto_reason(as.double(NA), 0.5),
    "stability statistics unavailable")
  expect_equal(.bp_veto_reason(0.9, as.double(NA)),
    "stability statistics unavailable")
})

test_that("S1: the CV argmin passing the veto is selected", {
  ev <- kcore_set_argmin(
    kcore_ev(2:5, c(300, 240, 260, 280), c(.9, .95, .85, .82),
      c(.35, .4, .3, .25)), 3L)
  out <- .bp_k_arbitrate(ev, cv = TRUE, mode = "explicit")
  expect_true(out$selected[2L])                # row for k = 3
  expect_equal(sum(out$selected), 1L)
})

test_that("S2: a vetoed argmin falls back to the first passing rank in CV-error order", {
  # k=3 is the argmin (cv 240) but unstable; passers k=4 (260), k=5 (280),
  # k=2 (300) -> fallback scan on CV-error ascending picks k=4
  ev <- kcore_set_argmin(
    kcore_ev(2:5, c(300, 240, 260, 280), c(.9, .5, .85, .82),
      c(.35, .1, .3, .25)), 3L)
  out <- .bp_k_arbitrate(ev, cv = TRUE, mode = "explicit")
  expect_false(out$selected[2L])
  expect_true(out$selected[3L])
})

test_that("S3: with layer 1 switched off the lowest passing rank wins", {
  ev <- kcore_ev(2:5, c(300, 240, 260, 280), c(.5, .5, .85, .82),
    c(.1, .1, .3, .25))
  out <- .bp_k_arbitrate(ev, cv = FALSE, mode = "explicit")
  expect_false(any(out$is_argmin))             # layer 1 absent erases the column
  expect_true(out$selected[3L])                # k=4, the lowest passing rank
})

test_that("S4: no rank passing raises in the explicit mode", {
  ev <- kcore_set_argmin(
    kcore_ev(2:5, c(300, 240, 260, 280), rep(.4, 4), rep(.1, 4)), 3L)
  expect_error(.bp_k_arbitrate(ev, cv = TRUE, mode = "explicit"),
    "undecidable")
})

test_that("S5: no rank passing degrades to the argmin with a warning in pipeline mode", {
  ev <- kcore_set_argmin(
    kcore_ev(2:5, c(300, 240, 260, 280), rep(.4, 4), rep(.1, 4)), 3L)
  expect_warning(out <- .bp_k_arbitrate(ev, cv = TRUE, mode = "degrade"),
    "undecidable")
  expect_true(out$selected[2L])                # the argmin (k = 3)
})

test_that("S6: the lowest tied rank carries the argmin mark into the core", {
  # all ranks pass -> the carried argmin (lowest of the cv tie at 240) wins
  ev <- kcore_set_argmin(
    kcore_ev(2:5, c(300, 240, 240, 280), c(.9, .95, .9, .82),
      c(.35, .4, .35, .25)), 3L)
  out <- .bp_k_arbitrate(ev, cv = TRUE, mode = "explicit")
  expect_true(out$selected[2L])                # k = 3
})

test_that("S7: fallback among passers breaks CV ties by the lower rank", {
  # argmin k=2 (240) vetoed; passers k=3 (240), k=4 (260), k=5 (280);
  # the scan takes k=3 first
  ev <- kcore_set_argmin(
    kcore_ev(2:5, c(240, 240, 260, 280), c(.4, .95, .9, .82),
      c(.1, .4, .35, .25)), 2L)
  out <- .bp_k_arbitrate(ev, cv = TRUE, mode = "explicit")
  expect_true(out$selected[2L])                # k = 3
  # vetoing k=3 leaves passers k=4 and k=5 tied at 260 -> the scan takes
  # the lower rank k=4
  ev2 <- kcore_set_argmin(
    kcore_ev(2:5, c(240, 240, 260, 260), c(.4, .4, .9, .85),
      c(.1, .4, .35, .3)), 2L)
  out2 <- .bp_k_arbitrate(ev2, cv = TRUE, mode = "explicit")
  expect_true(out2$selected[3L])               # k = 4
})

test_that("the argmin mark is computed on finite CV values only, lowest rank on ties", {
  # survey layer fills is_argmin; port law: among finite cv entries the
  # first minimum (which.min = lowest rank on ties); all-NA curve -> none
  expect_equal(.bp_cv_argmin(c(300, 240, 240, 280)), 2L)
  expect_equal(.bp_cv_argmin(rep(as.double(NA), 4)), integer(0))
  expect_equal(.bp_cv_argmin(c(as.double(NA), 240, 100, 280)), 3L)
  expect_equal(.bp_cv_argmin(c(240, as.double(NA), 240, 280)), 1L)
})

test_that("input validation fails fast with named messages", {
  good <- kcore_ev(2:3, c(300, 240), c(.9, .9), c(.3, .3))
  expect_error(.bp_k_arbitrate(as.data.frame(list(k = 1)), cv = TRUE,
    mode = "explicit"), "column")
  expect_error(.bp_k_arbitrate(transform(good, k = c(3L, 2L)), cv = TRUE,
    mode = "explicit"), "increasing")
  expect_error(.bp_k_arbitrate(transform(good, k = c(2L, 2L)), cv = TRUE,
    mode = "explicit"), "unique")
  expect_error(.bp_k_arbitrate(good, cv = NA, mode = "explicit"), "cv")
  expect_error(.bp_k_arbitrate(good, cv = TRUE, mode = "wat"), "should be one of")
  expect_error(.bp_k_arbitrate(transform(good, passes_veto = c("yes", "no")),
    cv = TRUE, mode = "explicit"), "column")
  # at most one row may carry the layer-1 argmin mark
  expect_error(.bp_k_arbitrate(transform(good, is_argmin = c(TRUE, TRUE)),
    cv = TRUE, mode = "explicit"), "argmin")
})
