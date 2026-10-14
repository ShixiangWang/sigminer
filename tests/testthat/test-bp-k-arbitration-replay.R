# Fixture replay for the bp K-arbitration core (parity plan P2).
#
# Consumes the committed reference capture
# fixtures/bp-k-arb/msuiter_core_cases.rds: seven hand-built evidence
# frames (S1..S7 of the parity design) run through the msuiter reference
# core (.ms_k_arbitrate, msuiter main@dad30440) on the controller host,
# recorded with their run parameters. Replaying them here asserts the
# ported core takes the identical decision on every recorded frame.
#
# The unit is msuiter-absent-safe: the .rds payload is plain data (data
# frames, scalars, NULL), loaded with readRDS only. The provenance
# generator next to the fixture re-produces both .rds files verbatim on
# the controller host; it is never sourced by R CMD check.

karb_fixture_path <- function(... ) {
  file.path(test_path(), "fixtures", "bp-k-arb", ...)
}

# Boundary validation of the consumed capture (input-validation law): a
# silently corrupt fixture must fail loudly here, not skew the replay.
validate_core_capture <- function(capture) {
  if (!is.list(capture) || !all(c("root_seed", "cases") %in% names(capture))) {
    stop("core capture must carry `root_seed` and `cases`", domain = NA)
  }
  for (nm in names(capture$cases)) {
    cs <- capture$cases[[nm]]
    need <- c("ev", "cv", "mode", "expect_error", "warned")
    if (!all(need %in% names(cs))) {
      stop(sprintf("core case `%s` is missing recorded field(s): %s", nm,
        paste0(dQuote(setdiff(need, names(cs))), collapse = ", ")),
        domain = NA)
    }
    if (!is.data.frame(cs$ev) ||
        !all(.BP_K_EVIDENCE_REQUIRED %in% names(cs$ev))) {
      stop(sprintf("core case `%s` carries an invalid evidence frame", nm),
        domain = NA)
    }
    if (!identical(cs$mode, "explicit") && !identical(cs$mode, "degrade")) {
      stop(sprintf("core case `%s` records an unknown mode", nm), domain = NA)
    }
    if (!isTRUE(cs$expect_error) && !isFALSE(cs$expect_error)) {
      stop(sprintf("core case `%s` has no recorded error outcome", nm),
        domain = NA)
    }
  }
  capture
}

capture <- validate_core_capture(readRDS(karb_fixture_path(
  "msuiter_core_cases.rds")))

test_that("the reference capture is the recorded seven-scenario battery", {
  expect_equal(capture$root_seed, 20261015L)
  expect_equal(names(capture$cases),
    c("S1_argmin_passes", "S2_argmin_vetoed_fallback",
      "S3_cv_off_lowest_passing", "S4_none_passes_explicit_error",
      "S5_none_passes_degrade_to_argmin", "S6_cv_tie_lowest_rank",
      "S7_vetoed_argmin_tie_fallback"))
})

test_that("every recorded core case replays to the identical reference decision", {
  for (nm in names(capture$cases)) {
    cs <- capture$cases[[nm]]
    if (cs$expect_error) {
      # layer-3 parity: the reference errored; the port must too, in the
      # explicit mode, with the undecidable diagnosis.
      expect_error(.bp_k_arbitrate(cs$ev, cv = cs$cv, mode = cs$mode),
        "undecidable", info = nm)
      next
    }
    if (cs$warned) {
      expect_warning(out <- .bp_k_arbitrate(cs$ev, cv = cs$cv,
        mode = cs$mode), "undecidable", info = nm)
    } else {
      out <- .bp_k_arbitrate(cs$ev, cv = cs$cv, mode = cs$mode)
    }
    expect_equal(out$selected, cs$selected, info = nm)
    expect_equal(out$is_argmin, cs$is_argmin_out, info = nm)
    expect_equal(sum(out$selected), 1L, info = nm)
  }
})

test_that("the replayed decisions are stable across repeated invocations", {
  # determinism cross-check: the core is a pure function, so two
  # invocations per case must agree bit-for-bit through serialization.
  for (nm in names(capture$cases)) {
    cs <- capture$cases[[nm]]
    if (cs$expect_error) next
    a <- suppressWarnings(.bp_k_arbitrate(cs$ev, cv = cs$cv, mode = cs$mode))
    b <- suppressWarnings(.bp_k_arbitrate(cs$ev, cv = cs$cv, mode = cs$mode))
    expect_equal(serialize(a, connection = NULL),
      serialize(b, connection = NULL), info = nm)
  }
})
