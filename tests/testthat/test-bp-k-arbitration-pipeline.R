# Pipeline-level acceptance units for the ported K-arbitration (parity plan
# P2): bp_select_k() -- the CV-survey + three-layer arbitration over a rank
# grid -- and its opt-in hookup into the bp extraction route.
#
# Load-bearing parity evidence: the expected verdicts are the committed
# msuiter reference captures fixtures/bp-k-arb/msuiter_verdicts.rds
# (produced on the controller host by the committed generator; provenance
# header there). The calibration cells themselves are rebuilt on the fly
# from the committed pure builder synthetic_cells.R under the pinned root
# seed, so tested cells and recorded cells are the same numbers by
# construction. Everything is SYNTHETIC (D3).

# ---- fixture plumbing -------------------------------------------------------

karb_cells_env <- new.env()
source(file.path(test_path(), "fixtures", "bp-k-arb", "synthetic_cells.R"),
  local = karb_cells_env)

karb_verdict_capture <- function() {
  cap <- readRDS(file.path(test_path(), "fixtures", "bp-k-arb",
    "msuiter_verdicts.rds"))
  # Boundary validation of the consumed capture (input-validation law).
  if (!is.list(cap) || !all(c("root_seed", "k_grid", "cells") %in% names(cap))) {
    stop("verdict capture must carry root_seed, k_grid and cells", domain = NA)
  }
  for (nm in names(cap$cells)) {
    cl <- cap$cells[[nm]]
    if (!all(c("name", "planted_k", "k", "dim", "evidence") %in% names(cl))) {
      stop(sprintf("verdict cell `%s` is missing recorded fields", nm),
        domain = NA)
    }
    if (!is.data.frame(cl$evidence) || nrow(cl$evidence) != length(cap$k_grid)) {
      stop(sprintf("verdict cell `%s` carries an invalid evidence table", nm),
        domain = NA)
    }
  }
  cap
}

# ---- T1: K-agreement with the committed reference verdicts -----------------

test_that("bp_select_k reaches the committed msuiter verdict on every tier", {
  cap <- karb_verdict_capture()
  cells <- karb_cells_env$cal_cells(karb_cells_env$CAL_ROOT_SEED)
  expect_equal(vapply(cells, function(cl) cl$name, character(1L)),
    names(cap$cells))
  for (cl in cells) {
    expect_equal(dim(cl$counts), cap$cells[[cl$name]]$dim, info = cl$name)
    res <- bp_select_k(cl$counts, k_grid = cap$k_grid,
      seed = cap$root_seed)
    expect_equal(res$k, cap$cells[[cl$name]]$k, info = cl$name)
    expect_equal(res$k, cl$planted_k, info = cl$name)
    ev <- res$k_evidence
    expect_equal(ev$k, as.integer(cap$k_grid), info = cl$name)
    expect_equal(sum(ev$selected), 1L, info = cl$name)
    sel <- ev[ev$selected, ]
    expect_true(sel$passes_veto, info = cl$name)
    expect_false(anyNA(ev$cv_test_deviance), info = cl$name)
    expect_false(anyNA(ev$passes_veto), info = cl$name)
  }
})

# ---- T2: single-root-seed law ----------------------------------------------

test_that("one root seed yields byte-identical verdicts and diagnostics", {
  cells <- karb_cells_env$cal_cells(karb_cells_env$CAL_ROOT_SEED)
  cl <- cells[[2L]]
  a <- bp_select_k(cl$counts, k_grid = 2:5, seed = 20261015L)
  b <- bp_select_k(cl$counts, k_grid = 2:5, seed = 20261015L)
  expect_equal(a$k, b$k)
  expect_identical(serialize(a$k_evidence, connection = NULL),
    serialize(b$k_evidence, connection = NULL))
  # A different root seed is a different stream: the diagnostics must
  # actually move (no frozen pseudoreplication); the decision columns stay
  # well-formed for any stream.
  d <- bp_select_k(cl$counts, k_grid = 2:5, seed = 99991L)
  expect_false(identical(serialize(a$k_evidence, connection = NULL),
    serialize(d$k_evidence, connection = NULL)))
  expect_equal(sum(d$k_evidence$selected), 1L)
})

# ---- T3: a sabotaged stability at the argmin flips the verdict -------------

test_that("an unstable CV argmin flips K to the best passing rank", {
  cells <- karb_cells_env$cal_cells(karb_cells_env$CAL_ROOT_SEED)
  cl <- cells[[3L]]                             # planted tier K = 4
  honest <- bp_select_k(cl$counts, k_grid = 2:5, seed = 20261015L)
  expect_equal(honest$k, 4L)                   # argmin lands on the planted tier

  real_stability <- get(".bp_ensemble_stability", envir = asNamespace(
    "sigminer"))
  on.exit(assignInNamespace(".bp_ensemble_stability", real_stability,
    envir = asNamespace("sigminer")), add = TRUE)
  assignInNamespace(".bp_ensemble_stability", function(counts, k, replicates,
      seed) {
    if (identical(as.integer(k), 4L)) {
      # the argmin rank is declared unstable, exactly the S2 shape
      list(avg = 0.3, min = 0.05, per_signature = c(0.3, 0.05))
    } else {
      real_stability(counts, k, replicates, seed)
    }
  }, envir = asNamespace("sigminer"))

  flipped <- bp_select_k(cl$counts, k_grid = 2:5, seed = 20261015L)
  ev <- flipped$k_evidence
  expect_false(identical(flipped$k, 4L))
  sel <- ev[ev$selected, ]
  expect_true(sel$passes_veto)
  expect_false(sel$is_argmin)                  # not the (vetoed) argmin row
  # and it is the best passing rank in CV order (fallback-scan semantics)
  cand <- ev[ev$passes_veto & !ev$is_argmin, ]
  ord <- order(cand$cv_test_deviance, cand$k, na.last = TRUE)
  expect_equal(flipped$k, cand$k[ord[1L]])
})

# ---- T4: diagnostics are observable on the returned object ------------------

test_that("the returned evidence table carries the full diagnostic schema", {
  cells <- karb_cells_env$cal_cells(karb_cells_env$CAL_ROOT_SEED)
  res <- bp_select_k(cells[[1L]]$counts, k_grid = 2:5, seed = 20261015L)
  ev <- res$k_evidence
  want <- c("k", "cv_test_deviance", "cv_test_mse", "cv_train_mse",
    "fold_test_deviance", "n_seeds_converged", "avg_stability",
    "min_cluster_stability", "per_signature_stability",
    "wilcoxon_p_vs_prev", "wilcoxon_l2_median_delta", "passes_veto",
    "is_argmin", "selected", "veto_reason")
  expect_equal(names(ev), want)
  expect_type(ev$veto_reason, "character")
  expect_type(ev$selected, "logical")
  expect_is.list(ev$fold_test_deviance)
  expect_equal(vapply(ev$fold_test_deviance, length, integer(1L)),
    rep(8L, 4L))
  expect_is.list(ev$per_signature_stability)
})

# ---- T5: cv = FALSE demotes to the veto-only scan ----------------------------

test_that("cv = FALSE arbitrates on the stability veto alone", {
  cells <- karb_cells_env$cal_cells(karb_cells_env$CAL_ROOT_SEED)
  res <- bp_select_k(cells[[1L]]$counts, k_grid = 2:5, cv = FALSE,
    seed = 20261015L)
  ev <- res$k_evidence
  expect_false(any(ev$is_argmin))              # layer 1 absent -> mark erased
  expect_true(all(is.na(ev$cv_test_deviance))) # CV columns demoted to NA
  expect_true(all(vapply(ev$fold_test_deviance, function(z)
    length(z) == 8L && all(is.na(z)))))
  cand <- ev[ev$passes_veto, ]
  expect_true(nrow(cand) > 0L)
  expect_equal(res$k, min(cand$k))             # lowest passing rank wins
})

# ---- T6: input validation fails fast at the boundary ------------------------

test_that("bp_select_k validates its inputs with named messages", {
  good <- karb_cells_env$cal_cells(karb_cells_env$CAL_ROOT_SEED)[[1L]]$counts
  expect_error(bp_select_k(matrix(1:4, 2), k_grid = 2:3, seed = 1L),
    "numeric matrix")
  bad <- good
  bad[1L, 1L] <- -1
  expect_error(bp_select_k(bad, k_grid = 2:3, seed = 1L), "non-negative")
  expect_error(bp_select_k(good, k_grid = c(2L, 0L), seed = 1L),
    "positive whole numbers")
  expect_error(bp_select_k(good, k_grid = 2:100L, seed = 1L),
    "must not exceed")
  expect_error(bp_select_k(good, k_grid = 2:3, seed = 1L, cv = NA),
    "single TRUE or FALSE")
  expect_error(bp_select_k(good, k_grid = 2:3, seed = 1.5), "whole number")
  expect_error(bp_select_k(good, k_grid = 2:3, seed = 1L, replicates = 1L),
    "replicates")
  expect_error(bp_select_k(good, k_grid = 2:3, seed = 1L, k_folds = 200L),
    "channel count")
})

# ---- T7: opt-in hook on the bp route, legacy preserved ---------------------

test_that("k_select = 'arbitrate' annotates the route; the default is untouched", {
  cells <- karb_cells_env$cal_cells(karb_cells_env$CAL_ROOT_SEED)
  nmf_matrix <- t(cells[[2L]]$counts)         # samples x components for the route
  legacy <- bp_extract_signatures(nmf_matrix, range = 2:3, n_bootstrap = 0L,
    n_nmf_run = 2L, cores = 1L, cores_solution = 1L, seed = 123456L)
  arbit <- bp_extract_signatures(nmf_matrix, range = 2:3, n_bootstrap = 0L,
    n_nmf_run = 2L, cores = 1L, cores_solution = 1L, seed = 123456L,
    k_select = "arbitrate")
  # Legacy preservation: every slot the legacy result has (except the
  # intentionally re-arbitrated `suggested`) is byte-identical.
  for (slot in setdiff(names(legacy), "suggested")) {
    expect_identical(serialize(legacy[[slot]], connection = NULL),
      serialize(arbit[[slot]], connection = NULL), info = slot)
  }
  expect_equal(arbit$suggested_legacy, legacy$suggested)
  expect_type(arbit$suggested, "integer")
  ev <- arbit$k_evidence
  expect_true(is.data.frame(ev))
  expect_equal(ev$k[ev$selected], arbit$suggested)
  expect_true(any(ev$selected))
})
