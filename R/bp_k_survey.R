# bp_k_survey.R: the per-rank CV-survey engine feeding the bp K-arbitration
# core (parity plan P2). Semantic port of the msuiter ms_select_k() survey
# (msuiter R/kselect-select.R, main @ dad30440): the SUITOR-style fold CV
# curve (layer-1 evidence), the replicate-ensemble consensus stability
# (layer-2 evidence) and the diagnostic-only Wilcoxon columns (layer 3;
# the reference's deliberate divergence -- diagnostics never arbitrate).
#
# Seed law (single root seed): every stochastic draw -- fold partition,
# NMF inits, ensemble replicates -- is derived from ONE root seed through
# .bp_derived_seed(); one root seed yields byte-identical verdicts and
# diagnostics, and there is no ambient-RNG dependence anywhere in the
# survey (every NMF call re-seeds internally from its derived stream).
#
# Divergences from the reference, recorded for the design digest: the
# reference kernel masks SUITOR cells inside its rust engine -- this port
# masks cells at the R level (train = masked-to-zero counts, held-out
# cells scored against the training basis through a fixed-round KL
# multiplicative projection); n_seeds_converged is computed here rather
# than left a placeholder; the stability scale is the documented
# agreement mapping 2C - 1 of the cophenetic co-assignment frequency.

# Survey budgets (named constants; the reference-side analogues are the
# ms_select_k defaults replicates/k_folds/n_seeds, cited as a reference,
# not recalibrated -- contract: changing these changes verdicts and
# requires re-running the agreement fixtures).
.BP_SURVEY_REPLICATES <- 4L        # ensemble fits per rank (stability)
.BP_SURVEY_FOLDS <- 8L             # SUITOR folds over the held-out cells
.BP_SURVEY_NSEEDS <- 2L            # best-of-n fits per (rank, fold)
.BP_SURVEY_PROJ_ITERS <- 40L       # KL-multiplicative projection rounds
.BP_SURVEY_NMF_METHOD <- "brunet"  # the house NMF fit method
.BP_K_TINY <- 1e-300               # log/quotient-domain guard

# The full contractual evidence schema (column order is contractual,
# mirroring the reference's .ms_k_evidence_columns): the 8 columns the
# arbitration core consumes (.BP_K_EVIDENCE_REQUIRED) plus the 7 survey
# diagnostic columns.
.BP_K_EVIDENCE_COLUMNS <- c("k", "cv_test_deviance", "cv_test_mse",
  "cv_train_mse", "fold_test_deviance", "n_seeds_converged",
  "avg_stability", "min_cluster_stability", "per_signature_stability",
  "wilcoxon_p_vs_prev", "wilcoxon_l2_median_delta", "passes_veto",
  "is_argmin", "selected", "veto_reason")

# One Park-Miller LCG (a = 16807, m = 2^31 - 1; mirror of the committed
# fixture builder tests/testthat/fixtures/bp-k-arb/synthetic_cells.R --
# all intermediates below 2^53, IEEE-exact, so streams are bit-replayable).
.bp_lcg_unit <- function(n, seed) {
  m <- 2147483647
  state <- abs(seed) %% (m - 1) + 1                # state in [1, m - 1]
  out <- numeric(n)
  for (i in seq_len(n)) {
    state <- (16807 * state) %% m
    out[i] <- state / m
  }
  out
}

# The single-root-seed stream derivation: fold each extra coordinate into
# the state with one LCG step, so every (root, purpose, rank, replicate)
# tuple owns an independent-but-reproducible stream. Character coordinates
# (the purpose labels) are folded through their character codes, so the
# derivation stays RNG-free and host-reproducible without a hash package.
.bp_derived_seed <- function(root, ...) {
  m <- 2147483647
  s <- abs(root) %% (m - 1) + 1
  coords <- list(...)
  for (x in coords) {
    if (is.character(x)) {
      codes <- as.numeric(utf8ToInt(paste0(x, collapse = "")))
      t <- sum(codes) * 131 + length(codes)
    } else {
      t <- as.numeric(x)
    }
    t <- abs(t) %% (m - 1) + 1
    s <- (16807 * s) %% m
    s <- (s + 16807 * t) %% m
  }
  as.integer(if (s < 1) 1 else s)
}

# Scalar whole-number gate for the survey budgets.
.bp_whole_scalar <- function(value, min, name) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value != floor(value) || value < min) {
    stop(sprintf("`%s` must be a single whole number >= %d", name, min),
      domain = NA)
  }
  as.integer(value)
}

# KL divergence sum(x log(x/mu) - x + mu) with the x = 0 convention (the
# log term drops out), optionally restricted to a cell mask. The survey's
# own deviance, so the CV curve is self-consistent across ranks without
# reaching into NMF internals.
.bp_kl_deviance <- function(x, mu, mask = NULL) {
  if (is.null(mask)) {
    mask <- matrix(TRUE, nrow = nrow(x), ncol = ncol(x))
  }
  mu <- pmax(mu, .BP_K_TINY)
  log.pos <- x > 0 & mask
  log.part <- x[log.pos] * log(x[log.pos] / mu[log.pos])
  sum(log.part) - sum(x[mask]) + sum(mu[mask])
}

# One deterministic NMF fit: the house route (check_nmf_matrix hygiene,
# method "brunet"), serial backend, single run -- pure in (data, k, seed).
# The information-level non-convergence chatter of the fit is silenced here
# (not swallowed): non-factorisations surface as a non-finite deviance at
# the convergence gate below.
.bp_nmf_fit <- function(x, k, seed) {
  fit <- suppressWarnings(NMF::nmf(check_nmf_matrix(x), rank = k,
    method = .BP_SURVEY_NMF_METHOD, nrun = 1L, seed = seed,
    .options = "nm", .pbackend = "serial"))
  if (!inherits(fit, "NMFfit")) {
    fit <- fit[[1L]]                              # nrun > 1 shapes
  }
  list(W = NMF::basis(fit), H = NMF::coef(fit))
}

# Best-of-n fit over derived per-seed streams: the factorisation with the
# smallest deviance on the observed cells (first wins on ties ->
# deterministic even against seed pathologies). n_converged is the honest
# convergence tally.
.bp_nmf_best_of_n <- function(x, k, n_seeds, seed_root, purpose,
                              mask = NULL) {
  best <- NULL
  best_dev <- Inf
  n_conv <- 0L
  for (i in seq_len(n_seeds)) {
    fit <- .bp_nmf_fit(x, k, .bp_derived_seed(seed_root, purpose, k, i))
    if (is.null(fit$W) || is.null(fit$H) ||
        !all(is.finite(fit$W)) || !all(is.finite(fit$H))) {
      next
    }
    mu <- fit$W %*% fit$H
    if (!all(is.finite(mu))) {
      next
    }
    dev <- .bp_kl_deviance(x, mu, mask)
    if (is.finite(dev) && dev < best_dev) {
      best_dev <- dev
      best <- fit
      n_conv <- n_conv + 1L
    }
  }
  list(fit = best, dev = if (is.null(best)) as.numeric(NA) else best_dev,
    n_converged = n_conv)
}

# SUITOR fold plan: a balanced, rank-independent deterministic partition of
# the COUNTS CELLS (the reference kernel holds out cells, "at least one
# row per column"), drawn from the root-seed stream; the permutation-cut
# construction makes every fold non-empty by design and keeps the fold
# curves of all ranks comparable (same partition across the grid).
.bp_cell_folds <- function(n_rows, n_cols, k_folds, root) {
  n <- n_rows * n_cols
  draws <- .bp_lcg_unit(n, .bp_derived_seed(root, "cells"))
  perm <- order(draws, seq_len(n))                # stable: index breaks ties
  fold <- integer(n)
  fold[perm] <- (seq_len(n) - 1L) %% k_folds + 1L
  matrix(fold, nrow = n_rows, ncol = n_cols)
}

# KL-multiplicative projection of the held-out cells onto the training
# basis: Lee-Seung KL update for H only. `mask` marks the held-out cells
# (TRUE); the update is driven by the observed cells (`observed = x *
# !mask`), whose reconstruction regularizes the denominator; held-out cells
# neither drive the update nor need their (masked) observation. Fixed
# round count, no RNG -> deterministic.
.bp_project_exposures <- function(x, W, mask) {
  k <- ncol(W)
  m <- ncol(x)
  observed <- x * (!mask)                          # held-out cells -> 0
  H <- matrix(rep(max(mean(observed), .BP_K_TINY), k * m), nrow = k,
    ncol = m)
  for (i in seq_len(.BP_SURVEY_PROJ_ITERS)) {
    WH <- W %*% H
    num <- t(W) %*% observed
    den <- t(W) %*% (WH * (!mask))
    H <- H * num / pmax(den, .BP_K_TINY)
  }
  H
}

# The fold-CV curve for one rank: the mean held-out KL deviance over the
# folds (layer-1 evidence), the mean training deviance of the fold fits
# and the per-fold held-out deviances. A fold whose fits all fail the
# convergence gate demotes its held-out entry to NA (finite-only curves
# downstream, matching the ported argmin law).
.bp_cv_rank_curve <- function(counts, k, foldmat, k_folds, root, n_seeds) {
  te <- rep(as.numeric(NA), k_folds)
  tr <- rep(as.numeric(NA), k_folds)
  conv <- 0L
  for (f in seq_len(k_folds)) {
    mask <- foldmat == f
    observed <- counts * (!mask)                  # train on observed cells
    fit <- .bp_nmf_best_of_n(observed, k, n_seeds, root, "cv", mask)
    conv <- conv + fit$n_converged
    if (is.null(fit$fit)) {
      next
    }
    H <- .bp_project_exposures(observed, fit$fit$W, mask)
    mu <- fit$fit$W %*% H
    te[f] <- .bp_kl_deviance(counts, mu, mask)
    tr[f] <- .bp_kl_deviance(observed, mu, !mask)
  }
  list(cv_test = if (all(is.na(te))) as.numeric(NA) else
      mean(te, na.rm = TRUE),
    cv_train = if (all(is.na(tr))) as.numeric(NA) else
      mean(tr, na.rm = TRUE),
    folds = te, converged = conv)
}

# Replicate-ensemble consensus stability for one rank: each replicate
# partition is the dominant-signature assignment of its fit; the
# cophenetic matrix C(i, j) is the across-replicate co-assignment
# frequency; average-linkage clustering on 1 - C cuts k clusters and each
# cluster enters as its mean pairwise agreement mapped to the [-1, 1]
# agreement scale (2C - 1): perfectly reproducible clusters score 1,
# anti-reproducible junk components score near -1 -- the reading the
# reference reports for avg_stability / min_cluster_stability /
# per_signature_stability. Also returns the per-sample L2 reconstruction
# residual of the ensemble's best full-data fit (layer-3 input).
.bp_ensemble_stability <- function(counts, k, replicates, seed) {
  m <- ncol(counts)
  C <- matrix(0, nrow = m, ncol = m)
  fitted <- 0L
  best <- NULL
  best_dev <- Inf
  for (r in seq_len(replicates)) {
    fit <- .bp_nmf_fit(counts, k, .bp_derived_seed(seed, "stab", k, r))
    if (is.null(fit$H) || anyNA(fit$H) || !all(is.finite(fit$H))) {
      next
    }
    a <- apply(fit$H, 2L, which.max)              # dominant signature
    if (anyNA(a)) {
      next
    }
    C <- C + outer(a, a, "==")
    fitted <- fitted + 1L
    dev <- .bp_kl_deviance(counts, fit$W %*% fit$H)
    if (is.finite(dev) && dev < best_dev) {
      best_dev <- dev
      best <- fit
    }
  }
  if (fitted < 2L) {
    # fewer than two usable partitions: an honest "unavailable" -- the
    # veto fails such ranks conservatively (stability-statistics law).
    return(list(avg = as.numeric(NA), min = as.numeric(NA),
      per_signature = numeric(0), resid = NULL))
  }
  C <- C / fitted
  dend <- stats::hclust(stats::as.dist(1 - C, diag = TRUE, upper = TRUE),
    method = "average")
  cl <- stats::cutree(dend, k = k)
  stabs <- vapply(split(seq_len(m), cl), function(idx) {
    2 * mean(C[idx, idx, drop = FALSE]) - 1
  }, numeric(1L))
  names(stabs) <- paste0("cluster", seq_along(stabs))
  resid <- NULL
  if (!is.null(best)) {
    mu <- best$W %*% best$H
    resid <- vapply(seq_len(m), function(j)
      sqrt(sum((counts[, j] - mu[, j])^2)), numeric(1L))
  }
  list(avg = mean(stabs), min = min(stabs), per_signature = stabs,
    resid = resid)
}

# Layer-3 diagnostics (diagnostic-only by the reference's ruling: the
# Wilcoxon layer never arbitrates). Paired signed-rank p of the per-sample
# L2 residuals against the previous grid rank; NA at the first rank and
# when the residual does not move (reference parity).
.bp_wilcoxon_pair <- function(resid_prev, resid) {
  if (is.null(resid_prev) || is.null(resid) ||
      length(resid_prev) != length(resid)) {
    return(as.numeric(NA))
  }
  d <- resid - resid_prev
  if (all(d == 0)) {
    return(as.numeric(NA))
  }
  suppressWarnings(stats::wilcox.test(d, correct = FALSE,
    exact = FALSE)$p.value)
}

# Boundary validation for the survey surface (fail fast, named messages).
.bp_survey_check <- function(counts, k_grid, replicates, k_folds, n_seeds,
                             cv, seed) {
  if (!is.matrix(counts) || !is.numeric(counts) || is.null(rownames(counts)) ||
      is.null(colnames(counts))) {
    stop("the K survey input `counts` must be a numeric matrix with row ",
      "and column names", domain = NA)
  }
  if (anyNA(counts) || !all(is.finite(counts)) || any(counts < 0)) {
    stop("the K survey `counts` must be finite and non-negative",
      domain = NA)
  }
  if (!is.numeric(k_grid) || length(k_grid) < 1L || anyNA(k_grid) ||
      !all(is.finite(k_grid)) || !all(k_grid == floor(k_grid)) ||
      any(k_grid < 1)) {
    stop("k_grid must be positive whole numbers", domain = NA)
  }
  bound <- min(nrow(counts), ncol(counts))
  if (any(k_grid > bound)) {
    stop(sprintf("every k in k_grid must not exceed min(channels, samples) = %d",
      bound), domain = NA)
  }
  if (!is.logical(cv) || length(cv) != 1L || is.na(cv)) {
    stop("cv must be a single TRUE or FALSE", domain = NA)
  }
  seed <- .bp_whole_scalar(seed, 0, "seed")
  replicates <- .bp_whole_scalar(replicates, 2L, "replicates")
  k_folds <- .bp_whole_scalar(k_folds, 2L, "k_folds")
  n_seeds <- .bp_whole_scalar(n_seeds, 1L, "n_seeds")
  if (k_folds > nrow(counts)) {
    stop(sprintf("k_folds must not exceed the channel count (%d)",
      nrow(counts)), domain = NA)
  }
  invisible(list(seed = seed, replicates = replicates, k_folds = k_folds,
    n_seeds = n_seeds))
}

# The survey over one grid: builds the full contractual evidence frame
# (column order .BP_K_EVIDENCE_COLUMNS). Pure over (counts, grid, budget,
# root seed); the stochastic layers all draw from that one seed.
.bp_k_survey <- function(counts, k_grid, replicates, k_folds, n_seeds,
                         cv, seed) {
  grid <- sort(unique(as.integer(k_grid)))
  n <- length(grid)
  n_cells <- nrow(counts) * ncol(counts)
  foldmat <- .bp_cell_folds(nrow(counts), ncol(counts), k_folds, seed)
  te <- rep(as.numeric(NA), n)
  tr <- rep(as.numeric(NA), n)
  folds_list <- vector("list", n)
  conv <- rep(0L, n)
  avg <- rep(as.numeric(NA), n)
  mn <- rep(as.numeric(NA), n)
  per_sig <- vector("list", n)
  wx_p <- rep(as.numeric(NA), n)
  wx_med <- rep(as.numeric(NA), n)
  prev_resid <- NULL
  for (i in seq_len(n)) {
    k <- grid[i]
    if (isTRUE(cv)) {
      curve <- .bp_cv_rank_curve(counts, k, foldmat, k_folds, seed, n_seeds)
      te[i] <- curve$cv_test
      tr[i] <- curve$cv_train
      folds_list[[i]] <- curve$folds
      conv[i] <- curve$converged
    } else {
      folds_list[[i]] <- rep(as.numeric(NA), k_folds)
    }
    stab <- .bp_ensemble_stability(counts, k, replicates, seed)
    avg[i] <- stab$avg
    mn[i] <- stab$min
    per_sig[[i]] <- stab$per_signature
    wx_med[i] <- if (is.null(stab$resid)) as.numeric(NA) else
      stats::median(stab$resid)
    wx_p[i] <- .bp_wilcoxon_pair(prev_resid, stab$resid)
    prev_resid <- stab$resid
  }
  ev <- data.frame(
    k = grid,
    cv_test_deviance = te,
    cv_test_mse = ifelse(is.na(te), as.numeric(NA), sqrt(te / n_cells)),
    cv_train_mse = ifelse(is.na(tr), as.numeric(NA),
      sqrt(tr / (n_cells * max(k_folds - 1L, 1L)))),
    stringsAsFactors = FALSE)
  ev$fold_test_deviance <- folds_list
  ev$n_seeds_converged <- conv
  ev$avg_stability <- avg
  ev$min_cluster_stability <- mn
  ev$per_signature_stability <- per_sig
  ev$wilcoxon_p_vs_prev <- wx_p
  ev$wilcoxon_l2_median_delta <- wx_med
  ev$veto_reason <- .bp_veto_reason(avg, mn)
  ev$passes_veto <- !nzchar(ev$veto_reason)
  ev$is_argmin <- rep(FALSE, n)
  if (isTRUE(cv)) {
    ai <- .bp_cv_argmin(te)
    if (length(ai)) {
      ev$is_argmin[ai] <- TRUE
    }
  }
  ev$selected <- rep(FALSE, n)
  ev[, .BP_K_EVIDENCE_COLUMNS]
}

#' Select the number of signatures for a channel-by-sample count matrix
#'
#' `bp_select_k()` runs the bp-route K survey over a rank grid and applies
#' the ported three-layer arbitration (parity plan P2, semantic port of
#' the msuiter `ms_select_k()` at main@dad30440): the SUITOR-style CV
#' argmin (layer 1, first minimum / lowest rank on ties; switchable off
#' with `cv = FALSE`), the consensus stability veto (layer 2, upstream
#' reference thresholds `avg >= 0.8`, `min_cluster >= 0.2`, sum `>= 1.0`)
#' and diagnostic-only Wilcoxon columns (layer 3 -- they never
#' arbitrate). All stochastic layers derive from the single `seed`
#' (seed-law: same seed, byte-identical verdict and diagnostics).
#'
#' @param counts a numeric channel-by-sample count matrix with row and
#'   column names, finite and non-negative (the extraction catalogue).
#' @param k_grid positive whole-number ranks to evaluate, each no larger
#'   than `min(nrow(counts), ncol(counts))`; evaluated ascending.
#' @param replicates single whole number `>= 2`: ensemble fits per rank
#'   feeding the consensus stability.
#' @param k_folds single whole number in `[2, nrow(counts)]`: SUITOR
#'   folds over the held-out cells.
#' @param n_seeds single whole number `>= 1`: best-of-n fits per fold.
#' @param cv single logical: run layer 1? `FALSE` demotes the CV columns
#'   to NA and arbitrates on the stability veto alone.
#' @param seed single whole number `>= 0`: the root seed of every
#'   derived stream in the survey.
#'
#' @return A named list with `k` (the selected rank, single positive
#'   integer) and `k_evidence` (the full-surface evidence table, one row
#'   per grid rank, column order is the contractual .BP_K_EVIDENCE_COLUMNS). An
#'   undecidable grid (no rank passes the veto) raises with the
#'   undecidable diagnosis -- the pipeline-integration analogue degrades
#'   instead (see `bp_extract_signatures(k_select = "arbitrate")`).
#'
#' @examples
#' # A planted two-signature catalogue (synthetic, deterministic): two
#' # blocks of channels with anti-phased sample loadings. The arbitration
#' # reads the planted K = 2 off the CV curve and the stability veto.
#' # The budget arguments are reduced below the production defaults so the
#' # example stays within CRAN runtime limits; the seeded verdict is
#' # unchanged under the full default budget.
#' m <- matrix(1, 60, 20,
#'   dimnames = list(paste0("ch", 1:60), paste0("smp", 1:20)))
#' m[1:6, 1:10] <- 40
#' m[31:36, 11:20] <- 40
#' m[31:36, 1:10] <- 8
#' m[1:6, 11:20] <- 8
#' res <- bp_select_k(m, k_grid = 2:3, seed = 7L,
#'   replicates = 2L, k_folds = 4L, n_seeds = 1L)
#' if (!identical(res$k, 2L)) stop("the survey must arbitrate K = 2")
#' if (!is.data.frame(res$k_evidence) || !any(res$k_evidence$selected))
#'   stop("the evidence table must carry exactly one selected row")
#'
#' @seealso `.bp_k_arbitrate()` for the private arbitration core,
#'   [bp_extract_signatures()] for the opt-in route integration.
#' @export
bp_select_k <- function(counts, k_grid = 2:5,
                        replicates = .BP_SURVEY_REPLICATES,
                        k_folds = .BP_SURVEY_FOLDS,
                        n_seeds = .BP_SURVEY_NSEEDS,
                        cv = TRUE, seed = 1L) {
  budget <- .bp_survey_check(counts, k_grid, replicates, k_folds, n_seeds,
    cv, seed)
  ev <- .bp_k_survey(counts, k_grid, budget$replicates, budget$k_folds,
    budget$n_seeds, cv, budget$seed)
  ev <- .bp_k_arbitrate(ev, cv = cv, mode = "explicit")
  list(k = as.integer(ev$k[ev$selected][1L]), k_evidence = ev)
}

# The route-side wrapper (called from bp_extract_signatures when
# k_select = "arbitrate"): survey + arbitration under the reference's
# default-pipeline DEGRADE semantics -- the annotation never errors the
# extraction (the standalone surface bp_select_k() raises instead).
.bp_route_arbitrate <- function(counts, grid, seed) {
  ev <- .bp_k_survey(counts, grid, .BP_SURVEY_REPLICATES, .BP_SURVEY_FOLDS,
    .BP_SURVEY_NSEEDS, cv = TRUE, seed = seed)
  .bp_k_arbitrate(ev, cv = TRUE, mode = "degrade")
}
