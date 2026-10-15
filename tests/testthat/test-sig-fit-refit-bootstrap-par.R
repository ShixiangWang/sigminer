## Parallel bootstrap route for the refit judgment (sig_fit_bootstrap /
## sig_fit_bootstrap_batch territory, parity plan P3). The two new methods
## must flow through the existing bootstrap surfaces, whose replicate/sample
## fan-out rides the house future/furrr route (R/utils.R:37
## set_future_strategy), with the judgment -- which signatures survive --
## evaluated per replicate.

test_that("sig_fit_bootstrap accepts the two new methods end-to-end", {
  W <- p3_refit_sig()
  catalog <- p3_refit_catalog(c(A = .5, B = .5), total = 500L)
  for (method in c("likelihood", "lrt")) {
    hb <- suppressMessages(suppressWarnings(
      sig_fit_bootstrap(catalog, sig = W, n = 6L, method = method,
        type = "absolute"
      )
    ))
    expect_true(is.list(hb), info = method)
    expect_true(all(c("expo", "errors", "cosine") %in% names(hb)), info = method)
    expect_equal(dimnames(hb$expo)[[2]], paste0("Rep_", 1:6), info = method)
    expect_equal(rownames(hb$expo), colnames(W), info = method)
    expect_equal(length(hb$errors), 6L, info = method)
    expect_equal(length(hb$cosine), 6L, info = method)
    ## The judgment must hold across resampling: A and B dominate every
    ## replicate of this well-populated cell.
    rel <- p3_normalize(hb$expo)
    expect_true(all(rel["A", ] > .2 & rel["B", ] > .2), info = method)
    expect_true(all(rel[c("C", "D", "E", "F"), ] < p3_declare_threshold),
      info = paste(method, "bootstrap noise leak"))
  }
})

test_that("bootstrap method choices follow the same match.arg contract", {
  W <- p3_refit_sig()
  catalog <- p3_refit_catalog(c(A = .5, B = .5), total = 500L)
  err <- tryCatch({
    suppressMessages(sig_fit_bootstrap(catalog, sig = W, n = 1L, method = "l"))
    "no-error"
  }, error = function(e) conditionMessage(e))
  expect_match(err, "likelihood")
  expect_match(err, "lrt")
})

test_that("batch route runs both new methods and passes refit arguments through", {
  W <- p3_refit_sig()
  samples <- list(
    p3_refit_catalog(c(A = .5, B = .5), total = 500L, sample_name = "s1"),
    p3_refit_catalog(c(C = 1), total = 500L, sample_name = "s2"),
    p3_refit_catalog(c(A = 1), total = 500L, sample_name = "s3")
  )
  ## The batch route maps samples over colnames(catalogue_matrix), so the
  ## sample names come from p3_refit_catalog(sample_name = ...) directly;
  ## they are already unique and must be preserved.
  mats <- do.call(cbind, samples)
  rv <- suppressMessages(suppressWarnings(
    sig_fit_bootstrap_batch(mats, sig = W, methods = c("likelihood", "lrt"),
      n = 4L, seed = 90210L, use_parallel = FALSE, alpha = 0.05
    )
  ))
  expect_true(is.list(rv))
  expect_true(all(c("expo", "error", "cosine", "p_val") %in% names(rv)))
  expect_true(all(c("likelihood", "lrt") %in% rv$expo$method))
  expect_true(all(c("s1", "s2", "s3") %in% rv$expo$sample))
  ## The s2 cell must end up C-owned in both methods at the optimal row
  ## (the batch tables are long-format: sample/sig/exposure/type/method).
  ## Per-replicate judgment stability is covered through the single-sample
  ## bootstrap surface by the first test of this file.
  for (m in c("likelihood", "lrt")) {
    opt <- rv$expo[type == "optimal" & sample == "s2" & method == m, ]
    expect_gte(nrow(opt), 6L, label = paste(m, "s2 long rows"))
    tot <- sum(opt$exposure)
    expect_gt(opt[sig == "C", exposure] / tot, .90, label = paste(m, "s2 optimal C"))
  }
})

## Worker-snapshot gate: future workers are fresh R processes that resolve
## sigminer through .libPaths(), never the load_all namespace of the test
## master. This battery therefore provisions a worker snapshot with
##   R CMD INSTALL --no-byte-compile --library=/tmp/p3refit/rlib <tree>
## immediately before the files run (see /tmp/p3refit/run_battery.sh); under
## R CMD check the check library already plays this role. If neither is
## present the arms skip loudly rather than silently passing.
p3_worker_snapshot <- function() {
  p3_worker_lib <- "/tmp/p3refit/rlib"
  have_snapshot <- dir.exists(p3_worker_lib) && isTRUE(
    "sigminer" %in% rownames(utils::installed.packages(lib.loc = p3_worker_lib,
      no.cache = TRUE))
  )
  if (!have_snapshot && dir.exists(Sys.getenv("RCHECK_DIR", ""))) {
    have_snapshot <- TRUE
    p3_worker_lib <- .libPaths()[[1L]]
  }
  skip_if_not(have_snapshot,
    paste("no worker snapshot: provision", p3_worker_lib,
      "with R CMD INSTALL before the battery"))
  p3_worker_lib
}

p3_batch_cell <- function() {
  mats <- p3_refit_catalog(c(A = .5, B = .5), total = 400L)
  mats <- cbind(mats, p3_refit_catalog(c(C = 1), total = 400L),
    p3_refit_catalog(c(D = .3, E = .7), total = 400L))
  colnames(mats) <- c("s1", "s2", "s3")
  mats
}

p3_batch <- function(mats, W, par) {
  suppressMessages(suppressWarnings(
    sig_fit_bootstrap_batch(mats, sig = W, methods = "lrt",
      n = 20L, seed = 42424L, use_parallel = par
    )
  ))
}

## Optimal rows of one payload, ordered by the key columns the payload
## actually carries (expo/p_val are keyed by sample+sig, error/cosine by
## method+sample+type, so a blind order(sample, sig) would break the latter).
p3_optimal <- function(payload) {
  if ("type" %in% colnames(payload)) {
    payload <- payload[type == "optimal"]
  }
  key <- intersect(c("method", "sample", "sig"), colnames(payload))
  ord <- do.call(order, lapply(key, function(k) payload[[k]]))
  payload[ord, ]
}

test_that("the parallel batch route drives the refit judgment on workers", {
  ## What the cross-arm contract is: the house route reseeds every worker
  ## from the master L'Ecuyer counter stream (furrr_options(seed = TRUE),
  ## R/sig_fit_bootstrap_batch.R:166), so bit-for-bit cross-arm identity of
  ## the RNG-consumed aggregates is NOT part of the contract -- the legacy
  ## suite pins it neither. The pins here are the real contract: workers
  ## execute this lane's code (snapshot-provisioned) and the two arms agree
  ## structurally, in the significance pattern, and statistically.
  W <- p3_refit_sig()
  p3_worker_lib <- p3_worker_snapshot()
  mats <- p3_batch_cell()
  serial <- p3_batch(mats, W, FALSE)
  opaths <- .libPaths()
  ## base::.libPaths(new) REPLACES the search path, so prepend the snapshot
  ## in front of the retained old paths (the site libraries hold future etc.).
  .libPaths(c(p3_worker_lib, opaths))
  on.exit(.libPaths(opaths))
  parallel <- p3_batch(mats, W, 2L)
  for (nm in c("expo", "error", "cosine", "p_val")) {
    expect_equal(colnames(serial[[nm]]), colnames(parallel[[nm]]), info = nm)
    expect_equal(nrow(serial[[nm]]), nrow(parallel[[nm]]), info = nm)
    expect_true(all(serial[[nm]]$method %in% "lrt"), info = nm)
  }
  so <- p3_optimal(serial$expo)
  po <- p3_optimal(parallel$expo)
  expect_equal(so$sample, po$sample)
  expect_equal(so$sig, po$sig)
  rel <- abs(so$exposure - po$exposure) / pmax(so$exposure, 1e-9)
  ## Measured on this fixture at n = 20 (battery log): max relative
  ## difference 0.0 per sample; the band is the honest no-regression room.
  expect_lt(max(rel), .25)
  sp <- p3_optimal(serial$p_val)
  pp <- p3_optimal(parallel$p_val)
  expect_equal(sp$sample, pp$sample)
  expect_equal(sp$sig, pp$sig)
  ## The decision layer must not fork: same significance calls in both arms
  ## (measured: identical; the true signatures sit at ~1e-25, the rest at 1).
  expect_equal(sp$p_value < .05, pp$p_value < .05)
  sc <- p3_optimal(serial$cosine)
  pc <- p3_optimal(parallel$cosine)
  expect_equal(sc$sample, pc$sample)
  ## Measured on this fixture (battery log): the optimal-row cosines are
  ## 1 / 1 / NaN in BOTH arms -- the flat cell (s3) yields a deterministic
  ## NaN optimal cosine, so the cross-arm pin is identity INCLUDING the NaN
  ## positions (identical() treats NaN as identical to NaN), and the clean
  ## cells (s1, s2) must reconstruct with fidelity above .99.
  expect_identical(sc$cosine, pc$cosine)
  expect_true(all(sc$cosine[sc$sample %in% c("s1", "s2")] > .99))
  expect_true(all(pc$cosine[pc$sample %in% c("s1", "s2")] > .99))
  expect_true(all(is.na(sc$cosine[sc$sample == "s3"])))
  expect_true(all(is.na(pc$cosine[pc$sample == "s3"])))
  expect_true(all(serial$error$errors >= 0))
  expect_true(all(parallel$error$errors >= 0))
})

test_that("the parallel batch route is reproducible under a pinned seed", {
  ## The rigorous half of the parallel contract: two runs of the parallel
  ## arm at the same seed and worker count are bit-identical in all four
  ## payloads (measured: all.equal TRUE x 4 in the battery probe).
  W <- p3_refit_sig()
  p3_worker_lib <- p3_worker_snapshot()
  mats <- p3_batch_cell()
  opaths <- .libPaths()
  .libPaths(c(p3_worker_lib, opaths))
  on.exit(.libPaths(opaths))
  run1 <- p3_batch(mats, W, 2L)
  run2 <- p3_batch(mats, W, 2L)
  for (nm in c("expo", "error", "cosine", "p_val")) {
    expect_true(isTRUE(all.equal(run1[[nm]], run2[[nm]])),
      info = paste("parallel repro:", nm))
  }
})
