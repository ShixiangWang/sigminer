# generate_verdicts.R -- provenance generator for the bp K-arbitration
# fixtures (parity plan P2). Run on the controller host to regenerate the
# two .rds files next to this script:
#
#   Rscript generate_verdicts.R <fixture-dir>
#
# NEVER run it from R CMD check / CI: it loads the msuiter reference via
# pkgload::load_all from a checkout path (the package is not on CRAN and
# must never appear in DESCRIPTION -- sigminer CI convention: the lockfile
# hard-fails on non-CRAN entries). The .rds payloads are plain numbers;
# the testthat units consume them with zero msuiter presence.
#
# Data law (D3): every numeric here is SYNTHETIC. The calibration cells
# come from the deterministic RNG-free builder below (Park-Miller LCG on a
# pinned root seed); no COSMIC or real-data-derived numbers are involved.
#
# Reference implementation: msuiter main@dad30440, R/kselect-select.R
# (.ms_k_arbitrate / .ms_veto_reason / ms_select_k). Reference host facts:
# R 4.5.2, controller host, captured 2026-10-15; the original recordings
# live in the lane logs /tmp/p2karb/logs/03-*.log and 04-*.log.

MSUITER_CHECKOUT <- "/Users/wsx/Documents/GitHub/msuiter"
ROOT_SEED <- 20261015L           # the single root seed of the whole campaign
K_GRID <- 2:5L                   # evaluated grid per calibration cell

# ---- calibration cells: the committed deterministic builder ----------------
# synthetic_cells.R next to this script is the single source of truth for
# the cells (the testthat units source it too, so the tested cells and the
# recorded reference cells are the same numbers by construction). Under the
# Rscript wrapper the script path arrives as the "--file=..." element of
# the full command line.
.cmdline <- commandArgs(trailing = FALSE)
.file.arg <- grep("^--file=", .cmdline)
if (!length(.file.arg)) {
  stop("cannot locate this script via --file=; invoke as 'Rscript generate_verdicts.R <fixture-dir>'", domain = NA)
}
.own.dir <- dirname(normalizePath(sub("^--file=", "", .cmdline[[.file.arg[1L]]])))
source(file.path(.own.dir, "synthetic_cells.R"))

# ---- the reference runs ---------------------------------------------------

main <- function(outdir) {
  suppressMessages(pkgload::load_all(MSUITER_CHECKOUT,
    export_all = FALSE, attach_testthat = FALSE))
  arb <- get(".ms_k_arbitrate", envir = asNamespace("msuiter"))
  veto <- get(".ms_veto_reason", envir = asNamespace("msuiter"))

  # (1) core-level cases: hand-built evidence frames through the reference
  # arbitration core, recorded with their run parameters for replay.
  mk_ev <- function(k, cv, avg, mn) {
    vr <- mapply(veto, avg, mn, USE.NAMES = FALSE)
    data.frame(
      k = as.integer(k),
      cv_test_deviance = cv,
      avg_stability = avg,
      min_cluster_stability = mn,
      passes_veto = vr == "",
      veto_reason = vr,
      is_argmin = rep(FALSE, length(k)),
      selected = rep(FALSE, length(k)),
      stringsAsFactors = FALSE)
  }
  with_argmin <- function(ev, k) {
    ev$is_argmin[ev$k == k] <- TRUE
    ev
  }
  cases <- list()
  add <- function(name, ev, cv, mode) {
    warned <- FALSE
    explicit <- identical(mode, "explicit")
    run <- function() arb(ev, cv = cv, explicit = explicit)
    got <- tryCatch(
      if (explicit) run() else withCallingHandlers(
        warning = function(w) { warned <<- TRUE; NULL },
        expr = run()),
      error = function(e) structure(list(), class = class(e)))
    is_err <- inherits(got, "error") && !is.data.frame(got)
    cases[[name]] <<- list(
      ev = ev, cv = cv, mode = mode,
      expect_error = is_err,
      error_class = if (is_err) class(got)[[1L]] else NA,
      warned = warned,
      selected = if (is_err) NULL else got$selected,
      is_argmin_out = if (is_err) NULL else got$is_argmin)
    cat(name, if (is_err) "-> reference error (undecidable)" else
      sprintf("-> reference selected k = %d",
        got$k[which(got$selected)[1L]]),
      if (warned) " [warned]" else "", "\n")
  }
  # the seven pinned scenarios (S1..S7) of the parity design
  add("S1_argmin_passes",
    with_argmin(mk_ev(2:5, c(300, 240, 260, 280), c(.9, .95, .85, .82),
      c(.35, .4, .3, .25)), 3L), TRUE, "explicit")
  add("S2_argmin_vetoed_fallback",
    with_argmin(mk_ev(2:5, c(300, 240, 260, 280), c(.9, .5, .85, .82),
      c(.35, .1, .3, .25)), 3L), TRUE, "explicit")
  add("S3_cv_off_lowest_passing",
    mk_ev(2:5, c(300, 240, 260, 280), c(.5, .5, .85, .82),
      c(.1, .1, .3, .25)), FALSE, "explicit")
  add("S4_none_passes_explicit_error",
    with_argmin(mk_ev(2:5, c(300, 240, 260, 280), rep(.4, 4), rep(.1, 4)),
      3L), TRUE, "explicit")
  add("S5_none_passes_degrade_to_argmin",
    with_argmin(mk_ev(2:5, c(300, 240, 260, 280), rep(.4, 4), rep(.1, 4)),
      3L), TRUE, "degrade")
  add("S6_cv_tie_lowest_rank",
    with_argmin(mk_ev(2:5, c(300, 240, 240, 280), c(.9, .95, .9, .82),
      c(.35, .4, .35, .25)), 3L), TRUE, "explicit")
  add("S7_vetoed_argmin_tie_fallback",
    with_argmin(mk_ev(2:5, c(240, 240, 260, 280), c(.4, .95, .9, .82),
      c(.1, .4, .35, .25)), 2L), TRUE, "explicit")
  saveRDS(list(root_seed = ROOT_SEED, cases = cases),
    file.path(outdir, "msuiter_core_cases.rds"))

  # (2) pipeline-level verdicts: the full reference arbiter over the
  # three-tier calibration cells.
  verdicts <- list()
  for (cell in cal_cells()) {
    cat(cell$name, "planted_k =", cell$planted_k,
      "dim:", nrow(cell$counts), "x", ncol(cell$counts), "\n")
    cat1 <- ms_catalog(cell$counts,
      list(name = paste0("SYN", cell$name), labels = rownames(cell$counts)),
      colnames(cell$counts), provenance = list(genome = "GRCh38"))
    res <- ms_select_k(cat1, k_grid = K_GRID, seed = ROOT_SEED)
    cat("  reference verdict k =", res$k, "\n")
    verdicts[[cell$name]] <- list(name = cell$name, planted_k = cell$planted_k,
      k = res$k, dim = dim(cell$counts),
      evidence = res$k_evidence[, c("k", "cv_test_deviance", "avg_stability",
        "min_cluster_stability", "passes_veto", "is_argmin", "selected")])
  }
  saveRDS(list(root_seed = ROOT_SEED, k_grid = K_GRID, cells = verdicts),
    file.path(outdir, "msuiter_verdicts.rds"))
  cat("fixtures written to ", outdir, "\n", sep = "")
}

args <- commandArgs(trailing = TRUE)
if (!length(args)) {
  stop("usage: Rscript generate_verdicts.R <fixture-dir>", domain = NA)
}
main(args[[1L]])
