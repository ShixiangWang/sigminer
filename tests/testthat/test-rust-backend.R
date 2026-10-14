# Rust-kernel backend tests for the P1 boundary (parity plan P1, closes #478).
#
# Product surface under test: the internal R/rust_backend.R layer that drives
# the certified msuiter KL-MU kernel through the frozen protocol-1 CLI
# contract, plus the backend = c("R", "rust") switches on sig_extract(),
# sig_estimate() and bp_extract_signatures().
#
# The single system2() contact point is .kernel_exec(); every scripted unit
# replaces it via assignInNamespace and restores it on exit (house idiom,
# see test-sigprofiler-roundtrip.R). Real-kernel units are env-gated on
# MSUITER_TEST_KERNEL_BIN and the bit-equivalence gate is additionally
# NMF-version-gated (>= 1.5.0).
#
# Floor law (DESCRIPTION: R >= 3.5): the base null-default operator is a
# 4.4.0 export only, so it must never appear in this battery or in the code
# it tests; exit-status handling goes through .rust_status().
#
# D3: every fixture here is synthetic, powers-of-two idiom (no COSMIC-derived
# numerics anywhere).

# ---- Shared fixture helpers ---------------------------------------------

# Samples x features orientation handed to sig_extract()/sig_estimate()
# (they transpose internally). Powers-of-two entries keep every derived
# quantity binary-exact.
rb_nmf_matrix <- function() {
  matrix(
    c(
      16, 8, 4, 2,   # Sample1 across the four features
      8, 16, 2, 4,   # Sample2
      4, 2, 16, 8,   # Sample3
      2, 4, 8, 16    # Sample4
    ),
    nrow = 4L, ncol = 4L,
    dimnames = list(paste0("Sample", 1:4),
                    c("A[C]A", "A[C]C", "A[T]A", "A[T]C"))
  )
}

# Kernel-orientation matrix (m rows x n cols) for direct .rust_fit_kl() use.
rb_live_v <- function(zero_column = FALSE) {
  rows <- list(
    c(16, 8, 4, 2, 1, 0, 4, 8),
    c( 8, 16, 2, 4, 0, 1, 8, 2),
    c( 4, 2, 16, 8, 2, 4, 0, 1),
    c( 2, 4, 8, 16, 4, 2, 1, 0),
    c( 1, 0, 2, 4, 16, 8, 2, 4),
    c( 0, 1, 4, 2, 8, 16, 1, 2)
  )
  v <- t(matrix(as.numeric(unlist(rows)), ncol = 6L, byrow = TRUE))
  if (zero_column) v[, 6L] <- 0
  v
}

rb_wired_attrs <- function(res) {
  attributes(res)[c("nrun", "method", "seed", "call_method")]
}

rb_catch_msg <- function(expr) {
  tryCatch({
    force(expr)
    NA_character_
  }, error = function(e) conditionMessage(e))
}

rb_expect_fails <- function(expr, pattern, fixed = FALSE) {
  msg <- rb_catch_msg(expr)
  expect_false(is.na(msg), info = "expected an error but the call succeeded")
  expect_true(grepl(pattern, msg, fixed = fixed),
               info = paste0("error message did not contain [", pattern, "]: ", msg))
  invisible(msg)
}

# ---- Mocking layer -------------------------------------------------------

# The frozen --version handshake line, verified live against the trunk binary
# (msuiter @ dad30440c2):
#   $ msuiter-kernel --version
#   msuiter-kernel protocol=1 engine="engine: pure numeric kernels, FFI-free" kl_eps=1e-12
rb_version_line <- paste0(
  'msuiter-kernel protocol=1 engine="engine: pure numeric kernels, FFI-free" kl_eps=1e-12'
)

rb_kernel_result <- function(status = 0L, stdout = character(), stderr = character()) {
  res <- 0L
  attr(res, "status") <- as.integer(status)
  attr(res, "stdout") <- as.character(stdout)
  attr(res, "stderr") <- as.character(stderr)
  res
}

# Parse the single-case fit-kl flag surface into a named list.
rb_parse_argv <- function(args) {
  flags <- c("--in", "--out", "--m", "--n", "--k", "--iters",
             "--init", "--seed", "--replicate")
  req <- list()
  i <- 1L
  while (i <= length(args)) {
    f <- args[i]
    if (f %in% flags) {
      req[[sub("^--", "", f)]] <- args[i + 1L]
      i <- i + 2L
    } else {
      i <- i + 1L
    }
  }
  req
}

rb_matrix_rows <- function(m) {
  apply(m, 1, function(r) paste(sprintf("%.17e", r), collapse = ","))
}

# Well-formed output set for a fit-kl request: W (m x k) and H (k x n) as
# %.17e comma grids over powers-of-two values, strictly decreasing objective
# trace of iters + 1 rows (row 0 at the initializer), and a meta.csv echoing
# the request with the frozen key order.
rb_default_files <- function(req, obj_final = 8) {
  m <- as.integer(req$m)
  n <- as.integer(req$n)
  k <- as.integer(req$k)
  iters <- as.integer(req$iters)
  wmat <- matrix(as.numeric(outer(seq_len(m), seq_len(k), function(i, s) 2^((i + s) %% 4L))),
                 nrow = m, ncol = k, byrow = TRUE)
  hmat <- matrix(as.numeric(outer(seq_len(k), seq_len(n), function(s, j) 2^((s + j + 1L) %% 4L))),
                 nrow = k, ncol = n, byrow = TRUE)
  len <- iters + 1L
  obj <- obj_final + (len - seq_len(len)) * (obj_final + 1)
  trace <- c("t,objective", sprintf("%d,%.17e", seq_len(len) - 1L, obj))
  init_val <- if (is.null(req$init)) "external" else req$init
  seed_val <- if (is.null(req$seed)) "na" else req$seed
  rep_val <- if (identical(init_val, "seeded")) {
    if (is.null(req$replicate)) "0" else req$replicate
  } else {
    "na"
  }
  meta <- c(
    "protocol,1", "command,fit-kl",
    paste0("init,", init_val),
    paste0("seed,", seed_val),
    paste0("replicate,", rep_val),
    paste0("m,", req$m), paste0("n,", req$n), paste0("k,", req$k),
    paste0("iters,", req$iters),
    paste0("iterations,", req$iters),
    paste0("trace_len,", as.character(iters + 1L)),
    paste0("obj_init,", sprintf("%.17e", obj[1L])),
    paste0("obj_final,", sprintf("%.17e", obj[len])),
    "eps,1e-12",
    "engine_crate,engine: pure numeric kernels, FFI-free"
  )
  list(
    "W.csv" = rb_matrix_rows(wmat),
    "H.csv" = rb_matrix_rows(hmat),
    "trace.csv" = trace,
    "meta.csv" = meta
  )
}

# Scripted .kernel_exec() responder. Returns a state environment carrying
# $fn (the closure installed as .kernel_exec()), $calls (recorded invocations)
# and $run. tamper(req, files) may edit or remove the output files before
# they are written (fail-closed battery); obj_finals varies the per-run
# objective tail so best-of-nrun picks are observable; status/stderr script a
# non-zero exit.
rb_responder <- function(tamper = NULL, obj_finals = NULL,
                         status = 0L, stderr = character()) {
  state <- new.env()
  state$run <- 0L
  state$calls <- list()
  state$fn <- function(bin, args) {
    state$calls[[length(state$calls) + 1L]] <- list(bin = bin, args = args)
    if (length(args) >= 1L && identical(args[[1L]], "--version")) {
      return(rb_kernel_result(0L, stdout = rb_version_line))
    }
    req <- rb_parse_argv(args)
    state$run <- state$run + 1L
    of <- if (is.null(obj_finals)) 8 else obj_finals[[state$run]]
    files <- rb_default_files(req, obj_final = of)
    if (!is.null(tamper)) files <- tamper(req, files)
    if (!is.null(files) && length(files) > 0L) {
      dir.create(req[["out"]], recursive = TRUE, showWarnings = FALSE)
      for (nm in names(files)) {
        writeLines(files[[nm]], con = file.path(req[["out"]], nm))
      }
    }
    rb_kernel_result(status, stderr = stderr)
  }
  state
}

# Install a responder (state env from rb_responder) or a bare responder
# function as .kernel_exec(); the returned state env is also the cleanup
# token for rb_cleanup_mocks().
rb_mock_exec <- function(responder) {
  if (is.function(responder)) {
    fn <- responder
    responder <- new.env()
    responder$fn <- fn
    responder$calls <- list()
  }
  ns <- asNamespace("sigminer")
  responder$ns <- ns
  responder$original <- get(".kernel_exec", envir = ns)
  assignInNamespace(".kernel_exec", responder$fn, ns)
  responder
}

rb_mock_internal <- function(name, value) {
  state <- new.env()
  state$ns <- asNamespace("sigminer")
  state$name <- name
  state$had <- exists(name, envir = state$ns, inherits = FALSE)
  state$original <- if (state$had) get(name, envir = state$ns) else NULL
  assignInNamespace(name, value, state$ns)
  state
}

rb_cleanup_mocks <- function(...) {
  for (state in list(...)) {
    if (is.null(state)) next
    if (!is.null(state$name)) {
      if (isTRUE(state$had)) {
        assignInNamespace(state$name, state$original, state$ns)
      } else {
        tryCatch(rm(list = state$name, envir = state$ns), error = function(e) NULL)
      }
    } else if (!is.null(state$ns) && !is.null(state$original)) {
      assignInNamespace(".kernel_exec", state$original, state$ns)
    }
  }
}

# Pin the resolution to a private dummy executable file and neutralise the
# runtime gate on the optional 'msuiter' package; both restored on exit.
rb_use_dummy_bin <- function() {
  state <- new.env()
  state$bin <- tempfile("rb-kernel-")
  file.create(state$bin)
  state$old_opt <- getOption("sigminer.rust_kernel")
  options(sigminer.rust_kernel = state$bin)
  state
}

rb_release_dummy_bin <- function(state) {
  options(sigminer.rust_kernel = state$old_opt)
  unlink(state$bin)
}

# The recorded fit-kl invocations (handshake calls excluded) and one flag.
rb_fit_calls <- function(mock) {
  Filter(function(cl) identical(cl$args[[1L]], "fit-kl"), mock$calls)
}

rb_fit_values <- function(mock, flag) {
  cl <- rb_fit_calls(mock)
  if (length(cl) == 0L) return(character())
  sapply(cl, function(x) {
    v <- rb_parse_argv(x$args)[[flag]]
    if (is.null(v)) NA_character_ else v
  })
}

# Locate R/<file> in the package source tree, walking up from getwd() so the
# unit works under both test_dir (getwd = tests/testthat) and R CMD check
# (getwd = tests).
rb_repo_r_path <- function(fname) {
  d <- normalizePath(getwd())
  for (i in 1:5) {
    cand <- file.path(d, "R", fname)
    if (file.exists(cand)) return(cand)
    parent <- normalizePath(file.path(d, ".."))
    if (identical(parent, d)) break
    d <- parent
  }
  stop("cannot locate R/", fname, " upward from ", getwd(), domain = NA)
}

# Package root (the directory holding DESCRIPTION), located by the same
# upward walk as rb_repo_r_path().
rb_pkg_root <- function() {
  d <- normalizePath(getwd())
  for (i in 1:5) {
    if (file.exists(file.path(d, "DESCRIPTION"))) return(d)
    parent <- normalizePath(file.path(d, ".."))
    if (identical(parent, d)) break
    d <- parent
  }
  stop("cannot locate the package root upward from ", getwd(), domain = NA)
}

# bp's parallel fan-out serializes call_nmf into separate future workers
# that re-resolve the sigminer namespace from the library path -- under
# devtools::load_all they would load the installed CRAN build and miss the
# P1 internals. Install the lane source once into a private library under
# tempdir() and let the unit prepend it to .libPaths(): future workers
# copy the parent's paths verbatim (the same lever future itself uses for
# covr in its worker bootstrap). Returns the library path, or NULL when the
# install is unavailable (unit then skips honestly).
rb_lane_worker_lib <- function() {
  lib <- file.path(tempdir(), "rb-worker-lib")
  if (file.exists(file.path(lib, "sigminer", "DESCRIPTION"))) return(lib)
  rbin <- file.path(R.home(), "bin", "R")
  out_file <- tempfile("rb-install-out-")
  err_file <- tempfile("rb-install-err-")
  res <- system2(c(rbin, "CMD", "INSTALL", "--library", lib,
                   "--no-byte-compile", "--no-build-vignettes",
                   rb_pkg_root()),
                 stdout = out_file, stderr = err_file, wait = TRUE)
  if (.rust_status(res) != 0L) {
    warn_out <- character()
    if (file.exists(err_file)) warn_out <- readLines(err_file, warn = FALSE)
    warn_out <- utils::tail(warn_out, 5L)
    warning("worker-side R CMD INSTALL failed: ",
            paste(warn_out, collapse = " | "), call. = FALSE)
    return(NULL)
  }
  lib
}

rb_live_kernel_available <- function() {
  bin <- Sys.getenv("MSUITER_TEST_KERNEL_BIN")
  nzchar(bin) && file.exists(bin)
}

# Snapshot preservation units must replay under the fresh-session RNG
# configuration the snapshot was captured in: a preceding suite file can
# load a namespace that re-tunes the RNG kind, and set.seed() alone does
# not reset kind/normal.kind/sample.kind. Pins the R defaults and hands
# back the previous law for restoration on exit.
rb_pin_default_rng <- function() {
  old <- suppressWarnings(RNGkind())
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))
  old
}

# ---- (i) R-path preservation ---------------------------------------------
#
# The baseline snapshot in fixtures/rust-backend/rpath-baseline.rds was
# captured from the UNMODIFIED tree (commit e5975f13, NMF 0.28) before any
# backend wiring existed; see /tmp/p1backend/logs/01-baseline-snapshot.log.

test_that("backend = \"R\" reproduces the pre-wiring baseline snapshot byte-for-byte", {
  snap <- readRDS(file.path(test_path(), "fixtures", "rust-backend",
                             "rpath-baseline.rds"))
  nmf_now <- as.character(packageVersion("NMF"))
  if (!identical(nmf_now, snap$meta$nmf_version)) {
    skip(paste0("engine-version gate: snapshot captured with NMF ",
                snap$meta$nmf_version, ", this host runs ", nmf_now,
                " (NMF-internal draw order may legitimately differ)"))
  }
  m <- rb_nmf_matrix()
  old_rng <- rb_pin_default_rng()
  # this R build rejects the vector form of RNGkind(x); pass the three
  # components positionally (documented deviation, probed live).
  on.exit(suppressWarnings(
    RNGkind(old_rng[1L], old_rng[2L], old_rng[3L])
  ), add = TRUE)
  res <- suppressMessages(suppressWarnings(
    sig_extract(m, 2, nrun = 2, seed = 9L, method = "brunet", cores = 1)
  ))
  expect_identical(res$Raw$W, snap$W)
  expect_identical(res$Raw$H, snap$H)
  expect_identical(res$Signature, snap$Signature)
  expect_identical(res$Exposure, snap$Exposure)
  expect_equal(as.character(attr(res, "nrun")), "2")
  expect_equal(as.character(attr(res, "method")), "brunet")
  expect_equal(as.character(attr(res, "seed")), "9")
  expect_identical(attr(res, "call_method"), "NMF")
})

test_that("the R path is deterministic within a session (same seed, same result)", {
  m <- rb_nmf_matrix()
  run <- function() {
    suppressMessages(suppressWarnings(
      sig_extract(m, 2, nrun = 2, seed = 9L, method = "brunet", cores = 1)
    ))
  }
  r1 <- run()
  r2 <- run()
  expect_identical(r1$Raw$W, r2$Raw$W)
  expect_identical(r1$Raw$H, r2$Raw$H)
  expect_identical(r1$Signature, r2$Signature)
  expect_identical(rb_wired_attrs(r1), rb_wired_attrs(r2))
})

# ---- (a) argv staging and pre-exec validation -----------------------------

test_that("fit-kl argv is the golden frozen flag order for a known call", {
  mock <- rb_mock_exec(rb_responder())
  bin <- rb_use_dummy_bin()
  gate <- rb_mock_internal(".require_msuiter", function() invisible(TRUE))
  on.exit({
    rb_cleanup_mocks(mock, gate)
    rb_release_dummy_bin(bin)
  }, add = TRUE)
  opts <- getOption("sigminer.rust_iters")
  options(sigminer.rust_iters = 13L)
  on.exit(options(sigminer.rust_iters = opts), add = TRUE)

  v <- rb_live_v()[1:4, 1:4]
  fit <- .rust_fit_kl(V = v, k = 2L, iters = 13L, seed = 7L, init = "nndsvd")
  expect_equal(length(mock$calls), 2L)  # handshake + fit
  argv <- mock$calls[[2L]]$args
  expect_equal(mock$calls[[2L]]$bin, normalizePath(bin$bin))
  expect_equal(argv[[1L]], "fit-kl")
  flag_order <- argv[grep("^--", argv)]
  expect_identical(
    flag_order,
    c("--in", "--out", "--m", "--n", "--k", "--iters", "--init", "--seed")
  )
  req <- rb_parse_argv(argv)
  expect_identical(req$m, "4")
  expect_identical(req$n, "4")
  expect_identical(req$k, "2")
  expect_identical(req$iters, "13")
  expect_identical(req$init, "nndsvd")
  expect_identical(req$seed, "7")
  expect_null(req$replicate)
  expect_true(grepl("sigminer-rust-", req[["in"]], fixed = TRUE))
  expect_true(grepl("sigminer-rust-", req[["out"]], fixed = TRUE))
  expect_true(is.matrix(fit$W))
})

test_that("hostile V entries are rejected before any exec (no spawn, no staging contact)", {
  mock <- rb_mock_exec(rb_responder())
  bin <- rb_use_dummy_bin()
  gate <- rb_mock_internal(".require_msuiter", function() invisible(TRUE))
  on.exit({
    rb_cleanup_mocks(mock, gate)
    rb_release_dummy_bin(bin)
  }, add = TRUE)

  good <- rb_live_v()[1:4, 1:4]
  bad <- list(
    "NaN entry" = { x <- good; x[2, 3] <- NaN; x },
    "Inf entry" = { x <- good; x[1, 1] <- Inf; x },
    "negative entry" = { x <- good; x[3, 2] <- -4; x },
    "NA entry" = { x <- good; x[4, 4] <- NA_real_; x },
    "character matrix" = matrix("4", nrow = 4, ncol = 4),
    "non-matrix numeric" = as.numeric(1:16),
    "empty matrix" = matrix(numeric(0), nrow = 0, ncol = 0)
  )
  for (nm in names(bad)) {
    msg <- rb_catch_msg(.rust_fit_kl(V = bad[[nm]], k = 2L, iters = 12L, seed = 1L))
    expect_false(is.na(msg), info = paste0("hostile V [", nm, "] was not rejected"))
  }
  # Not one of the hostile calls may have reached the exec seam.
  expect_equal(length(mock$calls), 0L)
})

test_that("hostile scalar arguments are rejected before any exec", {
  mock <- rb_mock_exec(rb_responder())
  bin <- rb_use_dummy_bin()
  gate <- rb_mock_internal(".require_msuiter", function() invisible(TRUE))
  on.exit({
    rb_cleanup_mocks(mock, gate)
    rb_release_dummy_bin(bin)
  }, add = TRUE)

  v <- rb_live_v()[1:4, 1:4]
  hostile <- list(
    "k = 0" = list(k = 0L),
    "k = 1.5" = list(k = 1.5),
    "k = \"two\"" = list(k = "two"),
    "k = -2" = list(k = -2L),
    "iters = 0" = list(iters = 0L),
    "iters = 2.5" = list(iters = 2.5),
    "iters = \"12\"" = list(iters = "12"),
    "iters = Inf" = list(iters = Inf),
    "seed = -1" = list(seed = -1L),
    "seed = Inf" = list(seed = Inf),
    "seed = \"1\"" = list(seed = "1"),
    "replicate = 1.5" = list(seed = 1L, init = "seeded", replicate = 1.5),
    "init not a character" = list(init = 4L)
  )
  for (nm in names(hostile)) {
    a <- hostile[[nm]]
    a$V <- v
    if (is.null(a$k)) a$k <- 2L
    if (is.null(a$iters)) a$iters <- 12L
    if (is.null(a$seed) && is.null(a$init)) a$seed <- 1L
    msg <- rb_catch_msg(do.call(.rust_fit_kl, a))
    expect_false(is.na(msg), info = paste0("hostile scalar [", nm, "] was not rejected"))
  }
  expect_equal(length(mock$calls), 0L)
})

test_that("the init/seed/replicate coherence of the kernel is enforced R-side fail-closed", {
  mock <- rb_mock_exec(rb_responder())
  bin <- rb_use_dummy_bin()
  gate <- rb_mock_internal(".require_msuiter", function() invisible(TRUE))
  on.exit({
    rb_cleanup_mocks(mock, gate)
    rb_release_dummy_bin(bin)
  }, add = TRUE)

  v <- rb_live_v()[1:4, 1:4]
  w0 <- matrix(2, nrow = 4L, ncol = 2L)
  h0 <- matrix(2, nrow = 2L, ncol = 4L)

  rb_expect_fails(
    .rust_fit_kl(V = v, k = 2L, iters = 12L, seed = 1L, init = "external", W0 = w0, H0 = h0),
    "seed"
  )
  rb_expect_fails(
    .rust_fit_kl(V = v, k = 2L, iters = 12L, seed = NULL, init = "nndsvd"),
    "seed"
  )
  rb_expect_fails(
    .rust_fit_kl(V = v, k = 2L, iters = 12L, seed = 1L, init = "nndsvd", replicate = 3L),
    "replicate"
  )
  rb_expect_fails(
    .rust_fit_kl(V = v, k = 2L, iters = 12L, seed = 1L, init = "gpl-kernel"),
    "init"
  )
  rb_expect_fails(
    .rust_fit_kl(V = v, k = 2L, iters = 12L, seed = NULL, init = "external"),
    "W0"
  )
  expect_equal(length(mock$calls), 0L)
})

test_that(".stage_inputs validates dims against m/n/k fail-closed", {
  dir <- tempfile("rb-stage-")
  v <- rb_live_v()[1:4, 1:4]
  rb_expect_fails(
    .stage_inputs(dir, V = v[, 1:3], W0 = NULL, H0 = NULL, m = 4L, n = 4L, k = 2L),
    "V"
  )
  rb_expect_fails(
    .stage_inputs(dir, V = v, W0 = matrix(1, 4, 3), H0 = matrix(1, 2, 4), m = 4L, n = 4L, k = 2L),
    "W0"
  )
  rb_expect_fails(
    .stage_inputs(dir, V = v, W0 = matrix(1, 4, 2), H0 = matrix(1, 3, 4), m = 4L, n = 4L, k = 2L),
    "H0"
  )
  unlink(dir, recursive = TRUE)
})

test_that(".stage_inputs writes the frozen CSV format: %.17e, row-major, comma, no header", {
  dir <- tempfile("rb-stage-")
  v <- rb_live_v()
  w0 <- matrix(as.numeric(outer(seq_len(nrow(v)), seq_len(2), function(i, s) 2^((i + s) %% 4L))),
               nrow = nrow(v), ncol = 2, byrow = TRUE)
  h0 <- matrix(as.numeric(outer(seq_len(2), seq_len(ncol(v)), function(s, j) 2^((s + j + 1L) %% 4L))),
               nrow = 2, ncol = ncol(v), byrow = TRUE)
  .stage_inputs(dir, V = v, W0 = w0, H0 = h0,
                 m = nrow(v), n = ncol(v), k = 2L)
  expect_true(file.exists(file.path(dir, "V.csv")))
  expect_true(file.exists(file.path(dir, "W0.csv")))
  expect_true(file.exists(file.path(dir, "H0.csv")))
  lines <- readLines(file.path(dir, "V.csv"), warn = FALSE)
  expect_equal(length(lines), nrow(v))
  fields <- strsplit(lines, ",", fixed = TRUE)
  expect_true(all(sapply(fields, length) == ncol(v)))
  # every field is plain %.17e scientific notation, headerless
  expect_true(all(grepl("^-?[0-9]\\.[0-9]{17}e[-+][0-9]+$", unlist(fields))))
  # bit-exact round trip for the powers-of-two fixture: fields are
  # row-major, so rebuild row-major with no transpose
  back <- matrix(as.numeric(unlist(sapply(fields, as.numeric,
    simplify = FALSE), use.names = FALSE)),
    nrow = nrow(v), ncol = ncol(v), byrow = TRUE)
  expect_equal(dim(back), dim(v))
  expect_equal(as.vector(back), as.vector(v))
  unlink(dir, recursive = TRUE)
})

# ---- (b) exit-status handling ---------------------------------------------

test_that(".rust_status extracts the status attribute and defaults absent attrs to 0", {
  with_attr <- rb_kernel_result(2L)
  expect_identical(.rust_status(with_attr), 2L)
  no_attr <- 7L
  expect_identical(.rust_status(no_attr), 0L)
})

test_that("a non-zero fit exit is fail-closed, naming the code and the stderr report", {
  mock <- rb_mock_exec(rb_responder(status = 2L,
                                    stderr = c("error: V.csv: expected 16 values, got 15")))
  bin <- rb_use_dummy_bin()
  gate <- rb_mock_internal(".require_msuiter", function() invisible(TRUE))
  on.exit({
    rb_cleanup_mocks(mock, gate)
    rb_release_dummy_bin(bin)
  }, add = TRUE)
  msg <- rb_expect_fails(
    .rust_fit_kl(V = rb_live_v()[1:4, 1:4], k = 2L, iters = 12L, seed = 1L),
    "status 2"
  )
  expect_true(grepl("expected 16 values", msg, fixed = TRUE),
              info = "the kernel stderr report should be quoted into the error")
})

# ---- (c) handshake ---------------------------------------------------------

test_that("the handshake rejects protocol mismatch, garbage, silence and bad exits", {
  cases <- list(
    "protocol 2" = rb_kernel_result(0L, stdout = paste0(
      'msuiter-kernel protocol=2 engine="x" kl_eps=1e-12'
    )),
    "garbage line" = rb_kernel_result(0L, stdout = "msuiter-kernel oops"),
    "silent" = rb_kernel_result(0L, stdout = character())
  )
  mock <- NULL
  on.exit(rb_cleanup_mocks(mock), add = TRUE)
  for (nm in names(cases)) {
    mock <- rb_mock_exec(function(bin, args) cases[[nm]])
    msg <- rb_catch_msg(.kernel_handshake("/unused/bin"))
    rb_cleanup_mocks(mock)
    expect_false(is.na(msg), info = paste0("handshake case [", nm, "] was accepted"))
    # protocol mismatch must quote the offending line
    if (identical(nm, "protocol 2")) {
      expect_true(grepl("protocol=2", msg, fixed = TRUE),
                  info = "the fail-closed must quote the line it saw")
    }
  }
  # a non-zero handshake exit is named with its code
  mock2 <- rb_mock_exec(function(bin, args) rb_kernel_result(3L, stdout = rb_version_line))
  on.exit(rb_cleanup_mocks(mock2), add = TRUE)
  rb_expect_fails(.kernel_handshake("/unused/bin"), "status 3")
})

test_that("a good handshake line parses into protocol, engine and kl_eps", {
  mock <- rb_mock_exec(function(bin, args) rb_kernel_result(0L, stdout = rb_version_line))
  on.exit(rb_cleanup_mocks(mock), add = TRUE)
  hs <- .kernel_handshake("/unused/bin")
  expect_identical(hs$protocol, 1L)
  expect_equal(hs$kl_eps, 1e-12)
  expect_true(grepl("pure numeric kernels", hs$engine, fixed = TRUE))
  expect_equal(hs$raw, rb_version_line)
})

# ---- (d) meta cross-checks --------------------------------------------------

test_that("meta mismatches are fail-closed and name the offending key", {
  specs <- list(
    list(t = function(req, files) {
      files$meta.csv <- sub("^seed,.*$", "seed,99", files$meta.csv)
      files
    }, want = "seed"),
    list(t = function(req, files) {
      files$meta.csv <- sub("^iters,.*$", "iters,99", files$meta.csv)
      files
    }, want = "iters"),
    list(t = function(req, files) {
      files$meta.csv <- files$meta.csv[!grepl("^obj_final,", files$meta.csv)]
      files
    }, want = "obj_final")
  )
  mock <- NULL
  bin <- NULL
  gate <- NULL
  # Block-level on.exit is the household guard (see roundtrip suite): it
  # restores the tokens even when an expectation aborts the loop mid-run.
  # Installs run at loop-body top level with plain <- -- <<- across the
  # tryCatch frame does not reliably land in the binding the loop reads
  # back (doTryCatch promise-frame interaction), which leaked mocks here.
  on.exit({
    rb_cleanup_mocks(mock, gate)
    rb_release_dummy_bin(bin)
  }, add = TRUE)
  for (spec in specs) {
    mock <- rb_mock_exec(rb_responder(tamper = spec$t))
    bin <- rb_use_dummy_bin()
    gate <- rb_mock_internal(".require_msuiter", function() invisible(TRUE))
    msg <- rb_catch_msg(
      .rust_fit_kl(V = rb_live_v()[1:4, 1:4], k = 2L, iters = 12L, seed = 1L)
    )
    rb_cleanup_mocks(mock, gate)
    rb_release_dummy_bin(bin)
    expect_false(is.na(msg), info = paste0("meta tampering [", spec$want, "] was accepted"))
    expect_true(grepl(spec$want, msg, fixed = TRUE),
                info = paste0("error did not name the key [", spec$want, "]: ", msg))
    expect_true(grepl("meta.csv", msg, fixed = TRUE),
                info = paste0("error did not name the file: ", msg))
  }
})

test_that("the comma-bearing engine_crate value survives the meta parse verbatim", {
  mock <- rb_mock_exec(rb_responder())
  bin <- rb_use_dummy_bin()
  gate <- rb_mock_internal(".require_msuiter", function() invisible(TRUE))
  on.exit({
    rb_cleanup_mocks(mock, gate)
    rb_release_dummy_bin(bin)
  }, add = TRUE)
  fit <- .rust_fit_kl(V = rb_live_v()[1:4, 1:4], k = 2L, iters = 12L, seed = 1L)
  expect_equal(fit$meta$engine_crate, "engine: pure numeric kernels, FFI-free")
  expect_equal(fit$meta$protocol, "1")
  expect_equal(fit$meta$eps, "1e-12")
  expect_equal(length(fit$objective), 13L)
})

# ---- (e) output shape fail-closed -------------------------------------------

test_that("missing and short outputs are fail-closed, naming the file", {
  specs <- list(
    list(t = function(req, files) {
      files["W.csv"] <- NULL
      files
    }),
    list(t = function(req, files) {
      files$W.csv <- files$W.csv[-length(files$W.csv)]
      files
    }),
    list(t = function(req, files) {
      files$W.csv[2L] <- sub(",[^,]*$", "", files$W.csv[2L])
      files
    })
  )
  mock <- NULL
  bin <- NULL
  gate <- NULL
  # Block-level on.exit is the household guard (see roundtrip suite): it
  # restores the tokens even when an expectation aborts the loop mid-run.
  # Installs run at loop-body top level with plain <- -- <<- across the
  # tryCatch frame does not reliably land in the binding the loop reads
  # back (doTryCatch promise-frame interaction), which leaked mocks here.
  on.exit({
    rb_cleanup_mocks(mock, gate)
    rb_release_dummy_bin(bin)
  }, add = TRUE)
  for (spec in specs) {
    mock <- rb_mock_exec(rb_responder(tamper = spec$t))
    bin <- rb_use_dummy_bin()
    gate <- rb_mock_internal(".require_msuiter", function() invisible(TRUE))
    msg <- rb_catch_msg(
      .rust_fit_kl(V = rb_live_v()[1:4, 1:4], k = 2L, iters = 12L, seed = 1L)
    )
    rb_cleanup_mocks(mock, gate)
    rb_release_dummy_bin(bin)
    expect_false(is.na(msg), info = "output tampering was accepted")
    expect_true(grepl("W.csv", msg, fixed = TRUE),
                info = paste0("error did not name the file: ", msg))
  }
})

test_that("negative H entries and broken traces are fail-closed, naming the file", {
  specs <- list(
    list(t = function(req, files) {
      files$H.csv[1L] <- sub("^[^,]*", "-2.50000000000000000e+00", files$H.csv[1L])
      files
    }, want = "H.csv"),
    list(t = function(req, files) {
      last <- length(files$trace.csv)
      files$trace.csv[last] <- sub(",.*$", ",NaN", files$trace.csv[last])
      files
    }, want = "trace.csv"),
    list(t = function(req, files) {
      files$trace.csv <- files$trace.csv[-length(files$trace.csv)]
      files
    }, want = "trace.csv"),
    list(t = function(req, files) {
      files$meta.csv <- sub("^trace_len,.*$", "trace_len,999", files$meta.csv)
      files
    }, want = "meta.csv")
  )
  mock <- NULL
  bin <- NULL
  gate <- NULL
  # Block-level on.exit is the household guard (see roundtrip suite): it
  # restores the tokens even when an expectation aborts the loop mid-run.
  # Installs run at loop-body top level with plain <- -- <<- across the
  # tryCatch frame does not reliably land in the binding the loop reads
  # back (doTryCatch promise-frame interaction), which leaked mocks here.
  on.exit({
    rb_cleanup_mocks(mock, gate)
    rb_release_dummy_bin(bin)
  }, add = TRUE)
  for (spec in specs) {
    mock <- rb_mock_exec(rb_responder(tamper = spec$t))
    bin <- rb_use_dummy_bin()
    gate <- rb_mock_internal(".require_msuiter", function() invisible(TRUE))
    msg <- rb_catch_msg(
      .rust_fit_kl(V = rb_live_v()[1:4, 1:4], k = 2L, iters = 12L, seed = 1L)
    )
    rb_cleanup_mocks(mock, gate)
    rb_release_dummy_bin(bin)
    expect_false(is.na(msg), info = paste0("tamper was accepted (want ", spec$want, ")"))
    expect_true(grepl(spec$want, msg, fixed = TRUE),
                info = paste0("error did not name [", spec$want, "]: ", msg))
  }
})

# ---- (f) backend idiom -------------------------------------------------------

test_that("match.arg backend idiom: exact, partial and case-sensitive rejection", {
  expect_identical(.r_assert_backend("R"), "R")
  expect_identical(.r_assert_backend("rust"), "rust")
  expect_identical(.r_assert_backend(c("R", "rust")), "R")   # untouched default
  expect_identical(.r_assert_backend("rus"), "rust")          # partial match, documented
  rb_expect_fails(.r_assert_backend("Rust"), ".")
  rb_expect_fails(.r_assert_backend("nimfa"), ".")
  rb_expect_fails(.r_assert_backend(4L), ".")
})

test_that("the exported surface validates backend and rejects pynmf x rust combinations early", {
  m <- rb_nmf_matrix()
  rb_expect_fails(sig_extract(m, 2, backend = "Rust"), ".")
  rb_expect_fails(sig_extract(m, 2, backend = "nimfa"), ".")
  rb_expect_fails(
    sig_extract(m, 2, backend = "rust", pynmf = TRUE),
    "mutually exclusive"
  )
  rb_expect_fails(
    bp_extract_signatures(m, range = 2, backend = "rust", pynmf = TRUE,
                          n_bootstrap = 20L, n_nmf_run = 1L),
    "mutually exclusive"
  )
  rb_expect_fails(sig_estimate(m, range = 2, backend = "qux"), ".")
})

test_that("the \"rust\" backend is gated at runtime on the optional msuiter package", {
  if (requireNamespace("msuiter", quietly = TRUE)) {
    expect_no_error(.require_msuiter())
  } else {
    msg <- rb_catch_msg(.require_msuiter())
    expect_false(is.na(msg))
    expect_true(grepl("msuiter", msg, fixed = TRUE))
    expect_true(grepl("rust", msg, fixed = TRUE),
                 info = "the gate message must name the backend it gates")
  }
})

# ---- (g) nrun seeding law -----------------------------------------------------

test_that("nrun sweeps advance --seed by one per run and best-of picks the minimum obj_final", {
  responder <- rb_responder(obj_finals = c(64, 8, 32))
  mock <- rb_mock_exec(responder)
  bin <- rb_use_dummy_bin()
  gate <- rb_mock_internal(".require_msuiter", function() invisible(TRUE))
  on.exit({
    rb_cleanup_mocks(mock, gate)
    rb_release_dummy_bin(bin)
  }, add = TRUE)
  opts <- getOption("sigminer.rust_iters")
  options(sigminer.rust_iters = 12L)
  on.exit(options(sigminer.rust_iters = opts), add = TRUE)

  v <- rb_live_v()[1:4, 1:4]
  fits <- .rust_runs(V = v, k = 2L, iters = 12L, seed = 9L, nrun = 3L,
                     init = "nndsvd")
  expect_equal(length(fits), 3L)
  expect_identical(rb_fit_values(mock, "seed"), c("9", "10", "11"))

  # reset the scripted objective sequence for the best-of phase
  mock$calls <- list()
  mock$run <- 0L
  best <- .rust_best_of_nrun(V = v, k = 2L, iters = 12L, seed = 9L,
                              init = "nndsvd", nrun = 3L)
  expect_equal(as.numeric(best$meta$obj_final), 8)
  expect_equal(tail(best$objective, 1L), 8)
  # the option value travelled verbatim to every spawn
  expect_true(all(rb_fit_values(mock, "iters") == "12"))
})

test_that("external init is refused for nrun > 1 and runs seedless for nrun = 1", {
  mock <- rb_mock_exec(rb_responder())
  bin <- rb_use_dummy_bin()
  gate <- rb_mock_internal(".require_msuiter", function() invisible(TRUE))
  on.exit({
    rb_cleanup_mocks(mock, gate)
    rb_release_dummy_bin(bin)
  }, add = TRUE)

  v <- rb_live_v()[1:4, 1:4]
  w0 <- matrix(2, nrow = 4L, ncol = 2L)
  h0 <- matrix(2, nrow = 2L, ncol = 4L)

  rb_expect_fails(
    .rust_best_of_nrun(V = v, k = 2L, iters = 12L, seed = NULL,
                       init = "external", nrun = 3L, W0 = w0, H0 = h0),
    "identical runs; use nrun = 1", fixed = TRUE
  )
  mock$calls <- list()

  fit <- .rust_best_of_nrun(V = v, k = 2L, iters = 12L, seed = NULL,
                            init = "external", nrun = 1L, W0 = w0, H0 = h0)
  argv <- mock$calls[[2L]]$args   # call 1 is the handshake
  expect_false("--seed" %in% argv, info = "external must not carry --seed")
  expect_equal(rb_parse_argv(argv)$init, "external")
  expect_equal(fit$meta$seed, "na")
})

# ---- (h) resolution precedence -------------------------------------------------

test_that(".kernel_resolve honours option > env > PATH, and fail-closes naming every source", {
  old_opt <- getOption("sigminer.rust_kernel")
  old_env <- Sys.getenv("MSUITER_KERNEL_BIN")
  old_path <- Sys.getenv("PATH")
  on.exit({
    options(sigminer.rust_kernel = old_opt)
    if (nzchar(old_env)) Sys.setenv(MSUITER_KERNEL_BIN = old_env) else Sys.unsetenv("MSUITER_KERNEL_BIN")
    Sys.setenv(PATH = old_path)
  }, add = TRUE)

  opt_bin <- tempfile("rb-opt-bin-")
  env_bin <- tempfile("rb-env-bin-")
  which_dir <- tempfile("rb-which-dir-")
  dir.create(which_dir, recursive = TRUE, showWarnings = FALSE)
  file.create(opt_bin)
  file.create(env_bin)
  file.create(file.path(which_dir, "msuiter-kernel"))
  Sys.chmod(file.path(which_dir, "msuiter-kernel"), mode = "0755")

  options(sigminer.rust_kernel = opt_bin)
  Sys.setenv(MSUITER_KERNEL_BIN = env_bin)
  Sys.setenv(PATH = paste0(which_dir, ":", old_path))
  expect_equal(.kernel_resolve(), normalizePath(opt_bin))

  options(sigminer.rust_kernel = NULL)
  expect_equal(.kernel_resolve(), normalizePath(env_bin))

  Sys.unsetenv("MSUITER_KERNEL_BIN")
  expect_equal(.kernel_resolve(), normalizePath(file.path(which_dir, "msuiter-kernel")))

  # all three sources missing: the message names each of them
  Sys.setenv(PATH = old_path)
  unlink(opt_bin)
  unlink(env_bin)
  unlink(which_dir, recursive = TRUE)
  ghost <- tempfile("rb-ghost-")
  options(sigminer.rust_kernel = ghost)
  msg <- rb_expect_fails(.kernel_resolve(), ".")
  expect_true(grepl("sigminer.rust_kernel", msg, fixed = TRUE))
  expect_true(grepl("MSUITER_KERNEL_BIN", msg, fixed = TRUE))
  expect_true(grepl("msuiter-kernel", msg, fixed = TRUE))
  # the ghost option path is named in the tried-sources report
  expect_true(grepl(ghost, msg, fixed = TRUE))
})

# ---- (b, cont.) floor-safety self-scan ----------------------------------------

test_that("the rust boundary uses the floor-safe status helper and never the base null-default operator", {
  src <- rb_repo_r_path("rust_backend.R")
  expect_true(grepl("rust_backend\\.R$", src))
  lines <- readLines(src, warn = FALSE)
  # Built from pieces so the forbidden token never appears literally here
  # either: the whole P1 change set stays at zero hits for a plain grep.
  forbidden <- paste0("%", "|", "|", "%")
  expect_equal(length(grep(forbidden, lines, fixed = TRUE)), 0L)
  expect_true(any(grepl("attr(res, \"status\")", lines, fixed = TRUE)),
              info = ".rust_status must read the status attribute directly")
})

# ---- sig_extract rust branch (attr stamps + wiring) ---------------------------

test_that("sig_extract backend = \"R\" explicit is identical to the untouched default", {
  m <- rb_nmf_matrix()
  default <- suppressMessages(suppressWarnings(
    sig_extract(m, 2, nrun = 2, seed = 9L, method = "brunet", cores = 1)
  ))
  explicit <- suppressMessages(suppressWarnings(
    sig_extract(m, 2, nrun = 2, seed = 9L, method = "brunet", cores = 1, backend = "R")
  ))
  expect_identical(default$Raw$W, explicit$Raw$W)
  expect_identical(default$Raw$H, explicit$Raw$H)
  expect_identical(rb_wired_attrs(default), rb_wired_attrs(explicit))
  expect_identical(attr(explicit, "call_method"), "NMF")
})

test_that("sig_extract backend = \"rust\" returns a Signature with the kernel attr stamps", {
  mock <- rb_mock_exec(rb_responder())
  bin <- rb_use_dummy_bin()
  gate <- rb_mock_internal(".require_msuiter", function() invisible(TRUE))
  on.exit({
    rb_cleanup_mocks(mock, gate)
    rb_release_dummy_bin(bin)
  }, add = TRUE)
  opts <- getOption("sigminer.rust_iters")
  options(sigminer.rust_iters = 12L)
  on.exit(options(sigminer.rust_iters = opts), add = TRUE)

  m <- rb_nmf_matrix()
  res <- suppressMessages(suppressWarnings(
    sig_extract(m, 2, nrun = 2, seed = 9L, backend = "rust")
  ))
  expect_s3_class(res, "Signature")
  expect_identical(attr(res, "call_method"), "msuiter-kernel fit-kl protocol=1")
  expect_equal(as.character(attr(res, "seed")), "9")
  expect_equal(as.character(attr(res, "nrun")), "2")
  expect_equal(dim(res$Raw$W), c(nrow(m), 2L))
  expect_equal(dim(res$Raw$H), c(2L, ncol(m)))
  expect_null(res$Raw$nmf_obj)
  # the nrun runs advanced the seed law
  expect_identical(rb_fit_values(mock, "seed"), c("9", "10"))
})

test_that("sig_extract rust branch refuses non-KL methods", {
  mock <- rb_mock_exec(rb_responder())
  bin <- rb_use_dummy_bin()
  gate <- rb_mock_internal(".require_msuiter", function() invisible(TRUE))
  on.exit({
    rb_cleanup_mocks(mock, gate)
    rb_release_dummy_bin(bin)
  }, add = TRUE)
  rb_expect_fails(
    sig_extract(rb_nmf_matrix(), 2, backend = "rust", method = "lee"),
    "KL"
  )
})

# ---- sig_estimate rust branch ---------------------------------------------------

test_that("sig_estimate backend = \"rust\" builds a Survey via the seed law and fail-closes on NMF-only arms", {
  mock <- rb_mock_exec(rb_responder())
  bin <- rb_use_dummy_bin()
  gate <- rb_mock_internal(".require_msuiter", function() invisible(TRUE))
  on.exit({
    rb_cleanup_mocks(mock, gate)
    rb_release_dummy_bin(bin)
  }, add = TRUE)
  opts <- getOption("sigminer.rust_iters")
  options(sigminer.rust_iters = 12L)
  on.exit(options(sigminer.rust_iters = opts), add = TRUE)

  m <- rb_nmf_matrix()
  res <- suppressMessages(suppressWarnings(
    sig_estimate(m, range = c(2L), nrun = 2L, seed = 5L, backend = "rust")
  ))
  expect_s3_class(res, "Survey")
  expect_false(is.null(res$survey))
  expect_true(all(c("rank", "rss", "evar", "cophenetic") %in% colnames(res$survey)))
  expect_equal(as.numeric(res$survey$rank), 2)
  expect_true(all(is.finite(res$survey$rss)))
  expect_true(res$survey$cophenetic >= 0 && res$survey$cophenetic <= 1)
  expect_identical(rb_fit_values(mock, "seed"), c("5", "6"))

  rb_expect_fails(
    sig_estimate(m, range = 2L, backend = "rust", use_random = TRUE), "use_random"
  )
  rb_expect_fails(
    sig_estimate(m, range = 2L, backend = "rust", save_plots = TRUE), "save_plots"
  )
  keep <- suppressMessages(suppressWarnings(
    sig_estimate(m, range = 2L, nrun = 1L, backend = "rust", keep_nmfObj = TRUE)
  ))
  expect_false(is.null(keep$nmfEstimate))
})

# ---- iters option ---------------------------------------------------------------

test_that("the sigminer.rust_iters option is honoured verbatim, defaults to 2000, rejects junk", {
  old_opt <- getOption("sigminer.rust_iters")
  on.exit(options(sigminer.rust_iters = old_opt), add = TRUE)
  options(sigminer.rust_iters = 7L)
  expect_identical(.rust_iters(), 7L)
  options(sigminer.rust_iters = NULL)
  expect_identical(.rust_iters(), 2000L)
  options(sigminer.rust_iters = -3)
  rb_expect_fails(.rust_iters(), "sigminer.rust_iters")
  options(sigminer.rust_iters = "banana")
  rb_expect_fails(.rust_iters(), "sigminer.rust_iters")
})

# ---- (j)-(n) real-kernel battery (env-gated) ------------------------------------

test_that("(j) the live kernel prints the frozen protocol-1 handshake", {
  if (!rb_live_kernel_available()) {
    skip("MSUITER_TEST_KERNEL_BIN not set: the real-kernel battery is opt-in")
  }
  bin <- Sys.getenv("MSUITER_TEST_KERNEL_BIN")
  res <- .kernel_exec(bin, c("--version"))
  expect_identical(.rust_status(res), 0L)
  lines <- attr(res, "stdout")
  expect_true(any(grepl("msuiter-kernel", lines)), info = "no handshake line on stdout")
  line <- lines[grepl("msuiter-kernel", lines)][1L]
  expect_equal(line, rb_version_line)
  hs <- .kernel_handshake(bin)
  expect_identical(hs$protocol, 1L)
  expect_equal(hs$kl_eps, 1e-12)
})

test_that("(k) a synthetic end-to-end fit returns dims-clean, finite, meta-verified results and is rerunnable", {
  if (!rb_live_kernel_available()) {
    skip("MSUITER_TEST_KERNEL_BIN not set: the real-kernel battery is opt-in")
  }
  old_opt <- getOption("sigminer.rust_kernel")
  options(sigminer.rust_kernel = Sys.getenv("MSUITER_TEST_KERNEL_BIN"))
  on.exit(options(sigminer.rust_kernel = old_opt), add = TRUE)
  gate <- rb_mock_internal(".require_msuiter", function() invisible(TRUE))
  on.exit(rb_cleanup_mocks(gate), add = TRUE)

  v <- rb_live_v()
  fit1 <- .rust_fit_kl(V = v, k = 2L, iters = 12L, seed = 1L, init = "nndsvd")
  expect_equal(dim(fit1$W), c(nrow(v), 2L))
  expect_equal(dim(fit1$H), c(2L, ncol(v)))
  expect_true(all(is.finite(fit1$W)) && all(fit1$W >= 0))
  expect_true(all(is.finite(fit1$H)) && all(fit1$H >= 0))
  expect_equal(length(fit1$objective), 13L)
  expect_true(all(is.finite(fit1$objective)))
  expect_equal(fit1$meta$protocol, "1")
  expect_equal(fit1$meta$init, "nndsvd")
  expect_equal(fit1$meta$seed, "1")
  expect_equal(fit1$meta$m, as.character(nrow(v)))
  expect_equal(fit1$meta$n, as.character(ncol(v)))
  expect_equal(fit1$meta$iters, "12")

  fit2 <- .rust_fit_kl(V = v, k = 2L, iters = 12L, seed = 1L, init = "nndsvd")
  expect_identical(fit1$W, fit2$W)
  expect_identical(fit1$H, fit2$H)
  expect_identical(fit1$objective, fit2$objective)
  expect_identical(fit1$meta, fit2$meta)
})

test_that("(l) bit-equivalence: the same external init drives both engines to the same W", {
  if (!rb_live_kernel_available()) {
    skip("MSUITER_TEST_KERNEL_BIN not set: the real-kernel battery is opt-in")
  }
  nmf_ver <- as.character(packageVersion("NMF"))
  if (nmf_ver < "1.5.0") {
    skip(paste0("NMF ", nmf_ver, " < 1.5.0: the modern nmf(mode/seed/normalize/init) surface is absent on this host; this gate runs where the engine ships >= 1.5.0 (logged in the lane report)"))
  }
  old_opt <- getOption("sigminer.rust_kernel")
  options(sigminer.rust_kernel = Sys.getenv("MSUITER_TEST_KERNEL_BIN"))
  on.exit(options(sigminer.rust_kernel = old_opt), add = TRUE)
  gate <- rb_mock_internal(".require_msuiter", function() invisible(TRUE))
  on.exit(rb_cleanup_mocks(gate), add = TRUE)

  v <- rb_live_v()
  k <- 2L
  i <- 12L
  s <- 3L
  w0 <- matrix(as.numeric(outer(seq_len(nrow(v)), seq_len(k), function(a, b) 2^((a + b) %% 3L))),
               nrow = nrow(v), ncol = k, byrow = TRUE)
  h0 <- matrix(as.numeric(outer(seq_len(k), seq_len(ncol(v)), function(a, b) 2^((a + b + 1L) %% 3L))),
               nrow = k, ncol = ncol(v), byrow = TRUE)

  fit <- .rust_fit_kl(V = v, k = k, iters = i, seed = NULL, init = "external",
                      W0 = w0, H0 = h0)
  suppressMessages(library(NMF))
  nmf_out <- suppressMessages(suppressWarnings(
    NMF::nmf(v, rank = k, nrun = 1, mode = "kl", seed = s, normalize = FALSE,
             init = list(W = w0, H = h0))
  ))
  wn <- unname(as.matrix(NMF::basis(nmf_out)))
  dimnames(wn) <- NULL

  expect_equal(dim(wn), c(nrow(v), k))
  relerr <- norm(wn - fit$W, "F") / norm(fit$W, "F")
  message(paste0(
    "P1(l) gate log | NMF ", nmf_ver,
    " | V nrow=", nrow(v), " ncol=", ncol(v), " k=", k, " iters=", i,
    " | kernel iterations=", fit$meta$iterations,
    " obj_final=", fit$meta$obj_final,
    " | relerr(W)=", format(relerr, scientific = TRUE)
  ))
  expect_lt(relerr, 1e-10)
})

test_that("(m) the rust bp path is thread-invariant: workers = 1 matches workers = 2", {
  if (!rb_live_kernel_available()) {
    skip("MSUITER_TEST_KERNEL_BIN not set: the real-kernel battery is opt-in")
  }
  if (!requireNamespace("synchronicity", quietly = TRUE) ||
      !requireNamespace("lpSolve", quietly = TRUE)) {
    skip("synchronicity / lpSolve not installed: bp pipeline is unavailable")
  }
  old_opt <- getOption("sigminer.rust_kernel")
  old_iters <- getOption("sigminer.rust_iters")
  options(sigminer.rust_kernel = Sys.getenv("MSUITER_TEST_KERNEL_BIN"),
          sigminer.rust_iters = 12L)
  on.exit({
    options(sigminer.rust_kernel = old_opt, sigminer.rust_iters = old_iters)
  }, add = TRUE)
  # The fan-out workers must resolve the lane namespace, not the installed
  # CRAN build (see rb_lane_worker_lib).
  lane_lib <- rb_lane_worker_lib()
  if (is.null(lane_lib)) {
    skip("no worker-side lane package install for the future fan-out")
  }
  old_libpaths <- .libPaths(c(lane_lib, .libPaths()))
  on.exit(.libPaths(old_libpaths), add = TRUE)
  gate <- rb_mock_internal(".require_msuiter", function() invisible(TRUE))
  on.exit(rb_cleanup_mocks(gate), add = TRUE)

  m <- rb_nmf_matrix()
  run_bp <- function(workers) {
    suppressMessages(suppressWarnings(
      bp_extract_signatures(
        m,
        range = 2L,
        n_bootstrap = 20L,          # >= 10 * n_nmf_run drives the future fan-out
        n_nmf_run = 2L,
        cores = workers,
        seed = 7L,
        handle_hyper_mutation = FALSE,
        only_core_stats = TRUE,
        keep_cache = FALSE,
        cache_dir = file.path(tempdir(), paste0("rb-bp-", workers)),
        backend = "rust"
      )
    ))
  }
  r1 <- run_bp(1)
  r2 <- run_bp(2)
  # Worker count must not touch the semantic payload of the ExtractionResult
  # (per-rank signature/exposure lists and the stats table).
  expect_identical(names(r1), names(r2))
  expect_identical(r1$signature, r2$signature)
  expect_identical(r1$exposure, r2$exposure)
  expect_identical(r1$stats, r2$stats)
})

test_that("(n) zero-lock regression: NNDSVD over a zero column keeps the objective strictly moving", {
  if (!rb_live_kernel_available()) {
    skip("MSUITER_TEST_KERNEL_BIN not set: the real-kernel battery is opt-in")
  }
  old_opt <- getOption("sigminer.rust_kernel")
  options(sigminer.rust_kernel = Sys.getenv("MSUITER_TEST_KERNEL_BIN"))
  on.exit(options(sigminer.rust_kernel = old_opt), add = TRUE)
  gate <- rb_mock_internal(".require_msuiter", function() invisible(TRUE))
  on.exit(rb_cleanup_mocks(gate), add = TRUE)

  v <- rb_live_v(zero_column = TRUE)
  expect_equal(as.vector(v[, 6L]), rep(0, nrow(v)))
  fit <- .rust_fit_kl(V = v, k = 2L, iters = 12L, seed = 1L, init = "nndsvd")
  first_steps <- diff(fit$objective[1:3])   # the first k = 2 update steps
  expect_true(all(first_steps < 0),
               info = paste0("objective froze in the first k steps: ",
                             paste(format(fit$objective[1:3]), collapse = " -> ")))
})
