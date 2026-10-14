# Internal interface to the certified msuiter KL-MU kernel (parity plan P1,
# closes #478). Everything in this file is internal: the user-facing surface
# is the `backend = c("R", "rust")` argument on sig_extract(), sig_estimate()
# and bp_extract_signatures().
#
# Design rulings that must not silently drift:
#  * The engine is reached ONLY as an isolated process through the frozen
#    protocol-1 CLI contract, through the single contact point
#    .kernel_exec() (system2 with a character-vector argv, never a shell).
#  * The boundary is Suggests-gated on 'msuiter', reticulate-idiom.
#  * Exit-status extraction is floor-safe: this package targets R >= 3.5 and
#    the base null-default operator is a 4.4.0 export, so the status
#    attribute is read directly in .rust_status() below.
#  * Seed law: every non-external invocation passes --seed explicitly and
#    run i of nrun uses seed + i - 1L (the pynmf path's silent seed drop is
#    the documented gap this must not repeat).
#  * All numerics crossing the boundary are validated fail-closed BEFORE
#    any spawn: no NaN/Inf/negatives/non-integer counts ever reach the
#    kernel; no un-validated output is ever returned.

.RUST_PROTOCOL <- 1L
.RUST_BACKENDS <- c("R", "rust")
.RUST_KERNEL_INITS <- c("external", "seeded", "nndsvd", "nndsvda")
.RUST_DEFAULT_ITERS <- 2000L
.RUST_CALL_METHOD <- "msuiter-kernel fit-kl protocol=1"
# Frozen meta.csv key order of the protocol-1 contract (bin source
# msuiter-kernel.rs lines 329-368, verified live).
.RUST_META_KEYS <- c("protocol", "command", "init", "seed", "replicate",
                     "m", "n", "k", "iters", "iterations", "trace_len",
                     "obj_init", "obj_final", "eps", "engine_crate")

# ---------------------------------------------------------------------------
# Status and process contact
# ---------------------------------------------------------------------------

# Floor-safe exit-status extraction. An absent "status" attribute (e.g. a
# mocked result) reads as success 0.
.rust_status <- function(res) {
  st <- attr(res, "status")
  if (is.null(st)) 0L else as.integer(st)
}

# The single exec seam. Character-vector argv, no shell; stdout/stderr are
# captured to files and attached as attributes so every caller (and the
# fail-closed reporting) sees them. Tests replace this function through
# assignInNamespace (house mock idiom).
.kernel_exec <- function(bin, args) {
  args <- as.character(args)
  out_f <- tempfile("rust-out-")
  err_f <- tempfile("rust-err-")
  on.exit({
    unlink(out_f)
    unlink(err_f)
  }, add = TRUE)
  res <- system2(bin, args, stdout = out_f, stderr = err_f)
  out <- character()
  err <- character()
  if (file.exists(out_f)) {
    out <- tryCatch(readLines(out_f, warn = FALSE), error = function(e) character())
  }
  if (file.exists(err_f)) {
    err <- tryCatch(readLines(err_f, warn = FALSE), error = function(e) character())
  }
  attr(res, "stdout") <- out
  attr(res, "stderr") <- err
  res
}

# ---------------------------------------------------------------------------
# Resolution and handshake
# ---------------------------------------------------------------------------

.rust_shown <- function(x) {
  if (is.null(x) || length(x) < 1L || is.na(x[[1L]]) || !nzchar(x[[1L]])) {
    "<unset>"
  } else {
    as.character(x[[1L]])
  }
}

# Resolution order: option sigminer.rust_kernel > environment variable
# MSUITER_KERNEL_BIN > PATH lookup for msuiter-kernel. Each candidate is
# probed with file.exists (never by branching on OS identity strings).
# Fail-closed and names every source it tried.
.kernel_resolve <- function() {
  opt <- getOption("sigminer.rust_kernel")
  env <- Sys.getenv("MSUITER_KERNEL_BIN")
  wh <- Sys.which("msuiter-kernel")
  wh_path <- if (length(wh) > 0L) as.character(wh[[1L]]) else NULL
  cands <- list(opt, env, wh_path)
  for (cand in cands) {
    if (!is.null(cand) && nzchar(cand) && file.exists(cand)) {
      return(normalizePath(cand))
    }
  }
  stop(paste0(
    "msuiter kernel binary not found. Tried: option sigminer.rust_kernel = ",
    .rust_shown(opt), "; environment variable MSUITER_KERNEL_BIN = ",
    .rust_shown(env), "; PATH lookup of 'msuiter-kernel' = ",
    .rust_shown(wh_path), ". Point one of these at the certified ",
    "msuiter-kernel binary (built from msuiter trunk, protocol=1)."
  ), domain = NA)
}

# Runs '<bin> --version' and validates the frozen handshake line:
#   msuiter-kernel protocol=1 engine="<engine>" kl_eps=<value>
# A protocol other than 1 is fail-closed and quotes the line it saw.
.kernel_handshake <- function(bin) {
  res <- .kernel_exec(bin, c("--version"))
  status <- .rust_status(res)
  if (!identical(status, 0L)) {
    stop(paste0("msuiter-kernel --version exited with status ", status,
                "; stderr: ", paste(as.character(attr(res, "stderr")),
                                    collapse = " | ")), domain = NA)
  }
  lines <- as.character(attr(res, "stdout"))
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0L) {
    stop("msuiter-kernel --version produced no output; cannot verify the ",
         "protocol-1 handshake.", domain = NA)
  }
  cand <- lines[grepl("^msuiter-kernel ", lines)]
  if (length(cand) == 0L) {
    stop(paste0("msuiter-kernel --version did not print its protocol line; saw: ",
                paste(lines, collapse = " / ")), domain = NA)
  }
  line <- cand[1L]
  shape <- paste0("^msuiter-kernel protocol=[0-9]+ engine=\"[^\"]*\" ",
                  "kl_eps=[-+0-9.eE]+$")
  if (!grepl(shape, line)) {
    stop(paste0("unexpected msuiter-kernel handshake line: ", line), domain = NA)
  }
  toks <- strsplit(line, " ", fixed = TRUE)[[1L]]
  proto_tok <- toks[grepl("^protocol=", toks)][1L]
  eps_tok <- toks[grepl("^kl_eps=", toks)][1L]
  if (is.na(proto_tok) || is.na(eps_tok)) {
    stop(paste0("unexpected msuiter-kernel handshake line: ", line), domain = NA)
  }
  protocol <- suppressWarnings(as.integer(sub("^protocol=", "", proto_tok)))
  if (is.na(protocol) || !identical(protocol, .RUST_PROTOCOL)) {
    stop(paste0("unsupported msuiter-kernel protocol: expected protocol=",
                .RUST_PROTOCOL, ", line said: ", line), domain = NA)
  }
  kl_eps <- suppressWarnings(as.numeric(sub("^kl_eps=", "", eps_tok)))
  if (!is.finite(kl_eps) || kl_eps < 0) {
    stop(paste0("unreadable kl_eps in handshake line: ", line), domain = NA)
  }
  engine <- sub("^.*engine=\"", "", line, )
  engine <- sub("\"[^\"]*$", "", engine)
  list(protocol = protocol, engine = engine, kl_eps = kl_eps, raw = line)
}

# Suggests gate in the reticulate idiom (see best_practice.R); the message
# names both the package and the backend it guards.
.require_msuiter <- function() {
  if (!requireNamespace("msuiter", quietly = TRUE)) {
    stop("Package 'msuiter' is required for the rust backend, please ",
         "install it firstly!", domain = NA)
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Argument validation (all fail-closed, all BEFORE any spawn)
# ---------------------------------------------------------------------------

# House match.arg idiom (sig_fit.R): validates against the exact choices and
# allows unambiguous partial matches ("rus" -> "rust").
.r_assert_backend <- function(backend) {
  if (!is.character(backend) || length(backend) < 1L || any(is.na(backend))) {
    stop("'backend' must be a character vector of \"R\" / \"rust\" choices ",
         "(the default is c(\"R\", \"rust\")).", domain = NA)
  }
  match.arg(backend, choices = .RUST_BACKENDS)
}

# Package option for the fixed iteration budget of the kernel (the engine
# has no convergence stop; NMF's brunet-KL runs until convergence, the
# kernel runs a fixed count -- this option maps one onto the other).
.rust_iters <- function() {
  v <- getOption("sigminer.rust_iters")
  if (is.null(v)) return(.RUST_DEFAULT_ITERS)
  if (!is.numeric(v) || length(v) != 1L || !is.finite(v) ||
      v != floor(v) || v < 1) {
    stop("option sigminer.rust_iters must be a single positive integer, ",
         "got: ", paste(as.character(v), collapse = ","), domain = NA)
  }
  as.integer(v)
}

# Single positive (or non-negative) whole count.
.rust_count <- function(x, name, allow.zero = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x != floor(x) || x < (if (allow.zero) 0 else 1)) {
    stop("'", name, "' must be a single ",
         if (allow.zero) "non-negative " else "positive ",
         "whole number, got: ", paste(as.character(x), collapse = ","),
         domain = NA)
  }
  as.integer(x)
}

# Numeric, finite, non-negative matrix.
.rust_validate_matrix <- function(x, name) {
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) < 1L || ncol(x) < 1L) {
    stop("'", name, "' must be a numeric matrix with at least one row and ",
         "one column.", domain = NA)
  }
  if (!all(is.finite(x))) {
    stop("'", name, "' must contain only finite values (no NaN/Inf/NA).",
         domain = NA)
  }
  if (any(x < 0)) {
    stop("'", name, "' must be non-negative (it is the NMF input).",
         domain = NA)
  }
  x
}

# ---------------------------------------------------------------------------
# Staging and parsing of the CSV contract
# ---------------------------------------------------------------------------

# Split a line at its first comma only (the meta value 'engine_crate'
# legitimately contains commas). Returns NULL when there is no comma.
.rust_split1 <- function(line) {
  chars <- strsplit(line, "")[[1L]]
  hit <- which(chars == ",")
  if (length(hit) == 0L) return(NULL)
  i <- hit[1L]
  list(left = substr(line, 1L, i - 1L),
       right = substr(line, i + 1L, length(chars)))
}

.rust_csv_rows <- function(mat) {
  apply(mat, 1L, function(r) paste(sprintf("%.17e", r), collapse = ","))
}

# Write the kernel's input grid, fail-closed on dimensions first.
.stage_inputs <- function(dir, V, W0, H0, m, n, k) {
  V <- .rust_validate_matrix(V, "V")
  if (!identical(nrow(V), as.integer(m)) || !identical(ncol(V), as.integer(n))) {
    stop("'V' must be ", m, " x ", n, ", got ", nrow(V), " x ", ncol(V), ".",
         domain = NA)
  }
  if (!is.null(W0)) {
    W0 <- .rust_validate_matrix(W0, "W0")
    if (!identical(dim(W0), c(as.integer(m), as.integer(k)))) {
      stop("'W0' must be ", m, " x ", k, ", got ", nrow(W0), " x ",
           ncol(W0), ".", domain = NA)
    }
  }
  if (!is.null(H0)) {
    H0 <- .rust_validate_matrix(H0, "H0")
    if (!identical(dim(H0), c(as.integer(k), as.integer(n)))) {
      stop("'H0' must be ", k, " x ", n, ", got ", nrow(H0), " x ",
           ncol(H0), ".", domain = NA)
    }
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(.rust_csv_rows(V), con = file.path(dir, "V.csv"))
  if (!is.null(W0)) {
    writeLines(.rust_csv_rows(W0), con = file.path(dir, "W0.csv"))
  }
  if (!is.null(H0)) {
    writeLines(.rust_csv_rows(H0), con = file.path(dir, "H0.csv"))
  }
  invisible(dir)
}

# Read one m x k %.17e comma grid back; every deviation is fail-closed and
# names the file.
.rust_read_matrix_csv <- function(path, nrows, ncols, name) {
  if (!file.exists(path)) {
    stop("msuiter-kernel did not produce ", name, " in the output directory.",
         domain = NA)
  }
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  if (!identical(length(lines), as.integer(nrows))) {
    stop(name, " must hold exactly ", nrows, " rows, got ", length(lines), ".",
         domain = NA)
  }
  fields <- strsplit(lines, ",", fixed = TRUE)
  bad <- which(sapply(fields, length) != as.integer(ncols))
  if (length(bad) > 0L) {
    stop(name, " row ", bad[1L], " has ", sapply(fields, length)[[bad[1L]]],
         " comma-separated fields, expected ", ncols, ".", domain = NA)
  }
  vals <- as.numeric(unlist(fields, use.names = FALSE))
  if (!all(is.finite(vals))) {
    stop(name, " contains a non-finite value.", domain = NA)
  }
  if (any(vals < 0)) {
    hit <- which(vals < 0)[1L]
    stop(name, " contains a negative entry (flattened position ", hit, ").",
         domain = NA)
  }
  matrix(vals, nrow = as.integer(nrows), ncol = as.integer(ncols),
         byrow = TRUE)
}

# meta.csv as key -> raw string map; duplicate or malformed lines are
# fail-closed.
.rust_read_meta_csv <- function(path) {
  if (!file.exists(path)) {
    stop("msuiter-kernel did not produce meta.csv in the output directory.",
         domain = NA)
  }
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  meta <- list()
  for (line in lines) {
    sp <- .rust_split1(line)
    if (is.null(sp)) {
      stop(paste0("meta.csv line is not a key,value pair: ", line),
           domain = NA)
    }
    if (!is.null(meta[[sp$left]])) {
      stop(paste0("meta.csv repeats the field '", sp$left, "'."), domain = NA)
    }
    meta[[sp$left]] <- sp$right
  }
  miss <- setdiff(.RUST_META_KEYS, names(meta))
  if (length(miss) > 0L) {
    stop(paste0("meta.csv is missing the required field '", miss[1L], "'."),
         domain = NA)
  }
  meta
}

# Read, validate and cross-check everything the kernel emitted. Every branch
# names the offending file (and, for meta, the field).
.parse_outputs <- function(dir, m, n, k, iters_expect, meta_expect) {
  meta <- .rust_read_meta_csv(file.path(dir, "meta.csv"))
  for (key in names(meta_expect)) {
    asked <- as.character(meta_expect[[key]])
    got <- as.character(meta[[key]])
    if (!identical(got, asked)) {
      stop(paste0("meta.csv field '", key, "' mismatch: asked ", asked,
                  ", kernel reported ", if (is.null(got)) "<absent>" else got, "."),
           domain = NA)
    }
  }
  iterations <- suppressWarnings(as.integer(meta$iterations))
  if (is.na(iterations) || !identical(iterations, .rust_count(iters_expect, "iters"))) {
    stop(paste0("meta.csv field 'iterations' must be ", iters_expect,
                ", kernel reported ", meta$iterations, "."), domain = NA)
  }
  w <- .rust_read_matrix_csv(file.path(dir, "W.csv"), m, k, "W.csv")
  h <- .rust_read_matrix_csv(file.path(dir, "H.csv"), k, n, "H.csv")

  trace_path <- file.path(dir, "trace.csv")
  if (!file.exists(trace_path)) {
    stop("msuiter-kernel did not produce trace.csv in the output directory.",
         domain = NA)
  }
  tlines <- readLines(trace_path, warn = FALSE)
  tlines <- tlines[nzchar(tlines)]
  if (!length(tlines) || !identical(tlines[1L], "t,objective")) {
    stop("trace.csv must start with the 't,objective' header.", domain = NA)
  }
  trows <- tlines[-1L]
  want_rows <- .rust_count(iters_expect, "iters") + 1L
  if (!identical(length(trows), want_rows)) {
    stop(paste0("trace.csv must hold iterations + 1 = ", want_rows,
                " objective rows, got ", length(trows), "."), domain = NA)
  }
  tvals <- sapply(trows, function(row) {
    sp <- .rust_split1(row)
    if (is.null(sp)) {
      stop(paste0("trace.csv row is not a t,objective pair: ", row),
           domain = NA)
    }
    list(t = suppressWarnings(as.integer(sp$left)),
         obj = suppressWarnings(as.numeric(sp$right)))
  }, simplify = FALSE)
  ts <- sapply(tvals, function(x) x$t)
  obj <- as.numeric(sapply(tvals, function(x) x$obj))
  if (!identical(unname(ts), seq_len(want_rows) - 1L)) {
    stop("trace.csv iteration column must count 0, 1, ... without gaps.",
         domain = NA)
  }
  if (!all(is.finite(obj))) {
    hit <- which(!is.finite(obj))[1L]
    stop(paste0("trace.csv objective is not finite at line ", hit + 1L, "."),
         domain = NA)
  }
  list(W = w, H = h, objective = obj, meta = meta)
}

# ---------------------------------------------------------------------------
# The fit
# ---------------------------------------------------------------------------

# One kernel fit-kl run through the whole pipeline: validate -> gate ->
# resolve -> handshake -> stage -> exec -> parse. Returns list(W, H,
# objective, meta). Numeric arguments are validated strictly before the
# process is spawned; every post-spawn deviation is fail-closed. `bin` may
# be supplied by callers that resolved the kernel earlier (the parallel bp
# fan-out runs in worker processes that cannot see the parent's option).
.rust_fit_kl <- function(V, k, iters, seed = NULL, init = "nndsvd",
                         replicate = NULL, W0 = NULL, H0 = NULL,
                         bin = NULL) {
  if (!is.character(init) || length(init) != 1L || !init %in% .RUST_KERNEL_INITS) {
    stop("'init' must be one of: ", paste(.RUST_KERNEL_INITS, collapse = ", "),
         ".", domain = NA)
  }
  V <- .rust_validate_matrix(V, "V")
  m <- nrow(V)
  n <- ncol(V)
  k <- .rust_count(k, "k")
  iters <- .rust_count(iters, "iters")

  if (identical(init, "external")) {
    if (!is.null(seed)) {
      stop("the kernel rejects --seed with init = \"external\" (the run is ",
           "fully determined by W0/H0); pass seed = NULL.", domain = NA)
    }
    W0 <- .rust_validate_matrix(W0, "W0")
    H0 <- .rust_validate_matrix(H0, "H0")
    if (!identical(dim(W0), c(m, k)) || !identical(dim(H0), c(k, n))) {
      stop("init = \"external\" requires W0 (", m, " x ", k, ") and H0 (",
           k, " x ", n, ").", domain = NA)
    }
  } else {
    seed <- .rust_count(seed, "seed", allow.zero = TRUE)
    if (!is.null(W0) || !is.null(H0)) {
      stop("W0/H0 are only meaningful with init = \"external\".", domain = NA)
    }
  }
  if (!identical(init, "seeded") && !is.null(replicate)) {
    stop("--replicate is only meaningful with init = \"seeded\".", domain = NA)
  }
  if (!is.null(replicate)) {
    replicate <- .rust_count(replicate, "replicate", allow.zero = TRUE)
  }

  .require_msuiter()
  if (is.null(bin)) {
    bin <- .kernel_resolve()
  } else if (!is.character(bin) || length(bin) != 1L || !nzchar(bin)) {
    stop("'bin' must be a single non-empty path.", domain = NA)
  }
  .kernel_handshake(bin)

  in_dir <- tempfile("sigminer-rust-in-")
  out_dir <- tempfile("sigminer-rust-out-")
  on.exit({
    unlink(in_dir, recursive = TRUE)
    unlink(out_dir, recursive = TRUE)
  }, add = TRUE)

  .stage_inputs(in_dir, V, W0, H0, m, n, k)

  argv <- c("fit-kl", "--in", in_dir, "--out", out_dir,
            "--m", as.character(m), "--n", as.character(n),
            "--k", as.character(k), "--iters", as.character(iters),
            "--init", init)
  if (!identical(init, "external")) {
    argv <- c(argv, "--seed", as.character(seed))
  }
  if (identical(init, "seeded") && !is.null(replicate)) {
    argv <- c(argv, "--replicate", as.character(replicate))
  }

  res <- .kernel_exec(bin, argv)
  status <- .rust_status(res)
  if (!identical(status, 0L)) {
    stop(paste0("msuiter-kernel fit-kl exited with status ", status,
                "; stderr: ", paste(as.character(attr(res, "stderr")),
                                    collapse = " | ")), domain = NA)
  }

  meta_expect <- list(
    protocol = as.character(.RUST_PROTOCOL),
    command = "fit-kl",
    init = init,
    seed = if (identical(init, "external")) "na" else as.character(seed),
    replicate = if (identical(init, "seeded"))
      as.character(if (is.null(replicate)) 0L else replicate) else "na",
    m = as.character(m), n = as.character(n), k = as.character(k),
    iters = as.character(iters),
    iterations = as.character(iters),
    trace_len = as.character(iters + 1L)
  )
  .parse_outputs(out_dir, m, n, k, iters, meta_expect)
}

# The nrun law: run i carries seed + i - 1L. external is deterministic, so
# more than one run of it is a user error, not a wasted spawn.
.rust_runs <- function(V, k, iters, seed, nrun = 1L, init = "nndsvd",
                       replicate = NULL, W0 = NULL, H0 = NULL,
                       bin = NULL) {
  nrun <- .rust_count(nrun, "nrun")
  if (identical(init, "external") && nrun > 1L) {
    stop("init = \"external\" is fully deterministic; nrun > 1 would launch ",
         "identical runs; use nrun = 1", domain = NA)
  }
  fits <- list()
  for (i in seq_len(nrun)) {
    run_seed <- if (identical(init, "external")) NULL else seed + i - 1L
    fits[[i]] <- .rust_fit_kl(V = V, k = k, iters = iters, seed = run_seed,
                              init = init, replicate = replicate,
                              W0 = W0, H0 = H0, bin = bin)
  }
  fits
}

.rust_best_of_nrun <- function(V, k, iters, seed, init = "nndsvd",
                               nrun = 1L, replicate = NULL,
                               W0 = NULL, H0 = NULL, bin = NULL) {
  fits <- .rust_runs(V = V, k = k, iters = iters, seed = seed,
                      nrun = nrun, init = init, replicate = replicate,
                      W0 = W0, H0 = H0, bin = bin)
  objs <- as.numeric(sapply(fits, function(f) f$meta$obj_final))
  if (!all(is.finite(objs))) {
    stop("kernel meta.obj_final is not finite for every run.", domain = NA)
  }
  fits[[which.min(objs)]]
}

# ---------------------------------------------------------------------------
# Rank-survey support for sig_estimate(backend = "rust")
# ---------------------------------------------------------------------------

# Cophenetic-style stability across the nrun fits: every run assigns each
# sample (a V column) to its argmax signature from H; the pairwise agreement
# vectors are correlated across runs; a degenerate (constant) assignment
# reads as perfect stability 1.
.rust_cophenetic <- function(fits) {
  nrun <- length(fits)
  if (nrun < 2L) return(1)
  labels <- lapply(fits, function(f) apply(f$H, 2L, which.max))
  agree <- lapply(labels, function(lm) as.vector(outer(lm, lm, `==`)) * 1)
  corrs <- c()
  for (i in seq_len(nrun - 1L)) {
    for (j in (i + 1L):nrun) {
      if (all(agree[[i]] == agree[[i]][1L]) ||
          all(agree[[j]] == agree[[j]][1L])) {
        r <- 1                       # degenerate clustering: fully stable
      } else {
        r <- suppressWarnings(cor(agree[[i]], agree[[j]]))
        if (is.na(r)) r <- 1
      }
      corrs <- c(corrs, r)
    }
  }
  vmean <- mean(corrs)
  min(max(vmean, 0), 1)
}

# Per-rank rust survey used by sig_estimate(): best-of-nrun fit, rss/evar
# against the input, cophenetic stability, and the kernel's obj_final.
.survey_rust <- function(mat, ranks, nrun, seed, iters, bin = NULL) {
  rows <- list()
  fits_by_rank <- list()
  for (rk in ranks) {
    k <- .rust_count(rk, "rank")
    fits <- .rust_runs(V = mat, k = k, iters = iters, seed = seed,
                       nrun = nrun, init = "nndsvd", bin = bin)
    best <- fits[[which.min(as.numeric(sapply(fits,
      function(f) as.numeric(f$meta$obj_final))))]]
    est <- best$W %*% best$H
    rss <- sum((mat - est)^2)
    tss <- sum((mat - mean(mat))^2)
    evar <- if (tss > 0) 1 - rss / tss else 1
    rows[[as.character(rk)]] <- data.frame(
      rank = as.numeric(rk),
      nrun = .rust_count(nrun, "nrun"),
      rss = rss,
      evar = evar,
      cophenetic = .rust_cophenetic(fits),
      obj_final = as.numeric(best$meta$obj_final),
      stringsAsFactors = FALSE
    )
    fits_by_rank[[as.character(rk)]] <- fits
  }
  list(survey = do.call("rbind", rows), fits = fits_by_rank)
}
