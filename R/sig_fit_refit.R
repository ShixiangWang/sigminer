## ---------------------------------------------------------------------------
## Internal refit decision layers for sig_fit(method = "likelihood") and
## sig_fit(method = "lrt").
##
## Clean-room law: the decision semantics below are transcribed from the
## primary-source captures in /tmp/p3refit/logs/00-primary-source-gate.md
## (MuSiCal Fig.5 bidirectional stepwise over the multinomial likelihood;
## PASA Step-1 presence tests, Step-2 stop-LRT and minimality at the alpha
## significance level, chi-squared statistic with |A1| - |A2| degrees of
## freedom). No external implementation code was consulted or ported.
##
## The refit engine (exposures given a selected set) deliberately reuses the
## package's own QP simplex route -- identical to what method = "QP" drives
## through quadprog::solve.QP -- so the decision layer is the only delta
## versus the legacy route. This reuse is an engineering choice and is NOT
## claimed to be MuSiCal/PASA implementation identity.
## ---------------------------------------------------------------------------

## Lane-calibrated acceptance constants. The MuSiCal paper states its
## likelihood thresholds are study-optimized parameters ("two likelihood
## thresholds were optimized", "optimal threshold was chosen through a grid
## search") and PASA keeps alpha as its significance-level parameter; the
## values here are therefore this lane's calibration (fixture evidence:
## /tmp/p3refit/logs/03b-nullcal-probe.log, 03c-stepwise-calibration.log),
## not paper-sourced constants, and every one is caller-overridable.
P3_LL_ADD <- 8.0
P3_LL_KEEP <- 8.0
P3_ALPHA <- 0.01
## Numerical floor for log() arguments in the multinomial likelihood.
P3_P_MIN <- 1e-12
## Termination guard for the bidirectional sweeps (see the hysteresis note in
## run_stepwise): cycles are impossible when ll_add >= ll_keep, the guard only
## bounds pathological caller-chosen threshold pairs.
P3_MAX_ROUNDS <- 64L
## Upper bound on the size of the Step-1 survivor set enumerated exhaustively
## by the LRT minimal-set search; larger survivor sets fall back to the
## greedy-argmax construction (documented engineering guard, no identity
## claim -- see find_minimal_set).
P3_MAX_V <- 12L

## log-likelihood + exposures of the best QP-simplex reconstruction of count
## vector d by the columns of Psel. Returns NULL when the QP solve fails.
## The empty column set reconstructs the uniform reference profile (documented
## implementation choice: the paper defines the nested-model comparisons, the
## empty-set reference profile is not paper-sourced).
.p3_fit_one <- function(d, Psel) {
  K <- length(d)
  if (is.null(Psel) || ncol(Psel) == 0) {
    p <- rep(1 / K, K)
    return(list(expo = numeric(0), logL = sum(d * log(pmax(p, P3_P_MIN)))))
  }
  total <- sum(d)
  m <- d / total
  N <- ncol(Psel)
  G <- t(Psel) %*% Psel
  Cmat <- cbind(rep(1, N), diag(N))
  b <- c(1, rep(0, N))
  dv <- t(m) %*% Psel
  out <- tryCatch(
    quadprog::solve.QP(Dmat = G, dvec = dv, Amat = Cmat, bvec = b, meq = 1),
    error = function(e) NULL
  )
  if (is.null(out)) {
    return(NULL)
  }
  e <- out$solution
  e[e < 0] <- 0
  if (sum(e) <= 0) {
    ## Degenerate solve: the model collapses onto the zero profile. Treat it
    ## as the uniform reference for comparison purposes.
    p <- rep(1 / K, K)
    return(list(expo = rep(0, N), logL = sum(d * log(pmax(p, P3_P_MIN)))))
  }
  e <- e / sum(e)
  p <- as.numeric(Psel %*% e)
  p <- p / sum(p)
  list(expo = e, logL = sum(d * log(pmax(p, P3_P_MIN))))
}

## Cached subset fitting over one sample: key is the sorted column-name
## string. Names of the reference matrix columns are unique by construction.
.p3_make_fitter <- function(d, P) {
  cache <- new.env()
  ## Control character that cannot occur in a column name, so the sentinel key
  ## of the empty set can never collide with a joined name of a real set.
  empty.key <- paste0("\r.empty\r")
  function(cols) {
    key <- if (length(cols) == 0) {
      ## The empty column set is a legitimate candidate (the uniform reference
      ## profile); the "+"-join of character(0) is not a length-1 name, so it
      ## needs a dedicated sentinel key.
      empty.key
    } else {
      paste0(cols[order(cols)], collapse = "+")
    }
    ## Floor-safe cache probe (base get() with inherits = FALSE since the
    ## dawn of R): cached fits are never NULL, so a failed lookup is the
    ## only NULL producer here and unambiguously means miss.
    hit <- tryCatch(get(key, envir = cache, inherits = FALSE),
      error = function(e) NULL)
    if (!is.null(hit)) {
      return(hit)
    }
    fit <- if (length(cols) == 0) {
      .p3_fit_one(d, NULL)
    } else {
      .p3_fit_one(d, P[, cols, drop = FALSE])
    }
    if (is.null(fit)) {
      ## Unfit-able subset (singular solver): poison-pill it so comparisons
      ## never select it.
      fit <- list(expo = rep(0, length(cols)), logL = -Inf)
    }
    assign(key, fit, envir = cache)
    fit
  }
}

## LRT p-value per the paper: lambda = -2(log P(D|A2) - log P(D|A1)) with A2
## a proper subset of A1, distributed chi-squared with |A1| - |A2| degrees of
## freedom. df here is a signature-COUNT difference (the paper's own rule).
.p3_lrt_p <- function(logL_big, logL_small, df_sig) {
  if (!is.finite(logL_big) || !is.finite(logL_small)) {
    return(1)
  }
  lambda <- 2 * (logL_big - logL_small)
  if (lambda <= 0 || df_sig <= 0) {
    return(1)
  }
  stats::pchisq(lambda, df = df_sig, lower.tail = FALSE)
}

## Shared tail: refit the finally selected set and hand off to return_expo.
## The returned vector stays UNNAMED on purpose: sig_fit assembles the sample
## results with dplyr::bind_rows(), which binds an unnamed per-sample vector
## as one column and a named one as one row (decompose_QP returns the same
## unnamed contract).
.p3_finish <- function(x, y, z, P, kept, type) {
  expo <- rep(0, ncol(P))
  if (length(kept) > 0) {
    fit <- .p3_fit_one(x, P[, kept, drop = FALSE])
    if (!is.null(fit)) {
      expo[match(kept, colnames(P))] <- fit$expo
    }
  }
  return_expo(expo = expo, y, type, total = sum(x))
}

## ---------------------------------------------------------------------------
## method = "likelihood": MuSiCal-branch semantics.
## "Matching and refitting with multinomial likelihood-based sparse NNLS",
## "bidirectional stepwise algorithm". Acceptance runs on the absolute gain of
## the multinomial log-likelihood; the two thresholds are parameters (paper:
## they are optimized per study), exposed as ll_add / ll_keep.
##
## Start from the full candidate set; alternate a backward sweep (drop the
## single member whose removal costs less than ll_keep, one at a time so each
## cost is re-measured against the current set) and a forward sweep (re-add
## the single best outsider when its gain reaches ll_add) until stable.
## The strict-complementary gates make the pair (ll_add = ll_keep) hysteresis
## stable; unequal caller pairs are bounded by P3_MAX_ROUNDS.
## ---------------------------------------------------------------------------
run_stepwise <- function(d, P, ll_add, ll_keep) {
  all_cols <- colnames(P)
  if (sum(d) == 0) {
    return(character(0))
  }
  fit <- .p3_make_fitter(d, P)
  kept <- all_cols
  rounds <- 0L
  repeat {
    rounds <- rounds + 1L
    if (rounds > P3_MAX_ROUNDS) {
      break()
    }
    changed <- FALSE
    ## Backward sweep: evaluate removal costs against the current set and
    ## drop the cheapest sub-threshold member, then re-measure.
    repeat {
      costs <- numeric(0)
      for (k in kept) {
        costs <- c(costs, fit(kept)$logL - fit(kept[kept != k])$logL)
      }
      if (length(costs) == 0) {
        break()
      }
      worst <- which.min(costs)
      if (length(worst) == 0 || is.na(costs[worst]) || costs[worst] >= ll_keep) {
        break()
      }
      kept <- kept[-worst]
      changed <- TRUE
    }
    ## Forward sweep: re-admit the best outsider when its improvement
    ## reaches ll_add.
    outsiders <- setdiff(all_cols, kept)
    if (length(outsiders) > 0) {
      base <- fit(kept)$logL
      gains <- numeric(0)
      for (k in outsiders) {
        gains <- c(gains, fit(c(kept, k))$logL - base)
      }
      j <- which.max(gains)
      if (length(j) > 0 && !is.na(gains[j]) && gains[j] >= ll_add) {
        kept <- c(kept, outsiders[j])
        changed <- TRUE
      }
    }
    if (!changed) {
      break()
    }
  }
  kept
}

## ---------------------------------------------------------------------------
## method = "lrt": PASA-branch semantics.
## Step 1: "the signature presence tests" -- every signature k of the full
## set is tested as (i) full-minus-k versus (ii) full; k survives when the
## LRT rejects the equality null, i.e. "the signatures in A1 - A2 are
## plausibly needed". The test is against all other signatures, so the
## survivor set "does not depend on the order of testing".
## Step 2: return a set A that is "good enough" against the full set
## (stop-LRT, p >= alpha) and minimal at alpha (removing any member drops
## significance, p < alpha). Satisfying the paper's own definition directly
## over the Step-1 survivor set is order-independent by construction; the
## search enumerates survivor subsets by increasing size (bounded by
## P3_MAX_V, greedy argmax fallback beyond -- a documented engineering
## guard, not a paper element).
## ---------------------------------------------------------------------------
find_minimal_set <- function(d, P, alpha) {
  all_cols <- colnames(P)
  if (sum(d) == 0) {
    return(character(0))
  }
  fit <- .p3_make_fitter(d, P)
  full <- fit(all_cols)
  ## Step 1: presence tests.
  V <- character(0)
  for (k in all_cols) {
    p <- .p3_lrt_p(full$logL, fit(all_cols[all_cols != k])$logL, df_sig = 1L)
    if (p < alpha) {
      V <- c(V, k)
    }
  }
  if (length(V) == 0) {
    return(character(0))
  }
  ## Step 2: minimal set inside V that is statistically as good as the full
  ## model. Enumerate subsets of V by increasing size, take the first one
  ## that qualifies.
  good <- NULL
  if (length(V) <= P3_MAX_V) {
    subs <- lapply(seq_len(2^length(V)), function(m) {
      ## Floor-safe (DESCRIPTION floor R >= 3.5) bitmask extraction: plain
      ## integer arithmetic, no dependence on the bit-operator base family
      ## that only exists in much later base versions.
      bits <- ((m - 1L) %/% 2^(seq_len(length(V)) - 1L)) %% 2 == 1L
      V[bits]
    })
  } else {
    ## Greedy fallback (guard only): grow from the empty set by the single
    ## best improvement, then prune. Documented engineering path, never taken
    ## by fixtures of this lane.
    subs <- list(character(0))
    cur <- character(0)
    for (k in V) {
      cand <- c(cur, k)
      if (fit(cand)$logL > fit(cur)$logL) {
        cur <- cand
      }
      subs <- c(subs, list(cur))
    }
  }
  subs <- subs[order(vapply(subs, length, 1L))]
  for (A in subs) {
    fitA <- fit(A)
    ## (a) stop condition: the current reconstruction is good enough against
    ## the full candidate set.
    if (.p3_lrt_p(full$logL, fitA$logL, df_sig = ncol(P) - length(A)) < alpha) {
      next
    }
    ## (b) minimality: removing any member must break significance.
    minimal <- TRUE
    for (k in A) {
      p <- .p3_lrt_p(fitA$logL, fit(A[A != k])$logL, df_sig = 1L)
      if (p >= alpha) {
        minimal <- FALSE
        break()
      }
    }
    if (!minimal) {
      next
    }
    ## Subsets are scanned by increasing size and (a) is monotone in A, so
    ## the first qualifying minimal set in that order is the answer. Ties
    ## inside a size class follow the enumeration order (= reference column
    ## order), which keeps the choice deterministic.
    good <- A
    break
  }
  if (is.null(good)) {
    ## No subset cleared both criteria (can only happen under the greedy
    ## guard path): the presence-test survivors are the honest answer.
    good <- V
  }
  good
}

## ---------------------------------------------------------------------------
## Per-sample entry points, same positional contract as decompose_QP:
##   x = catalogue column, y = rel_threshold, z = sample label, P = reference
## ---------------------------------------------------------------------------
decompose_LIKELIHOOD <- function(x, y, z, P, type = "absolute", auto_reduce = FALSE,
                                 ll_add = P3_LL_ADD, ll_keep = P3_LL_KEEP, ...) {
  if (auto_reduce) {
    send_stop(
      "Method 'likelihood' does not support 'auto_reduce': ",
      "the likelihood-based bidirectional stepwise decision layer already ",
      "selects the reference set."
    )
  }
  if (!is.na(z)) {
    send_info("Fitting sample: ", z)
  }
  if (ll_add < ll_keep) {
    send_warning(
      "'ll_add' <", ll_add, " is below 'll_keep' =", ll_keep,
      "; the stepwise hysteresis band is open, the termination guard applies."
    )
  }
  kept <- run_stepwise(x, P, ll_add = ll_add, ll_keep = ll_keep)
  .p3_finish(x, y, z, P, kept, type)
}

decompose_LRT <- function(x, y, z, P, type = "absolute", auto_reduce = FALSE,
                           alpha = P3_ALPHA, ...) {
  if (auto_reduce) {
    send_stop(
      "Method 'lrt' does not support 'auto_reduce': ",
      "the likelihood-ratio presence/stop tests already select the ",
      "reference set."
    )
  }
  if (!is.na(z)) {
    send_info("Fitting sample: ", z)
  }
  kept <- find_minimal_set(x, P, alpha = alpha)
  .p3_finish(x, y, z, P, kept, type)
}
