# synthetic_cells.R -- committed deterministic calibration-cell builder for
# the bp K-arbitration parity fixtures (parity plan P2).
#
# PURE MODULE: definitions only, no side effects, no reference-package
# dependency. Consumed by the testthat units (which build the cells on the
# fly under the pinned root seed and assert the committed msuiter verdicts
# in msuiter_verdicts.rds) and sourced by the provenance generator
# generate_verdicts.R next to this file.
#
# Data law (D3): every number this builder produces is SYNTHETIC -- Park-
# Miller LCG draws on pinned seeds, no COSMIC, no real-data-derived values.

CAL_ROOT_SEED <- 20261015L           # the single root seed of the campaign
CAL_K_GRID <- 2:5L                   # evaluated grid per calibration cell

CAL_CHANNELS <- 96L
CAL_PEAK_MASS <- 0.55
CAL_FLANK1_MASS <- 0.25
CAL_FLANK2_MASS <- 0.125
CAL_RESID_MASS <- 0.075
CAL_SIG_MASS <- c(0.85, 0.15)
CAL_BURDEN <- 100L
CAL_BLOCK <- 12L

# Park-Miller LCG (a = 16807, m = 2^31 - 1). Every intermediate stays
# below 2^53, so each product/remainder is exactly representable in an
# IEEE-754 double: bit-reproducible across hosts without integer-vector
# tricks.
CAL_LCG_M <- 2147483647
.cal_draw_unit <- function(n, seed) {
  state <- abs(seed) %% (CAL_LCG_M - 1L) + 1L        # state in [1, m - 1]
  out <- numeric(n)
  for (i in seq_len(n)) {
    state <- (16807L * state) %% CAL_LCG_M
    out[i] <- state / CAL_LCG_M
  }
  out
}

# Three-tier design: cell i plants planted[i] signatures. Signature s owns
# a contiguous CAL_BLOCK-channel window (peak + two flanking channels) plus
# an LCG-dithered residual spread over the rest; each sample loads on
# exactly two signatures with the CAL_SIG_MASS split, the pair rotating
# deterministically across samples. Above the planted K the extra basis
# explains noise only and yields unstable extra components (veto-
# detectable); below it the held-out deviance rises, so the CV argmin
# lands on the planted K.
cal_cells <- function(root_seed = CAL_ROOT_SEED) {
  planted <- c(2L, 3L, 4L)
  n_samples <- c(36L, 42L, 48L)
  cells <- vector("list", length(planted))
  for (i in seq_along(cells)) {
    k <- planted[i]
    seed <- root_seed + 1000L * i
    W <- matrix(0, nrow = CAL_CHANNELS, ncol = k)
    rownames(W) <- paste0("ch", seq_len(CAL_CHANNELS))
    offset <- seed %% 6L
    for (s in seq_len(k)) {
      start <- ((s - 1L) * CAL_BLOCK + offset * 3L) %% (CAL_CHANNELS - 2L) + 1L
      W[start, s] <- CAL_PEAK_MASS
      W[start + 1L, s] <- CAL_FLANK1_MASS
      W[start + 2L, s] <- CAL_FLANK2_MASS
      idx <- seq_len(CAL_CHANNELS)[-c(start, start + 1L, start + 2L)]
      resid <- .cal_draw_unit(CAL_CHANNELS - 3L, seed + 77L * s)
      W[idx, s] <- CAL_RESID_MASS * resid / sum(resid)
    }
    H <- matrix(0, nrow = k, ncol = n_samples[i])
    colnames(H) <- paste0("smp", seq_len(n_samples[i]))
    u <- .cal_draw_unit(2L * n_samples[i], seed + 991L)
    ua <- as.integer(floor(u[seq_len(n_samples[i])] * 1e6))
    ub <- as.integer(floor(u[seq_len(n_samples[i]) + n_samples[i]] * 1e6))
    for (j in seq_len(n_samples[i])) {
      a <- 1L + (ua[j] %% k)
      b <- 1L + (ub[j] %% (k - 1L))
      if (b >= a) b <- b + 1L
      H[a, j] <- CAL_BURDEN * CAL_SIG_MASS[1]
      H[b, j] <- CAL_BURDEN * CAL_SIG_MASS[2]
    }
    mu <- W %*% H
    dither <- matrix(.cal_draw_unit(CAL_CHANNELS * n_samples[i], seed + 313L),
      nrow = CAL_CHANNELS)
    counts <- floor(mu + dither)
    rownames(counts) <- rownames(W)
    colnames(counts) <- colnames(H)
    cells[[i]] <- list(name = paste0("tier_k", k), planted_k = k,
      counts = counts)
  }
  cells
}
