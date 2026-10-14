# bp_k_arbitration.R: the pure K-arbitration core of the bp de-novo
# survey route (parity plan P2). Semantic port of the msuiter
# ms_select_k() arbitration (msuiter R/kselect-select.R, main @
# dad30440): CV argmin (layer 1) -> stability veto (layer 2) ->
# diagnostics only (layer 3). The core is a pure function over the
# per-rank evidence frame: it neither fits nor samples, which is what
# makes its acceptance fixtures msuiter-absent-safe (the reference
# verdicts are committed as numbers, produced on the controller host by
# fixtures/bp-k-arb/generate_verdicts.R).

# Layer-2 veto thresholds: the documented upstream defaults
# (SigProfiler sigpro.py:352-354 semantics as referenced by the msuiter
# memo; a REFERENCE, not a recalibration claim -- any change here is a
# contract change and must re-run the agreement fixtures).
.BP_VETO_AVG_MIN <- 0.8
.BP_VETO_CLUSTER_MIN <- 0.2
.BP_VETO_SUM_MIN <- 1.0

# The k-evidence column schema consumed (and returned) by the core. The
# survey layer produces the full frame; the core reads these columns.
.BP_K_EVIDENCE_REQUIRED <- c("k", "cv_test_deviance", "avg_stability",
  "min_cluster_stability", "passes_veto", "veto_reason", "is_argmin",
  "selected")

# Layer-2 flag + reason for one rank from its consensus stability stats.
# NA stability (a degenerate consensus) fails the veto conservatively.
# Vectorised over the evidence frame; returns "" when the rank passes.
.bp_veto_reason <- function(avg_stability, min_cluster_stability) {
  if (length(avg_stability) != length(min_cluster_stability)) {
    stop("`avg_stability` and `min_cluster_stability` must have equal length",
      domain = NA)
  }
  if (!length(avg_stability)) {
    return(character(0))
  }
  reasons <- character(length(avg_stability))
  for (i in seq_along(reasons)) {
    avg <- avg_stability[i]
    mn <- min_cluster_stability[i]
    if (length(avg) != 1L || is.na(avg) || length(mn) != 1L || is.na(mn)) {
      reasons[i] <- "stability statistics unavailable"
      next
    }
    hit <- character(0)
    if (avg < .BP_VETO_AVG_MIN) {
      hit <- c(hit, sprintf("avg_stability %.4g < %.1f", avg, .BP_VETO_AVG_MIN))
    }
    if (mn < .BP_VETO_CLUSTER_MIN) {
      hit <- c(hit, sprintf("min_cluster_stability %.4g < %.1f", mn,
        .BP_VETO_CLUSTER_MIN))
    }
    if (avg + mn < .BP_VETO_SUM_MIN) {
      hit <- c(hit, sprintf("avg + min_cluster_stability %.4g < %.1f",
        avg + mn, .BP_VETO_SUM_MIN))
    }
    reasons[i] <- if (length(hit)) paste(hit, collapse = "; ") else ""
  }
  reasons
}

# Layer-1 argmin mark over the finite CV entries of one grid: the first
# minimum, i.e. the lowest rank on ties (which.min semantics, matching
# the reference implementation's documented intent). Non-finite entries
# (NA -- no finite fold -- or Inf) never carry the mark; an all-NON-finite
# curve marks nothing (the fallback scan then decides).
.bp_cv_argmin <- function(cv) {
  finite <- which(is.finite(cv))
  if (!length(finite)) {
    return(integer(0))
  }
  finite[which.min(cv[finite])]
}

# Validate the evidence frame at the boundary; KISS-friendly messages all
# name the offending part so a survey-layer bug is locatable from the
# error text alone.
.bp_k_evidence_check <- function(ev) {
  if (!is.data.frame(ev)) {
    stop("the K evidence must be a data.frame of per-rank rows", domain = NA)
  }
  missing <- setdiff(.BP_K_EVIDENCE_REQUIRED, names(ev))
  if (length(missing)) {
    stop(sprintf("the K evidence frame is missing required column(s): %s",
      paste0(dQuote(missing), collapse = ", ")), domain = NA)
  }
  if (nrow(ev) < 1L) {
    stop("the K evidence frame must carry at least one rank row", domain = NA)
  }
  if (!is.numeric(ev$k) || anyNA(ev$k) || !all(is.finite(ev$k)) ||
      !all(ev$k == floor(ev$k)) || !all(ev$k >= 1L)) {
    stop("the evidence `k` column must be positive whole numbers", domain = NA)
  }
  if (any(duplicated(ev$k))) {
    stop("the evidence `k` column must be unique (one row per rank)", domain = NA)
  }
  if (!identical(ev$k, sort(ev$k))) {
    stop("the evidence `k` column must be strictly increasing", domain = NA)
  }
  for (col in c("passes_veto", "is_argmin", "selected")) {
    if (!is.logical(ev[[col]]) || anyNA(ev[[col]])) {
      stop(sprintf("the evidence `%s` column must be logical (no NA)", col),
        domain = NA)
    }
  }
  if (!is.character(ev$veto_reason) || anyNA(ev$veto_reason)) {
    stop("the evidence `veto_reason` column must be character (no NA)",
      domain = NA)
  }
  if (sum(ev$is_argmin) > 1L) {
    stop("at most one evidence row may carry the `is_argmin` mark", domain = NA)
  }
  invisible(ev)
}

# The three-layer arbitration over a validated evidence frame.
#
#   cv     layer-1 switch: FALSE erases the argmin mark (layer 1 absent)
#          and the fallback scan degenerates to plain rank order.
#   mode   "explicit" -- no rank passing raises (the standalone surface;
#          full evidence frame is what the caller passed in and stays
#          available to it); "degrade" -- the default-pipeline shape:
#          fall back to the argmin (or the lowest rank) with a warning,
#          never an error, matching the reference pipeline behaviour where
#          the user's K is fixed and the arbitration only annotates.
#
# Returns a NEW frame (the input is never mutated) with is_argmin and
# selected normalised: exactly one selected row, the layer-2 flag respected
# throughout.
.bp_k_arbitrate <- function(ev, cv = TRUE, mode = c("explicit", "degrade")) {
  mode <- match.arg(mode, c("explicit", "degrade"))
  if (!is.logical(cv) || length(cv) != 1L || is.na(cv)) {
    stop("`cv` must be a single TRUE or FALSE", domain = NA)
  }
  .bp_k_evidence_check(ev)
  ev <- ev[seq_len(nrow(ev)), , drop = FALSE]   # fresh copy; never mutate input
  if (!isTRUE(cv)) {
    ev$is_argmin <- rep(FALSE, nrow(ev))
  }
  # Fallback scan order: CV-error ascending (NA = no finite fold last,
  # lower rank first on ties), or plain rank order when layer 1 is absent.
  ord <- if (isTRUE(cv)) {
    order(ev$cv_test_deviance, ev$k, na.last = TRUE)
  } else {
    order(ev$k)
  }
  ai <- which(ev$is_argmin)
  passing <- which(ev$passes_veto)
  sel_idx <- NULL
  if (length(ai) == 1L && ev$passes_veto[ai]) {
    sel_idx <- ai
  } else if (length(passing)) {
    in_ord <- passing[order(match(passing, ord))]
    sel_idx <- in_ord[1L]
  }
  if (is.null(sel_idx)) {
    reasons <- paste(sprintf("k=%d: %s", ev$k, ev$veto_reason), collapse = "; ")
    if (identical(mode, "explicit")) {
      stop(sprintf(paste0("K arbitration is undecidable: no rank in the grid ",
        "passes the stability veto (layer-2 rule: avg >= 0.8, min_cluster >= ",
        "0.2, sum >= 1.0). Veto reasons: ", "%s"), reasons), domain = NA)
    }
    fallback <- if (length(ai)) ai[1L] else ord[1L]
    warning(sprintf(paste0("K arbitration is undecidable: no rank passes the ",
      "stability veto; falling back to k = %d (evidence frame retained in ",
      "the survey result). Veto reasons: %s"), ev$k[fallback], reasons),
      domain = NA, call. = FALSE)
    sel_idx <- fallback
  }
  ev$selected <- rep(FALSE, nrow(ev))
  ev$selected[sel_idx] <- TRUE
  ev
}
