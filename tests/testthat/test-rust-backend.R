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

rb_wired_attrs <- function(res) {
  attributes(res)[c("nrun", "method", "seed", "call_method")]
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
