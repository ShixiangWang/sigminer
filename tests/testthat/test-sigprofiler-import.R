# Tests for the sigprofiler_import() reading surface over a committed,
# hand-authored fixture tree that mirrors the current SigProfilerExtractor
# output tree shape (Suggested_Solution/<MT>_De-Novo_Solution plus
# Suggested_Solution/COSMIC_<MT>_Decomposed_Solution plus All_Solutions).
#
# The refit leg pins the ruling from the P0 source-truth check (2026-10-13):
# the current upstream sigpro.py HEAD has no `refit_denovo_signatures` keyword
# (it lives only in a stale docstring), but decomposition runs unconditionally
# (`denovo_refit_option = True` at the spa_analyze call), so
# COSMIC_<MT>_Decomposed_Solution is always present and the
# `type = "refit"` import branch remains alive against the real tree.

sp_fixture_tree <- function() {
  file.path(testthat::test_path(), "fixtures", "sigprofiler")
}

sp_fixture_sbs96 <- function() {
  # The output directory passed to sigprofiler_import(); the extractor run
  # lives in the single sub-directory below it.
  file.path(sp_fixture_tree(), "SBS96")
}

test_that("sigprofiler_import reads the suggested (de-novo) solution from the fixture tree", {
  res <- sigprofiler_import(sp_fixture_sbs96(), type = "suggest")

  expect_equal(class(res$solution), "Signature")
  expect_equal(attr(res$solution, "call_method"), "SigProfiler")

  # Exposure is signatures x samples, renamed to Sig1..SigK.
  expect_equal(dim(res$solution$Exposure), c(2L, 3L))
  expect_equal(rownames(res$solution$Exposure), c("Sig1", "Sig2"))
  expect_equal(colnames(res$solution$Exposure), c("Sample1", "Sample2", "Sample3"))
  expect_equal(unname(res$solution$Exposure["Sig1", ]), c(1L, 2L, 4L))
  expect_equal(unname(res$solution$Exposure["Sig2", ]), c(3L, 2L, 0L))

  # Column-normalized exposures sum to one per sample.
  expect_equal(unname(colSums(res$solution$Exposure.norm)), c(1, 1, 1))

  # Signature.norm carries the raw fixture component columns; Signature is
  # the shipped rescale convention (column j times rowSums(Exposure)[j]).
  expect_equal(unname(res$solution$Signature.norm[, "Sig1"]), c(1L, 3L, 2L, 2L))
  expect_equal(unname(res$solution$Signature.norm[, "Sig2"]), c(1L, 1L, 1L, 1L))
  expect_equal(unname(res$solution$Signature[, "Sig1"]), c(7, 21, 14, 14))
  expect_equal(unname(res$solution$Signature[, "Sig2"]), c(5, 5, 5, 5))

  expect_equal(as.character(res$solution$Stats$signatures$Signatures),
               c("Sig1", "Sig2"))
  expect_equal(as.character(res$solution$Stats$samples$Samples),
               c("Sample1", "Sample2", "Sample3"))

  expect_equal(as.integer(res$all_stats$SignatureNumber), c(2L, 3L))
})

test_that("sigprofiler_import type = 'refit' reads the always-emitted COSMIC decomposed solution", {
  res <- sigprofiler_import(sp_fixture_sbs96(), type = "refit")

  # Ruling pin: the branch must keep working against the current upstream
  # output tree (COSMIC_<MT>_Decomposed_Solution, map_to_COSMIC stats file).
  expect_equal(class(res$solution), "Signature")
  expect_equal(attr(res$solution, "call_method"), "SigProfiler")
  expect_equal(dim(res$solution$Exposure), c(2L, 3L))
  expect_equal(as.character(res$solution$Stats$signatures$Signatures),
               c("Sig1", "Sig2"))
  # The map_to_COSMIC fallback carries the de-novo-to-COSMIC mapping rows,
  # one row per signature.
  expect_equal(nrow(res$solution$Stats$signatures), 2L)
})

test_that("sigprofiler_import type = 'all' labels solutions from the All_Solutions basenames", {
  res <- sigprofiler_import(sp_fixture_sbs96(), type = "all")

  # Solution labels follow the documented paste0("S", ...) convention.
  expect_equal(sort(names(res$solution_list)), c("SID83", "SSBS96"))
  for (nm in names(res$solution_list)) {
    expect_equal(class(res$solution_list[[nm]]), "Signature")
    expect_equal(dim(res$solution_list[[nm]]$Exposure), c(2L, 3L))
  }
  expect_equal(as.integer(res$all_stats$SignatureNumber), c(2L, 3L))
})

test_that("sigprofiler_import validates its input fail-fast", {
  missing_dir <- file.path(tempdir(), "no-such-sigprofiler-dir")
  expect_error(sigprofiler_import(missing_dir))

  empty_dir <- file.path(tempdir(), "sigprofiler-empty")
  if (!dir.exists(empty_dir)) dir.create(empty_dir)
  expect_error(sigprofiler_import(empty_dir), "Cannot find any result directory")

  two_dirs <- file.path(tempdir(), "sigprofiler-two-dirs")
  if (!dir.exists(two_dirs)) dir.create(two_dirs)
  dir.create(file.path(two_dirs, "A"))
  dir.create(file.path(two_dirs, "B"))
  expect_error(sigprofiler_import(two_dirs), "more than 1")

  expect_error(sigprofiler_import(sp_fixture_sbs96(), type = "bogus"))
})

test_that("sigprofiler_import leaves the fixture tree byte-identical (read-only import)", {
  if (!requireNamespace("digest", quietly = TRUE)) {
    skip("package 'digest' not installed")
  }
  files <- list.files(sp_fixture_tree(), recursive = TRUE, full.names = TRUE)
  digest_one <- function(f) digest::digest(file = f, algo = "sha256", serialize = FALSE)
  before <- vapply(files, digest_one, character(1L))

  invisible(sigprofiler_import(sp_fixture_sbs96(), type = "all"))
  invisible(sigprofiler_import(sp_fixture_sbs96(), type = "suggest"))

  after <- vapply(files, digest_one, character(1L))
  expect_equal(after, before)
})
