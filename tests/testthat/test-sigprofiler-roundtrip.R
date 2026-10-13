# Export -> import round-trip tests for the SigProfiler interop bridge.
#
# Two legs, both deterministic and python-free:
#  1. Export leg: sigprofiler_extract(output_matrix_only = TRUE) must write the
#     sigprofiler_input.txt hand-off file without touching the python /
#     reticulate / conda stack (the matrix-file mode is pure R).
#  2. Hand-off contract: the exported matrix and the committed result fixture
#     used by test-sigprofiler-import.R share one tiny universe of feature and
#     sample names and numerics, so export -> (SigProfiler run) -> import is
#     pinned end to end against known values.

sp_features <- c("A[C]A", "A[C]C", "A[T]A", "A[T]C")
sp_samples <- c("Sample1", "Sample2", "Sample3")

# features x samples catalogue (columns are the samples: Sample1 carries
# 1/3/2/2 across the four features, Sample2/3 carry 1/1/1/1, Sample4 is the
# null-catalogue column sigprofiler_extract() is documented to drop).
# sigprofiler_extract() takes the NMF orientation: SAMPLES as rows,
# mutational features as columns (it transposes internally, yielding the
# "Mutation Types" row-labelled export). Powers-of-two sums keep every
# derived quantity binary-exact.
sp_nmf_matrix <- function() {
  fm <- matrix(
    c(
      1, 3, 2, 2, # Sample1 across the four features
      1, 1, 1, 1, # Sample2
      1, 1, 1, 1, # Sample3
      0, 0, 0, 0  # Sample4: null catalogue, dropped on export
    ),
    nrow = 4L, ncol = 4L,
    dimnames = list(sp_features,
                    c(sp_samples, "Sample4"))
  )
  t(fm)
}

sp_fixture_activities_samples <- function() {
  act <- data.table::fread(
    file.path(testthat::test_path(), "fixtures", "sigprofiler", "SBS96",
              "run_001",
              "Suggested_Solution", "SBS96_De-Novo_Solution", "Activities.txt"),
    sep = "\t"
  )
  act[[1]]
}

test_that("output_matrix_only = TRUE writes the hand-off file without the python stack", {
  # The matrix-file mode must not reach env_install()/reticulate at all:
  # mock the python gate with a sentinel that aborts on contact.
  ns <- asNamespace("sigminer")
  original <- get("env_install", envir = ns)
  assignInNamespace("env_install",
                    function(...) stop("python machinery invoked by matrix-only export"),
                    ns)
  on.exit(assignInNamespace("env_install", original, ns), add = TRUE)

  out_dir <- file.path(tempdir(), "p0-export-out")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE)

  m <- sp_nmf_matrix()
  m_before <- m

  msgs <- capture.output(
    type = "message",
    expr = expect_no_error(
      res <- sigprofiler_extract(m, out_dir, output_matrix_only = TRUE)
    )
  )
  expect_null(res)

  in_file <- file.path(out_dir, "sigprofiler_input.txt")
  expect_true(file.exists(in_file))

  # Null-catalogue sample is reported and dropped.
  expect_match(paste0(msgs, collapse = "\n"), "dropped due to null catalogue")
  expect_match(paste0(msgs, collapse = "\n"), "Sample4")

  written <- data.table::fread(in_file, sep = "\t")
  expect_equal(colnames(written), c("Mutation Types", sp_samples))
  expect_equal(as.character(written[["Mutation Types"]]), sp_features)
  df <- as.data.frame(written)
  num <- matrix(as.numeric(as.matrix(df[, sp_samples, drop = FALSE])),
                nrow = length(sp_features))
  expect_equal(num, unname(as.matrix(t(m)[, sp_samples, drop = FALSE])))

  # Immutability: the user-passed matrix must come back untouched.
  expect_identical(m, m_before)
})

test_that("export hand-off shares one name universe with the import fixture tree", {
  out_dir <- file.path(tempdir(), "p0-export-contract")
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE)

  m <- sp_nmf_matrix()
  suppressMessages(sigprofiler_extract(m, out_dir, output_matrix_only = TRUE))

  written <- data.table::fread(file.path(out_dir, "sigprofiler_input.txt"),
                                sep = "\t")

  # Sample names as the fixture result tree carries them.
  fixture_samples <- sp_fixture_activities_samples()
  expect_equal(as.character(written[["Mutation Types"]]), sp_features)
  expect_equal(colnames(written)[-1L], as.character(fixture_samples))

  # Feature names match the fixture Signatures.txt components.
  sig_tbl <- data.table::fread(
    file.path(testthat::test_path(), "fixtures", "sigprofiler", "SBS96",
              "run_001",
              "Suggested_Solution", "SBS96_De-Novo_Solution", "Signatures.txt"),
    sep = "\t"
  )
  expect_equal(as.character(sig_tbl[[1]]), as.character(written[["Mutation Types"]]))
})

test_that("fixture import numerics line up with the exported matrix (deterministic round-trip)", {
  res <- suppressMessages(
    sigprofiler_import(
      file.path(testthat::test_path(), "fixtures", "sigprofiler", "SBS96"),
      type = "suggest"
    )
  )

  m <- sp_nmf_matrix()
  # Row (signature) totals equal the fixture activities row sums:
  # Sig1: 1 + 2 + 4 = 7, Sig2: 3 + 2 + 0 = 5.
  expect_equal(unname(rowSums(res$solution$Exposure)), c(7, 5))

  # Hand-off catalogue sample column totals: 1+3+2+2, 1+1+1+1, 1+1+1+1.
  col_tot <- colSums(as.matrix(t(m)[, sp_samples, drop = FALSE]))
  expect_equal(unname(col_tot), c(8, 4, 4))

  # Shipped rescale convention: Signature col j = Signature.norm col j times
  # rowSums(Exposure)[j]; fixture columns are (1,3,2,2) and (1,1,1,1).
  expect_equal(unname(colSums(res$solution$Signature)), c(8 * 7, 4 * 5))
})
