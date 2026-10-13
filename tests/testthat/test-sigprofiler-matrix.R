# Version-matrix tests for the SigProfiler interop bridge.
#
# The command/runner builders are deterministic string tests that always run.
# The actual container execution is local-only and double-gated: it needs
# SIGMINER_SP_MATRIX=1 in the environment AND a docker binary on PATH;
# otherwise the execution test skips cleanly so R CMD check stays green.

test_that("the matrix is parameterized over the ratified SigProfilerExtractor versions", {
  expect_equal(SIGMINER_SP_MATRIX_VERSIONS, c("1.1.3", "1.2.8", "1.4.0"))
})

test_that("image references are routed through the project registry mirror", {
  old_host <- Sys.getenv("SIGMINER_SP_MIRROR_HOST", "")
  old_set <- nzchar(old_host)
  if (old_set) Sys.unsetenv("SIGMINER_SP_MIRROR_HOST")
  on.exit({
    if (old_set) {
      Sys.setenv(SIGMINER_SP_MIRROR_HOST = old_host)
    } else {
      Sys.unsetenv("SIGMINER_SP_MIRROR_HOST")
    }
  }, add = TRUE)

  expect_equal(
    sigminer_sp_image_ref(),
    "9f6c5in3dh5vszgdz0.xuanyuan.run/library/python:3.11-slim"
  )

  Sys.setenv(SIGMINER_SP_MIRROR_HOST = "mirror.example")
  expect_equal(
    sigminer_sp_image_ref(),
    "mirror.example/library/python:3.11-slim"
  )
})

test_that("matrix command builders fail closed on shell-hostile inputs", {
  expect_error(sigminer_sp_validate_version("1.2.8; rm -rf /"))
  expect_error(sigminer_sp_validate_version("1.2.*"))
  expect_error(sigminer_sp_validate_token("a b", "workdir"))
  expect_error(sigminer_sp_validate_token("a$b", "workdir"))
  expect_error(sigminer_sp_validate_token("a`b", "workdir"))
  expect_error(sigminer_sp_validate_token("a'b", "workdir"))
  expect_error(sigminer_sp_validate_token("a\"b", "workdir"))
  expect_error(sigminer_sp_validate_token("100%", "workdir"))
  expect_error(sigminer_sp_validate_token("a\nb", "workdir"))
  expect_error(sigminer_sp_validate_token("", "workdir"))
  expect_error(
    sigminer_sp_matrix_command("1.2.8", workdir = "/tmp/has space")
  )
})

test_that("matrix command pins the pip install and mounts the shared workdir", {
  cmd <- sigminer_sp_matrix_command("1.2.8", workdir = "/tmp/sp-matrix")
  expect_match(cmd, "^docker run --rm -v /tmp/sp-matrix:/sp ")
  expect_match(cmd, "9f6c5in3dh5vszgdz0\\.xuanyuan\\.run/library/python:3\\.11-slim")
  expect_match(cmd, "SigProfilerExtractor==1\\.2\\.8")
  expect_match(cmd, "python /sp/run_cell\\.py\"$")
})

test_that("the in-container runner uses only fixed /sp/ literals", {
  wd <- file.path(tempdir(), "sp-runner-check")
  if (!dir.exists(wd)) dir.create(wd, recursive = TRUE)
  path <- sigminer_sp_write_runner(wd)
  expect_equal(path, file.path(wd, "run_cell.py"))
  lines <- readLines(path)
  joined <- paste0(lines, collapse = "\n")
  expect_match(joined, "from SigProfilerExtractor import sigpro as sig")
  expect_match(joined, "sig\\.sigProfilerExtractor\\(")
  expect_match(joined, "'/sp/out'")
  expect_match(joined, "'/sp/sigprofiler_input\\.txt'")
  # No host path ever leaks into the generated script.
  expect_equal(sum(grepl(wd, lines, fixed = TRUE)), 0L)
})

test_that("matrix execution runs only under the explicit opt-in flag", {
  if (!identical(Sys.getenv("SIGMINER_SP_MATRIX", ""), "1")) {
    skip(paste0(
      "SigProfilerExtractor version matrix is local-only: set ",
      "SIGMINER_SP_MATRIX=1 (see README-sigprofiler-matrix.md)"
    ))
  }
  docker <- Sys.which("docker")
  if (identical(docker, "")) skip("docker not available on PATH")

  for (version in SIGMINER_SP_MATRIX_VERSIONS) {
    wd <- file.path(tempdir(), paste0("sp-matrix-", gsub("\\.", "_", version)))
    if (dir.exists(wd)) unlink(wd, recursive = TRUE)
    dir.create(wd, recursive = TRUE)

    # Export the shared fixture input with the pure-R hand-off mode.
    m <- sp_nmf_matrix_export_helper()
    suppressMessages(
      sigprofiler_extract(m, wd, output_matrix_only = TRUE)
    )
    sigminer_sp_write_runner(wd)

    cmd <- sigminer_sp_matrix_command(version, workdir = wd)
    rc <- system(cmd, ignore.stderr = FALSE)
    expect_equal(rc, 0L, label = sprintf("SigProfilerExtractor %s ran clean", version))

    res <- suppressMessages(sigprofiler_import(file.path(wd, "out"),
                                               type = "suggest"))
    expect_equal(class(res$solution), "Signature",
                 label = sprintf("SigProfilerExtractor %s imports as Signature", version))
  }
})

# Helper for the opt-in leg: build the same tiny universe as
# test-sigprofiler-roundtrip.R (NMF orientation: samples as rows) without
# depending on that file's environment.
sp_nmf_matrix_export_helper <- function() {
  t(matrix(
    c(1, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1),
    nrow = 4L, ncol = 3L,
    dimnames = list(
      c("A[C]A", "A[C]C", "A[T]A", "A[T]C"),
      c("Sample1", "Sample2", "Sample3")
    )
  ))
}
