# Local-run-only version-matrix harness for the SigProfiler interop bridge.
#
# Parameterizes the interop surface over the ratified SigProfilerExtractor
# versions {1.1.3, 1.2.x -> 1.2.8, 1.4.0}. Executed ONLY when the
# environment variable SIGMINER_SP_MATRIX=1 is set; this is deliberately not
# wired into CI. See README-sigprofiler-matrix.md in this directory.
#
# The container coordinate follows the bench pattern: a mirrored python base
# image plus a pip pin of SigProfilerExtractor, with the run script written
# by R so the container never sees shell-embedded user paths (fail-closed:
# all in-container paths are fixed literals).

SIGMINER_SP_MATRIX_VERSIONS <- c("1.1.3", "1.2.8", "1.4.0")
SIGMINER_SP_BASE_IMAGE <- "python:3.11-slim"
SIGMINER_SP_DEFAULT_MIRROR <- "9f6c5in3dh5vszgdz0.xuanyuan.run"

sigminer_sp_mirror_host <- function(default = SIGMINER_SP_DEFAULT_MIRROR) {
  host <- Sys.getenv("SIGMINER_SP_MIRROR_HOST", "")
  if (nzchar(host)) host else default
}

# Docker Hub official images are namespaced under library/ behind the
# project's registry mirrors (bench/m7 policy: the mirror host spells the
# upstream registry, e.g. docker.io -> 9f6c5in3dh5vszgdz0.xuanyuan.run,
# ghcr.io -> 9f6c5in3dh5vszgdz0-ghcr.xuanyuan.run).
sigminer_sp_image_ref <- function(base = SIGMINER_SP_BASE_IMAGE,
                                  mirror = sigminer_sp_mirror_host()) {
  sigminer_sp_validate_token(base, "base image")
  sigminer_sp_validate_token(mirror, "mirror host")
  paste0(mirror, "/library/", base)
}

# Single shell token: no whitespace, quotes, or expansion characters.
# Anything else in a bare command position invites word splitting,
# globbing, or variable expansion -- fail closed instead.
sigminer_sp_validate_token <- function(value, what) {
  if (!is.character(value) || length(value) != 1L || anyNA(value) ||
      !nzchar(value)) {
    stop(sprintf("`%s` must be a single non-empty string", what),
         call. = FALSE)
  }
  bad <- c(" ", "\t", "\n", "\r", "'", "\"", "$", "`",
           "*", "?", "[", "]", ";", "|", "&", "<", ">", "\\", "%")
  chars <- strsplit(value, "")[[1L]]
  hit <- unique(chars[!is.na(match(chars, bad))])
  if (length(hit) > 0L) {
    stop(sprintf(
      paste0("`%s` contains shell-hostile character(s): %s. ",
             "Use a path without whitespace, quotes, or expansion syntax."),
      what, paste0("'", hit, collapse = "', '")),
      call. = FALSE)
  }
  invisible(value)
}

# Pin a released SigProfilerExtractor version: plain semver only.
sigminer_sp_validate_version <- function(version) {
  sigminer_sp_validate_token(version, "version")
  if (!grepl("^[0-9]+(\\.[0-9]+)*$", version)) {
    stop(sprintf("`version` must be plain semver like 1.2.8, got '%s'", version),
         call. = FALSE)
  }
  invisible(version)
}

# Write the in-container run script next to the exported matrix input file.
# All paths inside are fixed /sp/... literals baked into the mount, so the
# host path never re-enters a quoting layer.
sigminer_sp_write_runner <- function(workdir) {
  sigminer_sp_validate_token(workdir, "workdir")
  runner <- c(
    "from SigProfilerExtractor import sigpro as sig",
    "sig.sigProfilerExtractor(",
    "    'matrix',",
    "    '/sp/out',",
    "    '/sp/sigprofiler_input.txt',",
    "    reference_genome='GRCh37',",
    "    opportunity_genome='GRCh37',",
    "    minimum_signatures=2,",
    "    maximum_signatures=5,",
    "    nmf_replicates=10,",
    "    exome=False,",
    "    nmf_init='random',",
    "    make_decomposition_plots=False,",
    "    cpu=-1,",
    ")"
  )
  path <- file.path(workdir, "run_cell.py")
  writeLines(runner, con = path)
  path
}

sigminer_sp_matrix_command <- function(version, workdir,
                                       image = sigminer_sp_image_ref()) {
  sigminer_sp_validate_version(version)
  sigminer_sp_validate_token(workdir, "workdir")
  sigminer_sp_validate_token(image, "image")
  pip_index <- Sys.getenv("SIGMINER_SP_PIP_INDEX_URL", "")
  index_arg <- if (nzchar(pip_index)) {
    sigminer_sp_validate_token(pip_index, "pip index url")
    paste0("-i ", pip_index, " ")
  } else {
    ""
  }
  paste0(
    "docker run --rm -v ", workdir, ":/sp ", image,
    " bash -lc \"",
    "pip ", index_arg, "install --quiet SigProfilerExtractor==", version,
    " && python /sp/run_cell.py",
    "\""
  )
}
