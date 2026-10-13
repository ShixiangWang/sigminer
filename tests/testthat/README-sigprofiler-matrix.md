# SigProfilerExtractor version matrix (local-only)

This note documents the local-run capability provided by
`test-sigprofiler-matrix.R` and `helper-sigprofiler-matrix.R` in this
directory. The matrix is intentionally NOT wired into CI: every test in
`test-sigprofiler-matrix.R` skips cleanly unless the environment variable
`SIGMINER_SP_MATRIX=1` is set, and the container-executing test additionally
requires a `docker` binary on `PATH`. `R CMD check` therefore never pulls an
image.

## What is parameterized

The interop surface (`sigprofiler_extract()` export -> container run ->
`sigprofiler_import()`) is exercised against the ratified
`SigProfilerExtractor` version set:

| Pin | Meaning |
|---|---|
| `1.1.3` | the version sigminer's `sigprofiler_version` argument advertises today |
| `1.2.8` | the `1.2.x` line (latest 1.2 tag, adds the rn7/mm39 builds) |
| `1.4.0` | the latest tagged release line ratified into the matrix (a `v1.5.0` tag exists upstream and can be added once ratified) |

To extend the set (e.g. add `1.5.0` once tagged), edit
`SIGMINER_SP_MATRIX_VERSIONS` in `helper-sigprofiler-matrix.R`.

## How to run

From the package source root, with docker available:

```sh
SIGMINER_SP_MATRIX=1 Rscript -e 'devtools::test(filter = "matrix")'
```

or, without devtools, using the prebuilt-binary route described in the
development setup:

```sh
SIGMINER_SP_MATRIX=1 Rscript -e 'testthat::test_local(filter = "matrix")'
```

## Container coordinate and the registry mirror policy

Each matrix cell runs a one-shot container built from a mirrored python base
image, and pins `SigProfilerExtractor` through pip inside the container:

```
docker run --rm -v <workdir>:/sp <mirror>/library/python:3.11-slim \
  bash -lc "pip install --quiet SigProfilerExtractor==<pin> && python /sp/run_cell.py"
```

- `<mirror>` defaults to the project's Docker Hub mirror host
  `9f6c5in3dh5vszgdz0.xuanyuan.run` (the bench/G0 policy: base pulls spell
  the mirror host explicitly; direct upstream pulls are prohibited on this
  network and a 407 there is environmental evidence, not a bug). Override
  with `SIGMINER_SP_MIRROR_HOST`. The sibling mirrors for other registries
  follow the same family (for `ghcr.io` the host form is
  `9f6c5in3dh5vszgdz0-ghcr.xuanyuan.run`).
- If the default PyPI index is not reachable, point pip at a working index
  with `SIGMINER_SP_PIP_INDEX_URL` (passed through as `pip -i <url>`).
- The extraction run-script inside the container is written by R
  (`sigminer_sp_write_runner()`) and uses only fixed `/sp/...` literal
  paths, so no host path ever re-enters a quoting layer (same fail-closed
  discipline as the `#470` recall-command fix).

## Interpretation

A matrix cell failing is information, not necessarily a bug: older pins may
legitimately emit output trees `sigprofiler_import()` cannot read at the
current branch. The run is opt-in and its only contract is that whatever
fails, fails loudly per version pin.
