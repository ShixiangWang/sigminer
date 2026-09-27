## R CMD check results

0 errors | 0 warnings | 1 note

* The note is the standard CRAN incoming feasibility note ("Maintainer: ...")
  produced when checking from a local machine.
* This is a maintenance (bug fix) release:
  - `sig_estimate()`, `sig_extract()` and `bp_extract_signatures()` now honor
    core numbers > 2 by registering a foreach backend, because the built-in
    parallel backend of the NMF package caps requested workers at 2 (#479).
  - Fixed NULL row names of the NMF matrix returned by `sig_tally()` with the
    Wang method for single-sample inputs (#465).
  - Dropped the deprecated `refit_denovo_signatures` argument passed to
    SigProfilerExtractor for compatibility with its newer versions; the
    `refit` argument of `sigprofiler_extract()` is deprecated and ignored
    with a warning (#470).
  - Fixed "Unknown or uninitialised column" warnings of `read_maf_minimal()`
    when the input is a tibble (#461).

## Test environments

* Local: R 4.5.2, macOS (aarch64) -- 0 errors, 0 warnings
* GitHub Actions (r-lib/actions): R release on ubuntu-latest, macOS and
  Windows; R devel on ubuntu-latest -- all pass

## revdepcheck results

There are currently no reverse dependencies for this package on CRAN.
