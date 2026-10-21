## Resubmission

This is a resubmission of sigminer 2.4.0. The previous submission
(2026-10-20) was auto-rejected with 2 NOTEs on Debian / 1 NOTE on Windows.
Both are fixed in this version:

1. **NEWS.md parse NOTE** (`checking package subdirectories`):
   the leading `sigminer (unreleased)` section title carried no version
   number, so CRAN's NEWS parser could not extract version info from it.
   The empty section was removed; `NEWS.md` now opens with
   `sigminer 2.4.0`, the same format as all previously accepted releases.

2. **Example runtime NOTE** (`bp_select_k` example took ~8s, above the
   5s threshold): the example now passes an explicit reduced survey
   budget (`replicates = 2L, k_folds = 4L, n_seeds = 1L`) alongside the
   unchanged planted catalogue and pinned seed. The example runs in
   ~1.5s locally (well under the 5s limit even on slower machines), and
   the in-example assertions still verify the seeded K = 2 verdict and
   the single selected evidence row. The production defaults of
   `bp_select_k()` are unchanged.

No other code, API, or documentation changes relative to the previous
submission, apart from re-running `roxygen2` so three Rd files that had
gone stale relative to their R sources (`show_sig_bootstrap.Rd`,
`sigprofiler.Rd`, and the roxytest line reference in
`test-roxytest-testexamples-sig_fit.R`) match the code again.

## R CMD check results

Local `R CMD check --as-cran` on sigminer_2.4.0.tar.gz (R 4.5.2, macOS
aarch64) is clean for everything CRAN's pretest actually flagged:

* `checking package subdirectories ... OK` (the NEWS.md version-parse
  NOTE is gone).
* No example exceeds 5s (`bp_select_k` now runs in ~1.5s); the
  "Examples with CPU ... > 5s" NOTE is gone.

The remaining local findings are environment artifacts of the check
machine, not of the package (all were OK in CRAN's own pretest of the
previous tarball):

* `examples` ERROR: locally installed **parallelly** 1.38.0 predates
  the `max` argument of `availableCores()` (added in 1.41.0). Current
  CRAN parallelly is unaffected; CRAN's pretest ran the same examples
  fine.
* `PDF version of manual` WARNING/ERROR: no `pdflatex` on the local
  machine. Rd-to-LaTeX conversion itself succeeds.
* Other NOTEs: standard incoming-feasibility maintainer note, installed
  size 5.6Mb (macOS binaries), unverifiable local clock, R 4.5 HTML
  help validation cosmetics, and the leftover `sigminer-manual.tex`
  from the failed local PDF step.

## Test environments

* Local: R 4.5.2, macOS (aarch64).
* GitHub Actions (r-lib/actions): R release on ubuntu-latest, macOS and
  Windows; R devel on ubuntu-latest.

## revdepcheck results

There are currently no reverse dependencies for this package on CRAN.
