# sigminer 1.0.5

- Exported APOBEC result when the mode is 'ALL' in `sig_tally()`.
- Added batch bootstrap analysis feature (#158).
- Supported all common signature plotting.
- Added strand feature to signature profile.

# sigminer 1.0.4

- Added profile plot for DBS and INDEL.
- Fixed error for signature extraction in mode 'DBS' or 'ID'.
- Fixed method 'M' for CN tally cannot work when `cores > 1` (#161).

# sigminer 1.0.3

- Added multiple methods for `sig_fit()`.
- Added feature `sig_fit_bootstrap()` for bootstrap results.
- Added multiple classification method for SBS signature.
- Added strand bias enrichment analysis for SBS signature.
- Moved multiple packages from field `Imports` to `Suggests`.
- Added feature `report_bootstrap_p_value()` to report p values.
- Added common DBS and ID signature.
- Updated citation.

# sigminer 1.0.2

- Added merged transcript info for hg19 and hg38 build, this is availabe by `data()`.
- Added gene info for hg19 and hg38 build to extdata directory.

# sigminer 1.0.1

- Removed `fuzzyjoin` package from dependency.
- Moved `ggalluvial` package to field `suggsets`.

# sigminer 1.0.0

All users, this is a break-through version of **sigminer**,
most of functions have been modified, more features are implemented.
Please read the reference list to see the function groups and their
functionalities.

Please read the vignette for usage.

I Hope it helps your research work and makes a new contribution
to the scientific community.
