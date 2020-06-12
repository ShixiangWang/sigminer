# sigminer 1.0.7

- Removed useless options in `read_maf()`.
- Modify method 'LS' in `sig_fit()` to 'NNLS' and implement it with **pracma** package ([#216](https://github.com/ShixiangWang/sigminer/issues/216)).
- Made `use_all` option in `read_copynumber()` working correctly.
- Fixed potential problem raised by unordered copy number segments ([#217](https://github.com/ShixiangWang/sigminer/issues/217)).
- Fixed a typo, correct `MRSE` to `RMSE`.
- Added feature in `show_sig_bootstrap_*()` for plotting aggregated values.
- Fixed bug when use `get_groups()` for clustering.
- Fixed bug about using reference components from NatGen 2018 paper.
- Added option `highlight_size` for `show_sig_bootstrap_*()`.
- Fixed bug about signature profile plotting for method 'M'.

# sigminer 1.0.6

- Added "scatter" in `sig_fit()` function to better visualize a few samples.
- Added "highlight" option.
- `lsei` package was removed from CRAN, here I reset default method to 'QP' and tried best to keep the LS usage in sigminer ([#189](https://github.com/ShixiangWang/sigminer/issues/189)).
- Made consistent copy number labels in `show_sig_profile()` and added input checking for this function.
- Fixed unconsistent bootstrap when use `furrr`, solution is from <https://github.com/DavisVaughan/furrr/issues/107>.
- Properly handled null-count sample in `sig_fit()` for methods `QP` and `SA`.
- Supported boxplot or violin in `show_sig_fit()` and `show_sig_bootstrap_*` functions.
- Added job mode for `sig_fit_bootstrap_batch` for more useful in practice.
- Added `show_groups()` to show the signature contribution in each group from `get_groups()`.
- Expanded clustering in `get_groups()` to result of `sig_fit()`.
- Properly handled null-count samples in `sig_fit_bootstrap_batch()`.
- Added strand bias labeling for INDEL.
- Added COSMIC TSB signatures.

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
