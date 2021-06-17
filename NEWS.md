# sigminer 2.0.2

- Fixed `output_sig()` error in handling exposure plot with >9 signatures (#366).
- Added `limitsize = FALSE` for `ggsave()` or `ggsave2()` for handling big figure.

# sigminer 2.0.1

- Supported `mm9` genome build.
- Removed FTP link as CRAN suggested (#359).
- Updated README.

# sigminer 2.0.0

## BUG REPORTS

- Fixed the SigProfiler installation error due to Python version in conda environment.
- Fixed classification bug due to repeated function name `call_component`.
- Fixed the bug when `read_vcf()` with `##` commented VCF files.

## ENHANCEMENTS

- Added support for latest COSMIC v3.2 as reference signatures. You can obtain them by

```r
for (i in c("latest_SBS_GRCh37", "latest_DBS_GRCh37", "latest_ID_GRCh37",
            "latest_SBS_GRCh38", "latest_DBS_GRCh38",
            "latest_SBS_mm9", "latest_DBS_mm9",
            "latest_SBS_mm10", "latest_DBS_mm10",
            "latest_SBS_rn6", "latest_DBS_rn6")) {
  message(i)
  get_sig_db(i)
}
```

- Updated `keep_only_pass` to `FALSE` at default.
- Added RSS and unexplained variance calculation in `get_sig_rec_similarity()`.
- Added data check and filter in `output_tally()` and `show_catalogue()`.
- Enhanced `show_group_enrichment()` (#353) & added a new option to cluster rows.
- Removed unnecessary CN classifications code in recent development.

## NEW FUNCTIONS

## DEPRECATED

- Dropped copy number "M"" method to avoid misguiding user to use/read wrong signature profile and keep code simple.

# sigminer 1.2.5

## BUG REPORTS

## ENHANCEMENTS

- Modified the default visualization of `bp_show_survey()`.
- Enhanced `torch` check.

## NEW FUNCTIONS

- `read_sv_as_rs()` and `sig_tally.RS()` for simplified genome rearrangement classification matrix generation (experimental).

## DEPRECATED

# sigminer 1.2.4

## BUG REPORTS

- Fixed the assign problem about match pair in `bp_extract_signatures()` 
with `lpSolve` package instead of using my problematic code.

## ENHANCEMENTS

- Supported `mm10` in `read_vcf()`.
- Removed large data files and store them in Zenodo to reduce package size.
- Added cores check.
- Upgraded SP to v1.1.0 (need test).
- Tried installing Torch before SP (need test).

## NEW FUNCTIONS

## DEPRECATED

# sigminer 1.2.2

## BUG REPORTS

- Fixed bug in silhouette calculation in `bp_extract_signatures()` (#332).
PAY ATTENTION: this may affect results.
- Fixed bug using custom signature name in `show_sig_profile_loop()`.

## ENHANCEMENTS

- Subset signatures to plot is available by `sig_names` option.
- sigminer is available in bioconda channel: <https://anaconda.org/bioconda/r-sigminer/>
- Updated `ms` strategy in `sig_auto_extract()` by assigning each signature to its
best matched reference signatures.
- Added `get_shannon_diversity_index()` to get diversity index for signatures (#333).
- Added new method "S" (from Steele et al. 2019) for tallying copy number data (#329).
- Included new (RS) reference signatures (related to #331).
- Updated the internal code for getting relative activity in `get_sig_exposure()`.

## NEW FUNCTIONS

- `bp_get_clustered_sigs()` to get clustered mean signatures.

## DEPRECATED

# sigminer 1.2.1

- Updated author list.

## BUG REPORTS

## ENHANCEMENTS

- Added a quick start vignette.
- A new option `highlight` is added to `show_sig_number_survey()` and `bp_show_survey2()` to highlight a selected number.

## NEW FUNCTIONS

## DEPRECATED

# sigminer 1.2.0

## BUG REPORTS

## ENHANCEMENTS

- A new option `cut_p_value` is added to `show_group_enrichment()` to cut continous p values as binned regions.
- A Python backend for `sig_extract()` is provided.
- User now can directly use `sig_extract()` and `sig_auto_extract()` instead of loading NMF package firstly.
- Added benchmark results for different extraction approaches in README.
- The threshold for `auto_reduce` in `sig_fit()` is modified from 0.99 to 0.95 and similarity update threshold updated from `>0` to `>=0.01`.
- Removed `pConstant` option from `sig_extract()` and `sig_estimate()`. Now a
auto-check function is created for avoiding the error from NMF package due to
no contribution of a component in all samples.

## NEW FUNCTIONS

- `bp_show_survey2()` to plot a simplified version for signature number survey (#330).
- `read_xena_variants()` to read variant data from UCSC Xena as a `MAF` object for signature analysis.
- `get_sig_rec_similarity()` for getting reconstructed profile similarity for `Signature` object (#293).
- Added functions start with `bp_` which are combined to provide a best practice for extracting
signatures in cancer researches. See more details, run `?bp` in your R console.

## DEPRECATED

# sigminer 1.1.0

- Added data simulation.
- Suppressed `future` warnings.
- Fixed p value calculation in bootstrap analysis.
- Fixed typo in `show_cor()`, thanks to @Miachol.
- Added `y_tr` option in `show_sig_profile()` to transform y axis values.
- Optimized default behavior of `read_copynumber()`.
  - Support LOH records when user input minor allele copy number.
  - Set `complement = FALSE` as default.
  - Free dependencies between option `use_all` and `complement`.
- Added visualization support for genome rearrangement signatures (#300).
- Added four database for reference signatures from <https://doi.org/10.1038/s43018-020-0027-5> (#299).
- Added new measure 'CV' for `show_sig_bootstrap()` (#298).
- Added `group_enrichment()` and `show_group_enrichment()` (#277).
- Optimized signature profile visualization (#295).
- Updated `?sigminer` documentation.
- Added `ms` strategy to select optimal solution by maximizing cosine similarity
to reference signatures.
- Added `same_size_clustering()` for same size clustering.
- Added `show_cosmic()` to support reading COSMIC signatures in web browser (#288).
- Changed argument `rel_threshold` behavior in `sig_fit()` and `get_sig_exposure()`.
Made them more consistent and allowed un-assigned signature contribution (#285).
- Updated all COSMIC signatures to v3.1 and their aetiologies (#287).

# sigminer 1.0.19

- Added more specific reference signatures from SigProfiler, e.g. `SBS_mm9`.
- Supported `data.frame` as input object for `sig` in `get_sig_similarity()` and `sig_fit()`.
- Modified `g_label` option in `show_group_distribution()` to better control group names. 
- Added `test` option and variable checking in `show_cor()`.
- Updated `output_sig()` to output signature exposure distribution ([#280](https://github.com/ShixiangWang/sigminer/issues/280)).
- Added `show_cor()` for general association analysis.
- Added options in `show_group_distribution()` to control segments.

# sigminer 1.0.18

- Fixed bugs when outputing only 1 signatures.
- Fixed label inverse bug in `add_labels()`, thanks to TaoTao for reporting.

# sigminer 1.0.17

- Handled `,` seperated indices in show_cosmic_signatures.
- Added option `set_order` in `get_sig_similarity()` (#274).
- Outputed more stats information in `output_sig()`.
- Fixed default y axis title in `show_sig_bootstrap_error()`, now it is "Reconstruction error (L2 norm)"

# sigminer 1.0.16

- Added `auto_reduce` option in `sig_fit*` functions to improve signature fitting.
- Return cosine similarity for sample profile in `sig_fit()`.
- Set default strategy in `sig_auto_extract()` to 'optimal'.
- Supported search reference signature index in `get_sig_cancer_type_index()`.
- Outputed legacy COSMIC similarity for SBS signatures.
- Added new option in `sigprofiler_extract()` to reduce failure in when `refit` is enabled.
- Outputed both relative and absolute signature exposure in `output_sig()`.
- Updated background color in `show_group_distribution()`.
- Modified the default theme for signature profile in COSMIC style.
- Updated the copy number classification method.

# sigminer 1.0.15

- Handled null catalogue.
- Supported ordering the signatures for results from SigProfiler.
- Supported importing refit results from SigProfiler.
- Set `optimize` option in `sig_extract()` and `sig_auto_extract()`.

# sigminer 1.0.14

- Supported signature index separated by `,` in `sig_fit()` and `sig_fit_bootstrap*` functions.
- Added `output_*` functions from [sigflow](https://github.com/ShixiangWang/sigflow).
- Enhanced DBS search and error handling in `sig_tally()`.
- Added option `highlight_genes` in `show_cn_group_profile()` to show gene labels.
- Added `get_sig_cancer_type_index()` to get reference signature index.
- Added `show_group_distribution()` to show group distribution.
- Added options in `show_cn_profile()` to show specified ranges and add copy number value labels.
- Used package `nnls` instead of `pracma` for NNLS implementation in `sig_fit()`.

# sigminer 1.0.13

- Supported `BSgenome.Hsapiens.1000genomes.hs37d5` in `sig_tally()`.
- Remove changing `MT` to `M` in mutation data.
- Fixed bug in extract numeric signature names and signature orderings in `show_sig_exposure()`.
- Added `letter_colors` as an unexported discrete palette.

# sigminer 1.0.12

- Added `transform_seg_table()`.
- Added `show_cn_group_profile()`.
- Added `show_cn_freq_circos()`.
- `sig_orders` option in `show_sig_profile()` function now can select and order signatures to plot.
- Added `show_sig_profile_loop()` for better signature profile visualization.

# sigminer 1.0.11

- Added option to control the SigProfilerExtractor to avoid issue in docker image build.

# sigminer 1.0.10

- Some updates.
- Compatible with SigProfiler 1.0.15

# sigminer 1.0.9

- Tried to speed up joining adjacent segments in `read_copynumber()`, got 200% improvement.

# sigminer 1.0.8

- Tried to speed up joining adjacent segments in `read_copynumber()`, got 20% improvement.
- Added `cosine()` function.
- Added and exported `get_sig_db()` to let users directly load signature database.
- Added `sigprofiler_extract()` and `sigprofiler_import()` to call SigProfiler and import results.
- Added `read_vcf()` for simply reading VCF files.
- Implemented DBS-1248.
- Added `show_sig_profile_heatmap()`.
- Supported mouse genome 'mm10' ([#241](https://github.com/ShixiangWang/sigminer/issues/241)).
- Added `read_copynumber_seqz()` to read sequenza result directory.
- Speed up the annotation process in `read_copynumber()`.

# sigminer 1.0.7

- Fixed bug in OsCN feature calculation.
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
