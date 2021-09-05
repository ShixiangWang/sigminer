#' An Unified Interface to Extract Signatures
#'
#' This function provides an unified interface to signature extractor
#' implemented in **sigminer**. If you determine a specific `approach`,
#' please also read the documentation of corresponding extractor.
#' See "Arguments" part.
#'
#' @inheritParams sigprofiler
#' @param approach approach name.
#' - "repeated_nmf" - [sig_extract]
#' - "bayes_nmf" - [sig_auto_extract]
#' - "bootstrap_nmf" - [bp_extract_signatures]
#' - "sigprofiler" - [sigprofiler]
#' @param ... other parameters passing to signature extractor based
#' on the `approach` setting.
#'
#' @return Result dependent on the `approach` setting.
#' @export
#' @examples
#' \donttest{
#' load(system.file("extdata", "toy_copynumber_tally_W.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Extract signatures
#' # It is same as sig_extract(cn_tally_W$nmf_matrix, 2, nrun = 1)
#' res <- sig_unify_extract(cn_tally_W$nmf_matrix, 2,
#'   nrun = 1,
#'   approach = "repeated_nmf"
#' )
#' # Auto-extract signatures based on bayesian NMF
#' res2 <- sig_unify_extract(cn_tally_W$nmf_matrix,
#'   nrun = 1,
#'   approach = "bayes_nmf"
#' )
#' }
#' @testexamples
#' expect_s3_class(res, "Signature")
#' expect_s3_class(res2, "Signature")
#' @seealso [sig_extract], [sig_auto_extract], [bp_extract_signatures],
#' [sigprofiler]
sig_unify_extract <- function(nmf_matrix,
                              range = 2:5,
                              nrun = 10,
                              approach = c(
                                "bayes_nmf", "repeated_nmf",
                                "bootstrap_nmf", "sigprofiler"
                              ),
                              cores = 1L,
                              ...) {
  approach <- match.arg(approach)

  if (approach == "repeated_nmf") {
    if (length(range) > 1) {
      stop("For approach 'repeated_nmf', range parameter can only set to length-1 integer for specified signature number")
    }
  }
  if (approach == "bayes_nmf") {
    message("Set ", max(range), " as the initial signature number")
  }


  if (approach == "repeated_nmf") {
    sig_extract(
      nmf_matrix = nmf_matrix,
      n_sig = range,
      nrun = nrun,
      cores = cores,
      ...
    )
  } else if (approach == "bayes_nmf") {
    sig_auto_extract(
      nmf_matrix = nmf_matrix,
      K0 = max(range),
      nrun = nrun,
      cores = cores,
      ...
    )
  } else if (approach == "bootstrap_nmf") {
    bp_extract_signatures(
      nmf_matrix = nmf_matrix,
      range = range,
      n_nmf_run = nrun,
      cores = cores,
      ...
    )
  } else {
    sigprofiler_extract(
      nmf_matrix = nmf_matrix,
      range = range,
      nrun = nrun,
      cores = cores,
      ...
    )
  }
}
