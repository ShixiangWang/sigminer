# ========================================================
# test analysis functions, use test_ as start
# ========================================================


#' @title Test running status for getting components
#' @description This is used when user find it is hard to determine component number
#' by [get_components] function.
#'
#'
#' @inheritParams get_components
#' @param min_comp minimal number of components to fit, default is 2.
#' @param max_comp maximal number of components to fit, default is 10.
#' @param feature_name feature name to test.
#' @author Shixiang Wang
#' @return a `list` contain `flexmix` object of copy-number features.
#' @importClassesFrom flexmix FLXcontrol
#' @export
#' @examples
#' # Load copy number features
#' load(system.file("extdata", "toy_cn_features.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' 
#' # Test component number from 2 to 3 for feature 'sigment size'
#' test_run_components(cn_features, feature_name = "segsize", max_comp = 3)
#' @family test function series
test_run_components <- function(CN_features,
                                feature_name = c("segsize", "bp10MB", "osCN", "bpchrarm", "changepoint", "copynumber"),
                                seed = 123456,
                                min_comp = 2,
                                max_comp = 10,
                                min_prior = 0.001,
                                model_selection = "BIC",
                                nrep = 1,
                                niter = 1000) {
  feature_name <- match.arg(feature_name)

  if (feature_name == "segsize") {
    dat <- as.numeric(CN_features[["segsize"]][, 2])
    message("Fit feature: Segment size")
    segsize_mm <-
      fitComponent(
        dat,
        seed = seed,
        model_selection = model_selection,
        min_prior = min_prior,
        niter = niter,
        nrep = nrep,
        min_comp = min_comp,
        max_comp = max_comp
      )
  } else if (feature_name == "bp10MB") {
    dat <- as.numeric(CN_features[["bp10MB"]][, 2])
    message("Fit feature: Breakpoint count per 10 Mb")
    bp10MB_mm <-
      fitComponent(
        dat,
        dist = "pois",
        seed = seed,
        model_selection = model_selection,
        min_prior = min_prior,
        niter = niter,
        nrep = nrep,
        min_comp = min_comp,
        max_comp = max_comp
      )
  } else if (feature_name == "osCN") {
    dat <- as.numeric(CN_features[["osCN"]][, 2])
    message("Fit feature: Length of oscillating copy-number chain")
    osCN_mm <-
      fitComponent(
        dat,
        dist = "pois",
        seed = seed,
        model_selection = model_selection,
        min_prior = min_prior,
        niter = niter,
        nrep = nrep,
        min_comp = min_comp,
        max_comp = max_comp
      )
  } else if (feature_name == "bpchrarm") {
    dat <- as.numeric(CN_features[["bpchrarm"]][, 2])
    message("Fit feature: Breakpoint count per arm")
    bpchrarm_mm <-
      fitComponent(
        dat,
        dist = "pois",
        seed = seed,
        model_selection = model_selection,
        min_prior = min_prior,
        niter = niter,
        nrep = nrep,
        min_comp = min_comp,
        max_comp = max_comp
      )
  } else if (feature_name == "changepoint") {
    dat <- as.numeric(CN_features[["changepoint"]][, 2])
    message("Fit feature: Copy number change")
    changepoint_mm <-
      fitComponent(
        dat,
        seed = seed,
        model_selection = model_selection,
        min_prior = min_prior,
        niter = niter,
        nrep = nrep,
        min_comp = min_comp,
        max_comp = max_comp
      )
  } else if (feature_name == "copynumber") {
    dat <- as.numeric(CN_features[["copynumber"]][, 2])
    message("Fit feature: Absolute copy number")
    copynumber_mm <-
      fitComponent(
        dat,
        seed = seed,
        model_selection = model_selection,
        nrep = nrep,
        min_comp = min_comp,
        max_comp = max_comp,
        min_prior = min_prior,
        niter = niter
      )
  } else {
    stop("Bad feature_name!")
  }

  invisible(NULL)
}
