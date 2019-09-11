#' Extract Signatures through the Automatic Relevance Determination Technique
#'
#' A bayesian variant of NMF algorithm to enable optimal inferences for the
#' number of signatures through the automatic relevance determination technique.
#' This functions delevers highly interpretable and sparse representations for
#' both signature profiles and attributions at a balance between data fitting and
#' model complexity. See detail part and references for more.
#'
#' There are three methods available in this function: "L1W.L2H", "L1KL" and "L2KL".
#' They use different priors for the bayesian variant of NMF algorithm
#' (see `method` parameter) written by reference #1 and implemented in
#' [SignatureAnalyzer software](https://software.broadinstitute.org/cancer/cga/msp)
#' (reference #2).
#'
#' I copied source code for the three methods from Broad Institute and supplementary
#' files of reference #3, and wrote this higher function. It is more friendly for users
#' to extract, visualize and analyze signatures by combining with other powerful functions
#' in **sigminer** package. Besides, I implemented parallel computation to speed up
#' the calculation process and a similar input and output structure like [sig_extract()].
#'
#'
#' @inheritParams sig_estimate
#' @inheritParams derive
#' @param result_prefix prefix for result data files.
#' @param destdir path to save data runs, default is `tempdir()`.
#' @param method default is "L1W.L2H", which uses an exponential prior for W and
#' a half-normal prior for H (This method is used by PCAWG project, see reference #3).
#' You can also use "L1KL" to set expoential priors for both W and H, and "L2KL" to
#' set half-normal priors for both W and H. The latter two methods are originally
#' implemented by [SignatureAnalyzer software](https://software.broadinstitute.org/cancer/cga/msp).
#' @param K0 number of initial signatures.
#' @param nrun number of independent simulations.
#' @param tol tolerance for convergence.
#' @param recover if `TRUE`, try to recover result from previous runs based on input `result_prefix`,
#' `destdir` and `nrun`. This is pretty useful for reproducing result.
#' @author Shixiang Wang
#' @references
#' Tan, Vincent YF, and Cédric Févotte. "Automatic relevance determination in nonnegative matrix factorization with the/spl beta/-divergence."
#'  IEEE Transactions on Pattern Analysis and Machine Intelligence 35.7 (2012): 1592-1605.
#'
#' Kim, Jaegil, et al. "Somatic ERCC2 mutations are associated with a distinct genomic signature in urothelial tumors."
#'  Nature genetics 48.6 (2016): 600.
#'
#' Alexandrov, Ludmil, et al. "The repertoire of mutational signatures in human cancer." BioRxiv (2018): 322859.
#'
#' @return a `list` with `Signature` class.
#' @export
#'
#' @examples
#' # Load copy number prepare object
#' load(system.file("extdata", "toy_copynumber_prepare.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' res <- sig_auto_extract(cn_prepare$nmf_matrix, result_prefix = "Test_copynumber", nrun = 1)
#' # At default, all run files are stored in tempdir()
#' dir(tempdir())
#' @seealso [sig_extract] for extracting signatures using **NMF** package, [sig_estimate] for
#' estimating signature number for [sig_extract].
sig_auto_extract <- function(nmf_matrix = NULL,
                             result_prefix = "BayesNMF",
                             destdir = tempdir(),
                             method = c("L1W.L2H", "L1KL", "L2KL"),
                             K0 = 25,
                             nrun = 10,
                             niter = 2e5,
                             tol = 1e-07,
                             cores = 1,
                             recover = FALSE) {
  on.exit(invisible(gc())) # clean when exit
  method <- match.arg(method)

  filelist <- file.path(destdir, paste(result_prefix, seq_len(nrun), "rds", sep = "."))

  if (recover) {
    message("Recover mode is on, check if all files exist...")
    all_exist = all(file.exists(filelist))
    if (!all_exist) {
      stop("Recover failed, cannot find previous result files, please run with 'recover = FALSE'.")
    }
    message("Yup! Recovering...")
  }

  if (!recover) {
    nmf_matrix <- t(nmf_matrix) # rows for mutation types and columns for samples

    future::plan("multiprocess", workers = cores)
    furrr::future_map(seq_len(nrun), function(i, method, filelist) {
      if (method == "L1W.L2H") {
        # Apply BayesNMF - L1W.L2H for an exponential prior for W and a half-normal prior for H
        res <- BayesNMF.L1W.L2H(nmf_matrix, niter, 10, 5, tol, K0, K0, 1)
      } else if (method == "L1KL") {
        # Apply BayesNMF - L1KL for expoential priors
        res <- BayesNMF.L1KL(nmf_matrix, niter, 10, tol, K0, K0, 1)
      } else {
        # Apply BayesNMF - L2KL for half-normal priors
        res <- BayesNMF.L2KL(nmf_matrix, niter, 10, tol, K0, K0, 1)
      }
      saveRDS(res, file = filelist[i])
    }, method = method, filelist = filelist, .progress = TRUE)
  }

  summary.run <- purrr::map_df(seq_len(nrun), function(i) {
    res <- readRDS(filelist[i])
    K <- sum(colSums(res[[1]]) > 1)
    posterior <- -res[[4]] # res[[4]] = -log(posterior)
    dplyr::tibble(
      Run = i,
      K = K,
      posterior = posterior,
      file = filelist[i]
    )
  })

  summary.run <- summary.run %>%
    dplyr::arrange(dplyr::desc(.data$posterior))

  # select best solution
  best <- names(sort(table(summary.run$K), decreasing = TRUE))[1] %>%
    as.integer()

  best_row <- dplyr::filter(summary.run, .data$K == best) %>%
    head(1)

  message("Select Run ", best_row$Run, ", which K = ", best_row$K, " as best solution.")
  best_solution <- get_bayesian_result(best_row)

  res <- list(
    Signature = best_solution$Signature,
    Signature.norm = best_solution$Signature.norm,
    Exposure = best_solution$Exposure,
    Exposure.norm = best_solution$Exposure.norm,
    K = best_solution$K,
    Raw = list(
      summary_run = summary.run,
      W = best_solution$W,
      H = best_solution$H,
      best_run = best_row$Run
    )
  )
  class(res) <- "Signature"
  attr(res, "nrun") <- nrun
  attr(res, "method") <- method
  # attr(res, "seed") <- seed
  attr(res, "call_method") <- "BayesianNMF"

  res
}