#' Estimate Signature Number
#'
#' Use **NMF** package to evaluate the optimal number of signatures.
#' This is used along with [sig_extract].
#' Users should `library(NMF)` firstly.
#'
#' The most common approach is to choose the smallest rank for which cophenetic correlation coefficient
#' starts decreasing (Used by this function). Another approach is to choose the rank for which the plot
#' of the residual sum of squares (RSS) between the input matrix and its estimate shows an inflection point.
#' More custom features please directly use [NMF::nmfEstimateRank].
#'
#' @param nmf_matrix a `matrix` used for NMF decomposition with rows indicate samples and columns indicate components.
#' @param range a `numeric` vector containing the ranks of factorization to try. Note that duplicates are removed and values are sorted in increasing order. The results are notably returned in this order.
#' @param nrun a `numeric` giving the number of run to perform for each value in `range`, `nrun` set to 30~50 is enough to achieve robust result.
#' @param what a character vector whose elements partially match one of the following item, which correspond to the measures computed by summary on each – multi-run – NMF result: ‘all’, ‘cophenetic’, ‘rss’, ‘residuals’, ‘dispersion’, ‘evar’, ‘silhouette’ (and more specific \*.coef, \*.basis, \*.consensus), ‘sparseness’ (and more specific \*.coef, \*.basis). It specifies which measure must be plotted (what='all' plots all the measures).
#' @param cores number of cpu cores to run NMF.
#' @param seed specification of the starting point or seeding method, which will compute a starting point, usually using data from the target matrix in order to provide a good guess.
#' @param use_random Should generate random data from input to test measurements. Default is `TRUE`.
#' @param save_plots if `TRUE`, save plots to local machine. Of note, if pdf file consensus map has extra blank page, please call
#' `pdf(..., onefile=FALSE)`, it may fix this problem, see
#' <https://stackoverflow.com/questions/12481267/in-r-how-to-prevent-blank-page-in-pdf-when-using-gridbase-to-embed-subplot-insi>.
#' @param plot_basename when save plots, set custom basename for file path.
#' @param method specification of the NMF algorithm. Use 'brunet' as default.
#' Available methods for nmf decompositions are 'brunet', 'lee', 'ls-nmf', 'nsNMF', 'offset'.
#' @param pConstant A small positive value to add to the matrix. Use it ONLY if the functions throws an \code{non-conformable arrays} error.
#' @param verbose if `TRUE`, print extra message.
#' @author Shixiang Wang
#' @references Gaujoux, Renaud, and Cathal Seoighe. "A flexible R package for nonnegative matrix factorization." BMC bioinformatics 11.1 (2010): 367.
#' @importFrom grDevices pdf dev.off
#' @return a `list` contains information of NMF run and rank survey.
#' @import NMF
#' @export
#' @examples
#' \donttest{
#' # Load copy number prepare object
#' load(system.file("extdata", "toy_copynumber_prepare.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' library(NMF)
#' cn_estimate <- sig_estimate(cn_prepare$nmf_matrix,
#'   cores = 1, nrun = 5,
#'   verbose = TRUE
#' )
#' }
#' @seealso [sig_extract] for extracting signatures using **NMF** package, [sig_auto_extract] for
#' extracting signatures using automatic relevance determination technique.
sig_estimate <-
  function(nmf_matrix,
             range = 2:5,
             nrun = 10,
             what = "all",
             cores = 1,
             seed = 123456,
             use_random = TRUE,
             save_plots = FALSE,
             plot_basename = file.path(tempdir(), "nmf"),
             method = "brunet",
             pConstant = NULL,
             verbose = FALSE) {
    mat <- t(nmf_matrix)

    # To avoid error due to non-conformable arrays
    if (!is.null(pConstant)) {
      if (pConstant < 0 | pConstant == 0) {
        stop("pConstant must be > 0")
      }
      mat <- mat + pConstant
    }

    if (cores > 1) {
      estim.r <-
        NMF::nmfEstimateRank(
          mat,
          range,
          method = method,
          nrun = nrun,
          verbose = verbose,
          seed = seed,
          .opt = paste0("P", cores)
        )
    } else {
      estim.r <-
        NMF::nmfEstimateRank(
          mat,
          range,
          method = method,
          nrun = nrun,
          verbose = verbose,
          seed = seed
        )
    }

    nmf.sum <- NMF::summary(estim.r) # Get summary of estimates
    if (verbose) {
      message("Estimation of rank based on observed data.")
      print(nmf.sum)
    }

    if (use_random) {
      if (verbose) message("Generating random matrix and run NMF...")
      V.random <- NMF::randomize(mat)

      if (cores > 1) {
        estim.r.random <-
          NMF::nmfEstimateRank(
            V.random,
            range,
            method = method,
            nrun = nrun,
            verbose = verbose,
            seed = seed,
            .opt = paste0("P", cores)
          )
      } else {
        estim.r.random <-
          NMF::nmfEstimateRank(
            V.random,
            range,
            method = method,
            nrun = nrun,
            verbose = verbose,
            seed = seed
          )
      }

      nmf.sum.random <- NMF::summary(estim.r.random) # Get summary of estimates
      if (verbose) {
        message("Estimation of rank based on random data.")
        print(nmf.sum.random)
      }
    } else {
      estim.r.random <- NULL
      nmf.sum.random <- NULL
    }

    if (use_random) {
      p <- NMF::plot(
        estim.r,
        estim.r.random,
        what = what,
        xname = "Observed",
        yname = "Randomised",
        main = "NMF Rank Survey"
      )
    } else {
      p <- NMF::plot(
        estim.r,
        what = what,
        main = "NMF Rank Survey"
      )
    }

    print(p)

    if (save_plots) {
      pdf(
        paste0(plot_basename, "_survey.pdf"),
        bg = "white",
        pointsize = 9,
        width = 6,
        height = 6,
        paper = "special"
      )
      print(p)
      dev.off()
      if (verbose) message("Created ", paste0(plot_basename, "_survey.pdf"))
    }

    return(
      list(
        nmfEstimate = estim.r,
        nmfEstimate.random = estim.r.random,
        survey = nmf.sum,
        survey.random = nmf.sum.random,
        survey_plot = p
      )
    )
  }
