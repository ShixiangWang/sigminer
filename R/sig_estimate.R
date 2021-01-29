#' Estimate Signature Number
#'
#' Use **NMF** package to evaluate the optimal number of signatures.
#' This is used along with [sig_extract].
#' Users should `library(NMF)` firstly. If NMF objects are returned,
#' the result can be further visualized by NMF plot methods like
#' `NMF::consensusmap()` and `NMF::basismap()`.
#'
#' The most common approach is to choose the smallest rank for which cophenetic correlation coefficient
#' starts decreasing (Used by this function). Another approach is to choose the rank for which the plot
#' of the residual sum of squares (RSS) between the input matrix and its estimate shows an inflection point.
#' More custom features please directly use [NMF::nmfEstimateRank].
#'
#' @param nmf_matrix a `matrix` used for NMF decomposition with rows indicate samples and columns indicate components.
#' @param range a `numeric` vector containing the ranks of factorization to try. Note that duplicates are removed
#' and values are sorted in increasing order. The results are notably returned in this order.
#' @param keep_nmfObj default is `FALSE`, if `TRUE`, keep NMF objects from runs, and the result may be huge.
#' @param nrun a `numeric` giving the number of run to perform for each value in `range`, `nrun` set to 30~50 is
#' enough to achieve robust result.
#' @param what a character vector whose elements partially match one of the following item, which correspond to
#' the measures computed by summary on each multi-run NMF result: ‘all’, ‘cophenetic’, ‘rss’, ‘residuals’,
#'  ‘dispersion’, ‘evar’, ‘silhouette’ (and more specific .coef, .basis, .consensus), ‘sparseness’
#'  (and more specific .coef, .basis). It specifies which measure must be plotted
#'  (what='all' plots all the measures).
#'
#' @param cores number of cpu cores to run NMF.
#' @param seed specification of the starting point or seeding method, which will compute a starting point,
#'  usually using data from the target matrix in order to provide a good guess.
#' @param use_random Should generate random data from input to test measurements. Default is `TRUE`.
#' @param save_plots if `TRUE`, save signature number survey plot to local machine.
#' @param plot_basename when save plots, set custom basename for file path.
#' @param method specification of the NMF algorithm. Use 'brunet' as default.
#' Available methods for NMF decompositions are 'brunet', 'lee', 'ls-nmf', 'nsNMF', 'offset'.
#' @param verbose if `TRUE`, print extra message.
#' @author Shixiang Wang
#' @references Gaujoux, Renaud, and Cathal Seoighe. "A flexible R package for nonnegative matrix factorization." BMC bioinformatics 11.1 (2010): 367.
#' @return a `list` contains information of NMF run and rank survey.
#' @export
#' @examples
#' load(system.file("extdata", "toy_copynumber_tally_M.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' \dontrun{
#' cn_estimate <- sig_estimate(cn_tally_M$nmf_matrix,
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
           use_random = FALSE,
           method = "brunet",
           seed = 123456,
           cores = 1,
           keep_nmfObj = FALSE,
           save_plots = FALSE,
           plot_basename = file.path(tempdir(), "nmf"),
           what = "all",
           verbose = FALSE) {
    eval(parse(text = "suppressMessages(library('NMF'))"))
    if (cores > 1) cores <- min(cores, future::availableCores())
    mat <- t(nmf_matrix)

    ii <- colSums(mat) < 0.01
    if (any(ii)) {
      message(
        "The follow samples dropped due to null catalogue:\n\t",
        paste0(colnames(mat)[ii], collapse = ", ")
      )
      mat <- mat[, !ii, drop = FALSE]
    }

    # To avoid error due to NMF
    mat <- check_nmf_matrix(mat)

    if (cores > 1) {
      estim.r <-
        NMF::nmfEstimateRank(
          mat,
          range,
          method = method,
          nrun = nrun,
          verbose = verbose,
          seed = seed,
          .opt = paste0("p", cores)
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
            .opt = paste0("p", cores)
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

    if (save_plots) {
      if (use_random) {
        p <- NMF::plot(
          estim.r,
          estim.r.random,
          what = what,
          xname = "Observed",
          yname = "Randomised",
          xlab = "Number of signature",
          main = "Signature number survey using NMF package"
        )
      } else {
        p <- NMF::plot(
          estim.r,
          what = what,
          xlab = "Number of signature",
          main = "Signature number survey using NMF package"
        )
      }

      destdir <- dirname(plot_basename)
      if (!dir.exists(destdir)) dir.create(destdir, recursive = TRUE)
      pdf(
        paste0(plot_basename, "_survey.pdf"),
        bg = "white",
        pointsize = 9,
        width = 6 + abs((nrow(nmf.sum) - 6) / 3),
        height = 6,
        paper = "special"
      )
      print(p)
      dev.off()
      if (verbose) message("Created ", paste0(plot_basename, "_survey.pdf"))
    }

    if (keep_nmfObj) {
      res <- list(
        nmfEstimate = estim.r,
        nmfEstimate.random = estim.r.random,
        survey = nmf.sum,
        survey.random = nmf.sum.random
      )
    } else {
      res <- list(
        survey = nmf.sum,
        survey.random = nmf.sum.random
      )
    }

    class(res) <- "Survey"
    res
  }
