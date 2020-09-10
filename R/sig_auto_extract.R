#' Extract Signatures through the Automatic Relevance Determination Technique
#'
#' A bayesian variant of NMF algorithm to enable optimal inferences for the
#' number of signatures through the automatic relevance determination technique.
#' This functions delevers highly interpretable and sparse representations for
#' both signature profiles and attributions at a balance between data fitting and
#' model complexity (this method may introduce more signatures than expected,
#' especially for copy number signatures (thus **I don't recommend you to use this feature
#' to extract copy number signatures**)). See detail part and references for more.
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
#' @inheritParams sig_estimate
#' @inheritParams sig_tally
#' @inheritParams sig_extract
#' @param result_prefix prefix for result data files.
#' @param destdir path to save data runs, default is `tempdir()`.
#' @param method default is "L1W.L2H", which uses an exponential prior for W and
#' a half-normal prior for H (This method is used by PCAWG project, see reference #3).
#' You can also use "L1KL" to set expoential priors for both W and H, and "L2KL" to
#' set half-normal priors for both W and H. The latter two methods are originally
#' implemented by [SignatureAnalyzer software](https://software.broadinstitute.org/cancer/cga/msp).
#' @param strategy the selection strategy for returned data. Set 'stable' for getting optimal
#' result from the most frequent K. Set 'optimal' for getting optimal result from all Ks.
#' If you want select other solution, please check [get_bayesian_result].
#' @param K0 number of initial signatures.
#' @param nrun number of independent simulations.
#' @param niter the maximum number of iterations.
#' @param tol tolerance for convergence.
#' @param skip if `TRUE`, it will skip running a previous stored result. This can be used to
#' extend run times, e.g. you try running 10 times firstly and then you want to extend it to
#' 20 times.
#' @param recover if `TRUE`, try to recover result from previous runs based on input `result_prefix`,
#' `destdir` and `nrun`. This is pretty useful for reproducing result. Please use `skip` if you want
#' to recover an unfinished job.
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
#' load(system.file("extdata", "toy_copynumber_tally_M.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' res <- sig_auto_extract(cn_tally_M$nmf_matrix, result_prefix = "Test_copynumber", nrun = 1)
#' # At default, all run files are stored in tempdir()
#' dir(tempdir(), pattern = "Test_copynumber")
#' @testexamples
#' expect_s3_class(res, "Signature")
#' res <- sig_auto_extract(cn_tally_M$nmf_matrix, result_prefix = "test",
#'                         method = "L1KL", nrun = 1)
#' expect_s3_class(res, "Signature")
#' res <- sig_auto_extract(cn_tally_M$nmf_matrix, result_prefix = "test",
#'                         method = "L2KL", nrun = 1)
#' expect_s3_class(res, "Signature")
#' @seealso [sig_tally] for getting variation matrix,
#' [sig_extract] for extracting signatures using **NMF** package, [sig_estimate] for
#' estimating signature number for [sig_extract].
sig_auto_extract <- function(nmf_matrix = NULL,
                             result_prefix = "BayesNMF",
                             destdir = tempdir(),
                             method = c("L1W.L2H", "L1KL", "L2KL"),
                             strategy = c("optimal", "stable"),
                             K0 = 25,
                             nrun = 10,
                             niter = 2e5,
                             tol = 1e-07,
                             cores = 1,
                             optimize = FALSE,
                             skip = FALSE,
                             recover = FALSE) {
  on.exit(invisible(gc())) # clean when exit
  method <- match.arg(method)
  strategy <- match.arg(strategy)
  if (!dir.exists(destdir)) dir.create(destdir, recursive = TRUE)

  filelist <- file.path(destdir, paste(result_prefix, seq_len(nrun), "rds", sep = "."))

  if (recover) {
    message("Recover mode is on, check if all files exist...")
    all_exist <- all(file.exists(filelist))
    if (!all_exist) {
      stop("Recover failed, cannot find previous result files, please run with 'recover = FALSE'.",
        call. = FALSE
      )
    }
    message("Yup! Recovering...")
  }

  if (!recover) {
    nmf_matrix <- t(nmf_matrix) # rows for mutation types and columns for samples

    ii <- colSums(nmf_matrix) < 0.01
    if (any(ii)) {
      message(
        "The follow samples dropped due to null catalogue:\n\t",
        paste0(colnames(nmf_matrix)[ii], collapse = ", ")
      )
      nmf_matrix <- nmf_matrix[, !ii, drop = FALSE]
    }

    oplan <- future::plan()
    future::plan("multiprocess", workers = cores)
    on.exit(future::plan(oplan), add = TRUE)

    furrr::future_map(seq_len(nrun), function(i, method, filelist, skip) {
      if (skip & file.exists(filelist[i])) {
        message("Run #", i, " exists, skipping...")
      } else {
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
      }
    }, method = method, filelist = filelist, skip = skip, .progress = TRUE)
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

  if (strategy == "stable") {
    # Find stable K
    best <- names(sort(table(summary.run$K), decreasing = TRUE))[1] %>%
      as.integer()
    best <- max(1, best, na.rm = TRUE)
  } else {
    best <- max(1, summary.run$K[1], na.rm = TRUE)
  }

  best_row <- dplyr::filter(summary.run, .data$K == best) %>%
    head(1)

  message("Select Run ", best_row$Run, ", which K = ", best_row$K, " as best solution.")
  best_solution <- get_bayesian_result(best_row)

  has_cn <- grepl("^CN[^C]", rownames(best_solution$Signature)) | startsWith(rownames(best_solution$Signature), "copynumber")
  mat <- nmf_matrix
  if (optimize) {
    message("Refit the denovo signatures with QP.")
    ## Optimize signature exposure
    if (any(has_cn)) {
      mat_cn <- mat[has_cn, ]
      W_cn <- best_solution$Signature[has_cn, ]
      W_cn <- apply(W_cn, 2, function(x) x / sum(x))

      ## Call LCD
      best_solution$Exposure <- sig_fit(
        catalogue_matrix = mat_cn,
        sig = W_cn,
        method = "QP",
        mode = "copynumber"
      )
    } else {
      ## Call LCD
      best_solution$Exposure <- sig_fit(
        catalogue_matrix = mat,
        sig = apply(
          best_solution$Signature,
          2, function(x) x / sum(x)
        ),
        method = "QP"
      )
    }

    Exposure <- best_solution$Exposure
    # Handle hyper mutant samples
    hyper_index <- grepl("_\\[hyper\\]_", colnames(Exposure))
    if (sum(hyper_index) > 0) {
      H.hyper <- Exposure[, hyper_index, drop = FALSE]
      H.nonhyper <- Exposure[, !hyper_index, drop = FALSE]
      sample.hyper <- sapply(
        colnames(H.hyper),
        function(x) strsplit(x, "_\\[hyper\\]_")[[1]][[1]]
      )
      unique.hyper <- unique(sample.hyper)
      n.hyper <- length(unique.hyper)
      x.hyper <- array(0, dim = c(nrow(H.hyper), n.hyper))
      for (i in 1:n.hyper) {
        x.hyper[, i] <- rowSums(H.hyper[, sample.hyper %in% unique.hyper[i], drop = FALSE])
      }
      colnames(x.hyper) <- unique.hyper
      rownames(x.hyper) <- rownames(Exposure)
      Exposure <- cbind(H.nonhyper, x.hyper)
      best_solution$Exposure <- Exposure
    }

    best_solution$Exposure.norm <- apply(
      best_solution$Exposure, 2,
      function(x) x / sum(x, na.rm = TRUE)
    )
    # When only one signature
    if (!is.matrix(best_solution$Exposure.norm)) {
      best_solution$Exposure.norm <- matrix(best_solution$Exposure.norm,
        nrow = 1,
        dimnames = list(NULL, names(best_solution$Exposure.norm))
      )
    }
    ## Scale the result
    if (any(has_cn)) {
      best_solution$Signature <- best_solution$Signature.norm * sum(nmf_matrix, na.rm = TRUE)
    } else {
      best_solution$Signature <- best_solution$Signature.norm * sum(best_solution$Exposure, na.rm = TRUE)
    }
  }

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
