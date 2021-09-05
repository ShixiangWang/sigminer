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
#' Set 'ms' for getting result with maximum mean cosine similarity with provided reference
#' signatures. See `ref_sigs` option for details.
#' If you want select other solution, please check [get_bayesian_result].
#' @param ref_sigs A Signature object or matrix or string for specifying
#' reference signatures, only used when `strategy = 'ms'`.
#' See `Signature` and `sig_db` options in [get_sig_similarity] for details.
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
#' load(system.file("extdata", "toy_copynumber_tally_W.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' res <- sig_auto_extract(cn_tally_W$nmf_matrix, result_prefix = "Test_copynumber", nrun = 1)
#' # At default, all run files are stored in tempdir()
#' dir(tempdir(), pattern = "Test_copynumber")
#' \donttest{
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read_maf(maf = laml.maf)
#' mt_tally <- sig_tally(
#'   laml,
#'   ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
#'   use_syn = TRUE
#' )
#'
#' x <- sig_auto_extract(mt_tally$nmf_matrix,
#'   strategy = "ms", nrun = 3, ref_sigs = "legacy"
#' )
#' x
#' }
#' @testexamples
#' expect_s3_class(res, "Signature")
#' res <- sig_auto_extract(cn_tally_W$nmf_matrix, result_prefix = "test",
#'                         method = "L1KL", nrun = 1)
#' expect_s3_class(res, "Signature")
#' res <- sig_auto_extract(cn_tally_W$nmf_matrix, result_prefix = "test",
#'                         method = "L2KL", nrun = 1)
#' expect_s3_class(res, "Signature")
#' expect_s3_class(x, "Signature")
#' @seealso [sig_tally] for getting variation matrix,
#' [sig_extract] for extracting signatures using **NMF** package, [sig_estimate] for
#' estimating signature number for [sig_extract].
sig_auto_extract <- function(nmf_matrix = NULL,
                             result_prefix = "BayesNMF",
                             destdir = tempdir(),
                             method = c("L1W.L2H", "L1KL", "L2KL"),
                             strategy = c("stable", "optimal", "ms"),
                             ref_sigs = NULL,
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
  if (strategy == "ms") {
    if (is.null(ref_sigs)) {
      stop("When strategy set to 'ms', the ref_sigs cannot be NULL!")
    }
  }
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
    future::plan(set_future_strategy(), workers = cores, gc = TRUE)
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
    },
    method = method, filelist = filelist, skip = skip,
    .progress = TRUE,
    .options = furrr::furrr_options(seed = TRUE)
    )
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

  if (strategy != "ms") {
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
  } else {
    # Select the optimal solution by maximizing mean cosine similarity between
    # extracted signatures and reference signatures.

    sim_list <- lapply(1:nrow(summary.run), function(x, ref) {
      xz <- get_bayesian_result(summary.run[x, ])
      if (is.character(ref)) {
        sim <- suppressMessages(get_sig_similarity(xz$Signature.norm,
          sig_db = ref
        ))
      } else {
        sim <- suppressMessages(get_sig_similarity(xz$Signature.norm, ref))
      }
      sim <- sim$similarity
      yid <- apply(sim, 1, which.max)
      y <- apply(sim, 1, max)
      names(y) <- colnames(sim)[yid]
      y
    }, ref = ref_sigs)

    names(sim_list) <- paste0(
      "Run#", summary.run$Run, ":",
      "K#", summary.run$K
    )

    sims <- sapply(sim_list, mean)
    summary.run$similarity_to_ref <- as.numeric(sims)
    ind <- order(sims, decreasing = TRUE)
    message("Solutions ordered by mean cosine similarity to references:")
    print(sims[ind])
    message("Details:")
    print(sim_list)
    best_row <- summary.run[ind[1], ]
    summary.run <- summary.run %>%
      dplyr::arrange(dplyr::desc(.data$similarity_to_ref))
  }

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

    # Handle hyper mutant samples
    best_solution$Exposure <- collapse_hyper_records(best_solution$Exposure)

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
  attr(res, "call_method") <- "BayesianNMF"

  res
}
