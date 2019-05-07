#################################################
##### Signature Analysis Pipeline Series ########
#################################################


# S3 method for signature analysis prepare --------------------------------

#' Prepare variation signature analysis
#'
#' Generate a matrix for NMF de-composition.
#'
#' The result matrix generated further need to transpose before calling NMF
#' if user use [NMF::nmf] by hand.
#' @param object a [CopyNumber] object or [MAF] object or
#' [GenomicVariation] (not support for now) object.
#' @param ... custom setting for operating object. Detail see S3 method for
#' corresponding class (e.g. `CopyNumber`).
#' @return a `list` contains a `matrix` used for NMF de-composition.
#' @author Shixiang Wang
#' @export
#' @examples
#' \donttest{
#' # Load copy number object
#' load(system.file("extdata", "toy_copynumber.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Prepare copy number signature analysis
#' cn_prepare <- sig_prepare(cn)
#' }
#' @family signature analysis series function
sig_prepare <- function(object, ...) {
  UseMethod("sig_prepare")
}

#' @describeIn sig_prepare Signature analysis prepare for CopyNumber object
#' @inheritParams prepare_copynumber
#' @export
sig_prepare.CopyNumber <- function(object, reference_components = FALSE,
                                   cores = 1, seed = 123456,
                                   min_comp = 2, max_comp = 10,
                                   min_prior = 0.001,
                                   model_selection = "BIC",
                                   nrep = 1, niter = 1000, rowIter = 1000, ...) {
  prepare_copynumber(object,
    reference_components = reference_components,
    cores = cores, seed = seed,
    min_comp = min_comp, max_comp = max_comp,
    min_prior = min_prior,
    model_selection = model_selection,
    nrep = nrep, niter = niter, rowIter = rowIter
  )
}

#' @describeIn sig_prepare Signature analysis prepare for CopyNumber object
#' @inheritParams prepare_maf
#' @references Mayakonda, Anand, et al. "Maftools: efficient and comprehensive analysis of somatic variants in cancer." Genome research 28.11 (2018): 1747-1756.
#' @export
sig_prepare.MAF <- function(object, ref_genome = NULL, prefix = NULL,
                            add = TRUE, ignoreChr = NULL, useSyn = TRUE, ...) {
  res <- prepare_maf(object,
    ref_genome = ref_genome, prefix = prefix,
    add = add, ignoreChr = ignoreChr, useSyn = useSyn, fn = NULL
  )
}

#' @describeIn sig_prepare Signature analysis prepare for GenomicVariation object
sig_prepare.GenomicVariation <- function(object, ...) {
  print("Not support right now.")
}



# Signature number estimation ---------------------------------------------

#' Estimate signature number
#'
#' Use **NMF** package to evaluate the optimal number of signatures.
#' Users should `library(NMF)` firstly.
#'
#' The most common approach is to choose the smallest rank for which cophenetic correlation coefficient
#' starts decreasing (Used by this function). Another approach is to choose the rank for which the plot
#' of the residual sum of squares (RSS) between the input matrix and its estimate shows an inflection point.
#' More custom features please directly use [NMF::nmfEstimateRank].
#' @param nmf_matrix a `matrix` used for NMF decomposition (with rownames and colnames),
#'  generate from [sig_prepare] function.
#' @param range a `numeric` vector containing the ranks of factorization to try. Note that duplicates are removed and values are sorted in increasing order. The results are notably returned in this order.
#' @param nrun a `numeric` giving the number of run to perform for each value in `range`, `nrun` set to 30~50 is enough to achieve robust result.
#' @param what a character vector whose elements partially match one of the following item, which correspond to the measures computed by summary on each – multi-run – NMF result: ‘all’, ‘cophenetic’, ‘rss’, ‘residuals’, ‘dispersion’, ‘evar’, ‘silhouette’ (and more specific \*.coef, \*.basis, \*.consensus), ‘sparseness’ (and more specific \*.coef, \*.basis). It specifies which measure must be plotted (what='all' plots all the measures).
#' @param cores number of cpu cores to run NMF.
#' @param seed specification of the starting point or seeding method, which will compute a starting point, usually using data from the target matrix in order to provide a good guess.
#' @param use_random Should generate random data from input to test measurements. Default is `TRUE`.
#' @param save_plots if `TRUE`, save plots to local machine.
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
#' @family signature analysis series function
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

    if (save_plots) {
      pdf(
        paste0(plot_basename, "_consensus.pdf"),
        bg = "white",
        pointsize = 9,
        width = 12,
        height = 12,
        paper = "special",
        onefile = FALSE
      )
      NMF::consensusmap(estim.r)
      dev.off()
      if (verbose) message("Created ", paste0(plot_basename, "_consensus.pdf"))
    }

    nmf.sum <- NMF::summary(estim.r) # Get summary of estimates
    if (verbose) {
      message("Estimation of rank based on observed data.")
      print(nmf.sum)
    }
    # nmf.sum$diff = c(0, diff(nmf.sum$cophenetic))
    # bestFit = nmf.sum$rank[which(nmf.sum$diff < 0)][1]
    # #bestFit = nmf.sum[diff < 0, rank][1] #First point where cophenetic correlation coefficient starts decreasing
    # if (is.na(bestFit) | is.null(bestFit)) {
    #   stop("Cannot find appropriate rank, please reset your arguments.")
    # }
    #
    # # https://blog.csdn.net/YJJ18636810884/article/details/83214566
    # # https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-367
    # message(
    #   paste(
    #     'Using ',
    #     bestFit,
    #     ' as a best-fit rank based on decreasing cophenetic correlation coefficient.',
    #     sep = ''
    #   )
    # )
    # n = as.numeric(bestFit)

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


# Extract signatures ------------------------------------------------------

#' Extract variation signatures
#'
#' Do NMF de-composition and then extract signatures.
#' @inheritParams sig_estimate
#' @param n_sig number of signature. Please run [sig_prepare] to select a suitable value.
#' @param mode variation type to decompose, currently support "copynumber" or "mutation".
#' @author Shixiang Wang
#' @references Gaujoux, Renaud, and Cathal Seoighe. "A flexible R package for nonnegative matrix factorization." BMC bioinformatics 11.1 (2010): 367.
#' @references Mayakonda, Anand, et al. "Maftools: efficient and comprehensive analysis of somatic variants in cancer." Genome research 28.11 (2018): 1747-1756.
#' @return a `list` contains NMF object, signature matrix and activity matrix.
#' @import NMF
#' @export
#' @examples
#' \donttest{
#' # Load copy number prepare object
#' load(system.file("extdata", "toy_copynumber_prepare.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Extract copy number signatures
#' res <- sig_extract(cn_prepare$nmf_matrix, 2, mode = "copynumber", nrun = 1)
#' }
#' @family signature analysis series function
sig_extract <- function(nmf_matrix,
                        n_sig,
                        mode = c("copynumber", "mutation"),
                        nrun = 10,
                        cores = 1,
                        method = "brunet",
                        pConstant = NULL,
                        seed = 123456) {
  mode <- match.arg(mode)

  # transpose matrix
  mat <- t(nmf_matrix)

  # To avoid error due to non-conformable arrays
  if (!is.null(pConstant)) {
    if (pConstant < 0 | pConstant == 0) {
      stop("pConstant must be > 0")
    }
    mat <- mat + pConstant
  }

  nmf.res <- NMF::nmf(
    mat,
    n_sig,
    seed = seed,
    nrun = nrun,
    method = method,
    .opt = paste0("p", cores)
  )


  # Signatures
  w <- NMF::basis(nmf.res)
  w <- apply(w, 2, function(x) x / sum(x)) # Scale the signatures (basis)
  colnames(w) <- paste("Signature", 1:ncol(w), sep = "_")

  # Contribution
  h <- NMF::coef(nmf.res)
  colnames(h) <- colnames(mat) # correct colnames
  # For single signature, contribution will be 100% per sample
  if (n_sig == 1) {
    h <- h / h
    rownames(h) <- paste("Signature", "1", sep = "_")
  } else {
    h <- apply(h, 2, function(x) x / sum(x)) # Scale contributions (coefs)
    rownames(h) <- paste("Signature", 1:nrow(h), sep = "_")
  }

  if (mode == "mutation") {
    #- Copy from maftools, should update them all if something changed in maftools source
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< BEGIN
    sigs <- data.table::fread(
      input = system.file("extdata", "signatures.txt", package = "maftools"),
      stringsAsFactors = FALSE,
      data.table = FALSE
    )
    colnames(sigs) <- gsub(
      pattern = " ",
      replacement = "_",
      x = colnames(sigs)
    )
    rownames(sigs) <- sigs$Somatic_Mutation_Type
    sigs <- sigs[, -c(1:3)]
    # sigs = sigs[,1:22] #use only first 21 validated sigantures
    sigs <- sigs[rownames(w), ]

    aetiology <- structure(
      list(
        aetiology = c(
          "spontaneous deamination of 5-methylcytosine",
          "APOBEC Cytidine Deaminase (C>T)",
          "defects in DNA-DSB repair by HR",
          "exposure to tobacco (smoking) mutagens",
          "Unknown",
          "defective DNA mismatch repair",
          "UV exposure",
          "Unknown",
          "defects in polymerase-eta",
          "defects in polymerase POLE",
          "exposure to alkylating agents",
          "Unknown",
          "APOBEC Cytidine Deaminase (C>G)",
          "Unknown",
          "defective DNA mismatch repair",
          "Unknown",
          "Unknown",
          "Unknown",
          "Unknown",
          "defective DNA mismatch repair",
          "unknown",
          "exposure to aristolochic acid",
          "Unknown",
          "exposures to aflatoxin",
          "Unknown",
          "defective DNA mismatch repair",
          "Unknown",
          "Unknown",
          "exposure to tobacco (chewing) mutagens",
          "Unknown"
        )
      ),
      .Names = "aetiology",
      row.names = c(
        "Signature_1",
        "Signature_2",
        "Signature_3",
        "Signature_4",
        "Signature_5",
        "Signature_6",
        "Signature_7",
        "Signature_8",
        "Signature_9",
        "Signature_10",
        "Signature_11",
        "Signature_12",
        "Signature_13",
        "Signature_14",
        "Signature_15",
        "Signature_16",
        "Signature_17",
        "Signature_18",
        "Signature_19",
        "Signature_20",
        "Signature_21",
        "Signature_22",
        "Signature_23",
        "Signature_24",
        "Signature_25",
        "Signature_26",
        "Signature_27",
        "Signature_28",
        "Signature_29",
        "Signature_30"
      ),
      class = "data.frame"
    )

    message(
      "Comparing against experimentally validated 30 signatures.. (See http://cancer.sanger.ac.uk/cosmic/signatures for details.)"
    )
    # corMat = c()
    coSineMat <- c()
    for (i in 1:ncol(w)) {
      sig <- w[, i]
      coSineMat <- rbind(coSineMat, apply(sigs, 2, function(x) {
        round(crossprod(sig, x) / sqrt(crossprod(x) * crossprod(sig)), digits = 3) # Estimate cosine similarity against all 30 signatures
      }))
      # corMat = rbind(corMat, apply(sigs, 2, function(x) cor.test(x, sig)$estimate[[1]])) #Calulate correlation coeff.
    }
    # rownames(corMat) = colnames(w)
    rownames(coSineMat) <- colnames(w)

    for (i in 1:nrow(coSineMat)) {
      ae <- aetiology[names(which(coSineMat[i, ] == max(coSineMat[i, ]))), ]
      ae <- paste0(
        "Aetiology: ",
        ae,
        " [cosine-similarity: ",
        max(coSineMat[i, ]),
        "]"
      )
      message(
        "Found ",
        rownames(coSineMat)[i],
        " most similar to validated ",
        names(which(coSineMat[i, ] == max(coSineMat[i, ]))),
        ". ",
        ae,
        sep = " "
      )
    }

    return(list(
      signature = w,
      activity = h,
      coSineSimMat = coSineMat,
      nmfObj = nmf.res
    ))
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< END
  } else {
    return(list(
      signature = w,
      activity = h,
      nmfObj = nmf.res
    ))
  }
}


# Assign sample/features to groups/clusters ----------------------------------------

#' Return sample clustering from NMF run results
#'
#' One of key results from NMF decomposition is to cluster samples into different
#' groups. This function takes NMF result (a `NMF` object) as input
#' and return the membership in each cluster.
#'
#' X = W x H
#'
#' W is the feature matrix, H is the sample matrix
#' After NMF run, use this function to select import features and assign groups for the two matrix.
#'
#' More detail please see [NMF::predict()].
#' @param nmfObj a `NMF` result object which is an element return from [sig_extract]
#' or run results of **NMF** package.
#' @param type cluster type, could be 'consensus' or 'samples'.
#' @param matchConseOrder if `TRUE`, the result will match order as shown in consensus map
#' when type argument is 'consensus'.
#' @return a `data.table` object
#' @import NMF cluster
#' @export
#' @examples
#' # Load copy number signature
#' load(system.file("extdata", "toy_copynumber_signature.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Assign samples to clusters
#' subtypes <- sig_assign_samples(res$nmfObj, type = "samples")
#' @family signature analysis series function
sig_assign_samples <- function(nmfObj, type = "consensus", matchConseOrder = F) {
  data <- NULL

  # loadNamespace("cluster")
  # loadNamespace("NMF")

  if (type == "consensus") {
    predict.consensus <- predict(nmfObj, what = "consensus")
    silhouette.consensus <- silhouette(nmfObj, what = "consensus")
    # It turns out the factor levels is the NMF_assigned_groups from consensus matrix
    # that matches the original sampleNames(nmfObj) order
    # The attributes(a.predict.consensus)$iOrd is the idx order for it to match the
    # order of the samples in consensusmap(nmfObj). It is just for displaying
    # Therefore, the merged data frame sampleNames(nmfObj) + a.predict.consensus is the final
    # consensus results.
    data <- data.frame(
      Sample_ID = sampleNames(nmfObj),
      nmf_subtypes = predict.consensus,
      sil_width = signif(silhouette.consensus[, "sil_width"], 3)
    )
    # If we want to display as we see in consensusmap, we just need to reoder everything.
    # Now re-order data to match consensusmap sample order
    if (matchConseOrder) {
      sample.order <- attributes(predict.consensus)$iOrd
      data <- data[sample.order, ]
    }
  } else if (type == "samples") {
    predict.samples <- predict(nmfObj, what = "samples", prob = T)
    silhouette.samples <- silhouette(nmfObj, what = "samples")
    data <- data.frame(
      Sample_ID = names(predict.samples$predict),
      nmf_subtypes = predict.samples$predict,
      sil_width = signif(silhouette.samples[, "sil_width"], 3),
      prob = signif(predict.samples$prob, 3)
    )
  } else {
    stop(paste("Wrong type:", type, "Possible options are: 'consensus', 'samples' "))
  }

  data <- data.table::as.data.table(data)
  return(data)
}


# Get signature activity --------------------------------------------------

#' Get signature activity
#'
#' @inheritParams draw_sig_profile
#' @return a `list`
#' @export
#' @examples
#' # Load copy number signature
#' load(system.file("extdata", "toy_copynumber_signature.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Get activity of signatures
#' sig_activity <- sig_get_activity(res$nmfObj)
#' @family signature analysis series function
sig_get_activity <- function(nmfObj) {
  if (!inherits(nmfObj, "NMFfit")) {
    stop("'nmfObj' should be a NMFfit object.")
  }
  # Activity
  h <- NMF::coef(nmfObj)
  rownames(h) <- paste("Signature", 1:nrow(h), sep = "_")
  # colnames(h) = colnames(mat) #correct colnames
  # For single signature, contribution will be 100% per sample
  if (nrow(h) == 1) {
    h.norm <- h / h
    rownames(h.norm) <- paste("Signature", "1", sep = "_")
  } else {
    h.norm <- apply(h, 2, function(x) x / sum(x)) # Scale contributions (coefs)
    rownames(h.norm) <- paste("Signature", 1:nrow(h.norm), sep = "_")
  }

  h <- as.data.frame(h)
  h.norm <- as.data.frame(h.norm)

  list(
    relative = h.norm,
    absolute = h
  )
}



# Get correlation matrix between signature activities -------------------------

#' Get correlation matrix between signature activities
#'
#' Compute correlation matrix and corresponding statistical test values between signature activies.
#'
#' @param cn_activity activity of copy number signature, a `list`, obtain it from
#' [sig_get_activity] function.
#' @param snv_activity activity of mutational signature, a `list`, obtain it from
#' [sig_get_activity] function.
#' @param type one of 'absolute' and 'relative'.
#' @param ... other arguments pass to [corrplot::cor.mtest()].
#' @return a `list`.
#' @author Shixiang Wang
#' @export
#' @examples
#' # Load copy number signature
#' load(system.file("extdata", "toy_copynumber_signature.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Get activity of signatures
#' sig_activity <- sig_get_activity(res$nmfObj)
#' # Get correlation matrix between signature activities
#' sig_cor <- sig_get_correlation(sig_activity)
#' @family signature analysis series function
sig_get_correlation <- function(cn_activity = NULL, snv_activity = NULL,
                                type = c("absolute", "relative"),
                                ...) {
  if (is.null(cn_activity) & is.null(snv_activity)) {
    stop("At least one of 'cn_activity' and 'snv_activity' should be setted.")
  }

  type <- match.arg(type)

  gen_mat <- function(cn_activity, snv_activity, type) {
    if (!is.null(cn_activity) & !is.null(snv_activity)) {
      cn_mat <- t(cn_activity[[type]])
      snv_mat <- t(snv_activity[[type]])
      colnames(cn_mat) <- paste0("CN_", colnames(cn_mat))
      colnames(snv_mat) <- paste0("SNV_", colnames(snv_mat))

      samps <- base::intersect(rownames(cn_mat), rownames(snv_mat))

      samps.diff.cn <- base::setdiff(rownames(cn_mat), samps)
      if (length(samps.diff.cn) > 0) {
        message("Following samples removed from copy number signature activity matrix:")
        print(samps.diff.cn)
      }
      samps.diff.snv <- base::setdiff(rownames(snv_mat), samps)
      if (length(samps.diff.snv) > 0) {
        message("Following samples removed from mutation signature activity matrix:")
        print(samps.diff.snv)
      }

      mat <- cbind(cn_mat[samps, ], snv_mat[samps, ])
      return(mat)
    } else if (!is.null(cn_activity)) {
      mat <- t(cn_activity[[type]])
      return(mat)
    } else {
      mat <- t(snv_activity[[type]])
      return(mat)
    }
  }

  mat <- gen_mat(cn_activity, snv_activity, type)

  corr <- list(correlation = stats::cor(mat))

  c(
    corr,
    corrplot::cor.mtest(mat, ...)
  )
}



# Get similarity between signatures ---------------------------------------


#' Get similarity between signatures
#'
#' Obtain similarity for two signature matrix, their rows (features) must match.
#'
#' @param sig1 signature 1, result of [sig_extract] function or
#' a `data.frame` which rows are features and columns are signatures
#' (colnames are necessary).
#' @param sig2 signature 2, result of [sig_extract] function or
#' a `data.frame` which rows are features and columns are signatures
#' (colnames are necessary).
#' @param type could be "cos" (default, cosine similarity) or "cor" (correlation).
#' @author Shixiang Wang
#' @return a `matrix` which rownames from `sig1` and colnames from `sig2`.
#' @export
#' @family signature analysis series function
sig_get_similarity <- function(sig1, sig2, type = c("cos", "cor")) {
  if (is.list(sig1)) {
    mat1 <- sig1[["signature"]]
  } else {
    mat1 <- sig1
  }

  if (is.list(sig2)) {
    mat2 <- sig2[["signature"]]
  } else {
    mat2 <- sig2
  }

  type <- match.arg(type)

  ResMat <- c()
  for (i in 1:ncol(mat1)) {
    sig <- mat1[, i]
    if (type == "cos") {
      ResMat <- rbind(ResMat, apply(mat2, 2, function(x) {
        round(crossprod(sig, x) / sqrt(crossprod(x) * crossprod(sig)), digits = 3)
      }))
    } else {
      ResMat <- rbind(ResMat, apply(mat2, 2, function(x) cor.test(x, sig)$estimate[[1]])) # Calulate correlation coeff.
    }
  }

  rownames(ResMat) <- colnames(mat1)

  message("rownames come from sig1 and colnames come from sig2.")
  ResMat
}

# Get signature subtype summary -------------------------------------------

#' Get summary of signature subtypes
#'
#' Summarize genotypes/phenotypes based on signature subtypes. For categorical
#' type, calculate fisher p value (using [stats::fisher.test]) and count table.
#' For continuous type, calculate anova p value (using [stats::aov]),
#' summary table and Tukey Honest significant difference (using [stats::TukeyHSD]).
#' The result of this function can be plotted by [draw_subtypes_comparison()].
#'
#' @param data a `data.frame` contains signature subtypes and genotypes/phenotypes
#' (including categorical and continuous type data) want to analyze. User need to
#' construct this `data.frame` by him/herself.
#' @param col_subtype column name of signature subtypes.
#' @param cols_to_summary column names of genotypes/phenotypes want to summarize based on subtypes.
#' @param type a characater vector with length same as `cols_to_summary`,
#' 'ca' for categorical type and 'co' for continuous type.
#' @param verbose if `TRUE`, print extra information.
#' @author Shixiang Wang <w_shixiang@163.com>
#' @return a `list` contains data, summary, p value etc..
#' @export
#' @examples
#' # Load copy number signature
#' load(system.file("extdata", "toy_copynumber_signature.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Assign samples to clusters
#' subtypes <- sig_assign_samples(res$nmfObj, type = "samples")
#'
#' set.seed(1234)
#' # Add custom groups
#' subtypes$new_group <- sample(c("1", "2", "3", "4"), size = nrow(subtypes), replace = TRUE)
#' # Summarize subtypes
#' subtypes.sum <- sig_summarize_subtypes(subtypes[, -1],
#'   col_subtype = "nmf_subtypes",
#'   cols_to_summary = colnames(subtypes[, -1])[c(-1, -2)],
#'   type = c("co", "ca"), verbose = TRUE
#' )
#' @family signature analysis series function
sig_summarize_subtypes <- function(data, col_subtype, cols_to_summary,
                                   type = "ca", verbose = FALSE) {
  if (!all(type %in% c("ca", "co"))) {
    stop("all elements in 'type' must be 'ca' for 'categorical' and 'co' for 'continuous'.")
  }
  if (!is.data.frame(data)) stop("'data' must be a data.frame object.")
  data.table::setDT(data)

  # subset
  data <- data[, c(col_subtype, cols_to_summary), with = FALSE]
  colnames(data)[1] <- "subtype"
  data <- data[!is.na(data[["subtype"]])]

  do_summary <- function(col, type = c("ca", "co"),
                           verbose = FALSE) {
    type <- match.arg(type)

    df <- data[, c("subtype", col), with = FALSE]
    df <- df[!is.na(df[[col]])]

    if (type == "ca") {
      if (verbose) message("Treat ", col, " as categorical variable.")

      table_df <- table(df[["subtype"]], df[[col]])

      table_p <- tryCatch({
        test <- fisher.test(table_df)
        test[["p.value"]]
      }, error = function(e) {
        NA
      })

      list(data = df, table = table_df, p_value = table_p, type = "categorical", extra = NA)
    } else {
      if (verbose) message("Treat ", col, " as continuous variable.")

      table_df <- summary(df)

      fit <- tryCatch({
        stats::aov(as.formula(paste0(col, " ~ subtype")), data = df)
      }, error = function(e) {
        NA
      })

      if (inherits(fit, "aov")) {
        p_value <- summary(fit)[[1]][["Pr(>F)"]][1] # get anova p value
        extra <- stats::TukeyHSD(fit)[[1]]
      } else {
        p_value <- NA
        extra <- NA
      }
      list(data = df, table = table_df, p_value = p_value, type = "continuous", extra = extra)
    }
  }

  res <- Map(do_summary, cols_to_summary, type, verbose)
  names(res) <- cols_to_summary
  res
}
