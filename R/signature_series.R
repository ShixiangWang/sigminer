#################################################
##### Signature Analysis Pipeline Series ########
#################################################


# S3 method for signature analysis prepare --------------------------------

#' Variation signature analysis prepare generic
#'
#' Generate a matrix for NMF de-composition.
#'
#' The result matrix generated further need to transpose before calling NMF
#' if user use [NMF::nmf] by hand.
#' @param object a [CopyNumber] object or [MAF] object or [GenomicVariation] object.
#' @param ... custom setting for operating object.
#' @return a `matrix` used for NMF de-composition.
#' @author Shixiang Wang
#' @export
#' @family signature analysis series function
sig_prepare = function(object, ...){
  UseMethod("sig_prepare")
}

#' @describeIn sig_prepare Signature analysis prepare for CopyNumber object
#' @inheritParams prepare_copynumber
sig_prepare.CopyNumber = function(object, reference_components = FALSE,
                                  cores = 1, seed = 123456,
                                  min_comp = 2, max_comp = 10,
                                  min_prior = 0.001,
                                  model_selection = "BIC",
                                  nrep = 1, niter = 1000, rowIter = 1000){
  prepare_copynumber(object, reference_components = reference_components,
                     cores = cores, seed = seed,
                     min_comp = min_comp, max_comp = max_comp,
                     min_prior = min_prior,
                     model_selection = model_selection,
                     nrep = nrep, niter = niter, rowIter = rowIter)
}

#' @describeIn sig_prepare Signature analysis prepare for CopyNumber object
#' @inheritParams prepare_maf
sig_prepare.MAF = function(object, ref_genome = NULL, prefix = NULL,
                           add = TRUE, ignoreChr = NULL, useSyn = TRUE){
  res = prepare_maf(object, ref_genome = ref_genome, prefix = prefix,
                    add = add, ignoreChr = ignoreChr, useSyn = useSyn, fn = NULL)
}

# # sig_prepare Signature analysis prepare for GenomicVariation object
# sig_prepare.GenomicVariation = function(object) {
#   print("Not support right now.")
# }



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
#' @importFrom grDevices pdf dev.off
#' @return a `list` contains information of NMF run and rank survey.
#' @import NMF
#' @export
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
           plot_basename = "nmf",
           method = "brunet",
           pConstant = NULL,
           verbose = FALSE)
  {

    # loadNamespace("utils")
    # loadNamespace("registry")
    # loadNamespace("pkgmaker")
    # loadNamespace("rngtools")
    # loadNamespace("cluster")
    # loadNamespace("NMF")

    mat = t(nmf_matrix)

    #To avoid error due to non-conformable arrays
    if(!is.null(pConstant)){
      if(pConstant < 0 | pConstant == 0){
        stop("pConstant must be > 0")
      }
      mat = mat+pConstant
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
          seed = seed)
    }

    if (save_plots) {
      pdf(
        paste0(plot_basename, "_consensus.pdf"),
        bg = 'white',
        pointsize = 9,
        width = 12,
        height = 12,
        paper = "special"
      )
      NMF::consensusmap(estim.r)
      dev.off()
      if (verbose) message('created ', paste0(plot_basename, "_consensus.pdf"))
    }

    nmf.sum = NMF::summary(estim.r) # Get summary of estimates
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
            seed = seed)
      }

      nmf.sum.random = NMF::summary(estim.r.random) # Get summary of estimates
      if (verbose) {
        message("Estimation of rank based on random data.")
        print(nmf.sum.random)
      }
    } else {
      nmf.sum.random = NULL
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
        bg = 'white',
        pointsize = 9,
        width = 6,
        height = 6,
        paper = "special"
      )
      print(p)
      dev.off()
      if (verbose) message('created ', paste0(plot_basename, "_survey.pdf"))
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
#' @return a `list` contains NMF object, signature matrix and contribution matrix.
#' @import NMF
#' @export
#' @family signature analysis series function
sig_extract = function(nmf_matrix,
                       n_sig,
                       mode = c("copynumber", "mutation"),
                       nrun = 10,
                       cores = 1,
                       method = "brunet",
                       pConstant = NULL,
                       seed = 123456){

  mode = match.arg(mode)

  #transpose matrix
  mat = t(nmf_matrix)

  #To avoid error due to non-conformable arrays
  if(!is.null(pConstant)){
    if(pConstant < 0 | pConstant == 0){
      stop("pConstant must be > 0")
    }
    mat = mat+pConstant
  }

  nmf.res = NMF::nmf(
    mat,
    n_sig,
    seed = seed,
    nrun = nrun,
    method = method,
    .opt = paste0("p", cores)
  )


  #Signatures
  w = NMF::basis(nmf.res)
  w = apply(w, 2, function(x) x/sum(x)) #Scale the signatures (basis)
  colnames(w) = paste('Signature', 1:ncol(w),sep='_')

  #Contribution
  h = NMF::coef(nmf.res)
  colnames(h) = colnames(mat) #correct colnames
  #For single signature, contribution will be 100% per sample
  if(n_sig == 1){
    h = h/h
    rownames(h) = paste('Signature', '1', sep = '_')
  }else{
    h = apply(h, 2, function(x) x/sum(x)) #Scale contributions (coefs)
    rownames(h) = paste('Signature', 1:nrow(h),sep='_')
  }

  if (mode == "mutation") {
    #- Copy from maftools, should update them all if something changed in maftools source
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< BEGIN
    sigs = data.table::fread(
      input = system.file('extdata', 'signatures.txt', package = 'maftools'),
      stringsAsFactors = FALSE,
      data.table = FALSE
    )
    colnames(sigs) = gsub(pattern = ' ',
                          replacement = '_',
                          x = colnames(sigs))
    rownames(sigs) = sigs$Somatic_Mutation_Type
    sigs = sigs[, -c(1:3)]
    #sigs = sigs[,1:22] #use only first 21 validated sigantures
    sigs = sigs[rownames(w), ]

    aetiology = structure(
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
      'Comparing against experimentally validated 30 signatures.. (See http://cancer.sanger.ac.uk/cosmic/signatures for details.)'
    )
    #corMat = c()
    coSineMat = c()
    for (i in 1:ncol(w)) {
      sig = w[, i]
      coSineMat = rbind(coSineMat, apply(sigs, 2, function(x) {
        round(crossprod(sig, x) / sqrt(crossprod(x) * crossprod(sig)), digits = 3) #Estimate cosine similarity against all 30 signatures
      }))
      #corMat = rbind(corMat, apply(sigs, 2, function(x) cor.test(x, sig)$estimate[[1]])) #Calulate correlation coeff.
    }
    #rownames(corMat) = colnames(w)
    rownames(coSineMat) = colnames(w)

    for (i in 1:nrow(coSineMat)) {
      ae = aetiology[names(which(coSineMat[i, ] == max(coSineMat[i, ]))), ]
      ae = paste0("Aetiology: ",
                  ae,
                  " [cosine-similarity: ",
                  max(coSineMat[i, ]),
                  "]")
      message(
        'Found ',
        rownames(coSineMat)[i],
        ' most similar to validated ',
        names(which(coSineMat[i, ] == max(coSineMat[i, ]))),
        '. ',
        ae,
        sep = ' '
      )
    }

    return(list(
      signatures = w,
      contributions = h,
      coSineSimMat = coSineMat,
      nmfObj = nmf.res
    ))
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< END
  } else {
    return(list(
      signatures = w,
      contributions = h,
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
#' @param type cluster type, could be 'consensus' or 'sample'.
#' @param matchConseOrder if `TRUE`, the result will match order as shown in consensus map.
#' @return a `data.frame`
#' @import NMF cluster
#' @export
#' @family signature analysis series function
sig_assign_samples = function(nmfObj, type="consensus", matchConseOrder=F){

  data <- NULL

  #loadNamespace("cluster")
  #loadNamespace("NMF")

  if(type=="consensus"){
    predict.consensus <- predict(nmfObj, what="consensus")
    silhouette.consensus <- silhouette(nmfObj, what="consensus")
    # It turns out the factor levels is the NMF_assigned_groups from consensus matrix
    # that matches the original sampleNames(nmfObj) order
    # The attributes(a.predict.consensus)$iOrd is the idx order for it to match the
    # order of the samples in consensusmap(nmfObj). It is just for displaying
    # Therefore, the merged data frame sampleNames(nmfObj) + a.predict.consensus is the final
    # consensus results.
    data <- data.frame(Sample_ID=sampleNames(nmfObj),
                       nmf_subtypes = predict.consensus,
                       sil_width = signif(silhouette.consensus[, "sil_width"], 3))
    # If we want to display as we see in consensusmap, we just need to reoder everything.
    # Now re-order data to match consensusmap sample order
    if(matchConseOrder){
      sample.order <- attributes(predict.consensus)$iOrd
      data <- data[sample.order, ]
    }
  }else if(type=="samples"){
    predict.samples <- predict(nmfObj, what="samples", prob=T)
    silhouette.samples <- silhouette(nmfObj, what="samples")
    data <- data.frame(Sample_ID=names(predict.samples$predict),
                       nmf_subtypes = predict.samples$predict,
                       sil_width = signif(silhouette.samples[, "sil_width"], 3),
                       prob = signif(predict.samples$prob, 3))
  }else{
    stop(paste("Wrong type:", type, "Possible options are: 'consensus', 'samples' "))
  }
  return(data)
}


#' #' #------------------------------------
#' #' @title Choose optimal number of signatures
#' #' @description This function use **NMF** package to evaluate the optimal number of signatures.
#' #' The most common approach is to choose the smallest rank for which cophenetic correlation coefficient
#' #' starts decreasing (Used by this function). Another approach is to choose the rank for which the plot
#' #' of the residual sum of squares (RSS) between the input matrix and its estimate shows an inflection point.
#' #' @param nmf_matrix a `matrix` use for NMF de-composition, generate from [sig_prepare] function.
#' #' @param nTry the maximal tried number of signatures, default is 12.
#' #' Of note, this value should far less than number of features or samples.
#' #' @param nrun the number of run to perform for each value in range of 2 to `nTry`, default is 10.
#' #' According to **NMF** package documentation, `nrun` set to 30~50 is enough to achieve robust result.
#' #' @inheritParams get_features
#' #' @param seed seed number.
#' #' @param plot logical. If `TRUE`, plot rank survey.
#' #' @param consensusmap_name a character, basename of consensus map output path.
#' #' @param testRandom Should generate random data from input to test measurements. Default is `TRUE`.
#' #' @param nmfalg specification of the NMF algorithm.
#' #' @author Geoffrey Macintyre, Shixiang Wang
#' #' @importFrom grDevices pdf dev.off
#' #' @return a `list` contains information of NMF run and rank survey.
#' #' @export
#' #' @family CNV analysis functions
#' #' @examples
#' #' \dontrun{
#' #' ## load example copy-number data from tcga
#' #' load(system.file("inst/extdata", "example_cn_list.RData", package = "VSHunter"))
#' #' ## generate copy-number features
#' #' tcga_features = get_features(CN_data = tcga_segTabs, cores = 1, genome_build = "hg19")
#' #' ## fit mixture model  (this will take some time)
#' #' tcga_components = get_components(CN_features = tcga_features, cores = 1)
#' #' ## generate a sample-by-component matrix
#' #' tcga_sample_component_matrix = cnv_generateSbCMatrix(tcga_features, tcga_components, cores = 1)
#' #' ## optimal rank survey
#' #'  tcga_sig_choose = cnv_chooseSigNumber(tcga_sample_component_matrix,
#' #'  nrun = 10, cores = 1, plot = FALSE)
#' #' }
#' #'
#' cnv_chooseSigNumber <-
#'   function(nmf_matrix,
#'            nTry = 12,
#'            nrun = 10,
#'            cores = 1,
#'            seed = 123456,
#'            plot = TRUE,
#'            consensusmap_name = "nmf_consensus",
#'            testRandom = TRUE,
#'            nmfalg = "brunet")
#'   {
#'     message('Estimating best rank..')
#'     #nmfalg <- "brunet"
#'
#'     #suppressMessages(library(NMF))
#'     estim.r <-
#'       NMF::nmfEstimateRank(
#'         t(nmf_matrix),
#'         seq(2, nTry),
#'         seed = seed,
#'         nrun = nrun,
#'         verbose = TRUE,
#'         method = nmfalg,
#'         .opt = list(shared.memory = FALSE, paste0("p", cores))
#'       )
#'
#'     pdf(
#'       paste0(consensusmap_name, ".pdf"),
#'       bg = 'white',
#'       pointsize = 9,
#'       width = 12,
#'       height = 12,
#'       paper = "special"
#'     )
#'     NMF::consensusmap(estim.r)
#'     dev.off()
#'     message('created ', paste0(consensusmap_name, ".pdf"))
#'
#'     #--- copy from maftools and modified ---#
#'     nmf.sum = NMF::summary(estim.r) # Get summary of estimates
#'     print(nmf.sum)
#'     nmf.sum$diff = c(0, diff(nmf.sum$cophenetic))
#'     bestFit = nmf.sum$rank[which(nmf.sum$diff < 0)][1]
#'     #bestFit = nmf.sum[diff < 0, rank][1] #First point where cophenetic correlation coefficient starts decreasing
#'     if (is.na(bestFit) | is.null(bestFit)) {
#'       stop("Cannot find appropriate rank, please reset your arguments.")
#'     }
#'
#'     # https://blog.csdn.net/YJJ18636810884/article/details/83214566
#'     # https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-367
#'     message(
#'       paste(
#'         'Using ',
#'         bestFit,
#'         ' as a best-fit rank based on decreasing cophenetic correlation coefficient.',
#'         sep = ''
#'       )
#'     )
#'     n = as.numeric(bestFit)
#'
#'     if (testRandom) {
#'       message("Generating random matrix and run NMF...")
#'       V.random <- NMF::randomize(t(nmf_matrix))
#'       estim.r.random <-
#'         NMF::nmfEstimateRank(
#'           V.random,
#'           seq(2, nTry),
#'           seed = seed,
#'           nrun = nrun,
#'           verbose = TRUE,
#'           method = nmfalg,
#'           .opt = list(shared.memory = FALSE, paste0("p", cores))
#'         )
#'     }
#'
#'     if (plot) {
#'       message('Creating nmf rank survey plot...')
#'
#'       if (testRandom) {
#'         p <- NMF::plot(
#'           estim.r,
#'           estim.r.random,
#'           what = c("cophenetic",
#'                    "dispersion",
#'                    "sparseness",
#'                    #"silhouette",
#'                    #"residuals",
#'                    "rss"),
#'           xname = "Observed",
#'           yname = "Randomised",
#'           main = "NMF Rank Survey"
#'         )
#'       } else {
#'         p <- NMF::plot(
#'           estim.r,
#'           what = c("cophenetic",
#'                    "dispersion",
#'                    "sparseness",
#'                    # "silhouette",
#'                    # "residuals",
#'                    "rss"),
#'           main = "NMF Rank Survey"
#'         )
#'       }
#'
#'       print(p)
#'
#'     }
#'
#'     if (!plot)
#'       p = NULL
#'     if (!testRandom)
#'       estim.r.random = NULL
#'
#'     return(
#'       list(
#'         nmfEstimate = estim.r,
#'         nmfEstimate.random = estim.r.random,
#'         bestRank = n,
#'         survey = nmf.sum,
#'         survey_plot = p,
#'         seed = seed
#'       )
#'     )
#'   }
#'
#' #--------------------------
#' # extract signatures
#' #' @title Extract signature based on specified rank value
#' #' @inheritParams cnv_chooseSigNumber
#' #' @param nsig specification of the factorization rank.
#' #' @author Geoffrey Macintyre, Shixiang Wang
#' #' @return a object of \code{NMF} run.
#' #' @export
#' #' @aliases signatures
#' #' @keywords signature
#' #' @family CNV analysis functions
#' #' @seealso [cnv_plotSignatures()] for plot signatures and their contributions.
#' #' @examples
#' #' \dontrun{
#' #' ## load example copy-number data from tcga
#' #' load(system.file("inst/extdata", "example_cn_list.RData", package = "VSHunter"))
#' #' ## generate copy-number features
#' #' tcga_features = get_features(CN_data = tcga_segTabs, cores = 1, genome_build = "hg19")
#' #' ## fit mixture model  (this will take some time)
#' #' tcga_components = get_components(CN_features = tcga_features, cores = 1)
#' #' ## generate a sample-by-component matrix
#' #' tcga_sample_component_matrix = cnv_generateSbCMatrix(tcga_features, tcga_components, cores = 1)
#' #' ## optimal rank survey
#' #'  tcga_sig_choose = cnv_chooseSigNumber(tcga_sample_component_matrix, nrun = 10,
#' #'  cores = 1, plot = FALSE)
#' #'  tcga_signatures = cnv_extractSignatures(tcga_sample_component_matrix, nsig = 3, cores = 1)
#' #' }
#' cnv_extractSignatures <-
#'   function(nmf_matrix,
#'            nsig,
#'            seed = 123456,
#'            nmfalg = "brunet",
#'            cores = 1)
#'   {
#'     message("Running NMF based on specified rank...")
#'     #suppressMessages(library(NMF))
#'     NMF::nmf(
#'       t(nmf_matrix),
#'       nsig,
#'       seed = seed,
#'       nrun = 1000,
#'       method = nmfalg,
#'       .opt = paste0("p", cores)
#'     )
#'   }
#'
#'
#' #-------------------------------------------------------------------------------------
#' #' Auto-capture signature and coresponding exposure
#' #' @description  this is a wrapper of \code{cnv_chooseSigNumber}, \code{cnv_extractSignatures}
#' #' and \code{cnv_quantifySigExposure} these three functions.
#' #'
#' #' @inheritParams cnv_chooseSigNumber
#' #' @author Geoffrey Macintyre, Shixiang Wang
#' #' @return a `list` contains results of NMF best rank survey, run, signature matrix, exposure list etc..
#' #' @export
#' #' @inherit cnv_extractSignatures seealso
#' #' @family CNV analysis functions
#' #' @examples
#' #' \dontrun{
#' #' ## load example copy-number data from tcga
#' #' load(system.file("inst/extdata", "example_cn_list.RData", package = "VSHunter"))
#' #' ## generate copy-number features
#' #' tcga_features = get_features(CN_data = tcga_segTabs, cores = 1, genome_build = "hg19")
#' #' ## fit mixture model  (this will take some time)
#' #' tcga_components = get_components(CN_features = tcga_features, cores = 1)
#' #' ## generate a sample-by-component matrix
#' #' tcga_sample_component_matrix = cnv_generateSbCMatrix(tcga_features, tcga_components, cores = 1)
#' #' ## optimal rank survey13
#' #' tcga_results = cnv_autoCaptureSignatures(tcga_sample_component_matrix, nrun=10, cores = 1)
#' #' }
#' cnv_autoCaptureSignatures = function(nmf_matrix,
#'                                      nTry = 12,
#'                                      nrun = 10,
#'                                      cores = 1,
#'                                      seed = 123456,
#'                                      plot = TRUE,
#'                                      consensusmap_name = "nmf_consensus",
#'                                      testRandom = TRUE,
#'                                      nmfalg = "brunet") {
#'   choose_res = cnv_chooseSigNumber(
#'     nmf_matrix,
#'     nTry,
#'     nrun,
#'     cores,
#'     seed,
#'     plot = plot,
#'     consensusmap_name = consensusmap_name,
#'     testRandom = testRandom,
#'     nmfalg = nmfalg
#'   )
#'   NMF_res = cnv_extractSignatures(nmf_matrix,
#'                                   nsig = choose_res$bestRank,
#'                                   cores = cores, nmfalg = nmfalg)
#'   #-- Signatures
#'   w = NMF::basis(NMF_res)
#'   w = apply(w, 2, function(x)
#'     x / sum(x)) # Scale signatures (basis)
#'   colnames(w) = paste("Signature", 1:ncol(w), sep = "_")
#'
#'   #-- Contributions
#'   h = NMF::coef(NMF_res)
#'   h = apply(h, 2, function(x)
#'     x / sum(x)) # Scale contributions (coefs)
#'
#'   rownames(h) = paste("Signature", 1:ncol(w), sep = "_")
#'   message("Done.")
#'
#'   return(
#'     list(
#'       NMF = NMF_res,
#'       signatures = w,
#'       contributions = h,
#'       nmfEstimate = choose_res$nmfEstimate,
#'       nmfEstimate.random = choose_res$nmfEstimate.random,
#'       bestRank = choose_res$bestRank,
#'       survey = choose_res$survey,
#'       survey_plot = choose_res$survey_plot,
#'       seed = choose_res$seed
#'     )
#'   )
#' }
