#' #------------------------------------
#' #' @title Choose optimal number of signatures
#' #' @description This function use **NMF** package to evaluate the optimal number of signatures.
#' #' The most common approach is to choose the smallest rank for which cophenetic correlation coefficient
#' #' starts decreasing (Used by this function). Another approach is to choose the rank for which the plot
#' #' of the residual sum of squares (RSS) between the input matrix and its estimate shows an inflection point.
#' #' @param sample_by_component a sample-by-component `matrix`, generate from [cnv_generateSbCMatrix] function.
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
#'   function(sample_by_component,
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
#'         t(sample_by_component),
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
#'       V.random <- NMF::randomize(t(sample_by_component))
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
#'   function(sample_by_component,
#'            nsig,
#'            seed = 123456,
#'            nmfalg = "brunet",
#'            cores = 1)
#'   {
#'     message("Running NMF based on specified rank...")
#'     #suppressMessages(library(NMF))
#'     NMF::nmf(
#'       t(sample_by_component),
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
#' cnv_autoCaptureSignatures = function(sample_by_component,
#'                                      nTry = 12,
#'                                      nrun = 10,
#'                                      cores = 1,
#'                                      seed = 123456,
#'                                      plot = TRUE,
#'                                      consensusmap_name = "nmf_consensus",
#'                                      testRandom = TRUE,
#'                                      nmfalg = "brunet") {
#'   choose_res = cnv_chooseSigNumber(
#'     sample_by_component,
#'     nTry,
#'     nrun,
#'     cores,
#'     seed,
#'     plot = plot,
#'     consensusmap_name = consensusmap_name,
#'     testRandom = testRandom,
#'     nmfalg = nmfalg
#'   )
#'   NMF_res = cnv_extractSignatures(sample_by_component,
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
