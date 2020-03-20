#' Tally a Genomic Alteration Object
#'
#' Tally a variation object like [MAF] and return a matrix for NMF de-composition and more.
#' This is a generic function,
#' so it can be further extended to other mutation cases. Please read details
#' about how to set sex for identifying copy number signatures.
#'
#' For identifying copy number signatures, we have to derive copy number
#' features firstly. Due to the difference of copy number values in sex chromosomes
#' between male and female, we have to do an extra step **if we don't want to
#' ignore them**.
#'
#' I create two options to control this, the default values are shown as
#' the following, you can use the same way to set (per R session).
#'
#' `options(sigminer.sex = "female", sigminer.copynumber.max = NA_integer_)`
#'
#' - If your cohort are all females, you can totally ignore this.
#' - If your cohort are all males, set `sigminer.sex` to 'male' and
#' `sigminer.copynumber.max` to a proper value (the best is consistent
#' with [read_copynumber]).
#' - If your cohort contains bother males and females, set `sigminer.sex`
#' as a `data.frame` with two columns "sample" and "sex". And
#' set `sigminer.copynumber.max` to a proper value (the best is consistent
#' with [read_copynumber]).
#'
#' @param object a [CopyNumber] object or [MAF] object.
#' @param ... custom setting for operating object. Detail see S3 method for
#' corresponding class (e.g. `CopyNumber`).
#' @return a `list` contains a `matrix` used for NMF de-composition.
#' @author Shixiang Wang
#' @export
#' @seealso [sig_estimate] for estimating signature number for [sig_extract],
#' [sig_auto_extract] for extracting signatures using automatic relevance determination technique.
#' @examples
#' # Load copy number object
#' load(system.file("extdata", "toy_copynumber.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' \donttest{
#' # Use method designed by Wang, Shixiang et al.
#' cn_tally_W <- sig_tally(cn, method = "W")
#' # Use method designed by Macintyre et al.
#' cn_tally_M <- sig_tally(cn, method = "M")
#' }
#' @tests
#' ## Load copy number object
#' load(system.file("extdata", "toy_copynumber.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Use method designed by Macintyre et al.
#' cn_tally_M <- sig_tally(cn, method = "M")
#' # Use method designed by Wang, Shixiang et al.
#' cn_tally_W <- sig_tally(cn, method = "W")
#'
#' expect_equal(length(cn_tally_M), length(cn_tally_W))
#'
#' ## for SBS
#'
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read_maf(maf = laml.maf)
#' if (require("BSgenome.Hsapiens.UCSC.hg19")) {
#'   mt_tally <- sig_tally(
#'     laml,
#'     ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
#'     prefix = "chr", add = TRUE, useSyn = TRUE
#'   )
#'
#'   expect_equal(length(mt_tally), 2L)
#' } else {
#'   message("Please install package 'BSgenome.Hsapiens.UCSC.hg19' firstly!")
#' }
#'
sig_tally <- function(object, ...) {
  UseMethod("sig_tally")
}

#' @describeIn sig_tally Returns copy number features, components and component-by-sample matrix
#' @param method method for feature classfication, can be one of "Macintyre" ("M") and
#' "Wang" ("W").
#' @param feature_setting a `data.frame` used for classification.
#' Only used when method is "Wang" ("W").
#' Default is [CN.features]. Users can also set custom input with "feature",
#' "min" and "max" columns available. Valid features can be printed by
#' `unique(CN.features$feature)`.
#' @param type one of "probability", "count". Default is "probability", return a matrix
#' with the sum of posterior probabilities for each components. If set to 'count',
#' return a matrix with event count assigned to each components. The result for
#' both types should be close. Only used when method is "Macintyre".
#' @param reference_components default is `FALSE`, calculate mixture components
#' from [CopyNumber] object. Only used when method is "Macintyre".
#' @param cores number of compute cores to run this task.
#' You can use [future::availableCores()] function to check how
#' many cores you can use.
#' @param seed seed number. Only used when method is "Macintyre".
#' @param min_comp minimal number of components to fit, default is 2.
#' Can also be a vector with length 6, which apply to each feature.
#' Only used when method is "Macintyre".
#' @param max_comp maximal number of components to fit, default is 15.
#' Can also be a vector with length 6, which apply to each feature.
#' Only used when method is "Macintyre".
#' @param min_prior the minimum relative size of components, default is 0.001.
#' Details about custom setting please refer to **flexmix** package.
#' Only used when method is "Macintyre".
#' @param model_selection model selection strategy, default is 'BIC'.
#' Details about custom setting please refer to **flexmix** package.
#' Only used when method is "Macintyre".
#' @param threshold default is `0.1`. Sometimes, the result components
#' include adjacent distributions with similar mu
#' (two and more distribution are very close), we use this threshold
#' to obtain a more meaningful fit with less components.
#' Only used when method is "Macintyre".
#' @param nrep number of run times for each value of component,
#' keep only the solution with maximum likelihood.
#' Only used when method is "Macintyre".
#' @param niter the maximum number of iterations.
#' Only used when method is "Macintyre".
#' @param keep_only_matrix if `TRUE`, keep only matrix for signature extraction.
#' @references Macintyre, Geoff, et al. "Copy number signatures and mutational
#' processes in ovarian carcinoma." Nature genetics 50.9 (2018): 1262.
#' @export
sig_tally.CopyNumber <- function(object,
                                 method = "Wang",
                                 ignore_chrs = NULL,
                                 feature_setting = sigminer::CN.features,
                                 type = c("probability", "count"),
                                 reference_components = FALSE,
                                 cores = 1, seed = 123456,
                                 min_comp = 2, max_comp = 15,
                                 min_prior = 0.001,
                                 model_selection = "BIC",
                                 threshold = 0.1,
                                 nrep = 1, niter = 1000,
                                 keep_only_matrix = FALSE,
                                 ...) {
  stopifnot(is.logical(reference_components) | is.list(reference_components) | is.null(reference_components))
  method <- match.arg(method, choices = c("Macintyre", "M", "Wang", "W"))

  cn_list <- get_cnlist(object, ignore_chrs = ignore_chrs)

  if (startsWith(method, "M")) {
    # Method: Macintyre
    type <- match.arg(type)

    message("=> Step: getting copy number features")
    cn_features <- get_features(
      CN_data = cn_list, cores = cores,
      genome_build = object@genome_build
    )
    cn_features <- lapply(cn_features, function(x) as.data.frame(x))

    message("=> Step: fitting copy number components")
    if (is.null(reference_components) | is.list(reference_components)) {
      message("Detected reference components.")
      cn_components <- reference_components
    } else {
      cn_components <- get_components(
        CN_features = cn_features, seed = seed,
        min_comp = min_comp, max_comp = max_comp,
        min_prior = min_prior,
        model_selection = model_selection,
        threshold = threshold,
        nrep = nrep, niter = niter, cores = cores
      )
    }

    if (type == "count") {
      message("=> Step: calculating the sum of posterior probabilities")
      cn_matrix <- get_matrix(cn_features, cn_components,
        type = "count",
        cores = cores
      )
    } else {
      message("=> Step: calculating the sum of posterior probabilities")
      cn_matrix <- get_matrix(cn_features, cn_components,
        type = "probability",
        cores = cores
      )
    }
  } else {
    # Method: Wang Shixiang
    message("=> Step: getting copy number features")
    cn_features <- get_features_wang(
      CN_data = cn_list, cores = cores,
      genome_build = object@genome_build,
      feature_setting = feature_setting
    )
    # Make order as unique(feature_setting)$feature
    # cn_features <- cn_features[unique(feature_setting$feature)]

    message("=> Step: generating copy number components")
    # Chck feature setting
    if (!inherits(feature_setting, "sigminer.features")) {
      feature_setting <- get_feature_components(feature_setting)
    }

    cn_components <- purrr::map2(cn_features, names(cn_features),
      count_components_wrapper,
      feature_setting = feature_setting
    )

    ## Remove BoChr value is 0 in features
    if ("BoChr" %in% names(cn_features)) {
      cn_features$BoChr <- cn_features$BoChr[cn_features$BoChr$value != 0]
    }

    message("=> Step: generating components by sample matrix")
    cn_matrix <- data.table::rbindlist(cn_components, fill = TRUE, use.names = TRUE) %>%
      dplyr::as_tibble() %>%
      tibble::column_to_rownames(var = "component") %>%
      as.matrix()
    # Order the matrix as feature_setting
    cn_matrix <- cn_matrix[feature_setting$component, ] %>%
      t()

    if (any(is.na(cn_matrix))) {
      message("Warning: NA detected. There may be an issue, please contact the developer!")
      message("\tData will still returned, but please take case of it.")
    }
    # cn_matrix[is.na(cn_matrix)] <- 0L
    feature_setting$n_obs <- colSums(cn_matrix, na.rm = TRUE)
  }

  message("=> Done.")
  if (keep_only_matrix) {
    cn_matrix
  } else {
    if (startsWith(method, "M")) {
      para_df <- get_tidy_parameter(cn_components)
    } else {
      para_df <- feature_setting
    }

    list(
      features = cn_features,
      components = cn_components,
      parameters = para_df,
      nmf_matrix = cn_matrix
    )
  }
}

#' @describeIn sig_tally Returns SBS mutation component-by-sample matrix and APOBEC enrichment
#' @inheritParams maftools::trinucleotideMatrix
#' @param ignore_chrs Chromsomes to ignore from analysis. e.g. chrX and chrY.
#' @param use_syn Logical. Whether to include synonymous variants in analysis. Defaults to TRUE
#' @references Mayakonda, Anand, et al. "Maftools: efficient and comprehensive analysis of somatic variants in cancer." Genome research 28.11 (2018): 1747-1756.
#' @examples
#' \donttest{
#' # Prepare SBS signature analysis
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read_maf(maf = laml.maf)
#' if (require("BSgenome.Hsapiens.UCSC.hg19")) {
#'   mt_tally <- sig_tally(
#'     laml,
#'     ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
#'     prefix = "chr", add = TRUE, useSyn = TRUE
#'   )
#' } else {
#'   message("Please install package 'BSgenome.Hsapiens.UCSC.hg19' firstly!")
#' }
#' }
#' @export
sig_tally.MAF <- function(object, ref_genome = NULL, prefix = NULL,
                          add = TRUE, ignore_chrs = NULL, use_syn = TRUE,
                          keep_only_matrix = FALSE,
                          ...) {
  # TODO: Rewrite this function instead of using maftools
  # Make result consistent with result from sig_tally.CopyNumber
  res <- maftools::trinucleotideMatrix(
    object,
    ref_genome = ref_genome, prefix = prefix,
    add = add, ignoreChr = ignore_chrs, useSyn = use_syn, fn = NULL
  )

  if (keep_only_matrix) {
    res$nmf_matrix
  } else {
    res
  }
}
