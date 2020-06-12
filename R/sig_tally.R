#' Tally a Genomic Alteration Object
#'
#' Tally a variation object like [MAF], [CopyNumber] and return a matrix for NMF de-composition and more.
#' This is a generic function,
#' so it can be further extended to other mutation cases.
#' Please read details about how to set sex for identifying copy number signatures.
#' Please read <https://osf.io/s93d5/> for the generation of SBS, DBS and ID (INDEL)
#' components. **Of note, many options are designed for method "M" only, and they are highlighted
#' by bold fonts** (you can ignore them if you don't use "M" method).
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
#' - If your cohort contains both males and females, set `sigminer.sex`
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
#' \donttest{
#' # Prepare SBS signature analysis
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read_maf(maf = laml.maf)
#' if (require("BSgenome.Hsapiens.UCSC.hg19")) {
#'   mt_tally <- sig_tally(
#'     laml,
#'     ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
#'     use_syn = TRUE
#'   )
#'   mt_tally$nmf_matrix[1:5, 1:5]
#'
#'   ## Use strand bias categories
#'   mt_tally <- sig_tally(
#'     laml,
#'     ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
#'     use_syn = TRUE, add_trans_bias = TRUE
#'   )
#'   ## Test it by enrichment analysis
#'   enrich_component_strand_bias(mt_tally$nmf_matrix)
#'   enrich_component_strand_bias(mt_tally$all_matrices$SBS_24)
#' } else {
#'   message("Please install package 'BSgenome.Hsapiens.UCSC.hg19' firstly!")
#' }
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
#'     use_syn = TRUE
#'   )
#'
#'   expect_equal(length(mt_tally), 3L)
#'
#'   ## Use strand bias categories
#'   mt_tally <- sig_tally(
#'     laml,
#'     ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
#'     use_syn = TRUE, add_trans_bias = TRUE
#'   )
#'   ## Test it by enrichment analysis
#'   dt1 = enrich_component_strand_bias(mt_tally$nmf_matrix)
#'   dt2 = enrich_component_strand_bias(mt_tally$all_matrices$SBS_24)
#'
#'   expect_s3_class(dt1, "data.table")
#'   expect_s3_class(dt2, "data.table")
#' } else {
#'   message("Please install package 'BSgenome.Hsapiens.UCSC.hg19' firstly!")
#' }
#'
sig_tally <- function(object, ...) {
  timer <- Sys.time()
  send_info("Started.")
  on.exit(send_elapsed_time(timer))

  UseMethod("sig_tally")
}

#' @describeIn sig_tally Returns copy number features, components and component-by-sample matrix
#' @param method method for feature classfication, can be one of "Macintyre" ("M"),
#' "Wang" ("W") and "Tao & Wang" ("T").
#' @param feature_setting a `data.frame` used for classification.
#' **Only used when method is "Wang" ("W")**.
#' Default is [CN.features]. Users can also set custom input with "feature",
#' "min" and "max" columns available. Valid features can be printed by
#' `unique(CN.features$feature)`.
#' @param type one of "probability", "count". Default is "probability", return a matrix
#' with the sum of posterior probabilities for each components. If set to 'count',
#' return a matrix with event count assigned to each components. The result for
#' both types should be close. **Only used when method is "Macintyre"**.
#' @param reference_components default is `FALSE`, calculate mixture components
#' from [CopyNumber] object. **Only used when method is "Macintyre"**.
#' @param cores number of computer cores to run this task.
#' You can use [future::availableCores()] function to check how
#' many cores you can use.
#' @param seed seed number. **Only used when method is "Macintyre"**.
#' @param min_comp minimal number of components to fit, default is 2.
#' Can also be a vector with length 6, which apply to each feature.
#' **Only used when method is "Macintyre"**.
#' @param max_comp maximal number of components to fit, default is 15.
#' Can also be a vector with length 6, which apply to each feature.
#' **Only used when method is "Macintyre"**.
#' @param min_prior the minimum relative size of components, default is 0.001.
#' Details about custom setting please refer to **flexmix** package.
#' **Only used when method is "Macintyre"**.
#' @param model_selection model selection strategy, default is 'BIC'.
#' Details about custom setting please refer to **flexmix** package.
#' **Only used when method is "Macintyre"**.
#' @param threshold default is `0.1`. Sometimes, the result components
#' include adjacent distributions with similar mu
#' (two and more distribution are very close), we use this threshold
#' to obtain a more meaningful fit with less components.
#' **Only used when method is "Macintyre"**.
#' @param nrep number of run times for each value of component,
#' keep only the solution with maximum likelihood.
#' **Only used when method is "Macintyre"**.
#' @param niter the maximum number of iterations.
#' **Only used when method is "Macintyre"**.
#' @param keep_only_matrix if `TRUE`, keep only matrix for signature extraction.
#' For a `MAF` object, this will just return the most useful matrix.
#' @references
#' Macintyre, Geoff, et al. "Copy number signatures and mutational
#' processes in ovarian carcinoma." Nature genetics 50.9 (2018): 1262.
#'
#' Wang, Shixiang, et al. "Copy number signature analyses in prostate cancer reveal
#' distinct etiologies and clinical outcomes." medRxiv (2020).
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
  method <- match.arg(method, choices = c("Macintyre", "M", "Wang", "W", "Tao & Wang", "T"))

  if (startsWith(method, "T")) {
    ## Add segment index for method "T" so the segments can be easily joined or checked
    cn_list <- get_cnlist(object, ignore_chrs = ignore_chrs, add_index = TRUE)
  } else {
    cn_list <- get_cnlist(object, ignore_chrs = ignore_chrs)
  }

  if (startsWith(method, "M")) {
    if (!requireNamespace("flexmix", quietly = TRUE)) {
      send_stop("Please install 'flexmix' package firstly.")
    }

    # Method: Macintyre
    type <- match.arg(type)

    send_info("Step: getting copy number features.")
    cn_features <- get_features(
      CN_data = cn_list, cores = cores,
      genome_build = object@genome_build
    )
    cn_features <- lapply(cn_features, function(x) as.data.frame(x))
    send_success("Gotten.")

    send_info("Step: fitting copy number components.")
    if (is.null(reference_components) | is.list(reference_components)) {
      send_success("Reference components detected.")
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
      send_success("Components fitted.")
    }

    if (type == "count") {
      send_info("Step: calculating the sum of cluster count based on posterior probabilities.")
      cn_matrix <- get_matrix(cn_features, cn_components,
        type = "count",
        cores = cores
      )
    } else {
      send_info("Step: calculating the sum of posterior probabilities.")
      cn_matrix <- get_matrix(cn_features, cn_components,
        type = "probability",
        cores = cores
      )
    }
  } else if (startsWith(method, "W")) {
    # Method: Wang Shixiang, 'W'

    send_info("Step: getting copy number features.")
    cn_features <- get_features_wang(
      CN_data = cn_list, cores = cores,
      genome_build = object@genome_build,
      feature_setting = feature_setting
    )
    send_success("Gotten.")
    # Make order as unique(feature_setting)$feature
    # cn_features <- cn_features[unique(feature_setting$feature)]

    send_info("Step: generating copy number components.")
    # Check feature setting
    if (!inherits(feature_setting, "sigminer.features")) {
      feature_setting <- get_feature_components(feature_setting)
    }
    send_success("{.code feature_setting} checked.")

    send_info("Step: counting components.")
    cn_components <- purrr::map2(cn_features, names(cn_features),
      count_components_wrapper,
      feature_setting = feature_setting
    )
    send_success("Counted.")

    ## Remove BoChr value is 0 in features
    if ("BoChr" %in% names(cn_features)) {
      cn_features$BoChr <- cn_features$BoChr[cn_features$BoChr$value != 0]
    }

    send_info("Step: generating components by sample matrix.")
    cn_matrix <- data.table::rbindlist(cn_components, fill = TRUE, use.names = TRUE) %>%
      dplyr::as_tibble() %>%
      tibble::column_to_rownames(var = "component") %>%
      as.matrix()
    # Order the matrix as feature_setting
    cn_matrix <- cn_matrix[feature_setting$component, ] %>%
      t()

    if (any(is.na(cn_matrix))) {
      send_warning("{.code NA} detected. There may be an issue, please contact the developer!")
      send_warning("Data will still returned, but please take case of it.")
    }
    # cn_matrix[is.na(cn_matrix)] <- 0L
    feature_setting$n_obs <- colSums(cn_matrix, na.rm = TRUE)
  } else {
    # Method: Tao & Wang, 'T'

    ## TODO: should also set a feature_setting dataset??
    ## Is genome_build necessary here??
    send_info("Step: getting copy number features.")
    cn_features <- get_features_mutex(
      CN_data = cn_list, cores = cores,
      genome_build = object@genome_build,
      feature_setting = feature_setting
    )
    send_success("Gotten.")

    ## Curretly return features and explore how to combine them
    return(cn_features)
    # Make order as unique(feature_setting)$feature
    # cn_features <- cn_features[unique(feature_setting$feature)]

    # send_info("Step: generating copy number components based on combination.")
    # # Check feature setting
    # if (!inherits(feature_setting, "sigminer.features")) {
    #   feature_setting <- get_feature_components(feature_setting)
    # }
    # send_success("{.code feature_setting} checked.")
    #
    # send_info("Step: counting components.")
    # cn_components <- purrr::map2(cn_features, names(cn_features),
    #                              count_components_wrapper,
    #                              feature_setting = feature_setting
    # )
    # send_success("Counted.")
    #
    # ## Remove BoChr value is 0 in features
    # if ("BoChr" %in% names(cn_features)) {
    #   cn_features$BoChr <- cn_features$BoChr[cn_features$BoChr$value != 0]
    # }
    #
    # send_info("Step: generating components by sample matrix.")
    # cn_matrix <- data.table::rbindlist(cn_components, fill = TRUE, use.names = TRUE) %>%
    #   dplyr::as_tibble() %>%
    #   tibble::column_to_rownames(var = "component") %>%
    #   as.matrix()
    # # Order the matrix as feature_setting
    # cn_matrix <- cn_matrix[feature_setting$component, ] %>%
    #   t()
    #
    # if (any(is.na(cn_matrix))) {
    #   send_warning("{.code NA} detected. There may be an issue, please contact the developer!")
    #   send_warning("Data will still returned, but please take case of it.")
    # }
    # # cn_matrix[is.na(cn_matrix)] <- 0L
    # feature_setting$n_obs <- colSums(cn_matrix, na.rm = TRUE)
  }

  send_success("Matrix generated.")

  if (keep_only_matrix) {
    cn_matrix
  } else {
    if (startsWith(method, "M")) {
      if (is.null(cn_components)) {
        cn_components <- readRDS(file.path(tempdir(), "Nat_Gen_component_parameters.rds"))
      }
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

#' @describeIn sig_tally Returns SBS mutation sample-by-component matrix and APOBEC enrichment
#' @inheritParams maftools::trinucleotideMatrix
#' @param mode type of mutation matrix to extract, can be one of 'SBS', 'DBS' and 'ID'.
#' @param genome_build genome build 'hg19' or 'hg38', if not set, guess it by `ref_genome`.
#' @param add_trans_bias if `TRUE`, consider transcriptional bias categories.
#' 'T:' for Transcribed (the variant is on the transcribed strand);
#' 'U:' for Un-transcribed (the variant is on the untranscribed strand);
#' 'B:' for Bi-directional (the variant is on both strand and is transcribed either way);
#' 'N:' for Non-transcribed (the variant is in a non-coding region and is untranslated);
#' 'Q:' for Questionable.
#' **NOTE**: the result counts of 'B' and 'N' labels are a little different from
#' SigProfilerMatrixGenerator, the reason is unknown (may be caused by annotation file).
#' @param ignore_chrs Chromsomes to ignore from analysis. e.g. chrX and chrY.
#' @param use_syn Logical. Whether to include synonymous variants in analysis. Defaults to TRUE
#' @references Mayakonda, Anand, et al. "Maftools: efficient and comprehensive analysis of somatic variants in cancer." Genome research 28.11 (2018): 1747-1756.
#' @references Roberts SA, Lawrence MS, Klimczak LJ, et al. An APOBEC Cytidine Deaminase Mutagenesis Pattern is Widespread in Human Cancers. Nature genetics. 2013;45(9):970-976. doi:10.1038/ng.2702.
#' @references Bergstrom EN, Huang MN, Mahto U, Barnes M, Stratton MR, Rozen SG, Alexandrov LB: SigProfilerMatrixGenerator: a tool for visualizing and exploring patterns of small mutational events. BMC Genomics 2019, 20:685 https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-6041-2
#' @export
sig_tally.MAF <- function(object, mode = c("SBS", "DBS", "ID", "ALL"),
                          ref_genome = NULL,
                          genome_build = NULL,
                          add_trans_bias = FALSE,
                          ignore_chrs = NULL,
                          use_syn = TRUE,
                          keep_only_matrix = FALSE,
                          ...) {
  if (!requireNamespace("BSgenome", quietly = TRUE)) {
    send_stop("Please install 'BSgenome' package firstly.")
  }

  mode <- match.arg(mode)

  if (is.null(genome_build)) {
    if (grepl("hg19", ref_genome)) {
      genome_build <- "hg19"
    } else if (grepl("hg38", ref_genome)) {
      genome_build <- "hg38"
    } else {
      send_stop("Cannot guess the genome build, please set it by hand!")
    }
  }

  hsgs.installed <- BSgenome::installed.genomes(splitNameParts = TRUE)
  data.table::setDT(x = hsgs.installed)
  # hsgs.installed = hsgs.installed[organism %in% "Hsapiens"]

  if (nrow(hsgs.installed) == 0) {
    send_stop("Could not find any installed BSgenomes. Use {.code BSgenome::available.genomes()} for options.")
  }

  if (is.null(ref_genome)) {
    send_info("User did not set {.code ref_genome}.")
    send_success("Found following BSgenome installtions. Using first entry.\n")
    print(hsgs.installed)
    ref_genome <- hsgs.installed$pkgname[1]
  } else {
    if (!ref_genome %in% hsgs.installed$pkgname) {
      send_error("Could not find BSgenome {.code ", ref_genome, "}.")
      send_info("Found following BSgenome installtions. Correct {.code ref_genome} argument if necessary.")
      print(hsgs.installed)
      send_stop("Exit.")
    }
  }

  ref_genome <- BSgenome::getBSgenome(genome = ref_genome)
  send_success("Reference genome loaded.")

  query <- maftools::subsetMaf(
    maf = object,
    query = "Variant_Type %in% c('SNP', 'INS', 'DEL')", fields = "Chromosome",
    includeSyn = use_syn, mafObj = FALSE
  )
  send_success("Variants from MAF object queried.")

  # Remove unwanted contigs
  if (!is.null(ignore_chrs)) {
    query <- query[!query$Chromosome %in% ignore_chrs]
    send_success("Unwanted contigs removed.")
  }

  if (nrow(query) == 0) {
    send_stop("Zero variants to analyze!")
  }

  query$Chromosome <- sub(
    pattern = "chr",
    replacement = "chr",
    x = as.character(query$Chromosome),
    ignore.case = TRUE
  )
  ## Make sure all have prefix
  if (any(!grepl("chr", query$Chromosome))) {
    query$Chromosome[!grepl("chr", query$Chromosome)] <-
      paste0("chr", query$Chromosome[!grepl("chr", query$Chromosome)])
  }

  send_success("Chromosome names checked.")

  ## Handle non-autosomes
  query$Chromosome <- sub(
    pattern = "x",
    replacement = "X",
    x = as.character(query$Chromosome),
    ignore.case = TRUE
  )

  query$Chromosome <- sub(
    pattern = "y",
    replacement = "Y",
    x = as.character(query$Chromosome),
    ignore.case = TRUE
  )

  query$Chromosome <- sub(
    pattern = "MT",
    replacement = "M",
    x = as.character(query$Chromosome),
    ignore.case = TRUE
  )

  # detect and transform chromosome 23 to "X"
  query$Chromosome <- sub("23", "X", query$Chromosome)
  # detect and transform chromosome 24 to "Y"
  query$Chromosome <- sub("24", "Y", query$Chromosome)

  send_success("Sex chromosomes properly handled.")

  query$Start_Position <- as.numeric(as.character(query$Start_Position))
  query$End_Position <- as.numeric(as.character(query$End_Position))

  send_success("Variant start and end position checked.")

  query_seq_lvls <- query[, .N, Chromosome]
  ref_seqs_lvls <- BSgenome::seqnames(x = ref_genome)
  query_seq_lvls_missing <- query_seq_lvls[!Chromosome %in% ref_seqs_lvls]

  if (nrow(query_seq_lvls_missing) > 0) {
    send_warning(paste0(
      "Chromosome names in MAF must match chromosome names in reference genome.\nIgnorinig ",
      query_seq_lvls_missing[, sum(N)],
      " single nucleotide variants from missing chromosomes ",
      paste(query_seq_lvls_missing[, Chromosome], collapse = ", ")
    ))
  }

  query <- query[!Chromosome %in% query_seq_lvls_missing[, Chromosome]]

  send_success("Variant data for matrix generation preprocessed.")

  if (mode == "SBS") {
    res <- generate_matrix_SBS(query, ref_genome, genome_build = genome_build, add_trans_bias = add_trans_bias)
  } else if (mode == "DBS") {
    res <- generate_matrix_DBS(query, ref_genome, genome_build = genome_build, add_trans_bias = add_trans_bias)
  } else if (mode == "ID") {
    ## INDEL
    res <- generate_matrix_INDEL(query, ref_genome, genome_build = genome_build, add_trans_bias = add_trans_bias)
  } else {
    send_info("All types of matrices generation - start.")

    res_SBS <- tryCatch(
      generate_matrix_SBS(query, ref_genome, genome_build = genome_build, add_trans_bias = add_trans_bias),
      error = function(e) {
        NULL
      }
    )
    res_DBS <- tryCatch(
      generate_matrix_DBS(query, ref_genome, genome_build = genome_build, add_trans_bias = add_trans_bias),
      error = function(e) {
        NULL
      }
    )
    res_ID <- tryCatch(
      generate_matrix_INDEL(query, ref_genome, genome_build = genome_build, add_trans_bias = add_trans_bias),
      error = function(e) {
        NULL
      }
    )

    send_info("All types of matrices generation (APOBEC scores included) - end.")
    res <- c(res_SBS$all_matrices, res_DBS$all_matrices, res_ID$all_matrices)
    res$APOBEC_scores <- res_SBS$APOBEC_scores
  }

  send_success("Done.")

  if (keep_only_matrix) {
    if (mode == "ALL") {
      send_warning("Mode 'ALL' cannot return a single matrix.")
      return(res)
    }
    return(res$nmf_matrix)
  } else {
    return(res)
  }
}


utils::globalVariables(
  c(
    ".N",
    "Chromosome"
  )
)
