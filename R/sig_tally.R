#' Tally a Genomic Alteration Object
#'
#' Tally a variation object like [MAF], [CopyNumber] and return a matrix for NMF de-composition and more.
#' This is a generic function,
#' so it can be further extended to other mutation cases.
#' **Please read details about how to set sex for identifying copy number signatures**.
#' Please read <https://osf.io/s93d5/> for the generation of SBS, DBS and ID (INDEL)
#' components.
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
#' @param object a [CopyNumber] object or [MAF] object or SV object (from [read_sv_as_rs]).
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
#' }
#' # Use method designed by Steele et al.
#' # See example in read_copynumber
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
#' # Use method designed by Wang, Shixiang et al.
#' cn_tally_W <- sig_tally(cn, method = "W")
#' # Use method designed by Tao & Wang.
#' cn_tally_T <- sig_tally(cn, method = "T")
#'
#' expect_equal(length(cn_tally_T), 5L)
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
#' @param indices integer vector indicating segments to keep.
#' @param method method for feature classification, can be one of
#' "Wang" ("W"), "S" (for method described in Steele et al. 2019).
#' @param add_loh flag to add LOH classifications.
#' @param feature_setting a `data.frame` used for classification.
#' **Only used when method is "Wang" ("W")**.
#' Default is [CN.features]. Users can also set custom input with "feature",
#' "min" and "max" columns available. Valid features can be printed by
#' `unique(CN.features$feature)`.
#' @param cores number of computer cores to run this task.
#' You can use [future::availableCores()] function to check how
#' many cores you can use.
#' @param keep_only_matrix if `TRUE`, keep only matrix for signature extraction.
#' For a `MAF` object, this will just return the most useful matrix.
#' @references
#' Wang, Shixiang, et al. "Copy number signature analyses in prostate cancer reveal
#' distinct etiologies and clinical outcomes." medRxiv (2020).
#'
#' Steele, Christopher D., et al. "Undifferentiated sarcomas develop through
#' distinct evolutionary pathways." Cancer Cell 35.3 (2019): 441-456.
#' @export
sig_tally.CopyNumber <- function(object,
                                 method = "Wang",
                                 ignore_chrs = NULL,
                                 indices = NULL,
                                 add_loh = FALSE,
                                 feature_setting = sigminer::CN.features,
                                 cores = 1,
                                 keep_only_matrix = FALSE,
                                 ...) {
  method <- match.arg(method, choices = c("Wang", "W", "Tao & Wang", "T", "X", "S"))

  if (startsWith(method, "T") | method == "X") {
    send_warning("Currently, the method 'T' is in experimental stage, please don't use it for now!")
    ## Add segment index for method "T" so the segments can be easily joined or checked
    cn_list <- get_cnlist(object, ignore_chrs = ignore_chrs, add_index = TRUE)
  } else {
    cn_list <- get_cnlist(object, ignore_chrs = ignore_chrs)
  }

  if (startsWith(method, "W")) {
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
  } else if (startsWith(method, "S")) {
    # When use method "S", join_adj_seg should set to FALSE in read_copynumber
    send_info("When you use method 'S', please make sure you have set 'join_adj_seg' to FALSE and 'add_loh' to TRUE in 'read_copynumber() in the previous step!")
    mat_list <- get_matrix_mutex_sv(data.table::rbindlist(cn_list, idcol = "sample"))
    cn_features <- NULL
    cn_components <- mat_list$data
    cn_matrix <- mat_list$CN_40
  } else {
    # Method: Shixiang Wang, Ziyu Tao and Tao Wu, short with 'T'
    send_info("Step: getting copy number features.")
    cn_features <- get_features_mutex(
      CN_data = cn_list,
      add_loh = add_loh,
      # 'X' for final version
      XVersion = method == "X",
      cores = cores
    )
    send_success("Gotten.")

    send_info("Step: generating copy number components based on combination.")
    cn_components <- get_components_mutex(cn_features, XVersion = method == "X")
    send_success("Classified and combined.")

    send_info("Step: generating components by sample matrix.")
    if (method != "X") {
      cn_matrix_list <- get_matrix_mutex(cn_components,
        indices = indices
      )
    } else {
      cn_matrix_list <- get_matrix_mutex_xv(cn_components,
        indices = indices
      )
    }

    cn_matrix <- cn_matrix_list$ss_mat

    if (keep_only_matrix) {
      send_info("When keep_only_matrix is TRUE, only standard matrix kept.")
    }
  }

  send_success("Matrix generated.")

  if (keep_only_matrix) {
    cn_matrix
  } else {
    if (startsWith(method, "W")) {
      para_df <- feature_setting
    } else if (startsWith(method, "T")) {
      para_df <- "Message: No this info for method T."
    } else if (startsWith(method, "X")) {
      para_df <- "Message: No this info for method X."
    } else if (startsWith(method, "S")) {
      para_df <- "Message: No this info for method S."
    }

    if (startsWith(method, "T") | method == "X" | method == "S") {
      res_list <- list(
        features = cn_features,
        components = cn_components,
        parameters = para_df,
        nmf_matrix = cn_matrix,
        all_matrices = if (method == "X") {
          list(
            simplified_matrix = cn_matrix_list$ss_mat,
            standard_matrix = cn_matrix_list$s_mat
          )
        } else if (method == "S") {
          list(
            CN_40 = mat_list$CN_40,
            CN_48 = mat_list$CN_48
          )
        } else {
          list(
            simplified_matrix = cn_matrix_list$ss_mat,
            standard_matrix = cn_matrix_list$s_mat,
            complex_matrix = cn_matrix_list$c_mat
          )
        }
      )
    } else {
      res_list <- list(
        features = cn_features,
        components = cn_components,
        parameters = para_df,
        nmf_matrix = cn_matrix
      )
    }

    return(res_list)
  }
}

#' @describeIn sig_tally Returns genome rearrangement sample-by-component matrix
#' @export
sig_tally.RS <- function(object, keep_only_matrix = FALSE, ...) {
  svlist <- get_svlist(object)
  send_success("Successfully get RS list!")

  sv_features <- get_features_sv(svlist)
  send_success("Successfully get RS features!")

  sv_component <- get_components_sv(sv_features)
  send_success("Successfully get RS component!")

  sv_matrix_list <- get_matrix_sv(CN_components = sv_component)
  send_success("Successfully get RS matrix!")

  res_list <- list(
    features = sv_features,
    components = sv_component,
    nmf_matrix = sv_matrix_list$RS_32,
    all_matrices = sv_matrix_list
  )

  if (keep_only_matrix) {
    return(res_list$nmf_matrix)
  } else {
    return(res_list)
  }
}


#' @describeIn sig_tally Returns SBS mutation sample-by-component matrix and APOBEC enrichment
#' @param mode type of mutation matrix to extract, can be one of 'SBS', 'DBS' and 'ID'.
#' @param ref_genome 'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38',
#' 'BSgenome.Mmusculus.UCSC.mm10',  'BSgenome.Mmusculus.UCSC.mm9', etc.
#' @param genome_build genome build 'hg19', 'hg38', 'mm9' or "mm10", if not set, guess it by `ref_genome`.
#' @param add_trans_bias if `TRUE`, consider transcriptional bias categories.
#' 'T:' for Transcribed (the variant is on the transcribed strand);
#' 'U:' for Un-transcribed (the variant is on the untranscribed strand);
#' 'B:' for Bi-directional (the variant is on both strand and is transcribed either way);
#' 'N:' for Non-transcribed (the variant is in a non-coding region and is untranslated);
#' 'Q:' for Questionable.
#' **NOTE**: the result counts of 'B' and 'N' labels are a little different from
#' SigProfilerMatrixGenerator, the reason is unknown (may be caused by annotation file).
#' @param ignore_chrs Chromsomes to ignore from analysis. e.g. chrX and chrY.
#' @param use_syn Logical. If `TRUE`, include synonymous variants in analysis.
#' @references Mayakonda, Anand, et al. "Maftools: efficient and comprehensive analysis of somatic variants in cancer." Genome research 28.11 (2018): 1747-1756.
#' @references Roberts SA, Lawrence MS, Klimczak LJ, et al. An APOBEC Cytidine Deaminase Mutagenesis Pattern is Widespread in Human Cancers. Nature genetics. 2013;45(9):970-976. doi:10.1038/ng.2702.
#' @references Bergstrom EN, Huang MN, Mahto U, Barnes M, Stratton MR, Rozen SG, Alexandrov LB: SigProfilerMatrixGenerator: a tool for visualizing and exploring patterns of small mutational events. BMC Genomics 2019, 20:685 https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-6041-2
#' @export
sig_tally.MAF <- function(object, mode = c("SBS", "DBS", "ID", "ALL"),
                          ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
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

  hsgs.installed <- BSgenome::installed.genomes(splitNameParts = TRUE)
  data.table::setDT(x = hsgs.installed)

  if (nrow(hsgs.installed) == 0) {
    send_stop("Could not find any installed BSgenomes. Use {.code BSgenome::available.genomes()} for options.")
  }

  send_info("We would assume you marked all variants' position in + strand.")
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

  if (is.null(genome_build)) {
    if (grepl("hg19", ref_genome)) {
      genome_build <- "hg19"
    } else if (grepl("hg38", ref_genome)) {
      genome_build <- "hg38"
    } else if (grepl("mm10$", ref_genome)) {
      genome_build <- "mm10"
    } else if (grepl("mm9$", ref_genome)) {
      genome_build <- "mm9"
    } else {
      send_stop("Cannot guess the genome build, please set it by hand!")
    }
  }

  ref_genome <- BSgenome::getBSgenome(genome = ref_genome)
  send_success("Reference genome loaded.")

  query <- maftools::subsetMaf(
    maf = object,
    query = "Variant_Type %in% c('SNP', 'DNP', 'INS', 'DEL')", fields = "Chromosome",
    includeSyn = use_syn, mafObj = FALSE
  )
  # Check NA in Reference_Allele Tumor_Seq_Allele2
  query <- query[!is.na(query$Reference_Allele) & !is.na(query$Tumor_Seq_Allele2)]
  if (identical(query$Reference_Allele, query$Tumor_Seq_Allele2)) {
    send_stop("Tumor_Seq_Allele2 (mutated allele) should not equal to Reference_Allele!")
  }
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
  query$Chromosome <- ifelse(startsWith(query$Chromosome, "chr"),
    query$Chromosome,
    paste0("chr", query$Chromosome)
  )

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

  # detect and transform chromosome 23 to "X"
  query$Chromosome <- sub("23", "X", query$Chromosome)
  # detect and transform chromosome 24 to "Y"
  query$Chromosome <- sub("24", "Y", query$Chromosome)
  send_success("Sex chromosomes properly handled.")

  # only keep standard chromosomes
  query <- query[query$Chromosome %in% paste0("chr", c(1:22, "X", "Y", "M", "MT"))]
  send_success("Only variants located in standard chromosomes (1:22, X, Y, M/MT) are kept.")

  query$Start_Position <- as.numeric(as.character(query$Start_Position))
  query$End_Position <- as.numeric(as.character(query$End_Position))

  send_success("Variant start and end position checked.")

  query_seq_lvls <- query[, .N, Chromosome]
  ref_seqs_lvls <- BSgenome::seqnames(x = ref_genome)
  query_seq_lvls_missing <- query_seq_lvls[!Chromosome %in% ref_seqs_lvls]

  if (nrow(query_seq_lvls_missing) > 3) {
    ## Some reference genome builds have no 'chr' prefix
    send_warning("Too many chromosome names cannot match reference genome. Try dropping 'chr' prefix to fix it...")
    query$Chromosome <- sub(
      pattern = "chr",
      replacement = "",
      x = as.character(query$Chromosome),
      ignore.case = TRUE
    )
    query_seq_lvls <- query[, .N, Chromosome]
    query_seq_lvls_missing <- query_seq_lvls[!Chromosome %in% ref_seqs_lvls]
    send_info("Dropped.")
  }

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
        if (e$message == "") {
          NULL
        } else {
          send_error("Unexpected error occurred:")
          send_stop(e$message)
        }
      }
    )
    res_DBS <- tryCatch(
      generate_matrix_DBS(query, ref_genome, genome_build = genome_build, add_trans_bias = add_trans_bias),
      error = function(e) {
        if (e$message == "") {
          NULL
        } else {
          send_error("Unexpected error occurred:")
          send_stop(e$message)
        }
      }
    )
    res_ID <- tryCatch(
      generate_matrix_INDEL(query, ref_genome, genome_build = genome_build, add_trans_bias = add_trans_bias),
      error = function(e) {
        if (e$message == "") {
          NULL
        } else {
          send_error("Unexpected error occurred:")
          send_stop(e$message)
        }
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
