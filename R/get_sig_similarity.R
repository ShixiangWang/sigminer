#' Calculate Similarity between Identified Signatures and Reference Signatures
#'
#' The reference signatures can be either a `Signature` object specified by `Ref` argument
#' or known COSMIC signatures specified by `sig_db` argument.
#' Two COSMIC databases are used for comparisons - "legacy" which includes 30 signaures,
#' and "SBS" - which includes updated/refined 65 signatures. This function is modified
#' from `compareSignatures()` in **maftools** package.
#' **NOTE**: all reference signatures are generated from gold standard tool:
#' SigProfiler.
#'
#' @param Signature a `Signature` object or a component-by-signature matrix/`data.frame`
#' (sum of each column is 1) or a normalized component-by-sample matrix/`data.frame`
#' (sum of each column is 1).
#' More please see examples.
#' @param Ref default is `NULL`, can be a same object as `Signature`.
#' @param sig_db default 'legacy', it can be 'legacy' (for [COSMIC v2 'SBS'](https://cancer.sanger.ac.uk/cosmic/signatures_v2.tt/)),
#' 'SBS', 'DBS', 'ID' and 'TSB' (for [COSMIV v3.1 signatures](https://cancer.sanger.ac.uk/cosmic/signatures/))
#' for small scale mutations.
#' For more specific details, it can also be 'SBS_hg19', 'SBS_hg38',
#' 'SBS_mm9', 'SBS_mm10', 'DBS_hg19', 'DBS_hg38', 'DBS_mm9', 'DBS_mm10' to use
#' COSMIC v3 reference signatures from Alexandrov, Ludmil B., et al. (2020) (reference #1).
#' In addition, it can be one of "SBS_Nik_lab_Organ", "RS_Nik_lab_Organ",
#' "SBS_Nik_lab", "RS_Nik_lab" to refer reference signatures from
#' Degasperi, Andrea, et al. (2020) (reference #2);
#' "RS_BRCA560", "RS_USARC" to reference signatures from BRCA560 and USARC cohorts;
#' "CNS_USARC" (40 categories), "CNS_TCGA" (48 categories) to reference copy number signatures from USARC cohort and TCGA;
#' "CNS_TCGA176" (176 categories) and "CNS_PCAWG176" (176 categories) to reference copy number signatures from PCAWG and TCGA separately.
#' **UPDATE**, the latest version of reference version can be automatically
#' downloaded and loaded from <https://cancer.sanger.ac.uk/signatures/downloads/>
#' when a option with `latest_` prefix is specified (e.g. "latest_SBS_GRCh37").
#' **Note**: the signature profile for different genome builds are basically same.
#' And specific database (e.g. 'SBS_mm10') contains less signatures than all COSMIC
#' signatures (because some signatures are not detected from Alexandrov, Ludmil B., et al. (2020)).
#' For all available options, check the parameter setting.
#'
#' @param db_type only used when `sig_db` is enabled.
#' "" for keeping default, "human-exome" for transforming to exome frequency of component,
#' and "human-genome" for transforming to whole genome frequency of component.
#' Currently only works for 'SBS'.
#' @param method default is 'cosine' for cosine similarity.
#' @param normalize one of "row" and "feature". "row" is typically used
#' for common mutational signatures. "feature" is designed by me to use when input
#' are copy number signatures.
#' @param set_order if `TRUE`, order the return similarity matrix.
#' @param pattern_to_rm patterns for removing some features/components in similarity
#' calculation. A vector of component name is also accepted.
#' The remove operation will be done after normalization. Default is `NULL`.
#' @param verbose if `TRUE`, print extra info.
#' @inheritParams sig_tally
#' @author Shixiang Wang <w_shixiang@163.com>
#' @references
#' Alexandrov, Ludmil B., et al. "The repertoire of mutational signatures in human cancer." Nature 578.7793 (2020): 94-101.
#'
#' Degasperi, Andrea, et al. "A practical framework and online tool for mutational signature analyses show intertissue variation and driver dependencies." Nature cancer 1.2 (2020): 249-263.
#'
#' Steele, Christopher D., et al. "Undifferentiated sarcomas develop through distinct evolutionary pathways." Cancer Cell 35.3 (2019): 441-456.
#'
#' Nik-Zainal, Serena, et al. "Landscape of somatic mutations in 560 breast cancer whole-genome sequences." Nature 534.7605 (2016): 47-54.
#'
#' Steele, Christopher D., et al. "Signatures of copy number alterations in human cancer." Nature 606.7916 (2022): 984-991.
#' @return a `list` containing smilarities, aetiologies if available, best match and RSS.
#' @export
#'
#' @examples
#' # Load mutational signature
#' load(system.file("extdata", "toy_mutational_signature.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#'
#' s1 <- get_sig_similarity(sig2, Ref = sig2)
#' s1
#'
#' s2 <- get_sig_similarity(sig2)
#' s2
#' s3 <- get_sig_similarity(sig2, sig_db = "SBS")
#' s3
#'
#' # Set order for result similarity matrix
#' s4 <- get_sig_similarity(sig2, sig_db = "SBS", set_order = TRUE)
#' s4
#'
#' ## Remove some components
#' ## in similarity calculation
#' s5 <- get_sig_similarity(sig2,
#'   Ref = sig2,
#'   pattern_to_rm = c("T[T>G]C", "T[T>G]G", "T[T>G]T")
#' )
#' s5
#'
#' ## Same to DBS and ID signatures
#' x1 <- get_sig_db("DBS_hg19")
#' x2 <- get_sig_db("DBS_hg38")
#' s6 <- get_sig_similarity(x1$db, x2$db)
#' s6
#' @testexamples
#' expect_equal(length(s1), 4L)
#' expect_equal(length(s2), 4L)
#' expect_equal(length(s3), 4L)
#' expect_equal(length(s4), 4L)
#' expect_equal(length(s5), 4L)
#' expect_is(s6, "list")
get_sig_similarity <- function(Signature, Ref = NULL,
                               sig_db = c(
                                 "SBS", "legacy", "DBS", "ID", "TSB",
                                 "SBS_Nik_lab", "RS_Nik_lab",
                                 "RS_BRCA560", "RS_USARC",
                                 "CNS_USARC", "CNS_TCGA",
                                 "CNS_TCGA176","CNS_PCAWG176",
                                 "SBS_hg19", "SBS_hg38", "SBS_mm9", "SBS_mm10",
                                 "DBS_hg19", "DBS_hg38", "DBS_mm9", "DBS_mm10",
                                 "SBS_Nik_lab_Organ", "RS_Nik_lab_Organ",
                                 "latest_SBS_GRCh37", "latest_DBS_GRCh37", "latest_ID_GRCh37",
                                 "latest_SBS_GRCh38", "latest_DBS_GRCh38",
                                 "latest_SBS_mm9", "latest_DBS_mm9",
                                 "latest_SBS_mm10", "latest_DBS_mm10",
                                 "latest_SBS_rn6", "latest_DBS_rn6",
                                 "latest_CN_GRCh37", "latest_RNA-SBS_GRCh37", "latest_SV_GRCh38"
                               ),
                               db_type = c("", "human-exome", "human-genome"),
                               method = "cosine",
                               normalize = c("row", "feature"),
                               feature_setting = sigminer::CN.features,
                               set_order = TRUE,
                               pattern_to_rm = NULL,
                               verbose = TRUE) {
  if (inherits(Signature, "Signature")) {
    w <- Signature$Signature.norm
  } else if (is.matrix(Signature) | is.data.frame(Signature)) {
    w <- as.matrix(Signature)
    if (!all(round(colSums(w)) == 1)) {
      stop("If Signature is a matrix, sum of each column is 1!", call. = FALSE)
    }
  } else {
    stop("Invalid input for 'Signature'", call. = FALSE)
  }

  normalize <- match.arg(normalize)
  db_type <- match.arg(db_type)

  if (normalize == "feature") {
    use_W <- any(grepl("\\[.*\\]$", rownames(w)))
    use_M <- any(grepl("\\d+$", rownames(w)))

    if (use_W | use_M) {
      w <- helper_normalize_by_feature(w, feature_setting)
    } else {
      message("normalize='feature' only supports copy number signatures. Please use other value!")
      return(invisible())
    }
  }

  sig_db <- match.arg(arg = sig_db)
  method <- match.arg(arg = method, choices = c("cosine"))

  if (is.null(Ref)) {
    sigs_db <- get_sig_db(sig_db)
    sigs <- sigs_db$db
    aetiology <- sigs_db$aetiology

    ## Some extra processing
    if (sig_db == "legacy" & db_type == "human-genome") {
      ## v2 comes from Exome
      sigs <- sig_convert(sig = sigs, from = "human-exome", to = "human-genome")
    } else if (sig_db == "SBS" & db_type == "human-exome") {
      ## v3 comes from WGS (PCAWG)
      ## Should DBS and ID also handle such cases?
      sigs <- sig_convert(sig = sigs, from = "human-genome", to = "human-exome")
    }
  } else {
    if (inherits(Ref, "Signature")) {
      sigs <- Ref$Signature.norm
    } else if (is.matrix(Ref) | is.data.frame(Ref)) {
      sigs <- as.matrix(Ref)
      if (!all(round(colSums(sigs)) == 1)) {
        stop("If Ref is a matrix, sum of each column is 1!", call. = FALSE)
      }
    } else {
      stop("Invalid input for 'Ref'", call. = FALSE)
    }

    if (normalize == "feature") {
      use_W <- any(grepl("\\[.*\\]$", rownames(sigs)))
      use_M <- any(grepl("\\d+$", rownames(sigs)))

      if (use_W | use_M) {
        sigs <- helper_normalize_by_feature(sigs, feature_setting)
      } else {
        message("normalize='feature' only supports copy number signatures. Please use other value!")
        return(invisible())
      }
    }
  }

  if (!all(rownames(w) %in% rownames(sigs))) {
    stop(paste0(
      "The following components cannot be found in reference!\n",
      paste(rownames(sigs)[!rownames(w) %in% rownames(sigs)], collapse = " ")
    ))
  }

  # Match the components
  sigs <- sigs[rownames(w), ]

  # Remove if pattern is set
  if (!is.null(pattern_to_rm)) {
    index2rm <- sapply(pattern_to_rm, grepl, x = rownames(w), fixed = TRUE)

    if (is.matrix(index2rm)) {
      index2rm <- apply(index2rm, 1, any)
    }

    if (!any(index2rm)) {
      warning("Pattern to remove set, but no component found to remove.")
    } else {
      sigs <- sigs[!index2rm, ]
      w <- w[!index2rm, ]
    }
  }

  corMat <- round(cosineMatrix(w, sigs), digits = 3)
  rownames(corMat) <- colnames(w)
  colnames(corMat) <- colnames(sigs)

  RssMat <- c()
  for (i in 1:ncol(w)) {
    sig <- w[, i]
    RssMat <- rbind(RssMat, apply(sigs, 2, function(x) {
      round(sum((sig - x)^2), digits = 6)
    }))
  }
  rownames(RssMat) <- colnames(w)
  colnames(RssMat) <- colnames(sigs)

  if (!exists("aetiology")) {
    aetiology <- NULL
  }

  best_matches <- lapply(1:nrow(corMat), function(i) {
    if (!is.null(aetiology)) {
      ae <- aetiology[names(which(corMat[i, ] == max(corMat[i, ]))), ]
    } else {
      ae <- NA
    }
    max_cor <- names(which(corMat[i, ] == max(corMat[i, ])))
    list(aetiology = ae, best_match = paste0(
      "Best match: ",
      max_cor,
      " [similarity: ",
      max(corMat[i, ]), "]"
    ))
  })
  names(best_matches) <- rownames(corMat)
  if (verbose) {
    if (is.null(Ref)) {
      message("-Comparing against COSMIC signatures")
      message("------------------------------------")
    } else {
      message("-Comparing against Custom signatures/profiles")
      message("---------------------------------------------")
    }

    all_matches <- c()
    for (i in 1:nrow(corMat)) {
      .to_match <- rownames(corMat)[i]
      .be_match <- names(which(corMat[i, ] == max(corMat[i, ])))
      all_matches <- c(all_matches, .be_match)

      if (is.null(Ref)) {
        ae <- aetiology[.be_match, ]
        ae <- paste0(
          "Aetiology: ", ae, " [similarity: ",
          max(corMat[i, ]), "]"
        )
        message("--Found ", .to_match,
          " most similar to ", .be_match,
          sep = ""
        )
        message(paste0("   ", ae))
      } else {
        message("--Found ", .to_match,
          " most similar to ", .be_match,
          paste0(" [similarity: ", max(corMat[i, ]), "]"),
          sep = ""
        )
      }
    }
    message("------------------------------------")
    message("Return result invisiblely.")
  }

  if (set_order) {
    all_matches <- unique(all_matches)
    corMat <- corMat[, c(all_matches, setdiff(colnames(corMat), all_matches)), drop = FALSE]
  }

  res <- list(
    similarity = corMat,
    aetiology_db = ifelse(!is.null(aetiology), aetiology, NA),
    best_match = best_matches,
    rss = RssMat
  )
  class(res) <- c("similarity", class(res))

  invisible(res)
}
