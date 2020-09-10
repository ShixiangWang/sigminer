#' Calculate Similarity between Identified Signatures and Reference Signatures
#'
#' The reference signatures can be either a `Signature` object specified by `Ref` argument
#' or known COSMIC signatures specified by `sig_db` argument.
#' Two COSMIC databases are used for comparisons - "legacy" which includes 30 signaures,
#' and "SBS" - which includes updated/refined 65 signatures. This function is modified
#' from `compareSignatures()` in **maftools** package.
#'
#' @param Signature a `Signature` object or a component-by-signature matrix (sum of each column is 1)
#' or a normalized component-by-sample matirx (sum of each column is 1).
#' More please see examples.
#' @param Ref default is `NULL`, can be a same object as `Signature`.
#' @param sig_db can be 'legacy' (for COSMIC v2 'SBS'),
#' 'SBS', 'DBS', 'ID' and 'TSB' (for SBS transcriptional strand bias signatures).
#' Default 'legacy'.
#' @param db_type only used when `sig_db` is enabled.
#' "" for keeping default, "human-exome" for transforming to exome frequency of component,
#' and "human-genome" for transforming to whole genome frequency of component.
#' Currently only works for 'SBS'.
#' @param method default is 'cosine' for cosine similarity.
#' @param normalize one of "row" and "feature". "row" is typically used
#' for common mutational signatures. "feature" is designed by me to use when input
#' are copy number signatures.
#' @param pattern_to_rm patterns for removing some features/components in similarity
#' calculation. A vector of component name is also accepted.
#' The remove operation will be done after normalization. Default is `NULL`.
#' @param verbose if `TRUE`, print extra info.
#' @inheritParams sig_tally
#' @author Shixiang Wang <w_shixiang@163.com>
#'
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
#' ## Remove some components
#' ## in similarity calculation
#' s4 <- get_sig_similarity(sig2,
#'   Ref = sig2,
#'   pattern_to_rm = c("T[T>G]C", "T[T>G]G", "T[T>G]T")
#' )
#' s4
#'
#' ## Same to DBS and ID signatures
#' @testexamples
#' expect_equal(length(s1), 4L)
#' expect_equal(length(s2), 4L)
#' expect_equal(length(s3), 4L)
#' expect_equal(length(s4), 4L)
get_sig_similarity <- function(Signature, Ref = NULL,
                               sig_db = "legacy",
                               db_type = c("", "human-exome", "human-genome"),
                               method = "cosine",
                               normalize = c("row", "feature"),
                               feature_setting = sigminer::CN.features,
                               pattern_to_rm = NULL,
                               verbose = TRUE) {
  if (inherits(Signature, "Signature")) {
    w <- Signature$Signature.norm
  } else if (is.matrix(Signature)) {
    w <- Signature
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

  sig_db <- match.arg(arg = sig_db, choices = c("legacy", "SBS", "DBS", "ID", "TSB"))
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
    } else if (is.matrix(Ref)) {
      sigs <- Ref
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
    } else {
      message("-Comparing against Custom signatures")
    }
    message("------------------------------------")
    for (i in 1:nrow(corMat)) {
      if (is.null(Ref)) {
        ae <- aetiology[names(which(corMat[i, ] == max(corMat[i, ]))), ]
        ae <- paste0(
          "Aetiology: ", ae, " [similarity: ",
          max(corMat[i, ]), "]"
        )
        message("--Found ", rownames(corMat)[i], " most similar to ",
          names(which(corMat[i, ] == max(corMat[i, ]))),
          sep = ""
        )
        message(paste0("   ", ae))
      } else {
        message("--Found ", rownames(corMat)[i], " most similar to ",
          names(which(corMat[i, ] == max(corMat[i, ]))),
          paste0(" [similarity: ", max(corMat[i, ]), "]"),
          sep = ""
        )
      }
    }
    message("------------------------------------")
    message("Return result invisiblely.")
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


# Get Reference Signature Database ----------------------------------------

#' Obtain Reference Signatures
#'
#' The signatures and their aetiologies mainly obtained from COSMIC database and cleaned before saving into
#' **sigminer** package.
#'
#' @inheritParams get_sig_similarity
#'
#' @return a `list`.
#' @export
#' @seealso [get_sig_similarity], [sig_fit] and [show_cosmic_sig_profile].
#'
#' @examples
#' s1 <- get_sig_db()
#' s2 <- get_sig_db("DBS")
#' s1
#' s2
#' @testexamples
#' expect_is(s1, "list")
#' expect_is(s2, "list")
get_sig_db <- function(sig_db = "legacy") {
  db_file <- switch(
    sig_db,
    legacy = system.file("extdata", "legacy_signatures.RDs",
      package = "maftools", mustWork = TRUE
    ),
    SBS = system.file("extdata", "SBS_signatures.RDs",
      package = "maftools", mustWork = TRUE
    ),
    DBS = system.file("extdata", "DBS_signatures.rds",
      package = "sigminer", mustWork = TRUE
    ),
    ID = system.file("extdata", "ID_signatures.rds",
      package = "sigminer", mustWork = TRUE
    ),
    TSB = system.file("extdata", "TSB_signatures.rds",
      package = "sigminer", mustWork = TRUE
    )
  )
  sigs_db <- readRDS(file = db_file)
  ## Make sure column-sum is 1, i.e. normalized
  sigs_db$db <- apply(sigs_db$db, 2, function(x) x / sum(x))
  sigs_db
}
