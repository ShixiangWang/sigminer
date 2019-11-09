#' Calculate Similarity between Identified Signatures and Reference Signatures
#'
#' The reference signatures can be either a `Signature` object specified by `Ref` argument
#' or known COSMIC signatures specified by `sig_db` argument.
#' Two COSMIC databases are used for comparisons - "legacy" which includes 30 signaures,
#' and "SBS" - which includes updated/refined 65 signatures. This function is modified
#' from `compareSignatures()` in **maftools** package.
#'
#' @param Signature a `Signature` object or a component-by-signature matrix (sum of each column is 1).
#' More please see examples.
#' @param Ref default is `NULL`, can be a same object as `Sig`.
#' @param sig_db can be 'legacy' or 'SBS'. Default 'legacy'.
#' @param method default is 'cosine' for cosine similarity.
#' @param verbose if `TRUE`, print extra info.
#' @author Shixiang Wang <w_shixiang@163.com>
#'
#' @return a `list` containing smilarities, aetiologies if available, and best match.
#' @export
#'
#' @examples
#' # Load mutational signature
#' load(system.file("extdata", "toy_mutational_signature.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#'
#' get_sig_similarity(sig2, Ref = sig2)
#' \donttest{
#' get_sig_similarity(sig2)
#' get_sig_similarity(sig2, sig_db = "SBS")
#' }
get_sig_similarity <- function(Signature, Ref = NULL, sig_db = "legacy", method = "cosine", verbose = TRUE) {
  if (class(Signature) == "Signature") {
    w <- Signature$Signature.norm
  } else if (is.matrix(Signature)) {
    w <- Signature
    if (!all(colSums(w) == 1)) {
      stop("If Signature is a matrix, sum of each column is 1!", call. = FALSE)
    }
  } else {
    stop("Invalid input for 'Signature'", call. = FALSE)
  }

  sig_db <- match.arg(arg = sig_db, choices = c("legacy", "SBS"))
  method <- match.arg(arg = method, choices = c("cosine"))

  if (is.null(Ref)) {
    if (sig_db == "legacy") {
      sigs_db <- readRDS(file = system.file("extdata", "legacy_signatures.RDs",
        package = "maftools", mustWork = TRUE
      ))
      sigs <- sigs_db$db
      aetiology <- sigs_db$aetiology
    }
    else {
      sigs_db <- readRDS(file = system.file("extdata", "SBS_signatures.RDs",
        package = "maftools", mustWork = TRUE
      ))
      sigs <- sigs_db$db
      aetiology <- sigs_db$aetiology
    }
  } else {
    if (class(Ref) == "Signature") {
      sigs <- Ref$Signature.norm
    } else if (is.matrix(Ref)) {
      sigs <- Ref
      if (!all(colSums(sigs) == 1)) {
        stop("If Ref is a matrix, sum of each column is 1!", call. = FALSE)
      }
    } else {
      stop("Invalid input for 'Ref'", call. = FALSE)
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

  if (method == "cosine") {
    corMat <- c()
    for (i in 1:ncol(w)) {
      sig <- w[, i]
      corMat <- rbind(corMat, apply(sigs, 2, function(x) {
        round(crossprod(sig, x) / sqrt(crossprod(x) * crossprod(sig)),
          digits = 3
        )
      }))
    }
    rownames(corMat) <- colnames(w)
  } else {
    # Other methods
  }

  best_matches <- lapply(1:nrow(corMat), function(i) {
    if (exists("aetiology", envir = parent.frame())) {
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

  invisible(list(
    similarity = corMat,
    aetiology_db = ifelse(exists("aetiology"), aetiology, NA),
    best_match = best_matches
  ))
}
