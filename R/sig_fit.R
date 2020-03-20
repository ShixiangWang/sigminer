#' Fit Signature Exposure with Linear Combination Decomposition
#'
#' The function performs a signatures decomposition of a given mutational
#' catalogue `V` with known signatures `W` by solving the minimization problem
#' `min(||W*H - V||)` with additional constraints of non-negativity on H where W and V are known.
#'
#' This is a modification based on `LCD` function from YAPSA pakcage.
#'
#' @param catalogue_matrix a numeric matrix `V` with row representing components and
#' columns representing samples, typically you can get `nmf_matrix` from `sig_tally()` and
#' transpose it by `t()`.
#' @inheritParams show_sig_profile
#' @inheritParams show_cosmic_sig_profile
#' @inheritParams get_sig_exposure
#' @param return_class string, 'matrix' or 'data.table'.
#' @param rel_threshold numeric vector, a relative exposure lower than this value will be set to 0.
#'
#' @return The exposure result either in `matrix` or `data.table` format.
#' @export
#'
#' @examples
#' W <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
#' colnames(W) <- c("sig1", "sig2")
#' W <- apply(W, 2, function(x) x / sum(x))
#'
#' H <- matrix(c(2, 5, 3, 6, 1, 9, 1, 2), ncol = 4)
#' colnames(H) <- paste0("samp", 1:4)
#'
#' V <- W %*% H
#' V
#'
#' H_infer <- sig_fit(V, W)
#' H_infer
#' H
#'
#' H_dt <- sig_fit(V, W, return_class = "data.table")
#' H_dt
#' @testexamples
#' expect_is(H_infer, "matrix")
#' expect_is(H_dt, "data.table")
sig_fit <- function(catalogue_matrix,
                    sig,
                    sig_index = NULL,
                    sig_db = "legacy",
                    show_index = TRUE,
                    type = c("absolute", "relative"),
                    return_class = c("matrix", "data.table"),
                    rel_threshold = 0,
                    mode = c("SBS", "copynumber")) {
  stopifnot(is.matrix(catalogue_matrix))

  if (is.null(sig_index)) {
    if (class(sig) == "Signature") {
      sig_matrix <- sig$Signature
    } else if (is.matrix(sig)) {
      sig_matrix <- sig
    } else {
      stop("Invalid input for 'sig'", call. = FALSE)
    }
  } else {
    if (sig_db == "legacy") {
      sigs_db <- readRDS(file = system.file("extdata", "legacy_signatures.RDs",
        package = "maftools", mustWork = TRUE
      ))
      sigs <- sigs_db$db
      avail_index <- substring(colnames(sigs), 8)
    } else {
      sigs_db <- readRDS(file = system.file("extdata", "SBS_signatures.RDs",
        package = "maftools", mustWork = TRUE
      ))
      sigs <- sigs_db$db
      avail_index <- substring(colnames(sigs), 4)
    }

    if (show_index) {
      msg <- paste(
        paste0("\nValid index for db '", sig_db, "':"),
        paste0(avail_index, collapse = " "),
        sep = "\n"
      )
      message(msg)
    }

    if (!is.character(sig_index)) {
      sig_index <- as.character(sig_index)
    }

    if (!all(sig_index %in% avail_index)) {
      stop(msg)
    }

    index <- c()
    for (i in sig_index) {
      index <- c(index, colnames(sigs)[avail_index == i])
    }

    sig_matrix <- as.matrix(sigs[, index, drop = FALSE])
  }

  mode <- match.arg(mode)
  type <- match.arg(type)
  return_class <- match.arg(return_class)

  if (mode == "copynumber") {
    ## For copy number signature
    ## Only feature of copy number used for exposure quantification
    ## So the result exposure is estimated copy number segments
    has_cn_w <- grepl("^CN[^C]", rownames(sig_matrix)) | startsWith(rownames(sig_matrix), "copynumber")
    sig_matrix <- sig_matrix[has_cn_w, , drop = FALSE]
    catalogue_matrix <- catalogue_matrix[rownames(sig_matrix), , drop = FALSE]
  }

  ## Keep sum of signature to 1
  sig_matrix <- apply(sig_matrix, 2, function(x) x / sum(x))


  ## Check if V and W have same rows
  if (nrow(catalogue_matrix) != nrow(sig_matrix)) {
    stop("Catalogue matrix and Signature matrix should have same rows, please check!")
  }

  ## If V and W have row names, check the order
  cat_rowname <- rownames(catalogue_matrix)
  sig_rowname <- rownames(sig_matrix)
  if (!is.null(cat_rowname) & !is.null(sig_rowname)) {
    if (!all(sig_rowname == cat_rowname)) {
      message("Matrix V and W don't have same orders. Try reordering...")
      if (all(sort(cat_rowname) == sort(sig_rowname))) {
        ## Set catalogue matrix as signature matrix
        catalogue_matrix <- catalogue_matrix[sig_rowname, , drop = FALSE]
      } else {
        stop("The rownames of matrix are not identical, please check your input!")
      }
    }
  }

  # Set constraints x >= 0
  G <- diag(dim(sig_matrix)[2])
  H <- rep(0, dim(sig_matrix)[2])

  expo <- purrr::map2(as.data.frame(catalogue_matrix), rel_threshold, function(x, y, type = "absolute") {
    expo <- lsei::lsei(
      a = sig_matrix,
      b = x,
      e = G, f = H
    )
    rel_expo <- expo / sum(expo)
    expo[rel_expo < y] <- 0
    if (type == "relative") {
      expo <- expo / sum(expo)
    }
    expo <- round(expo, digits = 3)
    expo
  }, type = type)

  expo <- dplyr::bind_rows(expo) %>%
    as.matrix()
  rownames(expo) <- colnames(sig_matrix)
  colnames(expo) <- colnames(catalogue_matrix)

  if (return_class != "matrix") {
    expo <- expo %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Sig") %>%
      tidyr::pivot_longer(cols = -"Sig", names_to = "sample", values_to = "expo") %>%
      tidyr::pivot_wider(id_cols = "sample", names_from = "Sig", values_from = "expo") %>%
      data.table::as.data.table()
  }
  expo
}
