#' Obtain or Modify Signature Information
#'
#' @param sig a `Signature` object obtained either from [sig_extract] or [sig_auto_extract].
#' @name sig_operation
#' @return a `Signature` object or data.
#' @export
#'
#' @examples
#' ## Operate signature names
#' load(system.file("extdata", "toy_mutational_signature.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' sig_names(sig2)
#' cc <- sig_modify_names(sig2, new_names = c("Sig2", "Sig1", "Sig3"))
#' sig_names(cc)
#'
#' # The older names are stored in tags.
#' print(attr(cc, "tag"))
#' @testexamples
#' expect_is(cc, "Signature")
sig_names <- function(sig) {
  stopifnot(inherits(sig, "Signature"))
  colnames(sig$Signature)
}

#' @param new_names new signature names.
#' @rdname sig_operation
#' @export
sig_modify_names <- function(sig, new_names) {
  ns <- sig_names(sig)
  if (length(new_names) != length(ns)) {
    stop("The new names should have same length as old names!")
  }

  colnames(sig$Signature) <- colnames(sig$Signature.norm) <- colnames(sig$Raw$W) <- rownames(sig$Exposure) <- rownames(sig$Exposure.norm) <- rownames(sig$Raw$H) <- new_names

  attr(sig, "tag") <- paste0("Older names:", paste(ns, collapse = ","), ";")
  sig
}

#' @rdname sig_operation
#' @export
#' @examples
#' ## Get signature number
#' sig_number(sig2)
sig_number <- function(sig) {
  stopifnot(inherits(sig, "Signature"))
  sig$K
}

#' @rdname sig_operation
#' @export
#' @examples
#' ## Get signature attributes
#' sig_number(sig2)
sig_attrs <- function(sig) {
  stopifnot(inherits(sig, "Signature"))
  attributes(sig)[!names(attributes(sig)) %in% c("names", "class")]
}

#' @param normalize one of 'row', 'column', 'raw' and "feature", for row normalization (signature),
#' column normalization (component), raw data, row normalization by feature, respectively.
#' @rdname sig_operation
#' @export
#' @examples
#' ## Get signature matrix
#' z <- sig_signature(sig2)
#' z <- sig_signature(sig2, normalize = "raw")
sig_signature <- function(sig, normalize = c("row", "column", "raw", "feature")) {
  stopifnot(inherits(sig, "Signature"))

  normalize <- match.arg(normalize)
  sig <- sig$Signature

  if (normalize == "row") {
    sig <- apply(sig, 2, function(x) x / sum(x))
  } else if (normalize == "column") {
    sig <- t(apply(sig, 1, function(x) x / sum(x)))
  } else if (normalize == "feature") {
    has_cn <- grepl("^CN[^C]", rownames(sig)) | startsWith(rownames(sig), "copynumber")

    if (!any(has_cn)) {
      stop("normalize method 'feature' is only suitable for copy number signature!")
    }

    mat <- as.data.frame(sig)
    mat$context <- rownames(mat)

    if (any(grepl("^CN[^C]", rownames(sig)))) {
      mat$base <- sub("\\[.*\\]$", "", mat$context)
    } else {
      mat$base <- sub("\\d+$", "", mat$context)
    }

    mat <- tidyr::gather(mat, class, signature, -c("context", "base"))

    mat <- mat %>%
      dplyr::group_by(.data$base, .data$class) %>%
      tidyr::nest() %>%
      dplyr::mutate(
        context = purrr::map(.data$data, ~ .$context),
        signature = purrr::map(.data$data, ~ .$signature / sum(.$signature))
      ) %>%
      dplyr::select(-"data") %>%
      tidyr::unnest(cols = c("context", "signature")) %>%
      dplyr::ungroup()

    mat <- tidyr::pivot_wider(mat, names_from = "class", values_from = "signature") %>%
      dplyr::select(-"base") %>%
      tibble::column_to_rownames("context") %>%
      as.matrix()

    # Keep row order same
    sig <- mat[rownames(sig), ]
  }

  return(sig)
}

#' @param type one of 'absolute' and 'relative'.
#' @rdname sig_operation
#' @export
#' @examples
#' ## Get exposure matrix
#' ## Of note, this is different from get_sig_exposure()
#' ## it returns a matrix instead of data table.
#' z <- sig_exposure(sig2) # it is same as sig$Exposure
#' z <- sig_exposure(sig2, type = "relative") # it is same as sig2$Exposure.norm
sig_exposure <- function(sig, type = c("absolute", "relative")) {
  stopifnot(inherits(sig, "Signature"))
  type <- match.arg(type)
  if (type == "absolute") {
    return(sig$Exposure)
  } else {
    return(sig$Exposure.norm)
  }
}
