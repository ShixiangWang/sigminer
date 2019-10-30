#' Get Signature Exposure from 'Signature' Object
#'
#' @param Signature a `Signature` object obtained either from [sig_extract] or [sig_auto_extract],
#' or just a raw exposure matrix with column representing samples (patients) and row
#' representing signatures.
#' @param type 'absolute' for signature exposure and 'relative' for signature relative exposure.
#' @return a `tibble`
#' @export
#'
#' @examples
#' # Load mutational signature
#' load(system.file("extdata", "toy_mutational_signature.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Get signature exposure
#' get_sig_exposure(sig2)
#' get_sig_exposure(sig2, type = "relative")
get_sig_exposure <- function(Signature, type = c("absolute", "relative")) {
  if (class(Signature) == "Signature") {
    h <- Signature$Exposure
  } else if (is.matrix(Signature)) {
    if (!all(startsWith(rownames(Signature), "Sig"))) {
      stop("If Signature is a matrix, column names must start with 'Sig'!", call. = FALSE)
    }
    h <- Signature
  } else {
    stop("Invalid input for 'Signature'", call. = FALSE)
  }

  if (is.null(rownames(h)) | is.null(colnames(h))) {
    stop("Rownames or Colnames cannot be NULL!")
  }

  type <- match.arg(type)
  if (type == "absolute") {
    h <- t(h) %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "sample") %>%
      dplyr::as_tibble()
    return(h)
  } else {
    h.norm <- apply(h, 2, function(x) x / sum(x))
    h.norm <- t(h.norm) %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "sample") %>%
      dplyr::as_tibble()
    return(h.norm)
  }
}
