## Used in sig_extract and sig_auto_extract
##
## By scaling the signature loading matrix has all mutation burdens
## - each signature (column in W) now represents a number of mutations
## assigned to each signature.
##

helper_scale_nmf_matrix <- function(Signature, Exposure, K, handle_cn = FALSE) {
  W <- Signature
  H <- Exposure

  if (handle_cn) {
    ## If copy number exposure is not handled by caller
    ## then handle it here
    has_cn <- grepl("^CN[^C]", rownames(Signature)) | startsWith(rownames(Signature), "copynumber")
    if (any(has_cn)) {
      # Handle copy number exposure
      W <- W[has_cn, ]
    }
  }

  for (j in seq_len(K)) {
    Signature[, j] <- Signature[, j] * rowSums(H)[j]
    Exposure[j, ] <- Exposure[j, ] * colSums(W)[j]
  }

  return(list(Signature = Signature, Exposure = Exposure))
}
