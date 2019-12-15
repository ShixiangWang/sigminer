helper_sort_signature = function(sig) {
  stopifnot(is.matrix(sig), !is.null(rownames(sig)))
  # Component names are shown as rownames
  # Sort mutational signatures by C>T mutation type
  # Sort copy number signatures by BP10MB feature
  is_c2t = grepl("C>T", rownames(sig))
  if (any(is_c2t)) {
    mat_colsum <- colSums(sig[is_c2t, ])
  } else {
    use_M <- any(grepl("\\d+$", rownames(sig)))
    if (use_M) {
    mat_colsum <- colSums(sig[startsWith(rownames(sig), "bp10MB"), ])
    } else {
      # Method: Wang
      mat_colsum <- colSums(sig[startsWith(rownames(sig), "BP10MB"), ])
    }
  }
  sig_order <- rank(mat_colsum)
  return(sig_order)

}
