helper_sort_signature <- function(sig) {
  stopifnot(is.matrix(sig), !is.null(rownames(sig)))
  # Component names are shown as rownames
  # Sort mutational signatures by C>T mutation type
  # Sort copy number signatures by fraction of the first component of BP10MB feature for method 'M'
  # Sort copy number signatures by fraction of the copy number value 2 for method 'W'
  is_c2t <- grepl("C>T", rownames(sig))
  if (any(is_c2t)) {
    to_rank <- colSums(sig[is_c2t, ])
  } else {
    use_M <- any(grepl("\\d+$", rownames(sig)))
    if (use_M) {
      mat <- sig[startsWith(rownames(sig), "bp10MB"), , drop = FALSE]
      to_rank <- mat[1, ] / colSums(mat)
    } else {
      # Method: Wang
      mat <- sig[startsWith(rownames(sig), "CN"), , drop = FALSE]
      to_rank <- mat[3, ] / colSums(mat)
    }
  }
  sig_order <- order(to_rank)
  return(sig_order)
}
