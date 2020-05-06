helper_sort_signature <- function(sig) {
  stopifnot(is.matrix(sig), !is.null(rownames(sig)))
  # Component names are shown as rownames
  # Sort mutational signatures by C>T mutation type
  # Sort copy number signatures by segsize/SS feature for method 'M'
  is_c2t <- grepl("C>T", rownames(sig))
  is_cn <- grepl("SS", rownames(sig)) | grepl("segsize", rownames(sig))
  if (any(is_c2t)) {
    to_rank <- colSums(sig[is_c2t, ])
    sig_order <- order(to_rank)
  } else if (any(is_cn)) {
    use_M <- any(grepl("\\d+$", rownames(sig)))
    if (use_M) {
      mat <- sig[startsWith(rownames(sig), "segsize"), , drop = FALSE]
    } else {
      # Method: Wang
      mat <- sig[startsWith(rownames(sig), "SS"), , drop = FALSE]
    }
    sig_order <- get_segsize_order(mat)
  } else {
    ## Keep not change for now
    sig_order <- seq_len(ncol(sig))
  }
  return(sig_order)
}

# Matrix: feature-by-signature
get_segsize_order <- function(mat) {
  index <- apply(mat, 2, which.max)
  is_unique <- length(index) == length(unique(index))
  if (!is_unique) {
    # Compare value with same index
    mat %>%
      t() %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(
        N = dplyr::row_number(),
        index = index
      ) %>%
      dplyr::arrange(.data$index) %>%
      dplyr::group_by(.data$index) %>%
      tidyr::nest() %>%
      dplyr::mutate(data = purrr::map2(.data$data, .data$index, .f = function(x, y) {
        if (nrow(x) == 1) {
          return(x$N)
        } else {
          N <- x$N[order(x[[y]])]
          return(N)
        }
      })) %>%
      tidyr::unnest("data") %>%
      dplyr::pull("data")
  } else {
    return(order(index))
  }
}
