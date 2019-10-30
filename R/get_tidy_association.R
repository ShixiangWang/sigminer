#' Get Tidy Signature Association Results
#'
#' @param cor_res data returned by [get_sig_feature_association()]
#'
#' @return a `data.frame`
#' @export
#' @seealso [get_sig_feature_association]
get_tidy_association <- function(cor_res) {
  all_names <- names(cor_res[["corr_ca"]])
  common_names <- setdiff(all_names, "data")
  res <- list()
  tidy <- list()
  for (i in common_names) {
    if (!identical(rownames(cor_res[["corr_co"]][[i]]), rownames(cor_res[["corr_ca"]][[i]]))) {
      stop("Bad rowname order, please report this issue to developer with your sample data")
    }
    res[[i]] <- cbind(cor_res[["corr_co"]][[i]], cor_res[["corr_ca"]][[i]])
    tidy[[i]] <- res[[i]] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("signature") %>%
      tibble::as_tibble() %>%
      tidyr::gather_(key_col = "features", value_col = "value", gather_cols = colnames(res[[i]]))
  }

  # To avoid Inf,
  # set p == 0 to the same value
  # as p value among p != 0
  res[["data"]] <- dplyr::bind_rows(tidy, .id = "data_type") %>%
    tidyr::spread_(key_col = "data_type", value_col = "value") %>%
    dplyr::mutate(
      score = ifelse(p == 0, -log10(min(p[p != 0], na.rm = TRUE)), -log10(p)),
      score = ifelse(measure < 0, -score, score)
    )
  res
}
