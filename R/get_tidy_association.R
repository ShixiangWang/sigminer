#' Get Tidy Signature Association Results
#'
#' @param cor_res data returned by [get_sig_feature_association()]
#'
#' @return a `data.frame`
#' @export
#' @seealso [get_sig_feature_association]
get_tidy_association <- function(cor_res) {
  co_null <- identical(cor_res$corr_co, list())
  ca_null <- identical(cor_res$corr_ca, list())

  if (all(co_null, ca_null)) {
    message("No correlation result found, please check the previous step!")
    return(invisible())
  }

  if (co_null) {
    all_names <- names(cor_res$corr_ca)
  } else {
    all_names <- names(cor_res$corr_co)
  }

  common_names <- setdiff(all_names, "data")
  res <- list()
  tidy <- list()
  for (i in common_names) {
    res[[i]] <- cbind(cor_res[["corr_co"]][[i]], cor_res[["corr_ca"]][[i]])
    tidy[[i]] <- res[[i]] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("signature") %>%
      tibble::as_tibble() %>%
      tidyr::gather_(key_col = "feature", value_col = "value", gather_cols = colnames(res[[i]]))
  }

  data <- dplyr::bind_rows(tidy, .id = "data_type") %>%
    tidyr::spread_(key_col = "data_type", value_col = "value") %>%
    dplyr::mutate(
      type = ifelse(.data$feature %in% colnames(cor_res$corr_ca$measure), "ca", "co"),
      feature = factor(.data$feature, levels = c(
        colnames(cor_res$corr_co$measure),
        colnames(cor_res$corr_ca$measure)
      ))
    )

  data
}
