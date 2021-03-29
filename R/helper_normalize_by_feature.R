helper_normalize_by_feature <- function(w, feature_setting = sigminer::CN.features) {
  mat <- as.data.frame(w)
  mat$context <- rownames(mat)

  # Method "W"
  mat$base <- sub("\\[.*\\]$", "", mat$context)
  mat <- tidyr::gather(mat, class, signature, -c("context", "base"))

  if (!inherits(feature_setting, "sigminer.features")) {
    feature_setting <- get_feature_components(feature_setting)
  }
  avail_features <- unique(feature_setting$feature)

  # Make sure
  f_orders <- avail_features[avail_features %in% mat$base]
  comp_orders <- feature_setting$component[feature_setting$feature %in% f_orders]

  mat <- dplyr::mutate(mat,
    context = factor(.data$context,
      levels = comp_orders
    ),
    base = factor(.data$base, levels = f_orders),
    class = factor(class)
  )


  mat <- mat %>%
    dplyr::group_by(.data$base, .data$class) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      context = purrr::map(.data$data, ~ .$context),
      signature = purrr::map(.data$data, ~ .$signature / sum(.$signature))
    ) %>%
    dplyr::select(-"data") %>%
    tidyr::unnest(cols = c("context", "signature")) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(names_from = "class", values_from = "signature") %>%
    dplyr::select(-"base") %>%
    tibble::column_to_rownames(var = "context")

  w <- as.matrix(mat)
  return(w)
}
