helper_normalize_by_feature <- function(w, feature_setting = sigminer::CN.features) {
  # Detect the pattern of component names
  # use_W = any(grepl("\\[.*\\]$", rownames(w)))
  use_M <- any(grepl("\\d+$", rownames(w)))

  mat <- as.data.frame(w)
  mat$context <- rownames(mat)

  if (use_M) {
    # Method "M"
    mat$base <- sub("\\d+$", "", mat$context)
    mat <- tidyr::gather(mat, class, signature, -c("context", "base"))

    mat <- mat %>%
      dplyr::mutate(
        base = factor(.data$base, levels = c(
          "bp10MB", "bpchrarm",
          "copynumber", "changepoint",
          "osCN", "segsize"
        ))
      ) %>%
      dplyr::arrange(.data$base)

    mat <- dplyr::mutate(mat,
      context = factor(.data$context,
        levels = unique(mat[["context"]])
      ),
      class = factor(class)
    )
  } else {
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
  }

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
