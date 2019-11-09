#' Map Groups using Sankey
#'
#' This feature is designed for signature analysis. However, users can also use
#' it in other similar situations.
#'
#' @param data a `data.frame` containing signature group and other categorical groups.
#' @param sig_col length-1 character showing the colname of signature group.
#' @param map_cols character vector showing colnames of other groups.
#' @param include_sig default if `FALSE`, if `TRUE`, showing signature group.
#' @param title the title.
#' @param xlab label for x axis.
#' @param ylab label for y axis.
#' @param custom_theme theme for plotting, default is [cowplot::theme_minimal_hgrid()].
#'
#' @return a `ggplot` object
#' @export
#'
#' @examples
#' data <- dplyr::tibble(
#'   Group1 = rep(LETTERS[1:5], each = 10),
#'   Group2 = rep(LETTERS[6:15], each = 5),
#'   zzzz = c(rep("xx", 20), rep("yy", 20), rep(NA, 10))
#' )
#'
#' show_group_mapping(data, sig_col = "Group1", map_cols = colnames(data)[-1])
#' show_group_mapping(data, sig_col = "Group1", map_cols = colnames(data)[-1], include_sig = TRUE)
show_group_mapping <- function(data, sig_col, map_cols, include_sig = FALSE,
                               title = NULL, xlab = NULL, ylab = NULL,
                               custom_theme = cowplot::theme_minimal_hgrid()) {
  stopifnot(is.data.frame(data), length(sig_col) == 1)

  if (!requireNamespace("ggfittext", quietly = TRUE)) {
    stop("Please install package 'ggfittext' before using this feature!")
  }

  if (include_sig) {
    data <- data %>%
      dplyr::mutate(class_ = .data[[sig_col]]) %>%
      dplyr::select(c(sig_col, map_cols, "class_"))

    data_long <- ggalluvial::to_lodes_form(data, c(sig_col, map_cols))
  } else {
    data <- data %>%
      dplyr::mutate(class_ = .data[[sig_col]]) %>%
      dplyr::select(c(map_cols, "class_"))

    data_long <- ggalluvial::to_lodes_form(data, map_cols)
  }


  # StatStratum <- ggalluvial::StatStratum
  if (!"ggalluvial" %in% .packages()) {
    attachNamespace("ggalluvial")
  }

  ggplot(
    data = data_long,
    aes_string(
      x = "x", alluvium = "alluvium",
      stratum = "stratum", label = "stratum"
    )
  ) +
    scale_x_discrete(expand = c(.1, 0)) +
    ggalluvial::geom_flow(aes_string(fill = "class_"), width = 1 / 4) +
    ggalluvial::geom_stratum(alpha = .5, width = 1 / 4) +
    ggfittext::geom_fit_text(stat = "stratum", width = 1 / 4, min.size = 3) +
    custom_theme +
    theme(legend.position = "none") +
    labs(title = title, x = xlab, y = ylab)
}
