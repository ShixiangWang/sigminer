#' Map Groups using Sankey
#'
#' This feature is designed for signature analysis. However, users can also use
#' it in other similar situations.
#'
#' @param data a `data.frame` containing signature group and other categorical groups.
#' @param col_to_flow length-1 character showing the column to flow, typically a signature group.
#' @param cols_to_map character vector showing colnames of other groups.
#' @param include_sig default if `FALSE`, if `TRUE`, showing signature group.
#' @param fill_na length-1 string to fill NA, default is `FALSE`.
#' @param title the title.
#' @param xlab label for x axis.
#' @param ylab label for y axis.
#' @param custom_theme theme for plotting, default is `cowplot::theme_minimal_hgrid()`.
#'
#' @return a `ggplot` object
#' @export
#'
#' @examples
#' \donttest{
#' data <- dplyr::tibble(
#'   Group1 = rep(LETTERS[1:5], each = 10),
#'   Group2 = rep(LETTERS[6:15], each = 5),
#'   zzzz = c(rep("xx", 20), rep("yy", 20), rep(NA, 10))
#' )
#' p1 <- show_group_mapping(data, col_to_flow = "Group1", cols_to_map = colnames(data)[-1])
#' p1
#'
#' p2 <- show_group_mapping(data,
#'   col_to_flow = "Group1", cols_to_map = colnames(data)[-1],
#'   include_sig = TRUE
#' )
#' p2
#' }
#' @testexamples
#' if (!is.null(p1)) {
#' expect_s3_class(p1, "ggplot")
#' }
#' if (!is.null(p2)) {
#' expect_s3_class(p2, "ggplot")
#' }
show_group_mapping <- function(data, col_to_flow, cols_to_map, include_sig = FALSE,
                               fill_na = FALSE,
                               title = NULL, xlab = NULL, ylab = NULL,
                               custom_theme = cowplot::theme_minimal_hgrid()) {
  stopifnot(is.data.frame(data), length(col_to_flow) == 1)

  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    message("Please install package 'ggalluvial' before using this feature!")
    return(invisible(NULL))
  }

  if (!requireNamespace("ggfittext", quietly = TRUE)) {
    message("Please install package 'ggfittext' before using this feature!")
    return(invisible(NULL))
  }

  if (include_sig) {
    data <- data %>%
      dplyr::mutate(class_ = .data[[col_to_flow]]) %>%
      dplyr::select(c(col_to_flow, cols_to_map, "class_"))

    if (!isFALSE(fill_na)) {
      data <- data %>%
        dplyr::mutate_all(dplyr::funs(ifelse(is.na(.), tidyr::replace_na(., fill_na), .)))
    }

    data_long <- ggalluvial::to_lodes_form(data, c(col_to_flow, cols_to_map))
  } else {
    data <- data %>%
      dplyr::mutate(class_ = .data[[col_to_flow]]) %>%
      dplyr::select(c(cols_to_map, "class_"))

    if (!isFALSE(fill_na)) {
      data <- data %>%
        dplyr::mutate_all(dplyr::funs(ifelse(is.na(.), tidyr::replace_na(., fill_na), .)))
    }

    data_long <- ggalluvial::to_lodes_form(data, cols_to_map)
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
