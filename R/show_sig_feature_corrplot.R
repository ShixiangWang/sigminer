#' Draw Corrplot for Signature Exposures and Other Features
#'
#' @param tidy_cor data returned by [get_tidy_association].
#' @param features_list a character vector contains features want to be plotted.
#' @param breaks_count breaks for sample count.
#' @param breaks_score breaks for score.
#'
#' @return a `ggplot2` object
#' @export
#' @seealso [get_tidy_association] and [get_sig_feature_association]
show_sig_feature_corrplot <- function(tidy_cor, features_list,
                                      breaks_count = c(0L, 200L, 400L, 600L, 800L, 1020L),
                                      breaks_score = c(-Inf, -16, -10, -5, -1.3, 1.3, 5, 10, 16, Inf)) {
  data <- tidy_cor$data %>%
    dplyr::mutate(
      Samples = cut(.data$count, breaks = breaks_count),
      Score = cut(.data$score, breaks = breaks_score)
    ) %>%
    dplyr::filter(.data$features %in% features_list, !is.na(.data$score))

  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = signature,
      y = factor(.data$features, levels = rev(features_list))
    )
  ) +
    ggplot2::geom_point(ggplot2::aes(
      colour = .data$Score,
      size = .data$Samples
    )) +
    ggplot2::scale_color_manual(values = rev(RColorBrewer::brewer.pal(9, "RdBu")), drop = FALSE) +
    ggplot2::scale_size_discrete(drop = FALSE) +
    cowplot::theme_cowplot() + ggpubr::rotate_x_text(45) +
    ggplot2::labs(x = "Signatures", y = "Features")
}
