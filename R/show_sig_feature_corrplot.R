#' Draw Corrplot for Signature Exposures and Other Features
#'
#' @param tidy_cor data returned by [get_tidy_association].
#' @param feature_list a character vector contains features want to be plotted.
#' If missing, all features will be used.
#' @param sort_features default is `FALSE`, use feature order obtained from the previous
#' step. If `TRUE`, sort features as `feature_list`.
#' @param return_plotlist if `TRUE`, return as a list of `ggplot` objects.
#' @param p_val p value threshold.
#' @param xlab label for x axis.
#' @param ylab label for y axis.
#' @param gradient_colors gradient colors used to plot.
#' @param align one of "row" and "column".
#' @param plot_ratio a length-2 numeric vector to set the height/width ratio.
#' @param breaks_count breaks for sample count.
#'
#' @return a `ggplot2` object
#' @export
#' @seealso [get_tidy_association] and [get_sig_feature_association]
show_sig_feature_corrplot <- function(tidy_cor, feature_list,
                                       sort_features = FALSE,
                                       return_plotlist = FALSE,
                                       p_val = 0.05,
                                       xlab = "Signatures", ylab = "Features",
                                       gradient_colors = c("blue", "white", "red"),
                                       align = c("row", "column"), plot_ratio = "auto",
                                       breaks_count = c(0L, 200L, 400L, 600L, 800L, 1020L)) {
  align <- match.arg(align)
  if (!is.character(plot_ratio)) {
    if (length(plot_ratio) != 2 | !is.numeric(plot_ratio)) {
      stop("plot_ratio must be a length-2 numeric vector!")
    }
  }

  if (missing(feature_list)) {
    feature_list <- tidy_cor$feature %>% unique()
  }

  data <- tidy_cor %>%
    dplyr::mutate(
      Samples = cut(.data$count, breaks = breaks_count)
    ) %>%
    dplyr::filter(.data$feature %in% feature_list,
                  .data$p <= p_val)

  .plot_cor <- function(data) {
    if (sort_features) {
      p <- ggplot2::ggplot(
        data,
        ggplot2::aes(
          x = signature,
          y = factor(.data$feature, levels = rev(feature_list))
        )
      )
    } else {
      p <- ggplot2::ggplot(
        data,
        ggplot2::aes(
          x = signature,
          y = .data$feature
        )
      )
    }

    p + ggplot2::geom_point(ggplot2::aes_string(
      colour = "measure",
      size = "Samples"
    )) +
      ggplot2::scale_color_gradientn(colors = gradient_colors) +
      ggplot2::scale_size_discrete(drop = FALSE)
  }

  data2 <- data
  data <- data %>%
    dplyr::group_by(.data$type) %>%
    tidyr::nest() %>%
    dplyr::mutate(gg = purrr::map(.data$data, .plot_cor)) %>%
    dplyr::ungroup()


  gglist <- data$gg
  names(gglist) <- data$type

  gglist <- purrr::map2(gglist, names(gglist), function(p, type) {
    p <- p +
      cowplot::theme_cowplot() +
      ggpubr::rotate_x_text(45) +
      ggplot2::labs(x = xlab, y = ylab)

    if (type == "co") {
      p <- p + labs(color = "Correlation\ncoefficient")
    } else {
      p <- p + labs(color = "Difference\nin means of\nexposure")
    }
  })

  if (length(gglist) > 1) {
    if (return_plotlist) {
      return(gglist)
    }

    if (is.character(plot_ratio)) {
      count_df <- data2 %>%
        dplyr::mutate(type = factor(.data$type, levels = c("co", "ca"))) %>%
        dplyr::group_by(.data$type) %>%
        dplyr::summarise(count = length(unique(.data$feature)))


      if (align == "row") {
        heights <- count_df$count
        widths <- NULL
      } else {
        heights <- NULL
        wdiths <- count_df$count
      }

      co <- gglist$co + labs(x = NULL) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
      ca <- gglist$ca + guides(size = FALSE) + labs(y = NULL)

      co + ca +
        plot_layout(
          byrow = align == "row",
          heights = heights,
          widths = widths
        )
    } else {
      if (align == "row") {
        heights <- plot_ratio
        widths <- NULL
      } else {
        heights <- NULL
        wdiths <- plot_ratio
      }

      co <- gglist$co + labs(x = NULL) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
      ca <- gglist$ca + ggplot2::guides(size = FALSE) + labs(y = NULL)

      co + ca +
        plot_layout(
          byrow = align == "row",
          heights = heights,
          widths = widths
        )
    }
  } else {
    return(gglist[[1]])
  }
}
