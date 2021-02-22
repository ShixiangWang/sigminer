#' Show Copy Number Feature Distributions
#'
#' @param features a feature `list` generate from [sig_tally] function.
#' @param rm_outlier default is `FALSE`, if `TRUE`, remove outliers. Only
#' works when method is "Wang" ("W").
#' @param ylab lab of y axis.
#' @param log_y logical, if `TRUE`, show `log10` based y axis, only
#' works for input from "Wang" ("W") method.
#' @param return_plotlist if `TRUE`, return a list of ggplot objects but a combined plot.
#' @param ... other options pass to \code{\link[cowplot]{plot_grid}} function of `cowplot` package.
#' @inheritParams sig_tally
#' @inheritParams show_cn_distribution
#' @inheritParams cowplot::plot_grid
#' @return a `ggplot` object
#' @export
#'
show_cn_features <- function(features,
                             method = "Wang",
                             rm_outlier = FALSE,
                             ylab = NULL,
                             log_y = FALSE,
                             return_plotlist = FALSE,
                             base_size = 12,
                             nrow = 2,
                             align = "hv",
                             ...) {
  method <- match.arg(method, choices = c("Wang", "W", "Tao & Wang", "T"))

  # Method: Wang or Tao & Wang
  p_list <- purrr::map2(features, names(features), function(feature_df, f_name, rm_outlier = FALSE, log_y = FALSE) {
    if (rm_outlier) {
      outliers <- boxplot(feature_df$value, plot = FALSE)$out
      feature_df <- feature_df[!value %in% outliers]
    }

    p <- ggplot(data = feature_df, aes(x = value)) +
      geom_histogram(binwidth = 1, color = "white") +
      labs(x = f_name, y = ylab) +
      theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) +
      cowplot::theme_cowplot(font_size = base_size)

    if (log_y) {
      p <- p + scale_y_continuous(trans = "log10")
    }

    p
  }, rm_outlier = rm_outlier, log_y = log_y)

  if (return_plotlist) {
    return(p_list)
  } else {
    if (startsWith(method, "W") & length(p_list) == 8) {
      top_2_row <- cowplot::plot_grid(plotlist = p_list[1:6], nrow = 2, align = align)
      bot_row <- cowplot::plot_grid(plotlist = p_list[-(1:6)], nrow = 1, align = align, rel_widths = c(1, 2))
      p <- cowplot::plot_grid(top_2_row, bot_row, nrow = 2, align = align, rel_heights = c(2, 1), ...)
    } else {
      p <- cowplot::plot_grid(plotlist = p_list, nrow = nrow, align = align, ...)
    }
  }


  p
}

utils::globalVariables(
  "value"
)
