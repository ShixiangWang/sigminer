#' Show Copy Number Components
#'
#' Show classified components ("Wang" ("W") method) for copy number data.
#'
#' @inheritParams show_cn_features
#' @param parameters a `data.frame` contain parameter components, obtain this
#' from [sig_tally] function.
#' @param show_weights default is `TRUE`, show weights for each component.
#' Only used when method is "Macintyre".
#' @param ... other options pass to \code{\link[cowplot]{plot_grid}} function of **cowplot** package.
#' @inheritParams show_cn_distribution
#' @inheritParams sig_tally
#' @author Shixiang Wang <w_shixiang@163.com>
#' @return a `ggplot` object
#' @export
show_cn_components <- function(parameters,
                               method = "Wang",
                               show_weights = TRUE,
                               log_y = FALSE,
                               return_plotlist = FALSE,
                               base_size = 12,
                               nrow = 2,
                               align = "hv",
                               ...) {
  stopifnot(is.logical(return_plotlist))

  method <- match.arg(method, choices = c("Wang", "W"))

  plot_fun2 <- function(df, f_name, base_size) {
    df$component <- factor(df$component, levels = df$component)
    p <- ggplot(data = df, aes_string(x = "component", y = "n_obs")) +
      geom_bar(stat = "identity", color = "white") +
      labs(x = NULL, y = NULL) +
      theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) +
      cowplot::theme_cowplot(font_size = base_size) +
      theme(
        axis.text.x = element_text(
          angle = 60, vjust = 1,
          hjust = 1,
          color = "black"
        )
      )

    if (log_y) {
      p <- p + scale_y_continuous(trans = "log10")
    }

    p
  }

  p_df <- parameters %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(.data$feature) %>%
    tidyr::nest() %>%
    dplyr::mutate(p_list = purrr::map2(data, .data$feature, plot_fun2, base_size = base_size))
  p_list <- p_df$p_list
  names(p_list) <- p_df$feature

  if (return_plotlist) {
    return(p_list)
  } else {
    if (length(p_list) != 8) {
      p <- cowplot::plot_grid(plotlist = p_list, nrow = nrow, align = align, ...)
    } else {
      top_2_row <- cowplot::plot_grid(plotlist = p_list[1:6], nrow = 2, align = align)
      bot_row <- cowplot::plot_grid(plotlist = p_list[-(1:6)], nrow = 1, align = align, rel_widths = c(1, 2))
      p <- cowplot::plot_grid(top_2_row, bot_row, nrow = 2, align = align, rel_heights = c(2, 1), ...)
    }
    p
  }
}
