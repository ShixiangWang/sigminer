#' Show Groupped Variable Distribution
#'
#' This is a general function, it can be used in any proper analysis.
#'
#' @param data a `data.frame`.
#' @param gvar a group variable name/index.
#' @param dvar a distribution variable name/index.
#' @param g_label a named vector to set facet labels, default is `NULL`.
#' @param xlab title for x axis.
#' @param ylab title for y axis.
#' @param nrow number of row.
#' @param background_color background color for plot panel.
#'
#' @return a `ggplot` object.
#' @export
#'
#' @examples
#' set.seed(1234)
#' data <- data.frame(
#'   y = rnorm(120),
#'   gr = c(rep("A", 50), rep("B", 40), rep("C", 30))
#' )
#' p <- show_group_distribution(data, gvar = 2, dvar = 1, background_color = "grey")
#' p
#' @testexamples
#' expect_is(p, "ggplot")
show_group_distribution <- function(data, gvar, dvar,
                                    g_label = NULL,
                                    angle_label = 0,
                                    xlab = NULL,
                                    ylab = NULL,
                                    nrow = 1L,
                                    background_color = "#DBD7D2") {
  stopifnot(length(gvar) == 1L, length(dvar) == 1L)

  data$.gvar <- data[[gvar]]
  data$.dvar <- data[[dvar]]
  d <- data %>%
    dplyr::mutate(.order = dplyr::row_number()) %>%
    dplyr::group_by(.data$.gvar) %>%
    dplyr::arrange(.data$.dvar) %>%
    dplyr::mutate(x = seq_len(dplyr::n())) %>%
    dplyr::ungroup()

  ds <- d %>%
    dplyr::group_by(.data$.gvar) %>%
    dplyr::summarise(
      x_m = median(.data$x, na.rm = TRUE),
      y_m = median(.data$.dvar, na.rm = TRUE),
      n = sum(!is.na(.data$.dvar)),
      .groups = "drop"
    ) %>%
    dplyr::transmute(
      g = .data$.gvar,
      x = .data$x_m - round(.data$n / 3),
      xend = .data$x_m + round(.data$n / 3),
      y = .data$y_m,
      yend = .data$y_m,
      label = paste0(.data$g, "\n(n=", .data$n, ")")
    )

  g_label <- ds$label
  names(g_label) <- ds$g

  p <- ggplot(d, aes_string(x = "x", y = ".dvar")) +
    geom_point(alpha = 0.8) +
    geom_segment(aes_string(x = "x", xend = "xend", y = "y", yend = "yend"),
      data = ds,
      color = "red",
      size = 2
    ) +
    facet_wrap(~.gvar,
      nrow = nrow,
      scales = "free_x",
      labeller = labeller(g = g_label)
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      strip.background.x = element_rect(color = "white", fill = "white"),
      strip.text.x = element_text(
        color = "black",
        face = "bold"
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0, "line"),
      panel.background = element_rect(fill = background_color)
    ) +
    labs(x = xlab, y = ylab)
  p$sampleOrder <- d$.order
  p

}
