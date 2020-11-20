#' Show Groupped Variable Distribution
#'
#' This is a general function, it can be used in any proper analysis.
#'
#' @param data a `data.frame`.
#' @param gvar a group variable name/index.
#' @param dvar a distribution variable name/index.
#' @param fun a function to summarize, default is [stats::median], can also be [mean].
#' @param order_by_fun if `TRUE`, reorder the groups by summary measure computed
#' by argument `fun`.
#' @param alpha alpha for points, range from 0 to 1.
#' @param g_label a string 'label' (default) for labeling with sample size,
#' or 'norm' to show just group name, or a named vector to set facet labels.
#' @param g_angle angle for facet labels, default is `0`.
#' @param g_position position for facet labels, default is 'top', can also
#' be 'bottom'.
#' @param point_size size of point.
#' @param segment_size size of segment.
#' @param segment_color color of segment.
#' @param xlab title for x axis.
#' @param ylab title for y axis.
#' @param nrow number of row.
#' @param background_color background color for plot panel.
#'
#' @author Shixiang Wang <w_shixiang@163.com>
#' @return a `ggplot` object.
#' @export
#'
#' @examples
#' set.seed(1234)
#' data <- data.frame(
#'   yval = rnorm(120),
#'   gr = c(rep("A", 50), rep("B", 40), rep("C", 30))
#' )
#' p <- show_group_distribution(data,
#'   gvar = 2, dvar = 1,
#'   g_label = "norm",
#'   background_color = "grey"
#' )
#' p
#' p2 <- show_group_distribution(data,
#'   gvar = "gr", dvar = "yval",
#'   g_position = "bottom",
#'   order_by_fun = TRUE,
#'   alpha = 0.3
#' )
#' p2
#'
#' # Set custom group names
#' p3 <- show_group_distribution(data,
#'   gvar = 2, dvar = 1,
#'   g_label = c("A" = "X", "B" = "Y", "C" = "Z")
#' )
#' p3
#' @testexamples
#' expect_is(p, "ggplot")
#' expect_is(p2, "ggplot")
#' expect_is(p3, "ggplot")
show_group_distribution <- function(data, gvar, dvar,
                                    fun = stats::median,
                                    order_by_fun = FALSE,
                                    alpha = 0.8,
                                    g_label = "label",
                                    g_angle = 0,
                                    g_position = "top",
                                    point_size = 1L,
                                    segment_size = 1L,
                                    segment_color = "red",
                                    xlab = NULL,
                                    ylab = NULL,
                                    nrow = 1L,
                                    background_color = c("#DCDCDC", "#F5F5F5")) {
  stopifnot(length(gvar) == 1L, length(dvar) == 1L)

  data$.gvar <- data[[gvar]]
  data$.dvar <- data[[dvar]]
  d <- data %>%
    dplyr::mutate(.order = dplyr::row_number()) %>%
    dplyr::group_by(.data$.gvar) %>%
    dplyr::arrange(.data$.gvar, .data$.dvar) %>%
    dplyr::mutate(x = seq_len(dplyr::n())) %>%
    dplyr::ungroup()

  ds <- d %>%
    dplyr::group_by(.data$.gvar) %>%
    dplyr::summarise(
      x_m = median(.data$x, na.rm = TRUE),
      y_m = fun(.data$.dvar, na.rm = TRUE),
      n = sum(!is.na(.data$.dvar)),
      .groups = "drop"
    ) %>%
    dplyr::transmute(
      .gvar = .data$.gvar,
      x = ifelse(.data$n > 3, .data$x_m - round(.data$n / 3),
        .data$x_m - .data$n / 3
      ),
      xend = ifelse(.data$n > 3, .data$x_m + round(.data$n / 3),
        .data$x_m + .data$n / 3
      ),
      y = .data$y_m,
      yend = .data$y_m,
      label = paste0(.data$.gvar, "\n(n=", .data$n, ")")
    )

  ## Use rect fill colors to set panel colors
  dp <- d %>%
    dplyr::group_by(.data$.gvar) %>%
    dplyr::summarise(
      n = dplyr::n(),
      xmin = min(.data$x, na.rm = TRUE),
      xmax = max(.data$x, na.rm = TRUE),
      ymin = min(.data$.dvar, na.rm = TRUE),
      ymax = max(.data$.dvar, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      xmin = ifelse(.data$n > 3, .data$xmin - (.data$xmax - .data$xmin) * 0.05,
        .data$xmin - 0.5
      ),
      xmax = ifelse(.data$n > 3, .data$xmax + (.data$xmax - .data$xmin) * 0.05,
        .data$xmax + 0.5
      ),
      ymin = min(.data$ymin) - (.data$ymax - .data$ymin) * 0.015,
      ymax = max(.data$ymax) + (.data$ymax - .data$ymin) * 0.015,
      ymin = min(.data$ymin),
      ymax = max(.data$ymax)
    )

  if (order_by_fun) {
    ds <- ds %>%
      dplyr::arrange(.data$y) %>%
      dplyr::mutate(.gvar = factor(.data$.gvar, levels = .data$.gvar))
    d$.gvar <- factor(d$.gvar, levels = levels(ds$.gvar))
    dp$.gvar <- factor(dp$.gvar, levels = levels(ds$.gvar))
  }

  b_colors <- rep(background_color, ceiling(nrow(dp) / length(background_color)))

  if (is.factor(d$.gvar)) {
    dp <- dp %>% dplyr::arrange(.data$.gvar)
  }

  dp$b_colors <- b_colors[1:nrow(dp)]
  dp$b_colors <- factor(dp$b_colors, levels = background_color)

  if (identical(g_label, "label")) {
    g_label <- ds$label
    names(g_label) <- ds$.gvar
  } else if (identical(g_label, "norm")) {
    g_label <- ds$.gvar
    names(g_label) <- ds$.gvar
  }

  p <- ggplot(d) +
    geom_rect(aes_string(
      xmin = "xmin", xmax = "xmax",
      ymin = "ymin", ymax = "ymax",
      fill = "b_colors", color = "b_colors"
    ),
    data = dp
    ) +
    geom_point(aes_string(x = "x", y = ".dvar"),
      alpha = alpha, data = d, size = point_size
    ) +
    geom_segment(aes_string(x = "x", xend = "xend", y = "y", yend = "yend"),
      data = ds,
      color = segment_color,
      size = segment_size
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    scale_fill_manual(values = background_color) +
    scale_color_manual(values = background_color) +
    facet_wrap(~.gvar,
      nrow = nrow,
      scales = "free_x",
      strip.position = g_position,
      labeller = labeller(.gvar = g_label)
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "none",
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      strip.background.x = element_rect(color = "white", fill = "white"),
      strip.text.x = element_text(
        color = "black",
        face = "bold",
        angle = g_angle
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0, "line")
      # panel.background = element_rect(fill = background_color)
    ) +
    labs(x = xlab, y = ylab)
  p$sampleOrder <- d$.order
  p
}
