#' Add Horizontal Arrow with Text Label to a ggplot
#'
#' @param p a `ggplot`.
#' @param x position at x axis.
#' @param y position at y axis.
#' @param label text label.
#' @param space a small space between arrow and text.
#' @param vjust vertical adjustment, set to 0 to align with the bottom,
#' 0.5 for the middle, and 1 (the default) for the top.
#' @param seg_len length of the arrow segment.
#' @param arrow_len length of the arrow.
#' @param arrow_type type of the arrow.
#' @param font_size font size.
#' @param font_family font family.
#' @param font_face font face.
#'
#' @return a `ggplot` object.
#' @export
add_h_arrow <- function(p, x, y, label = "optimal number",
                        space = 0.01,
                        vjust = 0.3,
                        seg_len = 0.1,
                        arrow_len = unit(2, "mm"),
                        arrow_type = c("closed", "open"),
                        font_size = 5,
                        font_family = c("serif", "sans", "mono"),
                        font_face = c("plain", "bold", "italic")) {
  font_family <- match.arg(font_family)
  font_face <- match.arg(font_face)
  arrow_type <- match.arg(arrow_type)

  p + annotate(
    geom = "segment",
    x = x + seg_len + space,
    y = y,
    xend = x + space,
    yend = y,
    arrow = arrow(length = arrow_len, type = arrow_type)
  ) +
    annotate(
      geom = "text",
      size = font_size,
      x = x + seg_len + 2 * space,
      y = y,
      label = label,
      hjust = 0,
      vjust = vjust,
      family = font_family,
      fontface = font_face
    )
}
