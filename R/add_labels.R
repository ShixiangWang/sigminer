#' Add Text Labels to a ggplot
#'
#' Add text labels to a ggplot object, such as the result
#' from [show_sig_profile].
#'
#' @inheritParams add_h_arrow
#' @param y_end end position of y axis when `n_label` is set.
#' @param n_label the number of label, when this is set,
#' the position of labels at y axis is auto-generated
#' according to `y` and `y_end`.
#' @param labels text labels or a `similarity` object from [get_sig_similarity].
#' @param ... other parameters passing to [ggplot2::annotate].
#'
#' @return a `ggplot` object.
#' @export
#'
#' @examples
#' # Load mutational signature
#' load(system.file("extdata", "toy_mutational_signature.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Show signature profile
#' p <- show_sig_profile(sig2, mode = "SBS")
#'
#' # Method 1
#' p1 <- add_labels(p,
#'   x = 0.75, y = 0.3, y_end = 0.9, n_label = 3,
#'   labels = paste0("text", 1:3)
#' )
#' p1
#'
#' # Method 2
#' p2 <- add_labels(p,
#'   x = c(0.15, 0.6, 0.75), y = c(0.3, 0.6, 0.9),
#'   labels = paste0("text", 1:3)
#' )
#' p2
#'
#' # Method 3
#' sim <- get_sig_similarity(sig2)
#' p3 <- add_labels(p,
#'   x = c(0.15, 0.6, 0.75), y = c(0.3, 0.6, 0.9),
#'   labels = sim, font_size = 2
#' )
#' p3
#' @testexamples
#' expect_s3_class(p, "ggplot")
#' expect_s3_class(p1, "ggplot")
#' expect_s3_class(p2, "ggplot")
#' expect_s3_class(p3, "ggplot")
add_labels <- function(p, x, y, y_end = NULL,
                       n_label = NULL,
                       labels = NULL,
                       font_size = 5,
                       font_family = "serif",
                       font_face = c("plain", "bold", "italic"),
                       ...) {
  stopifnot(!is.null(labels))

  if (!is.null(n_label)) {
    if (is.null(y_end)) {
      stop("When n_label is set, please also set y_end.", call. = FALSE)
    }
    y <- seq(from = y, to = y_end, length.out = n_label)
  }

  if (inherits(labels, "similarity")) {
    labels <- dplyr::bind_rows(labels$best_match)
    labels <- paste(labels$best_match, labels$aetiology, sep = "\n")
  }

  font_face <- match.arg(font_face)

  p + annotate(
    geom = "text",
    size = font_size,
    x = x,
    y = y,
    label = labels,
    family = font_family,
    fontface = font_face,
    ...
  )
}
