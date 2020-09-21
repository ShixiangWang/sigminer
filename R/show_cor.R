#' A Simple and General Way for Association Analysis
#'
#' All variables must be continuous.
#' The matrix will be returned as an element of `ggplot` object.
#' This is basically a wrapper of R package
#' [ggcorrplot](https://github.com/kassambara/ggcorrplot).
#'
#' @param data a `data.frame`.
#' @param x_vars variables/column names shown in x axis.
#' @param y_vars variables/column names shown in y axis.
#' @param cor_method method for correlation, default is 'spearman'.
#' @param vis_method visualization method, default is 'square',
#' can also be 'circle'.
#' @param lab logical value. If TRUE, add correlation coefficient on the plot.
#' @param hc_order logical value. If `TRUE`,
#' correlation matrix will be hc.ordered using `hclust` function.
#' @param p_adj p adjust method, see [stats::p.adjust] for details.
#' @param ... other parameters passing to `ggcorrplot::ggcorrplot()`.
#'
#' @return a `ggplot` object
#' @export
#' @seealso [show_sig_feature_corrplot] for specific and more powerful
#' association analysis and visualization.
#'
#' @examples
#' data("mtcars")
#' p1 <- show_cor(mtcars)
#' p2 <- show_cor(mtcars,
#'                x_vars = colnames(mtcars)[1:4],
#'                y_vars = colnames(mtcars)[5:8])
#' p3 <- show_cor(mtcars, vis_method = "circle", p_adj = "fdr")
#' p1
#' p1$cor
#' p2
#' p3
#' @testexamples
#' expect_is(p1, "ggplot")
#' expect_is(p2, "ggplot")
#' expect_is(p3, "ggplot")
show_cor <- function(data, x_vars = colnames(data), y_vars = x_vars,
                     cor_method = "spearman",
                     vis_method = "square",
                     lab = TRUE,
                     hc_order = FALSE,
                     p_adj = NULL,
                     ...) {
  stopifnot(is.data.frame(data), is.character(x_vars), is.character(y_vars))
  if (!requireNamespace("ggcorrplot", quietly = TRUE)) {
    message("'ggcorrplot' package is required to use this feature.")
    return(invisible(NULL))
  }

  all_vars <- union(x_vars, y_vars)
  data <- data %>%
    dplyr::as_tibble() %>%
    dplyr::select(all_vars)

  corr <- round(stats::cor(data,
                           use = "pairwise.complete.obs",
                           method = cor_method), 2)
  p_mat <- ggcorrplot::cor_pmat(data)

  if (!is.null(p_adj)) {
    pa <- stats::p.adjust(p_mat, method = p_adj)
    p_mat <- matrix(pa, nrow = nrow(p_mat), ncol = ncol(p_mat),
                    dimnames = list(rownames(p_mat),
                                    colnames(p_mat)))
  }

  corr <- corr[x_vars, y_vars, drop = FALSE]
  p_mat <- p_mat[x_vars, y_vars, drop = FALSE]

  p <- ggcorrplot::ggcorrplot(
    corr,
    method = vis_method,
    hc.order = hc_order,
    lab = TRUE,
    p.mat = p_mat,
    ...
  )
  p$cor <- list()
  p$cor$cor_mat <- corr
  p$cor$p_mat <- p_mat

  return(p)
}
