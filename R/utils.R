
# Source code from reporttools (https://github.com/cran/reporttools/blob/master/R/pairwise.fisher.test.r)
pairwise.fisher.test <- function(x, g, p.adjust.method, ...) {
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
  g <- factor(g)

  compare.levels <- function(i, j) {
    ind <- (as.integer(g) %in% c(i, j)) & (is.na(x) == FALSE) & (is.na(g) == FALSE)
    xi <- factor(x[ind], exclude = NULL)
    xj <- factor(g[ind], exclude = NULL)
    tab <- table(xi, xj)
    nonzeromarginal <- (min(apply(tab, 1, sum)) * min(apply(tab, 2, sum)) > 0)
    size <- ((nrow(tab) > 1) * (ncol(tab) > 1) > 0)
    if ((nonzeromarginal == TRUE) & (size == TRUE)) {
      fisher.test(xi, xj, ...)$p.value
    } else {
      NA
    }
  }

  PVAL <- pairwise.table(compare.levels, levels(g), p.adjust.method)
  ans <- list(data.name = DNAME, p.value = PVAL, p.adjust.method = p.adjust.method)
  class(ans) <- "pairwise.htest"
  return(ans)
}


# https://stackoverflow.com/questions/8396577/check-if-character-value-is-a-valid-r-object-name
isValidAndUnreserved <- function(string) {
  make.names(string) == string
}

`%:::%` <- function(pkg, fun) {
  get(fun,
    envir = asNamespace(pkg),
    inherits = FALSE
  )
}

# From https://gist.github.com/mbannert/e9fcfa86de3b06068c83
rgb2hex <- function(r, g, b) grDevices::rgb(r, g, b, maxColorValue = 255)
col2hex <- function(col, alpha) grDevices::rgb(t(grDevices::col2rgb(col)), alpha = alpha, maxColorValue = 255)


#' Set Color Style for Plotting
#'
#' @param style one of 'default' and 'cosmic'.
#'
#' @return color values.
#' @export
#'
#' @examples
#' use_color_style("default")
#' use_color_style("cosmic")
use_color_style <- function(style) {
  #c("red", "cyan", "yellow", "blue", "magenta", "gray50", "orange", "darkgreen", "brown", "black", rainbow(10)[4:10])
  if (style == "default") {
    palette <- c(
      "cyan", "red", "brown", "purple",
      "green", "blue", "black", "magenta",
      "orange", "gray50", "darkgreen", "yellow",
      rainbow(10)[4:10]
    )
  } else {
    colors <- list(
      c(3, 189, 239),
      c(1, 1, 1),
      c(228, 41, 38),
      c(203, 202, 202),
      c(162, 207, 99),
      c(236, 199, 197)
    )

    # https://github.com/AlexandrovLab/SigProfilerPlotting/blob/a60c48c2b74f037b9df5a3f905cc6ada22623341/sigProfilerPlotting/sigProfilerPlotting.py#L96
    # colors = [[3/256,189/256,239/256], [1/256,1/256,1/256],
    #           [228/256,41/256,38/256], [203/256,202/256,202/256], [162/256,207/256,99/256], [236/256,199/256,197/256]]

    palette <- c(
      sapply(colors, FUN = function(x) rgb2hex(x[1], x[2], x[3])) %>% as.character(),
      sapply(c("purple", "brown", "orange"), FUN = col2hex) %>% as.character()
    )
  }
  palette
}
