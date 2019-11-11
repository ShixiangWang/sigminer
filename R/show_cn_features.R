#' Show Copy Number Feature Distributions
#'
#' @param features a feature `list` generate from [sig_derive] function.
#' @param ylab lab of y axis.
#' @param log_segsize default is `TRUE`, show `log10` based segsize.
#' @param return_plotlist if `TRUE`, return a list of ggplot objects but a combined plot.
#' @param ... other options pass to \code{\link[cowplot]{plot_grid}} function of `cowplot` package.
#' @inheritParams show_cn_distribution
#' @return a `ggplot` object
#' @examples
#' # Load copy number prepare object
#' load(system.file("extdata", "toy_copynumber_prepare.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' show_cn_features(cn_prepare$features)
#' @export
#'
show_cn_features <- function(features, ylab = NULL,
                             log_segsize = TRUE,
                             return_plotlist = FALSE,
                             base_size = 12, ...) {
  features <- lapply(features, function(x) {
    x[["value"]] <- as.numeric(x[["value"]])
    return(x)
  })

  # requireNamespace("cowplot")
  # cowplot::theme_cowplot()
  # ggplot2::theme_set(cowplot::theme_cowplot(font_size = 12))
  if (log_segsize) {
    p_1 <- ggplot(data = features$segsize, aes(x = log10(value))) +
      geom_line(stat = "density") + labs(x = "Segment size (log10 based)", y = ylab) +
      theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + cowplot::theme_cowplot(font_size = base_size)
  } else {
    p_1 <- ggplot(data = features$segsize, aes(x = value)) +
      geom_line(stat = "density") + labs(x = "Segment size", y = ylab) +
      theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + cowplot::theme_cowplot(font_size = base_size)
  }

  # p_1 = p_1 + scale_x_continuous(breaks = 10 ^c(7, 8),
  #                               labels = scales::trans_format("log10", scales::math_format(10^.x)))

  p_2 <- ggplot(data = features$copynumber, aes(x = value)) +
    geom_line(stat = "density") + labs(x = "Copy number", y = ylab) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + cowplot::theme_cowplot(font_size = base_size)

  p_3 <- ggplot(data = features$changepoint, aes(x = value)) +
    geom_line(stat = "density") + labs(x = "Copy number changepoint", y = ylab) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + cowplot::theme_cowplot(font_size = base_size)

  p_4 <- ggplot(data = features$bp10MB, aes(x = value)) +
    geom_bar(stat = "count") + labs(x = "Breakpoint count per 10MB", y = ylab) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + cowplot::theme_cowplot(font_size = base_size)
  p_5 <- ggplot(data = features$bpchrarm, aes(x = value)) +
    geom_bar(stat = "count") + labs(x = "Breakpoint count per chr arm", y = ylab) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + cowplot::theme_cowplot(font_size = base_size)

  p_6 <- ggplot(data = features$osCN, aes(x = value)) +
    geom_bar(stat = "count") + labs(x = "Oscilating CN chain length", y = ylab) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + cowplot::theme_cowplot(font_size = base_size)

  if (return_plotlist) {
    return(list(
      segsize = p_1, copynumber = p_2, changepoint = p_3,
      bp10MB = p_4, bpchrarm = p_5, osCN = p_6
    ))
  }

  p <- cowplot::plot_grid(p_1,
    p_2,
    p_3,
    p_4,
    p_5,
    p_6,
    nrow = 2,
    align = "hv",
    ...
  )
  p
}

utils::globalVariables(
  "value"
)
