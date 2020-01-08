#' Show Copy Number Feature Distributions
#'
#' @param features a feature `list` generate from [sig_tally] function.
#' @param rm_outlier default is `FALSE`, if `TRUE`, remove outliers. Only
#' works when method is "Wang" ("W").
#' @param ylab lab of y axis.
#' @param log_segsize default is `TRUE`, show `log10` based segsize, only
#' works for input from "Macintyre" ("M") method.
#' @param return_plotlist if `TRUE`, return a list of ggplot objects but a combined plot.
#' @param ... other options pass to \code{\link[cowplot]{plot_grid}} function of `cowplot` package.
#' @inheritParams sig_tally
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
show_cn_features <- function(features,
                             method = "Macintyre",
                             rm_outlier = FALSE,
                             ylab = NULL,
                             log_segsize = TRUE,
                             return_plotlist = FALSE,
                             base_size = 12, ...) {
  method <- match.arg(method, choices = c("Macintyre", "M", "Wang", "W"))

  if (startsWith(method, "M")) {
    features <- lapply(features, function(x) {
      x[["value"]] <- as.numeric(x[["value"]])
      return(x)
    })

    if (log_segsize) {
      p_1 <- ggplot(data = features$segsize, aes(x = log10(value))) +
        geom_line(stat = "density") + labs(x = "Segment size (log10 based)", y = ylab) +
        theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + cowplot::theme_cowplot(font_size = base_size)
    } else {
      p_1 <- ggplot(data = features$segsize, aes(x = value)) +
        geom_line(stat = "density") + labs(x = "Segment size", y = ylab) +
        theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + cowplot::theme_cowplot(font_size = base_size)
    }

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
  } else {
    # Method: Wang
    p_list <- purrr::map2(features, names(features), function(feature_df, f_name, rm_outlier = FALSE) {
      if (rm_outlier) {
        outliers <- boxplot(feature_df$value, plot = FALSE)$out
        feature_df <- feature_df[!value %in% outliers]
      }

      ggplot(data = feature_df, aes(x = value)) +
        geom_histogram(binwidth = 1, color = "white") +
        labs(x = f_name, y = ylab) +
        theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) +
        cowplot::theme_cowplot(font_size = base_size)
    }, rm_outlier = rm_outlier)

    if (return_plotlist) {
      return(p_list)
    } else {
      if (length(p_list) != 7) {
        p <- cowplot::plot_grid(plotlist = p_list, nrow = 2, align = "hv", ...)
      } else {
        top_row <- cowplot::plot_grid(plotlist = p_list[1:4], nrow = 1, align = "hv", ...)
        bot_row <- cowplot::plot_grid(plotlist = p_list[-(1:4)], nrow = 1, align = "hv", rel_widths = c(1, 1, 2), ...)
        p <- cowplot::plot_grid(top_row, bot_row, nrow = 2)
      }
      p
    }
  }

  p
}

utils::globalVariables(
  "value"
)
