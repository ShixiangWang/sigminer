# Plot copy number distribution -------------------------------------------


#' Plot copy number distribution either by length or chromosome
#'
#' Visually summarize copy number distribution either by copy number segment length
#' or chromosome. When input is a [CopyNumber] object, `genome_build` option will
#' read from `genome_build` slot of object instead of using argument set in function.
#'
#' @param data result from [get_LengthFraction] function result or a [CopyNumber] object.
#' @param rm_normal logical. Whether remove normal copy (i.e. "segVal" equals 2), default is `TRUE`.
#' @param mode either "ld" for distribution by CN length or "cd" for distribution by chromosome.
#' @param fill when `mode` is "cd" and `fill` is `TRUE`, plot percentage instead of count.
#' @param scale_chr logical. If `TRUE`, normalize count to per Megabase unit.
#' @inheritParams get_features
#' @author Shixiang Wang <w_shixiang@163.com>
#' @return a ggplot object
#' @import ggplot2
#' @export
#' @family copy number plot
draw_cn_distribution = function(data,
                                rm_normal = TRUE,
                                mode = c("ld", "cd"),
                                fill = FALSE,
                                scale_chr = TRUE,
                                genome_build = c("hg19", "hg38")) {
  stopifnot(is.logical(rm_normal),
            is.data.frame(data) | inherits(data, "CopyNumber"),
            is.logical(fill))
  mode = match.arg(mode)

  genome_build = match.arg(genome_build)

  requireNamespace("ggplot2")
  requireNamespace("cowplot")

  if (inherits(data, "CopyNumber")) {
    genome_build = data@genome_build
    data = data@annotation
  }

  # if remove normal copy number segments
  if (rm_normal) {
    if (!"segVal" %in% colnames(data)) {
      stop("'segVal' must be provided as a column.")
    }

    data = data[segVal != 2]
  }

  if (mode == "ld") {
    # plot length distribution
    if (!"fraction" %in% colnames(data)) {
      stop("'fraction' must be provided as a column.")
    }
    ggplot(data, aes(x = fraction, y = ..density..)) +
      geom_histogram(bins = 100) +
      labs(x = "Length of SCNA\n(normalized to chromosome arms)",
           y = "Percentage\n(as fraction of all SCNAs)")
  } else if (mode == "cd") {
    # plot chr distribution
    if (!all(c("chromosome", "location") %in% colnames(data))) {
      stop("'chromosome', 'location' must be provided as columns.")
    }

    if (is.character(data$chromosome[1])) {
      data$chromosome = sub("chr", "", data$chromosome, ignore.case = TRUE)
      data = subset(data, chromosome %in% c(1:22, "X"))
    } else if (is.integer(data$chromosome[1])) {
      data = subset(data, chromosome %in% c(1:22, "X"))
    }

    # only keep chr 1 to 22, X
    data$chromosome = factor(as.character(data$chromosome),
                             levels = c(1:22, "X"))

    # only keep p, q, pq
    data$location = factor(sub("[0-9X]*", "", data$location),
                           levels = c("p", "pq", "q"))
    # if (sum(!(data$location %in% c("p", "pq", "q"))) > 0){
    #   message("Discarding segments which located in centromere region...")
    #   data = subset(data, location %in% c("p", "pq", "q"))
    # }

    if (scale_chr) {
      if (genome_build == "hg19") {
        data("chromsize.hg19",
             package = "VSHunter",
             envir = environment())
        chrlen = chromsize.hg19
      } else {
        data("chromsize.hg38",
             package = "VSHunter",
             envir = environment())
        chrlen = chromsize.hg38
      }

      p = ggplot(data, aes(x = chromosome, fill = location)) +
        geom_bar() + xlab("Chromosome")


      q = ggplot_build(p)$data[[1]][, c("x", "count", "fill")]

      q$x = factor(q$x, levels = c(1:23), labels = c(1:22, "X"))
      q$fill = factor(q$fill, levels = c("#F8766D", "#00BA38", "#619CFF"))

      chrlen$chrom = gsub(
        pattern = "chr",
        replacement = "",
        x = chrlen$chrom
      )
      q = merge(q, chrlen, by.x = "x", by.y = "chrom")
      q[["count"]] = 1000000 * (q[["count"]] / q[["size"]])

      if (!fill) {
        ggplot(q, aes(x, y = count, fill = fill)) +
          geom_bar(stat = "identity") +
          scale_fill_discrete(name = "location",
                              labels = c("p", "pq", "q")) +
          labs(x = "Chromosome", y = "Normalized count (per Mb)")

      } else {
        ggplot(q, aes(x, y = count, fill = fill)) +
          geom_bar(stat = "identity", position = "fill") +
          scale_fill_discrete(name = "location",
                              labels = c("p", "pq", "q")) +
          labs(x = "Chromosome", y = "Percentage")
      }

    } else {
      # plot
      if (!fill) {
        ggplot(data, aes(x = chromosome, fill = location)) +
          geom_bar() + xlab("Chromosome")
      } else {
        ggplot(data, aes(x = chromosome, fill = location)) +
          geom_bar(position = "fill") + ylab("Percentage") + xlab("Chromosome")
      }
    }
  }
}



# Plot copy number features -----------------------------------------------



#' Plot copy number feature distribution
#'
#' @param features a `list` generate from [get_features] or [sig_prepare] function.
#' @param ylab lab of y axis.
#' @param ... other options pass to \code{\link[cowplot]{plot_grid}} function of `cowplot` package.
#'
#' @return a ggplot object
#' @import ggplot2
#' @family copy number plot
#' @export
#'
draw_cn_features = function(features, ylab = "", ...) {
  features = lapply(features, function(x) {
    x[["value"]] = as.numeric(x[["value"]])
    return(x)
  })

  requireNamespace("cowplot")
  #cowplot::theme_cowplot()
  #ggplot2::theme_set(cowplot::theme_cowplot(font_size = 12))
  p_1 = ggplot(data = features$segsize, aes(x = log10(value))) +
    geom_line(stat = "density") + labs(x = "Segment size (log10 based)", y = ylab) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + cowplot::theme_cowplot(font_size = 12)
  #p_1 = p_1 + scale_x_continuous(breaks = 10 ^c(7, 8),
  #                               labels = scales::trans_format("log10", scales::math_format(10^.x)))

  p_2 = ggplot(data = features$copynumber, aes(x = value)) +
    geom_line(stat = "density") + labs(x = "Copy number", y = ylab) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + cowplot::theme_cowplot(font_size = 12)
  p_3 = ggplot(data = features$changepoint, aes(x = value)) +
    geom_line(stat = "density") + labs(x = "Copy number changepoint", y = ylab) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + cowplot::theme_cowplot(font_size = 12)

  p_4 = ggplot(data = features$bp10MB, aes(x = value)) +
    geom_bar(stat = "count") + labs(x = "Breakpoint count per 10MB", y = ylab) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + cowplot::theme_cowplot(font_size = 12)
  p_5 = ggplot(data = features$bpchrarm, aes(x = value)) +
    geom_bar(stat = "count") + labs(x = "Breakpoint count per chr arm", y = ylab) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + cowplot::theme_cowplot(font_size = 12)

  p_6 = ggplot(data = features$osCN, aes(x = value)) +
    geom_bar(stat = "count") + labs(x = "Oscilating CN chain length", y = ylab) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + cowplot::theme_cowplot(font_size = 12)

  p = cowplot::plot_grid(p_1,
                         p_2,
                         p_3,
                         p_4,
                         p_5,
                         p_6,
                         nrow = 2,
                         align = "hv",
                         ...)
  p
}


# Plot fit components -----------------------------------------------------


#' Plot mixture fit model components
#'
#' @inheritParams draw_cn_features
#' @param components a `list` contain flexmix object of copy-number features, obtain this
#' from [get_components] function or use pre-compiled components data which come from CNV signature paper
#' <https://www.nature.com/articles/s41588-018-0179-8> (set this argument as `NULL`).
#' @param ... other options pass to \code{\link[cowplot]{plot_grid}} function of **cowplot** package.
#' @author Shixiang Wang <w_shixiang@163.com>
#' @return a ggplot object
#' @import ggplot2
#' @export
#' @family copy number plot
#'
draw_cn_components = function(features, components, ...) {

  requireNamespace("cowplot")
  cbPalette <-
    c(
      RColorBrewer::brewer.pal(8, "Dark2"),
      RColorBrewer::brewer.pal(9, "Set1"),
      "black"
    )
  plotNormDensity = function(value, matrix, xlab) {
    p = ggplot(data = data.frame(x = seq(
      min(value), max(value), length.out = 100
    )),
    aes(x)) + ylab("")

    for (i in 1:ncol(matrix)) {
      p = p + stat_function(
        fun = stats::dnorm,
        n = 1000,
        args = list(mean = matrix[1, i],
                    sd = matrix[2, i]),
        color = cbPalette[i]
      )
    }

    p = p + xlab(xlab) + cowplot::theme_cowplot(font_size = 12)
    p
  }

  plotPoisDensity = function(value, lambda, xlab, max_value = 10) {
    if (is.null(max_value)) {
      p = ggplot(data = data.frame(x = seq(
        min(value), max(value), length.out = 100
      )),
      aes(x)) + ylab("")
    } else {
      p = ggplot(data = data.frame(x = seq(
        min(value), max_value, length.out = 100
      )),
      aes(x)) + ylab("")
    }


    for (i in 1:length(lambda)) {
      p = p + stat_function(
        geom = "line",
        n = 11,
        fun = stats::dpois,
        args = list(lambda = lambda[i]),
        color = cbPalette[i]
      )
    }

    p = p + xlab(xlab) + cowplot::theme_cowplot(font_size = 12)
    p
  }



  features = lapply(features, function(x) {
    x$value = as.numeric(x$value)
    return(x)
  })

  # norm distribution
  comp_segsize = log10(flexmix::parameters(components[["segsize"]]))
  comp_copynumber = flexmix::parameters(components[["copynumber"]])
  comp_changepoint = flexmix::parameters(components[["changepoint"]])
  # pois distribution
  comp_bp10MB = flexmix::parameters(components[["bp10MB"]])
  comp_bpchrarm = flexmix::parameters(components[["bpchrarm"]])
  comp_osCN = flexmix::parameters(components[["osCN"]])

  # norm plots
  p_1 = plotNormDensity(log10(features[["segsize"]]$value), comp_segsize, xlab = "Segment size (log10 based)")
  # p_1 = p_1 + scale_x_continuous(breaks = 10 ^c(0, 7:9),
  #         labels = scales::trans_format("log10", scales::math_format(10^.x)))
  p_2 = plotNormDensity(features[["copynumber"]]$value, comp_copynumber, xlab = "Copy number")
  p_3 = plotNormDensity(features[["changepoint"]]$value, comp_changepoint, xlab = "Copy-number changepoint")

  # pois plots
  p_4 = plotPoisDensity(features[["bp10MB"]]$value,
                        comp_bp10MB,
                        xlab = "Breakpoint count per 10MB",
                        max_value = 10)
  p_5 = plotPoisDensity(features[["bpchrarm"]]$value,
                        comp_bpchrarm,
                        xlab = "Breakpoint count per arm",
                        max_value = 50)
  p_6 = plotPoisDensity(features[["osCN"]]$value, comp_osCN, xlab = "Oscilating CN chain length")

  p = cowplot::plot_grid(p_1,
                         p_2,
                         p_3,
                         p_4,
                         p_5,
                         p_6,
                         nrow = 2,
                         align = "hv",
                         ...)
  p
}





# Global variables --------------------------------------------------------

utils::globalVariables(
  c("segVal",
    "fraction",
    "..density..",
    "location",
    "x",
    "count",
    "value")
)
