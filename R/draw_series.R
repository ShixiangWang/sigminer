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

  ggplot2::theme_set(cowplot::theme_cowplot())
  #requireNamespace("ggplot2")
  #requireNamespace("cowplot")

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

  #requireNamespace("cowplot")
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



# Plot signature profile --------------------------------------------------

#' Plot signature profile
#'
#' Currently support copy number signatures and mutation signatures.
#'
#' @inheritParams sig_extract
#' @inheritParams sig_assign_samples
#' @param y_scale one of 'relative' or 'absolute', if choose 'relative',
#' signature columns will be scaled to sum as 1.
#' @param font_scale a number used to set font scale.
#' @author Shixiang Wang
#' @return a `ggplot` object
#' @import ggplot2
#' @export
#' @family signature plot
draw_sig_profile = function(nmfObj, mode = c("copynumber", "mutation"),
                            y_scale = c("relative", "absolute"), font_scale = 1) {

    mode = match.arg(mode)
    y_scale = match.arg(y_scale)

    #Signatures
    w = NMF::basis(nmfObj)
    if (y_scale == "relative") {
      w = apply(w, 2, function(x) x/sum(x)) #Scale the signatures (basis)
    }
    colnames(w) = paste('Signature', 1:ncol(w),sep='_')

    #>>>>>>>>>>>>>>>>> Setting theme
    scale <- font_scale

    .theme_ss <- theme_bw(base_size=12) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                       hjust = 1, size=10*scale, family = "mono"),
            axis.text.y = element_text(hjust = 0.5,size=12*scale, family = "mono"),
            axis.text = element_text(size = 14*scale, family = "mono"))
    #<<<<<<<<<<<<<<<<< Setting theme

    #>>>>>>>>>>>>>>>>> identify mode and do data transformation
    mat = as.data.frame(w)
    mat$context = rownames(mat)

    if (mode == "copynumber") {
      mat$base = sub("\\d+$", "",  mat$context)

      mat = tidyr::gather(mat, class, signature, dplyr::contains("Signature"))
      mat = dplyr::mutate(mat,
                          context = factor(context),
                          base = factor(base, levels = c("bp10MB", "copynumber",
                                                         "changepoint", "bpchrarm",
                                                         "osCN", "segsize")),
                          class = factor(class)
      )
    } else {
      mat$base = sub("[ACGT]\\[(.*)\\][ACGT]", "\\1", mat$context)
      mat$context = sub("(\\[.*\\])", "\\.", mat$context)

      mat = tidyr::gather(mat, class, signature, dplyr::contains("Signature"))
      mat = dplyr::mutate(mat,
                          context = factor(context),
                          base = factor(base, levels = c("C>A", "C>G",
                                                         "C>T", "T>A",
                                                         "T>C", "T>G")),
                          class = factor(class)
      )
    }

    #>>>>>>>>>>>>>>>>>>>>>>> Plot
    p = ggplot(mat) +
      geom_bar(aes_string(x="context",y="signature",fill="base"),
               stat="identity", position="identity", colour="gray50") +
      scale_fill_manual(values=c("cyan","red","yellow","purple",
                                 "green","blue","black","gray"))

    if (mode == "copynumber") {
      p = p + facet_grid(class ~ ., scales = "free")
    } else {
      p = p + facet_grid(class ~ base, scales = "free")
    }

    p = p +
      guides(fill=FALSE) + .theme_ss +
      theme(axis.title.x = element_text(face="bold",colour="black",size=14*scale)) +
      theme(axis.title.y = element_text(face="bold",colour="black",size=14*scale))

    if (mode == "copynumber") {
      p = p + xlab("Components") + ylab("Contributions")
    } else {
      p = p + xlab("Motifs") + ylab("Contributions")
    }

    return(p)
}


# Plot signature activity -------------------------------------------------

#' Plot signature activity
#'
#' Currently support copy number signatures and mutation signatures.
#' @inheritParams draw_sig_profile
#' @param hide_samps if `TRUE`, not show sample names.
#' @author Shixiang Wang
#' @return a `ggplot` object
#' @import ggplot2
#' @importFrom grDevices rainbow
#' @export
#' @family signature plot
draw_sig_activity = function(nmfObj, mode = c("copynumber", "mutation"),
                             font_scale = 1, hide_samps = TRUE) {
  mode = match.arg(mode)

  scale = font_scale
  .theme_ss <- theme_bw(base_size=12) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                     hjust = 1, size=10*scale, family = "mono"),
          axis.text.y = element_text(hjust = 0.5,size=12*scale, family = "mono"),
          axis.text = element_text(size = 8*scale, family = "mono"))

  activity = sig_get_activity(nmfObj)

  h = activity[["absolute"]]
  h.norm = activity[["relative"]]

  # chop Signature off
  rownames(h) = sub(".*_(\\d+)$", "\\1", rownames(h))
  rownames(h.norm) = sub(".*_(\\d+)$", "\\1", rownames(h.norm))

  ordering = order(colSums(h),decreasing=TRUE)
  h = h[, ordering]
  h.norm = h.norm[, ordering]

  sample.ordering = colnames(h)

  h$Signature = rownames(h)
  h.norm$Signature = rownames(h.norm)
  x1 = tidyr::gather(h, Sample, Activity, -Signature)
  x2 = tidyr::gather(h.norm, Sample, Activity, -Signature)
  if (mode == "copynumber") {
    x1$class0 = "Weights"
  } else {
    x1$class0 = "Counts"
  }
  x2$class0 = "Fractions"
  df = rbind(x1, x2)

  if (mode == "copynumber") {
    df$class0 <- factor(df$class0,c("Weights","Fractions"))
  } else {
    df$class0 <- factor(df$class0,c("Counts","Fractions"))
  }

  df$Sample <- factor(df$Sample,sample.ordering)

  p = ggplot(df,aes(x=Sample,y=Activity,fill=Signature))
  p = p + geom_bar(stat="identity",position='stack',color='black',alpha=0.9)
  p = p + scale_fill_manual(values=c("red","cyan","yellow","blue","magenta","gray50","orange","darkgreen","brown","black",rainbow(10)[4:10]))
  p = p + facet_grid(class0 ~ ., scales = "free_y")
  p = p + xlab("Samples") + ylab("Signature Activities")
  p = p + .theme_ss
  p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=14*scale))
  p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=14*scale))
  #p = p + theme(legend.title=element_blank())
  p = p + theme(legend.position="top")

  if (hide_samps) {
    p = p + theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
  }

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
    "value",
    "Activity",
    "Sample",
    "Signature",
    "base",
    "context")
)
