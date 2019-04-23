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
#' @examples
#' # Load copy number object
#' load(system.file("extdata", "toy_copynumber.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Plot distribution
#' draw_cn_distribution(cn)
#' draw_cn_distribution(cn, mode = "cd")
#' draw_cn_distribution(cn, mode = "cd", fill = TRUE)
#' @family copy number plot
draw_cn_distribution <- function(data,
                                 rm_normal = TRUE,
                                 mode = c("ld", "cd"),
                                 fill = FALSE,
                                 scale_chr = TRUE,
                                 genome_build = c("hg19", "hg38")) {
  stopifnot(
    is.logical(rm_normal),
    is.data.frame(data) | inherits(data, "CopyNumber"),
    is.logical(fill)
  )
  mode <- match.arg(mode)

  genome_build <- match.arg(genome_build)

  ggplot2::theme_set(cowplot::theme_cowplot())
  # requireNamespace("ggplot2")
  # requireNamespace("cowplot")

  if (inherits(data, "CopyNumber")) {
    genome_build <- data@genome_build
    data <- data@annotation
  }

  # if remove normal copy number segments
  if (rm_normal) {
    if (!"segVal" %in% colnames(data)) {
      stop("'segVal' must be provided as a column.")
    }

    data <- data[segVal != 2]
  }

  if (mode == "ld") {
    # plot length distribution
    if (!"fraction" %in% colnames(data)) {
      stop("'fraction' must be provided as a column.")
    }
    ggplot(data, aes(x = fraction, y = ..density..)) +
      geom_histogram(bins = 100) +
      labs(
        x = "Length of SCNA\n(normalized to chromosome arms)",
        y = "Percentage\n(as fraction of all SCNAs)"
      )
  } else if (mode == "cd") {
    # plot chr distribution
    if (!all(c("chromosome", "location") %in% colnames(data))) {
      stop("'chromosome', 'location' must be provided as columns.")
    }

    if (is.character(data$chromosome[1])) {
      data$chromosome <- sub("chr", "", data$chromosome, ignore.case = TRUE)
      data <- subset(data, chromosome %in% c(1:22, "X"))
    } else if (is.integer(data$chromosome[1])) {
      data <- subset(data, chromosome %in% c(1:22, "X"))
    }

    # only keep chr 1 to 22, X
    data$chromosome <- factor(as.character(data$chromosome),
      levels = c(1:22, "X")
    )

    # only keep p, q, pq
    data$location <- factor(sub("[0-9X]*", "", data$location),
      levels = c("p", "pq", "q")
    )
    # if (sum(!(data$location %in% c("p", "pq", "q"))) > 0){
    #   message("Discarding segments which located in centromere region...")
    #   data = subset(data, location %in% c("p", "pq", "q"))
    # }

    if (scale_chr) {
      if (genome_build == "hg19") {
        data("chromsize.hg19",
          package = "sigminer",
          envir = environment()
        )
        chrlen <- chromsize.hg19
      } else {
        data("chromsize.hg38",
          package = "sigminer",
          envir = environment()
        )
        chrlen <- chromsize.hg38
      }

      p <- ggplot(data, aes(x = chromosome, fill = location)) +
        geom_bar() + xlab("Chromosome")


      q <- ggplot_build(p)$data[[1]][, c("x", "count", "fill")]

      q$x <- factor(q$x, levels = c(1:23), labels = c(1:22, "X"))
      q$fill <- factor(q$fill, levels = c("#F8766D", "#00BA38", "#619CFF"))

      chrlen$chrom <- gsub(
        pattern = "chr",
        replacement = "",
        x = chrlen$chrom
      )
      q <- merge(q, chrlen, by.x = "x", by.y = "chrom")
      q[["count"]] <- 1000000 * (q[["count"]] / q[["size"]])

      if (!fill) {
        ggplot(q, aes(x, y = count, fill = fill)) +
          geom_bar(stat = "identity") +
          scale_fill_discrete(
            name = "location",
            labels = c("p", "pq", "q")
          ) +
          labs(x = "Chromosome", y = "Normalized count (per Mb)")
      } else {
        ggplot(q, aes(x, y = count, fill = fill)) +
          geom_bar(stat = "identity", position = "fill") +
          scale_fill_discrete(
            name = "location",
            labels = c("p", "pq", "q")
          ) +
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
#' @examples
#' # Load copy number prepare object
#' load(system.file("extdata", "toy_copynumber_prepare.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' draw_cn_features(cn_prepare$features)
#' @export
#'
draw_cn_features <- function(features, ylab = "", ...) {
  features <- lapply(features, function(x) {
    x[["value"]] <- as.numeric(x[["value"]])
    return(x)
  })

  # requireNamespace("cowplot")
  # cowplot::theme_cowplot()
  # ggplot2::theme_set(cowplot::theme_cowplot(font_size = 12))
  p_1 <- ggplot(data = features$segsize, aes(x = log10(value))) +
    geom_line(stat = "density") + labs(x = "Segment size (log10 based)", y = ylab) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + cowplot::theme_cowplot(font_size = 12)
  # p_1 = p_1 + scale_x_continuous(breaks = 10 ^c(7, 8),
  #                               labels = scales::trans_format("log10", scales::math_format(10^.x)))

  p_2 <- ggplot(data = features$copynumber, aes(x = value)) +
    geom_line(stat = "density") + labs(x = "Copy number", y = ylab) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + cowplot::theme_cowplot(font_size = 12)

  p_3 <- ggplot(data = features$changepoint, aes(x = value)) +
    geom_line(stat = "density") + labs(x = "Copy number changepoint", y = ylab) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + cowplot::theme_cowplot(font_size = 12)

  p_4 <- ggplot(data = features$bp10MB, aes(x = value)) +
    geom_bar(stat = "count") + labs(x = "Breakpoint count per 10MB", y = ylab) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + cowplot::theme_cowplot(font_size = 12)
  p_5 <- ggplot(data = features$bpchrarm, aes(x = value)) +
    geom_bar(stat = "count") + labs(x = "Breakpoint count per chr arm", y = ylab) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + cowplot::theme_cowplot(font_size = 12)

  p_6 <- ggplot(data = features$osCN, aes(x = value)) +
    geom_bar(stat = "count") + labs(x = "Oscilating CN chain length", y = ylab) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + cowplot::theme_cowplot(font_size = 12)

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
#' @examples
#' # Load copy number prepare object
#' load(system.file("extdata", "toy_copynumber_prepare.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' draw_cn_components(cn_prepare$features, cn_prepare$components)
#' @family copy number plot
#'
draw_cn_components <- function(features, components, ...) {
  requireNamespace("cowplot")
  cbPalette <-
    c(
      RColorBrewer::brewer.pal(8, "Dark2"),
      RColorBrewer::brewer.pal(9, "Set1"),
      "black"
    )
  plotNormDensity <- function(value, matrix, xlab) {
    p <- ggplot(
      data = data.frame(x = seq(
        min(value), max(value),
        length.out = 100
      )),
      aes(x)
    ) + ylab("")

    for (i in 1:ncol(matrix)) {
      p <- p + stat_function(
        fun = stats::dnorm,
        n = 1000,
        args = list(
          mean = matrix[1, i],
          sd = matrix[2, i]
        ),
        color = cbPalette[i]
      )
    }

    p <- p + xlab(xlab) + cowplot::theme_cowplot(font_size = 12)
    p
  }

  plotPoisDensity <- function(value, lambda, xlab, max_value = 10) {
    if (is.null(max_value)) {
      p <- ggplot(
        data = data.frame(x = seq(
          min(value), max(value),
          length.out = 100
        )),
        aes(x)
      ) + ylab("")
    } else {
      p <- ggplot(
        data = data.frame(x = seq(
          min(value), max_value,
          length.out = 100
        )),
        aes(x)
      ) + ylab("")
    }


    for (i in 1:length(lambda)) {
      p <- p + stat_function(
        geom = "line",
        n = 11,
        fun = stats::dpois,
        args = list(lambda = lambda[i]),
        color = cbPalette[i]
      )
    }

    p <- p + xlab(xlab) + cowplot::theme_cowplot(font_size = 12)
    p
  }



  features <- lapply(features, function(x) {
    x$value <- as.numeric(x$value)
    return(x)
  })

  # norm distribution
  comp_segsize <- log10(flexmix::parameters(components[["segsize"]]))
  comp_copynumber <- flexmix::parameters(components[["copynumber"]])
  comp_changepoint <- flexmix::parameters(components[["changepoint"]])
  # pois distribution
  comp_bp10MB <- flexmix::parameters(components[["bp10MB"]])
  comp_bpchrarm <- flexmix::parameters(components[["bpchrarm"]])
  comp_osCN <- flexmix::parameters(components[["osCN"]])

  # norm plots
  p_1 <- plotNormDensity(log10(features[["segsize"]]$value), comp_segsize, xlab = "Segment size (log10 based)")
  # p_1 = p_1 + scale_x_continuous(breaks = 10 ^c(0, 7:9),
  #         labels = scales::trans_format("log10", scales::math_format(10^.x)))
  p_2 <- plotNormDensity(features[["copynumber"]]$value, comp_copynumber, xlab = "Copy number")
  p_3 <- plotNormDensity(features[["changepoint"]]$value, comp_changepoint, xlab = "Copy-number changepoint")

  # pois plots
  p_4 <- plotPoisDensity(features[["bp10MB"]]$value,
    comp_bp10MB,
    xlab = "Breakpoint count per 10MB",
    max_value = 10
  )
  p_5 <- plotPoisDensity(features[["bpchrarm"]]$value,
    comp_bpchrarm,
    xlab = "Breakpoint count per arm",
    max_value = 50
  )
  p_6 <- plotPoisDensity(features[["osCN"]]$value, comp_osCN, xlab = "Oscilating CN chain length")

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
#' @examples
#' # Load copy number signature
#' load(system.file("extdata", "toy_copynumber_signature.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' draw_sig_profile(res$nmfObj)
#' @family signature plot
draw_sig_profile <- function(nmfObj, mode = c("copynumber", "mutation"),
                             y_scale = c("relative", "absolute"), font_scale = 1) {
  mode <- match.arg(mode)
  y_scale <- match.arg(y_scale)

  # Signatures
  w <- NMF::basis(nmfObj)
  if (y_scale == "relative") {
    w <- apply(w, 2, function(x) x / sum(x)) # Scale the signatures (basis)
  }
  colnames(w) <- paste("Signature", 1:ncol(w), sep = "_")

  # >>>>>>>>>>>>>>>>> Setting theme
  scale <- font_scale

  .theme_ss <- theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(
        angle = 90, vjust = 0.5,
        hjust = 1, size = 10 * scale, family = "mono"
      ),
      axis.text.y = element_text(hjust = 0.5, size = 12 * scale, family = "mono"),
      axis.text = element_text(size = 14 * scale, family = "mono")
    )
  # <<<<<<<<<<<<<<<<< Setting theme

  # >>>>>>>>>>>>>>>>> identify mode and do data transformation
  mat <- as.data.frame(w)
  mat$context <- rownames(mat)

  if (mode == "copynumber") {
    mat$base <- sub("\\d+$", "", mat$context)

    mat <- tidyr::gather(mat, class, signature, dplyr::contains("Signature"))
    mat <- dplyr::mutate(mat,
      context = factor(context),
      base = factor(base, levels = c(
        "bp10MB", "copynumber",
        "changepoint", "bpchrarm",
        "osCN", "segsize"
      )),
      class = factor(class)
    )
  } else {
    mat$base <- sub("[ACGT]\\[(.*)\\][ACGT]", "\\1", mat$context)
    mat$context <- sub("(\\[.*\\])", "\\.", mat$context)

    mat <- tidyr::gather(mat, class, signature, dplyr::contains("Signature"))
    mat <- dplyr::mutate(mat,
      context = factor(context),
      base = factor(base, levels = c(
        "C>A", "C>G",
        "C>T", "T>A",
        "T>C", "T>G"
      )),
      class = factor(class)
    )
  }

  # >>>>>>>>>>>>>>>>>>>>>>> Plot
  p <- ggplot(mat) +
    geom_bar(aes_string(x = "context", y = "signature", fill = "base"),
      stat = "identity", position = "identity", colour = "gray50"
    ) +
    scale_fill_manual(values = c(
      "cyan", "red", "yellow", "purple",
      "green", "blue", "black", "gray"
    ))

  if (mode == "copynumber") {
    p <- p + facet_grid(class ~ ., scales = "free")
  } else {
    p <- p + facet_grid(class ~ base, scales = "free")
  }

  p <- p +
    guides(fill = FALSE) + .theme_ss +
    theme(axis.title.x = element_text(face = "bold", colour = "black", size = 14 * scale)) +
    theme(axis.title.y = element_text(face = "bold", colour = "black", size = 14 * scale))

  p <- p + xlab("Components") + ylab("Contributions")
  return(p)
}


# Plot signature activity -------------------------------------------------

#' Plot signature activity
#'
#' Currently support copy number signatures and mutational signatures.
#' @inheritParams draw_sig_profile
#' @param hide_samps if `TRUE`, not show sample names.
#' @author Shixiang Wang
#' @return a `ggplot` object
#' @import ggplot2
#' @importFrom grDevices rainbow
#' @export
#' @examples
#' # Load copy number signature
#' load(system.file("extdata", "toy_copynumber_signature.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' draw_sig_activity(res$nmfObj)
#' @family signature plot
draw_sig_activity <- function(nmfObj, mode = c("copynumber", "mutation"),
                              font_scale = 1, hide_samps = TRUE) {
  mode <- match.arg(mode)

  scale <- font_scale
  .theme_ss <- theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(
        angle = 90, vjust = 0.5,
        hjust = 1, size = 10 * scale, family = "mono"
      ),
      axis.text.y = element_text(hjust = 0.5, size = 12 * scale, family = "mono"),
      axis.text = element_text(size = 8 * scale, family = "mono")
    )

  activity <- sig_get_activity(nmfObj)

  h <- activity[["absolute"]]
  h.norm <- activity[["relative"]]

  # chop Signature off
  rownames(h) <- sub(".*_(\\d+)$", "\\1", rownames(h))
  rownames(h.norm) <- sub(".*_(\\d+)$", "\\1", rownames(h.norm))

  ordering <- order(colSums(h), decreasing = TRUE)
  h <- h[, ordering]
  h.norm <- h.norm[, ordering]

  sample.ordering <- colnames(h)

  h$Signature <- rownames(h)
  h.norm$Signature <- rownames(h.norm)
  x1 <- tidyr::gather(h, Sample, Activity, -Signature)
  x2 <- tidyr::gather(h.norm, Sample, Activity, -Signature)
  # if (mode == "copynumber") {
  #   x1$class0 = "Weights"
  # } else {
  #   x1$class0 = "Counts"
  # }
  x1$class0 <- "Coefficients"
  x2$class0 <- "Fractions"
  df <- rbind(x1, x2)

  if (mode == "copynumber") {
    df$class0 <- factor(df$class0, c("Coefficients", "Fractions"))
  } else {
    df$class0 <- factor(df$class0, c("Coefficients", "Fractions"))
  }

  df$Sample <- factor(df$Sample, sample.ordering)

  p <- ggplot(df, aes(x = Sample, y = Activity, fill = Signature))
  p <- p + geom_bar(stat = "identity", position = "stack", color = "black", alpha = 0.9)
  p <- p + scale_fill_manual(values = c("red", "cyan", "yellow", "blue", "magenta", "gray50", "orange", "darkgreen", "brown", "black", rainbow(10)[4:10]))
  p <- p + facet_grid(class0 ~ ., scales = "free_y")
  p <- p + xlab("Samples") + ylab("Signature Activities")
  p <- p + .theme_ss
  p <- p + theme(axis.title.x = element_text(face = "bold", colour = "black", size = 14 * scale))
  p <- p + theme(axis.title.y = element_text(face = "bold", colour = "black", size = 14 * scale))
  # p = p + theme(legend.title=element_blank())
  p <- p + theme(legend.position = "top")

  if (hide_samps) {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  }

  p
}


# Plot correlation between signature activity -----------------------------

#' Plot correlation between signature activities
#'
#' @inheritParams corrplot::corrplot
#' @param mat_list a `list` contain correlation and p value matrix etc., obtain
#' it from [sig_get_correlation] function.
#' @param ... other arguments pass to [corrplot::corrplot()] function.
#' @author Shixiang Wang
#' @return NULL
#' @importFrom grDevices colorRampPalette
#' @export
#' @examples
#' \donttest{
#' # Load copy number signature
#' load(system.file("extdata", "toy_copynumber_signature.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' sig_activity <- sig_get_activity(res$nmfObj)
#' sig_cor <- sig_get_correlation(sig_activity)
#' draw_sig_corrplot(sig_cor)
#' }
#' @family signature plot
draw_sig_corrplot <- function(mat_list, order = "original", type = "lower",
                              sig.level = .05, ...) {
  col3 <- colorRampPalette(c("blue", "white", "red"))
  corrplot::corrplot(mat_list[["correlation"]],
    order = order, col = col3(20),
    type = type, tl.col = "black", tl.srt = 45,
    p.mat = mat_list[["p"]], sig.level = sig.level,
    ...
  )
}



# Plot comparison between signature subtypes ------------------------------

#' Plot comparison between signature subtypes
#'
#' Using result data from [sig_summarize_subtypes] function, this function plot
#' genotypes/phenotypes comparison between signature subtypes using **ggplot2** package and return
#' a list of `ggplot` object contains individual and combined plots. The combined
#' plot is easily saved to local using [cowplot::save_plot()].
#' @param subtype_summary a `list` from result of [sig_summarize_subtypes] function.
#' @param xlab lab name of x axis for all plots.
#' @param ylab_co lab name of y axis for plots of continuous type data. Of note,
#' this argument should be a character vector has same length as `subtype_summary`,
#' the location for categorical type data should mark with `NA`.
#' @param legend_title_ca legend title for plots of categorical type data.
#' Of note,
#' this argument should be a character vector has same length as `subtype_summary`,
#' the location for continuous type data should mark with `NA`.
#' @param show_pvalue if `TRUE`, show p value for comparison of continuous data types.
#' @param ... other paramters pass to [ggpubr::stat_compare_means()].
#' @author Shixiang Wang <w_shixiang@163.com>
#' @return a `list` of `ggplot` objects.
#' @import ggplot2
#' @export
#' @examples
#' # Load copy number signature
#' load(system.file("extdata", "toy_copynumber_signature.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Assign samples to clusters
#' subtypes <- sig_assign_samples(res$nmfObj, type = "samples")
#' 
#' set.seed(1234)
#' # Add custom groups
#' subtypes$new_group <- sample(c("1", "2", "3", "4"), size = nrow(subtypes), replace = TRUE)
#' # Summarize subtypes
#' subtypes.sum <- sig_summarize_subtypes(subtypes[, -1],
#'   col_subtype = "nmf_subtypes",
#'   cols_to_summary = colnames(subtypes[, -1])[c(-1, -2)],
#'   type = c("co", "ca"), verbose = TRUE
#' )
#' draw_subtypes_comparison(subtypes.sum)
#' @family signature plot
draw_subtypes_comparison <- function(subtype_summary,
                                     xlab = "subtype", ylab_co = NA,
                                     legend_title_ca = NA,
                                     show_pvalue = TRUE,
                                     ...) {
  # parameter with ca/co in the end need fill values

  # add parameters
  subtype_summary <- Map(function(x, xlab, legend_title_ca, ylab_co) {
    x[["xlab"]] <- xlab
    if (x[["type"]] == "categorical") x[["legend_title"]] <- legend_title_ca
    if (x[["type"]] == "continuous") x[["ylab"]] <- ylab_co
    x
  }, subtype_summary, xlab, legend_title_ca, ylab_co)

  # split input into two list according to element is categorical or continuous
  ca_index <- which(sapply(subtype_summary, function(x) x$type) == "categorical")
  co_index <- which(sapply(subtype_summary, function(x) x$type) == "continuous")

  n_left <- length(subtype_summary) - length(ca_index) - length(co_index)

  if (n_left > 0) {
    warning(n_left, " elements drop of because their type is not either categorical or continuous.")
  }

  if (length(ca_index) > 0) ca_list <- subtype_summary[ca_index]
  if (length(co_index) > 0) co_list <- subtype_summary[co_index]

  # library(ggplot2)
  # library(cowplot)
  if (length(ca_index) > 0) {
    # plot categorical data
    ca_res <- lapply(ca_list, function(df) {
      data <- df[["data"]]
      p <- ggplot(data, aes_string(x = "subtype", fill = colnames(data)[2])) +
        geom_bar(position = "fill") +
        xlab(df[["xlab"]]) +
        cowplot::theme_cowplot() +
        theme(axis.title.y = element_blank())

      if (!is.na(df[["legend_title"]])) {
        p <- p + scale_fill_discrete(name = df[["legend_title"]])
      }
      p
    })
    names(ca_res) <- names(ca_list)
  } else {
    ca_res <- NA
  }

  if (length(co_index) > 0) {
    # plot continuous data
    co_res <- lapply(co_list, function(df, ...) {
      data <- df[["data"]]
      my_comparisons <- combn(unique(as.character(data[["subtype"]])),
        2,
        simplify = FALSE
      )
      p <- ggplot(data, aes_string(x = "subtype", y = colnames(data)[2])) +
        geom_boxplot() + xlab(df[["xlab"]]) + cowplot::theme_cowplot()

      if (show_pvalue) {
        if (!requireNamespace("ggpubr")) {
          stop("'ggpubr' package is need for plotting p values.")
        }
        p <- p + ggpubr::stat_compare_means(
          comparisons = my_comparisons,
          ...
        )
      }
      p
    }, ...)
    names(co_res) <- names(co_list)
  } else {
    co_res <- NA
  }

  if (all(!is.na(ca_res))) {
    if (length(ca_res) <= 3) {
      ca_comb <- cowplot::plot_grid(plotlist = ca_res, align = "h")
    } else if (length(ca_res) == 4) {
      ca_comb <- cowplot::plot_grid(plotlist = ca_res, align = "h", ncol = 2)
    } else if (length(ca_res) <= 9) {
      ca_comb <- cowplot::plot_grid(plotlist = ca_res, align = "h", ncol = 3)
    } else {
      ca_comb <- cowplot::plot_grid(plotlist = ca_res, align = "h", ncol = 4)
    }
  } else {
    ca_comb <- NA
  }

  if (all(!is.na(co_res))) {
    if (length(co_res) <= 3) {
      co_comb <- cowplot::plot_grid(plotlist = co_res, align = "h")
    } else if (length(co_res) == 4) {
      co_comb <- cowplot::plot_grid(
        plotlist = co_res, align = "h",
        axis = "l", ncol = 2
      )
    } else if (length(co_res) <= 9) {
      co_comb <- cowplot::plot_grid(
        plotlist = co_res, align = "h",
        axis = "l", ncol = 3
      )
    } else {
      co_comb <- cowplot::plot_grid(
        plotlist = co_res, align = "h",
        axis = "l", ncol = 4
      )
    }
  } else {
    co_comb <- NA
  }

  list(
    ca = ca_res,
    ca_comb = ca_comb,
    co = co_res,
    co_comb = co_comb
  )
}

# Global variables --------------------------------------------------------

utils::globalVariables(
  c(
    "segVal",
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
    "context"
  )
)
