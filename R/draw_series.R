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

  # Sort params
  comp_segsize <- comp_segsize[, order(comp_segsize[1, ])]
  comp_copynumber <- comp_copynumber[, order(comp_copynumber[1, ])]
  comp_changepoint <- comp_changepoint[, order(comp_changepoint[1, ])]

  comp_bp10MB <- comp_bp10MB[order(comp_bp10MB)]
  comp_bpchrarm <- comp_bpchrarm[order(comp_bpchrarm)]
  comp_osCN <- comp_osCN[order(comp_osCN)]

  # output parameters
  parameters <- dplyr::tibble(
    components = c(
      paste0("segsize", 1:ncol(comp_segsize)),
      paste0("copynumber", 1:ncol(comp_copynumber)),
      paste0("changepoint", 1:ncol(comp_changepoint)),
      paste0("bp10MB", 1:length(comp_bp10MB)),
      paste0("bpchrarm", 1:length(comp_bpchrarm)),
      paste0("osCN", 1:length(comp_osCN))
    ),
    stats = c(
      comp_segsize[1, ], comp_copynumber[1, ],
      comp_changepoint[1, ],
      comp_bp10MB, comp_bpchrarm,
      comp_osCN
    ),
    dist = c(
      rep("norm", ncol(comp_segsize) + ncol(comp_copynumber) + ncol(comp_changepoint)),
      rep("pois", length(comp_bp10MB) + length(comp_bpchrarm) + length(comp_osCN))
    )
  )

  # norm plots
  p_1 <- plotNormDensity(log10(features[["segsize"]]$value), comp_segsize, xlab = "Segment size (log10 based)")
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
  p[["parameters"]] <- parameters
  p
}



# Plot signature profile --------------------------------------------------

#' Plot signature profile
#'
#' Currently support copy number signatures and mutation signatures.
#'
#' @inheritParams sig_extract
#' @inheritParams sig_assign_samples
#' @param params params data of components.
#' @param y_expand y expand height for plotting params of copy number signatures.
#' @param digits digits for plotting params of copy number signatures.
#' @param y_scale one of 'relative' or 'absolute', if choose 'relative',
#' signature columns will be scaled to sum as 1.
#' @param font_scale a number used to set font scale.
#' @param sig_names set name of signatures, can be a character vector.
#' Default is `NULL`, prefix 'Signature_' plus number is used.
#' @param sig_orders set order of signatures, can be a character vector.
#' Default is `NULL`, the signatures are ordered by alphabetical order.
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
draw_sig_profile <- function(nmfObj, mode = c("copynumber", "mutation"), params=NULL, y_expand=1,
                             digits = 1,
                             y_scale = c("relative", "absolute"), font_scale = 1,
                             sig_names = NULL, sig_orders = NULL) {
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
  # w = rbind(w, matrix(c(0.1, 0.2), ncol = 2, dimnames = list("segsize10")))
  mat <- as.data.frame(w)
  mat$context <- rownames(mat)

  if (mode == "copynumber") {
    mat$base <- sub("\\d+$", "", mat$context)

    mat <- tidyr::gather(mat, class, signature, dplyr::contains("Signature"))
    mat <- dplyr::mutate(mat,
      context = factor(context,
        levels = unique(mat[["context"]])
      ),
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

  # >>>>>>> Set signature name and order
  if (!is.null(sig_names)) {
    if (length(sig_names) != length(unique(mat[["class"]]))) {
      stop("The length of input signature names is not equal to signature number")
    }
    names(sig_names) <- paste0("Signature_", seq_along(sig_names))
    mat[["class"]] <- sig_names[mat[["class"]]]
  }

  if (!is.null(sig_orders)) {
    mat[["class"]] <- factor(mat[["class"]], levels = sig_orders)
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

    if (!is.null(params)) {
      params$class = levels(mat[["class"]])[1]
      p <- p + geom_text(aes(x=components, y=Inf,
                             label=ifelse(dist=="norm",
                                          paste0("\u03BC=",round(stats,digits)),
                                          paste0("\u03BB=",round(stats,digits)))),
                         data = params, angle=90,
                         hjust=-0.1, vjust=0.5) +
        coord_cartesian(clip = "off")
    }
    p <- p + facet_grid(class ~ ., scales = "free")

  } else {
    p <- p + facet_grid(class ~ base, scales = "free")
  }

  p <- p +
    guides(fill = FALSE) + .theme_ss +
    theme(axis.title.x = element_text(face = "bold", colour = "black", size = 14 * scale)) +
    theme(axis.title.y = element_text(face = "bold", colour = "black", size = 14 * scale))

  if (all(mode=="copynumber", !is.null(params))) {
    p <- p + theme(plot.margin = margin(30*y_expand, 2, 2, 2, unit = "pt")) # Add regions
  }

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
#' plot is easily saved to local using [cowplot::save_plot()]. Of note, default fisher
#' test p values are showed for categorical data and fdr values are showd for
#' continuous data.
#' @inheritParams get_adj_p
#' @param subtype_summary a `list` from result of [sig_summarize_subtypes] function.
#' @param xlab lab name of x axis for all plots. if it is `NA`, remove title for x axis.
#' @param ylab_co lab name of y axis for plots of continuous type data. Of note,
#' this argument should be a character vector has same length as `subtype_summary`,
#' the location for categorical type data should mark with `NA`.
#' @param legend_title_ca legend title for plots of categorical type data.
#' @param legend_position_ca legend position for plots of categorical type data.
#' Of note,
#' this argument should be a character vector has same length as `subtype_summary`,
#' the location for continuous type data should mark with `NA`.
#' @param show_pvalue if `TRUE`, show p values.
#' @param ... other paramters pass to [ggpubr::compare_means()] or [ggpubr::stat_compare_means()]
#' according to the specified `method`.
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
#'   col_subtype = "Subgroup",
#'   cols_to_summary = colnames(subtypes[, -1])[c(-1, -2)],
#'   type = c("co", "ca"), verbose = TRUE
#' )
#' draw_subtypes_comparison(subtypes.sum)
#'
#' # Another example
#' library(ggpubr)
#' tt <- sig_summarize_subtypes(ToothGrowth,
#'   col_subtype = "dose",
#'   cols_to_summary = "len",
#'   type = "co", verbose = TRUE
#' )
#' tt_plot <- draw_subtypes_comparison(tt)
#' tt_plot$co$len
#' @family signature plot
draw_subtypes_comparison <- function(subtype_summary,
                                     xlab = "subtype", ylab_co = NA,
                                     legend_title_ca = NA,
                                     legend_position_ca = "bottom",
                                     show_pvalue = TRUE,
                                     method = "wilcox.test", p.adjust.method = "fdr",
                                     ...) {
  # parameter with ca/co in the end need fill values

  # add parameters
  subtype_summary <- Map(function(x, xlab, legend_title_ca, ylab_co, legend_position_ca) {
    x[["xlab"]] <- xlab
    if (x[["type"]] == "categorical") {
      x[["legend_title"]] <- legend_title_ca
      x[["legend_position"]] <- legend_position_ca
    }
    if (x[["type"]] == "continuous") x[["ylab"]] <- ylab_co
    x
  }, subtype_summary, xlab, legend_title_ca, ylab_co, legend_position_ca)

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
      data_sum <- data %>% dplyr::count_("subtype")
      data_sum[["labels"]] <- paste(data_sum[["subtype"]], paste0("(n=", data_sum[["n"]], ")"), sep = "\n")

      var_name <- colnames(data)[2]
      var_name2 <- ifelse(isValidAndUnreserved(var_name), var_name, paste0("`", var_name, "`"))

      p <- ggplot(data, aes_string(x = "subtype", fill = var_name2)) +
        geom_bar(position = "fill") +
        cowplot::theme_cowplot() +
        theme(axis.title.y = element_blank()) +
        scale_x_discrete(
          breaks = data_sum[["subtype"]],
          labels = data_sum[["labels"]]
        )
      if (is.na(df[["xlab"]])) {
        p <- p + theme(axis.title.x = element_blank())
      } else {
        p <- p + xlab(df[["xlab"]])
      }

      if (!is.na(df[["legend_title"]])) {
        p <- p + scale_fill_discrete(name = df[["legend_title"]])
      }

      if (show_pvalue) {
        if (!is.na(df[["p_value"]])) {
          p <- p + labs(title = paste0("P=", signif(df[["p_value"]], 3)))
        }
      }

      p + theme(legend.position = df[["legend_position"]])
    })
    names(ca_res) <- names(ca_list)
  } else {
    ca_res <- NA
  }

  if (length(co_index) > 0) {
    # plot continuous data
    co_res <- lapply(co_list, function(df, ...) {
      data <- df[["data"]]
      data_sum <- data %>% dplyr::count_("subtype")
      data_sum[["labels"]] <- paste(data_sum[["subtype"]], paste0("(n=", data_sum[["n"]], ")"), sep = "\n")
      my_comparisons <- combn(unique(as.character(data[["subtype"]])),
        2,
        simplify = FALSE
      )

      var_name <- colnames(data)[2]
      var_name2 <- ifelse(isValidAndUnreserved(var_name), var_name, paste0("`", var_name, "`"))

      p <- ggplot(data, aes_string(x = "subtype", y = var_name2)) +
        geom_boxplot() + cowplot::theme_cowplot() +
        scale_x_discrete(
          breaks = data_sum[["subtype"]],
          labels = data_sum[["labels"]]
        )

      if (is.na(df[["xlab"]])) {
        p <- p + theme(axis.title.x = element_blank())
      } else {
        p <- p + xlab(df[["xlab"]])
      }

      if (!is.na(df[["ylab"]])) {
        p <- p + ylab(df[["ylab"]])
      }

      if (show_pvalue) {
        if (!requireNamespace("ggpubr")) {
          stop("'ggpubr' package is needed for plotting p values.")
        }

        if ((method != "wilcox.test") & (method != "t.test")) {
          p <- p + ggpubr::stat_compare_means(method = method, ...)
        } else {
          p_df <- get_adj_p(data,
            .col = var_name2, .grp = "subtype",
            method = method, p.adjust.method = p.adjust.method, ...
          )

          # p <- p + ggpubr::stat_compare_means(
          #   comparisons = my_comparisons,
          #   ...
          # )
          p <- p + ggpubr::stat_pvalue_manual(p_df, label = "p.adj")
        }
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
      ca_comb <- cowplot::plot_grid(plotlist = ca_res, align = "hv", ncol = 3)
    } else {
      ca_comb <- cowplot::plot_grid(plotlist = ca_res, align = "hv", ncol = 4)
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
        plotlist = co_res, align = "hv",
        axis = "l", ncol = 3
      )
    } else {
      co_comb <- cowplot::plot_grid(
        plotlist = co_res, align = "hv",
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

#' Plot signature enrichment result
#'
#' A modified version of [maftools::plotEnrichmentResults].
#'
#' @inheritParams maftools::plotEnrichmentResults
#' @param enrich_res results from [sig_group_enrichment],
#' or [maftools::clinicalEnrichment] or [maftools::signatureEnrichment].
#' @param addGeneDist Distance added to Gene label. Default is 0.
#' @param title title for plot.
#' @param showLegend if `TRUE`, show legend.
#' @return a base plot
#' @export
draw_sig_enrichment <-
  function(enrich_res, pVal = 0.05, cols = NULL, addGeneDist = 0, annoFontSize = 0.8,
             geneFontSize = 0.8, legendFontSize = 0.8, title = NULL, showTitle = TRUE, showLegend = TRUE) {
    # A modified version of maftools::plotEnrichmentResults
    binconf <- "maftools" %:::% "binconf"
    res <- enrich_res$groupwise_comparision
    plot.dat <- data.table::data.table(
      Hugo_Symbol = res$Hugo_Symbol,
      g1_muts = as.numeric(sapply(strsplit(
        x = res$n_mutated_group1,
        split = " of "
      ), "[[", 1)), g1_tot = as.numeric(sapply(strsplit(
        x = res$n_mutated_group1,
        split = " of "
      ), "[[", 2)), g2_muts = as.numeric(sapply(strsplit(
        x = res$n_mutated_group2,
        split = " of "
      ), "[[", 1)), g2_tot = as.numeric(sapply(strsplit(
        x = res$n_mutated_group2,
        split = " of "
      ), "[[", 2)), P_value = res$p_value,
      Group1 = res$Group1, Group2 = "Res"
    )
    plot.dat <- plot.dat[P_value < pVal]
    if (nrow(plot.dat) < 1) {
      stop(paste0(
        "No significant associations found at p-value < ",
        pVal
      ))
    }
    conf_int_g1 <- lapply(1:nrow(plot.dat), function(i) {
      as.data.frame(binconf(x = plot.dat[i, g1_muts], n = plot.dat[
        i,
        g1_tot
      ], alpha = pVal))
    })
    conf_int_g1 <- data.table::rbindlist(l = conf_int_g1)
    conf_int_g2 <- lapply(1:nrow(plot.dat), function(i) {
      as.data.frame(binconf(x = plot.dat[i, g2_muts], n = plot.dat[
        i,
        g2_tot
      ], alpha = pVal))
    })
    conf_int_g2 <- data.table::rbindlist(l = conf_int_g2)
    plot.dat$g1_muts_fract <- apply(plot.dat, 1, function(x) round(as.numeric(x[2]) / as.numeric(x[3]),
        digits = 3
      ))
    plot.dat$g2_muts_fract <- apply(plot.dat, 1, function(x) round(as.numeric(x[4]) / as.numeric(x[5]),
        digits = 3
      ))
    plot.dat[, `:=`(g1_title, paste0(g1_muts, "/", g1_tot))]
    plot.dat[, `:=`(g2_title, paste0(g2_muts, "/", g2_tot))]
    if (is.null(cols)) {
      cols <- c(
        RColorBrewer::brewer.pal(n = 9, name = "Set1"),
        RColorBrewer::brewer.pal(n = 8, name = "Dark2"),
        RColorBrewer::brewer.pal(n = 8, name = "Accent")
      )
      cols <- cols[1:length(unique(plot.dat$Group1))]
      names(cols) <- as.character(unique(plot.dat$Group1))
    }
    else {
      names(cols) <- as.character(enrich_res$cf_sizes$cf)
    }
    bar.cols <- cols[plot.dat[, Group1]]
    legend.cols <- cols[unique(plot.dat[, Group1])]
    add_legend <- function(...) {
      opar <- par(
        fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
        mar = c(0, 0, 0, 0), new = TRUE
      )
      on.exit(par(opar))
      plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
      legend(...)
    }
    yl_max <- max(rbind(conf_int_g1, conf_int_g2), na.rm = TRUE)
    if (yl_max <= 1) {
      yl_max <- 1
    }
    data.table::setDF(x = conf_int_g1)
    data.table::setDF(x = conf_int_g2)
    par(
      bty = "n", mgp = c(0.5, 0.5, 0), las = 1, tcl = -0.25,
      font.main = 4, xpd = TRUE, mar = c(4, 3, 3.5, 1)
    )
    b <- barplot(height = plot.dat$g1_muts_fract, ylim = c(
      -1.25,
      yl_max
    ), axes = FALSE, border = 0.1, col = bar.cols)
    axis(
      side = 2, at = seq(-1, 1, 0.25), labels = c(rev(seq(
        0,
        1, 0.25
      )), seq(0, 1, 0.25)[2:5]), lwd = 1.2, font.axis = 2,
      cex = 1.5, font = 1
    )
    for (i in 1:nrow(conf_int_g1)) {
      segments(x0 = b[i, 1], y0 = conf_int_g1[i, 2], x1 = b[
        i,
        1
      ], y1 = conf_int_g1[i, 3], lwd = 1.5)
    }
    text(b, conf_int_g1$Upper + 0.03, plot.dat$g1_title,
      cex = annoFontSize,
      las = 2, srt = 90, adj = 0, xpd = TRUE, font = 1
    )
    b2 <- barplot(
      height = -plot.dat$g2_muts_fract, add = TRUE,
      axes = FALSE, border = 0.1
    )
    for (i in 1:nrow(conf_int_g2)) {
      segments(x0 = b2[i, 1], y0 = -conf_int_g2[i, 2], x1 = b2[
        i,
        1
      ], y1 = -conf_int_g2[i, 3], lwd = 1.5)
    }
    text(b, -conf_int_g2$Upper - 0.03, plot.dat$g2_title,
      cex = annoFontSize,
      las = 2, srt = 90, adj = 1, xpd = TRUE, font = 1
    )
    text(b, -0.75 + (-addGeneDist), plot.dat$Hugo_Symbol,
      cex = geneFontSize,
      las = 2, srt = 90, adj = 1, xpd = TRUE, font = 3
    )
    if (length(legend.cols) <= 4) {
      n_col <- 1
    }
    else {
      n_col <- (length(legend.cols) %/% 4) + 1
    }

    if (showLegend) {
      legend(
        x = 0, y = -1.1, pt.lwd = 2, ncol = n_col, legend = c(
          names(legend.cols),
          "Rest"
        ), fill = c(legend.cols, "gray70"), bty = "n",
        cex = legendFontSize, border = NA, xpd = TRUE, text.font = 3
      )
    }
    if (showTitle) {
      if (is.null(title)) {
        title(
          main = enrich_res$clinicalFeature, adj = 0, cex.main = 1,
          outer = FALSE, font.main = 1
        )
      } else {
        title(
          main = title, adj = 0, cex.main = 1,
          outer = FALSE, font.main = 1
        )
      }
    }
  }

#' Filter enrichment result
#'
#' @param enrich results from [sig_group_enrichment],
#' or [maftools::clinicalEnrichment] or [maftools::signatureEnrichment].
#' @param subset logical expression used to subset data according to `type`.
#' @param type data type to subset, either 'groupwise' or 'pairwise'.
#' @return a filtered enrichment result
#' @export
filter_enrichment <- function(enrich, subset = TRUE, type = c("groupwise", "pairwise")) {
  type <- match.arg(type)
  enclos <- parent.frame()
  subset <- substitute(subset)
  if (type == "groupwise") {
    data <- enrich$groupwise_comparision
    row_selector <- eval(subset, data, enclos)
    data <- data[row_selector, ]
    enrich$groupwise_comparision <- data
  } else {
    data <- enrich$pairwise_comparision
    row_selector <- eval(subset, data, enclos)
    data <- data[row_selector, ]
    enrich$pairwise_comparision <- data
  }

  enrich
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
    "context",
    "P_value",
    "g1_muts",
    "g1_title",
    "g1_tot",
    "g2_muts",
    "g2_title",
    "g2_tot",
    "components",
    "stats"
  )
)
