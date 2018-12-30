# Visualization Part ------------------------------------------------------
#' Plot CNV signatures
#'
#' Plot CNV signature from multiple result objects from \code{VSHunter} package.
#'
#' @param Res a result \code{list} generated from CNV signature pipeline or a \code{NMFfitX} object from \code{NMF} package.
#' @param contributions if \code{TRUE}, plot contributions instead of signatures.
#' @param color user custom colors.
#' @param patient_order user custom patient order.
#' @param font_size font size.
#' @param show_title if \code{TRUE}, show title.
#' @param axis_lwd line width of axis.
#' @param title_size title size.
#' @param show_barcodes if \code{TRUE}, show barcode.
#' @param barcode_size barcode size when plot contribution.
#' @param yaxisLim limit of y axis.
#' @param ... other options pass to \code{barplot}.
#'
#' @return Nothing
#' @family CNV analysis functions
#' @export
#'
cnv_plotSignatures = function(Res = NULL,
                              contributions = FALSE,
                              color = NULL,
                              patient_order = NULL,
                              font_size = 1.2,
                              show_title = TRUE,
                              axis_lwd = 2,
                              title_size = 1.3,
                              show_barcodes = FALSE,
                              barcode_size = 0.5,
                              yaxisLim = 0.3,
                              ...) {
  # modify from https://github.com/PoisonAlien/maftools/blob/master/R/plotSignatures.R

  # input can be multipe objects from different functions
  # basically, we only use w and h matrix of NMF

  # assuming input is a basic NMF result
  if (inherits(Res, "NMFfitX1")) {
    #-- Signatures
    w = NMF::basis(Res)
    w = apply(w, 2, function(x)
      x / sum(x)) # Scale signatures (basis)
    colnames(w) = paste("Signature", 1:ncol(w), sep = "_")

    #-- Contributions
    h = NMF::coef(Res)
    h = apply(h, 2, function(x)
      x / sum(x)) # Scale contributions (coefs)

    rownames(h) = paste("Signature", 1:ncol(w), sep = "_")
  } else if (is.list(Res)) {
    if (!all(c("signature", "exposure") %in% names(Res)))
      stop("signature and exposure elements must in input list.")
    w = Res$signature
    w = apply(w, 2, function(x)
      x / sum(x)) # Scale signatures (basis)

    if (is.null(colnames(w))) {
      colnames(w) = paste("Signature", 1:ncol(w), sep = "_")
    }

    h = Res$exposure$relative_exposure
    if (is.null(rownames(h))) {
      rownames(h) = paste("Signature", 1:ncol(w), sep = "_")
    }
  }

  contrib = h

  if (contributions) {
    contribt = t(contrib)
    #calculate sd
    if (!is.null(patient_order)) {
      contribt = contribt[patient_order, ] #order on user-specified ordering of the genomes
    } else{
      contribt = contribt[order(contribt[, ncol(contribt)]), ] #order according to standard deviation
    }

    #contrib = t(contribt[,1:(ncol(contribt)-1)])
    contrib = t(contribt[, 1:(ncol(contribt))])

    cols = RColorBrewer::brewer.pal(n = 8, name = 'Set2')

    if (show_barcodes) {
      lo = layout(mat = matrix(data = c(1, 2), nrow = 2),
                  heights = c(6, 2))
      par(mar = c(6, 4, 2, 1))
      b = barplot(
        contrib,
        axes = FALSE,
        horiz = FALSE,
        col = cols,
        border = NA,
        names.arg = rep("", ncol(contrib))
      )
      axis(
        side = 1,
        at = b,
        labels = colnames(contrib),
        lwd = 2,
        cex.axis = barcode_size,
        las = 2,
        line = 1,
        hadj = 0.8,
        font = 2,
        tick = FALSE
      )
      axis(
        side = 2,
        at = seq(0, 1, 0.25),
        lwd = 3,
        font = 2,
        las = 2,
        cex.axis = 0.9
      )
      mtext(
        text = "Signature exposures",
        side = 2,
        font = 2,
        cex = 1,
        line = 2.8,
        srt = 90
      ) # strange srt
      par(mar = rep(2, 4))
      plot.new()
      #par(mar = c(2, 3, 0, 0))
      legend(
        x = "left",
        legend = rownames(contrib),
        col = cols[1:nrow(contrib)],
        border = NA,
        bty = "n",
        pch = 15,
        xpd = TRUE,
        ncol = 1,
        cex = 1.2,
        pt.cex = 1.5,
        horiz = TRUE
      )
    } else{
      lo = layout(mat = matrix(data = c(1, 2), nrow = 2),
                  heights = c(6, 2))
      par(mar = c(3, 4, 2, 1))
      b = barplot(
        contrib,
        axes = FALSE,
        horiz = FALSE,
        col = cols,
        border = NA,
        names.arg = rep("", ncol(contrib))
      )
      axis(
        side = 2,
        at = seq(0, 1, 0.25),
        lwd = 3,
        font = 2,
        las = 2,
        cex.axis = 0.9
      )
      mtext(
        text = "Signature exposure",
        side = 2,
        font = 2,
        cex = 1,
        line = 2.8,
        srt = 90
      )
      par(mar = rep(2, 4))
      plot.new()
      #par(mar = c(2, 3, 0, 0))
      legend(
        x = "left",
        legend = rownames(contrib),
        col = cols[1:nrow(contrib)],
        border = NA,
        bty = "n",
        pch = 15,
        xpd = TRUE,
        ncol = 1,
        cex = 1.2,
        pt.cex = 1.5,
        horiz = TRUE
      )
    }
  } else{
    plotData = as.data.frame(t(w))
    nsigs = nrow(plotData)

    if (is.null(color)) {
      #color = c("blue","black","red","gray","green","maroon")
      color = c(
        'coral4',
        'lightcyan4',
        'deeppink3',
        'lightsalmon1',
        'forestgreen',
        'cornflowerblue'
      )
    }

    cnames = colnames(plotData)

    # order the features
    plotData = plotData[, c(
      grep("bp10MB", cnames, value = TRUE),
      grep("copynumber", cnames, value = TRUE),
      grep("changepoint", cnames, value = TRUE),
      grep("bpchrarm", cnames, value = TRUE),
      grep("osCN", cnames, value = TRUE),
      grep("segsize", cnames, value = TRUE)
    )]


    len_c1 = sum(grepl("bp10MB", cnames))
    len_c2 = sum(grepl("copynumber", cnames))
    len_c3 = sum(grepl("changepoint", cnames))
    len_c4 = sum(grepl("bpchrarm", cnames))
    len_c5 = sum(grepl("osCN", cnames))
    len_c6 = sum(grepl("segsize", cnames))

    len_seqc = c(len_c1,
                 len_c2,
                 len_c3,
                 len_c4,
                 len_c5,
                 len_c6)

    colors = rep(color, times = len_seqc)

    len_seq = c(0L)
    for (i in 1:length(len_seqc)) {
      s = len_seqc[i] * 2 + len_seq[i]
      len_seq = c(len_seq, s)
    }

    par(
      mfrow = c(nsigs, 1),
      oma = c(5, 4, 0, 0) + 0.1,
      mar = c(0, 0, 2.5, 0) + 0.1,
      las = 1,
      tcl = -.25,
      font.main = 4,
      xpd = NA
    )

    for (i in 1:nsigs) {
      ae = rownames(plotData)[i]
      d = as.matrix(plotData[i, ])
      if (is.na(yaxisLim)) {
        bh = ceiling(max(d, na.rm = TRUE) * 10) / 10 #Bar height
      } else{
        bh = 0.3
      }

      barplot(
        d,
        xaxt = "n",
        yaxt = "n",
        col = colors,
        beside = TRUE,
        ylim = c(-0.1, bh),
        cex.main = 1,
        border = NA,
        font.axis = 2,
        font.lab = 2,
        adj = 0.25,
        ...
      )
      if (show_title) {
        title(
          main = ae,
          cex.main = title_size,
          line = 0,
          adj = 0.98
        )
      }

      #mtext(text = ae, side = 1, line = 2, font = 1, cex = 0.5, at = 0.3)
      axis(
        side = 2,
        at = seq(0, bh, 0.1),
        pos = -2,
        las = 2,
        lwd = axis_lwd,
        hadj = 1.1,
        font = 2,
        cex.axis = font_size
      )
      #abline(h = seq(0, 0.3, 0.1),lty=2,lwd=0.3, col = 'gray70')
      rect(
        xleft = len_seq,
        ybottom = -0.05,
        xright = ncol(plotData) * 2,
        ytop = -0.02,
        col = color,
        border = 'gray70'
      )
      if (i == nsigs) {
        text(
          labels = c(
            "bp10MB",
            "copynumber",
            "changepoint",
            "bpchrarm",
            "osCN",
            "segsize"
          ),
          y = rep(-0.1, 6),
          x = len_seq[2:7] - len_seqc,
          cex = font_size,
          font = 1.8,
          font.lab = 2,
          pos = 1.2,
          srt = 30
        )
      }
    }
  }


}


#' Plot copy number feature distribution
#'
#' @param features a \code{list} generate from \code{get_features} function.
#' @param ylab lab of y axis.
#' @param ... other options pass to \code{\link[cowplot]{plot_grid}} function of \code{cowplot} package.
#'
#' @return a ggplot object
#' @import cowplot
#' @family CNV analysis functions
#' @export
#'
cnv_plotFeatureDistribution = function(features, ylab = "", ...) {
  features = lapply(features, function(x) {
    x$value = as.numeric(x$value)
    return(x)
  })

  #requireNamespace("cowplot")
  #cowplot::theme_cowplot()
  #ggplot2::theme_set(cowplot::theme_cowplot(font_size = 12))
  p_1 = ggplot(data = features$segsize, aes(x = log10(value))) +
    geom_line(stat = "density") + labs(x = "Segment size (log10 based)", y = ylab) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + theme_cowplot(font_size = 12)
  #p_1 = p_1 + scale_x_continuous(breaks = 10 ^c(7, 8),
  #                               labels = scales::trans_format("log10", scales::math_format(10^.x)))

  p_2 = ggplot(data = features$copynumber, aes(x = value)) +
    geom_line(stat = "density") + labs(x = "Copy number", y = ylab) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + theme_cowplot(font_size = 12)
  p_3 = ggplot(data = features$changepoint, aes(x = value)) +
    geom_line(stat = "density") + labs(x = "Copy number changepoint", y = ylab) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + theme_cowplot(font_size = 12)

  p_4 = ggplot(data = features$bp10MB, aes(x = value)) +
    geom_bar(stat = "count") + labs(x = "Breakpoint count per 10MB", y = ylab) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + theme_cowplot(font_size = 12)
  p_5 = ggplot(data = features$bpchrarm, aes(x = value)) +
    geom_bar(stat = "count") + labs(x = "Breakpoint count per chr arm", y = ylab) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + theme_cowplot(font_size = 12)

  p_6 = ggplot(data = features$osCN, aes(x = value)) +
    geom_bar(stat = "count") + labs(x = "Oscilating CN chain length", y = ylab) +
    theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) + theme_cowplot(font_size = 12)
  # p_6 = ggplot(data = features$osCN, aes(x = as.factor(value))) +
  #     geom_bar(stat = "count") + labs(x = "Oscilating CN chain length", y = ylab) + theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"))

  # osCN_tab = length(table(features$osCN$value))
  #
  # if (osCN_tab > 15 & osCN_tab <= 20) {
  #     p_6 = p_6 + theme(axis.text.x = element_text(size = 9))
  # } else if (osCN_tab > 20 & osCN_tab <= 30) {
  #     p_6 = p_6 + theme(axis.text.x = element_text(size = 7))
  # } else if (osCN_tab > 30 & osCN_tab <= 40) {
  #     p_6 = p_6 + theme(axis.text.x = element_text(size = 4))
  # } else if (osCN_tab > 40) {
  #     p_6 = p_6 + theme(axis.text.x = element_text(size = 3))
  # }

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


#' Plot mixture fit model components
#'
#' @inheritParams cnv_plotFeatureDistribution
#' @param components a \code{list} contain \code{flexmix} object of copy-number features, obtain this
#' from \code{cnv_fitMixModels} function or use pre-compiled components data which come from CNV signature paper
#' https://www.nature.com/articles/s41588-018-0179-8 (set this argument as \code{NULL}).
#' @param ... other options pass to \code{\link[cowplot]{plot_grid}} function of \code{cowplot} package.
#'
#' @return a ggplot object
#' @import cowplot
#' @export
#' @family CNV analysis functions
#'
cnv_plotMixComponents = function(features, components, ...) {
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

    p = p + xlab(xlab) + theme_cowplot(font_size = 12)
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

    p = p + xlab(xlab) + theme_cowplot(font_size = 12)
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



#' Plot copy number distribution either by length or chromosome
#'
#' @param data a \code{data.frame}, result of \code{get_LengthFraction} function.
#' @param rm_normal logical. Whether removel normal copy (i.e. "segVal" equals 2), default is \code{TRUE}.
#' @param mode either "ld" for distribution by CN length or "cd" for distribution by chromosome.
#' @param fill when \code{mode} is "cd", if \code{fill} is \code{TRUE}, plot percentage instead of count.
#' @param scale_chr logical. If `TRUE`, normalize count to per Megabase unit.
#' @inheritParams get_features
#' @author Shixiang Wang <w_shixiang@163.com>
#'
#' @return a ggplot object
#' @import cowplot
#' @export
#' @family CNV analysis functions
#' @examples
#' \dontrun{
#' cnv_plotLengthSummary(data, mode = "cd")
#' }
cnv_plotDistributionProfile = function(data,
                                       rm_normal = TRUE,
                                       mode = c("ld", "cd"),
                                       fill = FALSE,
                                       scale_chr = TRUE,
                                       genome_build = c("hg19", "hg38")) {
  stopifnot(is.logical(rm_normal),
            is.data.frame(data),
            is.logical(fill))
  mode = match.arg(mode)
  genome_build = match.arg(genome_build)

  # if remove normal copy number segments
  if (rm_normal) {
    if (!"segVal" %in% colnames(data)) {
      stop("'segVal' must be provided as a column.")
    }
    data = base::subset(data, data[, "segVal"] != 2)
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
      data = subset(data, chromosome %in% as.character(1:22))
    } else if (is.integer(data$chromosome[1])) {
      data = subset(data, chromosome %in% 1:22)
    }

    # only keep chr 1 to 22
    data$chromosome = factor(as.character(data$chromosome), levels = as.character(1:22))

    # only keep p, q, pq
    data$location = factor(sub("[0-9]*", "", data$location),
                           levels = c("p", "pq", "q"))
    if (sum(!(data$location %in% c("p", "pq", "q"))) > 0){
      message("Discarding segments which located in centromere region...")
      data = subset(data, location %in% c("p", "pq", "q"))
    }

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

      q$x = factor(as.character(q$x), levels = as.character(1:22))
      q$fill = factor(q$fill, levels = c("#F8766D", "#00BA38", "#619CFF"))

      chrlen$chrom = gsub(
        pattern = "chr",
        replacement = "",
        x = chrlen$chrom
      )
      q = merge(q, chrlen, by.x = "x", by.y = "chrom")
      q[["count"]] = 1000000 * (q[["count"]] / q[["size"]])
      # data.table::setDT(data)
      # data.table::setDT(chrlen)
      # chrlen$chrom = gsub(pattern = "chr", replacement = "", x = chrlen$chrom)
      # data = merge(data, chrlen, by.x = "chromosome", by.y = "chrom")
      # data[, c(), by = chromosome]
      # plot
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





utils::globalVariables(
  c(
    ".x",
    "aes",
    "chrom",
    "chromosome",
    "element_text",
    "geom_bar",
    "geom_line",
    "ggplot",
    "labs",
    "scale_x_continuous",
    "stat_function",
    "theme",
    "theme_classic",
    "value",
    "x",
    "ylab",
    "temp",
    "fraction",
    "..density..",
    "geom_histogram",
    "location",
    "xlab",
    "unit",
    "count",
    "scale_fill_discrete",
    "ggplot_build"
  )
)
