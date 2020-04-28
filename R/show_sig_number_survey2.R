#' Show Comprehensive Signature Number Survey
#'
#' This function is modified from **NMF** package to
#' better help users to explore survey of signature number.
#'
#' @param x a `data.frame` or `NMF.rank` object obtained from [sig_estimate()].
#' @param y for random simulation,
#' a `data.frame` or `NMF.rank` object obtained from [sig_estimate()].
#' @inheritParams NMF::nmfEstimateRank
#'
#' @return a `ggplot` object
#' @export
#'
#' @examples
#' \donttest{
#' load(system.file("extdata", "toy_copynumber_tally_M.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' library(NMF)
#' cn_estimate <- sig_estimate(cn_tally_M$nmf_matrix,
#'   cores = 1, nrun = 5,
#'   verbose = TRUE,
#'   keep_nmfObj = TRUE
#' )
#'
#' # Show from data.frame obtained by sig_estimate()
#' show_sig_number_survey2(cn_estimate$survey)
#' show_sig_number_survey2(cn_estimate$survey, y = cn_estimate$survey.random)
#'
#' # Show directly from NMF.rank object
#' show_sig_number_survey2(cn_estimate$nmfEstimate)
#' show_sig_number_survey2(cn_estimate$nmfEstimate, y = cn_estimate$nmfEstimate.random)
#' }
show_sig_number_survey2 <- function(x, y = NULL, what = c(
                                      "all", "cophenetic", "rss", "residuals",
                                      "dispersion", "evar", "sparseness", "sparseness.basis", "sparseness.coef",
                                      "silhouette", "silhouette.coef", "silhouette.basis", "silhouette.consensus"
                                    ),
                                    na.rm = FALSE, xlab = "Number of signature",
                                    ylab = "", main = "Signature number survey using NMF package") {

  # Useless, just store it in case I need
  # to modify in the future
  xname <- "x"
  yname <- "y"

  if (is.character(y) && missing(what)) {
    what <- y
    y <- NULL
  }
  what <- match.arg(what, several.ok = TRUE)
  if ("all" %in% what) {
    what <- c(
      "cophenetic", "rss", "residuals", "dispersion",
      "evar", "sparseness", "silhouette"
    )
  }
  .getvals <- function(x, xname) {
    measures <- x
    iwhat <- unlist(lapply(
      paste("^", what, sep = ""), grep,
      colnames(measures)
    ))
    if (na.rm) {
      measures <- measures[apply(measures, 1, function(row) !any(is.na(row[iwhat]))), ]
    }
    vals <- measures[, iwhat, drop = FALSE]
    x <- as.numeric(measures$rank)
    xlim <- range(x)
    measure.type <- setNames(
      rep("Best fit", ncol(measures)),
      colnames(measures)
    )
    cons.measures <- c(
      "silhouette.consensus", "cophenetic",
      "cpu.all"
    )
    measure.type[match(cons.measures, names(measure.type))] <- "Consensus"
    measure.type[grep("\\.coef$", names(measure.type))] <- "Coefficients"
    measure.type[grep("\\.basis$", names(measure.type))] <- "Basis"
    measure.type <- factor(measure.type)
    pdata <- tidyr::pivot_longer(cbind(rank = x, vals),
      cols = colnames(vals),
      names_to = "variable"
    )
    pdata$Type <- measure.type[as.character(pdata$variable)]
    pdata$Measure <- gsub("^([^.]+).*", "\\1", pdata$variable)
    pdata$Data <- xname
    pdata
  }

  if (inherits(x, "NMF.rank")) {
    x <- x$measure
    pdata <- .getvals(x, xname)
  } else {
    pdata <- .getvals(x, xname)
  }

  if (!is.null(y)) {
    if (inherits(y, "NMF.rank")) {
      y <- y$measure
      pdata.y <- .getvals(y, yname)
    } else {
      pdata.y <- .getvals(y, yname)
    }

    pdata <- rbind(pdata, pdata.y)
  }

  p <- ggplot(pdata, aes_string(x = "rank", y = "value")) +
    geom_line(aes_string(linetype = "Data", colour = "Type")) +
    geom_point(size = 2, aes_string(shape = "Data", colour = "Type")) +
    theme_bw() +
    scale_x_continuous(xlab, breaks = unique(pdata$rank)) +
    scale_y_continuous(ylab) +
    ggtitle(main)
  if (!is(y, "NMF.rank")) {
    p <- p + scale_shape(guide = "none") + scale_linetype(guide = "none")
  }
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    send_stop("Please install 'RColorBrewer' package firstly.")
  }
  myColors <- RColorBrewer::brewer.pal(5, "Set1")
  names(myColors) <- levels(pdata$Type)
  p <- p + scale_colour_manual(name = "Measure type", values = myColors)
  p <- p + facet_wrap(~Measure, scales = "free")
  p
}
