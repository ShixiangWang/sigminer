#' Estimate Signature Number
#'
#' Use **NMF** package to evaluate the optimal number of signatures.
#' This is used along with [sig_extract].
#' Users should `library(NMF)` firstly. If NMF objects are returned,
#' the result can be further visualized by NMF plot methods like
#' `NMF::consensusmap()` and `NMF::basismap()`.
#'
#' The most common approach is to choose the smallest rank for which cophenetic correlation coefficient
#' starts decreasing (Used by this function). Another approach is to choose the rank for which the plot
#' of the residual sum of squares (RSS) between the input matrix and its estimate shows an inflection point.
#' More custom features please directly use [NMF::nmfEstimateRank].
#'
#' @name sig_estimate
#' @param nmf_matrix a `matrix` used for NMF decomposition with rows indicate samples and columns indicate components.
#' @param range a `numeric` vector containing the ranks of factorization to try. Note that duplicates are removed
#' and values are sorted in increasing order. The results are notably returned in this order.
#' @param keep_nmfObj default is `FALSE`, if `TRUE`, keep NMF objects from runs, and the result may be huge.
#' @param nrun a `numeric` giving the number of run to perform for each value in `range`, `nrun` set to 30~50 is
#' enough to achieve robust result.
#' @param what a character vector whose elements partially match one of the following item, which correspond to
#' the measures computed by summary on each multi-run NMF result: ‘all’, ‘cophenetic’, ‘rss’, ‘residuals’,
#'  ‘dispersion’, ‘evar’, ‘silhouette’ (and more specific .coef, .basis, .consensus), ‘sparseness’
#'  (and more specific .coef, .basis). It specifies which measure must be plotted
#'  (what='all' plots all the measures).
#'
#' @param cores number of cpu cores to run NMF.
#' @param seed specification of the starting point or seeding method, which will compute a starting point,
#'  usually using data from the target matrix in order to provide a good guess.
#' @param use_random Should generate random data from input to test measurements. Default is `TRUE`.
#' @param save_plots if `TRUE`, save signature number survey plot to local machine.
#' @param plot_basename when save plots, set custom basename for file path.
#' @param method specification of the NMF algorithm. Use 'brunet' as default.
#' Available methods for NMF decompositions are 'brunet', 'lee', 'ls-nmf', 'nsNMF', 'offset'.
#' @param verbose if `TRUE`, print extra message.
#' @author Shixiang Wang
#' @references Gaujoux, Renaud, and Cathal Seoighe. "A flexible R package for nonnegative matrix factorization." BMC bioinformatics 11.1 (2010): 367.
#' @return - sig_estimate: a `list` contains information of NMF run and rank survey.
#' @export
#' @seealso [sig_extract] for extracting signatures using **NMF** package, [sig_auto_extract] for
#' extracting signatures using automatic relevance determination technique.
sig_estimate <-
  function(nmf_matrix,
           range = 2:5,
           nrun = 10,
           use_random = FALSE,
           method = "brunet",
           seed = 123456,
           cores = 1,
           keep_nmfObj = FALSE,
           save_plots = FALSE,
           plot_basename = file.path(tempdir(), "nmf"),
           what = "all",
           verbose = FALSE) {
    eval(parse(text = "suppressMessages(library('NMF'))"))
    if (cores > 1) cores <- min(cores, future::availableCores())
    mat <- t(nmf_matrix)

    ii <- colSums(mat) < 0.01
    if (any(ii)) {
      message(
        "The follow samples dropped due to null catalogue:\n\t",
        paste0(colnames(mat)[ii], collapse = ", ")
      )
      mat <- mat[, !ii, drop = FALSE]
    }

    # To avoid error due to NMF
    mat <- check_nmf_matrix(mat)

    if (cores > 1) {
      estim.r <-
        NMF::nmfEstimateRank(
          mat,
          range,
          method = method,
          nrun = nrun,
          verbose = verbose,
          seed = seed,
          .opt = paste0("p", cores)
        )
    } else {
      estim.r <-
        NMF::nmfEstimateRank(
          mat,
          range,
          method = method,
          nrun = nrun,
          verbose = verbose,
          seed = seed
        )
    }

    nmf.sum <- NMF::summary(estim.r) # Get summary of estimates
    if (verbose) {
      message("Estimation of rank based on observed data.")
      print(nmf.sum)
    }

    if (use_random) {
      if (verbose) message("Generating random matrix and run NMF...")
      V.random <- NMF::randomize(mat)

      if (cores > 1) {
        estim.r.random <-
          NMF::nmfEstimateRank(
            V.random,
            range,
            method = method,
            nrun = nrun,
            verbose = verbose,
            seed = seed,
            .opt = paste0("p", cores)
          )
      } else {
        estim.r.random <-
          NMF::nmfEstimateRank(
            V.random,
            range,
            method = method,
            nrun = nrun,
            verbose = verbose,
            seed = seed
          )
      }

      nmf.sum.random <- NMF::summary(estim.r.random) # Get summary of estimates
      if (verbose) {
        message("Estimation of rank based on random data.")
        print(nmf.sum.random)
      }
    } else {
      estim.r.random <- NULL
      nmf.sum.random <- NULL
    }

    if (save_plots) {
      if (use_random) {
        p <- NMF::plot(
          estim.r,
          estim.r.random,
          what = what,
          xname = "Observed",
          yname = "Randomised",
          xlab = "Number of signature",
          main = "Signature number survey using NMF package"
        )
      } else {
        p <- NMF::plot(
          estim.r,
          what = what,
          xlab = "Number of signature",
          main = "Signature number survey using NMF package"
        )
      }

      destdir <- dirname(plot_basename)
      if (!dir.exists(destdir)) dir.create(destdir, recursive = TRUE)
      pdf(
        paste0(plot_basename, "_survey.pdf"),
        bg = "white",
        pointsize = 9,
        width = 6 + abs((nrow(nmf.sum) - 6) / 3),
        height = 6,
        paper = "special"
      )
      print(p)
      dev.off()
      if (verbose) message("Created ", paste0(plot_basename, "_survey.pdf"))
    }

    if (keep_nmfObj) {
      res <- list(
        nmfEstimate = estim.r,
        nmfEstimate.random = estim.r.random,
        survey = nmf.sum,
        survey.random = nmf.sum.random
      )
    } else {
      res <- list(
        survey = nmf.sum,
        survey.random = nmf.sum.random
      )
    }

    class(res) <- "Survey"
    res
  }

#' Show Simplified Signature Number Survey
#'
#' `sig_estimate()` shows comprehensive rank survey generated by
#' **NMF** package, sometimes
#' it is hard to consider all measures. `show_sig_number_survey()` provides a
#' one or two y-axis visualization method to help users determine
#' the optimal signature number (showing both
#' stability ("cophenetic") and error (RSS) at default).
#' Users can also set custom measures to show.
#'
#' @rdname sig_estimate
#' @param object a `Survey` object generated from [sig_estimate], or
#' a `data.frame` contains at least rank columns and columns for
#' one measure.
#' @param x column name for x axis.
#' @param left_y column name for left y axis.
#' @param right_y column name for right y axis.
#' @param left_name label name for left y axis.
#' @param right_name label name for right y axis.
#' @param left_color color for left axis.
#' @param right_color color for right axis.
#' @param left_shape,right_shape,shape_size shape setting.
#' @param highlight a `integer` to highlight a `x`.
#' @return - show_sig_number_survey: a `ggplot` object
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
#'   verbose = TRUE
#' )
#'
#' p <- show_sig_number_survey2(cn_estimate$survey)
#' p
#'
#' # Show two measures
#' show_sig_number_survey(cn_estimate)
#' # Show one measure
#' p1 <- show_sig_number_survey(cn_estimate, right_y = NULL)
#' p1
#' p2 <- add_h_arrow(p, x = 4.1, y = 0.953, label = "selected number")
#' p2
#'
#' # Show data from a data.frame
#' p3 <- show_sig_number_survey(cn_estimate$survey)
#' p3
#' # Show other measures
#' head(cn_estimate$survey)
#' p4 <- show_sig_number_survey(cn_estimate$survey,
#'   right_y = "dispersion",
#'   right_name = "dispersion"
#' )
#' p4
#' p5 <- show_sig_number_survey(cn_estimate$survey,
#'   right_y = "evar",
#'   right_name = "evar"
#' )
#' p5
#' }
#' @testexamples
#' expect_is(p, "ggplot")
#' expect_is(p1, "ggplot")
#' expect_is(p2, "ggplot")
#' expect_is(p3, "ggplot")
#' expect_is(p4, "ggplot")
#' expect_is(p5, "ggplot")
#' @seealso [sig_estimate] for estimating signature number for [sig_extract],
#' [show_sig_number_survey2] for more visualization method.
show_sig_number_survey <- function(object, x = "rank",
                                   left_y = "cophenetic", right_y = "rss",
                                   left_name = left_y, right_name = toupper(right_y),
                                   left_color = "black", right_color = "red",
                                   left_shape = 16, right_shape = 18,
                                   shape_size = 4, highlight = NULL) {
  stopifnot(is(object, "Survey") | is.data.frame(object))
  if (is(object, "Survey")) {
    survey <- object$survey
  } else {
    survey <- object
  }

  if (is.null(right_y)) {
    # Show one-axis plot
    p <- ggplot(data = survey) +
      geom_point(aes_string(x = x, y = left_y), color = left_color, shape = left_shape, size = shape_size) +
      geom_line(aes_string(x = x, y = left_y), color = left_color) +
      scale_x_continuous(breaks = unique(survey[[x]]), labels = unique(survey[[x]])) +
      scale_y_continuous(name = left_name) +
      cowplot::theme_cowplot() +
      xlab("Total signatures")
  } else {
    # Show two-axis plot
    survey$new_right <- norm2rg(
      survey[[right_y]],
      range(survey[[left_y]])
    )

    p <- ggplot(data = survey) +
      geom_point(aes_string(x = x, y = left_y), color = left_color, shape = left_shape, size = shape_size) +
      geom_point(aes_string(x = x, y = "new_right"), color = right_color, shape = right_shape, size = shape_size) +
      geom_line(aes_string(x = x, y = left_y), color = left_color) +
      geom_line(aes_string(x = x, y = "new_right"), color = right_color) +
      scale_x_continuous(breaks = unique(survey[[x]]), labels = unique(survey[[x]])) +
      scale_y_continuous(
        name = left_name,
        sec.axis = sec_axis(~ norm2rg(., range(survey[[right_y]])),
          name = right_name
        )
      ) +
      cowplot::theme_cowplot() +
      xlab("Total signatures") +
      theme(
        axis.title.y.left = element_text(color = left_color),
        axis.text.y.left = element_text(color = left_color),
        axis.title.y.right = element_text(color = right_color),
        axis.text.y.right = element_text(color = right_color),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)
      )
  }

  if (!is.null(highlight)) {
    p <- p + geom_vline(xintercept = highlight, linetype = 2, size = 1, color = "blue")
  }

  return(p)
}

#' Show Comprehensive Signature Number Survey
#'
#' `show_sig_number_survey2()` is modified from **NMF** package to
#' better help users to explore survey of signature number.
#'
#' @rdname sig_estimate
#' @param x a `data.frame` or `NMF.rank` object obtained from [sig_estimate()].
#' @param y for random simulation,
#' a `data.frame` or `NMF.rank` object obtained from [sig_estimate()].
#' @param what a character vector whose elements partially match one of the following item,
#' which correspond to the measures computed by `summary()` on each – multi-run – NMF result:
#' 'all', 'cophenetic', 'rss', 'residuals', 'dispersion', 'evar', 'silhouette'
#' (and more specific `*.coef`, `*.basis`, `*.consensus`), 'sparseness'
#' (and more specific `*.coef`, `*.basis`).
#' It specifies which measure must be plotted (what='all' plots all the measures).
#' @inheritParams NMF::nmfEstimateRank
#'
#' @return - show_sig_number_survey2: a `ggplot` object
#' @export
show_sig_number_survey2 <- function(x, y = NULL, what = c(
                                      "all", "cophenetic", "rss", "residuals",
                                      "dispersion", "evar", "sparseness", "sparseness.basis", "sparseness.coef",
                                      "silhouette", "silhouette.coef", "silhouette.basis", "silhouette.consensus"
                                    ),
                                    na.rm = FALSE, xlab = "Total signatures",
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
