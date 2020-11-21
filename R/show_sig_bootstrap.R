#' Show Signature Bootstrap Analysis Results
#'
#' See details for description.
#'
#' Functions:
#'
#' - [show_sig_bootstrap_exposure] - this function plots exposures from bootstrap samples with both dotted boxplot.
#' The optimal exposure (the exposure from original input) is shown as triangle point. **Only one sample can be plotted**.
#' - [show_sig_bootstrap_error] - this function plots decomposition errors from bootstrap samples with both dotted boxplot.
#' The error from optimal solution (the decomposition error from original input) is shown as triangle point. **Only one sample can be plotted**.
#' - [show_sig_bootstrap_stability] - this function plots the signature exposure instability for specified signatures. Currently,
#' the instability measure supports 3 types:
#'   - 'RMSE' for Mean Root Squared Error (default) of bootstrap exposures and original exposures for each sample.
#'   - 'CV' for  Coefficient of Variation (CV) based on RMSE (i.e. `RMSE / btExposure_mean`).
#'   - 'MAE' for Mean Absolute Error of bootstrap exposures and original exposures for each sample.
#'   - 'AbsDiff' for Absolute Difference between mean bootstram exposure and original exposure.
#'
#' @inheritParams ggpubr::ggboxplot
#' @inheritParams ggpubr::ggpar
#' @inheritParams ggplot2::geom_boxplot
#' @inheritParams sig_fit_bootstrap_batch
#' @param bt_result result object from [sig_fit_bootstrap_batch].
#' @param sample a sample id.
#' @param signatures signatures to show.
#' @param measure measure to estimate the exposure instability, can be one of 'RMSE', 'CV', 'MAE' and 'AbsDiff'.
#' @param dodge_width dodge width.
#' @param plot_fun set the plot function.
#' @param agg_fun set the aggregation function when `sample` is `NULL`.
#' @param highlight set the color for optimal solution. Default is "auto", which use the same color as
#' bootstrap results, you can set it to color like "red", "gold", etc.
#' @param highlight_size size for highlighting triangle, default is `4`.
#' @param ... other parameters passing to [ggpubr::ggboxplot] or [ggpubr::ggviolin].
#'
#' @name show_sig_bootstrap
#' @return a `ggplot` object
#' @seealso [sig_fit_bootstrap_batch], [sig_fit], [sig_fit_bootstrap]
#' @references Huang X, Wojtowicz D, Przytycka TM. Detecting presence of mutational signatures in cancer with confidence. Bioinformatics. 2018;34(2):330–337. doi:10.1093/bioinformatics/btx604
#' @examples
#' \donttest{
#' if (require("BSgenome.Hsapiens.UCSC.hg19")) {
#'   laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#'   laml <- read_maf(maf = laml.maf)
#'   mt_tally <- sig_tally(
#'     laml,
#'     ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
#'     use_syn = TRUE
#'   )
#'
#'   library(NMF)
#'   mt_sig <- sig_extract(mt_tally$nmf_matrix,
#'     n_sig = 3,
#'     nrun = 2,
#'     cores = 1
#'   )
#'
#'   mat <- t(mt_tally$nmf_matrix)
#'   mat <- mat[, colSums(mat) > 0]
#'   bt_result <- sig_fit_bootstrap_batch(mat, sig = mt_sig, n = 10)
#'   ## Parallel computation
#'   ## bt_result = sig_fit_bootstrap_batch(mat, sig = mt_sig, n = 10, use_parallel = TRUE)
#'
#'   ## At default, mean bootstrap exposure for each sample has been calculated
#'   p <- show_sig_bootstrap_exposure(bt_result, methods = c("QP"))
#'   ## Show bootstrap exposure (optimal exposure is shown as triangle)
#'   p1 <- show_sig_bootstrap_exposure(bt_result, methods = c("QP"), sample = "TCGA-AB-2802")
#'   p1
#'   p2 <- show_sig_bootstrap_exposure(bt_result,
#'     methods = c("QP"),
#'     sample = "TCGA-AB-3012",
#'     signatures = c("Sig1", "Sig2")
#'   )
#'   p2
#'
#'   ## Show bootstrap error
#'   ## Similar to exposure above
#'   p <- show_sig_bootstrap_error(bt_result, methods = c("QP"))
#'   p
#'   p3 <- show_sig_bootstrap_error(bt_result, methods = c("QP"), sample = "TCGA-AB-2802")
#'   p3
#'
#'   ## Show exposure (in)stability
#'   p4 <- show_sig_bootstrap_stability(bt_result, methods = c("QP"))
#'   p4
#'   p5 <- show_sig_bootstrap_stability(bt_result, methods = c("QP"), measure = "MAE")
#'   p5
#'   p6 <- show_sig_bootstrap_stability(bt_result, methods = c("QP"), measure = "AbsDiff")
#'   p6
#'   p7 <- show_sig_bootstrap_stability(bt_result, methods = c("QP"), measure = "CV")
#'   p7
#' } else {
#'   message("Please install package 'BSgenome.Hsapiens.UCSC.hg19' firstly!")
#' }
#' }
#' @testexamples
#' expect_s3_class(p1, "ggplot")
#' expect_s3_class(p2, "ggplot")
#' expect_s3_class(p3, "ggplot")
#' expect_s3_class(p4, "ggplot")
#' expect_s3_class(p5, "ggplot")
#' expect_s3_class(p6, "ggplot")
#' expect_s3_class(p7, "ggplot")
NULL

#' @rdname show_sig_bootstrap
#' @export
show_sig_bootstrap_exposure <- function(bt_result, sample = NULL, signatures = NULL,
                                        methods = "QP", plot_fun = c("boxplot", "violin"),
                                        agg_fun = c("mean", "median", "min", "max"),
                                        highlight = "auto", highlight_size = 4,
                                        palette = "aaas", title = NULL,
                                        xlab = FALSE, ylab = "Signature exposure", width = 0.3,
                                        dodge_width = 0.8, outlier.shape = NA,
                                        add = "jitter", add.params = list(alpha = 0.3),
                                        ...) {
  stopifnot(is.list(bt_result))
  plot_fun <- match.arg(plot_fun)
  agg_fun <- match.arg(agg_fun)
  agg_mode <- FALSE
  plot_fun <- switch(
    plot_fun,
    boxplot = ggpubr::ggboxplot,
    violin = ggpubr::ggviolin
  )

  timer <- Sys.time()
  send_info("Started.")
  on.exit(send_elapsed_time(timer))

  dat <- bt_result$expo
  rm(bt_result)
  dat <- dplyr::filter(dat, .data$method %in% methods)
  if (!is.null(sample)) {
    if (length(sample) > 1) {
      send_stop("Only one sample can be used!")
    }
    dat <- dplyr::filter(dat, .data$sample %in% .env$sample)
  } else {
    agg_mode <- TRUE
    send_success("NULL sample detected, aggregate mode enabled.")
    send_info("The summarized values will be stored in 'summary' element of result ggplot2 object.")
    agg_fun <- switch(agg_fun,
      mean = mean,
      median = median,
      min = min,
      max = max
    )
    dat <- dat %>%
      dplyr::mutate(type = ifelse(.data$type != "optimal", "bootstrap", .data$type)) %>%
      dplyr::group_by(.data$type, .data$method, .data$sample, .data$sig) %>%
      dplyr::summarise(exposure = agg_fun(.data$exposure, na.rm = TRUE)) %>%
      dplyr::ungroup()
  }
  if (!is.null(signatures)) {
    dat <- dplyr::filter(dat, .data$sig %in% signatures)
  }
  if (!nrow(dat) > 0) {
    send_stop("No data left to plot, could you check your input?")
  }

  if (is.null(title)) {
    if (agg_mode) {
      title <- "All samples"
    } else {
      title <- unique(dat$sample)
    }
  }

  send_info("Plotting.")
  ## Plotting
  p <- plot_fun(subset(dat, dat$type != "optimal"),
    x = "sig", y = "exposure", color = "method", outlier.shape = outlier.shape,
    palette = palette, width = width, add = add, add.params = add.params,
    title = title, xlab = xlab, ylab = ylab, ...
  )

  if (!agg_mode) {
    if (highlight == "auto") {
      p <- p + ggplot2::geom_point(
        data = subset(dat, dat$type == "optimal"),
        mapping = ggplot2::aes_string(x = "sig", y = "exposure", color = "method"),
        shape = 17, size = highlight_size,
        position = ggplot2::position_dodge2(width = dodge_width, preserve = "single")
      )
    } else {
      p <- p + ggplot2::geom_point(
        data = subset(dat, dat$type == "optimal"),
        mapping = ggplot2::aes_string(x = "sig", y = "exposure"),
        shape = 17, size = highlight_size, color = highlight,
        position = ggplot2::position_dodge2(width = dodge_width, preserve = "single")
      )
    }
  } else {
    p$summary <- dat
  }

  p
}


#' @rdname show_sig_bootstrap
#' @export
show_sig_bootstrap_error <- function(bt_result, sample = NULL,
                                     methods = "QP", plot_fun = c("boxplot", "violin"),
                                     agg_fun = c("mean", "median"),
                                     highlight = "auto", highlight_size = 4,
                                     palette = "aaas", title = NULL,
                                     xlab = FALSE, ylab = "Reconstruction error (L2 norm)", width = 0.3,
                                     dodge_width = 0.8, outlier.shape = NA,
                                     add = "jitter", add.params = list(alpha = 0.3),
                                     legend = "none",
                                     ...) {
  stopifnot(is.list(bt_result))
  plot_fun <- match.arg(plot_fun)
  agg_fun <- match.arg(agg_fun)
  agg_mode <- FALSE
  plot_fun <- switch(
    plot_fun,
    boxplot = ggpubr::ggboxplot,
    violin = ggpubr::ggviolin
  )

  timer <- Sys.time()
  send_info("Started.")
  on.exit(send_elapsed_time(timer))

  dat <- bt_result$error
  rm(bt_result)
  dat <- dplyr::filter(dat, .data$method %in% methods)
  if (!is.null(sample)) {
    if (length(sample) > 1) {
      send_stop("Only one sample can be used!")
    }
    dat <- dplyr::filter(dat, .data$sample %in% .env$sample)
  } else {
    agg_mode <- TRUE
    send_success("NULL sample detected, aggregate mode enabled.")
    send_info("The summarized values will be stored in 'summary' element of result ggplot2 object.")
    agg_fun <- switch(agg_fun,
      mean = mean,
      median = median,
      min = min,
      max = max
    )
    dat <- dat %>%
      dplyr::mutate(type = ifelse(.data$type != "optimal", "bootstrap", .data$type)) %>%
      dplyr::group_by(.data$type, .data$method, .data$sample) %>%
      dplyr::summarise(errors = agg_fun(.data$errors, na.rm = TRUE)) %>%
      dplyr::ungroup()
  }
  if (!nrow(dat) > 0) {
    send_stop("No data left to plot, could you check your input?")
  }

  if (is.null(title)) {
    if (agg_mode) {
      title <- "All samples"
    } else {
      title <- unique(dat$sample)
    }
  }

  send_info("Plotting.")
  ## Plotting
  p <- plot_fun(subset(dat, dat$type != "optimal"),
    x = "method", y = "errors", color = "method", outlier.shape = outlier.shape,
    palette = palette, width = width, add = add, add.params = list(alpha = 0.3),
    title = title, xlab = xlab, ylab = ylab, legend = legend, ...
  )

  if (!agg_mode) {
    if (highlight == "auto") {
      p <- p + ggplot2::geom_point(
        data = subset(dat, dat$type == "optimal"),
        mapping = ggplot2::aes_string(x = "method", y = "errors", color = "method"),
        shape = 17, size = highlight_size,
        position = ggplot2::position_dodge2(width = dodge_width, preserve = "single")
      )
    } else {
      p <- p + ggplot2::geom_point(
        data = subset(dat, dat$type == "optimal"),
        mapping = ggplot2::aes_string(x = "method", y = "errors"),
        shape = 17, size = highlight_size, color = highlight,
        position = ggplot2::position_dodge2(width = dodge_width, preserve = "single")
      )
    }
  } else {
    p$summary <- dat
  }
  p
}


#' @rdname show_sig_bootstrap
#' @export
show_sig_bootstrap_stability <- function(bt_result, signatures = NULL,
                                         measure = c("RMSE", "CV", "MAE", "AbsDiff"),
                                         methods = "QP", plot_fun = c("boxplot", "violin"),
                                         palette = "aaas", title = NULL,
                                         xlab = FALSE, ylab = "Signature instability",
                                         width = 0.3, outlier.shape = NA,
                                         add = "jitter", add.params = list(alpha = 0.3),
                                         ...) {
  stopifnot(is.list(bt_result))
  plot_fun <- match.arg(plot_fun)
  plot_fun <- switch(
    plot_fun,
    boxplot = ggpubr::ggboxplot,
    violin = ggpubr::ggviolin
  )
  measure <- match.arg(measure)

  timer <- Sys.time()
  send_info("Started.")
  on.exit(send_elapsed_time(timer))

  dat <- bt_result$expo
  rm(bt_result)
  dat <- dplyr::filter(dat, .data$method %in% methods)

  if (!is.null(signatures)) {
    dat <- dplyr::filter(dat, .data$sig %in% signatures)
  }
  if (!nrow(dat) > 0) {
    send_stop("No data left to plot, could you check your input?")
  }

  if (measure == "AbsDiff") {
    dat_opt <- dplyr::filter(dat, .data$type == "optimal")
    dat_bt <- dplyr::filter(dat, .data$type != "optimal")

    dat_bt <- dat_bt %>%
      dplyr::group_by(.data$method, .data$sample, .data$sig) %>%
      dplyr::summarise(exposure = mean(.data$exposure)) %>%
      dplyr::ungroup()
    dat_bt$type <- "bootstrap"

    dat <- dplyr::bind_rows(dat_opt, dat_bt)

    dat <- dat %>%
      tidyr::pivot_wider(names_from = "type", values_from = "exposure") %>%
      dplyr::group_by(.data$method, .data$sig, .data$sample) %>%
      dplyr::summarise(measure = abs(.data$optimal - .data$bootstrap)) %>%
      dplyr::ungroup()
  } else {
    ## Calculate RMSE（Root Mean Squared Error), CV or MAE (Mean Absolute Error)
    if (measure %in% c("RMSE", "CV")) {
      ## across solution: https://github.com/tidyverse/dplyr/issues/5230
      dat <- dat %>%
        tidyr::pivot_wider(names_from = "type", values_from = "exposure") %>%
        dplyr::mutate(mean_rep = dplyr::select(., dplyr::starts_with("Rep_")) %>% rowMeans()) %>%
        dplyr::mutate_at(dplyr::vars(dplyr::starts_with("Rep_", ignore.case = FALSE)), ~ (. - .data$optimal)^2) %>%
        dplyr::mutate(measure = dplyr::select(., -c("method", "sample", "optimal", "sig")) %>% rowMeans() %>% sqrt())


      if (measure == "CV") {
        # Calculate Coefficient of Variation (CV) with RMSE
        dat <- dat %>%
          dplyr::mutate(measure = measure / mean_rep)
      }

      dat <- dat %>%
        dplyr::select(c("method", "sample", "sig", "measure"))
    } else {
      ## MAE
      dat <- dat %>%
        tidyr::pivot_wider(names_from = "type", values_from = "exposure") %>%
        dplyr::mutate_at(dplyr::vars(dplyr::starts_with("Rep_")), ~ abs(. - .data$optimal)) %>%
        dplyr::mutate(measure = dplyr::select(., -c("method", "sample", "optimal", "sig")) %>% rowMeans()) %>%
        dplyr::select(c("method", "sample", "sig", "measure"))
    }
  }

  send_info("Plotting.")
  ## Plotting
  plot_fun(dat,
    x = "sig", y = "measure", color = "method", outlier.shape = outlier.shape,
    palette = palette, width = width, add = add, add.params = add.params,
    title = title, xlab = xlab, ylab = ylab, ...
  )
}
