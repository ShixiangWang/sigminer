#' Show Signature Fit Result
#'
#' See [sig_fit] for examples.
#'
#' @inheritParams ggpubr::ggboxplot
#' @inheritParams ggpubr::ggpar
#' @inheritParams ggplot2::geom_boxplot
#' @inheritParams sig_fit_bootstrap_batch
#' @inheritParams show_sig_bootstrap
#' @param fit_result result object from [sig_fit].
#' @param samples samples to show, if `NULL`, all samples are used.
#' @return a `ggplot` object.
#' @export
#' @seealso [sig_fit], [show_sig_bootstrap_exposure], [sig_fit_bootstrap], [sig_fit_bootstrap_batch]
show_sig_fit <- function(fit_result, samples = NULL, signatures = NULL,
                         plot_fun = c("boxplot", "violin", "scatter"),
                         palette = "aaas",
                         title = NULL,
                         xlab = FALSE, ylab = "Signature exposure", legend = "none",
                         width = 0.3, outlier.shape = NA,
                         add = "jitter", add.params = list(alpha = 0.3),
                         ...) {
  fun_setting <- plot_fun <- match.arg(plot_fun)
  plot_fun <- switch(plot_fun,
    boxplot = ggpubr::ggboxplot,
    violin = ggpubr::ggviolin,
    scatter = ggpubr::ggscatter
  )

  timer <- Sys.time()
  send_info("Started.")
  on.exit(send_elapsed_time(timer))

  send_info("Checking input format.")
  if (data.table::is.data.table(fit_result)) {
    dat <- data.table::copy(fit_result)
  } else if (is.list(fit_result)) {
    dat <- fit_result$expo
  } else if (is.matrix(fit_result)) {
    dat <- fit_result %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Sig") %>%
      tidyr::pivot_longer(cols = -"Sig", names_to = "sample", values_to = "expo") %>%
      tidyr::pivot_wider(id_cols = "sample", names_from = "Sig", values_from = "expo") %>%
      data.table::as.data.table()
  } else {
    send_stop("Bad input, could you take a check?")
  }

  send_success("Checked.")

  send_info("Checking filters.")
  if (!is.null(samples)) {
    dat <- dplyr::filter(dat, .data$sample %in% samples)
  }
  if (!is.null(signatures)) {
    dat <- dplyr::select(dat, c("sample", signatures))
  }
  if (!nrow(dat) > 0) {
    send_stop("No data left to plot, could you check your input?")
  }
  send_info("Checked.")

  dat <- dat %>%
    dplyr::as_tibble() %>%
    tidyr::pivot_longer(-"sample", names_to = "sig", values_to = "exposure")

  send_info("Plotting.")
  ## Plotting
  if (isFALSE(fun_setting == "scatter")) {
    plot_fun(dat,
      x = "sig", y = "exposure", color = "sig", outlier.shape = outlier.shape,
      palette = palette, width = width, add = add, add.params = add.params,
      title = title, xlab = xlab, ylab = ylab, legend = legend, ...
    )
  } else {
    if (legend == "none") {
      send_warning("When plot_fun='scatter', setting legend='top' is recommended.")
    }
    plot_fun(dat,
      x = "sig", y = "exposure", color = "sample", shape = "sample",
      palette = palette,
      title = title, xlab = xlab, ylab = ylab, legend = legend, ...
    )
  }
}
