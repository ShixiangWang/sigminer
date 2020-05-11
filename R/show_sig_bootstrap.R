show_sig_bootstrap_exposure <- function(bt_result, sample = NULL, signatures = NULL,
                                        method = "LS", palette = "aaas", title = NULL,
                                        xlab = FALSE, ylab = "Signature exposure", width = 0.3,
                                        dodge_width = 0.8,
                                        add = "jitter",
                                        ...) {
  stopifnot(is.list(bt_result))

  timer <- Sys.time()
  send_info("Started.")
  on.exit(send_elapsed_time(timer))

  dat = bt_result$expo
  rm(bt_result)
  dat = dplyr::filter(dat, .data$method %in% .env$method)
  if (!is.null(sample)) {
    if (length(sample) > 1) {
      send_stop("Only one sample can be used!")
    }
    dat = dplyr::filter(dat, .data$sample %in% .env$sample)
  } else {
    send_warning("Top 1 sample is selected due to no sample specified in argument {.code sample}.")
    dat = dplyr::filter(dat, .data$sample %in% .data$sample[1])
  }
  if (!is.null(signatures)) {
    dat = dplyr::filter(dat, .data$sig %in% signatures)
  }
  if (!nrow(dat) > 0) {
    send_stop("No data left to plot, could you check your input?")
  }

  if (is.null(title)) {
    title <- unique(dat$sample)
  }

  send_info("Plotting.")
  ## Plotting
  ggpubr::ggboxplot(subset(dat, dat$type != "optimal"), x = "sig", y = "exposure", color = "method",  outlier.shape = NA,
                    palette = palette, xlab = xlab, ylab = ylab, width = width, add = add, title = title, ...) +
    ggplot2::geom_point(data = subset(dat, dat$type == "optimal"),
                        mapping = aes_string(x = "sig", y = "exposure", color = "method"),
                        shape = 2, position = ggplot2::position_dodge2(width = dodge_width, preserve = "single"))

}


show_sig_bootstrap_error <- function(bt_result, sample = NULL,
                                     method = "LS", palette = "aaas", title = NULL,
                                     xlab = FALSE, ylab = "Decomposition error", width = 0.3,
                                     dodge_width = 0.8,
                                     add = "jitter",
                                     legend = "none",
                                     ...) {
  stopifnot(is.list(bt_result))

  timer <- Sys.time()
  send_info("Started.")
  on.exit(send_elapsed_time(timer))

  dat = bt_result$error
  rm(bt_result)
  dat = dplyr::filter(dat, .data$method %in% .env$method)
  if (!is.null(sample)) {
    if (length(sample) > 1) {
      send_stop("Only one sample can be used!")
    }
    dat = dplyr::filter(dat, .data$sample %in% .env$sample)
  } else {
    send_warning("Top 1 sample is selected due to no sample specified in argument {.code sample}.")
    dat = dplyr::filter(dat, .data$sample %in% .data$sample[1])
  }
  if (!nrow(dat) > 0) {
    send_stop("No data left to plot, could you check your input?")
  }

  if (is.null(title)) {
    title <- unique(dat$sample)
  }

  send_info("Plotting.")
  ## Plotting
  ggpubr::ggboxplot(subset(dat, dat$type != "optimal"), x = "method", y = "errors", color = "method",  outlier.shape = NA,
                    palette = palette, xlab = xlab, ylab = ylab, width = width, add = add, title = title,  legend = legend, ...) +
    ggplot2::geom_point(data = subset(dat, dat$type == "optimal"),
                        mapping = aes_string(x = "method", y = "errors", color = "method"),
                        shape = 2, position = ggplot2::position_dodge2(width = dodge_width, preserve = "single"))
}


show_sig_bootstrap_stability <- function(bt_result, signatures = NULL,
                                         method = "LS", palette = "aaas", title = NULL,
                                         xlab = FALSE, ylab = "Signature instability", width = 0.3,
                                         add = "jitter",
                                         ...) {
  stopifnot(is.list(bt_result))

  timer <- Sys.time()
  send_info("Started.")
  on.exit(send_elapsed_time(timer))

  dat = bt_result$expo
  rm(bt_result)
  dat = dplyr::filter(dat, .data$method %in% .env$method)

  if (!is.null(signatures)) {
    dat = dplyr::filter(dat, .data$sig %in% signatures)
  }
  if (!nrow(dat) > 0) {
    send_stop("No data left to plot, could you check your input?")
  }

  dat_opt = dplyr::filter(dat, .data$type == "optimal")
  dat_bt  = dplyr::filter(dat, .data$type != "optimal")

  dat_bt = dat_bt %>%
    dplyr::group_by(.data$method, .data$sample, .data$sig) %>%
    dplyr::summarise(exposure = mean(.data$exposure)) %>%
    dplyr::ungroup()
  dat_bt$type = "bootstrap"

  dat = dplyr::bind_rows(dat_opt, dat_bt)

  ## Calculate MSE（Mean Squared Error）
  dat = dat %>%
    tidyr::pivot_wider(names_from = "type", values_from = "exposure") %>%
    dplyr::group_by(.data$method, .data$sig, .data$sample) %>%
    dplyr::summarise(mse = mean((.data$optimal - .data$bootstrap)^2)) %>%
    dplyr::ungroup()

  send_info("Plotting.")
  ## Plotting
  ggpubr::ggboxplot(dat, x = "sig", y = "mse", color = "method", outlier.shape = NA,
                    palette = palette, xlab = xlab, ylab = ylab, width = width, add = add, title = title, ...)

}
