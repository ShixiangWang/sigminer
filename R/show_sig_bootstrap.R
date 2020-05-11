show_sig_bootstrap <- function(bt_result, sample = NULL, signatures = NULL, type = c("exposure", "error"),
                               palette = "aaas", xlab = FALSE, ylab = "Signature exposure", width = 0.3,
                               add = "jitter",
                               ...) {
  stopifnot(is.list(bt_result))
  type = match.arg(type)

  timer <- Sys.time()
  send_info("Started.")
  on.exit(send_elapsed_time(timer))

  if (type == "exposure") {
    dat = bt_result$expo
    rm(bt_result)
    if (!is.null(sample)) {
      dat = subset(dat, subset = dat$sample %in% sample)
    } else {
      send_warning("Top 1 sample is selected due to no sample specified in argument {.code sample}.")
    }
    if (!is.null(signatures)) {
      dat = subset(dat, subset = dat$sig %in% signatures)
    }
    if (!nrow(dat) > 0) {
      send_stop("No data left to plot, could you check your input?")
    }

    ## Plotting
    dat[, `:=`(method = factor(method),  sig = factor(sig))]
    ggpubr::ggboxplot(subset(dat, dat$type != "optimal"), x = "sig", y = "exposure", color = "method",
                      palette = palette, xlab = xlab, ylab = ylab, width = width, add = add, ...)
    + ggplot2::geom_point(data = subset(dat, dat$type == "optimal"),
                          mapping = aes_string(x = "sig", y = "exposure", color = "method"),
                          shape = 2, position = ggplot2::position_dodge2(width = 0.8))

  } else {
    dat = bt_result$error
    rm(bt_result)
    if (!is.null(sample)) {
      dat = subset(dat, subset = dat$sample %in% sample)
    } else {
      send_warning("Top 1 sample is selected due to no sample specified in argument {.code sample}.")
    }
    if (!is.null(signatures)) {
      send_warning("{.code type} is 'error', so the input {.code signatures} are ignored.")
    }
    if (!nrow(dat) > 0) {
      send_stop("No data left to plot, could you check your input?")
    }
  }
}
