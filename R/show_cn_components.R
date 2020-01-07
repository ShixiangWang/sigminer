#' Show Copy Number Components
#'
#' Show mixture fit model components ("Macintyre" ("M") method) or
#' standard classified components ("Wang" ("W") method) for copy number data.
#'
#' @inheritParams show_cn_features
#' @param parameters a `data.frame` contain parameter components, obtain this
#' from [sig_feed] function.
#' @param auto_transform default is `TRUE`, it will auto increase the `SD` for components
#' for showing them better in the plot. Only used when method is "Macintyre".
#' @param show_weights default is `TRUE`, show weights for each component.
#' Only used when method is "Macintyre".
#' @param ... other options pass to \code{\link[cowplot]{plot_grid}} function of **cowplot** package.
#' @inheritParams show_cn_distribution
#' @inheritParams sig_feed
#' @author Shixiang Wang <w_shixiang@163.com>
#' @return a `ggplot` object
#' @export
#' @examples
#' # Load copy number prepare object
#' load(system.file("extdata", "toy_copynumber_prepare.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' show_cn_components(cn_prepare$parameters)
#' show_cn_components(cn_prepare$parameters, show_weights = FALSE)
show_cn_components <- function(parameters,
                               method = "Macintyre",
                               show_weights = TRUE,
                               log_segsize = TRUE,
                               auto_transform = TRUE,
                               return_plotlist = FALSE,
                               base_size = 12, ...) {
  stopifnot(is.logical(return_plotlist))

  method <- match.arg(method, choices = c("Macintyre", "M", "Wang", "W"))

  if (startsWith(method, "M")) {
    dat <- parameters %>%
      dplyr::group_by(.data$feature) %>%
      dplyr::mutate(weights = .data$n_obs / sum(.data$n_obs)) %>%
      dplyr::ungroup()

    ft_order <- c("segsize", "copynumber", "changepoint", "bp10MB", "bpchrarm", "osCN")

    if (log_segsize) {
      dat$mean[dat$feature == "segsize"] <- log10(dat$mean[dat$feature == "segsize"])
      dat$sd[dat$feature == "segsize"] <- log10(dat$sd[dat$feature == "segsize"])
      ft_annot <- c(
        "Segment size (log10 based)", "Copy number", "Copy number changepoint", "Breakpoint count per 10MB",
        "Breakpoint count per chr arm", "Oscilating CN chain length"
      )
    } else {
      ft_annot <- c(
        "Segment size", "Copy number", "Copy number changepoint", "Breakpoint count per 10MB",
        "Breakpoint count per chr arm", "Oscilating CN chain length"
      )
    }

    if (auto_transform) {
      # Make components with very small SD and norm distribution
      # have two orders of magnitude
      # smaller than Mean
      .transform <- function(mean, sd) {
        tms <- mean / sd
        ifelse(
          tms > 1000,
          sd * 10^ceiling(log10(tms) - 2),
          sd
        )
      }

      dat <- dat %>%
        dplyr::mutate(
          sd = ifelse(.data$dist == "norm",
            .transform(.data$mean, .data$sd),
            .data$sd
          )
        )
    }

    annot_df <- dplyr::tibble(
      feature = ft_order,
      annotation = ft_annot
    )

    if (!show_weights) {
      dat$weights <- 1
    }

    requireNamespace("cowplot")

    cbPalette <-
      c(
        RColorBrewer::brewer.pal(8, "Dark2"),
        RColorBrewer::brewer.pal(9, "Set1"),
        "black"
      )

    plot_mix_comp_norm <- function(x, mean, sd, weight) {
      weight * stats::dnorm(x, mean, sd)
    }
    plot_mix_comp_pois <- function(x, mean, weight) {
      weight * stats::dpois(x, mean)
    }


    plot_fun <- function(df, x_label = NULL, base_size = 12) {
      row_index <- seq_len(nrow(df))
      set.seed(1234)
      if (df$dist[1] == "pois") {
        # Plot pois dist
        # Generate 1000 random data to get range of x-axis
        rdata <- purrr::map(df$mean, function(x) stats::rpois(1000, lambda = x)) %>%
          Reduce(c, .) %>%
          range(na.rm = TRUE)
        p <- ggplot(data.frame(x = rdata), aes(.data$x)) + ylab(NULL)
        for (i in row_index) {
          p <- p + stat_function(
            geom = "line",
            n = diff(rdata) + 1,
            fun = plot_mix_comp_pois,
            args = list(df$mean[i], df$weights[i]),
            color = cbPalette[i]
          )
        }
        p <- p + xlab(x_label) + cowplot::theme_cowplot(font_size = base_size)
        p
      } else {
        # Plot normal dist
        # Similar to pois dist
        rdata <- purrr::map2(df$mean, df$sd, function(x, y) stats::rnorm(1000, mean = x, sd = y)) %>%
          Reduce(c, .) %>%
          range(na.rm = TRUE)
        p <- ggplot(data.frame(x = rdata), aes(.data$x)) + ylab(NULL)
        for (i in row_index) {
          p <- p + stat_function(
            geom = "line",
            n = 1000,
            fun = plot_mix_comp_norm,
            args = list(df$mean[i], df$sd[i], df$weights[i]),
            color = cbPalette[i]
          )
        }
        p <- p + xlab(x_label) + cowplot::theme_cowplot(font_size = base_size)
        p
      }
    }

    zz <- dat %>%
      dplyr::group_by(.data$feature) %>%
      tidyr::nest() %>%
      dplyr::left_join(annot_df, by = "feature") %>%
      dplyr::mutate(
        density_plot = purrr::map2(.data$data, .data$annotation, plot_fun, base_size = base_size)
      )
    plotlist <- zz %>%
      dplyr::pull(.data$density_plot)
    names(plotlist) <- zz$feature
    rm(zz)
    # Keep same order as feature distribution
    plotlist <- plotlist[ft_order]

    if (return_plotlist) {
      return(plotlist)
    } else {
      p <- cowplot::plot_grid(
        plotlist = plotlist,
        nrow = 2,
        align = "hv",
        ...
      )
      p
    }
  } else {
    plot_fun2 <- function(df, f_name, base_size) {
      df$component <- factor(df$component, levels = df$component)
      ggplot(data = df, aes_string(x = "component", y = "n_obs")) +
        geom_bar(stat = "identity", color = "white") +
        labs(x = NULL, y = NULL) +
        theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm")) +
        cowplot::theme_cowplot(font_size = base_size) +
        theme(
          axis.text.x = element_text(
            angle = 60, vjust = 1,
            hjust = 1,
            color = "black"
          )
        )
    }

    p_df <- parameters %>%
      dplyr::as_tibble() %>%
      dplyr::group_by(.data$feature) %>%
      tidyr::nest() %>%
      dplyr::mutate(p_list = purrr::map2(data, .data$feature, plot_fun2, base_size = base_size))
    p_list <- p_df$p_list
    names(p_list) <- p_df$feature

    if (return_plotlist) {
      return(p_list)
    } else {
      if (length(p_list) != 7) {
        p <- cowplot::plot_grid(plotlist = p_list, nrow = 2, align = "hv", ...)
      } else {
        top_row <- cowplot::plot_grid(plotlist = p_list[1:4], nrow = 1, align = "h", ...)
        bot_row <- cowplot::plot_grid(plotlist = p_list[-(1:4)], nrow = 1, align = "h", rel_widths = c(1, 1, 2), ...)
        p <- cowplot::plot_grid(top_row, bot_row, nrow = 2)
      }
      p
    }
  }
}
