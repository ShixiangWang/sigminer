#' Show Mixture Fit Model Components
#'
#' Show mixture fit model components for copy number data.
#'
#' @inheritParams show_cn_features
#' @param parameters a `data.frame` contain parameter components, obtain this
#' from [derive] function.
#' @param show_weights default is `TRUE`, show weights for each component.
#' @param ... other options pass to \code{\link[cowplot]{plot_grid}} function of **cowplot** package.
#' @inheritParams show_cn_distribution
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
                               show_weights = TRUE,
                               return_plotlist = FALSE,
                               base_size = 12, ...) {
  stopifnot(is.logical(return_plotlist))
  dat <- parameters %>%
    dplyr::group_by(.data$feature) %>%
    dplyr::mutate(weights = .data$n_obs / sum(.data$n_obs)) %>%
    dplyr::ungroup()

  ft_order <- c("segsize", "copynumber", "changepoint", "bp10MB", "bpchrarm", "osCN")
  ft_annot <- c(
    "Segment size", "Copy number", "Copy number changepoint", "Breakpoint count per 10MB",
    "Breakpoint count per chr arm", "Oscilating CN chain length"
  )
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
}
