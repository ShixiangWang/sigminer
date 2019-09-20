#' Show Mixture Fit Model Components
#'
#' Show mixture fit model components for copy number data.
#'
#' @inheritParams show_cn_features
#' @param components a `list` contain flexmix object of copy-number features, obtain this
#' from [derive] function.
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
#' show_cn_components(cn_prepare$features, cn_prepare$components)
show_cn_components <- function(features, components, return_plotlist = FALSE,
                               base_size = 12, ...) {

  # //TODO: rewrite this function to make sure all components shown with their weights
  stopifnot(is.logical(return_plotlist))
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

    p <- p + xlab(xlab) + cowplot::theme_cowplot(font_size = base_size)
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

    p <- p + xlab(xlab) + cowplot::theme_cowplot(font_size = base_size)
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
    mu = c(
      10^comp_segsize[1, ], comp_copynumber[1, ],
      comp_changepoint[1, ],
      comp_bp10MB, comp_bpchrarm,
      comp_osCN
    ),
    sd = c(
      10^comp_segsize[2, ], comp_copynumber[2, ],
      comp_changepoint[2, ],
      rep(NA, length(comp_bp10MB)), rep(NA, length(comp_bpchrarm)),
      rep(NA, length(comp_osCN))
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

  if (return_plotlist) {
    return(list(
      segsize = p_1, copynumber = p_2, changepoint = p_3,
      bp10MB = p_4, bpchrarm = p_5, osCN = p_6, parameters = parameters
    ))
  }

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
