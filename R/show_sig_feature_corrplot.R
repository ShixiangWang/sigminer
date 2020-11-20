#' Draw Corrplot for Signature Exposures and Other Features
#'
#' This function is for association visualization. Of note,
#' the parameters `p_val` and `drop` will affect the visualization
#' of association results under p value threshold.
#'
#' @param tidy_cor data returned by [get_tidy_association].
#' @param feature_list a character vector contains features want to be plotted.
#' If missing, all features will be used.
#' @param sort_features default is `FALSE`, use feature order obtained from the previous
#' step. If `TRUE`, sort features as `feature_list`.
#' @param sig_orders signature levels for ordering.
#' @param return_plotlist if `TRUE`, return as a list of `ggplot` objects.
#' @param p_val p value threshold. If p value larger than this threshold,
#' the result becomes blank white.
#' @param drop if `TRUE`, when a feature has no association with all signatures
#' (p value larger than threshold set by `p_val`), this feature will be removed
#' from the plot. Otherwise, this feature (a row) will keep with all blank white.
#' @param xlab label for x axis.
#' @param ylab label for y axis.
#' @param co_gradient_colors a Scale object representing gradient colors used to plot for continuous features.
#' @param ca_gradient_colors a Scale object representing gradient colors used to plot for categorical features.
#' @param plot_ratio a length-2 numeric vector to set the height/width ratio.
#' @param breaks_count breaks for sample count. If set it to `NULL`,
#' ggplot `bin` scale will be used to automatically determine the
#' breaks. If set it to `NA`, `aes` for sample will be not used.
#' @return a `ggplot2` object
#' @export
#' @examples
#'
#' # The data is generated from Wang, Shixiang et al.
#' load(system.file("extdata", "asso_data.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#'
#' p <- show_sig_feature_corrplot(tidy_data.seqz.feature, p_val = 0.05)
#' p
#' @testexamples
#' expect_s3_class(p, "ggplot")
#' @seealso [get_tidy_association] and [get_sig_feature_association]
show_sig_feature_corrplot <- function(tidy_cor, feature_list,
                                      sort_features = FALSE,
                                      sig_orders = NULL,
                                      drop = TRUE,
                                      return_plotlist = FALSE,
                                      p_val = 0.05,
                                      xlab = "Signatures",
                                      ylab = "Features",
                                      co_gradient_colors = scale_color_gradient2(
                                        low = "blue",
                                        mid = "white", high = "red", midpoint = 0
                                      ),
                                      ca_gradient_colors = co_gradient_colors,
                                      plot_ratio = "auto",
                                      breaks_count = c(
                                        0L,
                                        200L, 400L, 600L, 800L, 1020L
                                      )) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    message("'patchwork' package is required to use this feature.")
    return(invisible(NULL))
  }

  if (!is.character(plot_ratio)) {
    if (length(plot_ratio) != 2 | !is.numeric(plot_ratio)) {
      stop("plot_ratio must be a length-2 numeric vector!")
    }
  }
  if (missing(feature_list)) {
    feature_list <- tidy_cor$feature %>% unique()
  }

  if (is.null(breaks_count)) {
    data <- tidy_cor %>%
      dplyr::mutate(
        Samples = .data$count,
        signature = factor(
          .data$signature,
          levels = if (is.null(sig_orders)) unique(.data$signature) else sig_orders
        )
      ) %>%
      dplyr::filter(.data$feature %in% feature_list)

    size_limits <- range(data$Samples)
  } else if (any(is.na(breaks_count))) {
    data <- tidy_cor %>%
      dplyr::mutate(
        signature = factor(
          .data$signature,
          levels = if (is.null(sig_orders)) unique(.data$signature) else sig_orders
        )
      ) %>%
      dplyr::filter(.data$feature %in% feature_list)
  } else {
    data <- tidy_cor %>%
      dplyr::mutate(
        Samples = cut(.data$count,
          breaks = breaks_count
        ),
        signature = factor(
          .data$signature,
          levels = if (is.null(sig_orders)) unique(.data$signature) else sig_orders
        )
      ) %>%
      dplyr::filter(.data$feature %in% feature_list)
  }


  if (drop) {
    data <- data %>%
      dplyr::filter(.data$p <= p_val)
  } else {
    # Fill measure with 0, so the feature is kept with blank
    data <- data %>%
      dplyr::mutate(
        measure = ifelse(.data$p > p_val | is.na(.data$p), 0, .data$measure),
      )

    if (is.numeric(breaks_count)) {
      size_levels <- levels(data$Samples)
      data <- data %>%
        dplyr::mutate(
          Samples = ifelse(is.na(.data$Samples),
            size_levels[1], .data$Samples %>% as.character()
          ),
          Samples = factor(.data$Samples, levels = size_levels)
        )
    }
  }

  .plot_cor <- function(data, type = "co") {
    if (sort_features) {
      p <- ggplot2::ggplot(data, ggplot2::aes(
        x = signature,
        y = factor(.data$feature, levels = rev(feature_list))
      ))
    }
    else {
      p <- ggplot2::ggplot(data, ggplot2::aes(
        x = signature,
        y = .data$feature
      ))
    }

    if (any(is.na(breaks_count))) {
      p <- p + ggplot2::geom_point(ggplot2::aes_string(
        color = "measure"
      ), size = 5)
    } else {
      p <- p + ggplot2::geom_point(ggplot2::aes_string(
        color = "measure",
        size = "Samples"
      ))

      if (is.null(breaks_count)) {
        p <- p +
          ggplot2::scale_size_binned(
            limits = size_limits,
            guide = ggplot2::guide_bins(show.limits = TRUE)
          )
      } else {
        p <- p + ggplot2::scale_size_discrete(drop = FALSE)
      }
    }

    if (type == "co") {
      p <- p + co_gradient_colors
    } else {
      p <- p + ca_gradient_colors
    }

    p +
      ggplot2::scale_x_discrete(drop = FALSE)
  }
  data2 <- data
  data <- data %>%
    dplyr::group_by(.data$type) %>%
    tidyr::nest() %>%
    dplyr::mutate(gg = purrr::map2(.data$data, .data$type, .plot_cor)) %>%
    dplyr::ungroup()
  gglist <- data$gg
  names(gglist) <- data$type
  gglist <- purrr::map2(gglist, names(gglist), function(p,
                                                        type) {
    p <- p + cowplot::theme_cowplot() +
      ggplot2::theme(axis.text.x = element_text(
        angle = 90, vjust = 0.5, hjust = 1
      )) +
      ggplot2::labs(x = xlab, y = ylab)
    if (type == "co") {
      p <- p + labs(color = "Correlation\ncoefficient")
    }
    else {
      p <- p + labs(color = "Difference\nin means of\nexposure")
    }
  })
  if (length(gglist) > 1) {
    if (return_plotlist) {
      return(gglist)
    }
    if (is.character(plot_ratio)) {
      count_df <- data2 %>%
        dplyr::mutate(type = factor(.data$type,
          levels = c("co", "ca")
        )) %>%
        dplyr::group_by(.data$type) %>%
        dplyr::summarise(count = length(unique(.data$feature)))
      heights <- count_df$count
      widths <- NULL
      co <- gglist$co + labs(x = NULL) + theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
      ca <- gglist$ca + guides(size = FALSE) +
        labs(y = NULL)
      co + ca + patchwork::plot_layout(
        byrow = TRUE, heights = heights,
        widths = widths
      )
    }
    else {
      heights <- plot_ratio
      widths <- NULL
      co <- gglist$co + labs(x = NULL) + theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
      ca <- gglist$ca + ggplot2::guides(size = FALSE) +
        labs(y = NULL)
      co + ca + patchwork::plot_layout(
        byrow = TRUE, heights = heights,
        widths = widths
      )
    }
  }
  else {
    return(gglist[[1]])
  }
}
