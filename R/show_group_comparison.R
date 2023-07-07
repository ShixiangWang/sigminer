#' Plot Group Comparison Result
#'
#' Using result data from [get_group_comparison], this function plots
#' genotypes/phenotypes comparison between signature groups using **ggplot2** package and return
#' a list of `ggplot` object contains individual and combined plots. The combined
#' plot is easily saved to local using [cowplot::save_plot()]. Of note, default fisher
#' test p values are shown for categorical data and fdr values are shown for
#' continuous data.
#' @inheritParams get_adj_p
#' @inheritParams show_cn_distribution
#' @param group_comparison a `list` from result of [get_group_comparison] function.
#' @param xlab lab name of x axis for all plots. if it is `NA`, remove title for x axis.
#' @param ylab_co lab name of y axis for plots of continuous type data. Of note,
#' this argument should be a character vector has same length as `group_comparison`,
#' the location for categorical type data should mark with `NA`.
#' @param legend_title_ca legend title for plots of categorical type data.
#' @param legend_position_ca legend position for plots of categorical type data.
#' Of note,
#' this argument should be a character vector has same length as `group_comparison`,
#' the location for continuous type data should mark with `NA`.
#' @param set_ca_sig_yaxis if `TRUE`, use y axis to show signature proportion instead of
#' variable proportion.
#' @param set_ca_custom_xlab only works when `set_ca_sig_yaxis` is `TRUE`. If
#' `TRUE`, set x labels using input `xlab`, otherwise variable names will be used.
#' @param show_pvalue if `TRUE`, show p values.
#' @param ca_p_threshold a p threshold for categorical variables, default is 0.01.
#' A p value less than 0.01 will be shown as `P < 0.01`.
#' @param font_size_x font size for x.
#' @param text_angle_x text angle for x.
#' @param text_hjust_x adjust x axis text
#' @param ... other paramters pass to [ggpubr::compare_means()] or [ggpubr::stat_compare_means()]
#' according to the specified `method`.
#' @author Shixiang Wang <w_shixiang@163.com>
#' @return a `list` of `ggplot` objects.
#' @import ggplot2
#' @export
#' @examples
#' \donttest{
#' load(system.file("extdata", "toy_copynumber_signature_by_W.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#'
#' # Assign samples to clusters
#' groups <- get_groups(sig, method = "k-means")
#'
#' set.seed(1234)
#'
#' groups$prob <- rnorm(10)
#' groups$new_group <- sample(c("1", "2", "3", "4", NA), size = nrow(groups), replace = TRUE)
#'
#' # Compare groups (filter NAs for categorical coloumns)
#' groups.cmp <- get_group_comparison(groups[, -1],
#'   col_group = "group",
#'   cols_to_compare = c("prob", "new_group"),
#'   type = c("co", "ca"), verbose = TRUE
#' )
#'
#' # Compare groups (Set NAs of categorical columns to 'Rest')
#' groups.cmp2 <- get_group_comparison(groups[, -1],
#'   col_group = "group",
#'   cols_to_compare = c("prob", "new_group"),
#'   type = c("co", "ca"), NAs = "Rest", verbose = TRUE
#' )
#'
#' show_group_comparison(groups.cmp)
#'
#' ggcomp <- show_group_comparison(groups.cmp2)
#' ggcomp$co_comb
#' ggcomp$ca_comb
#' }
#' @testexamples
#' expect_equal(length(groups.cmp) + 1L, length(groups.cmp) + 1)
#' expect_s3_class(ggcomp$ca$new_group, "ggplot")
show_group_comparison <- function(group_comparison,
                                  xlab = "group", ylab_co = NA,
                                  legend_title_ca = NA,
                                  legend_position_ca = "bottom",
                                  set_ca_sig_yaxis = FALSE,
                                  set_ca_custom_xlab = FALSE,
                                  show_pvalue = TRUE,
                                  ca_p_threshold = 0.01,
                                  method = "wilcox.test",
                                  p.adjust.method = "fdr",
                                  base_size = 12,
                                  font_size_x = 12,
                                  text_angle_x = 30,
                                  text_hjust_x = 0.2,
                                  ...) {
  # parameter with ca/co in the end need fill values

  # add parameters
  group_comparison <- Map(function(x, xlab, legend_title_ca, ylab_co, legend_position_ca) {
    x[["xlab"]] <- xlab
    if (x[["type"]] == "categorical") {
      x[["legend_title"]] <- legend_title_ca
      x[["legend_position"]] <- legend_position_ca
    }
    if (x[["type"]] == "continuous") x[["ylab"]] <- ylab_co
    x
  }, group_comparison, xlab, legend_title_ca, ylab_co, legend_position_ca)

  # split input into two list according to element is categorical or continuous
  ca_index <- which(sapply(group_comparison, function(x) x$type) == "categorical")
  co_index <- which(sapply(group_comparison, function(x) x$type) == "continuous")

  n_left <- length(group_comparison) - length(ca_index) - length(co_index)

  if (n_left > 0) {
    warning(n_left, " elements drop of because their type is not either categorical or continuous.")
  }

  if (length(ca_index) > 0) ca_list <- group_comparison[ca_index]
  if (length(co_index) > 0) co_list <- group_comparison[co_index]

  if (set_ca_custom_xlab) {
    ca_list <- lapply(ca_list, function(x) {
      x$extra <- TRUE
      return(x)
    })
  }
  if (length(ca_index) > 0) {
    # plot categorical data
    ca_res <- lapply(ca_list, function(df) {
      data <- df[["data"]] %>% dplyr::as_tibble()
      if (set_ca_sig_yaxis) {
        if (!isTRUE(df[["extra"]])) {
          df[["xlab"]] <- colnames(data)[2]
        }
        colnames(data) <- c("Dominant Signature", "group")
      }
      data_sum <- data %>% dplyr::count(.data$group)
      data_sum[["labels"]] <- paste(data_sum[["group"]], paste0("(n=", data_sum[["n"]], ")"), sep = "\n")

      var_name <- setdiff(colnames(data), "group")
      var_name2 <- ifelse(isValidAndUnreserved(var_name), var_name, paste0("`", var_name, "`"))

      p <- ggplot(data, aes_string(x = "group", fill = var_name2)) +
        geom_bar(position = "fill") +
        cowplot::theme_cowplot(font_size = base_size) +
        theme(axis.title.y = element_blank()) +
        scale_x_discrete(
          breaks = data_sum[["group"]],
          labels = data_sum[["labels"]]
        )
      if (is.na(df[["xlab"]])) {
        p <- p + theme(axis.title.x = element_blank())
      } else {
        p <- p + xlab(df[["xlab"]])
      }

      if (!is.na(df[["legend_title"]])) {
        p <- p + scale_fill_discrete(name = df[["legend_title"]])
      }

      if (show_pvalue) {
        if (!is.na(df[["p_value"]])) {
          p <- p + labs(title = ifelse(df$p_value < ca_p_threshold,
            paste("P <", signif(ca_p_threshold, 3)),
            paste("P =", signif(df$p_value, 3))
          ))
        }
      }

      p + theme(
        legend.position = df[["legend_position"]],
        axis.text.x = element_text(size = font_size_x, angle = text_angle_x, hjust = text_hjust_x)
      )
    })
    names(ca_res) <- names(ca_list)
  } else {
    ca_res <- NA
  }

  if (length(co_index) > 0) {
    # plot continuous data
    co_res <- lapply(co_list, function(df, ...) {
      data <- df[["data"]] %>% dplyr::as_tibble()
      data_sum <- data %>% dplyr::count(.data$group)
      data_sum[["labels"]] <- paste(data_sum[["group"]], paste0("(n=", data_sum[["n"]], ")"), sep = "\n")
      my_comparisons <- combn(unique(as.character(data[["group"]])),
        2,
        simplify = FALSE
      )

      var_name <- colnames(data)[2]
      var_name2 <- ifelse(isValidAndUnreserved(var_name), var_name, paste0("`", var_name, "`"))

      p <- ggplot(data, aes_string(x = "group", y = var_name2)) +
        geom_boxplot() +
        cowplot::theme_cowplot(font_size = base_size) +
        scale_x_discrete(
          breaks = data_sum[["group"]],
          labels = data_sum[["labels"]]
        )

      if (is.na(df[["xlab"]])) {
        p <- p + theme(axis.title.x = element_blank())
      } else {
        p <- p + xlab(df[["xlab"]])
      }

      if (!is.na(df[["ylab"]])) {
        p <- p + ylab(df[["ylab"]])
      }

      if (show_pvalue) {
        if (!requireNamespace("ggpubr", quietly = TRUE)) {
          stop("'ggpubr' package is needed for plotting p values.")
        }

        if ((method != "wilcox.test") & (method != "t.test")) {
          p <- p + ggpubr::stat_compare_means(method = method, ...)
        } else {
          p_df <- get_adj_p(data,
            .col = var_name2, .grp = "group",
            method = method, p.adjust.method = p.adjust.method, ...
          )

          # p <- p + ggpubr::stat_compare_means(
          #   comparisons = my_comparisons,
          #   ...
          # )
          p <- p + ggpubr::stat_pvalue_manual(p_df, label = "p.adj")
        }
      }
      p + theme(axis.text.x = element_text(size = font_size_x, angle = text_angle_x, hjust = text_hjust_x))
    }, ...)
    names(co_res) <- names(co_list)
  } else {
    co_res <- NA
  }

  if (all(!is.na(ca_res))) {
    if (length(ca_res) <= 3) {
      ca_comb <- cowplot::plot_grid(plotlist = ca_res, align = "h")
    } else if (length(ca_res) == 4) {
      ca_comb <- cowplot::plot_grid(plotlist = ca_res, align = "hv", ncol = 2)
    } else if (length(ca_res) <= 9) {
      ca_comb <- cowplot::plot_grid(plotlist = ca_res, align = "hv", ncol = 3)
    } else {
      ca_comb <- cowplot::plot_grid(plotlist = ca_res, align = "hv", ncol = 4)
    }
  } else {
    ca_comb <- NA
  }

  if (all(!is.na(co_res))) {
    if (length(co_res) <= 3) {
      co_comb <- cowplot::plot_grid(plotlist = co_res, align = "h")
    } else if (length(co_res) == 4) {
      co_comb <- cowplot::plot_grid(
        plotlist = co_res, align = "h",
        axis = "l", ncol = 2
      )
    } else if (length(co_res) <= 9) {
      co_comb <- cowplot::plot_grid(
        plotlist = co_res, align = "hv",
        axis = "l", ncol = 3
      )
    } else {
      co_comb <- cowplot::plot_grid(
        plotlist = co_res, align = "hv",
        axis = "l", ncol = 4
      )
    }
  } else {
    co_comb <- NA
  }

  list(
    ca = ca_res,
    ca_comb = ca_comb,
    co = co_res,
    co_comb = co_comb
  )
}
