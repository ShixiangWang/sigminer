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
#' @param show_pvalue if `TRUE`, show p values.
#' @param font_size_x font size for x.
#' @param text_angle_x text angle for x.
#' @param hjust adjust x axis text
#' @param ... other paramters pass to [ggpubr::compare_means()] or [ggpubr::stat_compare_means()]
#' according to the specified `method`.
#' @author Shixiang Wang <w_shixiang@163.com>
#' @return a `list` of `ggplot` objects.
#' @import ggplot2
#' @export
#' @examples
#' \donttest{
#' # Load copy number prepare object
#' load(system.file("extdata", "toy_copynumber_prepare.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Extract copy number signatures
#' library(NMF)
#' sig <- sig_extract(cn_prepare$nmf_matrix, 2, nrun = 10)
#'
#' # Assign samples to clusters
#' groups <- get_groups(sig, method = "samples")
#'
#' set.seed(1234)
#' # Add custom groups
#' groups$new_group <- sample(c("1", "2", "3", "4", NA), size = nrow(groups), replace = TRUE)
#'
#' # Compare groups (filter NAs for categorical coloumns)
#' groups.cmp <- get_group_comparison(groups[, -1],
#'   col_group = "group",
#'   cols_to_compare = colnames(groups[, -1])[c(-1, -2)],
#'   type = c("co", "ca"), verbose = TRUE
#' )
#'
#' # Compare groups (Set NAs of categorical columns to 'Rest')
#' groups.cmp2 <- get_group_comparison(groups[, -1],
#'   col_group = "group",
#'   cols_to_compare = colnames(groups[, -1])[c(-1, -2)],
#'   type = c("co", "ca"), NAs = "Rest", verbose = TRUE
#' )
#'
#' ggcomp <- show_group_comparison(groups.cmp2)
#' ggcomp$co_comb
#' ggcomp$ca_comb
#' }
show_group_comparison <- function(group_comparison,
                                  xlab = "group", ylab_co = NA,
                                  legend_title_ca = NA,
                                  legend_position_ca = "bottom",
                                  show_pvalue = TRUE,
                                  method = "wilcox.test",
                                  p.adjust.method = "fdr",
                                  base_size = 12,
                                  font_size_x = 12,
                                  text_angle_x = 30,
                                  hjust = 0.2,
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

  # library(ggplot2)
  # library(cowplot)
  if (length(ca_index) > 0) {
    # plot categorical data
    ca_res <- lapply(ca_list, function(df) {
      data <- df[["data"]]
      data_sum <- data %>% dplyr::count_("group")
      data_sum[["labels"]] <- paste(data_sum[["group"]], paste0("(n=", data_sum[["n"]], ")"), sep = "\n")

      var_name <- colnames(data)[2]
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
          p <- p + labs(title = paste0("P=", signif(df[["p_value"]], 3)))
        }
      }

      p + theme(
        legend.position = df[["legend_position"]],
        axis.text.x = element_text(size = font_size_x, angle = text_angle_x, hjust = hjust)
      )
    })
    names(ca_res) <- names(ca_list)
  } else {
    ca_res <- NA
  }

  if (length(co_index) > 0) {
    # plot continuous data
    co_res <- lapply(co_list, function(df, ...) {
      data <- df[["data"]]
      data_sum <- data %>% dplyr::count_("group")
      data_sum[["labels"]] <- paste(data_sum[["group"]], paste0("(n=", data_sum[["n"]], ")"), sep = "\n")
      my_comparisons <- combn(unique(as.character(data[["group"]])),
        2,
        simplify = FALSE
      )

      var_name <- colnames(data)[2]
      var_name2 <- ifelse(isValidAndUnreserved(var_name), var_name, paste0("`", var_name, "`"))

      p <- ggplot(data, aes_string(x = "group", y = var_name2)) +
        geom_boxplot() + cowplot::theme_cowplot(font_size = base_size) +
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
        if (!requireNamespace("ggpubr")) {
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
      p + theme(axis.text.x = element_text(size = font_size_x, angle = text_angle_x, hjust = hjust))
    }, ...)
    names(co_res) <- names(co_list)
  } else {
    co_res <- NA
  }

  if (all(!is.na(ca_res))) {
    if (length(ca_res) <= 3) {
      ca_comb <- cowplot::plot_grid(plotlist = ca_res, align = "h")
    } else if (length(ca_res) == 4) {
      ca_comb <- cowplot::plot_grid(plotlist = ca_res, align = "h", ncol = 2)
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