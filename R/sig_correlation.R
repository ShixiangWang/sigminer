#' Calculate association between signature activity and other features
#'
#' @param data a `data.frame` contains signature activity and other features
#' @param cols_to_sigs colnames for signature activity
#' @param cols_to_summary colnames for other features
#' @param type a character vector contains 'ca' for categorical variable and 'co' for continuous variable,
#' it must have same length as `cols_to_summary`.
#' @param method_co method for continuous variable, default is "spearman", could also be "pearson" and "kendall".
#' @param method_ca method for categorical variable, default is "wilcox.test"
#' @param min_n a minimal fraction (e.g. 0.01) or a integer number (e.g. 10) for filtering some variables with few positive events.
#' Default is 0.01.
#' @param verbose if `TRUE`, print extra message.
#' @param ... other arguments passing to
#'
#' @return a `list`
#' @export
#'
sig_summarize_correlation <- function(data, cols_to_sigs, cols_to_summary,
                                      type = "ca",
                                      method_co = c("spearman", "pearson", "kendall"),
                                      method_ca = stats::wilcox.test,
                                      min_n=0.01, verbose = FALSE, ...) {
  if (!all(type %in% c("ca", "co"))) {
    stop("all elements in 'type' must be 'ca' for 'categorical' variable and 'co' for 'continuous' variable.")
  }
  if (!is.data.frame(data)) stop("'data' must be a data.frame object.")
  method_co = match.arg(method_co)
  method_ca = match.fun(method_ca)

  if (!dplyr::is.tbl(data)) data = dplyr::as_tibble(data)

  data <- data[, c(cols_to_sigs, cols_to_summary)]


  cols_co = cols_to_summary[type == "co"]
  cols_ca = cols_to_summary[type == "ca"]

  if (verbose) message("-> Detecting and transforming possibly ordinal variables...")
  data <- data %>%
    dplyr::mutate_if(is.ordered, as.numeric)

  data_co = data[, c(cols_to_sigs, cols_co)]
  data_ca = data[, c(cols_to_sigs, cols_ca)]

  # copy from psych package
  pairwise_count = function (x, y = NULL, diagonal = TRUE)
  {
    if (is.null(y)) {
      n <- t(!is.na(x)) %*% (!is.na(x))
    }
    else {
      n <- t(!is.na(x)) %*% (!is.na(y))
    }
    if (!diagonal & is.null(y))
      diag(n) <- NA
    return(n)
  }

  matrix_cortest = function (mat, ...)
  {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat <- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tryCatch({
          tmp <- cor.test(x = mat[, i], y = mat[, j], ...)
          p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        }, error = function(e) {
          p.mat[i, j] <- NA
        })
      }
    }
    rownames(p.mat) <- colnames(mat)
    colnames(p.mat) <- colnames(mat)
    p.mat
  }

  cor.ca = function(df, cols_ca, cols_to_sigs, .fun=stats::t.test, min_n=0.01, verbose=FALSE, ...) {
    # Split data into two groups
    # and compare them with
    # statistical tests
    nest_df = lapply(cols_ca, function(x) {
      df[, c(cols_to_sigs, x)] %>%
        tidyr::gather_("signature", "exposure", cols_to_sigs) %>%
        dplyr::group_by(signature) %>%
        tidyr::nest()
    })
    names(nest_df) = cols_ca

    # 构建检验模型
    # min_n 可以是一个小于1的值，表示比例
    # 或者一个整数
    test_fun = function(df, .fun=stats::t.test, min_n=0.05, ...) {
      var_name = setdiff(colnames(df), "exposure")
      # Filter NAs
      df = df[!is.na(df[["exposure"]]) & !is.na(df[[var_name]]) ,]
      var_tb = table(df[[var_name]])
      if (min_n < 1) {
        if (any(prop.table(var_tb) < min_n)) return(NA)
      }
      if (min_n > 1) {
        if (any(var_tb < min_n)) return(NA)
      }
      var_name = ifelse(isValidAndUnreserved(var_name), var_name, paste0("`",var_name,"`"))
      f = as.formula(paste("exposure", "~", var_name))
      mod = .fun(f, data = df, ...)
      mod$counts=sum(var_tb, na.rm = TRUE)
      return(mod)
    }

    if (verbose) message("--> building test models, be patient...")
    tt = lapply(nest_df, function(x, .fun, min_n, ...){
      x %>%
        dplyr::transmute(signature, model = purrr::map(data, test_fun, .fun=.fun, min_n=min_n))
    }, .fun=.fun, min_n=min_n, ...)
    if (verbose) message("--> done")

    if (verbose) message("--> obtaining model results...")
    # measure在t检验里计算的是均值之差
    # 使用wilcox test没有estimate
    res = lapply(tt, function(x) {
      x %>%
        dplyr::mutate(measure=purrr::map_dbl(model, ~ifelse(all(is.na(.)), NA,
                                                            ifelse(is.null(.$estimate), NA, diff(.$estimate)))),
                      count=purrr::map_int(model, ~ifelse(all(is.na(.)), NA, .$counts)),
                      p=purrr::map_dbl(model, ~ifelse(all(is.na(.)), NA, .$p.value))
        )
    })
    if (verbose) message("--> done")

    if (verbose) message("--> collecting data...")
    # 提取数据组成矩阵
    res = dplyr::bind_rows(res, .id="features")

    res_list = list()
    res_list$measure = res[, c("signature", "features", "measure")] %>%
      tidyr::spread_(key_col = "features", value_col="measure") %>%
      dplyr::arrange(match(signature, cols_to_sigs)) %>%
      dplyr::select(c("signature", cols_ca)) %>%
      tibble::column_to_rownames(var="signature") %>% as.matrix() %>% signif(digits = 3)
    res_list$count = res[, c("signature", "features", "count")] %>%
      tidyr::spread_(key_col = "features", value_col="count") %>%
      dplyr::arrange(match(signature, cols_to_sigs)) %>%
      dplyr::select(c("signature", cols_ca)) %>%
      tibble::column_to_rownames(var="signature") %>% as.matrix()
    res_list$p = res[, c("signature", "features", "p")] %>%
      tidyr::spread_(key_col = "features", value_col="p") %>%
      dplyr::arrange(match(signature, cols_to_sigs)) %>%
      dplyr::select(c("signature", cols_ca)) %>%
      tibble::column_to_rownames(var="signature") %>% as.matrix() %>% signif(digits = 3)
    res_list$data = res

    if (verbose) message("--> done")
    res_list
  }

  if (verbose) message("-> Calculating correlations for continuous variables using ", method_co, "...")
  if (!all(sapply(data_co, is.numeric))) {
    stop("All continuous variables must be numeric (ordinal variable should be factor type)")
  }

  corr_co = list()
  corr_co$measure = cor(data_co, use = "pairwise.complete.obs",
                        method = method_co)[cols_to_sigs,
                                            c(cols_to_sigs, cols_co)] %>%
    signif(digits = 3)
  corr_co$count = pairwise_count(data_co)[cols_to_sigs, c(cols_to_sigs, cols_co)]
  corr_co$p = matrix_cortest(data_co, method = method_co)[cols_to_sigs,
                                                          c(cols_to_sigs, cols_co)] %>%
    signif(digits = 3)

  if (verbose) message("-> Calculating correlations for categorical variables using ",
                       as.list(method_ca)[[3]] %>% as.character() %>% .[[2]], "...")
  if (!all(sapply(data_ca[, cols_ca], is.logical))) {
    stop("All categorical variables must be logical, TRUE for positive events")
  }
  corr_ca = cor.ca(data_ca, cols_ca, cols_to_sigs, .fun=method_ca, min_n= min_n, verbose = verbose, ...)

  message("-> Done")
  list(corr_co=corr_co, corr_ca=corr_ca)
}


#' Tidy signature correlation results
#'
#' @param cor_res data returned by [sig_summarize_correlation()]
#'
#' @return a `data.frame`
#' @export
sig_tidy_correlation = function(cor_res) {
  all_names= names(cor_res[["corr_ca"]])
  common_names = setdiff(all_names, "data")
  res = list()
  tidy = list()
  for (i in common_names) {
    if (!identical(rownames(cor_res[["corr_co"]][[i]]), rownames(cor_res[["corr_ca"]][[i]]))) {
      stop("Bad rowname order, please report this issue to developer with your sample data")
    }
    res[[i]] = cbind(cor_res[["corr_co"]][[i]], cor_res[["corr_ca"]][[i]])
    tidy[[i]] = res[[i]] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("signature") %>%
      tibble::as_tibble() %>%
      tidyr::gather_(key_col = "features", value_col = "value", gather_cols = colnames(res[[i]]))
  }

  # To avoid Inf,
  # set p == 0 to the same value
  # as p value among p != 0
  res[["data"]] = dplyr::bind_rows(tidy, .id = "data_type") %>%
    tidyr::spread_(key_col = "data_type", value_col = "value") %>%
    dplyr::mutate(score = ifelse(p == 0, -log10(min(p[p!=0], na.rm = TRUE)), -log10(p)),
                  score = ifelse(measure<0, -score, score))
  res
}


#' Draw corrplot for signature activity and other features
#'
#' @param tidy_cor data returned by [sig_tidy_correlation()].
#' @param features_list a character vector contains features want to be plotted.
#' @param breaks_count breaks for sample count.
#' @param breaks_score breaks for score.
#'
#' @return a `ggplot2` object
#' @export
draw_features_corrplot = function(tidy_cor, features_list,
                                  breaks_count = c(0L, 200L, 400L, 600L, 800L, 1020L),
                                  breaks_score = c(-Inf, -16, -10, -5, -1.3, 1.3, 5, 10, 16, Inf)) {
  data = tidy_cor$data %>%
    dplyr::mutate(
      Samples = cut(count, breaks = breaks_count),
      Score = cut(score, breaks = breaks_score)
    ) %>%
    dplyr::filter(features %in% features_list, !is.na(score))

  ggplot2::ggplot(data,
                  ggplot2::aes(x = signature,
                               y = factor(features, levels = rev(features_list)))) +
    ggplot2::geom_point(ggplot2::aes(colour = Score,
                                     size = Samples))  +
    ggplot2::scale_color_manual(values = rev(RColorBrewer::brewer.pal(9,"RdBu")), drop=FALSE) +
    ggplot2::scale_size_discrete(drop=FALSE) +
    cowplot::theme_cowplot() + ggpubr::rotate_x_text(45) +
    ggplot2::labs(x="Signatures", y ="Features")
}

utils::globalVariables(
  c("features", "Score", "Samples", "model", "p", "measure")
)
