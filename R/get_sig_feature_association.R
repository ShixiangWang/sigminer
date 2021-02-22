#' Calculate Association between Signature Exposures and Other Features
#'
#' Association of signature exposures with other features will be performed using one of two procedures:
#' for a continuous association variable (including ordinal variable), correaltion is performed;
#' for a binary association variable, samples will be divided into two groups and Mann-Whitney U-test
#' is performed to test for differences in signature exposure medians between the two groups.
#' See [get_tidy_association] for cleaning association result.
#'
#' @param data a `data.frame` contains signature exposures and other features
#' @param cols_to_sigs colnames for signature exposure
#' @param cols_to_features colnames for other features
#' @param type a character vector containing 'ca' for categorical variable and 'co' for continuous variable,
#' it must have the same length as `cols_to_features`.
#' @param method_co method for continuous variable, default is "spearman", could also be "pearson" and "kendall".
#' @param method_ca method for categorical variable, default is "wilcox.test"
#' @param min_n a minimal fraction (e.g. 0.01) or a integer number (e.g. 10) for filtering some variables with few positive events.
#' Default is 0.01.
#' @param verbose if `TRUE`, print extra message.
#' @param ... other arguments passing to test functions, like `cor.test`.
#'
#' @return a `list`. For 'co' features, 'measure' means correlation coefficient.
#' For 'ca' features, 'measure' means difference in means of signature exposure.
#' @export
#' @seealso [get_tidy_association]
get_sig_feature_association <- function(data, cols_to_sigs, cols_to_features, type = "ca",
                                        method_co = c("spearman", "pearson", "kendall"), method_ca = stats::wilcox.test,
                                        min_n = 0.01, verbose = FALSE, ...) {
  if (!all(type %in% c("ca", "co"))) {
    stop("all elements in 'type' must be 'ca' for 'categorical' variable and 'co' for 'continuous' variable.")
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame object.")
  }
  method_co <- match.arg(method_co)
  method_ca <- match.fun(method_ca)
  if (!dplyr::is.tbl(data)) {
    data <- dplyr::as_tibble(data)
  }
  data <- data[, c(cols_to_sigs, cols_to_features)]
  cols_co <- cols_to_features[type == "co"]
  cols_ca <- cols_to_features[type == "ca"]
  if (verbose) {
    message("-> Detecting and transforming possibly ordinal variables...")
  }
  data <- data %>% dplyr::mutate_if(is.ordered, as.numeric)
  data_co <- data[, c(cols_to_sigs, cols_co)]
  data_ca <- data[, c(cols_to_sigs, cols_ca)]
  pairwise_count <- function(x, y = NULL, diagonal = TRUE) {
    if (is.null(y)) {
      n <- t(!is.na(x)) %*% (!is.na(x))
    }
    else {
      n <- t(!is.na(x)) %*% (!is.na(y))
    }
    if (!diagonal & is.null(y)) {
      diag(n) <- NA
    }
    return(n)
  }
  matrix_cortest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat <- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tryCatch(
          {
            suppressWarnings(tmp <- cor.test(x = mat[
              ,
              i
            ], y = mat[, j], ...))
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
          },
          error = function(e) {
            p.mat[i, j] <- NA
          }
        )
      }
    }
    rownames(p.mat) <- colnames(mat)
    colnames(p.mat) <- colnames(mat)
    p.mat
  }
  cor.ca <- function(df, cols_ca, cols_to_sigs, .fun = stats::t.test,
                     min_n = 0.01, verbose = FALSE, ...) {
    nest_df <- lapply(cols_ca, function(x) {
      df[, c(cols_to_sigs, x)] %>%
        tidyr::gather_(
          "signature",
          "exposure", cols_to_sigs
        ) %>%
        dplyr::group_by(signature) %>%
        tidyr::nest()
    })
    names(nest_df) <- cols_ca
    test_fun <- function(df, .fun = stats::t.test, min_n = 0.05,
                         ...) {
      var_name <- setdiff(colnames(df), "exposure")
      df <- df[!is.na(df[["exposure"]]) & !is.na(df[[var_name]]), ]
      var_tb <- table(df[[var_name]])
      if (min_n < 1) {
        if (any(prop.table(var_tb) < min_n)) {
          return(NA)
        }
      }
      if (min_n > 1) {
        if (any(var_tb < min_n)) {
          return(NA)
        }
      }
      var_name2 <- ifelse(isValidAndUnreserved(var_name),
        var_name, paste0("`", var_name, "`")
      )
      f <- as.formula(paste("exposure", "~", var_name2))
      mod <- .fun(f, data = df, ...)
      mod$counts <- sum(var_tb, na.rm = TRUE)

      mod$measure <- df %>%
        split(., .[[var_name]]) %>%
        purrr::map_dbl(., ~ mean(.$exposure, na.rm = TRUE))

      # if (grepl("wilcox", body(.fun)[[2]])) {
      #   mod$measure <- df %>% split(., .[[var_name]]) %>%
      #     purrr::map_dbl(., ~median(.$exposure, na.rm = TRUE))
      # }
      # else if (grepl("t.test", body(.fun)[[2]])) {
      #   mod$measure <- df %>% split(., .[[var_name]]) %>%
      #     purrr::map_dbl(., ~mean(.$exposure, na.rm = TRUE))
      # }
      return(mod)
    }
    if (verbose) {
      message("--> building test models, be patient...")
    }
    tt <- lapply(nest_df, function(x, .fun, min_n, ...) {
      x %>%
        dplyr::ungroup() %>%
        dplyr::transmute(signature,
          model = purrr::map(data, test_fun,
            .fun = .fun,
            min_n = min_n
          )
        )
    }, .fun = .fun, min_n = min_n, ...)
    if (verbose) {
      message("--> done")
    }
    if (verbose) {
      message("--> obtaining model results...")
    }
    res <- lapply(tt, function(x) {
      x %>% dplyr::mutate(measure = purrr::map_dbl(
        .data$model,
        ~ ifelse(all(is.na(.)), NA, ifelse(is.null(.$measure),
          NA, diff(.$measure)
        ))
      ), count = purrr::map_int(
        .data$model,
        ~ ifelse(all(is.na(.)), NA, .$counts)
      ), p = purrr::map_dbl(
        .data$model,
        ~ ifelse(all(is.na(.)), NA, .$p.value)
      ))
    })
    if (verbose) {
      message("--> done")
    }
    if (verbose) {
      message("--> collecting data...")
    }
    res <- dplyr::bind_rows(res, .id = "features")
    res_list <- list()
    res_list$measure <- res[, c(
      "signature", "features",
      "measure"
    )] %>%
      tidyr::spread_(
        key_col = "features",
        value_col = "measure"
      ) %>%
      dplyr::arrange(match(
        signature,
        cols_to_sigs
      )) %>%
      dplyr::select(c("signature", cols_ca)) %>%
      tibble::column_to_rownames(var = "signature") %>%
      as.matrix() %>%
      signif(digits = 3)
    res_list$count <- res[, c("signature", "features", "count")] %>%
      tidyr::spread_(key_col = "features", value_col = "count") %>%
      dplyr::arrange(match(signature, cols_to_sigs)) %>%
      dplyr::select(c("signature", cols_ca)) %>%
      tibble::column_to_rownames(var = "signature") %>%
      as.matrix()
    res_list$p <- res[, c("signature", "features", "p")] %>%
      tidyr::spread_(key_col = "features", value_col = "p") %>%
      dplyr::arrange(match(signature, cols_to_sigs)) %>%
      dplyr::select(c("signature", cols_ca)) %>%
      tibble::column_to_rownames(var = "signature") %>%
      as.matrix() %>%
      signif(digits = 3)
    res_list$data <- res
    if (verbose) {
      message("--> done")
    }
    res_list
  }


  corr_co <- list()
  corr_ca <- list()


  if (length(cols_co) > 0) {
    if (verbose) {
      message(
        "-> Calculating correlations for continuous variables using ",
        method_co, "..."
      )
    }
    if (!all(sapply(data_co, is.numeric))) {
      stop("All continuous variables must be numeric (ordinal variable should be factor type)")
    }

    corr_co$measure <- cor(data_co,
      use = "pairwise.complete.obs",
      method = method_co
    )[cols_to_sigs, c(cols_to_sigs, cols_co)] %>%
      signif(digits = 3)
    corr_co$count <- pairwise_count(data_co)[cols_to_sigs, c(
      cols_to_sigs,
      cols_co
    )]
    corr_co$p <- matrix_cortest(data_co, method = method_co)[
      cols_to_sigs,
      c(cols_to_sigs, cols_co)
    ] %>% signif(digits = 3)
  }


  if (length(cols_ca) > 0) {
    if (verbose) {
      message(
        "-> Calculating correlations for categorical variables using ",
        as.list(method_ca)[[3]] %>% as.character() %>% .[[2]],
        "..."
      )
    }
    if (!all(sapply(data_ca[, cols_ca], is.logical))) {
      stop("All categorical variables must be logical, TRUE for positive events")
    }
    corr_ca <- cor.ca(data_ca, cols_ca, cols_to_sigs,
      .fun = method_ca,
      min_n = min_n, verbose = verbose, ...
    )
  }

  message("-> Done")
  list(corr_co = corr_co, corr_ca = corr_ca)
}
