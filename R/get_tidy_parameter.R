#' Get Tidy Parameter from Flexmix Model
#'
#' When users derive copy number features, it is useful to know the parameters
#' of the fit components, including mean, sd and coefficient of variation.
#' This function is used by [sig_tally] function and
#' exported to users for extra usage.
#'
#' @param x a `flexmix` object or a `list` of `flexmix` objects.
#'
#' @return a `tibble`.
#' @export
#'
#' @examples
#' load(system.file("extdata", "toy_copynumber_tally_M.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Get all parameters
#' d1 <- get_tidy_parameter(cn_tally_M$components)
#' d1
#' # Get parameters for segsize feature
#' d2 <- get_tidy_parameter(cn_tally_M$components$segsize)
#' d2
#' @testexamples
#' expect_equal(nrow(d1), 20L)
#' expect_equal(nrow(d2), 10L)
get_tidy_parameter <- function(x) {
  if (is.list(x)) {
    purrr::map_df(x, .get_parameter, .id = "feature") %>%
      dplyr::arrange(.data$feature, .data$mean) %>%
      dplyr::group_by(.data$feature) %>%
      dplyr::mutate(components = paste0(.data$feature, dplyr::row_number())) %>%
      dplyr::select(.data$feature, .data$components, dplyr::everything()) %>%
      dplyr::ungroup()
  } else if (assert_class(x, "flexmix")) {
    .get_parameter(x) %>%
      dplyr::arrange(.data$mean)
  }
}

.get_parameter <- function(x) {
  assert_class(x, "flexmix")
  paras <- flexmix::parameters(x)
  # weight is how many events assigned
  # to a cluster (component)
  # i.e. number of observations
  #
  # Info from package author:
  # the cluster sizes indicate the number
  # of observations assigned to each of the
  # clusters according to the a-posteriori probabilities.
  .get_weight <- function(mean, x) {
    wt_tb <- flexmix::clusters(x) %>%
      table()
    wt <- as.numeric(wt_tb)
    if (length(wt) == length(mean)) {
      return(wt)
    } else {
      names(wt) <- names(wt_tb)
      all_names <- seq_along(mean) %>%
        as.character()
      wt[setdiff(all_names, names(wt))] <- 0
      wt[sort(names(wt))] %>% as.numeric()
    }
  }

  if (is.null(dim(paras))) {
    # Assume it is pois distribution
    dplyr::tibble(
      dist = "pois",
      mean = as.numeric(paras),
      sd = sqrt(.data$mean)
    ) %>%
      dplyr::mutate(
        cv = .data$sd / .data$mean, # Coefficient of variation
        n_obs = .get_weight(.data$mean, x)
      )
  } else {
    # Assume it is normal distribution
    dplyr::tibble(
      dist = "norm",
      mean = as.numeric(paras[1, ]),
      sd = as.numeric(paras[2, ])
    ) %>%
      dplyr::mutate(
        cv = .data$sd / .data$mean, # Coefficient of variation
        n_obs = .get_weight(.data$mean, x)
      )
  }
}
