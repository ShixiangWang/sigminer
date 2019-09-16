#' Get Tidy Parameter from Flexmix Model
#'
#' When users derive copy number features, it is useful to know the parameters
#' of the fit components. This function is used by [derive] function and
#' exported to users for extra usage.
#'
#' @param x a `flexmix` object or a `list` of `flexmix` objects.
#'
#' @return a `tibble`.
#' @export
#'
#' @examples
#' # Load copy number prepare object
#' load(system.file("extdata", "toy_copynumber_prepare.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Get all parameters
#' get_tidy_parameter(cn_prepare$components)
#' # Get parameters for segsize feature
#' get_tidy_parameter(cn_prepare$components$segsize)
get_tidy_parameter = function(x) {
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

.get_parameter = function(x) {
  assert_class(x, "flexmix")
  paras = flexmix::parameters(x)
  if (is.null(dim(paras))) {
    # Assume it is pois distribution
    dplyr::tibble(
      dist = "pois",
      mean = as.numeric(paras),
      sd = sqrt(.data$mean)
    )
  } else {
    # Assume it is normal distribution
    dplyr::tibble(
      dist = "norm",
      mean = as.numeric(paras[1, ]),
      sd = as.numeric(paras[2, ])
    )
  }
}
