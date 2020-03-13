#' Get Ploidy from Absolute Copy Number Profile
#'
#' @inheritParams show_cn_profile
#'
#' @return a value or a `data.table`
#' @export
#'
#' @examples
#' # Load copy number object
#' load(system.file("extdata", "toy_copynumber.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#'
#' df <- get_cn_ploidy(cn)
#' df
#' @testexamples
#' expect_equal(nrow(df), 10L)
get_cn_ploidy <- function(data) {
  stopifnot(is.data.frame(data) | inherits(data, "CopyNumber"))
  if (is.data.frame(data)) {
    nc_cols <- c("chromosome", "start", "end", "segVal")
    if (!all(nc_cols %in% colnames(data))) {
      stop("Invalid input, it must contain columns: ", paste(nc_cols, collapse = " "))
    }
  } else {
    data <- data@data
  }

  has_sample <- "sample" %in% colnames(data)
  data.table::setDT(data)
  data$segLen <- data$end - data$start + 1

  if (has_sample) {
    ploidy <- data %>%
      dplyr::as_tibble() %>%
      dplyr::group_by(.data$sample) %>%
      dplyr::summarise(
        ploidy = sum((.data$segLen / sum(.data$segLen)) * .data$segVal)
      ) %>%
      dplyr::mutate(
        ploidy = round(.data$ploidy, 2)
      ) %>%
      data.table::as.data.table()
  } else {
    ploidy <- sum((data$segLen / sum(data$segLen)) * data$segVal)
    ploidy <- round(ploidy, 2)
  }

  ploidy
}
