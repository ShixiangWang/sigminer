#' Read Copy Number Data from ASCAT Result Files
#'
#' Note, the result is not a `CopyNumber` object, you need to generate it
#' by yourself.
#'
#' @param x one or more `.rds` format files which contains `ASCAT` object from result of `ascat.runAscat()`
#' in [**ASCAT**](https://github.com/VanLoo-lab/ascat) package.
#'
#' @return a tidy `list`.
#' @export
read_copynumber_ascat <- function(x) {
  if (any(!file.exists(x))) {
    stop("file ", paste(x[!file.exists(x)], collapse = ", ")," could not be found!")
  }
  if (length(x) > 1) {
    return(read_ascat_cn_list(x))
  } else {
    return(read_ascat_cn(x))
  }
}

read_ascat_cn <- function(x) {
  message("reading ", x)
  x2 <- x
  x <- readRDS(x)

  if (object.size(x) < 2000) {
    warning(x2, " was failed in ASCAT calling.", immediate. = TRUE)
    return(invisible(NULL))
  }
  x <- x[c("segments", "aberrantcellfraction", "ploidy")]
  names(x) <- c("data", "purity", "ploidy")
  colnames(x$data) <- c("sample", "chromosome", "start", "end", "major_cn", "minor_cn")
  x$data <- x$data %>%
    dplyr::mutate(total_cn = .data$major_cn + .data$minor_cn) %>%
    dplyr::select(c("chromosome", "start", "end", "total_cn", "minor_cn", "sample"))
  x$data$source <- basename(x2) # track the source file
  return(x)
}

read_ascat_cn_list <- function(x_list) {
  message("reading file list:")
  x_list <- purrr::transpose(purrr::compact(lapply(x_list, read_ascat_cn)))
  y_list <- list()
  message("transforming data")
  y_list$data <- data.table::rbindlist(x_list$data)
  y_list$purity <- purrr::reduce(x_list$purity, c)
  y_list$ploidy <- purrr::reduce(x_list$ploidy, c)
  message("done")
  return(y_list)
}
