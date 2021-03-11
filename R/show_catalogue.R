#' Show Alteration Catalogue Profile
#'
#' @param catalogue result from [sig_tally] or a
#' matrix with row representing components (motifs) and
#' column representing samples
#' @inheritParams show_sig_profile
#' @param normalize normalize method.
#' @param style plot style, one of 'default' and 'cosmic'.
#' @param samples default is `NULL`, show sum of all samples in one row.
#' If not `NULL`, show specified samples.
#' @param samples_name set the sample names shown in plot.
#' @param ... other arguments passing to [show_sig_profile].
#'
#' @return a `ggplot` object
#' @export
#'
#' @examples
#' data("simulated_catalogs")
#' p <- show_catalogue(simulated_catalogs$set1, style = "cosmic")
#' p
#' @testexamples
#' expect_s3_class(p, "ggplot")
show_catalogue <- function(catalogue,
                           mode = c("SBS", "copynumber", "DBS", "ID", "RS"),
                           method = "Wang",
                           normalize = c("raw", "row", "feature"),
                           style = c("default", "cosmic"),
                           samples = NULL,
                           samples_name = NULL,
                           x_lab = "Components",
                           y_lab = "Counts",
                           ...) {
  if (is.list(catalogue)) {
    if (!is.null(catalogue$nmf_matrix)) {
      mat <- catalogue$nmf_matrix %>% t()
    } else {
      stop("A list input can only be the result from sig_tally!")
    }
  } else {
    mat <- as.matrix(catalogue)
  }

  if (is.null(samples)) {
    ## If samples is NULL, sum all values
    mat <- mat %>%
      rowSums(na.rm = TRUE) %>%
      as.matrix()
    colnames(mat) <- "Total"
  } else {
    mat <- mat[, samples, drop = FALSE]
    ii <- colSums(mat) < 0.01
    if (any(ii)) {
      message(
        "The follow samples dropped due to null catalogue:\n\t",
        paste0(colnames(mat)[ii], collapse = ", ")
      )
      mat <- mat[, !ii, drop = FALSE]
      if (ncol(mat) < 1) {
        warning("No sample left. Please check mutation number of your input!", immediate. = TRUE)
        return(invisible(NULL))
      }
    }
  }

  mode <- match.arg(mode)
  normalize <- match.arg(normalize)
  style <- match.arg(style)

  show_sig_profile(mat,
    check_sig_names = FALSE,
    set_gradient_color = FALSE,
    mode = mode,
    method = method,
    style = style,
    normalize = normalize,
    x_lab = x_lab,
    y_lab = y_lab,
    sig_names = samples_name,
    sig_orders = samples_name,
    ...
  )
}
