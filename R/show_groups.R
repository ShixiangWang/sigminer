#' Show Signature Contribution in Clusters
#'
#' See example section in [sig_fit()] for an examples.
#'
#' @param grp_dt a result `data.table` from [get_groups].
#' @param ... parameters passing to [legend()], e.g. `x = "topleft"`.
#'
#' @return nothing.
#' @export
#' @seealso [get_groups], [sig_fit].
show_groups <- function(grp_dt, ...) {
  x2 <- attr(grp_dt, "map_table")
  if (is.null(x2)) {
    send_stop("Wrong input, should be result from 'get_groups'.")
  }
  barplot(t(x2), legend.text = TRUE, args.legend = list(...))
}
