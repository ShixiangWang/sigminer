
#' Reset Signature Names
#'
#' @inheritParams show_sig_profile
#' @param names names to reset
#' @param order can be logical or a character vector.
#' If `TRUE`, signature names are sorted by [sort] function.
#' @param ... other parameters passing to [sort].
#'
#' @return a `Signature` object or a matrix.
#' @export
#'
#' @examples
#' print("yes")
reset_sig_names = function(Signature, names, order=TRUE) {
  if (inherits(Signature, "Signature")) {

  } else if (is.matrix(Signature)) {

  } else {
    stop("Invalid input!")
  }
}
