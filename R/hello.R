#' Say hello to users
#' @export
#' @examples
#'  hello()
hello <- function() {
  if (!requireNamespace("cowsay",  quietly = TRUE)) {
    stop("Package 'cowsay' is not installed.")
  } else {
    cowsay::say("Welcome to 'sigminer' package!", by = "spider")
  }
}
