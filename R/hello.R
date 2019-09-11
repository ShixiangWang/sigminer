#' Say Hello to Users
#' @export
#' @examples
#' hello()
hello <- function() {
  if (!requireNamespace("cowsay", quietly = TRUE)) {
    message("Please install package 'cowsay' firstly.")
  } else {
    msg <- paste0("Welcome to 'sigminer' package!
=========================================================================
Quick overview please type ?sigminer

Github page: https://github.com/ShixiangWang/sigminer
Bug report: https://github.com/ShixiangWang/sigminer/issues
Documentation: https://github.com/ShixiangWang/sigminer
=========================================================================
                 ")
    cowsay::say(msg, by = "random")
  }
}
