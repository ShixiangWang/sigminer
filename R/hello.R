#' Say Hello to Users
#' @export
#' @examples
#' hello()
#' @testexamples
#' expect_message(hello(), "sigminer")
hello <- function() {
  msg <- paste0("Thanks for using 'sigminer' package!
=========================================================================
Version: ", packageVersion("sigminer"), "
Run citation('sigminer') to see how to cite sigminer in publications.

Project home : https://github.com/ShixiangWang/sigminer
Bug report   : https://github.com/ShixiangWang/sigminer/issues
Documentation: https://shixiangwang.github.io/sigminer-doc/
=========================================================================
                 ")
  message(msg)
}
