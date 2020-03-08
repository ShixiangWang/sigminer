#' Say Hello to Users
#' @export
#' @examples
#' hello()
hello <- function() {
  msg <- paste0("Welcome to 'sigminer' package!
=========================================================================
Quick overview please type ?sigminer

Github page: https://github.com/ShixiangWang/sigminer
Bug report: https://github.com/ShixiangWang/sigminer/issues
Documentation: https://github.com/ShixiangWang/sigminer
=========================================================================
                 ")
  message(msg)
}
