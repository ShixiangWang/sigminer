#' Say hello to users
#' @export
#' @examples
#' hello()
hello <- function() {
  if (!requireNamespace("cowsay",  quietly = TRUE)) {
    stop("Package 'cowsay' is not installed.")
  } else {
    msg = paste0("Welcome to 'sigminer' package!
=========================================================================
Quick overview please type ?sigminer

Github page: https://github.com/ShixiangWang/sigminer
Bug report: https://github.com/ShixiangWang/sigminer/issues
Documentation: https://github.com/ShixiangWang/sigminer

If you use 'sigminer' in academic field, please cite:
Macintyre, Geoff, et al. \"Copy number signatures and mutational
    processes in ovarian carcinoma.\" Nature genetics 50.9 (2018): 1262.
Wang, Shixiang, et al. \"APOBEC3B and APOBEC mutational signature
    as potential predictive markers for immunotherapy
    response in non-small cell lung cancer.\" Oncogene (2018).
=========================================================================
                 ")
    cowsay::say(msg, by = "random")
  }
}
