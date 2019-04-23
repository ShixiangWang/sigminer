#' Say hello to users
#' @export
#' @examples
#' hello()
hello <- function() {
  if (!requireNamespace("cowsay", quietly = TRUE)) {
    stop("Package 'cowsay' is not installed.")
  } else {
    msg <- paste0("Welcome to 'sigminer' package!
=========================================================================
Quick overview please type ?sigminer

Github page: https://github.com/ShixiangWang/sigminer
Bug report: https://github.com/ShixiangWang/sigminer/issues
Documentation: https://github.com/ShixiangWang/sigminer

==
If you use 'sigminer' in academic field, please cite:

Wang, Shixiang, et al. \"The predictive power of tumor mutational burden
    in lung cancer immunotherapy response is influenced by patients' sex.\"
    International journal of cancer (2019).
==
Also cite the following if you extract copy number signatures

Macintyre, Geoff, et al. \"Copy number signatures and mutational
    processes in ovarian carcinoma.\" Nature genetics 50.9 (2018): 1262.
==
Also cite the following if you extract mutational signatures

Mayakonda, Anand, et al. \"Maftools: efficient and comprehensive analysis
    of somatic variants in cancer.\" Genome research 28.11 (2018): 1747-1756.
=========================================================================
                 ")
    cowsay::say(msg, by = "random")
  }
}
