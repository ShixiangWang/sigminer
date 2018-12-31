.onAttach = function(libname, pkgname) {
    version = packageDescription(pkgname, fields = "Version")

    msg = paste0("=========================================================================
", pkgname, " version ", version, "
Github page: https://github.com/ShixiangWang/sigminer
Documentation: https://github.com/ShixiangWang/sigminer
If you use it in published research, please cite:
Macintyre, Geoff, et al. \"Copy number signatures and mutational
    processes in ovarian carcinoma.\" Nature genetics 50.9 (2018): 1262.
Wang, Shixiang, et al. \"APOBEC3B and APOBEC mutational signature
    as potential predictive markers for immunotherapy
    response in non-small cell lung cancer.\" Oncogene (2018).
=========================================================================
                 ")
    packageStartupMessage(msg)
}

.onLoad <- function(libname, pkgname) {
    # for (i in c("NMF")) {
    #     suppressMessages(requireNamespace(i, quietly = TRUE))
    # }
    invisible(suppressPackageStartupMessages(
        sapply(c("NMF"),
               requireNamespace, quietly = TRUE)
    ))
}
