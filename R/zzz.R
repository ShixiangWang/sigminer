.onAttach <- function(libname, pkgname) {
  version <- utils::packageDescription(pkgname, fields = "Version")

  msg <- paste0("=========================================================================================
", pkgname, " version ", version, "
Project URL: https://github.com/ShixiangWang/sigminer
Usages: https://shixiangwang.github.io/sigminer-book/ or builtin vignettes
---
If you use it in published research, please cite any proper paper :
  (1) Wang, Shixiang, et al. Copy number signature analysis tool and its application in
      prostate cancer reveals distinct mutational processes and clinical outcomes.
      PLoS genetics 17.5 (2021): e1009557.
  (2) Wang, Shixiang, et al. Sigflow: an automated and comprehensive pipeline
      for cancer genome mutational signature analysis. Bioinformatics 37.11 (2021): 1590-1592.
---
Suppress this message by `suppressPackageStartupMessages(library(sigminer))`
=========================================================================================
                              --Enjoy it--")
  base::packageStartupMessage(msg)
}

.onLoad <- function(libname, pkgname) {
  invisible(suppressPackageStartupMessages(
    sapply(c("NMF"),
      requireNamespace,
      quietly = TRUE
    )
  ))
}
