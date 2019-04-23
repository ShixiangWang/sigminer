.onAttach <- function(libname, pkgname) {
  version <- packageDescription(pkgname, fields = "Version")

  msg <- paste0("Welcome to 'sigminer' package!
======================================================
", pkgname, " version ", version, "
Github page: https://github.com/ShixiangWang/sigminer

More info please call 'hello()' in console.
======================================================
                 ")
  packageStartupMessage(msg)
}

.onLoad <- function(libname, pkgname) {
  # for (i in c("NMF")) {
  #     suppressMessages(requireNamespace(i, quietly = TRUE))
  # }
  invisible(suppressPackageStartupMessages(
    sapply(c("NMF"),
      requireNamespace,
      quietly = TRUE
    )
  ))
}
