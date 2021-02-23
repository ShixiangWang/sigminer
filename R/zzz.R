.onAttach <- function(libname, pkgname) {
  version <- packageDescription(pkgname, fields = "Version")

  msg <- paste0(pkgname, " version ", version, "\n- Star me at https://github.com/ShixiangWang/sigminer\n- Run hello() to see usage and citation.")
  packageStartupMessage(msg)
}

.onLoad <- function(libname, pkgname) {
  invisible(suppressPackageStartupMessages(
    sapply(c("NMF"),
      requireNamespace,
      quietly = TRUE
    )
  ))
}
