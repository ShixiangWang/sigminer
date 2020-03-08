.onAttach <- function(libname, pkgname) {
  version <- packageDescription(pkgname, fields = "Version")

  msg <- paste0(pkgname, " version ", version, ", run hello() to see more.")
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
