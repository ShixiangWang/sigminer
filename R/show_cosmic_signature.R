#' Show Signature Information in Web Browser
#'
#' @param x a string indicating location
#' ("home" for COSMIC signature home, "legacy" for COSMIC v2 signatures,
#' "SBS" for COSMIC v3 SBS signatures, "DBS" for COSMIC v3 DBS signatures,
#' "ID" for COSMIC v3 INDEL signatures) or signature index (e.g.
#' "SBS1", "DBS2", "ID3").
#'
#' @return Nothing.
#' @export
#'
#' @examples
#' \dontrun{
#' show_cosmic()
#' show_cosmic("legacy")
#' show_cosmic("SBS")
#' show_cosmic("DBS")
#' show_cosmic("ID")
#' show_cosmic("SBS1")
#' show_cosmic("DBS2")
#' show_cosmic("ID3")
#' }
show_cosmic <- function(x = "home") {
  stopifnot(length(x) == 1L, is.character(x))
  url <- switch(x,
                home = "https://cancer.sanger.ac.uk/cosmic/signatures",
                legacy = "https://cancer.sanger.ac.uk/cosmic/signatures_v2.tt",
                SBS = "https://cancer.sanger.ac.uk/cosmic/signatures/SBS/index.tt",
                DBS = "https://cancer.sanger.ac.uk/cosmic/signatures/DBS/index.tt",
                ID = "https://cancer.sanger.ac.uk/cosmic/signatures/ID/index.tt",
                parse_cosmic(x))
  message("Opening your default web browser...")
  utils::browseURL(url)
}

parse_cosmic <- function(x) {
  if (startsWith(x, "SBS")) {
    paste0("https://cancer.sanger.ac.uk/cosmic/signatures/SBS/", x, ".tt")
  } else if (startsWith(x, "DBS")) {
    paste0("https://cancer.sanger.ac.uk/cosmic/signatures/DBS/", x, ".tt")
  } else if (startsWith(x, "ID")) {
    paste0("https://cancer.sanger.ac.uk/cosmic/signatures/ID/", x, ".tt")
  } else {
    stop("Signature name starts with 'SBS', 'DBS' or 'ID'!")
  }
}
