#' Plot COSMIC Signature Profile
#'
#' @param sig_index a vector for signature index.
#' @param show_index if `TRUE`, show valid index.
#' @inheritParams get_sig_similarity
#' @param ... other arguments passing to [show_sig_profile].
#'
#' @return a `ggplot` object
#' @export
#' @author Shixiang Wang <w_shixiang@163.com>
#' @examples
#' show_cosmic_sig_profile()
#' show_cosmic_sig_profile(sig_db = "SBS")
#' show_cosmic_sig_profile(sig_index = 1:5)
#' show_cosmic_sig_profile(sig_db = "SBS", sig_index = c("10a", "17a"))
show_cosmic_sig_profile = function(sig_index=NULL, show_index=TRUE, sig_db = "legacy", ...) {

  sig_db <- match.arg(arg = sig_db, choices = c("legacy", "SBS"))

  if (sig_db == "legacy") {
    sigs_db <- readRDS(file = system.file("extdata", "legacy_signatures.RDs",
                                          package = "maftools", mustWork = TRUE
    ))
    sigs <- sigs_db$db
    aetiology <- sigs_db$aetiology
    avail_index = substring(colnames(sigs), 8)
  } else {
    sigs_db <- readRDS(file = system.file("extdata", "SBS_signatures.RDs",
                                          package = "maftools", mustWork = TRUE
    ))
    sigs <- sigs_db$db
    aetiology <- sigs_db$aetiology
    avail_index = substring(colnames(sigs), 4)
  }

  if (show_index) {
    msg =  paste(
      paste0("\nValid index for db '", sig_db, "':"),
      paste0(avail_index, collapse = " "),
      sep = "\n"
    )
    message(msg)
  }

  if (is.null(sig_index)) {
    return(invisible())
  } else {
    if (!is.character(sig_index)) {
      sig_index = as.character(sig_index)
    }
  }

  if (!all(sig_index %in% avail_index)) {
    stop(msg)
  }

  index = colnames(sigs)[avail_index %in% sig_index]
  sig_mat = as.matrix(sigs[, index, drop = FALSE])
  gg = show_sig_profile(sig_mat, check_sig_names = FALSE,
                        mode = "mutation", normalize = "raw", ...)
  gg$aetiology = aetiology[index, ,drop = FALSE]
  gg
}
