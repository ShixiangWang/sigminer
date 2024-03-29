#' Plot Reference (Mainly COSMIC) Signature Profile
#'
#' @param sig_index a vector for signature index. "ALL" for all signatures.
#' @param show_index if `TRUE`, show valid indices.
#' @inheritParams get_sig_similarity
#' @param ... other arguments passing to [show_sig_profile].
#'
#' @return a `ggplot` object
#' @export
#' @author Shixiang Wang <w_shixiang@163.com>
#' @examples
#' \donttest{
#' show_cosmic_sig_profile()
#' show_cosmic_sig_profile(sig_db = "SBS")
#' show_cosmic_sig_profile(sig_index = 1:5)
#' show_cosmic_sig_profile(sig_db = "SBS", sig_index = c("10a", "17a"))
#'
#' gg <- show_cosmic_sig_profile(sig_index = 1:5)
#' gg$aetiology
#' }
#' @testexamples
#' expect_s3_class(gg, "ggplot")
show_cosmic_sig_profile <- function(sig_index = NULL, show_index = TRUE, sig_db = "legacy", ...) {
  if (packageVersion("maftools") < "2.2.0") {
    send_error("This feature requires maftools >=2.2.0, please install it firstly!")
    return(invisible())
  }

  sigs_db <- get_sig_db(sig_db)
  sigs <- sigs_db$db
  aetiology <- sigs_db$aetiology

  avail_index <- switch(sig_db,
    legacy = substring(colnames(sigs), 8),
    SBS = substring(colnames(sigs), 4),
    DBS = substring(colnames(sigs), 4),
    ID = substring(colnames(sigs), 3),
    TSB = substring(colnames(sigs), 4),
    {
      if (startsWith(colnames(sigs)[1], "Ref_Sig")) {
        substring(colnames(sigs), 9)
      } else {
        colnames(sigs)
      }
    }
  )

  msg <- paste(
    paste0("\nValid index for db '", sig_db, "':"),
    paste0(avail_index, collapse = " "),
    sep = "\n"
  )

  if (show_index) {
    message(msg)
  }

  if (is.null(sig_index)) {
    return(invisible())
  } else {
    if (!is.character(sig_index)) {
      sig_index <- as.character(sig_index)
    }

    if ("ALL" %in% sig_index) {
      sig_index <- avail_index
    }

    if (length(sig_index) == 1) {
      if (grepl(",", sig_index)) {
        sig_index <- unique(split_seq(sig_index))
      }
    }
  }

  if (!all(sig_index %in% avail_index)) {
    stop(msg)
  }


  index <- c()
  for (i in sig_index) {
    index <- c(index, colnames(sigs)[avail_index == i])
  }

  sig_mat <- as.matrix(sigs[, index, drop = FALSE])

  if (sig_db %in% c("legacy", "SBS", "TSB")) {
    mode <- "SBS"
  } else if (startsWith(sig_db, "SBS")) {
    mode <- "SBS"
  } else if (startsWith(sig_db, "RS")) {
    mode <- "RS"
  } else if (startsWith(sig_db, "CNS")) {
    mode <- "copynumber"
  } else {
    mode <- sig_db
  }

  gg <- show_sig_profile(sig_mat,
    check_sig_names = FALSE,
    mode = mode,
    normalize = "row",
    sig_orders = index,
    ...
  )
  gg$aetiology <- aetiology[index, , drop = FALSE]
  gg
}
