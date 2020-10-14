#' Show Signature Profile with Loop Way
#'
#' @inheritParams show_sig_profile
#' @inheritParams cowplot::plot_grid
#' @param ... other parameters but `sig_order` passing to [show_sig_profile].
#'
#' @return a `ggplot` result from `cowplot::plot_grid()`.
#' @export
#' @seealso [show_sig_profile]
#'
#' @examples
#'
#' load(system.file("extdata", "toy_mutational_signature.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Show signature profile
#' p1 <- show_sig_profile_loop(sig2, mode = "SBS")
#' p1
#' p2 <- show_sig_profile_loop(sig2, mode = "SBS", style = "cosmic", sig_names = c("A", "B", "C"))
#' p2
#' @testexamples
#' expect_s3_class(p1, "ggplot")
#' expect_s3_class(p2, "ggplot")
show_sig_profile_loop <- function(Signature, sig_names = NULL,
                                  ncol = 1,
                                  nrow = NULL,
                                  x_lab = "Components",
                                  ...) {
  if (inherits(Signature, "Signature")) {
    Sig <- Signature$Signature
  } else if (is.matrix(Signature)) {
    Sig <- Signature
  } else {
    stop("Invalid input for 'Signature'", call. = FALSE)
  }

  x_lab_bk <- x_lab
  nc <- ncol(Sig)
  plist <- list()

  for (i in seq_len(nc)) {
    if (i != nc) {
      x_lab <- NULL
    } else {
      x_lab <- x_lab_bk
    }

    plist[[i]] <- show_sig_profile(
      Signature = Sig,
      x_lab = x_lab,
      sig_names = if (is.null(sig_names)) NULL else sig_names[i],
      sig_orders = i,
      ...
    )
  }
  # cannot modify gglotify ggplot
  # plist[1:(nc-1)] <- lapply(plist[1:(nc-1)], function(x) x + labs(x = NULL))
  p <- cowplot::plot_grid(plotlist = plist,
                          ncol = ncol,
                          nrow = nrow)
                          #align = "hv",
                          #axis = "tblr")

  return(p)
}
