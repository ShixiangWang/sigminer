#' Show Signature Profile
#'
#' Who don't like to show a barplot for signature profile? This is for it.
#'
#' @inheritParams sig_extract
#' @inheritParams show_cn_distribution
#' @param Signature a `Signature` object obtained either from [sig_extract] or [sig_auto_extract],
#' or just a raw signature matrix with row representing components (motifs) and column
#' representing signatures (column names must start with 'Sig').
#' @param mode signature type for plotting, now supports 'copynumber' or 'mutation'.
#' @param normalize one of 'row', 'column' and 'raw', for row normalization (signature),
#' column normalization (component) and raw data, respectively.
#' @param params params `data.frame` of components, obtained from [derive].
#' @param params_label_size font size for params label.
#' @param params_label_angle font angle for params label.
#' @param y_expand y expand height for plotting params of copy number signatures.
#' @param digits digits for plotting params of copy number signatures.
#' @param font_scale a number used to set font scale.
#' @param sig_names set name of signatures, can be a character vector.
#' Default is `NULL`, prefix 'Signature_' plus number is used.
#' @param sig_orders set order of signatures, can be a character vector.
#' Default is `NULL`, the signatures are ordered by alphabetical order.
#' @author Shixiang Wang
#' @return a `ggplot` object
#' @export
#' @examples
#' \donttest{
#' # Load mutational signature
#' load(system.file("extdata", "toy_mutational_signature.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Show signature profile
#' show_sig_profile(sig2, mode = "mutation")
#'
#'
#' # Load copy number signature
#' load(system.file("extdata", "toy_copynumber_signature.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Show signature profile
#' show_sig_profile(sig)
#'
#' # Add params label
#' #=================
#' # Load copy number prepare object
#' load(system.file("extdata", "toy_copynumber_prepare.RData",
#'                  package = "sigminer", mustWork = TRUE
#'                  ))
#' show_sig_profile(sig, params = cn_prepare$parameters, y_expand = 2)
#'
#' }
show_sig_profile <- function(Signature, mode = c("copynumber", "mutation"),
                             normalize = c("row", "column", "raw"),
                             params = NULL, params_label_size = 3,
                             params_label_angle = 60, y_expand = 1,
                             digits = 2, base_size = 12, font_scale = 1,
                             sig_names = NULL, sig_orders = NULL) {
  if (class(Signature) == "Signature") {
    Sig = Signature$Signature
  } else if (is.matrix(Signature)) {
    if (!all(startsWith(colnames(Signature), "Sig"))) {
      stop("If Signature is a matrix, column names must start with 'Sig'!", call. = FALSE)
    }
    Sig = Signature
  } else {
    stop("Invalid input for 'Signature'", call. = FALSE)
  }

  mode <- match.arg(mode)
  normalize = match.arg(normalize)

  if (normalize == "row") {
    Sig <- apply(Sig, 2, function(x) x / sum(x))
  } else if (normalize == "column") {
    Sig <- t(apply(Sig, 1, function(x) x / sum(x)))
  }

  # >>>>>>>>>>>>>>>>> Setting theme
  scale <- font_scale

  .theme_ss <- theme_bw(base_size = base_size) +
    theme(
      axis.text.x = element_text(
        angle = 60, vjust = 1,
        hjust = 1, size = (base_size - 4) * scale
      ),
      axis.text.y = element_text(hjust = 0.5, size = base_size * scale)
    )
  # <<<<<<<<<<<<<<<<< Setting theme

  # >>>>>>>>>>>>>>>>> identify mode and do data transformation
  mat <- as.data.frame(Sig)
  mat$context <- rownames(mat)

  if (mode == "copynumber") {
    mat$base <- sub("\\d+$", "", mat$context)

    mat <- tidyr::gather(mat, class, signature, dplyr::contains("Sig"))
    mat <- dplyr::mutate(mat,
                         context = factor(.data$context,
                                          levels = unique(mat[["context"]])
                         ),
                         base = factor(.data$base, levels = c(
                           "bp10MB", "copynumber",
                           "changepoint", "bpchrarm",
                           "osCN", "segsize"
                         )),
                         class = factor(class)
    )
  } else {
    mat$base <- sub("[ACGT]\\[(.*)\\][ACGT]", "\\1", mat$context)
    mat$context <- sub("(\\[.*\\])", "\\-", mat$context)

    mat <- tidyr::gather(mat, class, signature, dplyr::contains("Sig"))
    mat <- dplyr::mutate(mat,
                         context = factor(.data$context),
                         base = factor(.data$base, levels = c(
                           "C>A", "C>G",
                           "C>T", "T>A",
                           "T>C", "T>G"
                         )),
                         class = factor(class)
    )
  }

  # >>>>>>> Set signature name and order
  if (!is.null(sig_names)) {
    if (length(sig_names) != length(unique(mat[["class"]]))) {
      stop("The length of input signature names is not equal to signature number")
    }
    names(sig_names) <- paste0("Signature_", seq_along(sig_names))
    mat[["class"]] <- sig_names[mat[["class"]]]
  }

  if (!is.null(sig_orders)) {
    mat[["class"]] <- factor(mat[["class"]], levels = sig_orders)
  } else {
    mat[["class"]] <- factor(mat[["class"]])
  }
  # >>>>>>>>>>>>>>>>>>>>>>> Plot
  p <- ggplot(mat) +
    geom_bar(aes_string(x = "context", y = "signature", fill = "base"),
             stat = "identity", position = "identity", colour = "gray50"
    ) +
    scale_fill_manual(values = c(
      "cyan", "red", "yellow", "purple",
      "green", "blue", "black", "gray"
    ))

  if (mode == "copynumber") {
    if (!is.null(params)) {
      params$class <- factor(levels(mat[["class"]])[1], levels = levels(mat[["class"]]))
      p <- p + geom_text(aes(
        x = .data$components, y = Inf,
        label = ifelse(dist == "norm",
                       paste0(" \u03BC=", signif(.data$mean, digits)),
                       paste0(" \u03BB=", signif(.data$mean, digits))
        )
      ),
      data = params,
      size = params_label_size,
      angle = params_label_angle,
      hjust = 0, vjust = 0.5
      ) +
        coord_cartesian(clip = "off")
    }
    p <- p + facet_grid(class ~ ., scales = "free")
  } else {
    p <- p + facet_grid(class ~ base, scales = "free")
  }

  p <- p +
    guides(fill = FALSE) + .theme_ss

  if (all(mode == "copynumber", !is.null(params))) {
    p <- p + theme(plot.margin = margin(30 * y_expand, 2, 2, 2, unit = "pt")) # Add regions
  }

  p <- p + xlab("Components") + ylab("Contributions")

  return(p)
}
