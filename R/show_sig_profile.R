#' Show Signature Profile
#'
#' Who don't like to show a barplot for signature profile? This is for it.
#'
#' @inheritParams sig_extract
#' @inheritParams show_cn_distribution
#' @param params params data of components.
#' @param params_label_size font size for params label.
#' @param y_expand y expand height for plotting params of copy number signatures.
#' @param digits digits for plotting params of copy number signatures.
#' @param y_scale one of 'relative' or 'absolute', if choose 'relative',
#' signature columns will be scaled to sum as 1.
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
#' # Load copy number signature
#' #// Remove this data and create a new one
#' load(system.file("extdata", "toy_copynumber_signature.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' show_sig_profile(res$nmfObj)
#' }
show_sig_profile <- function(nmfObj, mode = c("copynumber", "mutation"),
                             params = NULL, params_label_size = 3, y_expand = 1,
                             digits = 1, base_size = 12,
                             y_scale = c("relative", "absolute"), font_scale = 1,
                             sig_names = NULL, sig_orders = NULL) {
  mode <- match.arg(mode)
  y_scale <- match.arg(y_scale)

  # Signatures
  w <- NMF::basis(nmfObj)
  if (y_scale == "relative") {
    if (mode == "copynumber") {
      w <- t(apply(w, 1, function(x) x / sum(x))) # Scale the component weight
    } else if (mode == "mutation") {
      w <- apply(w, 2, function(x) x / sum(x)) # Scale the signatures
    }
  }
  colnames(w) <- paste("Signature", 1:ncol(w), sep = "_")

  # >>>>>>>>>>>>>>>>> Setting theme
  scale <- font_scale

  .theme_ss <- theme_bw(base_size = base_size) +
    theme(
      axis.text.x = element_text(
        angle = 90, vjust = 0.5,
        hjust = 1, size = (base_size - 2) * scale
      ),
      axis.text.y = element_text(hjust = 0.5, size = base_size * scale)
    )
  # <<<<<<<<<<<<<<<<< Setting theme

  # >>>>>>>>>>>>>>>>> identify mode and do data transformation
  # w = rbind(w, matrix(c(0.1, 0.2), ncol = 2, dimnames = list("segsize10")))
  mat <- as.data.frame(w)
  mat$context <- rownames(mat)

  if (mode == "copynumber") {
    mat$base <- sub("\\d+$", "", mat$context)

    mat <- tidyr::gather(mat, class, signature, dplyr::contains("Signature"))
    mat <- dplyr::mutate(mat,
                         context = factor(context,
                                          levels = unique(mat[["context"]])
                         ),
                         base = factor(base, levels = c(
                           "bp10MB", "copynumber",
                           "changepoint", "bpchrarm",
                           "osCN", "segsize"
                         )),
                         class = factor(class)
    )
  } else {
    mat$base <- sub("[ACGT]\\[(.*)\\][ACGT]", "\\1", mat$context)
    mat$context <- sub("(\\[.*\\])", "\\.", mat$context)

    mat <- tidyr::gather(mat, class, signature, dplyr::contains("Signature"))
    mat <- dplyr::mutate(mat,
                         context = factor(context),
                         base = factor(base, levels = c(
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
        x = components, y = Inf,
        label = ifelse(dist == "norm",
                       paste0(" \u03BC=", signif(mu, digits), ";\u03C3=", signif(sd, digits)),
                       paste0(" \u03BB=", signif(mu, digits))
        )
      ),
      data = params,
      size = params_label_size,
      angle = 90,
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

  p <- p + xlab("Components")

  if (mode == "copynumber") {
    p <- p + ylab("Weights")
  } else if (mode == "mutation") {
    p <- p + ylab("Contributions")
  }
  return(p)
}
