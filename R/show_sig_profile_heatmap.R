#' Show Signature Profile with Heatmap
#'
#' This is a complementary function to [show_sig_profile()], it is used for visualizing
#' some big signatures, i.e. SBS-1536, not all signatures are supported. See details for
#' current supported signatures.
#'
#' Support:
#'
#' - SBS-24
#' - SBS-96
#' - SBS-384
#' - SBS-1536
#' - SBS-6144
#' - DBS-78
#' - DBS-186
#'
#' @inheritParams show_sig_profile
#' @param mode one of "SBS" and "DBS".
#' @param x_lab x label.
#' @param y_lab y label.
#' @param legend_name name of figure legend.
#' @param palette color for value.
#' @param x_label_angle angle for x axis text.
#' @param x_label_vjust vjust for x axis text.
#' @param x_label_hjust hjust for x axis text.
#' @param y_label_angle angle for y axis text.
#' @param y_label_vjust vjust for y axis text.
#' @param y_label_hjust hjust for y axis text.
#' @param flip_xy if `TRUE`, flip x axis and y axis.
#
#' @return a `ggplot` object.
#' @export
#'
#' @examples
#' \donttest{
#' # Load SBS signature
#' load(system.file("extdata", "toy_mutational_signature.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Show signature profile
#' p1 <- show_sig_profile_heatmap(sig2, mode = "SBS")
#' p1
#' }
#' @testexamples
#' expect_s3_class(p1, "ggplot")
show_sig_profile_heatmap <- function(Signature, mode = c("SBS", "DBS"),
                                     normalize = c("row", "column", "raw"),
                                     filters = NULL,
                                     x_lab = NULL,
                                     y_lab = NULL,
                                     legend_name = "auto",
                                     palette = "red",
                                     x_label_angle = 90,
                                     x_label_vjust = 1,
                                     x_label_hjust = 0.5,
                                     y_label_angle = 0,
                                     y_label_vjust = 0.5,
                                     y_label_hjust = 1,
                                     flip_xy = FALSE,
                                     sig_names = NULL, sig_orders = NULL,
                                     check_sig_names = TRUE) {
  stopifnot(all(!is.na(sig_names)))

  if (inherits(Signature, "Signature")) {
    Sig <- Signature$Signature
  } else if (is.matrix(Signature)) {
    if (check_sig_names) {
      if (!all(startsWith(colnames(Signature), "Sig"))) {
        stop("If Signature is a matrix, column names must start with 'Sig'!", call. = FALSE)
      }
    }
    Sig <- Signature
  } else {
    stop("Invalid input for 'Signature'", call. = FALSE)
  }

  mode <- match.arg(mode)
  ## TODO: only support "Tao & Wang" for mode "copynumber"
  # method <- match.arg(method, choices = c("Wang", "W", "Tao & Wang", "T"))
  normalize <- match.arg(normalize)

  if (normalize == "row") {
    Sig <- apply(Sig, 2, function(x) x / sum(x))
  } else if (normalize == "column") {
    Sig <- t(apply(Sig, 1, function(x) x / sum(x)))
  }

  if (!is.null(filters)) {
    Sig <- Sig[sapply(rownames(Sig), function(x) {
      any(sapply(filters, function(y) {
        grepl(y, x)
      }))
    }), ]
  }

  ## Check if it has any transcription labels
  has_labels <- grepl("^[TUBNQ]:", rownames(Sig))
  if (any(has_labels)) {
    has_labels <- TRUE
  } else {
    has_labels <- FALSE
  }

  if (is.numeric(sig_orders)) {
    Sig <- Sig[, sig_orders, drop = FALSE]
  }

  # >>>>>>>>>>>>>>>>> identify mode and do data transformation
  mat <- as.data.frame(Sig)
  mat$context <- rownames(mat)

  if (has_labels) {
    mat$context <- substring(mat$context, 3)
    mat$label <- factor(substr(rownames(mat), 1, 1))
  }

  if (mode == "copynumber") {
    ## Tao & Wang
    mat$base <- sub("^([A-Z]:[A-Z]):[0-9]\\+?$", "\\1", mat$context)
    mat$context <- sub("^[A-Z]:[A-Z]:([0-9]\\+?)$", "\\1", mat$context)
    mat <- tidyr::gather(mat, class, signature, -c("context", "base"))

    mat <- dplyr::mutate(mat,
      context = factor(.data$context),
      base = factor(.data$base, levels = c(
        "O:S", "O:M", "O:L",
        "R:S", "R:M", "R:L",
        "F:S", "F:M", "F:L",
        "N:S", "N:M", "N:L"
      )),
      class = factor(class)
    )
  } else if (mode == "SBS") {
    if (any(grepl("^[ACGT]>[ACGT]$", rownames(Sig)))) {
      stop("Unsupported 6 compoents!", call. = FALSE)
    }
    mat$base <- sub("^[ACGT]{0,2}\\[([ACGT]>[ACGT])\\][ACGT]{0,2}", "\\1", mat$context)
    if (!any(grepl("C>A", mat$base))) {
      send_stop("You may choose wrong 'mode' option, this is designed for 'SBS'.")
    }
    mat$context <- sub("(\\[.*\\])", "\\-", mat$context)

    if (any(grepl(">", mat$context))) {
      mat$context <- mat$label
    } else if (all(nchar(mat$context) == 5L)) {
      if (has_labels) {
        # T:NN[N>N]NN
        mat$context2 <- substr(mat$context, 2L, 4L)
        mat$context <- paste0(
          mat$label, ":",
          substr(mat$context, 1, 1),
          "---",
          substr(mat$context, 5, 5)
        )
      } else {
        # NN[N>N]NN
        mat$context2 <- substr(mat$context, 2L, 4L)
        mat$context <- paste0(
          substr(mat$context, 1, 1),
          "---",
          substr(mat$context, 5, 5)
        )
      }
    } else {
      if (has_labels) {
        mat$context2 <- mat$context
        mat$context <- mat$label
      }
    }
    mat$label <- NULL

    if (!is.null(mat$context2)) {
      mat <- tidyr::gather(mat, class, signature, -c("context", "base", "context2"))
    } else {
      mat <- tidyr::gather(mat, class, signature, -c("context", "base"))
    }

    mat <- dplyr::mutate(mat,
      context = factor(.data$context),
      base = factor(.data$base, levels = c(
        "C>A", "C>G",
        "C>T", "T>A",
        "T>C", "T>G"
      )),
      class = factor(class)
    )
  } else if (mode == "DBS") {
    if (any(nchar(mat$context) <= 5)) {
      mat$base <- paste0(substr(mat$context, 1, 3), "NN")
      mat$context <- substr(mat$context, 4, 5)
    } else {
      mat$base <- paste0(substr(mat$context, 3, 5), "NN")
      mat$context2 <- substr(mat$context, 6, 7)
      mat$context <- paste0(
        substr(mat$context, 1, 1),
        "--",
        substr(mat$context, 9, 9)
      )
    }
    if (!any(grepl("[ACGT]{2}>[NN]{2}", mat$base))) {
      send_stop("You may choose wrong 'mode' option, this is designed for 'DBS'.")
    }


    mat$base <- substr(mat$base, 1, 2)
    if (has_labels) {
      ## T:NN>NN
      mat$base <- paste0(mat$label, ":", substr(mat$base, 1, 2))
    }
    mat$label <- NULL

    if (!is.null(mat$context2)) {
      mat <- tidyr::gather(mat, class, signature, -c("context", "base", "context2"))
    } else {
      mat <- tidyr::gather(mat, class, signature, -c("context", "base"))
    }


    mat <- dplyr::mutate(mat,
      context = factor(.data$context),
      base = factor(.data$base),
      class = factor(class)
    )
  }

  # >>>>>>> Set signature name and order
  if (!is.null(sig_names)) {
    if (length(sig_names) != length(unique(mat[["class"]]))) {
      stop("The length of input signature names is not equal to signature number!")
    }
    names(sig_names) <- paste0("Sig_", seq_along(sig_names))
    mat[["class"]] <- sig_names[mat[["class"]]]
  }

  if (!is.null(sig_orders) & !is.numeric(sig_orders)) {
    mat[["class"]] <- factor(mat[["class"]], levels = sig_orders)
  } else {
    mat[["class"]] <- factor(mat[["class"]])
  }
  # >>>>>>>>>>>>>>>>>>>>>>> Plot

  ## TODO: 处理转录链
  if (is.null(mat$context2)) {
    p <- ggplot(mat, aes_string(x = "context", y = "base", fill = "signature")) +
      geom_tile(color = "black")

    p <- p + facet_grid(class ~ ., scales = "free")
  } else {
    p <- ggplot(mat, aes_string(x = "context2", y = "context", fill = "signature")) +
      geom_tile(color = "black")

    p <- p + facet_grid(class ~ base, scales = "free")
  }


  if (is.character(palette)) {
    p <- p + scale_fill_gradient(low = "white", high = palette)
  } else if (inherits(palette, "Scale")) {
    p <- p + palette
  } else {
    stop("Unvalid input for 'palette'!", call. = FALSE)
  }

  p <- p + xlab(x_lab) + ylab(y_lab)

  if (legend_name == "auto") {
    if (normalize == "column") {
      p <- p + labs(fill = "Weights")
    } else if (normalize == "row" | normalize == "feature") {
      p <- p + labs(fill = "Contributions")
    } else {
      p <- p + labs(fill = "Estimated counts")
    }
  } else {
    p <- p + labs(fill = legend_name)
  }

  .theme_ss <- theme_minimal(
    base_family = "sans"
  ) +
    theme(
      axis.text.x = element_text(
        angle = x_label_angle,
        vjust = x_label_vjust,
        hjust = x_label_hjust,
        color = "black",
        family = "mono"
      ),
      axis.text.y = element_text(
        angle = y_label_angle,
        vjust = y_label_vjust,
        hjust = y_label_hjust,
        color = "black",
        family = "mono"
      ),
      legend.position = "bottom",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  p <- p + .theme_ss

  if (flip_xy) {
    p <- p + coord_flip()
  }

  return(p)
}
