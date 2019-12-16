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
#' @param method method for copy number feature classfication, can be one of "Macintyre" ("M") and
#' "Wang" ("W").
#' @param normalize one of 'row', 'column', 'raw' and "feature", for row normalization (signature),
#' column normalization (component), raw data, row normalization by feature, respectively.
#' Of note, 'feature' only works when the mode is 'copynumber'.
#' @param set_gradient_color default is `FALSE`, if `TRUE`, use gradient colors
#' to fill bars. **This is very useful when signatures are extracted from "Macintyre" method.**
#' @param x_label_angle font angle for x label.
#' @param params params `data.frame` of components, obtained from [sig_derive].
#' @param show_cv default is `FALSE`, if `TRUE`, show coefficient of variation when
#' `params` is not `NULL`.
#' @param params_label_size font size for params label.
#' @param params_label_angle font angle for params label.
#' @param y_expand y expand height for plotting params of copy number signatures.
#' @param digits digits for plotting params of copy number signatures.
#' @param font_scale a number used to set font scale.
#' @param sig_names set name of signatures, can be a character vector.
#' Default is `NULL`, prefix 'Sig_' plus number is used.
#' @param sig_orders set order of signatures, can be a character vector.
#' Default is `NULL`, the signatures are ordered by alphabetical order.
#' @param check_sig_names if `TRUE`, check signature names when input is
#' a matrix, i.e., all signatures (colnames) must start with 'Sig'.
#' @author Shixiang Wang
#' @inheritParams sig_derive
#' @return a `ggplot` object
#' @export
#' @examples
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
#' # =================
#' # Load copy number prepare object
#' load(system.file("extdata", "toy_copynumber_prepare.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' params <- get_tidy_parameter(cn_prepare$components)
#' show_sig_profile(sig, params = params, y_expand = 2)
show_sig_profile <- function(Signature, mode = c("copynumber", "mutation"),
                             method = "Macintyre",
                             normalize = c("row", "column", "raw", "feature"),
                             feature_setting = sigminer::CN.features,
                             set_gradient_color = FALSE,
                             x_label_angle = 60,
                             params = NULL, show_cv = FALSE,
                             params_label_size = 3,
                             params_label_angle = 60, y_expand = 1,
                             digits = 2, base_size = 12, font_scale = 1,
                             sig_names = NULL, sig_orders = NULL,
                             check_sig_names = TRUE) {
  if (class(Signature) == "Signature") {
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
  method <- match.arg(method, choices = c("Macintyre", "M", "Wang", "W"))
  normalize <- match.arg(normalize)

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
        angle = x_label_angle, vjust = 1,
        hjust = 1, size = (base_size - 4) * scale,
        color = "black"
      ),
      axis.text.y = element_text(
        hjust = 0.5,
        size = base_size * scale,
        color = "black"
      )
    )
  # <<<<<<<<<<<<<<<<< Setting theme

  # >>>>>>>>>>>>>>>>> identify mode and do data transformation
  mat <- as.data.frame(Sig)
  mat$context <- rownames(mat)

  if (mode == "copynumber") {
    if (startsWith(method, "M")) {
      mat$base <- sub("\\d+$", "", mat$context)
      mat <- tidyr::gather(mat, class, signature, -c("context", "base"))

      mat <- mat %>%
        dplyr::mutate(
          base = factor(.data$base, levels = c(
            "bp10MB", "bpchrarm",
            "copynumber", "changepoint",
            "osCN", "segsize"
          ))
        ) %>%
        dplyr::arrange(.data$base)

      mat <- dplyr::mutate(mat,
        context = factor(.data$context,
          levels = unique(mat[["context"]])
        ),
        class = factor(class)
      )
    } else {
      mat$base <- sub("\\[.*\\]$", "", mat$context)
      mat <- tidyr::gather(mat, class, signature, -c("context", "base"))

      if (!inherits(feature_setting, "sigminer.features")) {
        feature_setting <- get_feature_components(feature_setting)
      }
      avail_features <- unique(feature_setting$feature)

      # Make sure
      f_orders <- avail_features[avail_features %in% mat$base]
      comp_orders <- feature_setting$component[feature_setting$feature %in% f_orders]

      mat <- dplyr::mutate(mat,
        context = factor(.data$context,
          levels = comp_orders
        ),
        base = factor(.data$base, levels = f_orders),
        class = factor(class)
      )
    }
  } else {
    mat$base <- sub("[ACGT]\\[(.*)\\][ACGT]", "\\1", mat$context)
    mat$context <- sub("(\\[.*\\])", "\\-", mat$context)

    mat <- tidyr::gather(mat, class, signature, -c("context", "base"))
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

  if (normalize == "feature") {
    if (mode != "copynumber") {
      message("normalize='feature' only supports when the argument mode is 'copynumber'. Please use another value!")
      return(invisible())
    } else {
      mat <- mat %>%
        dplyr::group_by(.data$base, .data$class) %>%
        tidyr::nest() %>%
        dplyr::mutate(
          context = purrr::map(.data$data, ~ .$context),
          signature = purrr::map(.data$data, ~ .$signature / sum(.$signature))
        ) %>%
        dplyr::select(-"data") %>%
        tidyr::unnest(cols = c("context", "signature")) %>%
        dplyr::ungroup()
    }
  }

  # >>>>>>> Set signature name and order
  if (!is.null(sig_names)) {
    if (length(sig_names) != length(unique(mat[["class"]]))) {
      stop("The length of input signature names is not equal to signature number")
    }
    names(sig_names) <- paste0("Sig_", seq_along(sig_names))
    mat[["class"]] <- sig_names[mat[["class"]]]
  }

  if (!is.null(sig_orders)) {
    mat[["class"]] <- factor(mat[["class"]], levels = sig_orders)
  } else {
    mat[["class"]] <- factor(mat[["class"]])
  }
  # >>>>>>>>>>>>>>>>>>>>>>> Plot

  if (set_gradient_color) {
    if (mode == "mutation") {
      message("Gradient colors are not supported for mutational signatures!")
      message("Please set it to FALSE and re-run.")
      return(invisible())
    }

    col_df <- mat %>%
      dplyr::filter(.data$class == .data$class[1]) %>%
      dplyr::group_by(.data$base) %>%
      dplyr::summarise(N = dplyr::n())

    p <- ggplot(mat) +
      geom_bar(aes_string(x = "context", y = "signature", fill = "context"),
        stat = "identity", position = "identity", colour = "gray50"
      ) +
      scale_fill_manual(values = helper_create_colormap(col_df$base, col_df$N))
  } else {
    p <- ggplot(mat) +
      geom_bar(aes_string(x = "context", y = "signature", fill = "base"),
        stat = "identity", position = "identity", colour = "gray50"
      ) +
      scale_fill_manual(values = c(
        "cyan", "red", "yellow", "purple",
        "green", "blue", "black", "gray"
      ))
  }


  if (mode == "copynumber" & startsWith(method, "M")) {
    if (!is.null(params)) {
      params$class <- factor(levels(mat[["class"]])[1], levels = levels(mat[["class"]]))
      p <- p + geom_text(aes(
        x = .data$components, y = Inf,
        label = ifelse(.data$dist == "norm",
          paste0(
            " \u03BC=", signif(.data$mean, digits),
            ifelse(rep(show_cv, length(.data$cv)), paste0("; cv=", signif(.data$cv, digits)), "")
          ),
          paste0(
            " \u03BB=", signif(.data$mean, digits),
            ifelse(rep(show_cv, length(.data$cv)), paste0("; cv=", signif(.data$cv, digits)), "")
          )
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

  # Remove prefix to keep space
  if (startsWith(method, "W")) {
    p <- p + scale_x_discrete(
      breaks = mat$context,
      labels = sub(".*(\\[.*\\])$", "\\1", mat$context)
    )
  }

  p <- p +
    guides(fill = FALSE) + .theme_ss

  if (all(mode == "copynumber", !is.null(params))) {
    p <- p + theme(plot.margin = margin(30 * y_expand, 2, 2, 2, unit = "pt")) # Add regions
  }

  p <- p + xlab("Components")

  if (normalize == "column") {
    p <- p + ylab("Weights")
  } else if (normalize == "row" | normalize == "feature") {
    p <- p + ylab("Contributions")
  } else {
    if (mode == "copynumber") {
      p <- p + ylab("Estimated number of events")
    } else {
      p <- p + ylab("Estimated mutation counts")
    }
  }

  return(p)
}
