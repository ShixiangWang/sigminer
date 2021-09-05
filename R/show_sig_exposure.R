#' Plot Signature Exposure
#'
#' Currently support copy number signatures and mutational signatures.
#' @inheritParams show_sig_profile
#' @param Signature a `Signature` object obtained either from [sig_extract] or [sig_auto_extract],
#' or just a raw **absolute** exposure matrix with column representing samples (patients) and row
#' representing signatures (signature names must end with different digital numbers,
#' e.g. Sig1, Sig10, x12). If you named signatures with letters,
#' you can specify them by `sig_names` parameter.
#' @param sig_names set name of signatures, can be a character vector.
#' @param groups sample groups, default is `NULL`.
#' @param grp_order order of groups, default is `NULL`.
#' @param grp_size font size of groups.
#' @param cutoff a cutoff value to remove hyper-mutated samples.
#' @param palette palette used to plot, default use a built-in palette
#' according to parameter `style`.
#' @param rm_space default is `FALSE`. If `TRUE`, it will remove border color
#' and expand the bar width to 1. This is useful when the sample size is big.
#' @param hide_samps if `TRUE`, hide sample names.
#' @param legend_position position of legend, default is 'top'.
#' @author Shixiang Wang
#' @return a `ggplot` object
#' @import ggplot2
#' @importFrom grDevices rainbow
#' @export
#' @examples
#' # Load mutational signature
#' load(system.file("extdata", "toy_mutational_signature.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Show signature exposure
#' p1 <- show_sig_exposure(sig2)
#' p1
#'
#' # Load copy number signature
#' load(system.file("extdata", "toy_copynumber_signature_by_W.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Show signature exposure
#' p2 <- show_sig_exposure(sig)
#' p2
#' @testexamples
#' expect_s3_class(p1, "ggplot")
#' expect_s3_class(p2, "ggplot")
show_sig_exposure <- function(Signature,
                              sig_names = NULL,
                              groups = NULL,
                              grp_order = NULL,
                              grp_size = NULL,
                              cutoff = NULL,
                              style = c("default", "cosmic"),
                              palette = use_color_style(style),
                              base_size = 12,
                              font_scale = 1,
                              rm_space = FALSE,
                              rm_grid_line = TRUE,
                              rm_panel_border = FALSE,
                              hide_samps = TRUE,
                              legend_position = "top") {
  if (inherits(Signature, "Signature")) {
    h <- Signature$Exposure
  } else if (is.matrix(Signature)) {
    h <- Signature
  } else {
    stop("Invalid input for 'Signature'", call. = FALSE)
  }

  if (!is.null(cutoff)) {
    hyper_index <- which(colSums(h) > cutoff)
    h <- h[, -hyper_index, drop = FALSE]
  }

  style <- match.arg(style)

  h.norm <- apply(h, 2, function(x) x / sum(x))
  h.norm[is.na(h.norm)] <- 0

  if (nrow(h) == 1) {
    h.norm <- as.matrix(h.norm) %>% t()
    rownames(h.norm) <- rownames(h)
  }

  h <- as.data.frame(h)
  h.norm <- as.data.frame(h.norm)

  scale <- font_scale

  .theme_ss <- theme_bw(
    base_size = base_size,
    base_family = "sans"
  ) +
    theme(
      axis.text.x = element_text(
        angle = 90, vjust = 0.5, color = "black",
        hjust = 1, size = (base_size - 2) * scale,
        family = "mono"
      ),
      axis.text.y = element_text(
        hjust = 0.5,
        size = base_size * scale,
        color = "black"
      ),
      strip.text.y = element_text(face = "bold"),
      strip.text.x = element_text(size = grp_size)
    )

  if (style == "cosmic") {
    .theme_ss <- .theme_ss + theme(
      panel.spacing.x = unit(0, "line"),
      strip.background.x = element_rect(color = "white"),
      strip.text.x = element_text(
        color = "white",
        size = grp_size,
        face = "bold"
      )
    )
  }

  if (rm_panel_border) {
    .theme_ss <- .theme_ss + theme(
      panel.border = element_blank()
    )
  }

  if (rm_grid_line) {
    .theme_ss <- .theme_ss + theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  }

  if (is.null(sig_names)) {
    # chop Signature off
    rownames(h.norm) <- rownames(h) <- gsub("[^0-9]", "", rownames(h))
  } else {
    rownames(h.norm) <- rownames(h) <- sig_names
  }


  ordering <- order(colSums(h), decreasing = TRUE)
  h <- h[, ordering]
  h.norm <- h.norm[, ordering]

  sample.ordering <- colnames(h)

  h$Signature <- rownames(h)
  h.norm$Signature <- rownames(h.norm)
  x1 <- tidyr::gather(h, "Sample", "Exposure", -"Signature")
  x2 <- tidyr::gather(h.norm, "Sample", "Exposure", -"Signature")

  if (!is.null(sig_names)) {
    x1$Signature <- factor(x1$Signature, levels = sig_names)
    x2$Signature <- factor(x2$Signature, levels = sig_names)
  } else {
    x1$Signature <- factor(x1$Signature, levels = rownames(h))
    x2$Signature <- factor(x2$Signature, levels = rownames(h))
  }

  x1$class0 <- "Est_Counts"
  x2$class0 <- "Fraction"
  df <- rbind(x1, x2)

  df$class0 <- factor(df$class0, c("Est_Counts", "Fraction"))
  df$Sample <- factor(df$Sample, sample.ordering)

  if (!is.null(groups)) {
    if (is.character(groups)) {
      if (is.null(names(groups))) {
        names(groups) <- sample.ordering
      }
      group_df <- data.frame(
        Sample = names(groups),
        groups = as.character(groups)
      )

      if (!all(df$Sample %in% group_df$Sample)) {
        warning("Not all samples can be found in 'groups'!")
      }

      df <- merge(df, group_df, by = "Sample", all.x = TRUE)

      if (!is.null(grp_order)) {
        df$groups <- factor(df$groups, levels = grp_order)
        df <- df %>%
          dplyr::mutate(
            groups = factor(.data$groups, levels = grp_order)
          )
        # Order by groups
        df <- df %>%
          dplyr::group_by(.data$groups)
      }
    } else {
      stop("groups can only be a character vector!", call. = FALSE)
    }
  }

  p <- ggplot(df, aes_string(x = "Sample", y = "Exposure", fill = "Signature"))
  if (rm_space) {
    p <- p + geom_bar(stat = "identity", position = "stack", alpha = 0.9, width = 1)
  } else {
    p <- p + geom_bar(stat = "identity", position = "stack", color = "black", alpha = 0.9)
  }

  p <- p + scale_fill_manual(values = palette)
  if (!is.null(groups)) {
    p <- p + facet_grid(class0 ~ groups, scales = "free", space = "free_x")
  } else {
    p <- p + facet_grid(class0 ~ ., scales = "free_y")
  }
  p <- p + xlab("Samples") + ylab("Exposure")
  p <- p + .theme_ss
  p <- p + theme(legend.position = legend_position)

  if (hide_samps) {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  }

  if (style == "cosmic") {
    if (!requireNamespace("ggplotify", quietly = TRUE)) {
      send_stop("Please install 'ggplotify' package firstly.")
    }
    g <- ggplot_gtable(ggplot_build(p))

    strip_t <- which(grepl("strip-t", g$layout$name))
    k <- 1
    for (i in strip_t) {
      j <- which(grepl("rect", g$grobs[[i]]$grobs[[1]]$childrenOrder))
      g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- palette[k]
      k <- k + 1
    }
    p <- ggplotify::as.ggplot(g)
  }

  p
}
