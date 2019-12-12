#' Plot Signature Exposure
#'
#' Currently support copy number signatures and mutational signatures.
#' @inheritParams show_sig_profile
#' @param Signature a `Signature` object obtained either from [sig_extract] or [sig_auto_extract],
#' or just a raw exposure matrix with column representing samples (patients) and row
#' representing signatures (row names must start with 'Sig').
#' @param rm_space default is `FALSE`. If `TRUE`, it will remove border color
#' and expand the bar width to 1. This is useful when the sample size is big.
#' @param hide_samps if `TRUE`, hide sample names.
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
#' show_sig_exposure(sig2)
#'
#' # Load copy number signature
#' load(system.file("extdata", "toy_copynumber_signature.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Show signature exposure
#' show_sig_exposure(sig)
show_sig_exposure <- function(Signature,
                              base_size = 12,
                              font_scale = 1,
                              rm_space = FALSE,
                              hide_samps = TRUE) {
  if (class(Signature) == "Signature") {
    h <- Signature$Exposure
  } else if (is.matrix(Signature)) {
    if (!all(startsWith(rownames(Signature), "Sig"))) {
      stop("If Signature is a matrix, column names must start with 'Sig'!", call. = FALSE)
    }
    h <- Signature
  } else {
    stop("Invalid input for 'Signature'", call. = FALSE)
  }

  h.norm <- apply(h, 2, function(x) x / sum(x))
  h <- as.data.frame(h)
  h.norm <- as.data.frame(h.norm)

  scale <- font_scale

  .theme_ss <- theme_bw(base_size = base_size) +
    theme(
      axis.text.x = element_text(
        angle = 90, vjust = 0.5, color = "black",
        hjust = 1, size = (base_size - 2) * scale
      ),
      axis.text.y = element_text(
        hjust = 0.5,
        size = base_size * scale,
        color = "black"
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  # chop Signature off
  rownames(h) <- sub(".*[^\\d+](\\d+)$", "\\1", rownames(h))
  rownames(h.norm) <- sub(".*[^\\d+](\\d+)$", "\\1", rownames(h.norm))

  ordering <- order(colSums(h), decreasing = TRUE)
  h <- h[, ordering]
  h.norm <- h.norm[, ordering]

  sample.ordering <- colnames(h)

  h$Signature <- rownames(h)
  h.norm$Signature <- rownames(h.norm)
  x1 <- tidyr::gather(h, "Sample", "Exposure", -"Signature")
  x2 <- tidyr::gather(h.norm, "Sample", "Exposure", -"Signature")

  x1$class0 <- "Contribution"
  x2$class0 <- "Fraction"
  df <- rbind(x1, x2)

  df$class0 <- factor(df$class0, c("Contribution", "Fraction"))
  df$Sample <- factor(df$Sample, sample.ordering)

  p <- ggplot(df, aes_string(x = "Sample", y = "Exposure", fill = "Signature"))
  if (rm_space) {
    p <- p + geom_bar(stat = "identity", position = "stack", alpha = 0.9, width = 1)
  } else {
    p <- p + geom_bar(stat = "identity", position = "stack", color = "black", alpha = 0.9)
  }

  p <- p + scale_fill_manual(values = c("red", "cyan", "yellow", "blue", "magenta", "gray50", "orange", "darkgreen", "brown", "black", rainbow(10)[4:10]))
  p <- p + facet_grid(class0 ~ ., scales = "free_y")
  p <- p + xlab("Samples") + ylab("Exposure")
  p <- p + .theme_ss
  p <- p + theme(legend.position = "top")

  if (hide_samps) {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  }
  p
}
