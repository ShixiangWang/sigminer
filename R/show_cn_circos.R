#' Show Copy Number Profile in Circos
#'
#' Another visualization method for copy number profile like [show_cn_profile].
#'
#' @inheritParams show_cn_profile
#' @param show_title if `TRUE` (default), show title with sample ID.
#' @param col colors for the heatmaps. If it is `NULL`, set to
#' `circlize::colorRamp2(c(1, 2, 4), c("blue", "black", "red"))`.
#' @param samples default is `NULL`, can be a chracter vector representing multiple samples or
#' number of samples to show.
#' If data argument is a `data.frame`, a column called sample must exist.
#' @param side side of the heatmaps.
#' @param ... other parameters passing to [circlize::circos.genomicHeatmap].
#'
#' @return a circos plot
#' @export
#'
#' @examples
#' load(system.file("extdata", "toy_copynumber.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' \donttest{
#' show_cn_circos(cn, samples = 1)
#' show_cn_circos(cn, samples = "TCGA-99-7458-01A-11D-2035-01")
#'
#' ## Remove title
#' show_cn_circos(cn, samples = 1, show_title = FALSE)
#'
#' ## Subset chromosomes
#' show_cn_circos(cn, samples = 1, chrs = c("chr1", "chr2", "chr3"))
#'
#' ## Arrange plots
#' layout(matrix(1:4, 2, 2))
#' show_cn_circos(cn, samples = 4)
#'
#' layout(1) # reset layout
#' }
show_cn_circos <- function(data, samples = NULL,
                           show_title = TRUE,
                           chrs = paste0("chr", 1:22),
                           genome_build = c("hg19", "hg38", "T2T", "mm10", "mm9", "ce11"),
                           col = NULL,
                           side = "inside",
                           ...) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    send_stop("Please install 'circlize' package firstly.")
  }

  if (is.null(col)) {
    col <- circlize::colorRamp2(c(1, 2, 4), c("blue", "black", "red"))
  }

  ## Copy from show_cn_profile: start
  stopifnot(is.data.frame(data) | inherits(data, "CopyNumber"))
  if (is.data.frame(data)) {
    if (is.null(samples)) {
      nc_cols <- c("chromosome", "start", "end", "segVal")
    } else {
      nc_cols <- c("chromosome", "start", "end", "segVal", "sample")
    }
    if (!all(nc_cols %in% colnames(data))) {
      stop("Invalid input, it must contain columns: ", paste(nc_cols, collapse = " "))
    }
  }

  genome_build <- match.arg(genome_build)
  if (inherits(data, "CopyNumber")) {
    genome_build <- data@genome_build
    data <- data@data
  } else {
    data <- data.table::as.data.table(data)
  }

  # Filter data
  if (!is.null(samples)) {
    if (is.character(samples)) {
      data <- data[data$sample %in% samples]
      data$sample <- factor(data$sample, levels = samples)
    } else {
      # Assume it is integer
      samples <- unique(data$sample)[1:samples]
      data <- data[data$sample %in% samples]
      data$sample <- factor(data$sample, levels = samples)
    }
  }

  data$chromosome <- ifelse(startsWith(data$chromosome, prefix = "chr"),
    data$chromosome,
    paste0("chr", data$chromosome)
  )

  data <- data[data$chromosome %in% chrs]
  ## Copy from show_cn_profile: end

  colnames(data)[1] <- "chr"

  if (ncol(data) >= 5L) {
    if (ncol(data) > 5) {
      data = data[, 1:5]
    }
    data_bed <- split(data, data$sample)
    data_bed <- lapply(data_bed, function(x) {
      x$sample <- NULL
      return(x)
    })
    for (i in seq_along(data_bed)) {
      plot_cn_circos(data_bed[[i]],
        species = genome_build, chrs = chrs,
        col = col, side = side,
        title = ifelse(show_title,
          names(data_bed)[i],
          NA
        ),
        ...
      )
    }
  } else {
    plot_cn_circos(data,
      species = genome_build, chrs = chrs,
      col = col, side = side, ...
    )
  }
}


# Plot one sample
plot_cn_circos <- function(data_bed, species, chrs,
                           col, side, title = NA_character_, ...) {
  circlize::circos.initializeWithIdeogram(
    species = species,
    chromosome.index = chrs
  )

  circlize::circos.genomicHeatmap(data_bed,
    col = col,
    side = side,
    ...
  )
  if (!is.na(title)) {
    title(title)
  }
  circlize::circos.clear()
}
