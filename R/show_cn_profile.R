#' Show Sample Copy Number Profile
#'
#' Sometimes it is very useful to check details about copy number profile for one or multiple
#' samples. This function is designed to do this job and can be further modified by **ggplot2**
#' related packages.
#'
#' @param data a [CopyNumber] object or a `data.frame` containing at least 'chromosome', 'start',
#' 'end', 'segVal' these columns.
#' @param samples default is NULL, can be a chracter vector representing multiple samples. If `data` argument
#' is a `data.frame`, a column called `sample` must exist.
#' @param show_n number of samples to show, this is used for checking.
#' @param show_title if `TRUE`, show title for multiple samples.
#' @param show_labels one of `NULL`, "s" (for labelling short segments < 1e7)
#' or "a" (all segments).
#' @param chrs chromosomes start with 'chr'.
#' @param position a position range, e.g. `"chr1:3218923-116319008"`. Only data
#' overlaps with this range will be shown.
#' @param genome_build genome build version, used when `data` is a `data.frame`, should be 'hg19' or 'hg38'.
#' @param ylim limites for y axis.
#' @param nrow number of rows in the plot grid when multiple samples are selected.
#' @param ncol number of columns in the plot grid when multiple samples are selected.
#' @param return_plotlist default is `FALSE`, if `TRUE`, return a plot list instead of a combined plot.
#'
#' @return a `ggplot` object or a `list`
#' @export
#'
#' @examples
#' # Load copy number object
#' load(system.file("extdata", "toy_copynumber.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#'
#' p <- show_cn_profile(cn, nrow = 2, ncol = 1)
#' p
#' \donttest{
#' p2 <- show_cn_profile(cn,
#'   nrow = 2, ncol = 1,
#'   position = "chr1:3218923-116319008"
#' )
#' p2
#' }
#' @testexamples
#' expect_s3_class(p, "ggplot")
#' expect_s3_class(p2, "ggplot")
show_cn_profile <- function(data, samples = NULL,
                            show_n = NULL, show_title = FALSE,
                            show_labels = NULL,
                            chrs = paste0("chr", 1:22),
                            position = NULL,
                            genome_build = c("hg19", "hg38", "mm10"),
                            ylim = NULL,
                            nrow = NULL, ncol = NULL,
                            return_plotlist = FALSE) {
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
    data <- data[data$sample %in% samples]
    data$sample <- factor(data$sample, levels = samples)
  }

  data$chromosome <- ifelse(startsWith(data$chromosome, prefix = "chr"),
    data$chromosome,
    paste0("chr", data$chromosome)
  )

  data <- data[data$chromosome %in% chrs]

  if (!is.null(position)) {
    stopifnot(length(position) == 1L)
    position <- unlist(strsplit(position, split = "[:-]"))
    position_df <- data.table::data.table(
      chromosome = position[1],
      start = as.integer(position[2]),
      end = as.integer(position[3])
    )
    data.table::setkey(position_df, chromosome, start, end)
    data_range <- data.table::foverlaps(data, position_df)[!is.na(start)]
    if (nrow(data_range) < 1) {
      message("No data in specified position can be found.")
      return(invisible(NULL))
    } else {
      data <- data_range[, c(-2, -3)]
      colnames(data)[2:3] <- c("start", "end")
      chrs <- position[1]
    }
  }

  # Get plot data
  coord_df <- build_chrom_coordinate(genome_build, chrs) %>%
    dplyr::mutate(labels = sub("chr", "", .data$chrom))
  merge_df <- data %>%
    dplyr::as_tibble() %>%
    dplyr::left_join(coord_df, by = c("chromosome" = "chrom")) %>%
    dplyr::mutate(
      start = .data$x_start + .data$start - 1,
      end = .data$x_start + .data$end,
      segType = dplyr::case_when(
        .data$segVal > 2 ~ "Amp",
        .data$segVal == 2 ~ "Normal",
        .data$segVal < 2 ~ "Del"
      )
    )

  if (!"sample" %in% colnames(merge_df)) {
    # Plot as a single sample
    gg <- plot_cn_profile(merge_df, coord_df, show_labels = show_labels)
    if (!is.null(ylim)) {
      gg <- gg + ggplot2::ylim(ylim)
    }
  } else {
    # Plot mutiple samples
    gg_df <- merge_df %>%
      dplyr::group_by(.data$sample) %>%
      tidyr::nest() %>%
      dplyr::mutate(gg = purrr::map(.data$data,
        plot_cn_profile,
        coord_df = coord_df,
        show_labels = show_labels
      ))

    ## group_by does not maintain sample order
    ## https://github.com/tidyverse/dplyr/issues/3279
    if (!is.null(samples)) {
      gg_df <- gg_df %>%
        dplyr::arrange(.data$sample)
    }

    if (!is.null(show_n)) {
      gg_df <- gg_df %>%
        dplyr::ungroup() %>%
        dplyr::slice(1:show_n)
    }
    gg_list <- gg_df$gg
    if (show_title) {
      gg_list <- purrr::map2(gg_list, gg_df$sample, ~ .x + labs(title = .y))
    }
    if (!is.null(ylim)) {
      gg_list <- purrr::map2(gg_list, gg_df$sample, ~ .x + ylim(ylim))
    }
    names(gg_list) <- gg_df$sample
    if (return_plotlist) {
      return(gg_list)
    }
    gg <- cowplot::plot_grid(
      plotlist = gg_list, align = "hv",
      nrow = nrow, ncol = ncol
    )
  }
  return(gg)
}


# Build Chromosome Coordinates
build_chrom_coordinate <- function(genome_build, chrs) {
  chr_len <- get_genome_annotation(
    data_type = "chr_size",
    chrs = chrs,
    genome_build = genome_build
  )
  chr_len$x_end <- cumsum(as.numeric(chr_len$size))
  # chr_len$x_start = cumsum(as.numeric(c(0, chr_len$size[-1])))
  chr_len$x_start <- c(0, chr_len$x_end[-nrow(chr_len)])
  # Set lab location as middle of chromosome
  chr_len %>%
    dplyr::mutate(lab_loc = chr_len %>%
      dplyr::rowwise() %>%
      dplyr::do(lab_loc = mean(c(.$x_start, .$x_end))) %>%
      dplyr::summarise(lab_loc = round(.data$lab_loc)) %>%
      dplyr::pull(.data$lab_loc))
}

plot_cn_profile <- function(plot_df, coord_df, show_labels = NULL) {
  p <- ggplot() +
    geom_segment(aes(
      x = .data$start, xend = .data$end,
      y = .data$segVal, yend = .data$segVal,
      color = .data$segType
    ), data = plot_df)

  if (!is.null(show_labels)) {
    p <- p +
      ggrepel::geom_label_repel(
        aes(
          x = (.data$start + .data$end) / 2,
          y = .data$segVal + 0.05,
          label = .data$segVal
        ),
        size = 2,
        data = if (show_labels == "s") {
          plot_df %>%
            dplyr::filter(.data$end - .data$start < 1e7)
        } else {
          plot_df
        },
        label.padding = unit(0.1, "lines"),
        point.padding = NA,
        nudge_y = 0.3
      )
  }
  p +
    geom_vline(aes(xintercept = .data$x_start),
      alpha = 0.5,
      linetype = "dotted", data = coord_df
    ) +
    geom_vline(
      xintercept = coord_df$x_end[nrow(coord_df)],
      alpha = 0.5,
      linetype = "dotted"
    ) +
    scale_x_continuous(breaks = coord_df$lab_loc, labels = coord_df$labels) +
    scale_color_manual(values = c("Amp" = "red", "Normal" = "black", "Del" = "blue")) +
    labs(x = "Chromosome", y = "Copy number") +
    cowplot::theme_cowplot() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 60, hjust = 1, size = 9)
    )
}
