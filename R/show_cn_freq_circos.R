#' Show Copy Number Variation Frequency Profile with Circos
#'
#' @inheritParams show_cn_circos
#' @param data a `CopyNumber` object or a data.frame containing
#' at least 'chromosome', 'start', 'end', 'segVal', 'sample' these columns.
#' @param groups a named list or a column name for specifying groups.
#' @param cutoff copy number value cutoff for splitting data into AMP and DEL.
#' The values equal to cutoff are discarded. Default is `2`, you can also set
#' a length-2 vector, e.g. `c(2, 2)`.
#' @param resolution_factor an integer to control the resolution.
#' When it is `1` (default), compute frequency in each cytoband.
#' When it is `2`, use compute frequency in each half cytoband.
#' @param title length-2 titles for AMP and DEL.
#' @param cols length-2 colors for AMP and DEL.
#' @param plot_ideogram default is `TRUE`, show ideogram.
#' @param track_height track height in `mm` unit.
#' @param ideogram_height ideogram height in `mm` unit.
#' @param ... other parameters passing to [circlize::circos.genomicLines].
#'
#' @return Nothing.
#' @export
#'
#' @examples
#' \donttest{
#' load(system.file("extdata", "toy_copynumber.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#'
#' show_cn_freq_circos(cn)
#' ss <- unique(cn@data$sample)
#' show_cn_freq_circos(cn, groups = list(a = ss[1:5], b = ss[6:10]), cols = c("red", "green"))
#' }
show_cn_freq_circos <- function(data,
                                groups = NULL,
                                cutoff = 2L,
                                resolution_factor = 1L,
                                title = c("AMP", "DEL"),
                                chrs = paste0("chr", 1:22),
                                genome_build = c("hg19", "hg38", "mm10", "mm9"),
                                cols = NULL,
                                plot_ideogram = TRUE,
                                track_height = 0.5,
                                ideogram_height = 1,
                                ...) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    send_stop("Please install 'circlize' package firstly.")
  }

  if (is.null(cols)) {
    cols <- c("#FF000080", "#0000FF80")
  }

  stopifnot(is.data.frame(data) | inherits(data, "CopyNumber"))
  if (is.data.frame(data)) {
    nc_cols <- c("chromosome", "start", "end", "segVal", "sample")
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

  # The 'Groups' can be either a list of name index
  # or the column name
  if (is.list(groups)) {
    if (length(purrr::reduce(groups, intersect)) != 0) {
      stop("Duplicated sample names in group list!", call. = FALSE)
    }
    samples <- purrr::reduce(groups, c)
    data <- data[data$sample %in% samples]
    grp_df <- tibble::enframe(groups) %>%
      tidyr::unnest("value") %>%
      purrr::set_names(c("grp_name", "sample")) %>%
      data.table::as.data.table()

    data <- merge(data, grp_df, by = "sample")
  } else if (is.character(groups)) {
    if (length(groups) != 1) {
      stop("When 'groups' is character type, it represents one column name.",
        call. = FALSE
      )
    }
    samples <- unique(data$sample)
    colnames(data)[colnames(data) == groups] <- "grp_name"
  } else {
    samples <- unique(data$sample)
    data$grp_name <- "1"
  }

  data$sample <- factor(data$sample, levels = samples)

  data$chromosome <- ifelse(startsWith(data$chromosome, prefix = "chr"),
    data$chromosome,
    paste0("chr", data$chromosome)
  )

  data <- data[data$chromosome %in% chrs]

  ## Construct CNV frequency data
  ## for AMP/DEL and each group
  data_freq <- data[, get_cn_freq_table(.SD,
    genome_build = genome_build,
    cutoff = cutoff,
    resolution_factor = resolution_factor
  ), by = "grp_name"]
  data.table::setcolorder(data_freq, c(
    "chromosome", "start", "end",
    "freq_AMP", "freq_DEL"
  ))

  plot_cn_circos_density(split(data_freq, data_freq$grp_name),
    species = genome_build, chrs = chrs,
    cols = cols,
    title = title,
    plot_ideogram = plot_ideogram,
    track_height = track_height,
    ideogram_height = ideogram_height,
    ...
  )
}


# Plot density circles
plot_cn_circos_density <- function(bed_list, species, chrs,
                                   cols = c("#FF000080", "#0000FF80"),
                                   title = c("AMP", "DEL"),
                                   plot_ideogram = TRUE,
                                   track_height = 0.5,
                                   ideogram_height = 2,
                                   ...) {
  layout(matrix(1:2, 1, 2))
  on.exit(layout(1))
  on.exit(circlize::circos.clear())
  # on.exit(circlize::circos.par(RESET = TRUE))
  force(cols)

  bed_list1 <- lapply(bed_list, function(x) x[, c(1:4)])
  bed_list2 <- lapply(bed_list, function(x) x[, c(1:3, 5)])

  gap_after <- c(rep(2, length(chrs) - 1), 10)
  circlize::circos.par(gap.after = gap_after)

  message("Plotting AMP")
  circlize::circos.initializeWithIdeogram(
    species = species,
    chromosome.index = chrs,
    plotType = if (plot_ideogram) c("ideogram", "axis", "labels") else c("axis", "labels"),
    track.height = circlize::convert_height(track_height, "mm"),
    ideogram.height = circlize::convert_height(ideogram_height, "mm")
  )
  j <- 1
  for (bed in bed_list1) {
    message("Plotting group ", names(bed_list1)[j])
    circlize::circos.genomicTrack(
      bed,
      ylim = c(0, 1),
      panel.fun = function(region, value, ...) {
        circlize::circos.genomicLines(region,
          value,
          col = cols[1],
          area = TRUE, ...
        )
      }
    )
    if (j == 1) {
      circlize::circos.yaxis(
        labels.cex = 0.4,
        tick.length = circlize::convert_x(
          0.3, "mm",
          circlize::get.cell.meta.data("sector.index"),
          circlize::get.cell.meta.data("track.index")
        ),
        side = "right", sector.index = chrs[length(chrs)]
      )
    }
    j <- j + 1
  }
  if (!is.null(title)) {
    text(0, 0, title[1])
  }
  circlize::circos.clear()

  circlize::circos.par(gap.after = gap_after)

  message("Plotting DEL")
  circlize::circos.initializeWithIdeogram(
    species = species,
    chromosome.index = chrs,
    plotType = if (plot_ideogram) c("ideogram", "axis", "labels") else c("axis", "labels"),
    track.height = circlize::convert_height(track_height, "mm"),
    ideogram.height = circlize::convert_height(ideogram_height, "mm")
  )
  j <- 1
  for (bed in bed_list2) {
    message("Plotting group ", names(bed_list2)[j])
    circlize::circos.genomicTrack(
      bed,
      ylim = c(0, 1),
      panel.fun = function(region, value, ...) {
        circlize::circos.genomicLines(region,
          value,
          col = cols[2],
          area = TRUE, ...
        )
      }
    )
    if (j == 1) {
      circlize::circos.yaxis(
        labels.cex = 0.4,
        tick.length = circlize::convert_x(
          0.3, "mm",
          circlize::get.cell.meta.data("sector.index"),
          circlize::get.cell.meta.data("track.index")
        ),
        side = "right", sector.index = chrs[length(chrs)]
      )
    }
    j <- j + 1
  }
  if (!is.null(title)) {
    text(0, 0, title[2])
  }
  circlize::circos.clear()
}
