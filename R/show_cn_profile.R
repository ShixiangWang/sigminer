#' Show Sample Copy Number Profile
#'
#' Sometimes it is very useful to check details about copy number profile for one or multiple
#' sample. This function is designed to do this job and can be further modified by **ggplot2**
#' related packages.
#'
#' @param data a [CopyNumber] object or a `data.frame` containing at least 'chromosome', 'start',
#' 'end', 'segVal' these columns.
#' @param samples default is NULL, can be a chracter vector representing multiple samples. If `data` argument
#' is a `data.frame`, a column called `sample` must exist.
#' @param chrs chromosomes start with 'chr'.
#' @param genome_build genome build version, used when `data` is a `data.frame`, , should be 'hg19' or 'hg38'.
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
show_cn_profile = function(data, samples=NULL,
                           chrs = paste0("chr",c(1:22, "X")),
                           genome_build = c("hg19", "hg38"),
                           nrow=NULL, ncol=NULL, return_plotlist=FALSE) {
  stopifnot(is.data.frame(data) | inherits(data, "CopyNumber"))
  if (is.data.frame(data)) {
    if (is.null(samples)) {
      nc_cols = c("chromosome", "start", "end", "segVal")
    } else {
      nc_cols = c("chromosome", "start", "end", "segVal", "sample")
    }
    if (!all(nc_cols %in% data)) {
      stop("Invalid input, it must contain columns: ", paste(nc_cols, collapse = " "))
    }
  }

  genome_build = match.arg(genome_build)
  if (inherits(data, "CopyNumber")) {
    genome_build = data@genome_build
    data = data@data
  }
  data.table::setDT(data)

  # Filter data
  if (!is.null(samples)) {
    data = data[samples %in% data$sample]
  }
  data = data[data$chromosome %in% chrs]

  # Get plot data
  coord_df = build_chrom_coordinate(genome_build, chrs) %>%
    dplyr::mutate(labels = sub("chr", "", .data$chrom))
  merge_df = dplyr::left_join(data, coord_df, by=c("chromosome"="chrom")) %>%
    dplyr::mutate(start = .data$x_start + .data$start - 1,
                  end = .data$x_start + .data$end,
                  segType = dplyr::case_when(
                    .data$segVal > 2 ~ "Amp",
                    .data$segVal == 2 ~ "Normal",
                    .data$segVal < 2 ~ "Del"
                  ))

  plot_df = merge_df %>% dplyr::filter(sample == "TCGA-A8-A07S-01A-11D-A036-01")

  ggplot() +
    geom_segment(aes(x = .data$start, xend = .data$end,
                     y = .data$segVal, yend = .data$segVal,
                     color=.data$segType), data = plot_df) +
    geom_vline(aes(xintercept = .data$x_start), linetype="dotted", data = coord_df) +
    geom_vline(xintercept = coord_df$x_end[nrow(coord_df)], linetype="dotted") +
    scale_x_continuous(breaks = coord_df$lab_loc, labels = coord_df$labels) +
    scale_color_manual(values = c("Amp"="red", "Normal"="black", "Del"="blue"))
}


# Build Chromosome Coordinates
build_chrom_coordinate = function(genome_build, chrs){
  chr_len = get_genome_annotation(data_type = "chr_size",
                                  chrs = chrs,
                                  genome_build = genome_build)
  chr_len$x_end = cumsum(as.numeric(chr_len$size))
  # chr_len$x_start = cumsum(as.numeric(c(0, chr_len$size[-1])))
  chr_len$x_start = c(0, chr_len$x_end[-nrow(chr_len)])
  # Set lab location as middle of chromosome
  chr_len %>%
    dplyr::mutate(lab_loc = chr_len %>%
                    dplyr::rowwise() %>%
                    dplyr::do(lab_loc = mean(c(.$x_start, .$x_end))) %>%
                    dplyr::summarise(lab_loc = round(.data$lab_loc)) %>%
                    dplyr::pull(.data$lab_loc))

}
