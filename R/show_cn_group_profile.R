show_cn_group_profile = function(data,
                                 chrs = paste0("chr",c(1:22, "X")),
                                 genome_build = c("hg19", "hg38")) {
  stopifnot(is.data.frame(data) | inherits(data, "CopyNumber"))
  genome_build = match.arg(genome_build)

  if (is.data.frame(data)) {
    nc_cols = c("chromosome", "start", "end", "segVal")
    if (!all(nc_cols %in% data)) {
      stop("Invalid input, it must contain columns: ", paste(nc_cols, collapse = " "))
    }
  } else {
    genome_build = data@genome_build
    data = data@data
  }

  data_cnv = dplyr::filter(data, .data$segVal != 2)

  # Construct a GRanges object
  data_GR = GenomicRanges::GRanges(seqnames = data_cnv$chromosome,
                    ranges = IRanges::IRanges(start = data_cnv$start, end = data_cnv$end),
                    segVal = data_cnv$segVal)
  data_cnv_merge = IRanges::reduce(data_GR) %>% as.data.frame()
  data_cnv_merge %>%
    dplyr::group_by(.data$seqnames) %>%
    dplyr::do(
      locs = sort(unique(c(.$start, .$end)))
    ) -> tt

}
