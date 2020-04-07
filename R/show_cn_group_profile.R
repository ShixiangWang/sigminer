# Show a summary copy number profile for a group of samples
#
# TODO/NOTE:
# Implement it with foverlap and cytoband
# 对每个 cytoband 按密度设定进行分割
# 然后将病人的 CNV map 上去，计算 CNV 扩增或删除的比例
# 然后绘图

show_cn_group_profile <- function(data,
                                  chrs = paste0("chr", c(1:22, "X")),
                                  genome_build = c("hg19", "hg38")) {
  stopifnot(is.data.frame(data) | inherits(data, "CopyNumber"))
  genome_build <- match.arg(genome_build)

  if (is.data.frame(data)) {
    nc_cols <- c("chromosome", "start", "end", "segVal")
    if (!all(nc_cols %in% data)) {
      stop("Invalid input, it must contain columns: ", paste(nc_cols, collapse = " "))
    }
  } else {
    genome_build <- data@genome_build
    data <- data@data
  }

  # Keep CNV data
  data_cnv <- dplyr::filter(data, .data$segVal != 2)

  # Construct a GRanges object for all regions
  data_GR_total <- GenomicRanges::GRanges(
    seqnames = data$chromosome,
    ranges = IRanges::IRanges(start = data$start, end = data$end)
  ) %>%
    GenomicRanges::reduce()

  # Construct a GRanges object for CNV region
  data_GR_cnv <- GenomicRanges::GRanges(
    seqnames = data_cnv$chromosome,
    ranges = IRanges::IRanges(start = data_cnv$start, end = data_cnv$end)
  ) %>%
    GenomicRanges::reduce()

  # Get normal region
  data_GR_normal <- GenomicRanges::setdiff(data_GR_total, data_GR_cnv)
  data_normal <- data_GR_normal %>%
    as.data.frame() %>%
    dplyr::mutate(
      chromosome = as.character(.data$seqnames),
      segVal = 2
    ) %>%
    dplyr::select(-c("seqnames", "width", "strand"))

  # Construct regions
  # Example:
  # For a chromosome, two regions 1-100 and 20-120
  # split them into three: 1-20, 21-100, 101-120
  data_cnv_split <- purrr::map_df(unique(data_cnv$chromosome),
    function(x, df) {
      df <- dplyr::filter(df, .data$chromosome == x)
      loc <- sort(unique(c(df$start, df$end)))
      start <- loc[-length(loc)]
      start[-1] <- start[-1] + 1
      dplyr::tibble(
        chromosome = x,
        start = start,
        end = loc[-1]
      )
    },
    df = data_cnv
  )

  # data_join <- data.table::foverlaps(
  #   data_cnv_split, data_cnv,
  #   by = c("chromosome", "start", "end")
  # ) %>%
  #   dplyr::mutate(
  #     ID = paste(.data$chromosome.x, .data$start.x, .data$end.x, sep = "-"),
  #     segType = ifelse(.data$segVal > 2, "Amp", "Del")
  #   )
  data_join <- data_join %>%
    dplyr::group_by(.data$ID, .data$segType) %>%
    dplyr::summarise(segVal = mean(.data$segVal)) %>%
    tidyr::separate(col = "ID", into = c("chromosome", "start", "end")) %>%
    dplyr::mutate(
      start = as.numeric(.data$start),
      end = as.numeric(.data$end)
    ) %>%
    dplyr::select(-.data$segType) %>%
    dplyr::bind_rows(data_normal)

  show_cn_profile(data_join, chrs = chrs, genome_build = genome_build, .call = TRUE)
}
