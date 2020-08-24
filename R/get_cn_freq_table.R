#' Get CNV Frequency Table
#'
#' @inheritParams show_cn_freq_circos
#'
#' @return a `data.table`.
#' @export
get_cn_freq_table <- function(data, genome_build = "hg19", cutoff = 2L, resolution_factor = 1L) {
  stopifnot(is.data.frame(data) | inherits(data, "CopyNumber"), resolution_factor >= 1)
  if (is.data.frame(data)) {
    nc_cols <- c("chromosome", "start", "end", "segVal", "sample")
    if (!all(nc_cols %in% colnames(data))) {
      stop("Invalid input, it must contain columns: ", paste(nc_cols, collapse = " "))
    }
  }

  if (inherits(data, "CopyNumber")) {
    genome_build <- data@genome_build
    data <- data@data
  }
  data.table::setDT(data)

  annot <- get_genome_annotation(
    data_type = "cytobands",
    genome_build = genome_build
  )
  annot$start <- annot$start + 1L

  data.table::setDT(annot)
  ## Control the resolution
  if (resolution_factor > 1) {
    f <- function(x, y, n, chrom, band) {
      helper_create_chunks(x, y,
        n = n,
        chrom = chrom,
        band = paste(band, seq_len(n), sep = "-chunk-")
      )
    }
    annot <- purrr::pmap_df(
      data.frame(
        x = annot$start,
        y = annot$end,
        n = resolution_factor,
        chrom = annot$chrom,
        band = annot$band
      ),
      .f = f
    ) %>%
      data.table::as.data.table() %>%
      data.table::setcolorder(c("chrom", "start", "end", "band"))
  }
  data.table::setkey(annot, chrom, start, end)
  merge_dt <- data.table::foverlaps(data, annot,
    by.x = c("chromosome", "start", "end")
  )

  total_samps <- length(unique(data$sample))

  if (length(cutoff) == 1) {
    cutoff <- c(cutoff, cutoff)
  }
  merge_dt$AMP <- merge_dt$segVal > cutoff[2]
  merge_dt$DEL <- merge_dt$segVal < cutoff[1]

  res <- merge_dt[, list(
    freq_AMP = length(unique(sample[AMP])) / total_samps,
    freq_DEL = length(unique(sample[DEL])) / total_samps
  ), by = c("chromosome", "start", "end", "band")]

  colnames(annot)[1] <- "chromosome"

  match_cols <- c("chromosome", "start", "end", "band")
  res <- rbind(res,
    dplyr::setdiff(
      annot[, match_cols, with = FALSE],
      res[, match_cols, with = FALSE]
    ),
    fill = TRUE
  )

  res$freq_AMP <- ifelse(is.na(res$freq_AMP), 0, res$freq_AMP)
  res$freq_DEL <- ifelse(is.na(res$freq_DEL), 0, res$freq_DEL)
  res[order(chromosome, start)]
}
