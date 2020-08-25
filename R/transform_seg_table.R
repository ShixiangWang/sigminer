#' Transform Copy Number Table
#'
#' @inheritParams get_cn_freq_table
#' @inheritParams tidyr::pivot_wider
#' @param ref_type annotation data type used for constructing matrix.
#'
#' @return a `data.table`.
#' @export
#'
#' @examples
#' load(system.file("extdata", "toy_copynumber.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Compute the mean segVal in each cytoband
#' x <- transform_seg_table(cn, resolution_factor = 1)
#' x
#' # Compute the mean segVal in each half-cytoband
#' x2 <- transform_seg_table(cn, resolution_factor = 2)
#' x2
#' @testexamples
#' expect_is(x, "data.table")
#' expect_is(x2, "data.table")
transform_seg_table <- function(data,
                                genome_build = c("hg19", "hg38", "mm10"),
                                ref_type = c("cytoband", "gene"),
                                values_fill = NA,
                                values_fn = function(x, ...) {
                                  round(mean(x, ...))},
                                resolution_factor = 1L) {

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

  ref_type <- match.arg(ref_type)

  #data$sample <- factor(data$sample, levels = unique(data$sample))
  data$chromosome <- ifelse(startsWith(data$chromosome, prefix = "chr"),
                            data$chromosome,
                            paste0("chr", data$chromosome))

  if (ref_type == "cytoband") {
    annot <- get_genome_annotation(
      data_type = "cytobands",
      genome_build = genome_build
    )
    annot$start <- annot$start + 1L
  } else {
    if (genome_build == "mm10") {
      # Not support for now
      annot_file <- system.file("extdata", "mouse_mm10_gene_info.rds",
                                package = "sigminer", mustWork = TRUE)
    } else {
      annot_file <- system.file("extdata", paste0("human_", genome_build, "_gene_info.rds"),
                                package = "sigminer", mustWork = TRUE)
    }

    annot <- readRDS(annot_file)
    annot <- annot[, c("chrom", "start", "end", "gene_name", "gene_type")]
    colnames(annot)[4] <- "band"
  }


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
  merge_dt <- merge_dt %>%
    dplyr::as_tibble() %>%
    dplyr::select(-c("i.start", "i.end")) %>%
    na.omit() %>%
    tidyr::pivot_wider(names_from = "sample", values_from = "segVal",
                       values_fill = values_fill, values_fn = values_fn)
  colnames(merge_dt)[4] <- "label"
  merge_dt %>% data.table::as.data.table()
}
