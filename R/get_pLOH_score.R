#' Get proportions of pLOH score from Allele Specific Copy Number Profile
#'
#' pLOH score represents the genome that displayed LOH.
#'
#' @inheritParams read_copynumber
#' @param data a CopyNumber object or a `data.frame` containing at least
#' 'chromosome', 'start', 'end', 'segVal', "minor_cn", 'sample' these columns.
#' @param rm_chrs chromosomes to be removed in calculation. Default is sex
#' chromosomes (recommended).
#' @references
#' Steele, Christopher D., et al. "Signatures of copy number alterations in human cancer." bioRxiv (2021).
#'
#' @return A `data.frame`
#' @export
#'
#' @examples
#' # Load toy dataset of absolute copynumber profile
#' load(system.file("extdata", "toy_segTab.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#'
#' set.seed(1234)
#' segTabs$minor_cn <- sample(c(0, 1), size = nrow(segTabs), replace = TRUE)
#' cn <- read_copynumber(segTabs,
#'   seg_cols = c("chromosome", "start", "end", "segVal"),
#'   genome_measure = "wg", complement = TRUE, add_loh = TRUE
#' )
#'
#' df <- get_pLOH_score(cn)
#' df
#'
#' df2 <- get_pLOH_score(cn@data)
#' df2
#'
#' @testexamples
#' expect_equal(nrow(df), 10L)
#' expect_equal(nrow(df2), 10L)
get_pLOH_score <- function(data, rm_chrs = c("chrX", "chrY"), genome_build = "hg19") {
  stopifnot(is.data.frame(data) | inherits(data, "CopyNumber") | "chrY" %in% rm_chrs)
  if (!is.data.frame(data)) {
    data <- data@data
  }

  nc_cols <- c("chromosome", "start", "end", "segVal", "minor_cn", "sample")
  if (!all(nc_cols %in% colnames(data))) {
    stop("Invalid input, it must contain columns: ", paste(nc_cols, collapse = " "))
  }

  data <- data.table::as.data.table(data)[!chromosome %in% rm_chrs, nc_cols, with = FALSE]

  # Make sure the minor_cn is really minor copy
  data$minor_cn <- pmin(data$segVal - data$minor_cn, data$minor_cn)

  chr_size <- get_genome_annotation(genome_build = genome_build)
  chr_size <- chr_size[!chr_size$chrom %in% rm_chrs, ]
  chr_size <- sum(chr_size$size)

  data[segVal >= 1 & minor_cn == 0, list(pLOH = sum(end - start + 1L) / chr_size), by = "sample"]
}


utils::globalVariables("minor_cn")
