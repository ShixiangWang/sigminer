#' Get Genome Annotation
#'
#' @param data_type 'chr_size' for chromosome size,
#' 'centro_loc' for location of centromeres,
#' 'cytobands' for location of chromosome cytobands
#' and 'transcript' for location of transcripts.
#' @param genome_build one of 'hg19', 'hg38'
#' @param chrs chromosomes start with 'chr'
#' @return a `data.frame` containing annotation data
#' @examples
#' df1 <- get_genome_annotation()
#' df1
#'
#' df2 <- get_genome_annotation(genome_build = "hg38")
#' df2
#'
#' df3 <- get_genome_annotation(data_type = "centro_loc")
#' df3
#'
#' df4 <- get_genome_annotation(data_type = "centro_loc", genome_build = "hg38")
#' df4
#'
#' df5 <- get_genome_annotation(data_type = "cytobands")
#' df5
#'
#' df6 <- get_genome_annotation(data_type = "cytobands", genome_build = "hg38")
#' df6
#' @testexamples
#' expect_equal(nrow(df1), nrow(df2))
#' expect_equal(nrow(df3), nrow(df4))
#' expect_equal(nrow(df5), nrow(df6))
#'
#' expect_equal(identical(df1, df2), FALSE)
#' expect_equal(identical(df3, df4), FALSE)
#' expect_equal(identical(df5, df6), FALSE)
#' @export
get_genome_annotation <- function(data_type = c("chr_size", "centro_loc", "cytobands", "transcript"),
                                  chrs = paste0("chr", c(1:22, "X", "Y")),
                                  genome_build = c("hg19", "hg38", "mm10", "mm9")) {
  data_type <- match.arg(data_type)
  genome_build <- match.arg(genome_build)

  prefix <- switch(data_type,
    chr_size = "chromsize",
    centro_loc = "centromeres",
    cytobands = "cytobands",
    transcript = "transcript"
  )

  annot <- get(paste(prefix, genome_build, sep = "."), envir = as.environment("package:sigminer"))
  res <- annot %>%
    dplyr::filter(.data$chrom %in% chrs) %>%
    dplyr::arrange(factor(.data$chrom, chrs))

  res
}
