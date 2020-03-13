#' Get Genome Annotation
#'
#' @param data_type 'chr_size' for chromosome size,
#' 'centro_loc' for location of centromeres
#' @param genome_build one of 'hg19', 'hg38'
#' @param chrs chromosomes start with 'chr'
#' @return a `data.frame` containing annotation data
#' @examples
#' df1 <- get_genome_annotation()
#' df1
#'
#' df2 <- get_genome_annotation(genome_build = "hg38")
#' df2
#' @testexamples
#' expect_equal(nrow(df1), nrow(df2))
#' @export
get_genome_annotation <- function(data_type = c("chr_size", "centro_loc"),
                                  chrs = paste0("chr", c(1:22, "X", "Y")),
                                  genome_build = c("hg19", "hg38")) {
  data_type <- match.arg(data_type)
  genome_build <- match.arg(genome_build)

  if (data_type == "chr_size") {
    if (genome_build == "hg19") {
      data("chromsize.hg19",
        package = "sigminer",
        envir = environment()
      )
      res <- chromsize.hg19 %>%
        dplyr::filter(.data$chrom %in% chrs) %>%
        dplyr::arrange(factor(.data$chrom, chrs))
    } else {
      data("chromsize.hg38",
        package = "sigminer",
        envir = environment()
      )
      res <- chromsize.hg38 %>%
        dplyr::filter(.data$chrom %in% chrs) %>%
        dplyr::arrange(factor(.data$chrom, chrs))
    }
  } else {
    if (genome_build == "hg19") {
      data("centromeres.hg19",
        package = "sigminer",
        envir = environment()
      )
      res <- centromeres.hg19 %>%
        dplyr::filter(.data$chrom %in% chrs) %>%
        dplyr::arrange(factor(.data$chrom, chrs))
    } else {
      data("centromeres.hg38",
        package = "sigminer",
        envir = environment()
      )
      res <- centromeres.hg38 %>%
        dplyr::filter(.data$chrom %in% chrs) %>%
        dplyr::arrange(factor(.data$chrom, chrs))
    }
  }

  res
}


utils::globalVariables(
  c(
    "centromeres.hg19",
    "centromeres.hg38",
    "chromsize.hg19",
    "chromsize.hg38"
  )
)
