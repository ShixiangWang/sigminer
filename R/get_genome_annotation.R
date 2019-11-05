#' Get Genome Annotation
#'
#' @param data_type 'chr_size' for chromosome size,
#' 'centro_loc' for location of centromeres
#' @param genome_build one of 'hg19', 'hg38'
#' @return a `data.frame` containing annotation data
get_genome_annotation = function(data_type = c("chr_size", "centro_loc"),
                                 genome_build = c("hg19", "hg38")) {

  data_type = match.arg(data_type)
  genome_build = match.arg(genome_build)

  if (data_type == "chr_size") {
    if (genome_build == "hg19") {
      data("chromsize.hg19",
           package = "sigminer",
           envir = environment()
      )
      res <- chromsize.hg19
    } else {
      data("chromsize.hg38",
           package = "sigminer",
           envir = environment()
      )
      res <- chromsize.hg38
    }
  } else {
    if (genome_build == "hg19") {
      data("centromeres.hg19",
           package = "sigminer",
           envir = environment()
      )
      res <- chromsize.hg19
    } else {
      data("centromeres.hg38",
           package = "sigminer",
           envir = environment()
      )
      res <- chromsize.hg38
    }
  }

  res
}
