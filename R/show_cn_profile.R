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
#' @param genome_build genome build version, used when `data` is a `data.frame`, , should be 'hg19' or 'hg38'.
#' @param nrow number of rows in the plot grid when multiple samples are selected.
#' @param ncol number of columns in the plot grid when multiple samples are selected.
#' @param return_plotlist default is `FALSE`, if `TRUE`, return a plot list instead of a combined plot.
#'
#' @return a `ggplot` object or a `list`
#' @export
#'
#' @examples
#' print("hello")
show_cn_profile = function(data, samples=NULL, genome_build = c("hg19", "hg38"),
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

# Get chromosome size

  data.table::setDT(data)
  if (!is.null(samples)) {
    data = data[samples %in% data$sample]
  }

}

# Generate ggplot2 using "chromosome", "start", "end", "segVal"
show_one_cn_profile = function(df) {

}
