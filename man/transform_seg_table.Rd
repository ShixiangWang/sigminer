% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transform_seg_table.R
\name{transform_seg_table}
\alias{transform_seg_table}
\title{Transform Copy Number Table}
\usage{
transform_seg_table(
  data,
  genome_build = c("hg19", "hg38", "T2T", "mm10", "mm9", "ce11"),
  ref_type = c("cytoband", "gene"),
  values_fill = NA,
  values_fn = function(x, ...) {
     round(mean(x, ...))
 },
  resolution_factor = 1L
)
}
\arguments{
\item{data}{a \code{CopyNumber} object or a data.frame containing
at least 'chromosome', 'start', 'end', 'segVal', 'sample' these columns.}

\item{genome_build}{genome build version, used when \code{data} is a \code{data.frame}, should be 'hg19' or 'hg38'.}

\item{ref_type}{annotation data type used for constructing matrix.}

\item{values_fill}{Optionally, a (scalar) value that specifies what each
\code{value} should be filled in with when missing.

This can be a named list if you want to apply different fill values to
different value columns.}

\item{values_fn}{Optionally, a function applied to the value in each cell
in the output. You will typically use this when the combination of
\code{id_cols} and \code{names_from} columns does not uniquely identify an
observation.

This can be a named list if you want to apply different aggregations
to different \code{values_from} columns.}

\item{resolution_factor}{an integer to control the resolution.
When it is \code{1} (default), compute frequency in each cytoband.
When it is \code{2}, use compute frequency in each half cytoband.}
}
\value{
a \code{data.table}.
}
\description{
Transform Copy Number Table
}
\examples{
\donttest{
load(system.file("extdata", "toy_copynumber.RData",
  package = "sigminer", mustWork = TRUE
))
# Compute the mean segVal in each cytoband
x <- transform_seg_table(cn, resolution_factor = 1)
x
# Compute the mean segVal in each half-cytoband
x2 <- transform_seg_table(cn, resolution_factor = 2)
x2
}
}
