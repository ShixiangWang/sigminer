% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show_cn_group_profile.R
\name{show_cn_group_profile}
\alias{show_cn_group_profile}
\title{Show Summary Copy Number Profile for Sample Groups}
\usage{
show_cn_group_profile(
  data,
  groups = NULL,
  fill_area = TRUE,
  cols = NULL,
  chrs = paste0("chr", c(1:22, "X")),
  genome_build = c("hg19", "hg38", "T2T", "mm10", "mm9", "ce11"),
  cutoff = 2L,
  resolution_factor = 1L,
  force_y_limit = TRUE,
  highlight_genes = NULL,
  repel = FALSE,
  nrow = NULL,
  ncol = NULL,
  return_plotlist = FALSE
)
}
\arguments{
\item{data}{a \code{CopyNumber} object or a data.frame containing
at least 'chromosome', 'start', 'end', 'segVal', 'sample' these columns.}

\item{groups}{a named list or a column name for specifying groups.}

\item{fill_area}{default is \code{TRUE}, fill area with colors.}

\item{cols}{length-2 colors for AMP and DEL.}

\item{chrs}{chromosomes start with 'chr'.}

\item{genome_build}{genome build version, used when \code{data} is a \code{data.frame}, should be 'hg19' or 'hg38'.}

\item{cutoff}{copy number value cutoff for splitting data into AMP and DEL.
The values equal to cutoff are discarded. Default is \code{2}, you can also set
a length-2 vector, e.g. \code{c(2, 2)}.}

\item{resolution_factor}{an integer to control the resolution.
When it is \code{1} (default), compute frequency in each cytoband.
When it is \code{2}, use compute frequency in each half cytoband.}

\item{force_y_limit}{default is \code{TRUE}, force multiple plots}

\item{highlight_genes}{gene list to highlight.
have same y ranges. You can also set a length-2 numeric value.}

\item{repel}{if \code{TRUE} (default is \code{FALSE}), repel highlight genes to
avoid overlap.}

\item{nrow}{number of rows in the plot grid when multiple samples are selected.}

\item{ncol}{number of columns in the plot grid when multiple samples are selected.}

\item{return_plotlist}{default is \code{FALSE}, if \code{TRUE}, return a plot list instead of a combined plot.}
}
\value{
a (list of) \code{ggplot} object.
}
\description{
Show Summary Copy Number Profile for Sample Groups
}
\examples{
load(system.file("extdata", "toy_copynumber.RData",
  package = "sigminer", mustWork = TRUE
))

p1 <- show_cn_group_profile(cn)
p1
\donttest{
ss <- unique(cn@data$sample)
p2 <- show_cn_group_profile(cn, groups = list(a = ss[1:5], b = ss[6:10]))
p2
p3 <- show_cn_group_profile(cn,
  groups = list(g1 = ss[1:5], g2 = ss[6:10]),
  force_y_limit = c(-1, 1), nrow = 2
)
p3

## Set custom cutoff for custom data
data <- cn@data
data$segVal <- data$segVal - 2L
p4 <- show_cn_group_profile(data,
  groups = list(g1 = ss[1:5], g2 = ss[6:10]),
  force_y_limit = c(-1, 1), nrow = 2,
  cutoff = c(0, 0)
)
p4

## Add highlight gene
p5 <- show_cn_group_profile(cn, highlight_genes = c("TP53", "EGFR"))
p5
}
}
