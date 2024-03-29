% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_copynumber_ascat.R
\name{read_copynumber_ascat}
\alias{read_copynumber_ascat}
\title{Read Copy Number Data from ASCAT Result Files}
\usage{
read_copynumber_ascat(x)
}
\arguments{
\item{x}{one or more \code{.rds} format files which contains \code{ASCAT} object from result of \code{ascat.runAscat()}
in \href{https://github.com/VanLoo-lab/ascat}{\strong{ASCAT}} package.}
}
\value{
a tidy \code{list}.
}
\description{
Note, the result is not a \code{CopyNumber} object, you need to generate it
by yourself.
}
