
# Set Class ---------------------------------------------------------------

#' Class MAF
#' @description S4 class for storing summarized MAF.
#' @slot data data.table of MAF file containing all non-synonymous variants.
#' @slot variants.per.sample table containing variants per sample
#' @slot variant.type.summary table containing variant types per sample
#' @slot variant.classification.summary table containing variant classification per sample
#' @slot gene.summary table containing variant classification per gene
#' @slot summary table with basic MAF summary stats
#' @slot maf.silent subset of main MAF containing only silent variants
#' @slot clinical.data clinical data associated with each sample/Tumor_Sample_Barcode in MAF.
#' @import methods
#' @import data.table
#' @export MAF
#' @exportClass MAF

## MAF object
MAF <- setClass(
  Class = "MAF",
  slots = c(
    data = "data.table",
    variants.per.sample = "data.table",
    variant.type.summary = "data.table",
    variant.classification.summary = "data.table",
    gene.summary = "data.table",
    summary = "data.table",
    maf.silent = "data.table",
    clinical.data = "data.table"
  )
)

setMethod(
  f = "show",
  signature = "MAF",
  definition = function(object) {
    cat(paste("An object of class", class(object), "\n"))
    cat("======================\n")
    print(object@summary)
  }
)


#' Class CopyNumber
#' @description S4 class for storing summarized absolute copy number profile.
#' @slot data data.table of absolute copy number calling.
#' @slot summary.per.sample data.table of copy number variation summary per sample.
#' @slot genome_build genome build version, should be one of 'hg19' or 'hg38'.
#' @slot genome_measure Set 'called' will use autosomo called segments size to compute total size
#' for CNA burden calculation, this option is useful for WES and target sequencing.
#' Set 'wg' will autosome size from genome build, this option is useful for WGS, SNP etc..
#' @slot annotation data.table of annotation for copy number segments.
#' @slot dropoff.segs data.table of copy number segments dropped from raw input.
#' @slot clinical.data data associated with each sample in copy number profile.
#' @export CopyNumber
#' @exportClass CopyNumber
CopyNumber <- setClass(
  Class = "CopyNumber",
  slots = c(
    data = "data.table",
    summary.per.sample = "data.table",
    genome_build = "character",
    genome_measure = "character",
    annotation = "data.table",
    dropoff.segs = "data.table",
    clinical.data = "data.table"
  )
)


setMethod(
  f = "show",
  signature = "CopyNumber",
  definition = function(object) {
    cat(paste("An object of class", class(object), "\n"))
    cat("=============================\n")
    print(object@summary.per.sample)
  }
)


#' Class GenomicVariation
#' @description S4 class for storing summarized genomic variation profile.
#' @slot CopyNumber a `CopyNumber` object.
#' @slot MAF a `MAF` object.
#' @slot clinical.data data associated with each sample in copy number profile
#' and MAF. This slot is intersection of `clinical.data` slot in both `CopyNumber`
#' object and `MAF` object.
#' @export GenomicVariation
#' @exportClass GenomicVariation
GenomicVariation <- setClass(
  Class = "GenomicVariation",
  slots = c(
    CopyNumber = "CopyNumber", MAF = "MAF",
    clinical.data = "data.table"
  )
)

setMethod(
  f = "show",
  signature = "GenomicVariation",
  definition = function(object) {
    cat(paste("An object of class", class(object), "\n"))
    cat("===================================\n")
    str(object, max.level = 2)
  }
)



# Methods -----------------------------------------------------------------

#' Subsetting CopyNumber object
#'
#' Subset data slot of [CopyNumber] object, un-selected rows will move to
#' dropoff.segs slot, annotation slot will update in the same way.
#'
#' @param x a [CopyNumber] object to be subsetted.
#' @param subset logical expression indicating rows to keep.
#' @param ... further arguments to be passed to or from other methods.
#' Useless here.
#' @author Shixiang Wang
#' @return a [CopyNumber] object
#' @export
subset.CopyNumber <- function(x, subset = TRUE, ...) {
  data <- x@data
  enclos <- parent.frame()
  subset <- substitute(subset)
  row_selector <- eval(subset, data, enclos)

  x@dropoff.segs <- rbind(x@dropoff.segs, data[!row_selector])
  x@data <- data[row_selector]
  x@annotation <- x@annotation[row_selector]

  x@summary.per.sample <- get_cnsummary_sample(x@data,
    genome_build = x@genome_build,
    genome_measure = x@genome_measure
  )
  x
}
