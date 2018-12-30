
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
#' @export MAF
#' @exportClass MAF

## MAF object
MAF = setClass(
  Class = 'MAF',
  slots =  c(
    data = 'data.table',
    variants.per.sample = 'data.table',
    variant.type.summary = 'data.table',
    variant.classification.summary = 'data.table',
    gene.summary = 'data.table',
    summary = 'data.table',
    maf.silent = 'data.table',
    clinical.data = 'data.table'
  )
)

setMethod(
  f = 'show',
  signature = 'MAF',
  definition = function(object) {
    cat(paste('An object of class ', class(object), "\n"))
    print(object@summary)
  }
)


#' Class CopyNumber
#' @description S4 class for storing summarized absolute copy number profile.
#' @slot data data.table of absolute copy number calling.
#' @slot genome_build genome build version, should be one of 'hg19' or 'hg38'.
#' @slot annotation data.table of annotation for copy number segments.
#' @slot dropoff.segs data.table of copy number segments dropped from raw input.
#' @slot clinical.data data associated with each sample in copy number profile.
#' @export CopyNumber
#' @exportClass CopyNumber
CopyNumber = setClass(
  Class = "CopyNumber",
  slots = c(
    data = "data.table",
    genome_build = "character",
    annotation = "data.table",
    dropoff.segs = "data.table",
    clinical.data = "data.table"
  )
)


setMethod(
  f = 'show',
  signature = 'CopyNumber',
  definition = function(object) {
    cat(paste('An object of class ', class(object), "\n"))
    print(object@data)
  }
)


#' Class GenomicVariation
#' @description S4 class for storing summaarized genomic variation profile.
#' @slot CopyNumber a `CopyNumber` object.
#' @slot MAF a `MAF` object.
#' @slot clinical.data data associated with each sample in copy number profile
#' and MAF. This slot is intersection of `clinical.data` slot in both `CopyNumber`
#' object and `MAF` object.
#' @export GenomicVariation
#' @exportClass GenomicVariation
GenomicVariation = setClass(
  Class = "GenomicVariation",
  slots = c(CopyNumber = "CopyNumber",  MAF = "MAF",
  clinical.data = "data.table")
)

setMethod(
  f = 'show',
  signature = 'GenomicVariation',
  definition = function(object) {
    cat(paste('An object of class ', class(object), "\n"))
    print(str(object, max.level = 2))
  }
)
