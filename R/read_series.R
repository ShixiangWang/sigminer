#=====================================================
# Function: read data as S4 object
#=====================================================


# Read MAF data -----------------------------------------------------------

#' @inherit maftools::read.maf
#' @family read genomic variation data function series
#' @examples
#' \dontrun{
#' laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#' laml <- read_maf(maf = laml.maf)
#' }
#'

read_maf = function(
  maf, clinicalData = NULL, removeDuplicatedVariants = TRUE,
  useAll = TRUE, gisticAllLesionsFile = NULL,
  gisticAmpGenesFile = NULL, gisticDelGenesFile = NULL,
  gisticScoresFile = NULL, cnLevel = "all", cnTable = NULL,
  isTCGA = FALSE, vc_nonSyn = NULL, verbose = TRUE
) {
  maftools::read.maf(
    maf, clinicalData = clinicalData,
    removeDuplicatedVariants = removeDuplicatedVariants,
    useAll = useAll,
    gisticAllLesionsFile = gisticAllLesionsFile,
    gisticAmpGenesFile = gisticAmpGenesFile,
    gisticDelGenesFile = gisticDelGenesFile,
    gisticScoresFile = gisticScoresFile,
    cnLevel = cnLevel, cnTable = cnTable,
    isTCGA = isTCGA, vc_nonSyn = vc_nonSyn,
    verbose = verbose
  )
}


# Read absolute copy number profile ------------------------------------------------

#' @title  Read absolute copy number profile
#' @description Read **absolute** copy number profile for preparing CNV signature
#' analysis.
#' @param input a `data.frame` or a file or a directory contains copy number profile.
#' @param pattern an optional regular expression used to select part of files if
#' `input` is a directory, more detail please see [list.files] function.
#' @param ignore_case logical. Should pattern-matching be case-insensitive?
#' @param seg_cols four characters used to specify chromosome, start position,
#'  end position and copy number value in `input`, respectively.
#'  Default use names from ABSOLUTE calling result.
#' @param samp_col a character used to specify the sample column name. If `input`
#' is a directory and cannot find `samp_col`, sample names will use file names
#' (set this parameter to `NULL` is recommended in this case).
#' @param use_all default is `FALSE`. If `True`, use all columns from raw input.
#' @param min_segnum minimal number of copy number segments within a sample.
#' @param genome_build genome build version, should be 'hg19' or 'hg38'.
#' @param clinical_data a `data.frame` representing clinical data
#' associated with each sample in copy number profile.
#' @param verbose print extra messages.
#' @param ... other parameters pass to [data.table::fread()]
#' @author Shixiang Wang <w_shixiang@163.com>
#' @return a [CopyNumber] object
#' @export
#' @family read genomic variation data function series
read_copynumber = function(input,
                           pattern = NULL,
                           ignore_case = FALSE,
                           seg_cols = c("Chromosome", "Start.bp", "End.bp", "modal_cn"),
                           samp_col = "sample",
                           use_all = FALSE,
                           min_segnum = 22,
                           genome_build = c("hg19", "hg38"),
                           clinical_data = NULL,
                           verbose = FALSE,
                           ...) {
  stopifnot(
    is.character(samp_col),
    length(samp_col) == 1,
    min_segnum >= 0,
    is.logical(verbose)
  )

  if (dir.exists(input)) {
    if (verbose) message("Treat input as a directory...")
    if (length(input) != 1) {
      stop("Only can take one directory as input!")
    }
    # get files and exclude directories
    all.files <- list.files(
      path = input,
      pattern = pattern,
      all.files = FALSE,
      recursive = FALSE,
      ignore.case = ignore_case
    )
    files = all.files[!file.info(file.path(input, all.files))$isdir]
    if (length(files) == 0)
      stop("No files exist, please check!")
    files_path = file.path(input, files)
    data_list = list()
    dropoff_list = list()

    # read files
    for (i in seq_along(files_path)) {
      if (verbose) message("Reading ", files_path[i])
      temp = data.table::fread(files_path[i], ...)
      if (!all(seg_cols %in% colnames(temp)))
        stop("Not all seg_cols are in file, please check.")

      if (length(samp_col %in% colnames(temp))==0 | !(samp_col %in% colnames(temp))) {
        if (verbose) message("Select file names as sample names.")
        temp[, "sample"] = files[i]
        sample_col = "sample"

      }

      tempName = unique(temp[[samp_col]])
      if (length(tempName) > 1) {
        stop("When input is a directory, a file can only contain one sample.")
      }

      # set column order
      data.table::setcolorder(temp, neworder = c(seg_cols, samp_col))
      new_cols = c("chromosome", "start", "end", "segVal", "sample")
      colnames(temp)[1:5] = new_cols
      if (!use_all) temp = temp[, new_cols, with = FALSE]
      if (nrow(temp) < min_segnum) {
        dropoff_list[[tempName]] =  temp
      } else {
        data_list[[tempName]] = temp
      }

    }

    if (length(data_list) >= 1) {
      data_df = data.table::rbindlist(data_list, use.names = TRUE, fill = TRUE)
    } else {
      data_df = data.table::data.table()
    }

    if (length(dropoff_list) >= 1) {
      dropoff_df = data.table::rbindlist(dropoff_list, use.names = TRUE, fill = TRUE)
    } else {
      dropoff_df = data.table::data.table()
    }

  } else if (all(is.character(input)) | is.data.frame(input)) {

    if (!is.data.frame(input)) {
      if (verbose) message("Treat input as a file...")
      if (length(input) > 1) {
        stop("Muliple files are not a valid input, please use directory as input.")
      }
      if (!file.exists(input))
        stop("Input file not exists.")

      temp = data.table::fread(input, ...)
    } else {
      if (verbose) message("Treat input as a data frame...")
      temp = data.table::as.data.table(input)
    }

    #--- Check column names
    if (is.null(samp_col))
      stop("'samp_col' parameter must set!")
    if (!all(seg_cols %in% colnames(temp)))
      stop("Not all seg_cols are in file, please check.")
    if (!(samp_col %in% colnames(temp)))
      stop("Column ", samp_col, " does not exist.")

    #--- Set column order
    data.table::setcolorder(temp, neworder = c(seg_cols, samp_col))
    new_cols = c("chromosome", "start", "end", "segVal", "sample")
    colnames(temp)[1:5] = new_cols
    if (!use_all) temp = temp[, new_cols, with = FALSE]

    dropoff_samples = temp[, .N, by = .(sample)][N < min_segnum][["sample"]]
    keep_samples = base::setdiff(unique(temp[["sample"]]), dropoff_samples)

    data_df = temp[sample %in% keep_samples]
    dropoff_df = temp[sample %in% dropoff_samples]

    } else {
      stop("Invalid input.")
    }

  #--- match genome build
  genome_build = match.arg(genome_build)

  annot = get_LengthFraction(data_df,
                             genome_build = genome_build,
                             seg_cols = new_cols[1:4],
                             samp_col = new_cols[5])

  if (is.null(clinical_data)) {
    CopyNumber(
      data =  data_df,
      genome_build = genome_build,
      annotation = annot,
      dropoff.segs = dropoff_df,
      clinical.data = data.table::data.table()
    )
  } else {
    CopyNumber(
      data =  data_df,
      genome_build = genome_build,
      annotation = annot,
      dropoff.segs = dropoff_df,
      clinical.data = data.table::as.data.table(clinical_data)
    )
  }
}



# Read genomic variation --------------------------------------------------

#' Read genomic variation profile
#'
#' @description Read [CopyNumber] and [MAF] object as a new S4 object [GenomicVariation]
#' for uniform variation analysis. **The function is initialized to construct structure of sigminer, please dont use it for now**.
#' @param copynumber a [CopyNumber] object
#' @param maf a [MAF] object
#' @param clinical_data clinical.data data associated with each sample in copy number profile
#' and MAF.
#' @author Shixiang Wang <w_shixiang@163.com>
#' @return a [GenomicVariation] object
#' @export
#' @family read genomic variation data function series
read_variation = function(copynumber, maf, clinical_data = NULL) {

  if (is.null(clinical_data)) {
    clinical_data = data.table::data.table()
  } else if (!inherits(clinical_data, "data.table")) {
    clinical_data = data.table::as.data.table(clinical_data)
  }

  GenomicVariation(
    CopyNumber = copynumber,
    MAF = maf,
    clinical.data = clinical_data
  )
}
