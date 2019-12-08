# Read absolute copy number profile ------------------------------------------------

#' @title  Read Absolute Copy Number Profile
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
#' @param max_copynumber bigger copy number within a sample will be reset to this value.
#' @param genome_build genome build version, should be 'hg19' or 'hg38'.
#' @param genome_measure default is 'called', can be 'wg' or 'called'.
#' Set 'called' will use called segments size to compute total size for CNA burden calculation,
#' this option is useful for WES and target sequencing.
#' Set 'wg' will use autosome size from genome build, this option is useful for WGS, SNP etc..
#' @param complement if `TRUE`, complement chromosome (except 'Y') does not show in input data
#' with normal copy 2 and force `use_all` to `FALSE` (no matter what user input).
#' @param verbose print extra messages.
#' @param ... other parameters pass to [data.table::fread()]
#' @author Shixiang Wang <w_shixiang@163.com>
#' @return a [CopyNumber] object
#' @export
#' @examples
#' # Load toy dataset of absolute copynumber profile
#' load(system.file("extdata", "toy_segTab.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' cn <- read_copynumber(segTabs,
#'   seg_cols = c("chromosome", "start", "end", "segVal"),
#'   genome_build = "hg19", complement = FALSE, verbose = TRUE
#' )
#' @seealso [read_maf] for reading mutation data to [MAF] object.
read_copynumber <- function(input,
                            pattern = NULL,
                            ignore_case = FALSE,
                            seg_cols = c("Chromosome", "Start.bp", "End.bp", "modal_cn"),
                            samp_col = "sample",
                            use_all = FALSE,
                            min_segnum = 0L,
                            max_copynumber = 20L,
                            genome_build = c("hg19", "hg38"),
                            genome_measure = c("called", "wg"),
                            complement = TRUE,
                            verbose = FALSE,
                            ...) {
  stopifnot(
    is.character(samp_col),
    length(samp_col) == 1,
    min_segnum >= 0,
    is.logical(verbose)
  )

  #--- match genome build
  genome_build <- match.arg(genome_build)
  genome_measure <- match.arg(genome_measure)

  # get chromosome lengths
  valid_chr <- c(paste0("chr", 1:22), "chrX", "chrY")
  chrlen <- get_genome_annotation(
    data_type = "chr_size",
    chrs = valid_chr,
    genome_build = genome_build
  )
  data.table::setDT(chrlen)

  if (tryCatch(dir.exists(input), error = function(e) FALSE)) {
    if (verbose) message("Treating input as a directory...")
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
    files <- all.files[!file.info(file.path(input, all.files))$isdir]
    if (length(files) == 0) {
      stop("No files exist, please check!")
    }
    files_path <- file.path(input, files)
    data_list <- list()
    dropoff_list <- list()

    # read files
    for (i in seq_along(files_path)) {
      if (verbose) message("Reading ", files_path[i])
      temp <- data.table::fread(files_path[i], ...)
      if (!all(seg_cols %in% colnames(temp))) {
        stop("Not all seg_cols are in file, please check.")
      }

      if (length(samp_col %in% colnames(temp)) == 0 | !(samp_col %in% colnames(temp))) {
        if (verbose) message("Select file names as sample names.")
        temp[, "sample"] <- files[i]
        sample_col <- "sample"
      }

      tempName <- unique(temp[[samp_col]])
      if (length(tempName) > 1) {
        stop("When input is a directory, a file can only contain one sample.")
      }

      # set column order
      data.table::setcolorder(temp, neworder = c(seg_cols, samp_col))
      new_cols <- c("chromosome", "start", "end", "segVal", "sample")
      colnames(temp)[1:5] <- new_cols

      if (any(is.na(temp$segVal))) {
        temp <- temp[!is.na(temp$segVal)]
      }

      # unify chromosome column
      if (verbose) message("Checking chromosome names...")
      temp[, chromosome := sub(
        pattern = "chr",
        replacement = "chr",
        x = as.character(chromosome),
        ignore.case = TRUE
      )]
      if (any(!grepl("chr", temp$chromosome))) {
        temp$chromosome[!grepl("chr", temp$chromosome)] <-
          paste0("chr", temp$chromosome[!grepl("chr", temp$chromosome)])
      }
      temp[, chromosome := sub(
        pattern = "x",
        replacement = "X",
        x = chromosome
      )]
      temp[, chromosome := sub(
        pattern = "y",
        replacement = "Y",
        x = chromosome
      )]

      # detect and transform chromosome 23 to "X"
      temp[["chromosome"]] <- sub("23", "X", temp[["chromosome"]])
      # detect and transform chromosome 24 to "Y"
      temp[["chromosome"]] <- sub("24", "Y", temp[["chromosome"]])

      if (complement) {
        # complement value 2 (normal copy) to chromosome not called
        if (verbose) message("Fill value 2 (normal copy) to uncalled chromosomes.")
        miss_index <- !valid_chr %in% unique(temp[["chromosome"]])
        miss_index[length(miss_index)] <- FALSE # disable Y
        if (any(miss_index)) {
          comp_df <- temp[rep(1, sum(miss_index))]
          comp_df[, c("chromosome", "start", "end", "segVal") := .(
            chrlen[["chrom"]][miss_index],
            1,
            chrlen[["size"]][miss_index],
            2
          )]
          temp <- rbind(temp, comp_df)
        }
        if (verbose) message("'complement' option is TRUE, thus use_all automatically set to FALSE.")
        use_all <- FALSE
      }

      if (!use_all) temp <- temp[, new_cols, with = FALSE]
      if (nrow(temp) < min_segnum) {
        dropoff_list[[tempName]] <- temp
      } else {
        data_list[[tempName]] <- temp
      }
    }

    if (length(data_list) >= 1) {
      data_df <- data.table::rbindlist(data_list, use.names = TRUE, fill = TRUE)
    } else {
      data_df <- data.table::data.table()
    }

    if (length(dropoff_list) >= 1) {
      dropoff_df <- data.table::rbindlist(dropoff_list, use.names = TRUE, fill = TRUE)
    } else {
      dropoff_df <- data.table::data.table()
    }
  } else if (all(is.character(input)) | is.data.frame(input)) {
    if (!is.data.frame(input)) {
      if (verbose) message("Treating input as a file...")
      if (length(input) > 1) {
        stop("Muliple files are not a valid input, please use directory as input.")
      }
      if (!file.exists(input)) {
        stop("Input file not exists.")
      }

      temp <- data.table::fread(input, ...)
    } else {
      if (verbose) message("Treating input as a data frame...")
      temp <- data.table::as.data.table(input)
    }

    #--- Check column names
    if (is.null(samp_col)) {
      stop("'samp_col' parameter must set!")
    }
    if (!all(seg_cols %in% colnames(temp))) {
      stop("Not all seg_cols are in file, please check.")
    }
    if (!(samp_col %in% colnames(temp))) {
      stop("Column ", samp_col, " does not exist.")
    }

    #--- Set column order
    data.table::setcolorder(temp, neworder = c(seg_cols, samp_col))
    new_cols <- c("chromosome", "start", "end", "segVal", "sample")
    colnames(temp)[1:5] <- new_cols

    if (any(is.na(temp$segVal))) {
      temp <- temp[!is.na(temp$segVal)]
    }

    # unify chromosome column
    if (verbose) message("Checking chromosome names...")
    temp[, chromosome := sub(
      pattern = "chr",
      replacement = "chr",
      x = as.character(chromosome),
      ignore.case = TRUE
    )]
    if (any(!grepl("chr", temp$chromosome))) {
      temp$chromosome[!grepl("chr", temp$chromosome)] <-
        paste0("chr", temp$chromosome[!grepl("chr", temp$chromosome)])
    }
    temp[, chromosome := sub(
      pattern = "x",
      replacement = "X",
      x = chromosome
    )]
    temp[, chromosome := sub(
      pattern = "y",
      replacement = "Y",
      x = chromosome
    )]

    # detect and transform chromosome 23 to "X"
    temp[["chromosome"]] <- sub("23", "X", temp[["chromosome"]])
    # detect and transform chromosome 24 to "Y"
    temp[["chromosome"]] <- sub("24", "Y", temp[["chromosome"]])

    if (complement) {
      # complement value 2 (normal copy) to chromosome not called
      if (verbose) message("Fill value 2 (normal copy) to uncalled chromosomes.")

      comp <- data.table::data.table()
      for (i in unique(temp[["sample"]])) {
        tmp_sample <- temp[i, on = "sample"]
        miss_index <- !valid_chr %in% unique(tmp_sample[["chromosome"]])
        miss_index[length(miss_index)] <- FALSE # disable Y
        if (any(miss_index)) {
          comp_df <- tmp_sample[rep(1, sum(miss_index))]
          comp_df[, c("chromosome", "start", "end", "segVal") := .(
            chrlen[["chrom"]][miss_index],
            1,
            chrlen[["size"]][miss_index],
            2
          )]
          comp <- rbind(comp, comp_df)
        }
      }
      temp <- rbind(temp, comp)
      if (verbose) message("complement is TRUE, thus use_all automatically set to FALSE.")
      use_all <- FALSE
    }

    if (!use_all) temp <- temp[, new_cols, with = FALSE]

    dropoff_samples <- temp[, .N, by = .(sample)][N < min_segnum][["sample"]]
    keep_samples <- base::setdiff(unique(temp[["sample"]]), dropoff_samples)

    data_df <- temp[sample %in% keep_samples]
    dropoff_df <- temp[sample %in% dropoff_samples]
  } else {
    stop("Invalid input.")
  }

  if (!all(data_df$chromosome %in% valid_chr)) {
    if (verbose) message("Filter some invalid segments... (not as 1:22 and X, Y)")
    data_drop <- data_df[!chromosome %in% valid_chr]
    if (nrow(dropoff_df) >= 1) {
      dropoff_df <- base::rbind(dropoff_df, data_drop)
    } else {
      dropoff_df <- data_drop
    }

    data_df <- data_df[chromosome %in% valid_chr]
  }

  if (verbose) {
    message("Segments info:")
    message("    Keep - ", nrow(data_df))
    message("  Filter - ", nrow(dropoff_df))
  }

  # reset copy number for high copy number segments
  data_df$segVal[data_df$segVal > max_copynumber] <- max_copynumber
  # make sure seg value is integer
  data_df[["segVal"]] <- as.integer(round(data_df[["segVal"]]))
  # make sure position is numeric
  data_df$start = as.numeric(data_df$start)
  data_df$end = as.numeric(data_df$end)

  if (verbose) message("Anotating...")
  annot <- get_LengthFraction(data_df,
    genome_build = genome_build,
    seg_cols = new_cols[1:4],
    samp_col = new_cols[5]
  )
  if (verbose) message("\nSummarizing per sample...")
  sum_sample <- get_cnsummary_sample(data_df,
    genome_build = genome_build,
    genome_measure = genome_measure
  )

  if (verbose) message("Done!")

  res <- CopyNumber(
    data = data_df,
    summary.per.sample = sum_sample,
    genome_build = genome_build,
    genome_measure = genome_measure,
    annotation = annot,
    dropoff.segs = dropoff_df
  )

  res <- validate_segTab(res, verbose = verbose)
  res
}


# Global variables --------------------------------------------------------

utils::globalVariables(
  c(
    ".",
    "N",
    ".N"
  )
)
