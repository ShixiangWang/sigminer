# Read absolute copy number profile ------------------------------------------------

#' @title  Read Absolute Copy Number Profile
#' @description Read **absolute** copy number profile for preparing CNV signature
#' analysis. See detail part of [sig_tally()] to see how to handle sex to get correct
#' summary.
#' @param input a `data.frame` or a file or a directory contains copy number profile.
#' @param pattern an optional regular expression used to select part of files if
#' `input` is a directory, more detail please see [list.files()] function.
#' @param ignore_case logical. Should pattern-matching be case-insensitive?
#' @param seg_cols four strings used to specify chromosome, start position,
#'  end position and copy number value in `input`, respectively.
#'  Default use names from ABSOLUTE calling result.
#' @param samp_col a character used to specify the sample column name. If `input`
#' is a directory and cannot find `samp_col`, sample names will use file names
#' (set this parameter to `NULL` is recommended in this case).
#' @param add_loh if `TRUE`, add LOH labels to segments. **NOTE** a column
#' 'minor_cn' must exist to indicate minor allele copy number value.
#' Sex chromosome will not be labeled.
#' @param loh_min_len The length cut-off for labeling a segment as 'LOH'.
#' Default is `10Kb`.
#' @param loh_min_frac When `join_adj_seg` set to `TRUE`, only the length fraction
#' of LOH region is larger than this value will be labeled as 'LOH'.
#' Default is 30%.
#' @param join_adj_seg if `TRUE` (default), join adjacent segments with
#' same copy number value. This is helpful for precisely count the number of breakpoint.
#' When set `use_all=TRUE`, the mean function will be applied to extra numeric columns
#' and unique string columns will be pasted by comma for joined records.
#' @param skip_annotation if `TRUE`, skip annotation step, it may affect some analysis
#' and visualization functionality, but speed up reading data.
#' @param use_all default is `FALSE`. If `True`, use all columns from raw input.
#' @param min_segnum minimal number of copy number segments within a sample.
#' @param max_copynumber bigger copy number within a sample will be reset to this value.
#' @param genome_build genome build version, should be 'hg19', 'hg38', 'mm9' or 'mm10'.
#' @param genome_measure default is 'called', can be 'wg' or 'called'.
#' Set 'called' will use called segments size to compute total size for CNA burden calculation,
#' this option is useful for WES and target sequencing.
#' Set 'wg' will use autosome size from genome build, this option is useful for WGS, SNP etc..
#' @param complement if `TRUE`, complement chromosome (except 'Y') does not show in input data
#' with normal copy 2.
#' @param ... other parameters pass to [data.table::fread()]
#' @author Shixiang Wang <w_shixiang@163.com>
#' @return a [CopyNumber] object.
#' @export
#' @examples
#' # Load toy dataset of absolute copynumber profile
#' load(system.file("extdata", "toy_segTab.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#'
#' \donttest{
#' cn <- read_copynumber(segTabs,
#'   seg_cols = c("chromosome", "start", "end", "segVal"),
#'   genome_build = "hg19", complement = FALSE
#' )
#' cn
#' cn_subset <- subset(cn, sample == "TCGA-DF-A2KN-01A-11D-A17U-01")
#'
#' # Add LOH
#' set.seed(1234)
#' segTabs$minor_cn <- sample(c(0, 1), size = nrow(segTabs), replace = TRUE)
#' cn <- read_copynumber(segTabs,
#'   seg_cols = c("chromosome", "start", "end", "segVal"),
#'   genome_measure = "wg", complement = TRUE, add_loh = TRUE
#' )
#' # Use tally method "S" (Steele et al.)
#' tally_s <- sig_tally(cn, method = "S")
#'
#' tab_file <- system.file("extdata", "metastatic_tumor.segtab.txt",
#'   package = "sigminer", mustWork = TRUE
#' )
#' cn2 <- read_copynumber(tab_file)
#' cn2
#' }
#' @testexamples
#' expect_s4_class(cn, "CopyNumber")
#' expect_s4_class(cn_subset, "CopyNumber")
#' expect_s4_class(cn2, "CopyNumber")
#' expect_is(tally_s, "list")
#' @seealso [read_maf] for reading mutation data to [MAF] object.
read_copynumber <- function(input,
                            pattern = NULL,
                            ignore_case = FALSE,
                            seg_cols = c("Chromosome", "Start.bp", "End.bp", "modal_cn"),
                            samp_col = "sample",
                            add_loh = FALSE,
                            loh_min_len = 1e4,
                            loh_min_frac = 0.05,
                            join_adj_seg = TRUE,
                            skip_annotation = FALSE,
                            use_all = add_loh,
                            min_segnum = 0L,
                            max_copynumber = 20L,
                            genome_build = c("hg19", "hg38", "mm10", "mm9"),
                            genome_measure = c("called", "wg"),
                            complement = FALSE,
                            ...) {
  stopifnot(
    is.character(samp_col),
    length(samp_col) == 1,
    min_segnum >= 0
  )

  timer <- Sys.time()
  send_info("Started.")
  on.exit(send_elapsed_time(timer))

  #--- match genome build
  genome_build <- match.arg(genome_build)
  genome_measure <- match.arg(genome_measure)

  send_info("Genome build  : ", genome_build, ".")
  send_info("Genome measure: ", genome_measure, ".")

  if (add_loh) {
    use_all <- TRUE
    send_info("When add_loh is TRUE, use_all is forced to TRUE.
              Please drop columns you don't want to keep before reading.")
  }

  # get chromosome lengths
  if (genome_build %in% c("mm10", "mm9")) {
    valid_chr <- c(paste0("chr", 1:19), "chrX", "chrY")
  } else {
    valid_chr <- c(paste0("chr", 1:22), "chrX", "chrY")
  }
  chrlen <- get_genome_annotation(
    data_type = "chr_size",
    chrs = valid_chr,
    genome_build = genome_build
  )
  data.table::setDT(chrlen)
  send_success("Chromosome size database for build obtained.")

  send_info("Reading input.")
  if (tryCatch(dir.exists(input), error = function(e) FALSE)) {
    send_success("A directory as input detected.")
    if (length(input) != 1) {
      send_stop("Only can take one directory as input!")
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
      send_stop("No files exist, please check!")
    }
    files_path <- file.path(input, files)
    data_list <- list()
    dropoff_list <- list()

    # read files
    sb <- cli::cli_status("{symbol$arrow_right} About to read files.")
    Sys.sleep(0.5)
    for (i in seq_along(files_path)) {
      cli::cli_status_update(id = sb, "{symbol$arrow_right} Reading file {files_path[i]}.")

      temp <- data.table::fread(files_path[i], ...)
      if (!all(seg_cols %in% colnames(temp))) {
        send_stop("Not all seg_cols are in file, please check.")
      }

      if (length(samp_col %in% colnames(temp)) == 0 | !(samp_col %in% colnames(temp))) {
        cli::cli_status_update(id = sb, "{symbol$arrow_right} Select file names as sample names.")
        temp[, "sample"] <- files[i]
        sample_col <- "sample"
      }

      tempName <- unique(temp[[samp_col]])
      if (length(tempName) > 1) {
        send_stop("When input is a directory, a file can only contain one sample.")
      }

      # set column order
      data.table::setcolorder(temp, neworder = c(seg_cols, samp_col))
      new_cols <- c("chromosome", "start", "end", "segVal", "sample")
      colnames(temp)[1:5] <- new_cols

      if (any(is.na(temp$segVal))) {
        temp <- temp[!is.na(temp$segVal)]
      }

      # unify chromosome column
      cli::cli_status_update(id = sb, "{symbol$arrow_right} Checking chromosome names.")
      temp[, chromosome := sub(
        pattern = "chr",
        replacement = "chr",
        x = as.character(chromosome),
        ignore.case = TRUE
      )]
      temp$chromosome <- ifelse(startsWith(temp$chromosome, "chr"),
        temp$chromosome,
        paste0("chr", temp$chromosome)
      )
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
        cli::cli_status_update(
          id = sb,
          "{symbol$arrow_right} Fill value 2 (normal copy) to uncalled chromosomes."
        )
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
          comp_df[, setdiff(
            colnames(comp_df),
            c("chromosome", "start", "end", "segVal", "sample")
          ) := NA]
          temp <- rbind(temp, comp_df, fill = TRUE)
        }
      }

      if (!use_all) temp <- temp[, new_cols, with = FALSE]
      if (nrow(temp) < min_segnum) {
        dropoff_list[[tempName]] <- temp
      } else {
        data_list[[tempName]] <- temp
      }
    }
    cli::cli_status_clear(sb)

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
      send_success("A file as input detected.")
      if (length(input) > 1) {
        send_stop("Muliple files are not a valid input, please use directory as input.")
      }
      if (!file.exists(input)) {
        send_stop("Input file not exists.")
      }

      temp <- data.table::fread(input, ...)
    } else {
      send_success("A data frame as input detected.")
      temp <- data.table::as.data.table(input)
    }

    #--- Check column names
    if (is.null(samp_col)) {
      send_stop("'samp_col' parameter must set!")
    }
    if (!all(seg_cols %in% colnames(temp))) {
      send_stop("Not all seg_cols are in file, please check.")
    }
    if (!(samp_col %in% colnames(temp))) {
      send_stop("Column ", samp_col, " does not exist.")
    }
    send_success("Column names checked.")

    #--- Set column order
    data.table::setcolorder(temp, neworder = c(seg_cols, samp_col))
    new_cols <- c("chromosome", "start", "end", "segVal", "sample")
    colnames(temp)[1:5] <- new_cols

    send_success("Column order set.")

    if (is.factor(temp$sample)) {
      temp$sample <- as.character(temp$sample)
    }

    if (any(is.na(temp$segVal))) {
      temp <- temp[!is.na(temp$segVal)]
      send_success("Rows with NA copy number removed.")
    }

    # unify chromosome column
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

    send_success("Chromosomes unified.")

    if (complement) {
      # complement value 2 (normal copy) to chromosome not called
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
          comp <- rbind(comp, comp_df, fill = TRUE)
        }
      }
      comp[, setdiff(
        colnames(comp),
        c("chromosome", "start", "end", "segVal", "sample")
      ) := NA]
      temp <- rbind(temp, comp, fill = TRUE)
      send_success("Value 2 (normal copy) filled to uncalled chromosomes.")
    }

    if (!use_all) temp <- temp[, new_cols, with = FALSE]

    dropoff_samples <- temp[, .N, by = .(sample)][N < min_segnum][["sample"]]
    keep_samples <- base::setdiff(unique(temp[["sample"]]), dropoff_samples)

    data_df <- temp[sample %in% keep_samples]
    dropoff_df <- temp[sample %in% dropoff_samples]
  } else {
    send_stop("Invalid input.")
  }

  send_success("Data imported.")

  if (!all(data_df$chromosome %in% valid_chr)) {
    data_drop <- data_df[!chromosome %in% valid_chr]
    if (nrow(dropoff_df) >= 1) {
      dropoff_df <- base::rbind(dropoff_df, data_drop)
    } else {
      dropoff_df <- data_drop
    }

    data_df <- data_df[chromosome %in% valid_chr]
    send_success("Some invalid segments (not 1:22 and X, Y) dropped.")
  }

  send_info("Segments info:")
  send_info("    Keep - ", nrow(data_df))
  send_info("    Drop - ", nrow(dropoff_df))

  # reset copy number for high copy number segments
  data_df$segVal[data_df$segVal > max_copynumber] <- max_copynumber
  # make sure seg copy number value is integer
  data_df[["segVal"]] <- as.integer(round(data_df[["segVal"]]))
  # make sure position is numeric
  data_df$start <- as.numeric(data_df$start)
  data_df$end <- as.numeric(data_df$end)

  data.table::setorderv(data_df, c("sample", "chromosome", "start"))
  send_success("Segments sorted.")

  if (add_loh) {
    send_info("Adding LOH labels...")
    if (!"minor_cn" %in% colnames(data_df)) {
      send_stop("When you want to add LOH infor, a column named as 'minor_cn' should exist!")
    }
    # Make sure the minor_cn is really minor copy
    data_df$minor_cn <- pmin(
      data_df$segVal - data_df$minor_cn,
      data_df$minor_cn
    )
    # Determine a LOH label
    data_df$loh <- data_df$segVal >= 1 & data_df$minor_cn == 0 &
      (data_df$end - data_df$start > loh_min_len - 1)
    # We don't label sex chromosomes
    data_df[data_df$chromosome %in% c("chrX", "chrY")]$loh <- FALSE
  }

  if (join_adj_seg) {
    send_info("Joining adjacent segments with same copy number value. Be patient...")
    # When LOH regions have same total copy number values as adjacent
    # regions, only label the segments harbor LOH with minimal length fraction
    data_df <- helper_join_segments2(data_df,
      add_loh = add_loh,
      loh_min_frac = loh_min_frac
    )
    send_success(nrow(data_df), " segments left after joining.")
  } else {
    send_info("Skipped joining adjacent segments with same copy number value.")
  }
  # order by segment start position by each chromosome in each sample
  data.table::setorderv(data_df, c("sample", "chromosome", "start"))
  data.table::setcolorder(data_df, c("chromosome", "start", "end", "segVal", "sample"))

  if ("groups" %in% names(attributes(data_df))) {
    attr(data_df, "groups") <- NULL
  }

  send_success("Segmental table cleaned.")

  if (skip_annotation) {
    annot <- data.table::data.table()
    send_info("Annotation skipped.")
  } else {
    send_info("Annotating.")
    annot <- get_LengthFraction(data_df,
      genome_build = genome_build,
      seg_cols = new_cols[1:4],
      samp_col = new_cols[5]
    )
    send_success("Annotation done.")
  }


  send_info("Summarizing per sample.")
  sum_sample <- get_cnsummary_sample(data_df,
    genome_build = genome_build,
    genome_measure = genome_measure
  )
  send_success("Summarized.")

  send_info("Generating CopyNumber object.")
  res <- CopyNumber(
    data = data_df,
    summary.per.sample = sum_sample,
    genome_build = genome_build,
    genome_measure = genome_measure,
    annotation = annot,
    dropoff.segs = dropoff_df
  )
  send_success("Generated.")

  send_info("Validating object.")
  res <- validate_segTab(res)
  send_success("Done.")
  res
}


# Global variables --------------------------------------------------------

utils::globalVariables(
  c(
    ".",
    "N",
    ".N",
    ".SD",
    "flag",
    "p_start",
    "p_end",
    "q_start",
    "q_end",
    "total_size"
  )
)
