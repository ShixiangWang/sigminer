#========================================================
# Function: analysis helper functions, use get_ as start
#========================================================

#' Calculate length fraction profile of copy number
#'
#' @inheritParams cnv_derivefeatures
#' @inheritParams read_copynumber
#' @author Shixiang Wang <w_shixiang@163.com>
#' @return a data frame
#' @export
#' @family CNV analysis functions
get_LengthFraction = function(CN_data,
                              genome_build = c("hg19", "hg38"),
                              seg_cols = c("Chromosome", "Start.bp", "End.bp", "modal_cn"),
                              samp_col = "sample") {
  stopifnot(is.list(CN_data) | is.data.frame(CN_data))
  genome_build = match.arg(genome_build)

  if (inherits(CN_data, "list")) {
    segTab = base::Reduce(rbind, CN_data)
    segTab$sample = base::rep(x = names(CN_data),
                              times = sapply(CN_data, function(x)
                                nrow(x)))
  } else {
    segTab = CN_data[, c(seg_cols, samp_col)]
    if (ncol(segTab) == 5) {
      colnames(segTab) = c("chromosome", "start", "end", "segVal", "sample")
    } else if (ncol(segTab) == 4) {
      colnames(segTab) = c("chromosome", "start", "end", "sample")
    } else {
      stop(
        "If input is a data.frame, must have 4 necessary columns (chr, start, end, sample) and 1 optional column (segVal)."
      )
    }
  }

  # unify chromosome column
  segTab$chromosome = as.character(segTab$chromosome)
  segTab$chromosome = sub(
    pattern = "chr",
    replacement = "chr",
    x = segTab$chromosome,
    ignore.case = TRUE
  )
  if (any(!grepl("chr", segTab$chromosome))) {
    segTab$chromosome[!grepl("chr", segTab$chromosome)] = paste0("chr", segTab$chromosome[!grepl("chr", segTab$chromosome)])
  }
  if (any(grepl("chr23", segTab$chromosome))) {
    warning("'23' is not a supported chromosome, related rows will be discarded.")
    segTab = segTab[!grepl("chr23", segTab$chromosome),]
  }

  valid_chr = c(paste0("chr", 1:22), "chrX", "chrY")
  if (!all(segTab$chromosome %in% valid_chr)) {
    message("Filter some invalid segments... (not as 1:22 and X, Y)")
    segTab = base::subset(segTab, chromosome %in% valid_chr)
  }

  arm_data = get_ArmLocation(genome_build)
  # logical operation

  assign_df = data.frame(
    location = vector("character", nrow(segTab)),
    annotation = vector("character", nrow(segTab)),
    fraction = vector("numeric", nrow(segTab)),
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(segTab)) {
    x = segTab[i,]

    # locate chromosome
    arm_loc = base::subset(arm_data, chrom == x$chromosome)

    y = c(x$start, x$end)
    if (y[2] <= arm_loc$p_end & y[1] >= arm_loc$p_start) {
      location = paste0(sub("chr", "", arm_loc$chrom), "p")
      annotation = "short arm"
      fraction = (y[2] - y[1] + 1) / (arm_loc$p_end - arm_loc$p_start + 1)
    } else if (y[2] <= arm_loc$q_end &
               y[1] >= arm_loc$q_start) {
      location = paste0(sub("chr", "", arm_loc$chrom), "q")
      annotation = "long arm"
      fraction = (y[2] - y[1] + 1) / (arm_loc$q_end - arm_loc$q_start + 1)
    } else if (y[1] >= arm_loc$p_start &
               y[1] <= arm_loc$p_end &
               y[2] >= arm_loc$q_start & y[2] <= arm_loc$q_end) {
      location = paste0(sub("chr", "", arm_loc$chrom), "pq") # across p and q arm
      annotation = "across short and long arm"
      fraction = 2 * ((y[2] - y[1] + 1) / arm_loc$total_size)
    } else if (y[1] < arm_loc$p_end & y[2] < arm_loc$q_start) {
      location = paste0(sub("chr", "", arm_loc$chrom), "p")
      annotation = "short arm intersect with centromere region"
      # only calculate region does not intersect
      fraction = (y[2] - y[1] + 1 - (y[2] - arm_loc$p_end)) / (arm_loc$p_end - arm_loc$p_start + 1)
    } else if (y[1] > arm_loc$p_end &
               y[1] < arm_loc$q_start & y[2] > arm_loc$q_start) {
      location = paste0(sub("chr", "", arm_loc$chrom), "q")
      annotation = "long arm intersect with centromere region"
      # only calculate region does not intersect
      fraction = (y[2] - y[1] + 1 - (y[1] - arm_loc$q_start)) / (arm_loc$q_end - arm_loc$q_start + 1)
    } else {
      location = "Unknown"
      annotation = "Unknow segment, locate in centromere region"
      fraction = NA_real_
    }

    assign_df[i,] = c(location, annotation, fraction)
    #assign_df = base::rbind(assign_df, c(location, annotation, percentage))
  }
  res = base::cbind(segTab, assign_df)
  res$fraction = as.numeric(res$fraction)
  res
}


get_ArmLocation = function(genome_build = c("hg19", "hg38")) {
  # get chromosome lengths and centromere locations
  if (genome_build == "hg19") {
    data("chromsize.hg19",
         package = "VSHunter",
         envir = environment())
    data("centromeres.hg19",
         package = "VSHunter",
         envir = environment())
    chrlen = chromsize.hg19
    centromeres = centromeres.hg19
  } else {
    data("chromsize.hg38",
         package = "VSHunter",
         envir = environment())
    data("centromeres.hg38",
         package = "VSHunter",
         envir = environment())
    chrlen = chromsize.hg38
    centromeres = centromeres.hg38
  }

  # only keep 1:22 and x, y
  chrlen = chrlen[chrlen$chrom %in% centromeres$chrom,]

  # sort
  chrlen = chrlen[order(chrlen$chrom),]
  centromeres = centromeres[order(centromeres$chrom),]

  # compute and get results
  res = data.frame(
    chrom = vector(mode = "character", length = 24),
    p_start = vector("integer", length = 24),
    p_end = vector("integer", length = 24),
    p_length = vector("integer", length = 24),
    q_start = vector("integer", length = 24),
    q_end = vector("integer", length = 24),
    q_length = vector("integer", length = 24),
    total_size = vector("integer", length = 24),
    stringsAsFactors = FALSE
  )

  i = 1
  for (chr in chrlen$chrom) {
    chrom = chr

    # p
    p_start = 1
    p_end = centromeres$left.base[centromeres$chrom == chr]
    p_length = p_end - p_start + 1

    # q
    q_start = centromeres$right.base[centromeres$chrom == chr]
    q_end = chrlen$size[chrlen$chrom == chr]
    q_length = q_end - q_start + 1

    total_size = chrlen$size[chrlen$chrom == chr]

    res[i, 1] = as.character(chrom)
    res[i, 2:8] = c(p_start,
                    p_end,
                    p_length,
                    q_start,
                    q_end,
                    q_length,
                    total_size)
    i = i + 1
  }

  res
}
