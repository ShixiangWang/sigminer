# ========================================================
# internal analysis functions, use get_ as start
# ========================================================

# Get copy number list ----------------------------------------------------

#' Extract copy number profile as list from CopyNumber object
#' @author Shixiang Wang
#' @param CopyNumber a [CopyNumber] object.
#' @return a `list`
#' @export
#' @examples
#' extdata_dir <- system.file("extdata", package = "sigminer", mustWork = TRUE)
#' cp <- read_copynumber(extdata_dir, pattern = "txt", genome_build = "hg19")
#' cn_list <- get_cnlist(cp)
#' @family internal calculation function series

get_cnlist <- function(CopyNumber) {
  if (!inherits(CopyNumber, "CopyNumber")) {
    stop("Input must be a CopyNumber object!")
  }
  res <- split(CopyNumber@data, by = "sample")
  res
}

# Get feature distributions -----------------------------------------------

#' @title  Derive copy number feature distributions
#' @description This function summarise each copy-number profile using a number of different
#' feature distributions: sigment size, breakpoint number (per ten megabase), change-point copy-number,
#' segment copy-number, breakpoint number (per chromosome arm), length of segments with oscilating
#' copy-number.
#' @param CN_data a `list` contains multiple `data.frame`s (recommended),
#' each `data.frame` stores copy-number profile for one sample with 'chromosome', 'start', 'end' and
#' 'segVal' these four necessary columns. Of note, 'segVal' column shoule be absolute copy number values.
#' @param cores number of compute cores to run this task.
#' You can use [parallel::detectCores()] function to check how
#' many cores you can use.
#' @inheritParams read_copynumber
#' @author Geoffrey Macintyre, Shixiang Wang
#' @return a `list` contains six copy number feature distributions.
#' @import foreach
#' @export
#' @examples
#' # Load copy number list
#' load(system.file("extdata", "toy_cnlist.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' \donttest{
#' cn_features <- get_features(cn_list, cores = 1)
#' }
#' @family internal calculation function series

get_features <- function(CN_data,
                         cores = 1,
                         genome_build = c("hg19", "hg38")) {
  genome_build <- match.arg(genome_build)
  # get chromosome lengths and centromere locations
  if (genome_build == "hg19") {
    data("chromsize.hg19",
      package = "sigminer",
      envir = environment()
    )
    data("centromeres.hg19",
      package = "sigminer",
      envir = environment()
    )
    chrlen <- chromsize.hg19
    centromeres <- centromeres.hg19
  } else {
    data("chromsize.hg38",
      package = "sigminer",
      envir = environment()
    )
    data("centromeres.hg38",
      package = "sigminer",
      envir = environment()
    )
    chrlen <- chromsize.hg38
    centromeres <- centromeres.hg38
  }

  # only keep 1:22 and x, y
  chrlen <- chrlen[chrlen$chrom %in% centromeres$chrom, ]
  if (cores > 1) {

    # attachNamespace("foreach")
    # attachNamespace("doParallel")
    doParallel::registerDoParallel(cores = cores)

    temp_list <- foreach::foreach(i = 1:6) %dopar% {
      if (i == 1) {
        list(bp10MB = getBPnum(CN_data, chrlen))
      } else if (i == 2) {
        list(copynumber = getCN(CN_data))
      } else if (i == 3) {
        list(changepoint = getChangepointCN(CN_data))
      } else if (i == 4) {
        list(bpchrarm = getCentromereDistCounts(CN_data, centromeres, chrlen))
      } else if (i == 5) {
        list(osCN = getOscilation(CN_data))
      } else {
        list(segsize = getSegsize(CN_data))
      }
    }
    unlist(temp_list, recursive = FALSE)
  } else {
    bp10MB <- getBPnum(CN_data, chrlen)
    copynumber <- getCN(CN_data)
    changepoint <- getChangepointCN(CN_data)
    bpchrarm <-
      getCentromereDistCounts(CN_data, centromeres, chrlen)
    osCN <- getOscilation(CN_data)
    segsize <- getSegsize(CN_data)
    list(
      segsize = segsize,
      bp10MB = bp10MB,
      osCN = osCN,
      bpchrarm = bpchrarm,
      changepoint = changepoint,
      copynumber = copynumber
    )
  }
}




# Get mixture model components --------------------------------------------

#' @title Fit optimal number of mixture model components
#' @description Apply mixture modelling to breakdown each feature distribution into mixtures
#' of Gaussian or mixtures of Poison distributions using **flexmix** package. The order of
#' features is 'Segment size', 'Breakpoint count per 10 Mb', 'Length of oscillating copy-number chain',
#' 'Breakpoint count per arm',  'Copy number change', 'Absolute copy number'.
#'
#'
#' @param CN_features a `list` generate from [get_features()] function.
#' @param seed seed number.
#' @param min_comp minimal number of components to fit, default is 2.
#' Can also be a vector with length 6, which apply to each feature.
#' @param max_comp maximal number of components to fit, default is 10.
#' Can also be a vector with length 6, which apply to each feature.
#' @param min_prior minimal prior value, default is 0.001.
#' Details about custom setting please refer to **flexmix** package.
#' @param model_selection model selection strategy, default is 'BIC'.
#' Details about custom setting please refer to **flexmix** package.
#' @param nrep number of run times for each value of component,
#' keep only the solution with maximum likelihood.
#' @param niter maximal number of iteration to achive converge.
#' @author Geoffrey Macintyre, Shixiang Wang
#' @return a `list` contain `flexmix` object of copy-number features.
#' @importClassesFrom flexmix FLXcontrol
#' @export
#' @examples
#' \donttest{
#' # Load copy number features
#' load(system.file("extdata", "toy_cn_features.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' cn_components <- get_components(cn_features)
#' }
#' @family internal calculation function series
get_components <- function(CN_features,
                           seed = 123456,
                           min_comp = 2,
                           max_comp = 10,
                           min_prior = 0.001,
                           model_selection = "BIC",
                           nrep = 1,
                           niter = 1000) {
  flag_min <- FALSE
  flag_max <- FALSE

  if (length(min_comp) > 1) {
    if (length(min_comp) != 6) {
      stop("If you use more than 1 values, length must be 6, each value for corresponding feature!")
    }
    flag_min <- TRUE
  }

  if (length(max_comp) > 1) {
    if (length(max_comp) != 6) {
      stop("If you use more than 1 values, length must be 6, each value for corresponding feature!")
    }
    flag_max <- TRUE
  }

  dat <- as.numeric(CN_features[["segsize"]][, 2])
  message("Fit feature: Segment size")
  segsize_mm <-
    fitComponent(
      dat,
      seed = seed,
      model_selection = model_selection,
      min_prior = min_prior,
      niter = niter,
      nrep = nrep,
      min_comp = if (flag_min) min_comp[1] else min_comp,
      max_comp = if (flag_max) max_comp[1] else max_comp
    )

  dat <- as.numeric(CN_features[["bp10MB"]][, 2])
  message("Fit feature: Breakpoint count per 10 Mb")
  bp10MB_mm <-
    fitComponent(
      dat,
      dist = "pois",
      seed = seed,
      model_selection = model_selection,
      min_prior = min_prior,
      niter = niter,
      nrep = nrep,
      min_comp = if (flag_min) min_comp[2] else min_comp,
      max_comp = if (flag_max) max_comp[2] else max_comp
    )

  dat <- as.numeric(CN_features[["osCN"]][, 2])
  message("Fit feature: Length of oscillating copy-number chain")
  osCN_mm <-
    fitComponent(
      dat,
      dist = "pois",
      seed = seed,
      model_selection = model_selection,
      min_prior = min_prior,
      niter = niter,
      nrep = nrep,
      min_comp = if (flag_min) min_comp[3] else min_comp,
      max_comp = if (flag_max) max_comp[3] else max_comp
    )

  dat <- as.numeric(CN_features[["bpchrarm"]][, 2])
  message("Fit feature: Breakpoint count per arm")
  bpchrarm_mm <-
    fitComponent(
      dat,
      dist = "pois",
      seed = seed,
      model_selection = model_selection,
      min_prior = min_prior,
      niter = niter,
      nrep = nrep,
      min_comp = if (flag_min) min_comp[4] else min_comp,
      max_comp = if (flag_max) max_comp[4] else max_comp
    )

  dat <- as.numeric(CN_features[["changepoint"]][, 2])
  message("Fit feature: Copy number change")
  changepoint_mm <-
    fitComponent(
      dat,
      seed = seed,
      model_selection = model_selection,
      min_prior = min_prior,
      niter = niter,
      nrep = nrep,
      min_comp = if (flag_min) min_comp[5] else min_comp,
      max_comp = if (flag_max) max_comp[5] else max_comp
    )

  dat <- as.numeric(CN_features[["copynumber"]][, 2])
  message("Fit feature: Absolute copy number")
  copynumber_mm <-
    fitComponent(
      dat,
      seed = seed,
      model_selection = model_selection,
      nrep = nrep,
      min_comp = if (flag_min) min_comp[6] else min_comp,
      max_comp = if (flag_max) max_comp[6] else max_comp,
      min_prior = min_prior,
      niter = niter
    )

  list(
    segsize = segsize_mm,
    bp10MB = bp10MB_mm,
    osCN = osCN_mm,
    bpchrarm = bpchrarm_mm,
    changepoint = changepoint_mm,
    copynumber = copynumber_mm
  )
}



# Get sample-by-component matrix ------------------------------------------

#' @title Generate a sample-by-component matrix
#' @description Generate a sample-by-component matrix representing the sum of
#' posterior probabilities of each copy-number event being assigned to each component.
#' @param CN_features a `list` contains six copy number feature distributions,
#' obtain this from [get_features()] function.
#' @param all_components a `list` contain `flexmix` object of copy-number features, obtain this
#' from [get_components] function or use pre-compiled components data which come from CNV signature paper
#' https://www.nature.com/articles/s41588-018-0179-8 (set this parameter as `NULL`).
#' @inheritParams get_features
#' @param rowIter step size of iteration for rows of ech CNV feature.
#' @author Geoffrey Macintyre, Shixiang Wang
#' @return a numeric sample-by-component `matrix`
#' @importFrom utils data download.file str
#' @export
#' @examples
#' # Load copy number components
#' load(system.file("extdata", "toy_cn_components.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' # Load copy number features
#' load(system.file("extdata", "toy_cn_features.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' 
#' cn_matrix <- get_matrix(cn_features, cn_components)
#' @family internal calculation function series
get_matrix <- function(CN_features,
                       all_components = NULL,
                       cores = 1,
                       rowIter = 1000) {
  if (is.null(all_components)) {
    message(
      "About reference components\n   more detail please see https://github.com/ShixiangWang/absoluteCNVdata"
    )
    if (!file.exists("Nat_Gen_component_parameters.rds")) {
      message(
        "Nat_Gen_component_parameters.rds doesnot exist, will download reference components."
      )
      download.file(
        url = "https://github.com/ShixiangWang/absoluteCNVdata/raw/master/component_parameters.rds",
        destfile = file.path(tempdir(), "Nat_Gen_component_parameters.rds")
      )
    }
    all_components <-
      readRDS(file.path(tempdir(), "Nat_Gen_component_parameters.rds"))
  }

  full_mat <- cbind(
    calculateSumOfPosteriors(CN_features[["bp10MB"]],
      all_components[["bp10MB"]],
      "bp10MB",
      cores = cores
    ),
    calculateSumOfPosteriors(CN_features[["copynumber"]],
      all_components[["copynumber"]],
      "copynumber",
      cores = cores
    ),
    calculateSumOfPosteriors(CN_features[["changepoint"]],
      all_components[["changepoint"]],
      "changepoint",
      cores = cores
    ),
    calculateSumOfPosteriors(CN_features[["bpchrarm"]],
      all_components[["bpchrarm"]],
      "bpchrarm",
      cores = cores
    ),
    calculateSumOfPosteriors(CN_features[["osCN"]],
      all_components[["osCN"]],
      "osCN",
      cores = cores
    ),
    calculateSumOfPosteriors(CN_features[["segsize"]],
      all_components[["segsize"]],
      "segsize",
      cores = cores
    )
  )

  rownames(full_mat) <- unique(CN_features[["segsize"]][, 1])
  full_mat[is.na(full_mat)] <- 0
  full_mat
}



# Get copy number length profile ------------------------------------------

#' Calculate length fraction profile of copy number
#'
#' @inheritParams get_features
#' @inheritParams read_copynumber
#' @param CN_data a `data.frame` with 'chromosome', 'start', 'end' and 'segVal'
#' (optinal) and 'sample' these five columns
#' (specify column names using `seg_cols` and `samp_cols` options)
#' or a `list` contains multiple data.frames,
#' each `data.frame` stores copy-number profile for one sample with
#' 'chromosome', 'start', 'end' and 'segVal' (optional) these four columns.
#' If 'sample' column is not specified, will try using name of each `data.frame`.
#' @param samp_col a character used to specify the sample column name.
#' @author Shixiang Wang <w_shixiang@163.com>
#' @return a data table
#' @importFrom dplyr progress_estimated
#' @importFrom purrr pmap_df
#' @export
#' @examples
#' # Load copy number list
#' load(system.file("extdata", "toy_cnlist.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' annot <- get_LengthFraction(cn_list, seg_cols = c("chromosome", "start", "end", "segVal"))
#' @family internal calculation function series
get_LengthFraction <- function(CN_data,
                               genome_build = c("hg19", "hg38"),
                               seg_cols = c("Chromosome", "Start.bp", "End.bp", "modal_cn"),
                               samp_col = "sample") {
  stopifnot(is.list(CN_data) | is.data.frame(CN_data))
  genome_build <- match.arg(genome_build)

  if (inherits(CN_data, "list")) {
    segTab <- data.table::rbindlist(CN_data, use.names = TRUE, fill = TRUE)
    # segTab = base::Reduce(rbind, CN_data)
    if (!samp_col %in% colnames(segTab)) {
      segTab$sample <- base::rep(
        x = names(CN_data),
        times = sapply(CN_data, function(x)
          nrow(x))
      )
      samp_col <- "sample"
    }
  } else {
    segTab <- CN_data
  }

  if (inherits(segTab, "data.table")) {
    segTab <- segTab[, c(seg_cols, samp_col), with = FALSE]
  } else {
    segTab <- segTab[, c(seg_cols, samp_col)]
  }

  if (ncol(segTab) == 5) {
    colnames(segTab) <- c("chromosome", "start", "end", "segVal", "sample")
  } else if (ncol(segTab) == 4) {
    colnames(segTab) <- c("chromosome", "start", "end", "sample")
  } else {
    stop(
      "If input is a data.frame, must have 4 necessary columns (chr, start, end, sample) and 1 optional column (segVal)."
    )
  }

  data.table::setDT(segTab)
  segTab$start <- as.integer(segTab$start)
  segTab$end <- as.integer(segTab$end)
  # unify chromosome column
  segTab$chromosome <- as.character(segTab$chromosome)
  segTab$chromosome <- sub(
    pattern = "chr",
    replacement = "chr",
    x = segTab$chromosome,
    ignore.case = TRUE
  )
  if (any(!grepl("chr", segTab$chromosome))) {
    segTab$chromosome[!grepl("chr", segTab$chromosome)] <- paste0("chr", segTab$chromosome[!grepl("chr", segTab$chromosome)])
  }

  valid_chr <- c(paste0("chr", 1:22), "chrX", "chrY")
  if (!all(segTab$chromosome %in% valid_chr)) {
    message("Filter some invalid segments... (not as 1:22 and X, Y)")

    segTab <- segTab[valid_chr, on = "chromosome"]
  }

  arm_data <- get_ArmLocation(genome_build)
  data.table::setDT(arm_data)

  segTab <- dplyr::left_join(segTab, arm_data, by = c("chromosome" = "chrom"))

  .annot_fun <- function(chrom, start, end, p_start, p_end, p_length, q_start, q_end, q_length, total_size) {
    if (end <= p_end & start >= p_start) {
      location <- paste0(sub("chr", "", chrom), "p")
      annotation <- "short arm"
      fraction <- (end - start + 1) / (p_end - p_start + 1)
    } else if (end <= q_end &
      start >= q_start) {
      location <- paste0(sub("chr", "", chrom), "q")
      annotation <- "long arm"
      fraction <- (end - start + 1) / (q_end - q_start + 1)
    } else if (start >= p_start &
      start <= p_end &
      end >= q_start & end <= q_end) {
      location <- paste0(sub("chr", "", chrom), "pq") # across p and q arm
      annotation <- "across short and long arm"
      fraction <- 2 * ((end - start + 1) / total_size)
    } else if (start < p_end & end < q_start) {
      location <- paste0(sub("chr", "", chrom), "p")
      annotation <- "short arm intersect with centromere region"
      # only calculate region does not intersect
      fraction <- (end - start + 1 - (end - p_end)) / (p_end - p_start + 1)
    } else if (start > p_end &
      start < q_start & end > q_start) {
      location <- paste0(sub("chr", "", chrom), "q")
      annotation <- "long arm intersect with centromere region"
      # only calculate region does not intersect
      fraction <- (end - start + 1 - (start - q_start)) / (q_end - q_start + 1)
    } else {
      location <- paste0(sub("chr", "", chrom), "pq") # suppose as pq
      annotation <- "segment locate in centromere region"
      fraction <- 2 * ((end - start + 1) / total_size)
    }

    dplyr::tibble(location = location, annotation = annotation, fraction = fraction)
  }

  annot_fun <- function(chrom, start, end, p_start, p_end, p_length, q_start,
                          q_end, q_length, total_size, .pb = NULL) {
    if (.pb$i < .pb$n) .pb$tick()$print()
    .annot_fun(
      chrom, start, end, p_start, p_end, p_length, q_start,
      q_end, q_length, total_size
    )
  }

  pb <- progress_estimated(nrow(segTab), 0)

  annot <- purrr::pmap_df(
    list(
      chrom = segTab$chromosome,
      start = segTab$start,
      end = segTab$end,
      p_start = segTab$p_start,
      p_end = segTab$p_end,
      p_length = segTab$p_length,
      q_start = segTab$q_start,
      q_end = segTab$q_end,
      q_length = segTab$q_length,
      total_size = segTab$total_size
    ), annot_fun,
    .pb = pb
  )


  # location = vector("character", nrow(segTab))
  # annotation = vector("character", nrow(segTab))
  # fraction = vector("numeric", nrow(segTab))
  #
  # p = dplyr::progress_estimated(nrow(segTab))
  #
  # for (i in 1:nrow(segTab)) {
  #   # locate chromosome
  #   arm_loc = arm_data[segTab[["chromosome"]][i], on = "chrom"]
  #
  #   y = c(segTab$start[i], segTab$end[i])
  #   if (y[2] <= arm_loc$p_end & y[1] >= arm_loc$p_start) {
  #     location[i] = paste0(sub("chr", "", arm_loc$chrom), "p")
  #     annotation[i] = "short arm"
  #     fraction[i] = (y[2] - y[1] + 1) / (arm_loc$p_end - arm_loc$p_start + 1)
  #   } else if (y[2] <= arm_loc$q_end &
  #              y[1] >= arm_loc$q_start) {
  #     location[i] = paste0(sub("chr", "", arm_loc$chrom), "q")
  #     annotation[i] = "long arm"
  #     fraction[i] = (y[2] - y[1] + 1) / (arm_loc$q_end - arm_loc$q_start + 1)
  #   } else if (y[1] >= arm_loc$p_start &
  #              y[1] <= arm_loc$p_end &
  #              y[2] >= arm_loc$q_start & y[2] <= arm_loc$q_end) {
  #     location[i] = paste0(sub("chr", "", arm_loc$chrom), "pq") # across p and q arm
  #     annotation[i] = "across short and long arm"
  #     fraction[i] = 2 * ((y[2] - y[1] + 1) / arm_loc$total_size)
  #   } else if (y[1] < arm_loc$p_end & y[2] < arm_loc$q_start) {
  #     location[i] = paste0(sub("chr", "", arm_loc$chrom), "p")
  #     annotation[i] = "short arm intersect with centromere region"
  #     # only calculate region does not intersect
  #     fraction[i] = (y[2] - y[1] + 1 - (y[2] - arm_loc$p_end)) / (arm_loc$p_end - arm_loc$p_start + 1)
  #   } else if (y[1] > arm_loc$p_end &
  #              y[1] < arm_loc$q_start & y[2] > arm_loc$q_start) {
  #     location[i] = paste0(sub("chr", "", arm_loc$chrom), "q")
  #     annotation[i] = "long arm intersect with centromere region"
  #     # only calculate region does not intersect
  #     fraction[i] = (y[2] - y[1] + 1 - (y[1] - arm_loc$q_start)) / (arm_loc$q_end - arm_loc$q_start + 1)
  #   } else {
  #     location[i] = paste0(sub("chr", "", arm_loc$chrom), "pq") # suppose as pq
  #     annotation[i] = "segment locate in centromere region"
  #     fraction[i] = 2 * ((y[2] - y[1] + 1) / arm_loc$total_size)
  #   }
  #
  #   p$tick()$print()
  # }

  cbind(data.table::setDT(segTab)[, colnames(arm_data)[-1] := NULL], data.table::setDT(annot))
}


# Get arm location --------------------------------------------------------

#' Get chromosome arm location
#' @inheritParams read_copynumber
#' @export
#' @examples
#' hg19_arm <- get_ArmLocation("hg19")
#' hg38_arm <- get_ArmLocation("hg38")
#' @family internal calculation function series

get_ArmLocation <- function(genome_build = c("hg19", "hg38")) {
  genome_build <- match.arg(genome_build)
  # get chromosome lengths and centromere locations
  if (genome_build == "hg19") {
    data("chromsize.hg19",
      package = "sigminer",
      envir = environment()
    )
    data("centromeres.hg19",
      package = "sigminer",
      envir = environment()
    )
    chrlen <- chromsize.hg19
    centromeres <- centromeres.hg19
  } else {
    data("chromsize.hg38",
      package = "sigminer",
      envir = environment()
    )
    data("centromeres.hg38",
      package = "sigminer",
      envir = environment()
    )
    chrlen <- chromsize.hg38
    centromeres <- centromeres.hg38
  }

  # only keep 1:22 and x, y
  chrlen <- chrlen[chrlen$chrom %in% centromeres$chrom, ]

  # sort
  chrlen <- chrlen[order(chrlen$chrom), ]
  centromeres <- centromeres[order(centromeres$chrom), ]

  # compute and get results
  res <- data.frame(
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

  i <- 1
  for (chr in chrlen$chrom) {
    chrom <- chr

    # p
    p_start <- 1
    p_end <- centromeres$left.base[centromeres$chrom == chr]
    p_length <- p_end - p_start + 1

    # q
    q_start <- centromeres$right.base[centromeres$chrom == chr]
    q_end <- chrlen$size[chrlen$chrom == chr]
    q_length <- q_end - q_start + 1

    total_size <- chrlen$size[chrlen$chrom == chr]

    res[i, 1] <- as.character(chrom)
    res[i, 2:8] <- c(
      p_start,
      p_end,
      p_length,
      q_start,
      q_end,
      q_length,
      total_size
    )
    i <- i + 1
  }

  res
}



# Get summary of copy number variation per sample ------------------------------------

#' Get summary of copy number variation per sample
#'
#' Include number of CNV segments, CNA burden, number of CNV amplification segments,
#' number of CNV deletion segments etc..
#'
#' CNA burden, a simple metric of CNA level defined as the percent of the
#' autosomal tumor genome bearing CNAs, could be used as an informative measure of CNA.
#'
#' @inheritParams read_copynumber
#' @param segTab a `data.frame` with 'chromosome', 'start', 'end' and 'segVal'
#' and 'sample' these five ordered columns. 'chromosome' should have prefix "chr".
#' @param min_seg_len minimal length of CNV segment for CNA burden calculation,
#' default is 1000. (!NOT implement NOW!)
#' @references Hieronymus, Haley, et al. "Copy number alteration burden predicts prostate cancer relapse." Proceedings of the National Academy of Sciences 111.30 (2014): 11139-11144.
#' @author Shixiang Wang <w_shixiang@163.com>
#' @return a data table
#' @export
#' @examples
#' \donttest{
#' load(system.file("extdata", "example_cn_list.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#' segTabs <- data.table::rbindlist(tcga_segTabs, idcol = "sample")
#' segTabs$chromosome <- paste0("chr", segTabs$chromosome)
#' samp_sum <- get_cnsummary_sample(segTabs[, c(2:5, 1)])
#' }
#' @family internal calculation function series
get_cnsummary_sample <- function(segTab, genome_build = c("hg19", "hg38"),
                                 genome_measure = c("called", "wg"),
                                 min_seg_len = 1000L) {
  genome_build <- match.arg(genome_build)
  genome_measure <- match.arg(genome_measure)

  segTab <- segTab[, 1:5]
  if (ncol(segTab) == 5) {
    colnames(segTab) <- c("chromosome", "start", "end", "segVal", "sample")
  } else {
    stop(
      "Input must have 5 ordered columns (chr, start, end, segVal, sample)."
    )
  }

  data.table::setDT(segTab)
  segTab$start <- as.integer(segTab$start)
  segTab$end <- as.integer(segTab$end)

  autosome <- paste0("chr", 1:22)

  if (genome_measure == "wg") {
    if (genome_build == "hg19") {
      data("chromsize.hg19",
        package = "sigminer",
        envir = environment()
      )
      chrlen <- chromsize.hg19
    } else {
      data("chromsize.hg38",
        package = "sigminer",
        envir = environment()
      )
      chrlen <- chromsize.hg38
    }

    chrlen <- chrlen[chrlen[["chrom"]] %in% paste0("chr", 1:22), ]
    total_size <- sum(chrlen[["size"]])

    seg_summary <- segTab %>%
      dplyr::group_by(sample) %>%
      dplyr::summarise(
        n_of_cnv = sum(segVal != 2),
        n_of_amp = sum(segVal > 2),
        n_of_del = sum(segVal < 2),
        cna_burden = sum(end[(segVal != 2) & (chromosome %in% autosome)] - start[(segVal != 2) & (chromosome %in% autosome)] + 1) / total_size
      ) %>%
      data.table::as.data.table()
  } else {
    seg_summary <- segTab %>%
      dplyr::group_by(sample) %>%
      dplyr::summarise(
        n_of_cnv = sum(segVal != 2),
        n_of_amp = sum(segVal > 2),
        n_of_del = sum(segVal < 2),
        cna_burden = sum(end[(segVal != 2) & (chromosome %in% autosome)] - start[(segVal != 2) & (chromosome %in% autosome)] + 1) / sum(end[chromosome %in% autosome] - start[chromosome %in% autosome] + 1)
      ) %>%
      data.table::as.data.table()
  }
  seg_summary
}



# Global variables --------------------------------------------------------


utils::globalVariables(
  c(
    "centromeres.hg19",
    "centromeres.hg38",
    "chromsize.hg19",
    "chromsize.hg38",
    ".",
    "i",
    "N",
    "chrom",
    "chromosome"
  )
)
