# ========================================================
# internal analysis functions, use get_ as start
# ========================================================

# Get copy number list ----------------------------------------------------

get_cnlist <- function(CopyNumber) {
  if (!inherits(CopyNumber, "CopyNumber")) {
    stop("Input must be a CopyNumber object!")
  }
  res <- split(CopyNumber@data, by = "sample")
  res
}

# Get feature distributions -----------------------------------------------

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

  oplan <- future::plan()
  future::plan("multiprocess", workers = cores)
  on.exit(future::plan(oplan), add = TRUE)


  features <- c(
    "segsize", "bp10MB", "osCN", "bpchrarm",
    "changepoint", "copynumber"
  )

  .get_feature <- function(i) {
    if (i == "segsize") {
      message("Getting segsize...")
      getSegsize(CN_data)
    } else if (i == "bp10MB") {
      message("Getting bp10MB...")
      getBPnum(CN_data, chrlen)
    } else if (i == "osCN") {
      message("Getting osCN...")
      getOscilation(CN_data)
    } else if (i == "bpchrarm") {
      message("Getting bpchrarm...")
      getCentromereDistCounts(CN_data, centromeres, chrlen)
    } else if (i == "changepoint") {
      message("Getting changepoint...")
      getChangepointCN(CN_data)
    } else {
      message("Getting copynumber...")
      getCN(CN_data)
    }
  }

  res <- furrr::future_map(features, .get_feature,
    .progress = TRUE
  )
  res <- res %>% setNames(features)
  res
}


# Get mixture model components --------------------------------------------

get_components <- function(CN_features,
                           seed = 123456,
                           min_comp = 2,
                           max_comp = 15,
                           min_prior = 0.001,
                           model_selection = "BIC",
                           threshold = 0.1,
                           nrep = 3,
                           niter = 1000,
                           cores = 1) {
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

  dist_map <- c("norm", "pois", "pois", "pois", "norm", "norm")
  names(dist_map) <- c("segsize", "bp10MB", "osCN", "bpchrarm", "changepoint", "copynumber")
  apply_fitComponent <- function(CN_feature,
                                   feature_name,
                                   dist_map,
                                   seed = 123456,
                                   min_comp = 2,
                                   max_comp = 10,
                                   min_prior = 0.001,
                                   model_selection = "BIC",
                                   threshold = 0.1,
                                   nrep = 1,
                                   niter = 1000,
                                   cores = 1) {
    dat <- as.numeric(CN_feature[, 2])
    message("Fitting feature: ", feature_name)
    fit <-
      fitComponent(
        dat,
        dist = dist_map[feature_name],
        seed = seed,
        model_selection = model_selection,
        threshold = threshold,
        min_prior = min_prior,
        niter = niter,
        nrep = nrep,
        min_comp = min_comp,
        max_comp = max_comp,
        cores = cores
      )
    fit
  }

  purrr::pmap(list(
    CN_feature = CN_features,
    feature_name = names(CN_features),
    min_comp = if (flag_min) min_comp else rep(min_comp, 6),
    max_comp = if (flag_max) max_comp else rep(max_comp, 6)
  ), apply_fitComponent,
  seed = seed,
  min_prior = min_prior,
  model_selection = model_selection,
  threshold = threshold,
  nrep = nrep,
  niter = niter,
  dist_map = dist_map,
  cores = cores
  )
}



# Get sample-by-component matrix ------------------------------------------

get_matrix <- function(CN_features,
                       all_components = NULL,
                       type = c("probability", "count"),
                       cores = 1) {
  type <- match.arg(type)

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

  feature_orders <- c("bp10MB", "copynumber", "changepoint", "bpchrarm", "osCN", "segsize")
  # Make sure have same order
  CN_features <- CN_features[feature_orders]
  all_components <- all_components[feature_orders]

  oplan <- future::plan()
  future::plan("multiprocess", workers = cores)
  on.exit(future::plan(oplan), add = TRUE)

  full_mat <- furrr::future_pmap(
    list(
      feature = CN_features,
      component = all_components,
      name = feature_orders
    ),
    .f = calculateSumOfEvents,
    type = type,
    .progress = TRUE
  )

  full_mat <- purrr::reduce(full_mat, .f = base::cbind)

  # rownames(full_mat) <- unique(CN_features[["segsize"]][, 1])
  # full_mat[is.na(full_mat)] <- 0
  full_mat
}



# Get copy number length profile ------------------------------------------

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

  cbind(data.table::setDT(segTab)[, colnames(arm_data)[-1] := NULL], data.table::setDT(annot))
}


# Get arm location --------------------------------------------------------

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
    "i",
    "chrom",
    "chromosome",
    "segVal"
  )
)
