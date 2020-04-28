# ========================================================
# internal analysis functions, use get_ as start
# ========================================================

# Get copy number list ----------------------------------------------------

get_cnlist <- function(CopyNumber, ignore_chrs = NULL) {
  if (!inherits(CopyNumber, "CopyNumber")) {
    stop("Input must be a CopyNumber object!")
  }
  data <- data.table::copy(CopyNumber@data)
  if (!is.null(ignore_chrs)) {
    chrs_exist <- ignore_chrs %in% unique(data$chromosome)
    if (!any(chrs_exist)) {
      send_info(
        "No chromosome names called ",
        paste(ignore_chrs, collapse = ","), " found in data, skipping filter."
      )
    } else {
      send_info("Filtering out segments in ", paste(ignore_chrs[chrs_exist], collapse = ","))
      data <- data[!chromosome %in% ignore_chrs[chrs_exist]]
    }
  }
  res <- split(data, by = "sample")
  res
}

# Get feature distributions -----------------------------------------------

get_features <- function(CN_data,
                         cores = 1,
                         genome_build = c("hg19", "hg38")) {
  genome_build <- match.arg(genome_build)
  # get chromosome lengths and centromere locations
  chrlen <- get_genome_annotation(data_type = "chr_size", genome_build = genome_build)
  centromeres <- get_genome_annotation(data_type = "centro_loc", genome_build = genome_build)

  oplan <- future::plan()
  future::plan("multiprocess", workers = cores)
  on.exit(future::plan(oplan), add = TRUE)


  features <- c(
    "segsize", "bp10MB", "osCN", "bpchrarm",
    "changepoint", "copynumber"
  )

  .get_feature <- function(i) {
    if (i == "segsize") {
      send_info("Getting segment size...")
      getSegsize(CN_data)
    } else if (i == "bp10MB") {
      send_info("Getting breakpoint count per 10 Mb...")
      getBPnum(CN_data, chrlen)
    } else if (i == "osCN") {
      send_info("Getting length of chains of oscillating copy number...")
      getOscilation(CN_data)
    } else if (i == "bpchrarm") {
      send_info("Getting breakpoint count per chromosome arm...")
      getCentromereDistCounts(CN_data, centromeres)
    } else if (i == "changepoint") {
      send_info("Getting change-point copy number change...")
      getChangepointCN(CN_data)
    } else {
      send_info("Getting copy number...")
      getCN(CN_data)
    }
  }

  res <- furrr::future_map(features, .get_feature,
    .progress = TRUE
  )
  res <- res %>% setNames(features)
  res
}



get_features_wang <- function(CN_data,
                              cores = 1,
                              genome_build = c("hg19", "hg38"),
                              feature_setting = sigminer::CN.features) {
  genome_build <- match.arg(genome_build)
  # get chromosome lengths and centromere locations
  chrlen <- get_genome_annotation(data_type = "chr_size", genome_build = genome_build)
  centromeres <- get_genome_annotation(data_type = "centro_loc", genome_build = genome_build)

  oplan <- future::plan()
  future::plan("multiprocess", workers = cores)
  on.exit(future::plan(oplan), add = TRUE)

  features <- unique(feature_setting$feature)

  .get_feature <- function(i) {
    if (i == "SS") {
      send_info("Getting (log10 based) segment size...")
      zz <- getSegsize(CN_data)
      zz$value <- log10(zz$value)
      zz
    } else if (i == "BP10MB") {
      send_info("Getting breakpoint count per 10 Mb...")
      getBPnum(CN_data, chrlen)
    } else if (i == "OsCN") {
      send_info("Getting length of chains of oscillating copy number...")
      getOscilation(CN_data)
    } else if (i == "BPArm") {
      send_info("Getting breakpoint count per chromosome arm...")
      getCentromereDistCounts(CN_data, centromeres)
    } else if (i == "CNCP") {
      send_info("Getting change-point copy number change...")
      getChangepointCN(CN_data)
    } else if (i == "CN") {
      send_info("Getting copy number...")
      getCN(CN_data)
    } else if (i == "BoChr") {
      send_info("Getting burden of chromosome...")
      getBoChr(CN_data, genome_build)
    } else if (i == "NChrV") {
      send_info("Getting number of autosome with CNV...")
      getNChrV(CN_data, genome_build)
    } else if (i == "NC50") {
      send_info("Getting the minimal number of chromosome with 50% CNV...")
      getNC50(CN_data, genome_build)
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
    dat <- CN_feature$value
    send_info("Fitting feature: ", feature_name, ".")
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
    send_info(
      "More details about reference components please see {.url https://github.com/ShixiangWang/absoluteCNVdata}"
    )
    if (!file.exists("Nat_Gen_component_parameters.rds")) {
      send_info(
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
        times = sapply(CN_data, function(x) {
          nrow(x)
        })
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
    send_stop(
      "If input is a data.frame, must have 4 necessary columns (chr, start, end, sample) and 1 optional column (segVal)."
    )
  }

  data.table::setDT(segTab)
  segTab$start <- as.numeric(segTab$start)
  segTab$end <- as.numeric(segTab$end)
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
    segTab <- segTab[valid_chr, on = "chromosome"]
    send_success("Some invalid segments (not as 1:22 and X, Y) dropped.")
  }

  arm_data <- get_ArmLocation(genome_build)
  data.table::setDT(arm_data)

  segTab <- data.table::merge.data.table(segTab, arm_data, by.x = "chromosome", by.y = "chrom", all.x = TRUE)

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

  cbind(
    data.table::as.data.table(segTab)[, colnames(arm_data)[-1] := NULL],
    data.table::as.data.table(annot)
  )
}


# Get arm location --------------------------------------------------------

get_ArmLocation <- function(genome_build = c("hg19", "hg38")) {
  genome_build <- match.arg(genome_build)
  # get chromosome lengths and centromere locations
  chrlen <- get_genome_annotation(data_type = "chr_size", genome_build = genome_build)
  centromeres <- get_genome_annotation(data_type = "centro_loc", genome_build = genome_build)

  # compute and get results
  res <- data.frame(
    chrom = vector(mode = "character", length = 24),
    p_start = vector("numeric", length = 24),
    p_end = vector("numeric", length = 24),
    p_length = vector("numeric", length = 24),
    q_start = vector("numeric", length = 24),
    q_end = vector("numeric", length = 24),
    q_length = vector("numeric", length = 24),
    total_size = vector("numeric", length = 24),
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
                                 genome_measure = c("called", "wg")) {
  genome_build <- match.arg(genome_build)
  genome_measure <- match.arg(genome_measure)

  if (ncol(segTab) > 5) {
    segTab <- segTab[, 1:5]
  }
  if (ncol(segTab) == 5) {
    colnames(segTab) <- c("chromosome", "start", "end", "segVal", "sample")
  } else {
    stop(
      "Input must have 5 ordered columns (chr, start, end, segVal, sample)."
    )
  }

  # Handle sex
  sex <- getOption("sigminer.sex", default = "female")
  if (is.character(sex)) {
    if (sex == "male") {
      segTab <- segTab %>%
        dplyr::as_tibble() %>%
        dplyr::mutate(
          segVal = ifelse(.data$chromosome %in% c("chrX", "chrY"), 2L * .data$segVal, .data$segVal)
        ) %>%
        data.table::as.data.table()
    }
  }
  if (is.data.frame(sex)) {
    segTab <- segTab %>%
      dplyr::as_tibble() %>%
      dplyr::left_join(sex, by = "sample") %>%
      dplyr::mutate(
        segVal = dplyr::case_when(
          .data$sex == "male" & .data$chromosome %in% c("chrX", "chrY") ~ 2L * .data$segVal,
          TRUE ~ .data$segVal
        )
      ) %>%
      dplyr::select(-"sex") %>%
      dplyr::select(c("chromosome", "start", "end", "segVal", "sample")) %>%
      data.table::as.data.table()
  }

  data.table::setDT(segTab)
  segTab$start <- as.numeric(segTab$start)
  segTab$end <- as.numeric(segTab$end)

  autosome <- paste0("chr", 1:22)

  if (genome_measure == "wg") {
    chrlen <- get_genome_annotation(
      data_type = "chr_size",
      chrs = autosome,
      genome_build = genome_build
    )
    total_size <- sum(chrlen[["size"]])

    seg_sum1 <- segTab %>%
      dplyr::as_tibble() %>%
      dplyr::group_by(.data$sample) %>%
      dplyr::summarise(
        n_of_seg = sum(!is.na(segVal)),
        n_of_cnv = sum(segVal != 2),
        n_of_amp = sum(segVal > 2),
        n_of_del = sum(segVal < 2)
      )

    seg_sum2 <- segTab %>%
      dplyr::as_tibble() %>%
      dplyr::filter(chromosome %in% autosome) %>%
      dplyr::group_by(.data$sample) %>%
      dplyr::summarise(
        n_of_vchr = length(unique(chromosome[segVal != 2])),
        cna_burden = sum(end[segVal != 2] - start[segVal != 2] + 1) / total_size
      )

    seg_summary <- dplyr::full_join(seg_sum1, seg_sum2, by = "sample") %>%
      data.table::as.data.table()
  } else {
    seg_sum1 <- segTab %>%
      dplyr::as_tibble() %>%
      dplyr::group_by(.data$sample) %>%
      dplyr::summarise(
        n_of_seg = sum(!is.na(segVal)),
        n_of_cnv = sum(segVal != 2),
        n_of_amp = sum(segVal > 2),
        n_of_del = sum(segVal < 2)
      )

    seg_sum2 <- segTab %>%
      dplyr::as_tibble() %>%
      dplyr::filter(chromosome %in% autosome) %>%
      dplyr::group_by(.data$sample) %>%
      dplyr::summarise(
        n_of_vchr = length(unique(chromosome[segVal != 2])),
        cna_burden = sum(end[segVal != 2] - start[segVal != 2] + 1) / sum(end - start + 1)
      )

    seg_summary <- dplyr::full_join(seg_sum1, seg_sum2, by = "sample") %>%
      data.table::as.data.table()
  }

  seg_summary <- seg_summary[order(seg_summary$n_of_cnv)]
  seg_summary$cna_burden <- round(seg_summary$cna_burden, 3)
  attr(seg_summary, "Note") <- "n_of_vchr is the number of autosome with CNV.\n cna_burden also calculates from autosome."
  seg_summary
}

# Classify copy number features into components ------------------------------------
# Automatically handle and support official/custom copy number feature classfication setting

get_feature_components <- function(x) {
  stopifnot(is.data.frame(x), all(c("feature", "min", "max") %in% colnames(x)))

  x <- x %>%
    dplyr::mutate(
      component = dplyr::case_when(
        .data$min == .data$max ~ sprintf("%s[%s]", .data$feature, .data$min),
        is.infinite(.data$min) ~ sprintf("%s[<=%s]", .data$feature, .data$max),
        is.infinite(.data$max) ~ sprintf("%s[>%s]", .data$feature, .data$min),
        .data$max > .data$min ~ sprintf("%s[>%s & <=%s]", .data$feature, .data$min, .data$max),
        TRUE ~ NA_character_
      ),
      label = ifelse(.data$min == .data$max, "point", "range")
    ) %>%
    data.table::as.data.table()

  x <- x[, c("feature", "component", "label", "min", "max")]

  if (any(is.na(x$component))) {
    stop("Result componet column contain NAs, please check your input!")
  } else {
    class(x) <- c(class(x), "sigminer.features")
    return(x)
  }
}

# Global variables --------------------------------------------------------

utils::globalVariables(
  c(
    "i",
    "chrom",
    "chromosome",
    "segVal"
  )
)
