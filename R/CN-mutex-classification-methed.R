## This is a new copy number segment classification method
## It focuses on **segment** and segment's context, in both concept and implementation,
## it is same as SBS, DBS, INDEL etc.
##
## Firstly, we collect copy number (variation) features from segmental copy number data
## similar to previous work but here we focus on each **segment**.
## Secondly, we classified all segments into mutually exclusive types based on features.
get_features_mutex <- function(CN_data,
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

  #features <- unique(feature_setting$feature)
  features <- c("CN", "SS")
  # c("BP10MB", "CN", "SS", "CNCP-L", "CNCP-R", "CNCP-M", "OsCN") # more?

  .get_feature <- function(i) {
    if (i == "SS") {
      send_info("Getting floor(log10 based segment size)...")
      zz <- getSegsize_v2(CN_data)
      zz$value <- floor(log10(zz$value))
      ## TODO: set <= 1 to 2
      zz
    } else if (i == "CN") {
      send_info("Getting copy number...")
      getCN_v2(CN_data)
    }
  }

  res <- furrr::future_map(features, .get_feature,
                           .progress = TRUE
  )
  res <- res %>% setNames(features)
  res
}


## Make sure the index is also returned
getSegsize_v2 <- function(abs_profiles) {
  segsize <- purrr::map_df(abs_profiles, function(x) {
    x$segsize <- x$end - x$start + 1
    x[, c("sample", "segsize", "Index"), with = FALSE]
  })
  colnames(segsize) <- c("sample", "value", "Index")
  segsize
}

getCN_v2 <- function(abs_profiles) {
  abs_profiles <- handle_sex(abs_profiles)

  cn <- purrr::map_df(abs_profiles, function(x) {
    x[, c("sample", "segVal", "Index"), with = FALSE]
  })
  colnames(cn) <- c("sample", "value", "Index")
  cn
}


getBP10MB_v2 <- function(abs_profiles) {
  y = purrr::map_df(abs_profiles, function(x){

    x_cp = data.table::copy(x)
    x_cp$region_start <- x$start - 5000000L
    x_cp$region_end <- x$start + 5000000L
    x_cp$region_start <- ifelse(x_cp$region_start <= 0, 1L, x_cp$region_start)
    x_cp[, c("start", "end", "segVal") := NULL]

    data.table::setkey(x_cp, "chromosome", "region_start", "region_end")
    x_overlap <- data.table::foverlaps(x, x_cp,
                                       by.x = c("chromosome", "start", "end"),
                                       type = "any")
    x %>%
      dplyr::as_tibble() %>%
      dplyr::full_join(
        x_overlap %>%
          dplyr::as_tibble() %>%
          dplyr::group_by(.data$Index) %>%
          dplyr::summarise(
            value = dplyr::n() - 1L
          ),
        by = c("Index")
      ) %>%
      dplyr::select(c("sample", "value", "Index"))
  }) %>%
    data.table::as.data.table()

  y[order(y$Index)]

}

getCNCP_Left <- function(abs_profiles) {
  y = purrr::map_df(abs_profiles, function(x){
    x %>%
      dplyr::as_tibble() %>%
      dplyr::group_by(.data$chromosome) %>%
      dplyr::mutate(value = abs(c(0L, diff(.data$segVal)))) %>%
      dplyr::ungroup() %>%
      dplyr::select(c("sample", "value", "Index"))
  }) %>%
    data.table::as.data.table()

  y[order(y$Index)]
}

getCNCP_Right <- function(abs_profiles) {
  y = purrr::map_df(abs_profiles, function(x){
    x %>%
      dplyr::as_tibble() %>%
      dplyr::group_by(.data$chromosome) %>%
      dplyr::mutate(value = rev(abs(c(0L, diff(rev(.data$segVal)))))) %>%
      dplyr::ungroup() %>%
      dplyr::select(c("sample", "value", "Index"))
  }) %>%
    data.table::as.data.table()

  y[order(y$Index)]
}

getCNCP_Max <- function(abs_profiles){
  df_left <- getCNCP_Left(abs_profiles) %>%
    dplyr::as_tibble() %>%
    dplyr::rename(left_value = .data$value)
  df_right <- getCNCP_Right(abs_profiles) %>%
    dplyr::as_tibble() %>%
    dplyr::rename(right_value = .data$value)

  y <- dplyr::full_join(df_left, df_right %>% dplyr::select(-"sample"), by = "Index") %>%
    dplyr::mutate(value = pmax(.data$left_value, .data$right_value)) %>%
    dplyr::select(c("sample", "value", "Index")) %>%
    data.table::as.data.table()

  y[order(y$Index)]
}
