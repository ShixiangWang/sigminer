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
