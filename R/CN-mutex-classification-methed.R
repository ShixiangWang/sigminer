## This is a new copy number segment classification method
## It focuses on **segment** and segment's context, in both concept and implementation,
## it is same as SBS, DBS, INDEL etc.
##
## Firstly, we collect copy number (variation) features from segmental copy number data
## similar to previous work but here we focus on each **segment**.
## Secondly, we classified all segments into mutually exclusive types based on features.
get_features_mutex <- function(CN_data,
                               cores = 1) {

  oplan <- future::plan()
  future::plan("multiprocess", workers = cores)
  on.exit(future::plan(oplan), add = TRUE)

  #features <- unique(feature_setting$feature)
  features <- c("BP10MB", "CN", "SS", "CNCP_M", "OsCN", "StepRising", "StepFalling") # more?

  send_info("NOTE: this method derives features for each segment.")

  .get_feature <- function(i) {
    if (i == "SS") {
      send_info("Getting floor(log10 based segment size) of eash segment...")
      zz <- getSegsize_v2(CN_data)
      zz$value <- floor(log10(zz$value))
      ## TODO: set <= 1 to 2
      zz
    } else if (i == "CN") {
      send_info("Getting absolute copy number value of each segment...")
      getCN_v2(CN_data)
    } else if (i == "BP10MB") {
      send_info("Getting breakpoint number within upstream/downstream 5Mb-flank region of each segment...")
      getBP10MB_v2(CN_data)
    } else if (i == "OsCN") {
      send_info("Getting maximum length of chains of oscillating copy number for each segment...")
      getOscilation(CN_data, use_index = TRUE)
    } else if (i == "CNCP_M") {
      send_info("Getting maximum change-point amplitude at both side of each segment...")
      getCNCP_Max(CN_data)
    } else if (i == "StepRising") {
      send_info("Getting rising step size for each segment...")
      getSteps(CN_data, rising = TRUE)
    } else if (i == "StepFalling") {
      send_info("Getting falling step size for each segment...")
      getSteps(CN_data, rising = FALSE)
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
  segsize[order(segsize$Index)]
}

getCN_v2 <- function(abs_profiles) {
  abs_profiles <- handle_sex(abs_profiles)

  cn <- purrr::map_df(abs_profiles, function(x) {
    x[, c("sample", "segVal", "Index"), with = FALSE]
  })
  colnames(cn) <- c("sample", "value", "Index")
  cn[order(cn$Index)]
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

## Examples:
## Rising 2-3-4-5 (size:4)
## Falling 5-4-3-2 (size:4)
## 2-3 (rising size:2)
## 2 (size:0)
getSteps <- function(abs_profiles, rising = TRUE) {
  ## rising: up
  ## falling: down
  res <- purrr::map_df(abs_profiles, function(x) {
    x <- x %>%
      dplyr::as_tibble() %>%
      dplyr::group_by(.data$chromosome) %>%
      tidyr::nest() %>%
      dplyr::mutate(value = purrr::map(
        data, function(df) {
          currseg <- df$segVal
          len_seg <- length(currseg)
          if (len_seg >= 2) {
            if (rising) {
              target <- which(diff(currseg) > 0) + 1L
            } else {
              target <- which(diff(currseg) < 0) + 1L
            }
            steps <- vector("integer", length = len_seg)
            modify_list <- split(
              target,
              findInterval(target, target[which(diff(target) > 1)] + 1L)  # same if +2L
            )
            for (i in modify_list) {
              steps[c(i[1] - 1L, i)] <- length(i) + 1L
            }
            return(steps)
          } else {
            return(0L)
          }
        }
      )) %>%
      dplyr::ungroup()

    x %>% tidyr::unnest(c("data", "value")) %>%
      dplyr::select(c("value", "Index")) %>%
      data.table::as.data.table()
  }, .id = "sample")

  res[order(res$Index)]
}



# Get components ----------------------------------------------------------

## Use two classification systems:
## Standard system (S): keep simpler
## Complex system (S): keep comprehensive
get_components_mutex <- function(CN_features) {
  feature_names <- names(CN_features)

  purrr::map2(CN_features, feature_names, call_component)
}

call_component <- function(f_dt, f_name) {
  f_dt <- data.table::copy(f_dt)

  if (f_name == "BP10MB") {
    f_dt$C_BP10MB = cut(f_dt$value,
                        breaks = c(-Inf, 2L, 4L, 6L, 8L, 10L, Inf),
                        labels = c("2-", "4-", "6-", "8-", "10-", "11+"))
  } else if (f_name == "CN") {
    f_dt$S_CN = cut(f_dt$value,
                        breaks = c(-Inf, 0:8, Inf),
                        labels = c(as.character(0:8), "9+"))
    f_dt$L_CN = f_dt$S_CN
  } else if (f_name == "SS") {
    f_dt$S_CN = cut(f_dt$value,
                    breaks = c(-Inf, 50000, 5000000, Inf),
                    labels = c("S", "M", "L"))
    f_dt$L_CN = cut(f_dt$value,
                    breaks = c(-Inf, 2:7, Inf),
                    labels = c("2-", as.character(3:7), "8+"))
  } else if (f_name == "CNCP_M") {

  } else if (f_name == "OsCN") {
    f_dt$S_OsCN = cut(f_dt$value,
                    breaks = c(-Inf, 2L, Inf),
                    labels = c("N", "O"))
    f_dt$L_OsCN = cut(f_dt$value,
                    breaks = c(-Inf, 2L, 5L, Inf),
                    labels = c("N", "SO", "LO"))

  } else if (f_name == "StepRising") {
    f_dt$S_StepR = cut(f_dt$value,
                      breaks = c(-Inf, 2L, Inf),
                      labels = c("N", "R"))
    f_dt$L_StepR = cut(f_dt$value,
                      breaks = c(-Inf, 2L, 3L, Inf),
                      labels = c("N", "SR", "LR"))

  } else if (f_name == "StepFalling") {
    f_dt$S_StepF = cut(f_dt$value,
                       breaks = c(-Inf, 2L, Inf),
                       labels = c("N", "F"))
    f_dt$L_StepF = cut(f_dt$value,
                       breaks = c(-Inf, 2L, 3L, Inf),
                       labels = c("N", "SF", "LF"))
  }
}
