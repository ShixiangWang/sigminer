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

  # features <- unique(feature_setting$feature)
  features <- c("BP10MB", "CN", "SS", "CS", "AB")

  send_info("NOTE: this method derives features for each segment. Be patient...")

  .get_feature <- function(i) {
    if (i == "SS") {
      send_info("Getting floor(log10 based segment size) of eash segment...")
      zz <- getSegsize_v2(CN_data)
      # zz$value <- floor(log10(zz$value))  # do it in get_components
      zz
    } else if (i == "CN") {
      send_info("Getting absolute copy number value of each segment...")
      getCN_v2(CN_data)
    } else if (i == "BP10MB") {
      send_info("Getting breakpoint number within upstream/downstream 5Mb-flank region of each segment...")
      getBP10MB_v2(CN_data)
    } else if (i == "CS") {
      send_info("Getting context change shape based on left and right sides of each segment...")
      getCS(CN_data)
    } else if (i == "AB") {
      send_info("Getting change extent on left and right sides of each segment...")
      getAB(CN_data)
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
    x$segsize <- x$end - x$start + 1L
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
  y <- purrr::map_df(abs_profiles, function(x) {
    x_cp <- data.table::copy(x)
    x_cp$region_start <- x$start - 5000000L
    x_cp$region_end <- x$start + 5000000L
    x_cp$region_start <- ifelse(x_cp$region_start <= 0, 1L, x_cp$region_start)
    x_cp[, c("start", "end", "segVal") := NULL]

    data.table::setkey(x_cp, "chromosome", "region_start", "region_end")
    x_overlap <- data.table::foverlaps(x, x_cp,
      by.x = c("chromosome", "start", "end"),
      type = "any"
    )
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

## Get context shape (left and right size):
## HH, HL, LH, LL (H: high, L: low)
getCS <- function(abs_profiles) {
  y <- purrr::map_df(abs_profiles, function(x) {
    x %>%
      dplyr::as_tibble() %>%
      dplyr::group_by(.data$chromosome) %>%
      dplyr::mutate(
        lv = diff(c(2L, .data$segVal)),
        rv = -diff(c(.data$segVal, 2L))
      ) %>%
      dplyr::mutate(
        value = dplyr::case_when(
          .data$lv <= 0 & .data$rv <= 0 ~ "HH",
          .data$lv < 0 & .data$rv > 0 ~ "HL",
          .data$lv > 0 & .data$rv < 0 ~ "LH",
          .data$lv >= 0 & .data$rv >= 0 ~ "LL"
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(c("sample", "value", "Index"))
  }) %>%
    data.table::as.data.table()

  y[order(y$Index)]
}

## Get left and right CN change cut off
## (left&right):AA, AB, BA, BB (A <= 2, B >2)
# abs_profiles <- sigminer::get_cnlist(cn, add_index = TRUE)
getAB <- function(abs_profiles) {
  y <- purrr::map_df(abs_profiles, function(x) {
    x %>%
      dplyr::as_tibble() %>%
      dplyr::group_by(.data$chromosome) %>%
      dplyr::mutate(
        lv = abs(diff(c(2L, .data$segVal))),
        rv = abs(-diff(c(.data$segVal, 2L)))
      ) %>%
      dplyr::mutate(
        value = dplyr::case_when(
          .data$lv <= 2 & .data$rv <= 2 ~ "AA",
          .data$lv <= 2 & .data$rv > 2 ~ "AB",
          .data$lv > 2 & .data$rv <= 2 ~ "BA",
          .data$lv > 2 & .data$rv > 2 ~ "BB"
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(c("sample", "value", "Index"))
  }) %>%
    data.table::as.data.table()

  y[order(y$Index)]
}
# Get components ----------------------------------------------------------

## Use two classification systems:
## Standard system (S): keep simpler
## Complex system (C): keep comprehensive
get_components_mutex <- function(CN_features) {
  feature_names <- names(CN_features)

  purrr::map2(CN_features, feature_names, call_component)
}

call_component <- function(f_dt, f_name) {
  f_dt <- data.table::copy(f_dt)

  if (f_name == "BP10MB") {
    f_dt$C_BP10MB <- cut(f_dt$value,
      breaks = c(-Inf, 2L, 4L, 6L, 8L, 10L, Inf),
      labels = c("2-", "4-", "6-", "8-", "10-", "11+")
    )
  } else if (f_name == "CN") {
    f_dt$S_CN <- cut(f_dt$value,
      breaks = c(-Inf, 0:4, Inf),
      labels = c(as.character(0:4), "5+")
    )
    f_dt$C_CN <- cut(f_dt$value,
      breaks = c(-Inf, 0:8, Inf),
      labels = c(as.character(0:8), "9+")
    )
  } else if (f_name == "SS") {
    f_dt$S_SS <- cut(f_dt$value,
      breaks = c(-Inf, 50000L, 500000L, 5000000L, Inf),
      labels = c("S", "M", "L", "E")
    )
    f_dt$C_SS <- cut(floor(log10(f_dt$value)),
      breaks = c(-Inf, 2:7, Inf),
      labels = c("2-", as.character(3:7), "8+")
    )
  } else if (f_name == "CS") {
    f_dt$S_CS <- f_dt$C_CS <- factor(f_dt$value, levels = c("HH", "HL", "LH", "LL"))
  } else if (f_name == "AB") {
    f_dt$S_AB <- f_dt$C_AB <- factor(f_dt$value, levels = c("AA", "AB", "BA", "BB"))
  }
  f_dt$value <- NULL
  f_dt
}


# Get matrix --------------------------------------------------------------
# CN_components <- cn_tally$components
get_matrix_mutex <- function(CN_components, indices = NULL) {
  merged_dt <- purrr::reduce(CN_components, merge, by = c("sample", "Index"), all = TRUE)

  if (!is.null(indices)) {
    merged_dt <- merged_dt[merged_dt$Index %in% indices]
  }

  ## Standard Classificiations
  ## Complex Classifications
  dt_s <- merged_dt[, colnames(merged_dt) == "sample" | startsWith(colnames(merged_dt), "S_"), with = FALSE]
  dt_c <- merged_dt[, colnames(merged_dt) == "sample" | startsWith(colnames(merged_dt), "C_"), with = FALSE]

  ## 1. handle standard way
  s_class_levels <- vector_to_combination(levels(dt_s$S_SS), levels(dt_s$S_CS),
    levels(dt_s$S_CN), levels(dt_s$S_AB),
    c_string = ":"
  )


  dt_s$s_class <- paste(dt_s$S_SS, dt_s$S_CS, dt_s$S_CN, dt_s$S_AB, sep = ":")
  dt_s$s_class <- factor(dt_s$s_class, levels = s_class_levels)
  s_mat <- classDT2Matrix(dt_s, samp_col = "sample", component_col = "s_class") %>%
    as.data.frame()

  ## Code to combine catagories with very few counts
  ## for AB|BA|BB
  TM_set <- c(
    "S:HH:0:",
    "S:HL:1:", "S:HL:2:", "S:HL:3:", "S:HL:4:",
    "S:LH:1:", "S:LH:2:", "S:LH:3:", "S:LH:4:",
    "S:LL:3:", "S:LL:4:",
    "M:HH:0:", "M:HH:1:",
    "M:HL:1:", "M:HL:2:", "M:HL:3:", "M:HL:4:",
    "M:LH:1:", "M:LH:2:", "M:LH:3:", "M:LH:4:",
    "M:LL:3:", "M:LL:4:",
    "L:HH:0:", "L:HH:1:", "L:HH:2:", "L:HH:3:", "L:HH:4:", "L:HH:5+:",
    "L:HL:1:", "L:HL:2:", "L:HL:3:", "L:HL:4:",
    "L:LH:1:", "L:LH:2:", "L:LH:3:", "L:LH:4:",
    "L:LL:3:", "L:LL:4:",
    "E:HH:0:", "E:HH:1:", "E:HH:2:", "E:HH:3:", "E:HH:4:", "E:HH:5+:",
    "E:HL:1:", "E:HL:2:", "E:HL:3:", "E:HL:4:",
    "E:LH:1:", "E:LH:2:", "E:LH:3:", "E:LH:4:",
    "E:LL:3:", "E:LL:4:"
  )

  for (i in TM_set) {
    s_mat[[paste0(i, "3C")]] <- s_mat[[paste0(i, "AB")]] + s_mat[[paste0(i, "BA")]] + s_mat[[paste0(i, "BB")]]
    s_mat[[paste0(i, "AB")]] <- s_mat[[paste0(i, "BA")]] <- s_mat[[paste0(i, "BB")]] <- NULL
  }

  ## for AB|BA
  DM_set <- c(
    "S:HH:1:", "S:HH:3:", "S:HH:4:", "S:HH:5+:",
    "S:LL:5+:",
    "M:HH:2:", "M:HH:3:", "M:HH:4:", "M:HH:5+:",
    "M:LL:5+:"
  )

  for (i in DM_set) {
    s_mat[[paste0(i, "2C")]] <- s_mat[[paste0(i, "AB")]] + s_mat[[paste0(i, "BA")]]
    s_mat[[paste0(i, "AB")]] <- s_mat[[paste0(i, "BA")]] <- NULL
  }

  ## Delete 0 count classifications
  ## some classes have already been deleted in the previous step
  s_mat[
    ,
    c(
      "E:HL:0:AA", "E:HL:0:AB", "E:HL:0:BA", "E:HL:0:BB", "E:HL:1:AB",
      "E:HL:1:BB", "E:HL:2:AB", "E:HL:2:BB", "E:LH:0:AA", "E:LH:0:AB",
      "E:LH:0:BA", "E:LH:0:BB", "E:LH:1:BA", "E:LH:1:BB", "E:LH:2:BA",
      "E:LH:2:BB", "E:LL:0:AA", "E:LL:0:AB", "E:LL:0:BA", "E:LL:0:BB",
      "E:LL:1:AB", "E:LL:1:BA", "E:LL:1:BB", "E:LL:2:AB", "E:LL:2:BA",
      "E:LL:2:BB", "L:HL:0:AA", "L:HL:0:AB", "L:HL:0:BA", "L:HL:0:BB",
      "L:HL:1:AB", "L:HL:1:BB", "L:HL:2:AB", "L:HL:2:BB", "L:LH:0:AA",
      "L:LH:0:AB", "L:LH:0:BA", "L:LH:0:BB", "L:LH:1:BA", "L:LH:1:BB",
      "L:LH:2:BA", "L:LH:2:BB", "L:LL:0:AA", "L:LL:0:AB", "L:LL:0:BA",
      "L:LL:0:BB", "L:LL:1:AB", "L:LL:1:BA", "L:LL:1:BB", "L:LL:2:AB",
      "L:LL:2:BA", "L:LL:2:BB", "M:HL:0:AA", "M:HL:0:AB", "M:HL:0:BA",
      "M:HL:0:BB", "M:HL:1:AB", "M:HL:1:BB", "M:HL:2:AB", "M:HL:2:BB",
      "M:LH:0:AA", "M:LH:0:AB", "M:LH:0:BA", "M:LH:0:BB", "M:LH:1:BA",
      "M:LH:1:BB", "M:LH:2:BA", "M:LH:2:BB", "M:LL:0:AA", "M:LL:0:AB",
      "M:LL:0:BA", "M:LL:0:BB", "M:LL:1:AB", "M:LL:1:BA", "M:LL:1:BB",
      "M:LL:2:AB", "M:LL:2:BA", "M:LL:2:BB", "S:HL:0:AA", "S:HL:0:AB",
      "S:HL:0:BA", "S:HL:0:BB", "S:HL:1:AB", "S:HL:1:BB", "S:HL:2:AB",
      "S:HL:2:BB", "S:HL:3:BB", "S:LH:0:AA", "S:LH:0:AB", "S:LH:0:BA",
      "S:LH:0:BB", "S:LH:1:BA", "S:LH:1:BB", "S:LH:2:BA", "S:LH:2:BB",
      "S:LL:0:AA", "S:LL:0:AB", "S:LL:0:BA", "S:LL:0:BB", "S:LL:1:AB",
      "S:LL:1:BA", "S:LL:1:BB", "S:LL:2:AB", "S:LL:2:BA", "S:LL:2:BB"
    )
  ] <- NULL

  s_mat <- as.matrix(s_mat[, sort(colnames(s_mat))])

  ## 2. hanlde complex way
  c_class_levels <- vector_to_combination(levels(dt_c$C_SS), levels(dt_c$C_CS), levels(dt_c$C_CN),
    paste0("P", levels(dt_c$C_BP10MB)),
    c_string = ":"
  )
  dt_c$c_class <- paste(dt_c$C_SS, dt_c$C_CS, dt_c$C_CN,
    paste0("P", dt_c$C_BP10MB),
    sep = ":"
  )
  dt_c$c_class <- factor(dt_c$c_class, levels = c_class_levels)
  c_mat <- classDT2Matrix(dt_c, samp_col = "sample", component_col = "c_class")

  return(list(s_mat = s_mat, c_mat = c_mat))
}


classDT2Matrix <- function(dt, samp_col, component_col) {
  dt.summary <- dt[, .N, by = c(samp_col, component_col)]
  mat <- as.data.frame(data.table::dcast(dt.summary,
    formula = as.formula(paste(samp_col, "~", component_col)),
    fill = 0, value.var = "N", drop = FALSE
  ))

  rownames(mat) <- mat[, 1]
  mat <- mat[, -1]
  mat %>% as.matrix()
}

utils::globalVariables(
  c(".N")
)
