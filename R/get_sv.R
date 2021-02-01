# Handling SV classes -------------------------------------------------------

# read sv -----------------------------------------------------------------
read_sv <- function(input) {
  # read sv file and convert it to data.table
  if (is.data.frame(x = input)) {
    input <- data.table::as.data.table(input)
  }
  else {
    input <- data.table::fread(
      file = input,
      data.table = TRUE, header = TRUE
    )
  }

  # check necessary name
  necessary.fields <- c(
    "sample",
    "chr1", "start1", "end1",
    "chr2", "start2", "end2",
    "strand1", "strand2",
    "svclass"
  )

  colnames(input)  = tolower(colnames(input))
  idx = necessary.fields %in% colnames(input)
  if (!all(idx)) {
    stop(
      "Missing required fields from SV: ",
      paste(necessary.fields[!idx], collapse = ", ")
    )
  }

  # drop unnecessary fields
  input <- subset(input, select = necessary.fields)

  message("Succesfully read SV!")
  class(input) = c("SV", class(input))
  return(input)
}


# split by sample and collect in a list : svlist --------------------------

get_svlist <- function(data) {
  index <- seq(1, nrow(data))
  data$Index <- index
  res <- split(data, by = "sample")
}


# get size
getRearrSize_v1 <- function(sv_profiles) {
  # sv_profiles <- svlist
  rearrsize <- purrr::map_df(sv_profiles, function(x) {
    # x <- sv_profiles$PD26851a[1]
    if (x$svclass == "translocation") x$rearrsize <- NA
    length <- x$end2 - x$start1
    if (length < 1000) x$rearrsize <- "<1kb"
    if (length >= 1000 & length < 10000) x$rearrsize <- "1-10kb"
    if (length >= 10000 & length < 100000) x$rearrsize <- "10-100kb"
    if (length >= 100000 & length < 1000000) x$rearrsize <- "100kb-1Mb"
    if (length >= 1000000 & length <= 10000000) x$rearrsize <- "1Mb-10Mb"
    if (length > 10000000) x$rearrsize <- ">10Mb"
    x[, c("sample", "rearrsize", "Index"), with = FALSE]
    # x[, c("sample", "rearrsize")]
  })
  colnames(rearrsize) <- c("sample", "value", "Index")
  rearrsize <- rearrsize %>%
    data.table::as.data.table()
  rearrsize[order(rearrsize$Index)]
}




# copy from UCL-Research-Department-of-Pathology/RESIN --------------------
getDists <- function(chrom1, pos1, chrom2, pos2, doPCF = FALSE) {
  pos1 <- as.numeric(pos1)
  pos2 <- as.numeric(pos2)
  # get segments per chromosome
  dists <- sapply(unique(c(chrom1, chrom2)), FUN = function(x) {
    bp <- matrix(NA, ncol = 2, nrow = 0)
    # positions separate for each member of breakpoint
    pos <- c()
    index1 <- which(chrom1 == x)
    if (length(index1) > 0) {
      pos <- c(pos, pos1[index1])
      bp <- rbind(bp, cbind(index1, 1))
    }
    index2 <- which(chrom2 == x)
    if (length(index2) > 0) {
      pos <- c(pos, pos2[index2])
      bp <- rbind(bp, cbind(index2, 2))
    }
    # distances between breakpoint
    dists <- diff(sort(pos))
    # segment
    forCN <- data.frame(chrom = x, pos = sort(pos), dist = dists[order(pos)])
    if (!doPCF) {
      return(forCN[, 3])
    }
    return(list(info = forCN, seg = suppressMessages(copynumber::pcf(forCN, gamma = 25, kmin = 10))))
  }, simplify = FALSE)
  return(dists)
}



# get clustered -----------------------------------------------------------
getClustered_v1 <- function(sv_profiles, threshold = NULL) {
  # sv_profiles <- svlist
  # x <- sv_profiles$PD26879a
  # x <- sv_profiles$PD26851a
  # threshold = NULL
  # each list each row apply function
  clustered <- purrr::map(sv_profiles, function(x) {
    # get pos and chrom
    pos1 <- x$start1
    pos2 <- x$start2
    chrom1 <- x$chr1
    chrom2 <- x$chr2
    # get segments per chromosome
    dists <- getDists(chrom1, pos1, chrom2, pos2, doPCF = TRUE)
    if (is.null(threshold)) threshold <- 0.1 * mean(unlist(sapply(dists, FUN = function(x) x$info[, 3])), na.rm = TRUE)
    # which segments are below threshold
    regions <- do.call(rbind, sapply(dists, FUN = function(x) {
      x$seg[which(x$seg$mean < threshold), ]
    }, simplify = FALSE))
    if (nrow(regions) == 0) {
      clustered <- rep("non-clustered", length(chrom1))
      # return(rep("unclustered",length(chrom1)))
    } else {
      # which rearrangements are in clustered regions
      # print(paste0(regions$chrom,":",regions$start.pos,"-",regions$end.pos))
      regionGR <- as(paste0(regions$chrom, ":", regions$start.pos, "-", regions$end.pos), "GRanges")
      clustered <- rep("non-clustered", length = length(pos1))
      # print(paste0(chrom1,":",pos1,"-",pos1))
      index1 <- IRanges::findOverlaps(as(paste0(chrom1, ":", pos1, "-", pos1), "GRanges"), regionGR)@from
      if (length(index1) > 0) clustered[index1] <- "clustered"
      # print(paste0(chrom2,":",pos2,"-",pos2))
      index2 <- IRanges::findOverlaps(as(paste0(chrom2, ":", pos2, "-", pos2), "GRanges"), regionGR)@from
      if (length(index2) > 0) clustered[index2] <- "clustered"
    }
    x$clustered <- clustered
    x[, c("sample", "clustered", "Index"), with = FALSE]
  })

  clustered_dt <- plyr::ldply(clustered, data.frame, .id = NULL) %>%
    .[, c("sample", "clustered", "Index")] %>%
    data.table::as.data.table()
  colnames(clustered_dt) <- c("sample", "value", "Index")
  clustered_dt[order(clustered_dt$Index)]
}


# get type ----------------------------------------------------------------
getType_v1 <- function(sv_profiles) {
  # sv_profiles <- svlist
  type <- purrr::map_df(sv_profiles, function(x) {
    # x <- sv_profiles$PD26851a[1]
    if (x$svclass == "deletion") x$type <- "del"
    if (x$svclass == "inversion") x$type <- "inv"
    if (x$svclass == "tandem-duplication") x$type <- "tds"
    if (x$svclass == "translocation") x$type <- "trans"
    x[, c("sample", "type", "Index"), with = FALSE]
  })
  colnames(type) <- c("sample", "value", "Index")
  type <- type %>%
    data.table::as.data.table()
  type[order(type$Index)]
}




# get features list -------------------------------------------------------
# imitate : get_features_wang
get_features_sv <- function(sv_data) {
  field <- c(
    "clustered",
    "type",
    "size"
  )
  .get_feature <- function(i) {
    if (i == "clustered") {
      print("Getting clustered info...")
      zz <- getClustered_v1(sv_data)
      zz
    }
    else if (i == "type") {
      print("Getting type of segment ...")
      getType_v1(sv_data)
    }
    else if (i == "size") {
      print("Getting distance of two rearrange segments ...")
      getRearrSize_v1(sv_data)
    }
  }
  res <- furrr::future_map(field, .get_feature,
                           .progress = TRUE,
                           .options = furrr::furrr_options(seed = TRUE)
  )
  res <- res %>% setNames(field)
  res
}


# get cn component
# input : cn_feature, output : get_component_sv
get_components_sv <- function(CN_features) {
  # CN_features = test_feature
  feature_names <- names(CN_features)
  # 【pre】feature_names <- setdiff(names(CN_features), "LOH")

  purrr::map2(CN_features[feature_names], feature_names,
              .f = call_component
              # ,
              # extra = if ("LOH" %in% names(CN_features)) {
              #   CN_features$LOH
              # } else {
              #   NULL
              # },
              # XVersion = XVersion
  )
}

call_component <- function(f_dt, f_name) {
  # f_dt <- CN_features[feature_names]$clustered
  # f_name <- feature_names[1]
  f_dt <- data.table::copy(f_dt)
  if (f_name == "clustered") {
    f_dt$C_clustered <- factor(f_dt$value, levels = c("clustered", "non-clustered"))
  }
  if (f_name == "type") {
    f_dt$C_type <- factor(f_dt$value, levels = c("del", "inv", "tds", "trans"))
  }
  if (f_name == "size") {
    f_dt$C_size <- factor(f_dt$value, levels = c("<1kb", ">10Mb", "1-10kb", "10-100kb", "100kb-1Mb", "1Mb-10Mb"))
  }
  f_dt$value <- NULL
  f_dt
}




# get sv matrix
get_matrix_sv <- function(CN_components, indices = NULL) {
  # CN_components <- test_component

  merged_dt <- purrr::reduce(CN_components, merge, by = c("sample", "Index"), all = TRUE)
  dt_mg <- merged_dt
  sv_class_levels <- sigminer:::vector_to_combination(
    levels(dt_mg$C_clustered),
    levels(dt_mg$C_type),
    levels(dt_mg$C_size),
    c_string = ":"
  ) %>%
    ifelse(grepl("trans", .),
           sub("trans.+", "trans", .),
           .
    ) %>%
    unique()

  dt_mg$sv_class <- paste(dt_mg$C_clustered, dt_mg$C_type, dt_mg$C_size, sep = ":") %>%
    ifelse(grepl("trans", .),
           sub("trans.+", "trans", .),
           .
    )
  dt_mg$sv_class <- factor(dt_mg$sv_class, levels = sv_class_levels)

  sv_mat <- sigminer:::classDT2Matrix(dt_mg, samp_col = "sample", component_col = "sv_class") %>%
    as.data.frame()

  sv_mat <- as.matrix(sv_mat[, sort(colnames(sv_mat))])
  return(list(sv_mat = sv_mat))
}



# sig_tally_sv
sig_tally.SV <- function(object) {
  svlist <- get_svlist(object)
  print("successfully get svlist!")

  sv_features <- get_features_sv(svlist)
  print("successfully get sv features!")

  sv_component <- get_components_sv(sv_features)
  print("successfully get sv component!")

  sv_matrix <- get_matrix_sv(CN_components = sv_component)
  print("successfully get sv matrix!")

  res_list <- list(
    features = sv_features,
    components = sv_component,
    nmf_matrix = sv_matrix
  )
  return(res_list)
}
