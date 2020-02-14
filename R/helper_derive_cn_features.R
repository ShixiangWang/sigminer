# sigminer.sex can be a length-1 character indicating 'female' or 'male',
# also a data.frame with two columns 'sample' and 'sex'
#
# sigminer.copynumber.max can be an integer
options(sigminer.sex = "female", sigminer.copynumber.max = NA_integer_)

getSegsize <- function(abs_profiles) {
  segsize <- purrr::map_df(abs_profiles, function(x) {
    x$segsize <- x$end - x$start + 1
    x[, c("sample", "segsize"), with = FALSE]
  })
  colnames(segsize) <- c("ID", "value")
  segsize
}

getBPnum <- function(abs_profiles, chrlen) {
  res <- purrr::map_df(abs_profiles, function(x, chrlen) {
    calcBPnum <- function(df, c, chrlen) {
      intervals <-
        seq(1, chrlen[chrlen[, 1] == c, 2] + 10000000, 10000000)
      y <- tryCatch(
        hist(df$end[-nrow(df)],
          breaks = intervals,
          plot = FALSE
        )$counts,
        error = function(e) {
          stop(
            "Stop due to the following reason. Please check if your genome build is right.",
            "\n", e$message
          )
        }
      )
      y
    }

    x <- x %>%
      dplyr::as_tibble() %>%
      dplyr::group_by(.data$chromosome) %>%
      tidyr::nest() %>%
      dplyr::mutate(value = purrr::map2(
        data, .data$chromosome, calcBPnum,
        chrlen = chrlen
      )) %>%
      dplyr::ungroup() %>%
      dplyr::select(c("value"))
    data.table::data.table(
      value = purrr::reduce(x$value, c)
    )
  }, chrlen = chrlen, .id = "ID")

  return(res)
}

getOscilation <- function(abs_profiles) {
  oscCounts <- purrr::map_df(abs_profiles, function(x) {
    x <- x %>%
      dplyr::as_tibble() %>%
      dplyr::group_by(.data$chromosome) %>%
      tidyr::nest() %>%
      dplyr::mutate(value = purrr::map(
        data, function(df) {
          currseg <- df$segVal
          oscCounts <- c()
          if (length(currseg) > 3) {
            prevval <- currseg[1]
            count <- 0
            for (j in 3:length(currseg))
            {
              if (currseg[j] == prevval & currseg[j] != currseg[j - 1]) {
                count <- count + 1
              } else {
                oscCounts <- c(oscCounts, count)
                count <- 0
              }
              prevval <- currseg[j - 1]
            }
            return(oscCounts)
          } else {
            return(0L)
          }
        }
      )) %>%
      dplyr::ungroup() %>%
      dplyr::select(c("value"))
    data.table::data.table(
      value = purrr::reduce(x$value, c)
    )
  }, .id = "ID")

  oscCounts
}

getCentromereDistCounts <-
  function(abs_profiles, centromeres) {
    calcArmBP <- function(df, c, centromeres) {
      centstart <- centromeres[centromeres$chrom == c, 2]
      centend <- centromeres[centromeres$chrom == c, 3]

      # Remove segments located in centromere
      df <- df[!(df$start >= centstart & df$end <= centend), ]

      if (nrow(df) > 1) {
        starts <- df$start
        ends <- df$end

        p_count <- sum(ends < centend)
        q_count <- sum(ends > centend)

        if (!any(starts < centstart & ends > centend)) {
          # If there is no segments across p and q arms
          if (p_count > 0L) p_count <- p_count - 1L
        }
        if (q_count > 0L) q_count <- q_count - 1L

        return(c(p_count, q_count))
      } else {
        return(rep(0L, 2))
      }
    }

    res <- purrr::map_df(abs_profiles, function(x, centromeres) {
      x <- x %>%
        dplyr::as_tibble() %>%
        dplyr::group_by(.data$chromosome) %>%
        tidyr::nest() %>%
        dplyr::mutate(value = purrr::map2(
          data, .data$chromosome, calcArmBP,
          centromeres = centromeres
        )) %>%
        dplyr::ungroup() %>%
        dplyr::select(c("value"))
      data.table::data.table(
        value = purrr::reduce(x$value, c)
      )
    }, centromeres = centromeres, .id = "ID")

    return(res)
  }

getChangepointCN <- function(abs_profiles) {
  abs_profiles <- handle_sex(abs_profiles)

  cp <- purrr::map_df(abs_profiles, function(x) {
    x <- x %>%
      dplyr::as_tibble() %>%
      dplyr::group_by(.data$chromosome) %>%
      tidyr::nest() %>%
      dplyr::mutate(value = purrr::map(
        data, function(df) {
          if (nrow(df) <= 1) {
            return(0L)
          } else {
            return(as.integer(abs(diff(df$segVal))))
          }
        }
      )) %>%
      dplyr::ungroup() %>%
      dplyr::select(c("value"))
    data.table::data.table(
      value = purrr::reduce(x$value, c)
    )
  }, .id = "ID")

  return(cp)
}

getCN <- function(abs_profiles) {
  abs_profiles <- handle_sex(abs_profiles)

  cn <- purrr::map_df(abs_profiles, function(x) {
    x[, c("sample", "segVal"), with = FALSE]
  })
  colnames(cn) <- c("ID", "value")
  cn
}

# Number of Chromosome with CNV
getNChrV <- function(abs_profiles, genome_build = "hg38") {
  genome_build <- match.arg(genome_build, choices = c("hg19", "hg38"))
  abs_profiles <- handle_sex(abs_profiles)

  if (genome_build %in% c("hg19", "hg38")) {
    chrs <- c(paste0("chr", 1:22), "chrX")
  }

  cn <- purrr::map_df(abs_profiles, function(x) {
    x <- x[x$chromosome %in% chrs, c("sample", "chromosome", "segVal"), with = FALSE]
    x_cnv <- x[x$segVal != 2]
    value <- ifelse(nrow(x_cnv) == 0, 0L, length(unique(x_cnv$chromosome)))
    data.table::data.table(
      ID = x$sample[1],
      value = value
    )
  })
  cn
}

# The chromosome sequences (using integer as index) with copy number variation
# The count of this result represents the burden (contribution) of chromosome
getBoChr <- function(abs_profiles, genome_build = "hg38") {
  genome_build <- match.arg(genome_build, choices = c("hg19", "hg38"))
  abs_profiles <- handle_sex(abs_profiles)

  if (genome_build %in% c("hg19", "hg38")) {
    chrs <- c(paste0("chr", 1:22), "chrX")
    # Create a dict for mapping
    chrs_map <- as.character(1:23)
    names(chrs_map) <- chrs
  }

  cn <- purrr::map_df(abs_profiles, function(x) {
    x <- x[x$chromosome %in% chrs, c("sample", "chromosome", "segVal"), with = FALSE]
    x_cnv <- x[x$segVal != 2]
    if (nrow(x_cnv) == 0) {
      value <- 0L
    } else {
      value <- as.integer(chrs_map[x_cnv$chromosome])
    }
    data.table::data.table(
      ID = x$sample[1],
      value = value
    )
  })

  cn
}

# The minimal number of chromosome with 50% CNV
getNC50 <- function(abs_profiles, genome_build = "hg38") {
  genome_build <- match.arg(genome_build, choices = c("hg19", "hg38"))
  abs_profiles <- handle_sex(abs_profiles)

  if (genome_build %in% c("hg19", "hg38")) {
    chrs <- c(paste0("chr", 1:22), "chrX")
  }

  cn <- purrr::map_df(abs_profiles, function(x) {
    x <- x[x$chromosome %in% chrs, c("sample", "chromosome", "segVal"), with = FALSE]
    x_cnv <- x[x$segVal != 2]
    if (nrow(x_cnv) == 0) {
      value <- 0
    } else {
      value <- x_cnv %>%
        dplyr::as_tibble() %>%
        dplyr::count(.data$chromosome, sort = TRUE) %>%
        dplyr::mutate(
          n = cumsum(.data$n) / sum(.data$n),
          index = dplyr::row_number()
        ) %>%
        dplyr::filter(.data$n >= 0.5) %>%
        head(1) %>%
        dplyr::pull(.data$index)
    }
    data.table::data.table(
      ID = x$sample[1],
      value = value
    )
  })
  cn
}

count_components <- function(min, max, label, feature) {
  # Count all samples
  if (label == "point") {
    feature %>%
      dplyr::group_by(.data$ID) %>%
      dplyr::summarise(count = sum(.data$value == min, na.rm = TRUE)) %>%
      dplyr::pull("count")
  } else if (label == "range") {
    feature %>%
      dplyr::group_by(.data$ID) %>%
      dplyr::summarise(count = sum(.data$value > min & .data$value <= max, na.rm = TRUE)) %>%
      dplyr::pull("count")
  } else {
    stop("Bad labels for feature setting, can only be 'point' and 'range'!")
  }
}

count_components_wrapper <- function(feature_df, f_name, feature_setting) {
  samps <- unique(feature_df$ID)
  feature_df <- feature_df %>% dplyr::as_tibble()
  # Make sure sample order is consistent
  feature_df$ID <- factor(feature_df$ID, levels = samps)

  specific_f <- feature_setting[feature_setting$feature == f_name]

  # https://www.jianshu.com/p/24bbf44e4fa2
  specific_f %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(.data$component) %>%
    dplyr::summarize(
      count = list(count_components(.data$min, .data$max, .data$label, feature_df)),
      sample = list(samps)
    ) %>%
    tidyr::unnest(cols = c("count", "sample")) %>%
    tidyr::pivot_wider(names_from = "sample", values_from = "count", values_fill = list(count = 0L)) %>% # Should not have NA value, but take case of it with 0
    data.table::as.data.table()

  # Same result as above
  # But I don't know which is more efficient
  #
  # specific_f %>%
  #   dplyr::as_tibble() %>%
  #   dplyr::group_by(.data$component) %>%
  #   dplyr::summarize(
  #     count = paste(count_components(.data$min, .data$max, .data$label, feature_df), collapse = ","),
  #   ) %>%
  #   tidyr::separate(.data$count, samps, sep = ",", convert = TRUE) %>%
  #   data.table::as.data.table()
}

handle_sex <- function(abs_profiles) {
  # only works for feature 'CNCP' and 'CN' and "NChrV" and "BoChr"
  sex <- getOption("sigminer.sex", default = "female")
  cn_max <- getOption("sigminer.copynumber.max", default = NA_integer_)
  stopifnot(is.character(sex) | is.data.frame(sex), is.na(cn_max) | is.numeric(cn_max))

  if (is.character(sex)) {
    if (sex == "male") {
      data <- data.table::rbindlist(abs_profiles)
      data <- data %>%
        dplyr::as_tibble() %>%
        dplyr::mutate(
          segVal = ifelse(.data$chromosome %in% c("chrX", "chrY"), 2L * .data$segVal, .data$segVal)
        ) %>%
        data.table::as.data.table()
    } else {
      return(abs_profiles)
    }
  } else {
    sex <- sex %>%
      dplyr::select(c("sample", "sex"))

    data <- data.table::rbindlist(abs_profiles)
    data <- data %>%
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

  if (!is.na(cn_max)) {
    data$segVal[data$segVal > cn_max] <- as.integer(cn_max)
  }

  res <- split(data, by = "sample")
  return(res)
}
