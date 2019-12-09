getSegsize <- function(abs_profiles) {
  segsize = purrr::map_df(abs_profiles, function(x) {
    x$segsize = x$end - x$start + 1
    x[, c("sample", "segsize"), with = FALSE]
  })
  colnames(segsize) = c("ID", "value")
  segsize
}

getBPnum <- function(abs_profiles, chrlen) {
  res = purrr::map_df(abs_profiles, function(x, chrlen) {

    calcBPnum = function(df, c, chrlen) {
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

    x = x %>%
      dplyr::as_tibble() %>%
      dplyr::group_by(.data$chromosome) %>%
      tidyr::nest() %>%
      dplyr::mutate(value = purrr::map2(
        data, .data$chromosome, calcBPnum, chrlen = chrlen)) %>%
      dplyr::ungroup() %>%
      dplyr::select(c("value"))
    data.table::data.table(
      value = purrr::reduce(x$value, c)
    )
  }, chrlen = chrlen, .id = "ID")

  return(res)
}

getOscilation <- function(abs_profiles) {
  oscCounts = purrr::map_df(abs_profiles, function(x) {
    x = x %>%
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
            return(NULL)
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
  function(abs_profiles, centromeres, chrlen) {

    calcArmBP = function(df, c, centromeres, chrlen) {
      all_dists <- c()

      if (nrow(df) > 1) {
        starts <- df$start[-1]
        segstart <- df$start[1]
        ends <- df$end
        segend <- ends[length(ends)]
        ends <- ends[-length(ends)]

        centstart <- centromeres[centromeres$chrom == c, 2]
        centend <- centromeres[centromeres$chrom == c, 3]
        chrend <- chrlen[chrlen$chrom == c, 2]

        ndist <-
          cbind(rep(NA, length(starts)), rep(NA, length(starts)))
        ndist[starts <= centstart, 1] <-
          (centstart - starts[starts <= centstart]) / (centstart - segstart) * -1
        ndist[starts >= centend, 1] <-
          (starts[starts >= centend] - centend) / (segend - centend)
        ndist[ends <= centstart, 2] <-
          (centstart - ends[ends <= centstart]) / (centstart - segstart) * -1
        ndist[ends >= centend, 2] <-
          (ends[ends >= centend] - centend) / (segend - centend)
        ndist <- apply(ndist, 1, min)

        all_dists <- rbind(all_dists, sum(ndist > 0))
        all_dists <- rbind(all_dists, sum(ndist <= 0))

        if (nrow(all_dists) > 0) {
          return(as.integer(all_dists[, 1]))
        } else {
          return(rep(0L, 2))
        }

      } else {
        return(rep(0L, 2))
      }
    }

    res = purrr::map_df(abs_profiles, function(x, centromeres, chrlen) {

      x = x %>%
        dplyr::as_tibble() %>%
        dplyr::group_by(.data$chromosome) %>%
        tidyr::nest() %>%
        dplyr::mutate(value = purrr::map2(
          data, .data$chromosome, calcArmBP, centromeres = centromeres, chrlen = chrlen)) %>%
        dplyr::ungroup() %>%
        dplyr::select(c("value"))
      data.table::data.table(
        value = purrr::reduce(x$value, c)
      )
    }, centromeres = centromeres, chrlen = chrlen, .id = "ID")

    return(res)
}

getChangepointCN <- function(abs_profiles) {

  cp = purrr::map_df(abs_profiles, function(x) {
    x = x %>%
      dplyr::as_tibble() %>%
      dplyr::group_by(.data$chromosome) %>%
      tidyr::nest() %>%
      dplyr::mutate(value = purrr::map(
        data, function(df) {
          if (nrow(df) <= 1) {
            return(NULL)
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
  cn = purrr::map_df(abs_profiles, function(x) {
    x[, c("sample", "segVal"), with = FALSE]
  })
  colnames(cn) = c("ID", "value")
  cn
}

# Number of Chromosome with CNV
getNChrV <- function(abs_profiles, genome_build="hg38") {
  genome_build = match.arg(genome_build, choices = c("hg19", "hg38"))

  if (genome_build %in% c("hg19", "hg38")) {
    autosome = paste0("chr", 1:22)
  }

  cn = purrr::map_df(abs_profiles, function(x) {
    x = x[x$chromosome %in% autosome, c("sample", "chromosome", "segVal"), with = FALSE][x$segVal != 2]
    data.table::data.table(ID = x$sample[1], value = length(unique(x$chromosome)))
  })
  cn
}
