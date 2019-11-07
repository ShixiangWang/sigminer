getSegsize <- function(abs_profiles) {
  out <- c()
  samps <- names(abs_profiles)
  for (i in samps)
  {
    segTab <- abs_profiles[[i]]
    segTab$segVal[segTab$segVal < 0] <- 0
    seglen <- segTab$end - segTab$start
    seglen <- seglen[seglen > 0]
    out <-
      rbind(out, cbind(ID = rep(i, length(seglen)), value = seglen))
  }
  rownames(out) <- NULL
  data.frame(out, stringsAsFactors = F)
}

getBPnum <- function(abs_profiles, chrlen) {
  out <- c()
  samps <- names(abs_profiles)
  for (i in samps)
  {
    segTab <- abs_profiles[[i]]

    chrs <- unique(segTab$chromosome)

    allBPnum <- c()
    for (c in chrs)
    {
      currseg <- segTab[chromosome == c, ]
      intervals <-
        seq(1, chrlen[chrlen[, 1] == c, 2] + 10000000, 10000000)
      res <- tryCatch(
        hist(currseg$end[-nrow(currseg)],
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
      res <-

        allBPnum <- c(allBPnum, res)
    }
    out <-
      rbind(out, cbind(ID = rep(i, length(allBPnum)), value = allBPnum))
  }
  rownames(out) <- NULL
  data.frame(out, stringsAsFactors = F)
}

getOscilation <- function(abs_profiles) {
  out <- c()
  samps <- names(abs_profiles)
  for (i in samps)
  {
    segTab <- abs_profiles[[i]]

    chrs <- unique(segTab$chromosome)

    oscCounts <- c()
    for (c in chrs)
    {
      currseg <- segTab[chromosome == c][["segVal"]]
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
      }
    }
    out <-
      rbind(out, cbind(ID = rep(i, length(oscCounts)), value = oscCounts))
    if (length(oscCounts) == 0) {
      out <- rbind(out, cbind(ID = i, value = 0))
    }
  }
  rownames(out) <- NULL
  data.frame(out, stringsAsFactors = F)
}

getCentromereDistCounts <-
  function(abs_profiles, centromeres, chrlen) {
    out <- c()
    samps <- names(abs_profiles)
    for (i in samps)
    {
      segTab <- abs_profiles[[i]]

      chrs <- unique(segTab$chromosome)

      all_dists <- c()
      for (c in chrs)
      {
        if (nrow(segTab) > 1) {
          starts <- segTab[chromosome == c][["start"]][-1]
          segstart <- segTab[chromosome == c][["start"]][1]
          ends <- segTab[chromosome == c][["end"]]
          segend <- ends[length(ends)]
          ends <- ends[-length(ends)]
          # centstart <-
          #     as.numeric(centromeres[substr(centromeres[, 2], 4, 5) == c, 3])
          # centend <-
          #     as.numeric(centromeres[substr(centromeres[, 2], 4, 5) == c, 4])
          # chrend <- chrlen[substr(chrlen[, 1], 4, 5) == c, 2]

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
        }
      }
      if (nrow(all_dists) > 0) {
        out <- rbind(out, cbind(ID = i, value = all_dists[, 1]))
      }
    }
    rownames(out) <- NULL
    out <- data.frame(out, stringsAsFactors = F)
    out <- stats::na.omit(out)
    out
  }


getChangepointCN <- function(abs_profiles) {
  out <- c()
  samps <- names(abs_profiles)
  for (i in samps)
  {
    segTab <- abs_profiles[[i]]

    segTab$segVal[segTab$segVal < 0] <- 0
    chrs <- unique(segTab$chromosome)
    allcp <- c()
    for (c in chrs)
    {
      currseg <- segTab[chromosome == c][["segVal"]]
      allcp <-
        c(allcp, abs(currseg[-1] - currseg[-length(currseg)]))
    }
    if (length(allcp) == 0) {
      allcp <- 0 # if there are no changepoints
    }
    out <-
      rbind(out, cbind(ID = rep(i, length(allcp)), value = allcp))
  }
  rownames(out) <- NULL
  data.frame(out, stringsAsFactors = F)
}

getCN <- function(abs_profiles) {
  out <- c()
  samps <- names(abs_profiles)
  for (i in samps)
  {
    segTab <- abs_profiles[[i]]
    segTab$segVal[segTab$segVal < 0] <- 0
    cn <- segTab$segVal
    out <-
      rbind(out, cbind(ID = rep(i, length(cn)), value = cn))
  }
  rownames(out) <- NULL
  data.frame(out, stringsAsFactors = F)
}
