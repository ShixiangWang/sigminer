#========================================================
# internal analysis functions, use get_ as start
#========================================================

# License Part ------------------------------------------------------------

# MIT License
#
# Copyright (c) [2018] [Geoffrey Macintyre]
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
#     The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.



# Get feature distributions -----------------------------------------------

#' @title  Derive copy number feature distributions
#' @description This function summarise each copy-number profile using a number of different
#' feature distributions: sigment size, breakpoint number (per ten megabase), change-point copy-number,
#' segment copy-number, breakpoint number (per chromosome arm), length of segments with oscilating
#' copy-number.
#' @param CN_data a `list` contains multiple `data.frame`s (recommended),
#' each `data.frame` stores copy-number profile for one sample with 'chromosome', 'start', 'end' and
#' 'segVal' these four necessary columns. Of note, 'segVal' column shoule be absolute copy number values.
#' @param cores number of compute cores to run this task.
#' You can use [parallel::detectCores()] function to check how
#' many cores you can use.
#' @inheritParams read_copynumber
#' @author Geoffrey Macintyre, Shixiang Wang
#' @return a `list` contains six copy number feature distributions.
#' @export
#' @family internal calculation function series

get_features = function(CN_data,
                        cores = 1,
                        genome_build = c("hg19", "hg38")) {
  genome_build = match.arg(genome_build)
  # get chromosome lengths and centromere locations
  if (genome_build == "hg19") {
    data("chromsize.hg19",
         package = "VSHunter",
         envir = environment())
    data("centromeres.hg19",
         package = "VSHunter",
         envir = environment())
    chrlen = chromsize.hg19
    centromeres = centromeres.hg19
  } else {
    data("chromsize.hg38",
         package = "VSHunter",
         envir = environment())
    data("centromeres.hg38",
         package = "VSHunter",
         envir = environment())
    chrlen = chromsize.hg38
    centromeres = centromeres.hg38
  }

  # only keep 1:22 and x, y
  chrlen = chrlen[chrlen$chrom %in% centromeres$chrom,]
  if (cores > 1) {
    #require(foreach)
    requireNamespace("foreach", quietly = TRUE)
    #doMC::registerDoMC(cores)
    doParallel::registerDoParallel(cores = cores)

    temp_list = foreach::foreach(i = 1:6) %dopar% {
      if (i == 1) {
        list(segsize = getSegsize(CN_data))
      } else if (i == 2) {
        list(bp10MB = getBPnum(CN_data, chrlen))
      } else if (i == 3) {
        list(osCN = getOscilation(CN_data))
      } else if (i == 4) {
        list(bpchrarm = getCentromereDistCounts(CN_data, centromeres, chrlen))
      } else if (i == 5) {
        list(changepoint = getChangepointCN(CN_data))
      } else {
        list(copynumber = getCN(CN_data))
      }

    }
    unlist(temp_list, recursive = FALSE)
  } else {
    segsize <- getSegsize(CN_data)
    bp10MB <- getBPnum(CN_data, chrlen)
    osCN <- getOscilation(CN_data)
    bpchrarm <-
      getCentromereDistCounts(CN_data, centromeres, chrlen)
    changepoint <- getChangepointCN(CN_data)
    copynumber <- getCN(CN_data)

    list(
      segsize = segsize,
      bp10MB = bp10MB,
      osCN = osCN,
      bpchrarm = bpchrarm,
      changepoint = changepoint,
      copynumber = copynumber
    )
  }

}




# Get mixture model components --------------------------------------------

#' @title Fit optimal number of mixture model components
#' @description Apply mixture modelling to breakdown each feature distribution into mixtures
#' of Gaussian or mixtures of Poison distributions using **flexmix** package.
#'
#' @param CN_features a `list` generate from [get_features()] function.
#' @param seed seed number.
#' @param min_comp minimal number of components to fit, default is 2.
#' @param max_comp maximal number of components to fit, default is 10.
#' @param min_prior minimal prior value, default is 0.001.
#' Details about custom setting please refer to **flexmix** package.
#' @param model_selection model selection strategy, default is 'BIC'.
#' Details about custom setting please refer to **flexmix** package.
#' @param nrep number of run times for each value of component,
#' keep only the solution with maximum likelihood.
#' @param niter maximal number of iteration to achive converge.
#' @inheritParams get_features
#' @author Geoffrey Macintyre, Shixiang Wang
#' @return a `list` contain `flexmix` object of copy-number features.
#' @export
#' @family internal calculation function series
get_components = function(CN_features,
                            seed = 123456,
                            min_comp = 2,
                            max_comp = 10,
                            min_prior = 0.001,
                            model_selection = "BIC",
                            nrep = 1,
                            niter = 1000,
                            cores = 1) {

  featsToFit = seq(1, 6)

  if (cores > 1) {
    #require(foreach)
    requireNamespace("foreach", quietly = TRUE)
    #doMC::registerDoMC(cores)
    doParallel::registerDoParallel(cores = cores)

    temp_list = foreach(i = 1:6) %dopar% {
      if (i == 1 & i %in% featsToFit) {
        message("Fit feature: Segment size")
        dat <- as.numeric(CN_features[["segsize"]][, 2])
        list(
          segsize = fitComponent(
            dat,
            seed = seed,
            model_selection = model_selection,
            min_prior = min_prior,
            niter = niter,
            nrep = nrep,
            min_comp = min_comp,
            max_comp = max_comp
          )
        )

      } else if (i == 2 & i %in% featsToFit) {
        message("Fit feature: Breakpoint count per 10 Mb")
        dat <- as.numeric(CN_features[["bp10MB"]][, 2])
        list(
          bp10MB = fitComponent(
            dat,
            dist = "pois",
            seed = seed,
            model_selection = model_selection,
            min_prior = min_prior,
            niter = niter,
            nrep = nrep,
            min_comp = min_comp,
            max_comp = max_comp
          )
        )

      } else if (i == 3 & i %in% featsToFit) {
        message("Fit feature: Length of oscillating copy-number chain")
        dat <- as.numeric(CN_features[["osCN"]][, 2])
        list(
          osCN = fitComponent(
            dat,
            dist = "pois",
            seed = seed,
            model_selection = model_selection,
            min_prior = min_prior,
            niter = niter,
            nrep = nrep,
            min_comp = min_comp,
            max_comp = max_comp
          )
        )

      } else if (i == 4 & i %in% featsToFit) {
        message("Fit feature: Breakpoint count per arm")
        dat <- as.numeric(CN_features[["bpchrarm"]][, 2])
        list(
          bpchrarm = fitComponent(
            dat,
            dist = "pois",
            seed = seed,
            model_selection = model_selection,
            min_prior = min_prior,
            niter = niter,
            nrep = nrep,
            min_comp = min_comp,
            max_comp = max_comp
          )
        )

      } else if (i == 5 & i %in% featsToFit) {
        message("Fit feature: Copy number change")
        dat <- as.numeric(CN_features[["changepoint"]][, 2])
        list(
          changepoint = fitComponent(
            dat,
            seed = seed,
            model_selection = model_selection,
            min_prior = min_prior,
            niter = niter,
            nrep = nrep,
            min_comp = min_comp,
            max_comp = max_comp
          )
        )

      } else if (i == 6 & i %in% featsToFit) {
        message("Fit feature: Absolute copy number")
        dat <- as.numeric(CN_features[["copynumber"]][, 2])
        list(
          copynumber = fitComponent(
            dat,
            seed = seed,
            model_selection = model_selection,
            nrep = nrep,
            min_comp = min_comp,
            max_comp = max_comp,
            min_prior = 0.005,
            niter = 2000
          )
        )

      }

    }
    unlist(temp_list, recursive = FALSE)
  } else {
    dat <- as.numeric(CN_features[["segsize"]][, 2])
    message("Fit feature: Segment size")
    segsize_mm <-
      fitComponent(
        dat,
        seed = seed,
        model_selection = model_selection,
        min_prior = min_prior,
        niter = niter,
        nrep = nrep,
        min_comp = min_comp,
        max_comp = max_comp
      )

    dat <- as.numeric(CN_features[["bp10MB"]][, 2])
    message("Fit feature: Breakpoint count per 10 Mb")
    bp10MB_mm <-
      fitComponent(
        dat,
        dist = "pois",
        seed = seed,
        model_selection = model_selection,
        min_prior = min_prior,
        niter = niter,
        nrep = nrep,
        min_comp = min_comp,
        max_comp = max_comp
      )

    dat <- as.numeric(CN_features[["osCN"]][, 2])
    message("Fit feature: Length of oscillating copy-number chain")
    osCN_mm <-
      fitComponent(
        dat,
        dist = "pois",
        seed = seed,
        model_selection = model_selection,
        min_prior = min_prior,
        niter = niter,
        nrep = nrep,
        min_comp = min_comp,
        max_comp = max_comp
      )

    dat <- as.numeric(CN_features[["bpchrarm"]][, 2])
    message("Fit feature: Breakpoint count per arm")
    bpchrarm_mm <-
      fitComponent(
        dat,
        dist = "pois",
        seed = seed,
        model_selection = model_selection,
        min_prior = min_prior,
        niter = niter,
        nrep = nrep,
        min_comp = min_comp,
        max_comp = max_comp
      )

    dat <- as.numeric(CN_features[["changepoint"]][, 2])
    message("Fit feature: Copy number change")
    changepoint_mm <-
      fitComponent(
        dat,
        seed = seed,
        model_selection = model_selection,
        min_prior = min_prior,
        niter = niter,
        nrep = nrep,
        min_comp = min_comp,
        max_comp = max_comp
      )

    dat <- as.numeric(CN_features[["copynumber"]][, 2])
    message("Fit feature: Absolute copy number")
    copynumber_mm <-
      fitComponent(
        dat,
        seed = seed,
        model_selection = model_selection,
        nrep = nrep,
        min_comp = min_comp,
        max_comp = max_comp,
        min_prior = 0.005,
        niter = 2000
      )

    list(
      segsize = segsize_mm,
      bp10MB = bp10MB_mm,
      osCN = osCN_mm,
      bpchrarm = bpchrarm_mm,
      changepoint = changepoint_mm,
      copynumber = copynumber_mm
    )
  }
}



# Get sample-by-component matrix ------------------------------------------

#' @title Generate a sample-by-component matrix
#' @description Generate a sample-by-component matrix representing the sum of
#' posterior probabilities of each copy-number event being assigned to each component.
#' @param CN_features a `list` contains six copy number feature distributions,
#' obtain this from [get_features()] function.
#' @param all_components a `list` contain `flexmix` object of copy-number features, obtain this
#' from [get_components] function or use pre-compiled components data which come from CNV signature paper
#' https://www.nature.com/articles/s41588-018-0179-8 (set this parameter as `NULL`).
#' @inheritParams get_features
#' @param rowIter step size of iteration for rows of ech CNV feature.
#' @author Geoffrey Macintyre, Shixiang Wang
#' @return a numeric sample-by-component `matrix`
#' @importFrom utils data download.file str
#' @export
#' @family internal calculation function seriess
get_matrix = function(CN_features,
                      all_components = NULL,
                      cores = 1,
                      rowIter = 1000)
{
  if (is.null(all_components))
  {
    message(
      "About reference components\n   more detail please see https://github.com/ShixiangWang/absoluteCNVdata"
    )
    if (!file.exists("Nat_Gen_component_parameters.rds")) {
      message(
        "Nat_Gen_component_parameters.rds doesnot exist, will download reference components."
      )
      download.file(url = "https://github.com/ShixiangWang/absoluteCNVdata/raw/master/component_parameters.rds",
                    destfile = "Nat_Gen_component_parameters.rds")
    }
    all_components <-
      readRDS("Nat_Gen_component_parameters.rds")
  }

  full_mat <- cbind(
    calculateSumOfPosteriors(CN_features[["segsize"]],
                             all_components[["segsize"]],
                             "segsize",
                             cores = cores),
    calculateSumOfPosteriors(CN_features[["bp10MB"]],
                             all_components[["bp10MB"]],
                             "bp10MB",
                             cores = cores),
    calculateSumOfPosteriors(CN_features[["osCN"]],
                             all_components[["osCN"]],
                             "osCN",
                             cores = cores),
    calculateSumOfPosteriors(CN_features[["changepoint"]],
                             all_components[["changepoint"]],
                             "changepoint",
                             cores = cores),
    calculateSumOfPosteriors(CN_features[["copynumber"]],
                             all_components[["copynumber"]],
                             "copynumber",
                             cores = cores),
    calculateSumOfPosteriors(CN_features[["bpchrarm"]],
                             all_components[["bpchrarm"]],
                             "bpchrarm",
                             cores = cores)
  )


  rownames(full_mat) <- unique(CN_features[["segsize"]][, 1])
  full_mat[is.na(full_mat)] <- 0
  full_mat
}

#' Calculate length fraction profile of copy number
#'
#' @inheritParams get_features
#' @inheritParams read_copynumber
#' @author Shixiang Wang <w_shixiang@163.com>
#' @return a data frame
#' @export
#' @family internal calculation function series
get_LengthFraction = function(CN_data,
                              genome_build = c("hg19", "hg38"),
                              seg_cols = c("Chromosome", "Start.bp", "End.bp", "modal_cn"),
                              samp_col = "sample") {
  stopifnot(is.list(CN_data) | is.data.frame(CN_data))
  genome_build = match.arg(genome_build)

  if (inherits(CN_data, "list")) {
    segTab = base::Reduce(rbind, CN_data)
    segTab$sample = base::rep(x = names(CN_data),
                              times = sapply(CN_data, function(x)
                                nrow(x)))
  } else {
    segTab = CN_data[, c(seg_cols, samp_col)]
    if (ncol(segTab) == 5) {
      colnames(segTab) = c("chromosome", "start", "end", "segVal", "sample")
    } else if (ncol(segTab) == 4) {
      colnames(segTab) = c("chromosome", "start", "end", "sample")
    } else {
      stop(
        "If input is a data.frame, must have 4 necessary columns (chr, start, end, sample) and 1 optional column (segVal)."
      )
    }
  }

  # unify chromosome column
  segTab$chromosome = as.character(segTab$chromosome)
  segTab$chromosome = sub(
    pattern = "chr",
    replacement = "chr",
    x = segTab$chromosome,
    ignore.case = TRUE
  )
  if (any(!grepl("chr", segTab$chromosome))) {
    segTab$chromosome[!grepl("chr", segTab$chromosome)] = paste0("chr", segTab$chromosome[!grepl("chr", segTab$chromosome)])
  }
  if (any(grepl("chr23", segTab$chromosome))) {
    warning("'23' is not a supported chromosome, related rows will be discarded.")
    segTab = segTab[!grepl("chr23", segTab$chromosome),]
  }

  valid_chr = c(paste0("chr", 1:22), "chrX", "chrY")
  if (!all(segTab$chromosome %in% valid_chr)) {
    message("Filter some invalid segments... (not as 1:22 and X, Y)")
    segTab = base::subset(segTab, chromosome %in% valid_chr)
  }

  arm_data = get_ArmLocation(genome_build)
  # logical operation

  assign_df = data.frame(
    location = vector("character", nrow(segTab)),
    annotation = vector("character", nrow(segTab)),
    fraction = vector("numeric", nrow(segTab)),
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(segTab)) {
    x = segTab[i,]

    # locate chromosome
    arm_loc = base::subset(arm_data, chrom == x$chromosome)

    y = c(x$start, x$end)
    if (y[2] <= arm_loc$p_end & y[1] >= arm_loc$p_start) {
      location = paste0(sub("chr", "", arm_loc$chrom), "p")
      annotation = "short arm"
      fraction = (y[2] - y[1] + 1) / (arm_loc$p_end - arm_loc$p_start + 1)
    } else if (y[2] <= arm_loc$q_end &
               y[1] >= arm_loc$q_start) {
      location = paste0(sub("chr", "", arm_loc$chrom), "q")
      annotation = "long arm"
      fraction = (y[2] - y[1] + 1) / (arm_loc$q_end - arm_loc$q_start + 1)
    } else if (y[1] >= arm_loc$p_start &
               y[1] <= arm_loc$p_end &
               y[2] >= arm_loc$q_start & y[2] <= arm_loc$q_end) {
      location = paste0(sub("chr", "", arm_loc$chrom), "pq") # across p and q arm
      annotation = "across short and long arm"
      fraction = 2 * ((y[2] - y[1] + 1) / arm_loc$total_size)
    } else if (y[1] < arm_loc$p_end & y[2] < arm_loc$q_start) {
      location = paste0(sub("chr", "", arm_loc$chrom), "p")
      annotation = "short arm intersect with centromere region"
      # only calculate region does not intersect
      fraction = (y[2] - y[1] + 1 - (y[2] - arm_loc$p_end)) / (arm_loc$p_end - arm_loc$p_start + 1)
    } else if (y[1] > arm_loc$p_end &
               y[1] < arm_loc$q_start & y[2] > arm_loc$q_start) {
      location = paste0(sub("chr", "", arm_loc$chrom), "q")
      annotation = "long arm intersect with centromere region"
      # only calculate region does not intersect
      fraction = (y[2] - y[1] + 1 - (y[1] - arm_loc$q_start)) / (arm_loc$q_end - arm_loc$q_start + 1)
    } else {
      location = "Unknown"
      annotation = "Unknow segment, locate in centromere region"
      fraction = NA_real_
    }

    assign_df[i,] = c(location, annotation, fraction)
    #assign_df = base::rbind(assign_df, c(location, annotation, percentage))
  }
  res = base::cbind(segTab, assign_df)
  res$fraction = as.numeric(res$fraction)
  res
}


#' Get chromosome arm location
#' @inheritParams read_copynumber
#' @export
#' @family internal calculation function series

get_ArmLocation = function(genome_build = c("hg19", "hg38")) {
  # get chromosome lengths and centromere locations
  if (genome_build == "hg19") {
    data("chromsize.hg19",
         package = "VSHunter",
         envir = environment())
    data("centromeres.hg19",
         package = "VSHunter",
         envir = environment())
    chrlen = chromsize.hg19
    centromeres = centromeres.hg19
  } else {
    data("chromsize.hg38",
         package = "VSHunter",
         envir = environment())
    data("centromeres.hg38",
         package = "VSHunter",
         envir = environment())
    chrlen = chromsize.hg38
    centromeres = centromeres.hg38
  }

  # only keep 1:22 and x, y
  chrlen = chrlen[chrlen$chrom %in% centromeres$chrom,]

  # sort
  chrlen = chrlen[order(chrlen$chrom),]
  centromeres = centromeres[order(centromeres$chrom),]

  # compute and get results
  res = data.frame(
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

  i = 1
  for (chr in chrlen$chrom) {
    chrom = chr

    # p
    p_start = 1
    p_end = centromeres$left.base[centromeres$chrom == chr]
    p_length = p_end - p_start + 1

    # q
    q_start = centromeres$right.base[centromeres$chrom == chr]
    q_end = chrlen$size[chrlen$chrom == chr]
    q_length = q_end - q_start + 1

    total_size = chrlen$size[chrlen$chrom == chr]

    res[i, 1] = as.character(chrom)
    res[i, 2:8] = c(p_start,
                    p_end,
                    p_length,
                    q_start,
                    q_end,
                    q_length,
                    total_size)
    i = i + 1
  }

  res
}




# Global variables --------------------------------------------------------


# utils::globalVariables(
#   c(
#     "centromeres.hg19",
#     "centromeres.hg38",
#     "chromsize.hg19",
#     "chromsize.hg38",
#     "feat",
#     "i"
#   )
# )
