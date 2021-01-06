
# Extraction helpers ------------------------------------------------------

process_solution <- function(slist, catalogue_matrix,
                             report_integer_exposure = FALSE,
                             only_core_stats = FALSE) {
  on.exit(invisible(gc()), add = TRUE)
  send_info("Normalizing solutions to get Signature and Exposure.")
  out <- purrr::map(slist, .f = normalize_solution) %>%
    setNames(paste0("Run", seq_along(slist)))
  sn <- ncol(out[[1]]$Signature)
  send_info("Start processing ", sn, " signatures.")

  # If there are hyper-mutant records, collapse the samples into true samples
  if (any(grepl("_\\[hyper\\]_", colnames(catalogue_matrix)))) {
    send_info("Hyper-mutant records detected, revert them to true samples.")
    catalogue_matrix <- collapse_hyper_records(catalogue_matrix)
    purrr::map(seq_along(out), function(i) {
      out[[i]]$Exposure <<- collapse_hyper_records(out[[i]]$Exposure)
    })
    purrr::map(seq_along(out), function(i) {
      out[[i]]$H <<- collapse_hyper_records(out[[i]]$H)
    })
    send_success("Reverted.")
  }

  # If just one run, skip the following steps
  if (length(out) >= 2) {
    send_info("Running clustering with match algorithm for signature assignment in different NMF runs.")
    run_pairs <- combn(names(out), 2, simplify = FALSE)
    run_pairs_name <- purrr::map_chr(run_pairs, ~ paste(., collapse = "_"))
    send_success("Run pairs obtained.")
    # Get similarity distance
    pair_dist <- purrr::map(
      run_pairs,
      ~ get_cosine_distance(out[[.[1]]]$Signature, out[[.[2]]]$Signature)
    ) %>%
      setNames(run_pairs_name)
    send_success("Run pairwise similarity distance calculated.")
    pair_dist_mean <- purrr::map_dbl(pair_dist, mean) %>%
      setNames(run_pairs_name)
    match_list <- purrr::map2(pair_dist, run_pairs, get_matches)
    send_success("Match list obtained.")

    # Order the result by mean distance
    res_orders <- order(pair_dist_mean)
    run_pairs <- run_pairs[res_orders]
    pair_dist <- pair_dist[res_orders]
    pair_dist_mean <- pair_dist_mean[res_orders]
    match_list <- match_list[res_orders]
    send_success("Match list ordered by distance.")

    # Do clustering with match
    clusters <- clustering_with_match(match_list, n = length(out))
    send_success("Clustering done.")

    # 按 match cluster 排序
    # 样本排序也确保对齐
    send_info("Ordering data based on the clustering result.")
    samp_order <- colnames(out$Run1$Exposure)
    out <- purrr::map2(clusters, out, function(x, y, samp_order) {
      y$Signature <- y$Signature[, x, drop = FALSE]
      y$Exposure <- y$Exposure[x, samp_order, drop = FALSE]
      y
    }, samp_order = samp_order)
    send_success("Done.")
  } else {
    run_pairs <- NA
    pair_dist <- NA
    pair_dist_mean <- NA
    match_list <- NA
    clusters <- data.table::data.table(
      Run1 = paste0("S", seq_len(ncol(out$Run1$Signature)))
    )
  }

  send_info("Calculating stats for signatures and samples.")
  # 生成统计量
  # Remember the new order above
  # 分为 signature 和 样本两种，取每个度量的最大、最小值、平均值以及 SD
  # 重构相似性，轮廓系数，
  # 聚类平均相似距离, 平均错误，Exposure 相关性等
  send_info("Getting signature stats.")
  stat_sigs <- get_stat_sigs(out)
  send_info("Getting sample stats.")
  stat_samps <- get_stat_samps(out, mat = catalogue_matrix, only_core_stats = only_core_stats)
  send_success("Done.")

  send_info("Outputing extraction result and corresponding stats.")

  Signature <- stat_sigs$signature
  Exposure <- stat_samps$exposure
  stats_signature <- stat_sigs$stats
  stats_sample <- stat_samps$stats

  # Order by contribution and set signature names
  new_order <- order(colSums(Signature$signature_mean), decreasing = TRUE)

  Signature <- lapply(Signature, function(s) {
    s <- s[, new_order, drop = FALSE]
    rownames(s) <- rownames(catalogue_matrix)
    colnames(s) <- paste0("Sig", seq_along(new_order))
    s
  })

  Exposure <- lapply(Exposure, function(e) {
    e <- e[new_order, , drop = FALSE]
    colnames(e) <- colnames(catalogue_matrix)
    rownames(e) <- paste0("Sig", seq_along(new_order))
    e
  })

  # Update the order for stats_signature and set signature names
  stats_signature <- stats_signature[new_order, ]
  stats_signature$signature <- paste0("Sig", seq_along(new_order))

  # Merge the stats and generate measures for selecting signature number
  stats <- data.frame(
    signature_number = nrow(stats_signature),
    silhouette = mean(stats_signature$silhouette),
    sample_cosine_distance = mean(stats_sample$cosine_distance_mean),
    L1_error = mean(stats_sample$L1_mean),
    L2_error = mean(stats_sample$L2_mean),
    exposure_positive_correlation = mean(stats_signature$expo_pos_cor_mean),
    signature_similarity_within_cluster = mean(stats_signature$similarity_mean),
    signature_similarity_across_cluster = mean(stats_signature$cross_similarity_mean),
    silhouette_sample = if (isTRUE(only_core_stats)) NA else mean(stats_sample$silhouette),
    # 不同 runs 同一样本看作一个 cluster
    # 展示的是样本间的区分度
    stringsAsFactors = FALSE
  )
  send_success("Data.frame for stats generated.")

  # Generate Signature object
  object <- tf_signature(
    Signature$signature_mean,
    Exposure$exposure_mean,
    used_runs = length(out),
    catalogue_matrix = if (report_integer_exposure) catalogue_matrix else NULL
  )
  send_success("Signature object generated.")

  list(
    object = object,
    stats = stats,
    stats_signature = stats_signature,
    stats_sample = stats_sample,
    signature = Signature,
    exposure = Exposure
  )
}

normalize_solution <- function(solution) {
  # solution is a NMF fit result
  on.exit(invisible(gc()), add = TRUE)

  out <- c(helper_scale_nmf_matrix(solution$W, solution$H, solution$K, handle_cn = FALSE),
    list(W = solution$W, H = solution$H),
    KLD = solution$KLD
  )
  colnames(out$Signature) <- colnames(out$W) <- rownames(out$Exposure) <- rownames(out$H) <- paste0("S", seq_len(solution$K))
  out
}

# Transform data into Signature object
tf_signature <- function(s, e, used_runs, catalogue_matrix = NULL) {
  on.exit(invisible(gc()), add = TRUE)

  # If total_records is not NULL
  # generate integer counts based on resampling
  s.norm <- apply(s, 2, function(x) x / sum(x, na.rm = TRUE))
  e.norm <- apply(e, 2, function(x) x / sum(x, na.rm = TRUE))
  # When only one signature
  if (!is.matrix(e.norm)) {
    e.norm <- matrix(e.norm,
      nrow = 1,
      dimnames = list(rownames(e), names(e.norm))
    )
  }

  if (!is.null(catalogue_matrix)) {
    set.seed(123, kind = "L'Ecuyer-CMRG")
    s2 <- purrr::map2(
      as.data.frame(s),
      round((colSums(s) / sum(colSums(s))) * sum(catalogue_matrix)),
      simulate_catalogue
    ) %>%
      dplyr::as_tibble() %>%
      as.matrix()
    rownames(s2) <- rownames(s)
    colnames(s2) <- colnames(s)
    s <- s2

    if (nrow(e) < 2) {
      e2 <- matrix(colSums(catalogue_matrix), nrow = 1)
    } else {
      set.seed(123, kind = "L'Ecuyer-CMRG")
      e2 <- purrr::map2(
        as.data.frame(e),
        colSums(catalogue_matrix),
        simulate_catalogue
      ) %>%
        dplyr::as_tibble() %>%
        as.matrix()
    }
    rownames(e2) <- rownames(e)
    colnames(e2) <- colnames(e)
    e <- e2
  }

  obj <- list(
    Signature = s,
    Signature.norm = s.norm,
    Exposure = e,
    Exposure.norm = e.norm,
    K = nrow(e)
    # Raw = list()
  )

  class(obj) <- "Signature"
  attr(obj, "used_runs") <- used_runs
  attr(obj, "method") <- "brunet"
  attr(obj, "call_method") <- "NMF with best practice"

  obj
}

get_stat_sigs <- function(runs) {
  on.exit(invisible(gc()), add = TRUE)

  send_info("\t getting exposure array.")
  sig_list <- purrr::map(runs, "Signature")
  dm <- dim(sig_list[[1]])
  l <- length(sig_list)
  sig_array <- array(unlist(sig_list), dim = c(dm, l))
  rm(sig_list)
  invisible(gc())
  # signatures
  sig <- data.frame(
    signature_number = rep(dm[2], dm[2]),
    signature = seq_len(dm[2]) # Rename it outside the function
  )

  send_info("\t summarizing signatures from different runs.")
  s <- get_3d_array_stats(
    sig_array,
    c(
      "signature_mean", "signature_sd",
      "signature_min", "signature_max"
    )
  )

  send_info("\t summarizing sample profile similarity.")
  sim <- get_similarity_stats(
    sig_array,
    n = dm[2],
    ns = c(
      "similarity_mean", "similarity_sd",
      "similarity_min", "similarity_max"
    )
  )

  # cluster silhouette
  send_info("\t summarizing signature-wise profile similarity.")
  cross_sim <- get_similarity_stats(
    sig_array,
    n = dm[2],
    ns = c(
      "cross_similarity_mean", "cross_similarity_sd",
      "cross_similarity_min", "cross_similarity_max"
    ),
    type = "between-cluster"
  )
  b <- 1 - cross_sim$cross_similarity_mean
  a <- 1 - sim$similarity_mean
  sil_width <- data.frame(
    silhouette = (b - a) / pmax(a, b)
  )
  # The minimum KLD in all kept runs
  KLD <- data.frame(
    KLD = rep(min(purrr::map_dbl(runs, "KLD")), length(a))
  )

  # Get exposure correlation between different signatures
  expo_list <- purrr::map(runs, "Exposure")
  dm <- dim(expo_list[[1]])
  expo_array <- array(unlist(expo_list), dim = c(dm, l))
  rm(sig_array, expo_list, runs)
  invisible(gc())
  expo_cor <- get_expo_corr_stat(expo_array)

  list(
    signature = s,
    stats = cbind(sig, KLD, sil_width, expo_cor, sim, cross_sim)
  )
}

get_stat_samps <- function(runs, mat, only_core_stats = FALSE) {
  on.exit(invisible(gc()), add = TRUE)

  send_info("\t getting exposure array.")
  expo_list <- purrr::map(runs, "Exposure")
  dm <- dim(expo_list[[1]]) # the second value of dm indicates how many samples
  l <- length(expo_list)
  expo_array <- array(unlist(expo_list), dim = c(dm, l))
  rm(expo_list)
  invisible(gc())
  # exposures
  send_info("\t summarizing exposures from different runs.")
  e <- get_3d_array_stats(
    expo_array,
    c(
      "exposure_mean", "exposure_sd",
      "exposure_min", "exposure_max"
    )
  )
  # get catalog array
  send_info("\t getting catalog array.")
  W_list <- purrr::map(runs, "W")
  H_list <- purrr::map(runs, "H")
  catalog_list <- purrr::map2(W_list, H_list, ~ .x %*% .y)
  dm2 <- dim(catalog_list[[1]])
  catalog_array <- array(unlist(catalog_list), dim = c(dm2, l))
  rm(runs, expo_array, W_list, H_list, catalog_list)
  invisible(gc())

  if (isFALSE(only_core_stats)) {
    send_info("\t summarizing sample profile similarity.")
    sim <- get_similarity_stats(
      catalog_array,
      n = dm2[2],
      ns = c(
        "similarity_mean", "similarity_sd",
        "similarity_min", "similarity_max"
      )
    )

    # cluster silhouette
    send_info("\t summarizing sample-wise profile similarity.")
    cross_sim <- get_similarity_stats(
      catalog_array,
      n = dm2[2],
      ns = c(
        "cross_similarity_mean", "cross_similarity_sd",
        "cross_similarity_min", "cross_similarity_max"
      ),
      type = "between-cluster"
    )
    b <- 1 - cross_sim$cross_similarity_mean
    a <- 1 - sim$similarity_mean
    sil_width <- data.frame(
      silhouette = (b - a) / pmax(a, b)
    )
  }

  samp <- data.frame(
    signature_number = rep(dm[1], ncol(mat)),
    sample = colnames(mat)
  )

  # 重构 error：包括 重构相似度
  send_info("\t summarizing sample profile reconstruction stats.")
  error <- get_error_stats(
    catalog_array,
    mat
  )

  list(
    exposure = e,
    stats = if (only_core_stats) cbind(samp, error) else cbind(samp, sil_width, error)
  )
}

get_3d_array_stats <- function(x, ns = NULL) {
  on.exit(invisible(gc()), add = TRUE)

  r <- list(
    mean = apply(x, c(1, 2), mean),
    sd = apply(x, c(1, 2), sd),
    min = apply(x, c(1, 2), min),
    max = apply(x, c(1, 2), max)
  )
  if (!is.null(ns)) {
    names(r) <- ns
  }
  r
}

get_similarity_stats <- function(x,
                                 n,
                                 ns = NULL,
                                 type = "within-cluster") {
  on.exit(invisible(gc()), add = TRUE)

  if (type == "within-cluster") {
    d <- lapply(seq_len(n), function(i) {
      x2 <- x[, i, ]
      if (is.null(dim(x2))) {
        x2 <- matrix(x2, ncol = 1)
      }
      mat <- cosineMatrix(x2, x2)
      if (ncol(mat) > 1) {
        mat[upper.tri(mat)]
      } else {
        mat
      }
    })
  } else if (type == "between-cluster") {
    d <- lapply(seq_len(n), function(i) {
      if (dim(x)[2] >= 2) {
        x1 <- x[, i, ]
        x2 <- x[, -i, ]

        if (is.null(dim(x1))) {
          x1 <- matrix(x1, ncol = 1)
        }
        if (is.null(dim(x2))) {
          x2 <- matrix(x2, ncol = 1)
        }

        dim(x2) <- c(dim(x1)[1], prod(dim(x2)) / dim(x1)[1])
        cosineMatrix(x1, x2)
      } else {
        NA
      }
    })
  }
  r <- data.frame(
    mean = sapply(d, mean),
    sd = sapply(d, sd),
    min = sapply(d, min),
    max = sapply(d, max),
    stringsAsFactors = FALSE
  )
  if (!is.null(ns)) {
    colnames(r) <- ns
  }
  r
}

# Quantify the exposure correlation between different signatures with
# Pearson coefficient
get_expo_corr_stat <- function(x) {
  on.exit(invisible(gc()), add = TRUE)

  n <- dim(x)[1] # n signatures
  r <- dim(x)[3] # r runs

  if (n > 1) {
    d <- lapply(seq_len(n), function(i) {
      x1 <- x[i, , , drop = FALSE]
      x2 <- x[-i, , , drop = FALSE]

      rows <- dim(x2)[1]
      corr <- vector("numeric", rows * r)
      # calculate exposure corr in each run
      j <- 1L
      for (row in seq_len(rows)) {
        for (rn in seq_len(r)) {
          corr[j] <- stats::cor(
            x1[1, , rn],
            x2[row, , rn],
            method = "spearman"
          )
          j <- j + 1L
        }
      }

      # 仅关注正相关
      corr <- corr[corr > 0]
      # 乘以一个数量权重
      if (length(corr) == 0) {
        NA
      } else {
        corr <- corr * (length(corr) / (rows * r))
        corr
      }
    })
  } else {
    d <- NA
  }

  r <- data.frame(
    expo_pos_cor_mean = sapply(d, mean),
    expo_pos_cor_sd = sapply(d, sd),
    expo_pos_cor_min = sapply(d, min),
    expo_pos_cor_max = sapply(d, max)
  )
  purrr::map_df(r, ~ ifelse(is.na(.), 0, .)) %>%
    as.data.frame()
}

# The difference between reconstructed catalogs and the original catalog
# 用相似性距离、L1、L2范数
get_error_stats <- function(x, mat) {
  on.exit(invisible(gc()), add = TRUE)

  n <- dim(mat)[2] # n samples
  r <- dim(x)[3] # r runs

  # similarity distance
  d <- lapply(seq_len(n), function(i) {
    x2 <- x[, i, ]
    if (is.null(dim(x2))) {
      x2 <- matrix(x2, ncol = 1)
    }
    1 - cosineMatrix(x2, mat[, i, drop = FALSE])
  })
  # L1 and L2
  # Target at each column (i.e. sample)
  l1 <- lapply(seq_len(n), function(i) {
    colSums(abs(x[, i, ] - mat[, rep(i, r), drop = FALSE]))
  })
  l2 <- lapply(seq_len(n), function(i) {
    sqrt(colSums((x[, i, ] - mat[, rep(i, r), drop = FALSE])^2))
  })

  r <- data.frame(
    cosine_distance_mean = sapply(d, mean),
    cosine_distance_sd = sapply(d, sd),
    cosine_distance_min = sapply(d, min),
    cosine_distance_max = sapply(d, max),
    L1_mean = sapply(l1, mean),
    L1_sd = sapply(l1, sd),
    L1_min = sapply(l1, min),
    L1_max = sapply(l1, max),
    L2_mean = sapply(l2, mean),
    L2_sd = sapply(l2, sd),
    L2_min = sapply(l2, min),
    L2_max = sapply(l2, max),
    stringsAsFactors = FALSE
  )
  r
}

# Column representing signatures
get_cosine_distance <- function(x, y) {
  # rows for x and cols for y
  out <- 1 - cosineMatrix(x, y)
  rownames(out) <- paste0("S", seq_len(ncol(x)))
  colnames(out) <- paste0("S", seq_len(ncol(y)))
  out
}

# test_mat <- matrix(
#   c(0.6, 0.1, 0.5, 0.3, 0.1, 0.7, 0.2, 0.4, 0.9),
#   nrow = 3
# )
# rownames(test_mat) <- colnames(test_mat) <- c("S1", "S2", "S3")
get_matches <- function(mat, runs) {
  # Find the column index of the second run responds to the first run
  match_index <- apply(mat, 1, which.min) %>% as.integer()
  index_uniq <- unique(match_index)
  if (length(match_index) != length(index_uniq)) {
    # Some signatures in the first run is matched by more than once
    # Order the columns by cum minimum distance
    min_orders <- get_min_orders(mat)
    # Assign the match by the order
    index_uniq <- vector("integer", length = length(match_index))
    for (i in seq_along(min_orders)) {
      values_order <- order(mat[, min_orders[i]])

      for (j in seq_along(values_order)) {
        # possible index
        pi <- which(values_order == j)
        if (index_uniq[pi] == 0L) {
          index_uniq[pi] <- min_orders[i]
          break()
        }
      }
    }
  }
  out <- data.table::data.table(
    v1 = rownames(mat),
    v2 = colnames(mat)[index_uniq]
  )
  colnames(out) <- runs
  out
}

get_min_orders <- function(mat) {
  cummins <- apply(mat, 2, function(x) {
    cumsum(sort(x))
  }) %>% t()
  colnames(cummins) <- paste0("ord", seq_len(ncol(cummins)))
  cummins %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(.i = dplyr::row_number()) %>%
    dplyr::arrange_at(colnames(cummins)) %>%
    dplyr::pull(".i")
}

# My implementation of clustering with match algorithm proposed
# by Nature Cancer paper by following description in supplementary material
# Steps:
# 1. 得到不同 runs 的 signature 结果列表，进行编号
# 2. 一对一配对计算相似性，并得到距离矩阵
# 3. 每个距离矩阵计算最小平均距离
# 4. 按平均距离对配对 run 进行排序，得到排序好的列表
# 5. 初始化排序列表（步骤4的子集），按顺序利用算法逐步合并
clustering_with_match <- function(match_list, n) {
  on.exit(invisible(gc()), add = TRUE)
  # Input: a match list ordered by mean distance
  # This function implements the core step:
  #   collapse the match list into one data.table step by step
  len <- length(match_list)

  if (len >= 2) {
    reduceG <- function(G) {
      # Reduce elements of G if at least two elements
      # contain common column names
      # G >= 2 elements here
      if (length(G) < 2) {
        return(G)
      }
      cnames <- purrr::map(G, colnames)
      check_list <- combn(seq_along(cnames), 2, simplify = FALSE)
      common <- purrr::map(check_list, ~ intersect(cnames[[.[1]]], cnames[[.[2]]]))

      # Index to reduce
      ri <- purrr::map_lgl(common, ~ length(.) != 0)
      if (any(ri)) {
        purrr::map2(check_list[ri], common[ri], .f = function(x, y) {
          if (!is.na(G[x[1]]) & !is.na(G[x[2]])) {
            # Update global G in reduceG
            G[[min(x)]] <<- merge(G[[x[1]]], G[[x[2]]], by = y)
            # to make sure the data is removed and the index
            # is kept to avoid "subscript out of bounds" error
            G[[max(x)]] <<- NA
          }
        })
        # Remove elements set to NA
        G <- G[!is.na(G)]
        return(reduceG(G))
      } else {
        return(G)
      }
    }

    updateG <- function(G, m) {
      # Add a match m to the set G
      # Check if m can be merged into G
      # if not, add m as a new element of G
      # if yes, merge the m and go further
      # check if the G can be reduced
      flag <- FALSE
      for (i in seq_along(G)) {
        nm <- colnames(G[[i]])
        byi <- colnames(m) %in% nm
        if (any(byi)) {
          # Can be merged
          G[[i]] <- merge(G[[i]], m, by = colnames(m)[byi])
          flag <- TRUE
          break()
        }
      }
      if (isFALSE(flag)) {
        # Set as a new member of G
        G[[length(G) + 1]] <- m
      } else if (length(G) > 1) {
        # Go further check if G can be reduced
        # 这是个递归操作
        G <- reduceG(G)
      }

      return(G)
    }

    G <- match_list[1]
    pair_list <- purrr::map(match_list, colnames)

    # Loop for integration
    # 找到另外 n - 2 个需要合并的 match
    # 边查找，边合并
    for (i in seq_along(pair_list)) {
      if (i == 1L) next()
      if (any(purrr::map_lgl(G, ~ all(pair_list[[i]] %in% colnames(.))))) {
        # Do nothing because the result has been included
        next()
      } else {
        # Include the result
        G <- updateG(G, match_list[[i]])
      }
      if (length(G) == 1L & ncol(G[[1]]) == n) break()
    }
    G <- G[[1]]
    # Check the G set
    if (any(ncol(G) != n, nrow(G) != nrow(match_list[[1]]))) {
      stop("Clustering failed! There are some bugs in code the developer not found, please report it!")
    }

    G[, paste0("Run", seq_len(ncol(G))), with = FALSE]
  } else {
    match_list[[1]]
  }
}

rank_solutions <- function(stats) {
  # 同时迭代 signature number, 指标 data.frame 和对应的排序方法
  # stats 必须按 signature_number 从小到大 order，不然一些定量会有问题
  # 如果 rank 得分相同，后面应该选择 signature 数目较大的
  stats <- stats[order(stats$signature_number), ]

  measures <- c(
    "silhouette", "sample_cosine_distance", "L2_error",
    "exposure_positive_correlation", "signature_similarity_within_cluster"
  )
  send_info(
    "Calculating integrated rank score based on the measures: ",
    paste(measures, collapse = ", ")
  )
  weights <- c(0.25, 0.2, 0.2, 0.1, 0.25)
  send_info(
    "Corresponding weights for obtaining the aggregated score are: ",
    paste(weights, collapse = ", ")
  )
  types <- c("diff", "increase", "increase", "increase", "diff")

  rk <- purrr::map2_df(
    .x = stats[, measures],
    .y = types,
    .f = function(x, y) {
      if (y == "increase") {
        # Smaller is better
        100 - (rank(x) - 1) * 5
      } else if (y == "decrease") {
        # Larger is better
        100 - (length(x) - rank(x, na.last = FALSE)) * 5
      } else {
        # diff type
        df <- c(diff(x), 0)
        100 - (rank(df) - 1) * 5
      }
    }
  )

  rk$aggregated_score <- as.numeric(as.matrix(rk) %*% weights)
  rk <- cbind(data.frame(signature_number = stats$signature_number), rk)

  rk
}

# Activity helpers --------------------------------------------------------

# 构建一个 signature by sample 逻辑矩阵，
# 先指示某个 signature 是否在某个样本中存在
# 根据已知的 exposure 进行初始化，设定<0.05 相对贡献不存在
# 0 表示不存在 signature，1 表示存在，2 表示全局 signature
# 一开始的矩阵只可能是一个1和2构成的矩阵，后续的操作会进行填0操作
construct_sig_exist_matrix <- function(expo, sample_class = NULL, cutoff = 0.05) {
  out <- expo
  if (is.null(sample_class)) {
    # 没有群组标签（即1组），那么所有 signature 都有可能
    out[, ] <- 1L
  } else {
    # 只有 1 组标签也是如此
    grps <- unique(sample_class)
    if (length(grps) == 1L) {
      out[, ] <- 1L
    } else {
      out <- ifelse(out > cutoff, 1L, 0L) # 先初始化

      # Check sample_class
      sample_class <- sample_class[colnames(expo)]

      if (length(sample_class) != ncol(expo)) {
        stop("Not all samples are assigned to a class (subtype).")
      }

      # 用 group 标签覆盖样本标签
      lapply(grps, function(grp) {
        idx <- names(sample_class[sample_class == grp])
        ex <- sum(out[, idx] == 1L) > 5
        if (ex) {
          out[, idx] <<- 1L
        }
      })
      out[out == 0L] <- 2L # 还存在 0 的地方标记为全局 signature
    }
  }
  return(out)
}

## Bootstrap approach
optimize_exposure_in_one_sample_bt <- function(catalog,
                                               flag_vector, # a set of 1L or 2L
                                               sample,
                                               sig_matrix,
                                               tmp_dir,
                                               bt_use_prop = FALSE) {
  tmpf_expo <- file.path(
    tmp_dir,
    paste0(sample, "_expo.rds")
  )
  tmpf_sim <- file.path(
    tmp_dir,
    paste0(sample, "_sim.rds")
  )

  if (all(file.exists(tmpf_expo), file.exists(tmpf_sim))) {
    message("ALL result file exists, skipping...")
  } else {
    on.exit(invisible(gc()), add = TRUE)
    force(flag_vector)
    expo2 <- vector("numeric", length(flag_vector))
    flag1 <- which(flag_vector == 1L)

    message("Handling sample: ", sample)

    catalog_mat <- matrix(
      catalog,
      ncol = 1,
      dimnames = list(rownames(sig_matrix), sample)
    )

    out <- sig_fit_bootstrap(
      catalog_mat,
      sig_matrix[, flag1, drop = FALSE],
      n = 100,
      type = "absolute"
    )

    # Use median
    expo <- rowMedians(out$expo, na.rm = TRUE)
    sim <- median(out$cosine, na.rm = TRUE)

    th <- sum(catalog) * 0.01
    reset_idx <- if (bt_use_prop) {
      apply(out$expo, 1, function(x) {
        mean(x < th, na.rm = TRUE) > 0.05
      })
    } else {
      apply(out$expo, 1, function(x) {
        p <- my.t.test.p.value(x, mu = th, alternative = "greater")
        if (is.na(p)) {
          message("NA result detected from t.test, use empirical p value (proportion).")
          p <- mean(x < th, na.rm = TRUE)
        }

        p > 0.05
      })
    }

    if (any(reset_idx)) {
      expo[reset_idx] <- 0
    }

    expo2[flag1] <- expo

    message("Save expo result to temp file ", tmpf_expo)
    saveRDS(expo2, file = tmpf_expo)

    message("Save similarity result to temp file ", tmpf_sim)
    saveRDS(sim, file = tmpf_sim)
  }

  NULL
}


## Stepwise approach
optimize_exposure_in_one_sample <- function(catalog,
                                            flag_vector, # a set of 1L or 2L
                                            sample,
                                            sig_matrix,
                                            tmp_dir) {

  tmpf_expo <- file.path(
    tmp_dir,
    paste0(sample, "_expo.rds")
  )
  tmpf_sim <- file.path(
    tmp_dir,
    paste0(sample, "_sim.rds")
  )

  if (all(file.exists(tmpf_expo), file.exists(tmpf_sim))) {
    message("ALL result file exists, skipping...")
  } else {
    on.exit(invisible(gc()), add = TRUE)
    force(flag_vector)
    flag1_bk <- flag1 <- which(flag_vector == 1L)
    flag2 <- which(flag_vector == 2L)

    expo <- vector("numeric", length(flag_vector))
    catalog_mat <- matrix(
      catalog,
      ncol = 1,
      dimnames = list(rownames(sig_matrix), sample)
    )

    # 先处理得到一个 baseline similarity 值
    message("Processing sample: ", sample)
    message("\t getting baseline similarity.")
    baseline <- .get_one_catalog_similarity(
      catalog_mat, flag1, sig_matrix,
      return_all = TRUE
    )
    message("\t\t", baseline$sim, " based on ", length(flag1), " signatures.")

    # 先处理标记 1，如果相似性降低小于 0.01，移除
    message("\t getting updated similarity by removing one signature in batch.")
    sim_rm <- purrr::map_dbl(
      flag1,
      ~ .get_one_catalog_similarity(catalog_mat, setdiff(flag1, .), sig_matrix)
    )

    rm_id <- sim_rm - baseline$sim >= -0.01 ## 会不会出现全都可以扔掉？待观测

    if (any(rm_id)) {
      flag1 <- flag1[!rm_id]

      if (length(flag1) > 0) {
        baseline <- .get_one_catalog_similarity(
          catalog_mat, flag1, sig_matrix,
          return_all = TRUE
        )
        message("\t\t", baseline$sim, " with ", length(flag1), " signatures left.")
      }
    } else {
      message("\t no signature need to be removed.")
    }

    # 然后再处理可能的标记 2，如果相似性增加大于 0.05，则加入
    if (length(flag2) > 0) {
      message("\t getting updated similarity by adding one global signature in batch.")
      sim_add <- purrr::map_dbl(
        flag2,
        ~ .get_one_catalog_similarity(catalog_mat, c(flag1, .), sig_matrix)
      )

      add_id <- sim_add - baseline$sim > 0.05

      if (any(add_id)) {
        flag1 <- sort(c(flag1, flag2[add_id]))
        baseline <- .get_one_catalog_similarity(
          catalog_mat, flag1, sig_matrix,
          return_all = TRUE
        )
        message("\t\t", baseline$sim, " with ", length(flag1), " signatures left.")
      } else {
        message("\t no global signature need to be added.")
      }
    }

    # Output the result
    if (length(flag1) == 0L) {
      message("\t no signatures left in the whole procedure seems due to low similarity, just report baseline exposure.")
      flag1 <- flag1_bk
    }
    expo[flag1] <- baseline$expo[, 1]

    message("Save expo result to temp file ", tmpf_expo)
    saveRDS(expo, file = tmpf_expo)

    message("Save similarity result to temp file ", tmpf_sim)
    saveRDS(baseline$sim, file = tmpf_sim)
  }

  NULL
}

.get_one_catalog_similarity <- function(catalog_mat, flag, sig_matrix,
                                        return_all = FALSE) {
  expo <- suppressMessages(
    sig_fit(
      catalogue_matrix = catalog_mat,
      sig = sig_matrix[, flag, drop = FALSE],
      method = "QP",
      type = "absolute",
      return_class = "matrix"
    )
  )
  rec_catalog <- sig_matrix[, flag, drop = FALSE] %*% expo
  if (return_all) {
    list(
      sim = cosineMatrix(rec_catalog, catalog_mat)[1],
      expo = expo
    )
  } else {
    cosineMatrix(rec_catalog, catalog_mat)[1]
  }
}
