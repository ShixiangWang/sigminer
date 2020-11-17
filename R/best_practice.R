# A Best Practice for Signature Extraction and Activity (Exposure) Attribution
# 先基于已知的 signature 进行拟合，然后再提取额外的未知signature也是一种可行方式
# 以后还要处理不同平台数据混合导致的差异性！

# 当两个随机分布相同时，它们的相对熵为零，当两个随机分布的差别增
# 大时，它们的相对熵也会增大。
# deviance(res) # KLD
bp_extract_signatures <- function(nmf_matrix,
                                  range = 2:10,
                                  n_bootstrap = 20L,
                                  n_nmf_run = 50,
                                  RTOL = 1e-3, min_contribution = 0,
                                  cores = 1L, seed = 123456L,
                                  handle_hyper_mutation = TRUE,
                                  mpi_platform = FALSE,
                                  verbose = FALSE) {
  stopifnot(
    is.matrix(nmf_matrix),
    !is.null(rownames(nmf_matrix)), !is.null(colnames(nmf_matrix)),
    is.numeric(range), is.numeric(n_bootstrap), is.numeric(n_nmf_run),
    is.numeric(RTOL), RTOL > 0,
    is.numeric(min_contribution),
    min_contribution >= 0, min_contribution <= 0.1,
    is.numeric(seed),
    is.logical(handle_hyper_mutation), is.logical(mpi_platform)
  )
  seed <- as.integer(seed)

  timer <- Sys.time()
  send_info("Best practice for signature extraction started.")
  send_info("NOTE: the input should be a sample-by-component matrix.")
  on.exit(send_elapsed_time(timer))

  # Input: a matrix used for NMF decomposition with rows indicate samples and columns indicate components.
  input <- nmf_matrix

  # Dimension reduction: 去掉总贡献小的 components
  if (min_contribution != 0) {
    contris <- colSums(nmf_matrix) / sum(nmf_matrix)
    contris_index <- contris <= min_contribution
    if (any(contris_index)) {
      send_info(
        "Dropping the components with very few contribution: ",
        paste(names(contris[contris_index]), collapse = ", ")
      )
      nmf_matrix <- nmf_matrix[, !contris_index, drop = FALSE]
    }
  }

  # 超突变处理
  if (handle_hyper_mutation) {
    nmf_matrix <- handle_hyper_mutation(nmf_matrix)
  }

  seeds_bt <- seq(seed, length = n_bootstrap)
  # Generate bootstrapped catalogs based on n_bootstrap
  bt_catalog_list <- purrr::map(seeds_bt, function(x) {
    set.seed(x)
    simulate_catalogue_matrix(nmf_matrix)
  })
  nmf_matrix <- t(nmf_matrix)

  # NMF with brunet method
  if (isFALSE(mpi_platform)) {
    if (!requireNamespace("doFuture", quietly = TRUE)) {
      send_info("{.pkg doFuture} is recommended to install for improving computation.")
    } else {
      doFuture::registerDoFuture()
      suppressWarnings(future::plan("multiprocess", workers = cores))
    }
  }

  seeds <- seq(seed, length = n_bootstrap * n_nmf_run)
  solutions <- list()
  for (k in seq_along(range)) {
    send_info(
      "Extracting ", range[k], " signatures on ",
      n_bootstrap, " bootstrapped catalogs with ",
      n_nmf_run, " NMF runs for each."
    )
    solution_list <- suppressWarnings({
      foreach(
        s = seeds,
        bt_matrix = bt_catalog_list[rep(seq_len(n_bootstrap), each = n_nmf_run)],
        .packages = "NMF",
        .export = c("k", "range", "verbose"),
        .verbose = FALSE
      ) %dopar% {
        if (verbose) {
          message("Extracting ", range[k], " signatures with seed: ", s)
          print(t(bt_matrix))
        }
        NMF::nmf(t(bt_matrix), rank = range[k], method = "brunet", seed = s, nrun = 1L)
      }
    })

    # Filter solutions with RTOL threshold in each bootstrap chunk
    solutions[[paste0("K", range[k])]] <- purrr::map(
      chunk2(solution_list, n_bootstrap),
      .f = function(s) {
        KLD_list <- sapply(s, NMF::deviance)
        s <- s[KLD_list <= min(KLD_list) * (1 + RTOL)]
        if (length(s) > 10) s <- s[1:10] # Limits 10 best runs
        s
      }
    ) %>% unlist()
  }

  # Collect solutions
  # 先将所有 solution 标准化处理，得到 signature 和 activity
  # 然后针对 signature 使用 clustering with match 算法进行聚类
  # 聚类：使用 1 - cosine 相似性作为距离指标
  solutions <- purrr::map(solutions,
    .f = process_solution,
    nmf_matrix = nmf_matrix
  )
  solutions
}

process_solution <- function(slist, nmf_matrix) {
  out <- purrr::map(slist, .f = normalize_solution) %>%
    setNames(paste0("Run", seq_along(slist)))

  # If just one run, skip the following steps
  if (length(out) >= 2) {
    run_pairs <- combn(names(out), 2, simplify = FALSE)
    run_pairs_name <- purrr::map_chr(run_pairs, ~ paste(., collapse = "_"))
    # Get similarity distance
    pair_dist <- purrr::map(
      run_pairs,
      ~ get_cosine_distance(out[[.[1]]]$Signature, out[[.[2]]]$Signature)
    ) %>%
      setNames(run_pairs_name)
    pair_dist_mean <- purrr::map_dbl(pair_dist, mean) %>%
      setNames(run_pairs_name)
    match_list <- purrr::map2(pair_dist, run_pairs, get_matches)

    # Order the result by mean distance
    res_orders <- order(pair_dist_mean)
    run_pairs <- run_pairs[res_orders]
    pair_dist <- pair_dist[res_orders]
    pair_dist_mean <- pair_dist_mean[res_orders]
    match_list <- match_list[res_orders]

    # Do clustering with match
    clusters <- clustering_with_match(match_list, n = length(out))

    # 按 match cluster 排序
    # 样本排序也确保对齐
    samp_order <- colnames(out$Run1$Exposure)
    out <- purrr::map2(clusters, out, function(x, y, samp_order) {
      y$Signature <- y$Signature[, x, drop = FALSE]
      y$Exposure <- y$Exposure[x, samp_order, drop = FALSE]
      y
    }, samp_order = samp_order)
    Signature <- purrr::reduce(purrr::map(out, "Signature"), `+`) / length(out)
    Exposure <- purrr::reduce(purrr::map(out, "Exposure"), `+`) / length(out)
    KLD <- purrr::reduce(purrr::map(out, "KLD"), `+`) / length(out)
  } else {
    run_pairs <- NA
    pair_dist <- NA
    pair_dist_mean <- NA
    match_list <- NA
    clusters <- data.table::data.table(
      Run1 = paste0("S", seq_len(ncol(out$Run1$Signature)))
    )

    Signature <- out$Run1$Signature
    Exposure <- out$Run1$Exposure
    KLD <- out$Run1$KLD
  }

  # Order by contribution and rename signatures
  new_order <- order(colSums(Signature), decreasing = TRUE)
  Signature <- Signature[, new_order, drop = FALSE]
  Exposure <- Exposure[new_order, , drop = FALSE]
  colnames(Signature) <- rownames(Exposure) <- paste0("Sig", seq_along(new_order))

  # 生成统计量
  # Remember the new order above
  # 分为 signature 和 样本两种，取每个度量的最大、最小值、平均值以及 SD
  # 重构相似性，cophenetic，轮廓系数，
  # 聚类平均相似距离，RSS, 平均错误，Exposure 相关性等
  stat_sigs <- get_stat_sigs(out)
  stat_samps <- get_stat_samps(out, mat = nmf_matrix)

  list(
    runs = out,
    run_pairs = run_pairs,
    pair_dist = pair_dist,
    pair_dist_mean = pair_dist_mean,
    match_list = match_list,
    clusters = clusters,
    stat_sigs = stat_sigs,
    stat_samps = stat_samps
  )
}

get_3d_array_stats <- function(x, ns = NULL) {
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
  if (type == "within-cluster") {
    d <- lapply(seq_len(n), function(i) {
      mat <- cosineMatrix(x[, i, ], x[, i, ])
      mat[upper.tri(mat)]
    })
  } else if (type == "between-cluster") {
    d <- lapply(seq_len(n), function(i) {
      x1 <- x[, i, ]
      x2 <- x[, -i, ]
      dim(x2) <- c(dim(x1)[1], prod(dim(x2)) / dim(x1)[1])
      cosineMatrix(x1, x2)
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
  n <- dim(x)[1] # n signatures
  r <- dim(x)[3] # r runs

  get_cor <- function(x1, x2) {
    nc <- ncol(x1)
    sapply(seq_len(nc), function(i) {
      stats::cor(x1[, i], x2[, i])
    })
  }

  d <- lapply(seq_len(n), function(i) {
    x1 <- x[i, , ]
    x2 <- x[-i, , ]
    if (length(dim(x2)) < 3) {
      get_cor(x1, x2)
    } else {
      apply(x2, 1, function(m) {
        get_cor(x1, m)
      }) %>% rowMeans()
    }
  })

  r <- list(
    expo_cor_mean = sapply(d, mean),
    expo_cor_sd = sapply(d, sd),
    expo_cor_min = sapply(d, min),
    expo_cor_max = sapply(d, max)
  )
  r
}

# The difference between reconstructed catalogs and the original catalog
# 用相似性、L1、L2范数
get_error_stats <- function(x, mat) {
  n <- dim(mat)[2] # n samples
  r <- dim(x)[3] # n runs

  # similarity distance
  d <- lapply(seq_len(n), function(i) {
    1 - cosineMatrix(x[, i, ], mat[, i, drop = FALSE])
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

get_stat_sigs <- function(runs) {
  sig_list <- purrr::map(runs, "Signature")
  dm <- dim(sig_list[[1]])
  l <- length(sig_list)
  sig_array <- array(unlist(sig_list), dim = c(dm, l))
  # signatures
  sig <- data.frame(
    signature = seq_len(dm[2]) # Rename it outside the function
  )
  s <- get_3d_array_stats(
    sig_array,
    c(
      "signature_mean", "signature_sd",
      "signature_min", "signature_max"
    )
  )
  sim <- get_similarity_stats(
    sig_array,
    n = dm[2],
    ns = c(
      "similarity_mean", "similarity_sd",
      "similarity_min", "similarity_max"
    )
  )

  # cluster silhouette
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
  dm <- dim(expo_list[[2]])
  expo_array <- array(unlist(expo_list), dim = c(dm, l))
  expo_cor <- get_expo_corr_stat(expo_array)

  list(
    signature = s,
    stats = cbind(sig, KLD, sil_width, expo_cor, sim, cross_sim)
  )
}

get_stat_samps <- function(runs, mat) {
  expo_list <- purrr::map(runs, "Exposure")
  dm <- dim(expo_list[[2]]) # the second value of dm indicates how many samples
  l <- length(expo_list)
  expo_array <- array(unlist(expo_list), dim = c(dm, l))
  # exposures
  e <- get_3d_array_stats(
    expo_array,
    c(
      "exposure_mean", "exposure_sd",
      "exposure_min", "exposure_max"
    )
  )
  # get catalog array
  sig_list <- purrr::map(runs, "Signature")
  catalog_list <- purrr::map2(sig_list, expo_list, ~ .x %*% .y)
  dm2 <- dim(catalog_list[[1]])
  catalog_array <- array(unlist(catalog_list), dim = c(dm2, l))

  sim <- get_similarity_stats(
    catalog_array,
    n = dm2[2],
    ns = c(
      "similarity_mean", "similarity_sd",
      "similarity_min", "similarity_max"
    )
  )

  # cluster silhouette
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

  samp <- data.frame(
    sample = colnames(mat)
  )

  # error
  error <- get_error_stats(
    catalog_array,
    mat
  )

  list(
    exposure = e,
    stats = cbind(samp, sil_width, error)
  )
}

normalize_solution <- function(solution) {
  # solution is a NMF fit result

  W <- NMF::basis(solution)
  H <- NMF::coef(solution)
  K <- ncol(W)
  KLD <- NMF::deviance(solution)

  out <- c(helper_scale_nmf_matrix(W, H, K, handle_cn = FALSE), KLD = KLD)
  colnames(out$Signature) <- rownames(out$Exposure) <- paste0("S", seq_len(K))
  out
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

# 获取一些指定的信息
bp_get <- function() {

}

bp_attribute_activity <- function() {

}
