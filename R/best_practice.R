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
      future::plan("multiprocess", workers = cores)
    }
  }

  seeds <- seq(seed, length = n_bootstrap * n_nmf_run)
  solutions <- list()
  for (k in seq_along(range)) {
    send_info("Extracting ", range[k], " signatures on ",
              n_bootstrap, " bootstrapped catalogs with ",
              n_nmf_run, " NMF runs for each.")
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
      })
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

get_stat_sigs <- function(runs) {
  sig_list <- purrr::map(runs, "Signature")
  dm <- dim(sig_list[[1]])
  l <- length(sig_list)
  sig_array <- array(unlist(sig_list), dim = c(dm, l))
  # signatures
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

  list(
    signature = s,
    stats = cbind(KLD, sim, sil_width, cross_sim)
  )
}

get_stat_samps <- function(runs, mat) {
  expo_list <- purrr::map(runs, "Exposure")
  dm <- dim(expo_list[[2]]) # the second value indicates how many samples
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
  catalog_list <- purrr::map2(sig_list, expo_list, ~.x%*%.y)
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

  list(
    exposure = e,
    stats = cbind(sim, sil_width, cross_sim)
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
# 5. 初始化排序列表（步骤4的子集），按顺序利用算法逐步合并（左连接）
clustering_with_match <- function(match_list, n) {
  # Input: a match list ordered by mean distance
  # This function implements the core step:
  #   collapse the match list into one step by step

  # Init
  len <- length(match_list)

  if (len >= 2) {
    # L = match_list[seq(2L, len)]
    # G = match_list[1]
    pair_list <- purrr::map(match_list, colnames)

    G <- vector("character", 2 * (n - 1))
    G[1:2] <- pair_list[[1]]
    G_index <- vector("integer", n - 1)
    G_index[1] <- 1L

    # Loop for integration
    # 找到另外 n - 2 个需要合并的 match 的索引
    j <- 1L
    for (i in seq_along(pair_list)) {
      if (i == 1L) next()
      if (all(pair_list[[i]] %in% G)) {
        # Do nothing because the result has been included
        next()
      } else {
        # Include the result
        G[seq(2 * j + 1, 2 * j + 2)] <- pair_list[[i]]
        G_index[j + 1] <- i
        j <- j + 1L
      }
      if (G_index[n - 1] != 0) break()
    }

    integrated_list <- match_list[G_index]
    # Make sure the list can be merged correctly
    # res <- tryCatch(
    #   purrr::reduce(integrated_list, merge),
    #   error = function(e) {
    #     to_join <- integrated_list[[1]]
    #     be_join <- integrated_list[-1]
    #     while (length(be_join) > 0) {
    #       col_exist <- colnames(be_join[[1]]) %in% colnames(to_join)
    #       if (any(col_exist)) {
    #         to_join <- merge(to_join, be_join[[1]], by = colnames(be_join[[1]])[col_exist])
    #         be_join[[1]] <- NULL
    #       } else {
    #         be_join <- shifter(be_join)
    #       }
    #     }
    #   }
    # )
    res <- purrr::reduce(
      purrr::map(integrated_list, as.data.frame),
      merge
    ) %>% data.table::as.data.table()

    res[, paste0("Run", seq_len(ncol(res))), with = FALSE]
  } else {
    match_list[[1]]
  }
}

# 获取一些指定的信息
bp_get <- function() {

}

bp_attribute_activity <- function() {

}
