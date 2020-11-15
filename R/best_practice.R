# A Best Practice for Signature Extraction and Activity (Exposure) Attribution
# 先基于已知的 signature 进行拟合，然后再提取额外的未知signature也是一种可行方式
# 以后还要处理不同平台数据混合导致的差异性！

# 当两个随机分布相同时，它们的相对熵为零，当两个随机分布的差别增
# 大时，它们的相对熵也会增大。
# deviance(res) # KLD
bp_extract_signatures <- function(nmf_matrix,
                                  range = 2:10, nrun = 100L,
                                  RTOL = 1e-3, min_contribution = 0,
                                  cores = 1L, seed = 123456L,
                                  handle_hyper_mutation = TRUE,
                                  mpi_platform = FALSE,
                                  verbose = FALSE) {
  stopifnot(is.matrix(nmf_matrix),
            !is.null(rownames(nmf_matrix)), !is.null(colnames(nmf_matrix)),
            is.numeric(range), is.numeric(nrun), # nrun >= 100,
            is.numeric(RTOL), RTOL > 0,
            is.numeric(min_contribution),
            min_contribution >= 0, min_contribution <= 0.1,
            is.numeric(seed),
            is.logical(handle_hyper_mutation), is.logical(mpi_platform))
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
      send_info("Dropping the components with very few contribution: ",
                paste(names(contris[contris_index]), collapse = ", "))
      nmf_matrix <- nmf_matrix[, !contris_index, drop = FALSE]
    }
  }

  # 超突变处理
  if (handle_hyper_mutation) {
    nmf_matrix <- handle_hyper_mutation(nmf_matrix)
  }

  seeds <- seq(seed, length = nrun)
  # Generate bootstrapped catalogs based on nrun
  bt_catalog_list <- purrr::map(seeds, function(x) {
    set.seed(x)
    simulate_catalogue_matrix(nmf_matrix)
  })

  # NMF with brunet method
  if (isFALSE(mpi_platform)) {
    if (!requireNamespace("doFuture", quietly = TRUE)) {
      send_info("{.pkg doFuture} is recommended to install for improving computation.")
    } else {
      doFuture::registerDoFuture()
      future::plan("multiprocess", workers = cores)
    }
  }

  solutions <- list()
  for (k in seq_along(range)) {
    send_info("Extracting ", range[k], " signatures on bootstrapped catalogs...")
    solution_list <- suppressWarnings(
      {
        foreach(
          s = seeds,
          bt_matrix = bt_catalog_list,
          .packages = "NMF",
          .export = c("k", "range", "verbose"),
          .verbose = FALSE
        ) %dopar% {
          if (verbose) {
            message("Extracting ", range[k], " signatures with seed: ", s)
            print(bt_matrix)
          }
          NMF::nmf(bt_matrix, rank = range[k], method = "brunet", seed = s, nrun = 1L)
        }
      }
    )
    KLD_list <- sapply(solution_list, NMF::deviance)

    # Filter solutions with RTOL threshold
    solutions[[paste0("K", range[k])]] <- solution_list[KLD_list <= min(KLD_list) * (1 + RTOL)]
  }

  # Collect solutions
  # 先将所有 solution 标准化处理，得到 signature 和 activity
  # 然后针对 signature 使用 clustering with match 算法进行聚类
  solutions <- purrr::map(solutions, .f = function(solution_list) {
    out <- purrr::map(solution_list, .f = normalize_solution) %>%
      setNames(paste0("Run", seq_along(solution_list)))
    # To do: Do clustering with match
    run_pairs <- combn(names(out), 2, simplify = FALSE)
    # Get similarity distance
    out
  })
  solutions
  # 聚类：使用 cosine 或相关性作为距离指标

  # 生成统计量
  # 重构相似性，cophenetic，轮廓系数，
  # 聚类平均相似距离，RSS, 平均错误，Exposure 相关性

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

# My implementation of clustering with match algorithm proposed
# by Nature Cancer paper by following description in supplementary material
# Steps:
# 1. 得到不同 runs 的 signature 结果列表，进行编号
# 2. 一对一配对计算相似性，并得到距离矩阵
# 3. 每个距离矩阵计算最小平均距离
# 4. 按平均距离对配对 run 进行排序，得到排序好的列表
# 5. 初始化排序列表（步骤4的子集），按顺序利用算法逐步合并（左连接）
clustering_with_match <- function() {

}

# 获取一些指定的信息
bp_get <- function() {

}

bp_attribute_activity <- function() {

}
