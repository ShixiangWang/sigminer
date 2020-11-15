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
          .export = c("k", "verbose"),
          .verbose = FALSE
        ) %dopar% {
          if (verbose) {
            message("Extracting ", k, " signatures with seed: ", s)
            print(bt_matrix)
          }
          NMF::nmf(bt_matrix, rank = k, method = "brunet", seed = s, nrun = 1L)
        }
      }
    )
    KLD_list <- sapply(solution_list, NMF::deviance)

    # Filter solutions with RTOL threshold
    solutions[[paste0("K", range[k])]] <- solution_list[KLD_list <= min(KLD_list) * (1 + RTOL)]
  }

  # Collect solutions
  solutions
  # 聚类：使用 cosine 或相关性作为距离指标

  # 生成统计量
  # 重构相似性，cophenetic，轮廓系数，
  # 聚类平均相似距离，RSS, 平均错误，Exposure 相关性

}

# 获取一些指定的信息
bp_get <- function() {

}

bp_attribute_activity <- function() {

}
