# A Best Practice for Signature Extraction and Exposure (Activity) Attribution
# 先基于已知的 signature 进行拟合，然后再提取额外的未知signature也是一种可行方式
# 以后还要处理不同平台数据混合导致的差异性！

#' A Best Practice for Signature Extraction and Exposure (Activity) Attribution
#'
#' - `bp_extract_signatures()` for extracting signatures.
#'
#' @inheritParams sig_estimate
#' @param n_bootstrap number of bootstrapped (resampling) catalogs used.
#' When it is `0`, the original (input) mutation catalog is used for NMF decomposition.
#' @param n_nmf_run number of NMF runs for each bootstrapped or original catalog.
#' At default, in total n_bootstrap x n_nmf_run (i.e. 1000) NMF runs are used
#' for the task.
#' @param RTOL a threshold proposed by Nature Cancer paper to control how to
#' filter solutions of NMF. Default is `0.1%`, only NMF solutions with KLD
#' (KL deviance) <= `100.1%` minimal KLD are kept.
#' @param min_contribution a component contribution threshold to filer out small
#' contributed components.
#' @param seed a random seed to make reproducible result.
#' @param handle_hyper_mutation default is `TRUE`, handle hyper-mutant samples.
#' @param report_integer_exposure default is `TRUE`, report integer signature
#' exposure by bootstrapping technique.
#' @param mpi_platform keep for implementation.
#'
#' @return It depends on the called function.
#' @name bp
NULL

#' @rdname bp
#' @export
bp_extract_signatures <- function(nmf_matrix,
                                  range = 2:5,
                                  n_bootstrap = 20L,
                                  n_nmf_run = 50,
                                  RTOL = 1e-3, min_contribution = 0,
                                  cores = min(4L, future::availableCores()),
                                  seed = 123456L,
                                  handle_hyper_mutation = TRUE,
                                  report_integer_exposure = TRUE,
                                  mpi_platform = FALSE) {
  stopifnot(
    is.matrix(nmf_matrix),
    !is.null(rownames(nmf_matrix)), !is.null(colnames(nmf_matrix)),
    is.numeric(range),
    is.numeric(n_bootstrap), is.numeric(n_nmf_run),
    n_bootstrap >= 0, n_nmf_run > 0,
    is.numeric(RTOL), RTOL > 0,
    is.numeric(min_contribution),
    min_contribution >= 0, min_contribution <= 0.1,
    is.numeric(seed),
    is.logical(handle_hyper_mutation), is.logical(mpi_platform)
  )
  seed <- as.integer(seed)
  range <- sort(unique(range))

  ii <- rowSums(nmf_matrix) < 0.01
  if (any(ii)) {
    message(
      "The follow samples dropped due to null catalogue:\n\t",
      paste0(rownames(nmf_matrix)[ii], collapse = ", ")
    )
    nmf_matrix <- nmf_matrix[!ii, , drop = FALSE]
  }

  timer <- Sys.time()
  send_info("Best practice for signature extraction started.")
  send_info("NOTE: the input should be a sample-by-component matrix.")
  on.exit(send_elapsed_time(timer))
  on.exit(invisible(gc()), add = TRUE)

  # Input: a matrix used for NMF decomposition with rows indicate samples and columns indicate components.
  raw_catalogue_matrix <- t(nmf_matrix)
  send_success("Input matrix kept.")

  # Dimension reduction: 去掉总贡献小的 components
  if (min_contribution != 0) {
    send_info("Checking contribution of components.")
    contris <- colSums(nmf_matrix) / sum(nmf_matrix)
    contris_index <- contris <= min_contribution
    if (any(contris_index)) {
      send_info(
        "Dropping the components with very few contribution: ",
        paste(names(contris[contris_index]), collapse = ", ")
      )
      nmf_matrix <- nmf_matrix[, !contris_index, drop = FALSE]
    }
    if (ncol(nmf_matrix) < 3) {
      send_stop("Too few components (<3) left!")
    }
    send_success("Checked.")
  }

  # 超突变处理
  if (handle_hyper_mutation) {
    send_info("Spliting samples if it is hyper-mutant.")
    nmf_matrix <- handle_hyper_mutation(nmf_matrix)
    send_success("Done.")
  }

  send_info("Generating data for inputing NMF.")
  if (n_bootstrap == 0) {
    send_success("Resampling is disabled.")
    bt_catalog_list <- list(nmf_matrix)
    n_bootstrap <- 1L
    bt_flag <- FALSE
  } else {
    send_success("Starting resampling (get bootstrapped catalog matrix).")
    seeds_bt <- seq(seed, length = n_bootstrap)
    # Generate bootstrapped catalogs based on n_bootstrap
    bt_catalog_list <- purrr::map(seeds_bt, function(x) {
      set.seed(x)
      simulate_catalogue_matrix(nmf_matrix)
    })
    bt_flag <- TRUE
  }

  # 有必要的话添加一个极小的数，解决 NMF 包可能存在的异常报错问题
  # 一个 component 之和不能为 0，还有其他一些可能引发的异常
  bt_catalog_list <- purrr::map(bt_catalog_list, check_nmf_matrix)

  catalogue_matrix <- t(nmf_matrix)
  send_success("Done.")

  send_info("Running NMF with brunet method (Lee-KLD).")
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
  send_success("Seeds generated for reproducible research.")

  solutions <- list()
  for (k in seq_along(range)) {
    if (bt_flag) {
      send_info(
        "Extracting ", range[k], " signatures on ",
        n_bootstrap, " bootstrapped catalogs with ",
        n_nmf_run, " NMF runs for each."
      )
    } else {
      send_info(
        "Extracting ", range[k], " signatures on ",
        "the original catalog with ",
        n_nmf_run, " NMF runs in a parallel chunk."
      )
    }

    solution_list <- suppressWarnings({
      foreach(
        s = seeds,
        bt_matrix = bt_catalog_list[rep(seq_len(n_bootstrap), each = n_nmf_run)],
        .packages = "NMF",
        .export = c("k", "range"),
        .verbose = FALSE
      ) %dopar% {
        NMF::nmf(t(bt_matrix), rank = range[k], method = "brunet", seed = s, nrun = 1L)
      }
    })

    # Filter solutions with RTOL threshold
    if (bt_flag) {
      send_info("Keeping at most 10 NMF solutions with KLD within (1+RTOL) range of the best.")
    } else {
      send_info("Keeping at most 100 best NMF solutions.")
    }
    solutions[[paste0("K", range[k])]] <- purrr::map(
      chunk2(solution_list, n_bootstrap),
      .f = function(s) {
        KLD_list <- sapply(s, NMF::deviance)
        if (bt_flag) {
          ki <- KLD_list <= min(KLD_list) * (1 + RTOL)
          s <- s[ki]
          if (length(s) > 10) {
            # Limits 10 best runs
            KLD_list <- KLD_list[ki]
            s <- s[order(KLD_list)[1:10]]
          }
        } else if (length(s) > 100 & !bt_flag) {
          s <- s[order(KLD_list)[1:100]]
        }
        s
      }
    ) %>% unlist()
    send_success(
      "NMF done for this solution. Current memory size used: ",
      round(mem_used() / 2^20), "MB"
    )
  }

  send_info("Starting process the solution list.")
  # Collect solutions
  # 先将所有 solution 标准化处理，得到 signature 和 activity
  # 然后针对 signature 使用 clustering with match 算法进行聚类
  # 聚类：使用 1 - cosine 相似性作为距离指标
  if (sum(sapply(solutions, length)) < 200L & length(solutions) < 4) {
    solutions <- purrr::map(
      solutions,
      .f = process_solution,
      catalogue_matrix = catalogue_matrix,
      report_integer_exposure = report_integer_exposure
    )
  } else {
    oplan <- future::plan()
    future::plan("multiprocess", workers = cores)
    on.exit(future::plan(oplan), add = TRUE)
    solutions <- furrr::future_map(
      solutions,
      .f = process_solution,
      catalogue_matrix = catalogue_matrix,
      report_integer_exposure = report_integer_exposure,
      .progress = TRUE,
      .options = furrr::furrr_options(seed = TRUE)
    )
  }
  send_success(
    "Solution list processed. Current memory size used: ",
    round(mem_used() / 2^20), "MB"
  )

  send_info("Merging and checking the solution data.")
  # 合并 solutions
  solutions <- purrr::transpose(solutions)
  solutions$stats <- purrr::reduce(solutions$stats, rbind)
  solutions$stats_signature <- purrr::reduce(solutions$stats_signature, rbind)
  solutions$stats_sample <- purrr::reduce(solutions$stats_sample, rbind)
  # 追加属性
  solutions$object <- purrr::map(solutions$object, .f = function(obj) {
    attr(obj, "nrun") <- n_bootstrap * n_nmf_run
    attr(obj, "seed") <- seed
    obj
  })
  # 如果发现缺少 components，利用 0 进行回补
  kept_comps <- rownames(solutions$signature[[1]]$signature_mean)
  raw_comps <- rownames(raw_catalogue_matrix)
  if (length(kept_comps) < length(raw_comps)) {
    to_add <- setdiff(raw_comps, kept_comps)
    solutions$object <- purrr::map(solutions$object, .f = function(obj, to_add) {
      n <- ncol(obj$Signature)
      m <- length(to_add)

      mat_add <- matrix(rep(0, n * m),
        nrow = m,
        dimnames = list(to_add)
      )
      obj$Signature <- rbind(obj$Signature, mat_add)
      obj$Signature.norm <- rbind(obj$Signature.norm, mat_add)
      obj
    }, to_add = to_add)
  }
  send_success(
    "Merged. Current memory size used: ",
    round(mem_used() / 2^20), "MB"
  )

  if (nrow(solutions$stats) > 1) {
    send_info("Ranking signature extraction solutions.")
    rank_score <- rank_solutions(solutions$stats)
    suggested <- rank_score[order(
      rank_score$aggregated_score,
      rank_score$signature_number,
      decreasing = TRUE
    ), ][1, ]
    send_success(
      "The solution of ",
      suggested$signature_number,
      " signatures (aggregated score = ",
      suggested$aggregated_score,
      ") is the suggested result."
    )
    solutions$rank_score <- rank_score
    solutions$suggested <- as.integer(suggested$signature_number)
  } else {
    solutions$rank_score <- NA
    solutions$suggested <- NA_integer_
  }
  class(solutions) <- "ExtractionResult"

  send_success("Extraction procedure run successfully.")
  solutions
}

# 获取一些指定的信息
bp_get_sig_obj <- function(obj, signum = NULL) {
  assert_class(obj, "ExtractionResult")
  if (is.null(signum)) {
    message("When signum is NULL, output all signature objects as a list.")
    obj$object
  } else {
    signum <- paste0("K", signum)
    if (length(signum) > 1) {
      obj$object[signum]
    } else {
      obj$object[[signum]]
    }
  }
}

bp_get_stats <- function(obj) {
  assert_class(obj, "ExtractionResult")
  obj[c("stats", "stats_signature", "stats_sample")]
}

bp_get_rank_score <- function(obj) {
  assert_class(obj, "ExtractionResult")
  obj[["rank_score"]]
}

bp_show_survey <- function(obj, scales = c("free_y", "free"), fixed_ratio = TRUE) {
  assert_class(obj, "ExtractionResult")
  scales <- match.arg(scales)

  if (!is.data.frame(obj$rank_score)) {
    message("Show survey cannot work when only one signature extracted.")
    return(invisible())
  }

  cols <- c(
    "signature_number",
    "aggregated_score",
    "silhouette",
    "sample_cosine_distance",
    "L2_error",
    "exposure_positive_correlation",
    "signature_similarity_within_cluster"
  )

  rs <- obj$rank_score[, cols]
  nsig <- nrow(rs)
  plot_df <- merge(
    obj$stats,
    rs[, c("signature_number", "aggregated_score")],
    by = "signature_number"
  ) %>%
    dplyr::select(cols)
  colnames(plot_df) <- colnames(rs) <- c("sn", "as", "sil", "cos", "err", "corr", "sim")
  cn <- c("score", "silhouette", "distance", "error", "pos cor", "similarity")
  names(cn) <- c("as", "sil", "cos", "err", "corr", "sim")

  plot_df <- plot_df %>%
    tidyr::pivot_longer(
      cols = -"sn",
      names_to = "type",
      values_to = "measure"
    )
  rs <- rs %>%
    tidyr::pivot_longer(
      cols = -"sn",
      names_to = "type",
      values_to = "measure"
    )
  rs <- rs %>%
    dplyr::group_by(.data$type) %>%
    dplyr::mutate(rk = rank(.data$measure, na.last = FALSE)) %>%
    dplyr::select(-"measure")

  plot_df <- dplyr::left_join(plot_df, rs, by = c("sn", "type")) %>%
    dplyr::mutate(
      type = cn[.data$type],
      type = factor(.data$type, levels = cn)
    )

  p <- ggplot(plot_df, aes_string(x = "sn", y = "measure")) +
    geom_line() +
    geom_point() +
    geom_point(
      data = dplyr::filter(
        plot_df,
        .data$rk == nsig & .data$type != "score"
      ),
      color = "orange"
    ) +
    geom_point(
      data = dplyr::filter(
        plot_df,
        .data$rk == nsig & .data$type == "score"
      ),
      color = "red"
    ) +
    facet_wrap(~type, nrow = 2, scales = scales) +
    cowplot::theme_cowplot() +
    labs(x = NULL, y = NULL)

  if (fixed_ratio) {
    p <- p + theme(aspect.ratio = 1)
  }

  p
}

bp_attribute_activity <- function() {

}
