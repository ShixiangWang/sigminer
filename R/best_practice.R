# A Best Practice for Signature Extraction and Exposure (Activity) Attribution
# 先基于已知的 signature 进行拟合，然后再提取额外的未知signature也是一种可行方式
# 以后还要处理不同平台数据混合导致的差异性！

#' A Best Practice for Signature Extraction and Exposure (Activity) Attribution
#'
#' These functions are combined to provide a best practice for optimally
#' identifying mutational signatures and attributing their activities (exposures)
#' in tumor samples. They are listed in order to use.
#' - `bp_extract_signatures()` for extracting signatures.
#' - `bp_show_survey()` for showing measures change under different
#' signature numbers to help user select optimal signature number.
#' At default, an aggregated score (named score) is generated from 5 measures to
#' suggest the best solution. See section "Measure Explanation in Survey plot"
#' for more explanation.
#' - `bp_get_sig_obj()` for get a (list of) `Signature` object which is common
#' used in **sigminer** for analysis and visualization.
#' - `bp_attribute_activity()` for optimizing signature activities (exposures).
#' - `bp_extract_signatures_iter()` for extracting signature in a iteration way.
#' - `bp_cluster_iter_list()` for clustering iterated signatures to help collapse
#' multiple signatures into one. The result cluster can be visualized by
#' `plot()` or `factoextra::fviz_dend()`.
#' - Extra: `bp_get_stats`() for obtaining stats for signatures and samples of a solution.
#' These stats are aggregated (averaged) as the stats for a solution
#' (specific signature number).
#' - Extra: `bp_get_rank_score()` for obtaining rank score for all signature numbers.
#'
#' @details
#' The signature extraction approach is adopted from reference #1, #2, and
#' the whole best practice is adopted from the pipeline used by reference #3.
#' I implement the whole procedure with R code based on the method description
#' of papers. The code is well organized, tested and documented so user will
#' find it pretty simple and useful. Besides, the structure of the results is
#' very clear to see and also visualize like other approaches provided by **sigminer**.
#'
#' @section Measure Explanation in Survey Plot:
#' The survey plot provides a pretty good way to facilitate the signature number
#' selection. A `score` measure is calculated as the weighted mean of 5 important
#' measures and visualized as the first sub-plot. The optimal number is highlighted
#' with red color dot and the best values for 5 measures to be weighted are also
#' highlighted with orange color dots. The detail of 6 measures shown in plot are
#' explained as below.
#' - `score` - an aggregated score based on rank scores from 5 measures below.
#' The higher, the better. When two signature numbers have the same score,
#' the larger signature number is preferred (this is a rare situation, you
#' have to double check other measures).
#' - `silhouette` - the average silhouette width for signatures, also named as ASW in reference #2.
#' The signature number with silhouette decreases sharply is preferred.
#' - `distance` - the average sample reconstructed cosine distance, the lower value is better.
#' - `error` - the average sample reconstructed error calculated with L2 formula
#' (i.e. L2 error). This lower value is better. This measure represents a
#' similar concept like `distance` above, they are all used to quantify how well
#' sample mutation profiles can be reconstructed from signatures, but `distance`
#' cares the whole mutation profile similarity while `error` here cares value difference.
#' - `pos cor` - the average positive signature exposure correlation coefficient.
#' The lower value is better. This measure is constructed based on my understanding
#' about signatures: mutational signatures are typically treated as independent
#' recurrent patterns, so their activities are less correlated.
#' - `similarity` - the average similarity within in a signature cluster.
#' Like `silhouette`, the point decreases sharply is preferred.
#' In the practice, results from multiple NMF runs are clustered
#' with "clustering with match" algorithm proposed by reference #2. This value
#' indicates if the signature profiles extracted from different NMF runs are similar.
#' @inheritParams sig_estimate
#' @param n_bootstrap number of bootstrapped (resampling) catalogs used.
#' When it is `0`, the original (input) mutation catalog is used for NMF decomposition,
#' this is not recommended, just for testing, user should not set it to `0`.
#' @param n_nmf_run number of NMF runs for each bootstrapped or original catalog.
#' At default, in total n_bootstrap x n_nmf_run (i.e. 1000) NMF runs are used
#' for the task.
#' @param RTOL a threshold proposed by Nature Cancer paper to control how to
#' filter solutions of NMF. Default is `0.1%` (from reference #2),
#' only NMF solutions with KLD (KL deviance) <= `100.1%` minimal KLD are kept.
#' @param min_contribution a component contribution threshold to filer out small
#' contributed components.
#' @param seed a random seed to make reproducible result.
#' @param handle_hyper_mutation default is `TRUE`, handle hyper-mutant samples.
#' @param report_integer_exposure default is `TRUE`, report integer signature
#' exposure by bootstrapping technique.
#' @return It depends on the called function.
#' @name bp
#' @author Shixiang Wang <w_shixiang@163.com>
#' @references
#' Alexandrov, Ludmil B., et al. "Deciphering signatures of mutational processes operative in human cancer." Cell reports 3.1 (2013): 246-259.
#'
#' Degasperi, Andrea, et al. "A practical framework and online tool for mutational signature analyses show intertissue variation and driver dependencies." Nature cancer 1.2 (2020): 249-263.
#'
#' Alexandrov, Ludmil B., et al. “The repertoire of mutational signatures in human cancer.” Nature 578.7793 (2020): 94-101.
#' @seealso See [sig_estimate], [sig_extract], [sig_auto_extract],
#' [sigprofiler_extract] for other approaches.
#' @examples
#' data("simulated_catalogs")
#' \donttest{
#' # Here I reduce the values for n_bootstrap and n_nmf_run
#' # for reducing the run time.
#' # In practice, you should keep default or increase the values
#' # for better estimation.
#' #
#' # The input data here is simulated from 10 mutational signatures
#' e1 <- bp_extract_signatures(
#'   t(simulated_catalogs$set1),
#'   range = 8:12,
#'   n_bootstrap = 5,
#'   n_nmf_run = 10
#' )
#'
#' # See the survey for different signature numbers
#' # The suggested solution is marked as red dot
#' # with highest integrated score.
#' p1 <- bp_show_survey(e1)
#' p1
#' # You can also exclude plotting and highlighting the score
#' p2 <- bp_show_survey(e1, add_score = FALSE)
#' p2
#'
#' # Obtain the suggested solution from extraction result
#' obj_suggested <- bp_get_sig_obj(e1, e1$suggested)
#' obj_suggested
#' # If you think the suggested signature number is not right
#' # Just pick up the solution you want
#' obj_s8 <- bp_get_sig_obj(e1, 8)
#'
#' # After extraction, you can assign the signatures
#' # to reference COSMIC signatures
#' # More see ?get_sig_similarity
#' sim <- get_sig_similarity(obj_suggested)
#' # Visualize the match result
#' if (require(pheatmap)) {
#'   pheatmap::pheatmap(sim$similarity)
#' }
#'
#' # You already got the activities of signatures
#' # in obj_suggested, however, you can still
#' # try to optimize the result by removing a signature
#' # in a sample if it contributes <0.01 similarity
#' # For more controls, see its documentation
#' expo <- bp_attribute_activity(e1, return_class = "data.table")
#' expo$abs_activity
#'
#' # Iterative extraction:
#' # This procedure will rerun extraction step
#' # for those samples with reconstructed catalog similarity
#' # lower than a threshold (default is 0.95)
#' e2 <- bp_extract_signatures_iter(
#'   t(simulated_catalogs$set1),
#'   range = 9:11,
#'   n_bootstrap = 5,
#'   n_nmf_run = 5,
#'   sim_threshold = 0.99
#' )
#' e2
#' # When the procedure run multiple rounds
#' # you can cluster the signatures from different rounds by
#' # the following command
#' # bp_cluster_iter_list(e2)
#'
#' ## Extra utilities
#' rank_score <- bp_get_rank_score(e1)
#' rank_score
#' stats <- bp_get_stats(e2$iter1)
#' # Get the mean reconstructed similarity
#' 1 - stats$stats_sample$cosine_distance_mean
#' }
#' @testexamples
#' expect_is(e1, "ExtractionResult")
#' expect_is(e2, "ExtractionResultList")
#' expect_is(p1, "ggplot")
#' expect_is(p2, "ggplot")
#' expect_is(obj_suggested, "Signature")
#' expect_is(obj_s8, "Signature")
#' expect_is(sim, "list")
#' expect_is(expo, "list")
#' expect_is(rank_score, "data.frame")
#' expect_is(stats, "list")
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
                                  report_integer_exposure = TRUE) {
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
    is.logical(handle_hyper_mutation)
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
  if (!requireNamespace("doFuture", quietly = TRUE)) {
    send_info("{.pkg doFuture} is recommended to install for improving computation.")
  } else {
    doFuture::registerDoFuture()
    suppressWarnings(future::plan("multiprocess", workers = cores))
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
    rm(solution_list); invisible(gc(verbose = FALSE))
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

  solutions$catalog_matrix <- raw_catalogue_matrix
  class(solutions) <- "ExtractionResult"

  send_success("Extraction procedure run successfully.")
  solutions
}

#' @param sim_threshold a similarity threshold for selecting samples to auto-rerun
#' the extraction procedure (i.e. `bp_extract_signatures()`), default is `0.95`.
#' @param max_iter the maximum iteration size, default is 10, i.e., at most run
#' the extraction procedure 10 times.
#' @param cache_dir a directory for storing result for each round. If you
#' don't want to use it, set it to `FALSE`.
#' @rdname bp
#' @export
bp_extract_signatures_iter <- function(nmf_matrix,
                                       range = 2:5,
                                       sim_threshold = 0.95,
                                       max_iter = 10L,
                                       n_bootstrap = 20L,
                                       n_nmf_run = 50,
                                       RTOL = 1e-3, min_contribution = 0,
                                       cores = min(4L, future::availableCores()),
                                       seed = 123456L,
                                       handle_hyper_mutation = TRUE,
                                       report_integer_exposure = TRUE,
                                       cache_dir = getwd()) {
  iter_list <- list()
  for (i in seq_len(max_iter)) {
    message("Round #", i)
    message("===============================")
    bp <- bp_extract_signatures(
      nmf_matrix = nmf_matrix,
      range = range,
      n_bootstrap = n_bootstrap,
      n_nmf_run = n_nmf_run,
      RTOL = RTOL,
      min_contribution = min_contribution,
      cores = cores,
      seed = seed,
      handle_hyper_mutation = handle_hyper_mutation,
      report_integer_exposure = report_integer_exposure
    )
    # 检查寻找需要重新运行的样本，修改 nmf_matrix
    iter_list[[paste0("iter", i)]] <- bp

    if (!isFALSE(cache_dir)) {
      if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
      j <- 1L
      check_file <- file.path(cache_dir, paste0("job", j, "_round_1.rds"))
      while (file.exists(check_file)) {
        message("Cache file ", check_file, " exists, changing name...")
        j <- j + 1L
        check_file <- file.path(cache_dir, paste0("job", j, "_round_1.rds"))
      }
      cache_file <- file.path(cache_dir, paste0("job", j, "_round_", i, ".rds"))
      message("Save round #", i, " result to ", cache_file)
      saveRDS(bp, file = cache_file)
    }

    samp2rerun <- bp$stats_sample %>%
      dplyr::filter(.data$signature_number == bp$suggested) %>%
      dplyr::filter(.data$cosine_distance_mean > 1 - sim_threshold) %>%
      dplyr::pull("sample")

    nsamp <- dim(bp$catalog_matrix)[2]
    if (length(samp2rerun) < 2L) {
      if (length(samp2rerun) == 0) {
        message("All samples passed the rerun threshold in round #", i)
      } else {
        message("Only one sample did not pass the rerun threshold in round #", i, ". Stop here.")
      }
      message("Return.")
      break()
    } else if (length(samp2rerun) == nsamp) {
      message("All samples cannot pass rerun threshold in round #", i)
      if (i == 1L) message("It is the first round, maybe your should lower the value.")
      message("Return.")
      break()
    } else {
      # Rerun
      nmf_matrix <- t(bp$catalog_matrix)
      nmf_matrix <- nmf_matrix[samp2rerun, ]
    }
  }
  class(iter_list) <- "ExtractionResultList"
  iter_list
}

#' @param x result from [bp_extract_signatures_iter()] or a list of
#' `Signature` objects.
#' @param include_final_iteration if `FALSE`, exclude final iteration result
#' from clustering for input from [bp_extract_signatures_iter()], not applied
#' if input is a list of `Signature` objects.
#' @rdname bp
#' @export
bp_cluster_iter_list <- function(x, include_final_iteration = TRUE) {
  if (length(x) < 2) {
    stop("No need to cluster length-1 result list.")
  }
  if (inherits(x, "ExtractionResultList")) {
    x <- purrr::map(x, ~ bp_get_sig_obj(., .$suggested))
    if (isFALSE(include_final_iteration)) {
      x <- x[-length(x)]
    }
  }
  if (!inherits(x[[1]], "Signature")) {
    stop("The list element should be a Signature object.")
  }

  sig_list <- purrr::map(x, "Signature.norm")
  sigmat <- purrr::imap(sig_list, function(x, y) {
    colnames(x) <- paste(y, colnames(x), sep = ":")
    x
  }) %>% purrr::reduce(cbind)
  cosdist <- 1 - cosineMatrix(sigmat, sigmat)
  rownames(cosdist) <- colnames(cosdist) <- colnames(sigmat)
  # Do clustering
  cls <- stats::hclust(stats::as.dist(cosdist))
  r <- list(
    cluster = cls,
    distance = cosdist,
    sigmat = sigmat
  )
  class(r) <- "SignatureListClusters"
  r
}

#' @param obj a `ExtractionResult` object from [bp_extract_signatures()].
#' @param signum a integer vector to extract the corresponding `Signature` object(s).
#' If it is `NULL` (default), all will be returned.
#' @rdname bp
#' @export
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

#' @rdname bp
#' @export
bp_get_stats <- function(obj) {
  assert_class(obj, "ExtractionResult")
  obj[c("stats", "stats_signature", "stats_sample")]
}

#' @rdname bp
#' @export
bp_get_rank_score <- function(obj) {
  assert_class(obj, "ExtractionResult")
  obj[["rank_score"]]
}

#' @rdname bp
#' @param add_score if `FALSE`, don't show score and label optimal points by
#' rank score.
#' @param scales one of "free_y" (default) and "free" to control the scales
#' of plot facet.
#' @param fixed_ratio if `TRUE` (default), make the x/y axis ratio fixed.
#' @export
bp_show_survey <- function(obj,
                           add_score = TRUE,
                           scales = c("free_y", "free"),
                           fixed_ratio = TRUE) {
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
    "signature_similarity_within_cluster",
    "exposure_positive_correlation"
  )

  rs <- obj$rank_score[, cols]
  nsig <- nrow(rs)
  plot_df <- merge(
    obj$stats,
    rs[, c("signature_number", "aggregated_score")],
    by = "signature_number"
  ) %>%
    dplyr::select(cols)
  colnames(plot_df) <- colnames(rs) <- c("sn", "as", "sil", "cos", "err", "sim", "corr")
  cn <- c("score", "silhouette", "distance", "error", "similarity", "pos cor")
  names(cn) <- c("as", "sil", "cos", "err", "sim", "corr")

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
    dplyr::mutate(rk = rank(.data$measure, na.last = FALSE, ties.method = "first")) %>%
    dplyr::select(-"measure")

  plot_df <- dplyr::left_join(plot_df, rs, by = c("sn", "type")) %>%
    dplyr::mutate(
      type = cn[.data$type],
      type = factor(.data$type, levels = cn)
    )

  if (add_score) {
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
      )
  } else {
    p <- ggplot(
      plot_df %>%
        dplyr::filter(.data$type != "score"),
      aes_string(x = "sn", y = "measure")
    ) +
      geom_line() +
      geom_point()
  }

  p <- p +
    facet_wrap(~type, nrow = 2, scales = scales) +
    cowplot::theme_cowplot() +
    labs(x = NULL, y = NULL)
  if (fixed_ratio) {
    p <- p + theme(aspect.ratio = 1)
  }

  p
}

#' @param input result from [bp_extract_signatures()] or a Signature object.
#' @param sample_class a named string vector whose names are sample names
#' and values are class labels (i.e. cancer subtype). If it is `NULL` (the default),
#' treat all samples as one group.
#' @inheritParams sig_fit
#' @rdname bp
#' @export
bp_attribute_activity <- function(input,
                                  sample_class = NULL,
                                  nmf_matrix = NULL,
                                  return_class = c("matrix", "data.table")) {
  # logical: excludes class specific signatures if it contributes <0.01 similarity
  #          while include global signatures if it add >0.05 similarity

  return_class <- match.arg(return_class)
  if (inherits(input, "ExtractionResult")) {
    if (is.null(nmf_matrix)) {
      nmf_matrix <- t(input$catalog_matrix)
    }

    if (is.na(input$suggested)) {
      # Only one Signature object
      input <- input$object[[1]]
    } else {
      input <- bp_get_sig_obj(input, signum = input$suggested)
    }
  }

  if (inherits(input, "Signature")) {
    sig <- input$Signature.norm
    expo <- input$Exposure.norm
  } else {
    # 其他的输入情况，待定
    stop("Invalid input!")
  }

  if (is.null(nmf_matrix)) {
    stop("nmf_matrix cannot be NULL!")
  }

  # catalog matrix 的 component 顺序必须和 signature profile 矩阵保持一致
  # 存在矩阵和 signature profile 矩阵中 signature 顺序必须一致
  # 存在矩阵与 catalog matrix 中的 sample 顺序也必须一致
  sig_order <- colnames(sig)
  samp_order <- colnames(expo)
  catalog_df <- as.data.frame(t(nmf_matrix))
  catalog_df <- catalog_df[rownames(sig), samp_order, drop = FALSE]
  exist_df <- construct_sig_exist_matrix(expo, sample_class) %>%
    as.data.frame()
  exist_df <- exist_df[sig_order, , drop = FALSE]

  # Handle samples one by one (by columns)
  out <- purrr::pmap(
    .l = list(
      catalog = catalog_df,
      flag_vector = exist_df,
      sample = samp_order
    ),
    .f = optimize_exposure_in_one_sample,
    sig_matrix = sig
  )
  out <- purrr::transpose(out)
  expo <- purrr::reduce(out$expo, cbind)
  rel_expo <- apply(expo, 2, function(x) x / sum(x, na.rm = TRUE))
  sim <- purrr::reduce(out$similarity, c)

  if (return_class == "data.table") {
    expo <- .mat2dt(expo)
    rel_expo <- .mat2dt(rel_expo)
  }

  list(
    abs_activity = expo,
    rel_activity = rel_expo,
    similarity = sim
  )
}

.mat2dt <- function(x) {
  x %>%
    t() %>%
    data.table::as.data.table(keep.rownames = "sample")
}
