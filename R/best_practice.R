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
#' At default, an aggregated score (named score) is generated to
#' suggest the best solution.
#' - `bp_show_survey2()` for showing simplified signature number survey like
#' [show_sig_number_survey()].
#' - `bp_get_sig_obj()` for get a (list of) `Signature` object which is common
#' used in **sigminer** for analysis and visualization.
#' - `bp_attribute_activity()` for optimizing signature activities (exposures).
#' NOTE: the activities from extraction step may be better!
#' You can also use [sig_extract] to get optimal NMF result from multiple NMF runs.
#' Besides, you can use [sig_fit] to quantify exposures based on signatures extracted
#' from `bp_extract_signatures()`.
#' - `bp_extract_signatures_iter()` for extracting signature in a iteration way.
#' - `bp_cluster_iter_list()` for clustering (`hclust` with average linkage)
#' iterated signatures to help collapse
#' multiple signatures into one. The result cluster can be visualized by
#' `plot()` or `factoextra::fviz_dend()`.
#' - `bp_get_clustered_sigs()` for getting clustered (grouped) mean signatures from signature clusters.
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
#' selection. A `score` measure is calculated as the weighted mean of selected
#' measures and visualized as the first sub-plot. The optimal number is highlighted
#' with red color dot and the best values for each measures are also
#' highlighted with orange color dots. The detail of 6 measures shown in plot are
#' explained as below.
#' - `score` - an aggregated score based on rank scores from selected measures below.
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
#' @param cores_solution cores for processing solutions, default is equal to argument `cores`.
#' @param seed a random seed to make reproducible result.
#' @param handle_hyper_mutation default is `TRUE`, handle hyper-mutant samples.
#' @param report_integer_exposure if `TRUE`, report integer signature
#' exposure by bootstrapping technique.
#' @param only_core_stats if `TRUE`, only calculate the core stats for signatures and samples.
#' @param cache_dir a directory for storing intermediate result, also avoid
#' repeated computation for same data.
#' @param keep_cache default is `FALSE`, if `TRUE`, keep cache data.
#' For small input data, it is not necessary.
#' @param pynmf if `TRUE`, use Python NMF driver [Nimfa](http://nimfa.biolab.si/index.html).
#' The seed currently is not used by this implementation, so the only way to reproduce
#' your result is setting `keep_cache = TRUE`.
#' @param use_conda if `TRUE`, create an independent conda environment to run NMF.
#' @param py_path path to Python executable file, e.g. '/Users/wsx/anaconda3/bin/python'. In my
#' test, it is more stable than `use_conda=TRUE`. You can install the Nimfa package by yourself
#' or set `use_conda` to `TRUE` to install required Python environment, and then set this option.
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
#' # You can also plot a simplified version
#' p3 <- bp_show_survey2(e1, highlight = 10)
#' p3
#'
#' # Obtain the suggested solution from extraction result
#' obj_suggested <- bp_get_sig_obj(e1, e1$suggested)
#' obj_suggested
#' # If you think the suggested signature number is not right
#' # Just pick up the solution you want
#' obj_s8 <- bp_get_sig_obj(e1, 8)
#'
#' # Track the reconstructed profile similarity
#' rec_sim <- get_sig_rec_similarity(obj_s8, t(simulated_catalogs$set1))
#' rec_sim
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
#' # try to optimize the result.
#' # NOTE: the optimization step may not truly optimize the result!
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
#' expect_is(p2, "ggplot")
#' expect_is(obj_suggested, "Signature")
#' expect_is(obj_s8, "Signature")
#' expect_is(sim, "list")
#' expect_is(rec_sim, "data.table")
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
                                  cores_solution = min(cores, length(range)),
                                  seed = 123456L,
                                  handle_hyper_mutation = TRUE,
                                  report_integer_exposure = FALSE,
                                  only_core_stats = nrow(nmf_matrix) > 100,
                                  cache_dir = file.path(tempdir(), "sigminer_bp"),
                                  keep_cache = FALSE,
                                  pynmf = FALSE,
                                  use_conda = TRUE,
                                  py_path = "/Users/wsx/anaconda3/bin/python") {
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
    is.logical(handle_hyper_mutation),
    is.numeric(cores), is.numeric(cores_solution)
  )
  seed <- as.integer(seed)
  cores <- as.integer(cores)
  cores_solution <- as.integer(cores_solution)
  range <- sort(unique(range))

  ii <- rowSums(nmf_matrix) < 0.01
  if (any(ii)) {
    message(
      "The follow samples dropped due to null catalogue:\n\t",
      paste0(rownames(nmf_matrix)[ii], collapse = ", ")
    )
    nmf_matrix <- nmf_matrix[!ii, , drop = FALSE]
  }

  if (!requireNamespace("synchronicity")) {
    install.packages("synchronicity")
  }
  if (!requireNamespace("lpSolve")) {
    install.packages("lpSolve")
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
  if (pynmf) {
    bt_catalog_list <- purrr::map(bt_catalog_list, t)
  } else {
    bt_catalog_list <- purrr::map(bt_catalog_list, ~ t(check_nmf_matrix(.)))
  }
  catalogue_matrix <- t(nmf_matrix)
  send_success("Done.")

  # Construct cache file names
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
  if (pynmf) {
    cache_files <- file.path(
      cache_dir,
      paste0(
        digest::digest(list(catalogue_matrix, seed)),
        "_pyNMF_K",
        range,
        ".rds"
      )
    )
  } else {
    cache_files <- file.path(
      cache_dir,
      paste0(
        digest::digest(list(catalogue_matrix, seed)),
        "_NMF_K",
        rep(range, each = length(bt_catalog_list)),
        if (bt_flag) "_bt_catalog_" else "_observed_catalog_",
        seq_len(n_bootstrap),
        ".rds"
      )
    )
  }

  extract_nmf <- function(solution) {
    W <- NMF::basis(solution)
    H <- NMF::coef(solution)
    K <- ncol(W)
    KLD <- NMF::deviance(solution)
    return(list(
      W = W,
      H = H,
      K = K,
      KLD = KLD
    ))
  }

  filter_nmf <- function(s, bt_flag, RTOL) {
    KLD_list <- purrr::map_dbl(s, "KLD")
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

  solutions <- list()
  cache_list <- chunk2(cache_files, n = length(range))
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

    fl <- cache_list[[k]]

    if (pynmf) {
      # Use Python nimfa package
      env_install(use_conda, py_path, pkg = "nimfa", pkg_version = "1.4.0")
      send_info("Calling Python as backend to run NMF.")
      if (!file.exists(fl[1])) {
        slist <- call_pynmf(bt_catalog_list, range[k], n_nmf_run, cores)
        # Filter solutions with RTOL threshold
        if (bt_flag) {
          send_info("Keeping at most 10 NMF solutions with KLD within (1+RTOL) range of the best.")
        } else {
          send_info("Keeping at most 100 best NMF solutions.")
        }
        slist <- purrr::map(chunk2(slist, length(bt_catalog_list)),
          filter_nmf,
          bt_flag = bt_flag, RTOL = RTOL
        ) %>%
          purrr::flatten()

        solutions[[paste0("K", range[k])]] <- slist
        if (keep_cache) {
          send_info("Saving cache file to ", fl[1])
          saveRDS(slist, file = fl[1])
        }

        rm(slist)
        invisible(gc())
      } else {
        send_info("Cache run file ", fl[1], " already exists, skip calling.")
        send_info("Reading cache file...")
        solutions[[paste0("K", range[k])]] <- readRDS(fl[1])
      }
    } else {
      # Use R NMF package
      eval(parse(text = "suppressMessages(library('NMF'))"))

      call_nmf <- function(i, cores = 1L) {
        if (!file.exists(fl[i])) {
          send_info(sprintf(
            "Running signature number %s for %s catalog %d.",
            range[k],
            if (bt_flag) "bootstrapped" else "observed",
            i
          ))

          r <- NMF::nmf(
            bt_catalog_list[[i]],
            rank = range[k],
            method = "brunet",
            seed = seed,
            nrun = n_nmf_run,
            .options = paste0("v", if (ncol(bt_catalog_list[[i]]) > 100) 4 else 1, "mkp", cores)
          )
          if (inherits(r, "NMFfit")) {
            r <- list(r)
          } else {
            r <- r@.Data
          }
          r <- purrr::map(r, extract_nmf)

          # Filter solutions with RTOL threshold
          if (bt_flag) {
            send_info("Keeping at most 10 NMF solutions with KLD within (1+RTOL) range of the best.")
          } else {
            send_info("Keeping at most 100 best NMF solutions.")
          }
          r <- filter_nmf(r, bt_flag, RTOL)
          saveRDS(r, file = fl[i])
          rm(r)
          invisible(gc())
        } else {
          send_info("Cache run file ", fl[i], " already exists, skip.")
        }
      }

      if (n_bootstrap >= 10 * n_nmf_run) {
        oplan <- future::plan()
        future::plan(set_future_strategy(), workers = cores, gc = TRUE)
        on.exit(future::plan(oplan), add = TRUE, after = FALSE)
        furrr::future_map(seq_len(n_bootstrap), call_nmf,
          cores = 1L,
          .progress = TRUE,
          .options = furrr::furrr_options(
            seed = TRUE,
            stdout = FALSE
          )
        )
      } else {
        purrr::map(seq_len(n_bootstrap), call_nmf, cores = cores)
      }

      send_info("Reading NMF run files...")
      solutions[[paste0("K", range[k])]] <- purrr::map(cache_list[[k]], readRDS) %>%
        purrr::flatten()
      send_success("Read successfully.")
      if (!keep_cache) {
        signal <- file.remove(cache_list[[k]])
        if (all(signal)) {
          send_success("NMF run files deleted.")
        } else {
          send_warning("Delete some or all NMF run files failed.")
        }
      }
    }

    send_success(
      "Done for this signature number. Current memory size used: ",
      round(mem_used() / 2^20), "MB"
    )
    invisible(gc(verbose = FALSE))
  }

  send_info("Starting process the solution list.")
  # Collect solutions
  # 先将所有 solution 标准化处理，得到 signature 和 activity
  # 然后针对 signature 使用 clustering with match 算法进行聚类
  # 聚类：使用 1 - cosine 相似性作为距离指标
  if ((sum(sapply(solutions, length)) < 200L & length(solutions) < 4) | cores_solution == 1) {
    solutions <- purrr::map(
      solutions,
      .f = process_solution,
      catalogue_matrix = catalogue_matrix,
      report_integer_exposure = report_integer_exposure,
      only_core_stats = only_core_stats
    )
  } else {
    cores <- cores_solution
    cores <- min(cores, length(solutions))
    send_info(cores, " cores set for processing solutions.")
    oplan <- future::plan()
    future::plan(set_future_strategy(), workers = cores, gc = TRUE, .skip = TRUE)
    on.exit(future::plan(oplan), add = TRUE, after = FALSE)
    solutions <- furrr::future_map(
      solutions,
      .f = process_solution,
      catalogue_matrix = catalogue_matrix,
      report_integer_exposure = report_integer_exposure,
      only_core_stats = only_core_stats,
      .progress = TRUE,
      .options = furrr::furrr_options(
        seed = seed,
        stdout = FALSE
      )
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
                                       cores_solution = min(cores, length(range)),
                                       seed = 123456L,
                                       handle_hyper_mutation = TRUE,
                                       report_integer_exposure = FALSE,
                                       only_core_stats = nrow(nmf_matrix) > 100,
                                       cache_dir = file.path(tempdir(), "sigminer_bp"),
                                       keep_cache = FALSE,
                                       pynmf = FALSE,
                                       use_conda = FALSE,
                                       py_path = "/Users/wsx/anaconda3/bin/python") {
  force(only_core_stats)
  iter_list <- list()
  cache_file_list <- c()
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
      cores_solution = cores_solution,
      seed = seed,
      handle_hyper_mutation = handle_hyper_mutation,
      report_integer_exposure = report_integer_exposure,
      only_core_stats = only_core_stats,
      cache_dir = cache_dir,
      keep_cache = keep_cache,
      pynmf = pynmf,
      use_conda = use_conda,
      py_path = py_path
    )
    # 检查寻找需要重新运行的样本，修改 nmf_matrix
    iter_list[[paste0("iter", i)]] <- bp

    if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
    cache_file <- file.path(
      cache_dir,
      paste0(digest::digest(nmf_matrix), "_round_", i, ".rds")
    )
    message("Save round #", i, " result to ", cache_file)
    saveRDS(bp, file = cache_file)
    cache_file_list <- c(cache_file_list, cache_file)

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

  if (!keep_cache) {
    on.exit(unlink(cache_file_list, force = TRUE), add = TRUE)
  }

  message("Done.")
  class(iter_list) <- "ExtractionResultList"
  iter_list
}

#' @param x result from [bp_extract_signatures_iter()] or a list of
#' `Signature` objects.
#' @param k an integer sequence specifying the cluster number to get silhouette.
#' @param include_final_iteration if `FALSE`, exclude final iteration result
#' from clustering for input from [bp_extract_signatures_iter()], not applied
#' if input is a list of `Signature` objects.
#' @rdname bp
#' @export
bp_cluster_iter_list <- function(x, k = NULL, include_final_iteration = TRUE) {
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
  distobj <- stats::as.dist(cosdist)
  cls <- stats::hclust(distobj, method = "average")
  if (is.null(k)) {
    # Set default k sequence
    k <- seq(2, max(sapply(sig_list, ncol)))
  }
  sil_df <- purrr::map_df(
    k,
    function(x) {
      cls <- cluster::silhouette(cutree(cls, k = x), distobj)
      cbind(
        data.frame(
          signame = rownames(cosdist),
          k = x
        ),
        data.frame(cls[, 1:3])
      )
    }
  )
  sil_summary <- sil_df %>%
    dplyr::group_by(.data$k) %>%
    dplyr::summarise(
      min = min(.data$sil_width, na.rm = TRUE),
      mean = mean(.data$sil_width, na.rm = TRUE),
      max = max(.data$sil_width, na.rm = TRUE),
      sd = sd(.data$sil_width, na.rm = TRUE)
    ) %>%
    as.data.frame()
  r <- list(
    cluster = cls,
    distance = cosdist,
    sil_df = sil_df,
    sil_summary = sil_summary,
    sigmat = sigmat
  )
  class(r) <- "SignatureListClusters"
  r
}

bp_get_cluster_index_list <- function(x) {
  rg <- range(x)
  sq <- seq(rg[1], rg[2])
  y <- purrr::map(sq, ~ which(x == .))
  names(y) <- as.character(sq)
  y
}


#' @param SigClusters result from [bp_cluster_iter_list()].
#' @param cluster_label cluster labels for a specified cluster number, obtain it
#' from `SigClusters$sil_df`.
#' @rdname bp
#' @export
bp_get_clustered_sigs <- function(SigClusters, cluster_label) {
  sig_idx <- bp_get_cluster_index_list(cluster_label)
  sig_map <- purrr::map(sig_idx, ~ colnames(SigClusters$sigmat)[.])
  names(sig_map) <- paste0("Sig", names(sig_map))
  grp_sigs <- purrr::reduce(
    purrr::map(sig_idx, ~ t(t(rowMeans(SigClusters$sigmat[, ., drop = FALSE])))),
    cbind
  )
  colnames(grp_sigs) <- names(sig_map)

  sim_sig_to_grp_mean <- purrr::map(
    sig_idx,
    ~ as.numeric(
      cosine(
        SigClusters$sigmat[, ., drop = FALSE],
        t(t(rowMeans(SigClusters$sigmat[, ., drop = FALSE])))
      )
    )
  )
  sim_sig_to_grp_mean <- purrr::map2_df(
    sim_sig_to_grp_mean,
    sig_idx,
    ~ data.frame(sig = colnames(SigClusters$sigmat)[.y], sim = .x),
    .id = "grp_sig"
  )
  sim_sig_to_grp_mean$grp_sig <- paste0("Sig", sim_sig_to_grp_mean$grp_sig)

  return(
    list(
      grp_sigs = grp_sigs,
      sim_sig_to_grp_mean = sim_sig_to_grp_mean,
      sig_map = sig_map
    )
  )
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
#' @inheritParams show_sig_number_survey
#' @export
bp_show_survey2 <- function(obj,
                            x = "signature_number",
                            left_y = "silhouette", right_y = "L2_error",
                            left_name = left_y, right_name = right_y,
                            left_color = "black", right_color = "red",
                            left_shape = 16, right_shape = 18,
                            shape_size = 4, highlight = NULL) {
  stopifnot(class(obj) == "ExtractionResult" | is.data.frame(obj))
  if (class(obj) == "ExtractionResult") {
    obj <- obj$stats
  }

  message("Variables can be used: ", paste(colnames(obj), collapse = ", "))

  show_sig_number_survey(
    object = obj,
    x = x,
    left_y = left_y, right_y = right_y,
    left_name = left_name, right_name = right_name,
    left_color = left_color, right_color = right_color,
    left_shape = left_shape, right_shape = right_shape,
    shape_size = shape_size, highlight = highlight
  )
}

#' @rdname bp
#' @param add_score if `FALSE`, don't show score and label optimal points by
#' rank score.
#' @param scales one of "free_y" (default) and "free" to control the scales
#' of plot facet.
#' @param fixed_ratio if `TRUE` (default), make the x/y axis ratio fixed.
#' @export
bp_show_survey <- function(obj,
                           add_score = FALSE,
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

  if (!add_score) {
    plot_df <- plot_df %>%
      dplyr::filter(!type %in% c("silhouette", "error"))
  }

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
    facet_wrap(~type, nrow = if (length(unique(plot_df$type)) > 3) 2 else 1, scales = scales) +
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
#' @param method one of 'bt' (use bootstrap exposure median, from reference #2,
#' **the most recommended way in my personal view**) or stepwise'
#' (stepwise reduce and update signatures then do signature fitting
#' with last signature sets, from reference #2, the result tends to assign
#' the contribution of removed signatures to the remaining signatures,
#' **maybe I misunderstand the paper method? PAY ATTENTION**).
#' @param bt_use_prop this parameter is only used for `bt` method to reset
#' low contributing signature activity (relative activity `<0.01`). If `TRUE`,
#' use empirical P value calculation way (i.e. proportion, used by reference `#2`),
#' otherwise a `t.test` is applied.
#' @param cache_dir a directory for keep temp result files.
#' @param keep_cache if `TRUE`, keep cache results.
#' @inheritParams sig_fit
#' @inheritParams sig_fit_bootstrap_batch
#' @rdname bp
#' @export
bp_attribute_activity <- function(input,
                                  sample_class = NULL,
                                  nmf_matrix = NULL,
                                  method = c("bt", "stepwise"),
                                  bt_use_prop = FALSE,
                                  return_class = c("matrix", "data.table"),
                                  use_parallel = FALSE,
                                  cache_dir = file.path(tempdir(), "sigminer_attribute_activity"),
                                  keep_cache = FALSE) {
  # logical for bt: get median exposure value from signature fitting,
  #          we set exposures of a signature to zero if
  #          (1) bootstrapped exposure distribution is significantly less
  #              than the threshold of 1% of total mutations with t.test.
  #              If a t.test is impractical, then (2) is applied.
  #              In common sense, program will never run (2).
  #          (2) more than 5% of
  #          the bootstrapped exposure values (empirical P = 0.05) are
  #          below the threshold of 1% of total mutations in the sample.
  # logical for stepwise: excludes class specific signatures if it contributes <0.01 similarity
  #          while include global signatures if it add >0.05 similarity

  set.seed(123456, kind = "L'Ecuyer-CMRG")

  method <- match.arg(method)
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

  if (!is.null(sample_class)) {
    samps <- names(sample_class)
    if (any(duplicated(samps))) {
      stop("A sample cannot be assigned to two groups!")
    }
    idx <- colnames(expo) %in% samps
    if (length(samps) != ncol(expo) | !all(idx)) {
      n_total <- ncol(expo)
      samps2 <- colnames(expo)[idx]
      if (n_total - length(samps2) > 0) {
        warning(
          n_total - length(samps2),
          " samples cannot be found in 'sample_class', they will be removed!",
          immediate. = TRUE
        )
        expo <- expo[, samps2, drop = FALSE]
        nmf_matrix <- nmf_matrix[samps2, ]
      }
      sample_class <- sample_class[samps2]
    }
  }

  # catalog matrix 的 component 顺序必须和 signature profile 矩阵保持一致
  # 存在矩阵和 signature profile 矩阵中 signature 顺序必须一致
  # 存在矩阵与 catalog matrix 中的 sample 顺序也必须一致
  sig_order <- colnames(sig)
  samp_order <- colnames(expo)
  catalog_df <- as.data.frame(t(nmf_matrix))
  catalog_samp <- colnames(catalog_df)
  ix <- samp_order %in% catalog_samp
  if (!all(ix)) {
    message("The following samples have no original catalogs, they will be removed.")
    message("\t", paste(samp_order[!ix], collapse = ", "))
    samp_order <- samp_order[ix]
  }
  catalog_df <- catalog_df[rownames(sig), samp_order, drop = FALSE]
  exist_df <- as.data.frame(construct_sig_exist_matrix(expo, sample_class))
  exist_df <- exist_df[sig_order, samp_order, drop = FALSE]

  act_tmp_dir <- cache_dir
  if (!dir.exists(act_tmp_dir)) {
    dir.create(act_tmp_dir, recursive = TRUE)
  }
  message("Set cache dir: ", act_tmp_dir)

  if (use_parallel) {
    oplan <- future::plan()
    future::plan(set_future_strategy(),
      workers = if (is.logical(use_parallel)) {
        future::availableCores()
      } else {
        use_parallel
      }
    )
    on.exit(future::plan(oplan), add = TRUE)
  }

  if (method == "bt") {
    if (use_parallel) {
      suppressMessages(
        furrr::future_pmap(
          .l = list(
            catalog = catalog_df,
            flag_vector = exist_df,
            sample = samp_order
          ),
          .f = optimize_exposure_in_one_sample_bt,
          sig_matrix = sig,
          tmp_dir = act_tmp_dir,
          bt_use_prop = bt_use_prop,
          .progress = TRUE,
          .options = furrr::furrr_options(
            seed = TRUE,
            stdout = FALSE
          )
        )
      )
    } else {
      purrr::pmap(
        .l = list(
          catalog = catalog_df,
          flag_vector = exist_df,
          sample = samp_order
        ),
        .f = optimize_exposure_in_one_sample_bt,
        sig_matrix = sig,
        tmp_dir = act_tmp_dir,
        bt_use_prop = bt_use_prop
      )
    }
  } else {
    # Handle samples one by one (by columns)
    if (use_parallel) {
      suppressMessages(
        furrr::future_pmap(
          .l = list(
            catalog = catalog_df,
            flag_vector = exist_df,
            sample = samp_order
          ),
          .f = optimize_exposure_in_one_sample,
          sig_matrix = sig,
          tmp_dir = act_tmp_dir,
          .progress = TRUE,
          .options = furrr::furrr_options(
            seed = TRUE,
            stdout = FALSE
          )
        )
      )
    } else {
      purrr::pmap(
        .l = list(
          catalog = catalog_df,
          flag_vector = exist_df,
          sample = samp_order
        ),
        .f = optimize_exposure_in_one_sample,
        sig_matrix = sig,
        tmp_dir = act_tmp_dir
      )
    }
  }

  message()
  message("Reading result files...")
  expo_files <- file.path(
    act_tmp_dir,
    paste0(samp_order, "_expo.rds")
  )
  sim_files <- file.path(
    act_tmp_dir,
    paste0(samp_order, "_sim.rds")
  )
  expo <- purrr::map(expo_files, readRDS)
  sim <- purrr::map_dbl(sim_files, readRDS)
  expo <- purrr::reduce(expo, c)
  expo <- matrix(expo, ncol = length(samp_order))
  rownames(expo) <- sig_order
  colnames(expo) <- samp_order

  rel_expo <- apply(expo, 2, function(x) x / sum(x, na.rm = TRUE))

  if (!keep_cache) {
    message("Clean cache files when exits.")
    on.exit(unlink(act_tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)
  }

  if (return_class == "data.table") {
    expo <- .mat2dt(expo)
    rel_expo <- .mat2dt(rel_expo)
  }

  message("Done.")

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


# Python NMF ------------------------------------------------------------------

call_pynmf <- function(V_list, rank, nrun = 1L, cores = 1L) {
  on.exit(invisible(gc()))
  rank <- as.integer(rank)
  nrun <- as.integer(nrun)
  cores <- as.integer(cores)
  reticulate::source_python(system.file("py", "nmf.py", package = "sigminer", mustWork = TRUE))
  res <- MultiVNMF(V_list, rank, nrun, cores) %>%
    purrr::flatten() %>%
    purrr::flatten()
  res <- purrr::map(res, function(x) {
    rownames(x$W) <- rownames(V_list[[1]])
    colnames(x$H) <- colnames(V_list[[1]])
    x
  })
  res
}

# Install python environment
env_install <- function(use_conda, py_path, pkg, pkg_version) {
  if (!requireNamespace("reticulate")) {
    stop("Package 'reticulate' is required, please install it firstly!")
  }

  if (use_conda) {
    tryCatch(reticulate::conda_binary(),
      error = function(e) {
        message("Cannot find conda binary, installing miniconda...")
        reticulate::install_miniconda()
      }
    )

    ## Prepare conda environment and packages
    tryCatch(
      reticulate::use_condaenv("sigminer_sigprofiler", required = TRUE),
      error = function(e) {
        if (grepl("failed to initialize requested version of Python", e$message)) {
          stop(e$message, ". Please restart R and try again.")
        }
        message("Conda environment not detected, creat it and install required packages.")
        message("======")
        reticulate::conda_create("sigminer_sigprofiler")
        message("Installing packages, be patient...")
        message("======")
        reticulate::conda_install("sigminer_sigprofiler",
          packages = paste0(pkg, "==", pkg_version),
          pip = TRUE
        )
        reticulate::use_condaenv("sigminer_sigprofiler", required = TRUE)
      }
    )

    message("Python and conda environment configuration.")
    message("====================")
    print(reticulate::py_config())
  } else {
    if (is.null(py_path)) {
      if (!reticulate::py_available(initialize = TRUE)) {
        stop("Python is not available!")
      }
      config <- reticulate::py_config()
      reticulate::use_python(config$python, required = TRUE)
    } else {
      reticulate::use_python(py_path, required = TRUE)
    }

    message("Python environment configuration.")
    message("====================")
    print(reticulate::py_config())
  }

  if (pkg == "SigProfilerExtractor" & !reticulate::py_module_available("torch")) {
    message("torch not found, try installing it...")
    tryCatch(
      reticulate::py_install("torch==1.5.1", pip = TRUE),
      error = function(e) {
        message("Cannot install torch just with pip, try again with source from https://download.pytorch.org/whl/torch_stable.html")
        reticulate::py_install("torch==1.5.1", pip = TRUE, pip_options = "-f https://download.pytorch.org/whl/torch_stable.html")
      }
    )
  }
  if (!reticulate::py_module_available(pkg)) {
    message("Python module ", pkg, " not found, try installing it...")
    reticulate::py_install(paste0(pkg, "==", pkg_version), pip = TRUE)
  }

  if (pkg != "SigProfilerExtractor") {
    if (!reticulate::py_module_available("multiprocess")) {
      message("Python module multiprocess not found, try installing it...")
      reticulate::py_install("multiprocess", pip = TRUE)
    }
  }
}
