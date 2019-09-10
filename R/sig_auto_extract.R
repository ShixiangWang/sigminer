sig_auto_extract <- function(nmf_matrix,
                             result_prefix = "BayesNMF",
                             destdir = tempdir(),
                             K0 = 25,
                             nrun = 10,
                             n_iter = 2e5,
                             tol = 1e-07,
                             cores = 1) {

  on.exit(invisible(gc()))  # clean when exit
  nmf_matrix = t(nmf_matrix)  # rows for mutation types and columns for samples

  filelist = file.path(destdir, paste(result_prefix, seq_len(nrun), "rds", sep = "."))
  future::plan("multiprocess", workers = cores)
  furrr::future_map(seq_len(nrun), function(i) {
    # Apply BayesNMF - L1W.L2H for an exponential prior for W and a half-normal prior for H
    res <- BayesNMF.L1W.L2H(nmf_matrix, n_iter, 10, 5, tol, K0, K0, 1)
    saveRDS(res, file = filelist[i])
  }, .progress = TRUE)

  summary.run = purrr::map_df(seq_len(nrun), function(i) {
    res = readRDS(filelist[i])
    K = sum(colSums(res[[1]]) > 1)
    posterior = -res[[4]]  # res[[4]] = -log(posterior)
    dplyr::tibble(
     Run = i,
     K = K,
     posterior = posterior,
     file = filelist[i]
    )
  })

  summary.run = summary.run %>%
    dplyr::arrange(dplyr::desc(.data$posterior))

  # select best solution
  best = names(sort(table(summary.run$K), decreasing = TRUE))[1] %>%
    as.integer()

  best_row = dplyr::filter(summary.run, .data$K == best) %>%
    head(1)

  message("Select Run ", best_row$Run, ", which K = ", best_row$K, " as best solution.")
  best_solution = get_bayesian_result(best_row)

  res <- list(
    Signature = best_solution$Signature,
    Signature.norm = best_solution$Signature.norm,
    Exposure = best_solution$Exposure,
    Exposure.norm = best_solution$Exposure.norm,
    K = best_solution$K,
    Raw = list(
      summary_run = summary.run,
      W = best_solution$W,
      H = best_solution$H,
      best_run = best_row$Run
    )
  )
  class(res) <- "Signature"
  attr(res, "nrun") <- nrun
  #attr(res, "seed") <- seed
  attr(res, "call_method") <- "BayesianNMF"

  res
}
