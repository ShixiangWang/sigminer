#' Exposure Instability Analysis of Signature Exposures with Bootstrapping
#'
#' Read [sig_fit_bootstrap] for more option setting.
#'
#' @inheritParams sig_fit
#' @inheritParams sig_fit_bootstrap
#' @param methods a subset of `c("NNLS", "QP", "SA")`.
#' @param min_count minimal exposure in a sample, default is 1. Any patient has total exposure less
#' than this value will be filtered out.
#' @param p_val_thresholds a vector of relative exposure threshold for calculating p values.
#' @param use_parallel if `TRUE`, use parallel computation based on **furrr** package.
#' @param seed random seed to reproduce the result.
#' @param job_id a job ID, default is `NULL`, can be a string. When not `NULL`, all bootstrapped results
#' will be saved to local machine location defined by `result_dir`. This is very useful for running
#' more than 10 times for more than 100 samples.
#' @param result_dir see above, default is temp directory defined by R.
#' @param ... other common parameters passing to [sig_fit_bootstrap], including
#'  `sig`, `sig_index`, `sig_db`, `db_type`, `mode`, `auto_reduce` etc.
#'
#' @return a `list` of `data.table`.
#' @export
#' @seealso [sig_fit], [sig_fit_bootstrap]
#'
#' @examples
#' W <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
#' colnames(W) <- c("sig1", "sig2")
#' W <- apply(W, 2, function(x) x / sum(x))
#'
#' H <- matrix(c(2, 5, 3, 6, 1, 9, 1, 2), ncol = 4)
#' colnames(H) <- paste0("samp", 1:4)
#'
#' V <- W %*% H
#' V
#'
#' if (requireNamespace("quadprog")) {
#'   z10 <- sig_fit_bootstrap_batch(V, sig = W, n = 10)
#'   z10
#' }
#' @testexamples
#' expect_is(z10, "list")
#' z2 <- sig_fit_bootstrap_batch(V, sig = W, n = 2, use_parallel = TRUE)
#' expect_is(z2, "list")
sig_fit_bootstrap_batch <- function(catalogue_matrix,
                                    methods = c("QP"),
                                    n = 100L,
                                    min_count = 1L,
                                    p_val_thresholds = c(0.05),
                                    use_parallel = FALSE,
                                    seed = 123456L,
                                    job_id = NULL,
                                    result_dir = tempdir(),
                                    ...) {
  stopifnot(is.matrix(catalogue_matrix))

  set.seed(seed, kind = "L'Ecuyer-CMRG")
  methods <- match.arg(methods, choices = c("NNLS", "QP", "SA"), several.ok = TRUE)

  timer <- Sys.time()
  send_info("Batch Bootstrap Signature Exposure Analysis Started.")
  on.exit(send_elapsed_time(timer, "Total "))

  if (!is.null(job_id)) {
    if (!dir.exists(result_dir)) {
      dir.create(result_dir, recursive = TRUE)
    }
    send_success("Job mode is enabled. All bootstrapped results will be saved to ", result_dir)
  }

  samp_index <- colSums(catalogue_matrix) > min_count
  if (sum(samp_index) != nrow(catalogue_matrix)) {
    send_info(paste(
      "Samples to be filtered out:",
      paste(colnames(catalogue_matrix)[!samp_index], collapse = ",")
    ))
    catalogue_matrix <- catalogue_matrix[, samp_index, drop = FALSE]
  }

  if (ncol(catalogue_matrix) < 1) {
    send_stop("No sample left after filtering, please check your input!")
  }

  ## Get optimal exposures with different methods
  send_info("Finding optimal exposures (&errors) for different methods.")
  optimal_list <- list()
  for (m in methods) {
    send_info("Calling method {.code ", m, "}.")
    expo_list <- sig_fit(catalogue_matrix,
      method = m,
      return_error = TRUE,
      return_class = "data.table",
      ...
    )
    expo_list$expo <- data.table::melt(expo_list$expo, id.vars = "sample", variable.name = "sig", value.name = "exposure")
    expo_list$errors <- data.table::data.table(
      sample = names(expo_list$errors),
      errors = as.numeric(expo_list$errors)
    )
    expo_list$cosine <- data.table::data.table(
      sample = names(expo_list$cosine),
      cosine = as.numeric(expo_list$cosine)
    )
    optimal_list[[m]] <- expo_list
  }
  optimal_list <- purrr::transpose(optimal_list)
  optimal_list <- purrr::map(optimal_list, data.table::rbindlist, fill = TRUE, idcol = "method")

  ## Get bootstrap exposures with different methods
  send_info("Getting bootstrap exposures (&errors/similarity) for different methods.")
  send_info("This step is time consuming, please be patient.")
  call_bt <- function(x, sample, y, methods, n = 1000, ...) {
    if (!is.null(job_id)) {
      fpath <- file.path(result_dir, paste0(job_id, "_", sample, ".rds"))
      if (file.exists(fpath)) {
        send_info("Sample '", sample, "' has been processed.")
        return(fpath)
      }
    }
    on.exit(gc(verbose = FALSE))
    names(x) <- y
    out_list <- list()
    send_info("Processing sample {.code ", sample, "}.")
    for (m in methods) {
      out <- sig_fit_bootstrap(x, n = n, method = m, ...)
      out_list[[m]] <- out
    }
    if (!is.null(job_id)) {
      saveRDS(out_list, file = fpath)
      return(fpath)
    } else {
      return(out_list)
    }
  }

  if (use_parallel) {
    oplan <- future::plan()
    future::plan("multiprocess", workers = future::availableCores())
    # furrr::future_options(seed = TRUE)
    on.exit(future::plan(oplan), add = TRUE)
    bt_list <- furrr::future_map2(as.data.frame(catalogue_matrix), colnames(catalogue_matrix), call_bt,
      y = rownames(catalogue_matrix), methods = methods, n = n, ..., .progress = TRUE,
      .options = furrr::future_options(seed = seed) # set options
    )
  } else {
    bt_list <- purrr::map2(as.data.frame(catalogue_matrix), colnames(catalogue_matrix), call_bt,
      y = rownames(catalogue_matrix), methods = methods, n = n, ...
    )
  }

  if (is.character(bt_list[[1]])) {
    ## Assume file paths
    bt_list <- purrr::map(bt_list, function(x) {
      readRDS(x)
    })
  }

  send_success("Gotten.")

  send_info("Reporting p values...")
  p_val <- purrr::map(bt_list, report_bootstrap_p_value, thresholds = p_val_thresholds)
  send_success("Done.")

  send_info("Cleaning results...")

  ## Clean exposures and errors
  bt_list <- purrr::map(bt_list, function(x) {
    x <- purrr::map(x, function(out) {
      out$expo <- out$expo %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "type") %>%
        tidyr::pivot_longer(-"type", names_to = "sig", values_to = "exposure") %>%
        data.table::as.data.table()
      out$errors <- data.table::data.table(
        type = names(out$errors),
        errors = as.numeric(out$errors)
      )
      out$cosine <- data.table::data.table(
        type = names(out$cosine),
        cosine = as.numeric(out$cosine)
      )
      out
    })
    x <- purrr::transpose(x)
    x <- purrr::map(x, data.table::rbindlist, fill = TRUE, idcol = "method")
    x
  })

  bt_list <- purrr::transpose(bt_list)
  bt_list <- purrr::map(bt_list, data.table::rbindlist, fill = TRUE, idcol = "sample")

  ## Clean p values
  p_val <- purrr::map(p_val, function(x) {
    x <- purrr::map(x, function(out) {
      out <- out %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "threshold") %>%
        tidyr::pivot_longer(-"threshold", names_to = "sig", values_to = "p_value") %>%
        data.table::as.data.table()
      out
    })
    x <- data.table::rbindlist(x, fill = TRUE, idcol = "method")
    x
  })
  p_val <- data.table::rbindlist(p_val, fill = TRUE, idcol = "sample")
  p_val$threshold <- sub("threshold_", "", p_val$threshold)

  ## Outputing
  send_success("Outputing.")
  optimal_list$expo$type <- "optimal"
  optimal_list$errors$type <- "optimal"
  optimal_list$cosine$type <- "optimal"

  expos <- rbind(optimal_list$expo, bt_list$expo, fill = TRUE)
  errors <- rbind(optimal_list$errors, bt_list$errors, fill = TRUE)
  cosine <- rbind(optimal_list$cosine, bt_list$cosine, fill = TRUE)

  result <- list(
    expo = expos,
    error = errors,
    cosine = cosine,
    p_val = p_val
  )
  return(result)
}
