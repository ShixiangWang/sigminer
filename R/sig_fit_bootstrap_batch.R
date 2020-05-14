#' Exposure Instability Analysis of Signature Exposures with Bootstrap
#'
#' @inheritParams sig_fit
#' @inheritParams sig_fit_bootstrap
#' @param methods a subset of `c("LS", "QP", "SA")`.
#' @param p_val_thresholds a vector of relative exposure threshold for calculating p values.
#' @param use_parallel if `TRUE`, use parallel computation based on **furrr** package,
#' not implemented yet.
#' @param seed random seed to reproduce the result.
#' @param ... other common parameters passing to [sig_fit_bootstrap], including `sig`, `sig_index`,
#' `sig_db`, `db_type`, `mode`, etc.
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
#' if (requireNamespace("lsei") & requireNamespace("quadprog")) {
#'   z <- sig_fit_bootstrap_batch(V, sig = W, n = 2)
#'   z
#'   z10 <- sig_fit_bootstrap_batch(V, sig = W, n = 10)
#' }
#' @testexamples
#' expect_is(z, "list")
#' z2 <- sig_fit_bootstrap_batch(V, sig = W, n = 2, use_parallel = TRUE)
#' expect_is(z2, "list")
sig_fit_bootstrap_batch <- function(catalogue_matrix, methods = c("LS", "QP"), n = 100L,
                                    p_val_thresholds = c(0.05),
                                    use_parallel = FALSE,
                                    seed = 123456L,
                                    ...) {
  stopifnot(is.matrix(catalogue_matrix))

  set.seed(seed, kind = "L'Ecuyer-CMRG")
  methods <- match.arg(methods, choices = c("LS", "QP", "SA"), several.ok = TRUE)

  timer <- Sys.time()
  send_info("Batch Bootstrap Signature Exposure Analysis Started.")
  on.exit(send_elapsed_time(timer, "Total "))

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
    optimal_list[[m]] <- expo_list
  }
  optimal_list <- purrr::transpose(optimal_list)
  optimal_list <- purrr::map(optimal_list, data.table::rbindlist, fill = TRUE, idcol = "method")

  ## Get bootstrap exposures with different methods
  send_info("Getting bootstrap exposures (&errors) for different methods.")
  send_info("This step is time consuming, please be patient.")
  call_bt <- function(x, sample, y, methods, n = 1000, ...) {
    names(x) <- y
    out_list <- list()
    send_info("Processing sample {.code ", sample, "}.")
    for (m in methods) {
      out <- sig_fit_bootstrap(x, n = n, method = m, ...)
      out_list[[m]] <- out
    }
    return(out_list)
  }

  if (use_parallel) {
    oplan <- future::plan()
    future::plan("multiprocess", workers = future::availableCores())
    furrr::future_options(seed = TRUE)
    on.exit(future::plan(oplan), add = TRUE)
    bt_list <- furrr::future_map2(as.data.frame(catalogue_matrix), colnames(catalogue_matrix), call_bt,
      y = rownames(catalogue_matrix), methods = methods, n = n, ..., .progress = TRUE
    )
  } else {
    bt_list <- purrr::map2(as.data.frame(catalogue_matrix), colnames(catalogue_matrix), call_bt,
      y = rownames(catalogue_matrix), methods = methods, n = n, ...
    )
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

  expos <- rbind(optimal_list$expo, bt_list$expo, fill = TRUE)
  errors <- rbind(optimal_list$errors, bt_list$errors, fill = TRUE)

  result <- list(
    expo = expos,
    error = errors,
    p_val = p_val
  )
  return(result)
}
