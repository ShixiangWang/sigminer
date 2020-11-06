#' Obtain Bootstrap Distribution of Signature Exposures of a Certain Tumor Sample
#'
#' This can be used to obtain the confidence of signature exposures or search
#' the suboptimal decomposition solution.
#'
#' @inheritParams sig_fit
#' @param catalog a named numeric vector or a numeric matrix with dimension Nx1.
#' N is the number of component, 1 is the sample.
#' @param n the number of bootstrap replicates.
#' @param SA_not_bootstrap if `TRUE`, directly run 'SA' multiple times with original input instead of
#' bootstrap samples.
#' @param find_suboptimal logical, if `TRUE`, find suboptimal decomposition with
#' slightly higher error than the optimal solution by method 'SA'. This is useful
#' to explore hidden dependencies between signatures. More see reference.
#' @param suboptimal_ref_error baseline error used for finding suboptimal solution.
#' if it is `NULL`, then use 'SA' method to obtain the optimal error.
#' @param suboptimal_factor suboptimal factor to get suboptimal error, default is `1.05`,
#' i.e., suboptimal error is `1.05` times baseline error.
#' @references Huang X, Wojtowicz D, Przytycka TM. Detecting presence of mutational signatures in cancer with confidence. Bioinformatics. 2018;34(2):330–337. doi:10.1093/bioinformatics/btx604
#' @return a `list`
#' @export
#' @keywords bootstrap
#' @seealso [report_bootstrap_p_value], [sig_fit], [sig_fit_bootstrap_batch]
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
#' if (requireNamespace("quadprog", quietly = TRUE)) {
#'   H_bootstrap <- sig_fit_bootstrap(V[, 1], W, n = 10, type = "absolute")
#'   ## Typically, you have to run many times to get close to the answer
#'   boxplot(t(H_bootstrap$expo))
#'   H[, 1]
#'
#'   ## Return P values
#'   ## In practice, run times >= 100
#'   ## is recommended
#'   report_bootstrap_p_value(H_bootstrap)
#'   ## For multiple samples
#'   ## Input a list
#'   report_bootstrap_p_value(list(samp1 = H_bootstrap, samp2 = H_bootstrap))
#'
#   ## Find suboptimal decomposition
#   H_suboptimal <- sig_fit_bootstrap(V[, 1], W,
#     n = 10,
#     type = "absolute",
#     method = "SA",
#     find_suboptimal = TRUE
#   )
#' }
#' @testexamples
#' expect_is(H_bootstrap, "list")
# expect_is(H_suboptimal, "list")

sig_fit_bootstrap <- function(catalog,
                              sig,
                              n = 100L,
                              sig_index = NULL,
                              sig_db = "legacy",
                              db_type = c("", "human-exome", "human-genome"),
                              show_index = TRUE,
                              method = c("QP", "NNLS", "SA"),
                              auto_reduce = FALSE,
                              SA_not_bootstrap = FALSE,
                              type = c("absolute", "relative"),
                              rel_threshold = 0,
                              mode = c("SBS", "DBS", "ID", "copynumber"),
                              find_suboptimal = FALSE,
                              suboptimal_ref_error = NULL,
                              suboptimal_factor = 1.05,
                              ...) {
  method <- match.arg(method)
  mode <- match.arg(mode)

  timer <- Sys.time()
  send_info("Started.")
  on.exit(send_elapsed_time(timer))

  if (find_suboptimal) {
    send_info("Suboptimal mode enabled.")
    if (method != "SA") {
      send_warning("Only method 'SA' can be used to find suboptimal decomposition. Resetting method.")
      method <- "SA"
    }
  }

  send_info("Checking catalog.")

  if (is.matrix(catalog)) {
    if (ncol(catalog) != 1) {
      stop("'catalog' can only have one sample!")
    }
    temp <- catalog[, 1]
    names(temp) <- rownames(catalog)
    catalog <- temp
  } else if (!is.numeric(catalog)) {
    stop("Invalid input for 'catalog'.")
  }

  if (mode == "copynumber") {
    ## For copy number signature
    ## Only feature of copy number used for exposure quantification
    ## So the result exposure is estimated copy number segments
    has_cn_w <- grepl("^CN[^C]", names(catalog)) | startsWith(names(catalog), "copynumber")
    catalog <- catalog[has_cn_w]
  }

  send_success("Done.")

  ## To here, catalog is a named vector
  total_count <- sum(catalog)
  K <- length(catalog) # number of mutation types

  if (find_suboptimal & is.null(suboptimal_ref_error)) {
    send_info("Running 'SA' method to getting optimal error as baseline error.")
    catalog_mat <- matrix(catalog, ncol = 1)
    rownames(catalog_mat) <- names(catalog)
    optimal_res <- sig_fit(
      catalogue_matrix = catalog_mat,
      sig = sig,
      sig_index = sig_index,
      sig_db = sig_db,
      db_type = db_type,
      show_index = show_index,
      method = "SA",
      type = type,
      return_class = "matrix",
      return_error = TRUE,
      rel_threshold = rel_threshold,
      mode = mode,
      true_catalog = catalog,
      ...
    )
    suboptimal_ref_error <- optimal_res$errors
    if (suboptimal_ref_error == 0) {
      send_info("The optimal error is 0, resetting it to 1% of total exposure.")
      suboptimal_ref_error <- 0.01 * sum(catalog)
    }
    send_success("Got baseline error.")
  }

  Args_add <- list(...)
  if (!is.null(sig_index)) {
    sig <- NA
  }
  if (find_suboptimal) {
    threshold.stop <- suboptimal_factor * suboptimal_ref_error
  } else {
    threshold.stop <- NULL
  }

  Args <- list(
    catalogue_matrix = NA,
    sig = sig,
    sig_index = sig_index,
    sig_db = sig_db,
    db_type = db_type,
    show_index = show_index,
    method = method,
    auto_reduce = auto_reduce,
    type = type,
    return_class = "matrix",
    return_error = TRUE,
    rel_threshold = rel_threshold,
    mode = mode,
    true_catalog = catalog,
    threshold.stop = threshold.stop
  )
  Args <- c(Args, Args_add)

  send_info("About to start bootstrap.")
  sb <- cli::cli_status("{symbol$arrow_right} Bootstrapping {n} times.")

  if (SA_not_bootstrap && method == "SA") {
    send_info("'SA' method is detected & 'SA_not_bootstrap' is TRUE, run 'SA' multiple times with original catalogues.")
    catalog_mat <- matrix(catalog, ncol = 1)
    rownames(catalog_mat) <- names(catalog)
  }

  res <- sapply(1:n, function(i) {
    if (any(!SA_not_bootstrap, method != "SA")) {
      sampled <- sample(seq(K), total_count, replace = TRUE, prob = catalog / sum(catalog))
      catalog_mat <- as.integer(table(factor(sampled, levels = seq(K))))
      catalog_mat <- matrix(catalog_mat, ncol = 1)
      rownames(catalog_mat) <- names(catalog)
    }

    Args$catalogue_matrix <- catalog_mat

    cli::cli_status_update(id = sb, "{symbol$arrow_right} Total {n} times, starting no.{i}.")

    suppressMessages(
      base::do.call("sig_fit", args = Args)
    )
  })

  cli::cli_status_clear(sb)
  send_success("Bootstrap done.")

  expo <- res[1, ]
  expo <- sapply(expo, cbind)
  rownames(expo) <- rownames(res[1, ][[1]])

  send_success("Signature exposures collected.")

  errors <- sapply(res[2, ], c)
  cosine <- sapply(res[3, ], c)
  names(cosine) <- names(errors) <- colnames(expo) <- paste0("Rep_", seq(n))
  send_success("Errors and similarity collected.")

  send_success("Done.")

  return(list(expo = expo, errors = errors, cosine = cosine))
}


## Handle one or more samples
#' Report P Values from bootstrap Results
#'
#' See examples in [sig_fit_bootstrap].
#'
#' @param x a (list of) result from [sig_fit_bootstrap].
#' @param thresholds a vector of relative exposure threshold for calculating p values.
#'
#' @return a (list of) `matrix`
#' @export
report_bootstrap_p_value <- function(x, thresholds = c(0.01, 0.05, 0.1)) {
  stopifnot(is.list(x), is.numeric(thresholds))

  timer <- Sys.time()
  send_info("Started.")
  on.exit(send_elapsed_time(timer))

  if (!is.null(x$expo)) {
    send_success("Single sample mode enabled.")
    expo <- x$expo
    y <- bootstrap_p_value(expo, thresholds)
  } else {
    send_success("Batch mode enabled.")
    y <- lapply(x, function(xi, y) {
      bootstrap_p_value(xi$expo, y)
    }, y = thresholds)
  }

  send_success("Done.")
  return(y)
}

## Handle one sample
bootstrap_p_value <- function(x, y) {
  # if (any(is.null(rownames(x)), is.null(colnames(x)))) {
  #   send_stop("Exposure matrix should have both rownames and colnames.")
  # }
  x <- x / colSums(x)
  y_mat <- sapply(y, function(yi) {
    apply(x, 1, function(xi) {
      p <- my.t.test.p.value(xi, mu = yi, alternative = "greater")
      if (is.na(p)) {
        send_warning("NA result detected from t.test, reporting proportion as p value.")
        p <- mean(xi <= yi)
      }
      p
    })

  })
  colnames(y_mat) <- paste0("threshold_", y)
  return(y_mat)
}

my.t.test.p.value <- function(...) {
  obj <- try(t.test(...), silent = TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}
