#' Obtain Bootstrap Distribution of Signature Exposures of a Certain Tumor Sample
#'
#' This can be used to obtain the confidence of signature exposures or search
#' the suboptimal decomposition solution.
#'
#' @inheritParams sig_fit
#' @param catalog a named numeric vector or a numeric matrix with dimension Nx1.
#' N is the number of component, 1 is the sample.
#' @param n the number of bootstrap replicates.
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
#' if (requireNamespace("lsei", quietly = TRUE)) {
#'   H_bootstrap <- sig_fit_bootstrap(V[, 1], W, n = 100, type = "absolute")
#'   ## Typically, you have to run many times to get close to the answer
#'   boxplot(t(H_bootstrap$expo))
#'   H[, 1]
#' }
#' @testexamples
#' expect_is(H_bootstrap, "list")

sig_fit_bootstrap <- function(catalog,
                              sig,
                              n = 100L,
                              sig_index = NULL,
                              sig_db = "legacy",
                              db_type = c("", "human-exome", "human-genome"),
                              show_index = TRUE,
                              method = c("LS", "QP", "SA"),
                              type = c("absolute", "relative"),
                              rel_threshold = 0,
                              mode = c("SBS", "copynumber"),
                              find_suboptimal = FALSE,
                              suboptimal_ref_error = NULL,
                              suboptimal_factor = 1.05,
                              ...) {
  method <- match.arg(method)
  mode <- match.arg(mode)

  if (find_suboptimal) {
    if (method != "SA") {
      warning("Only method 'SA' can be used to find suboptimal decomposition.\n Resetting method.", immediate. = TRUE)
      method <- "SA"
    }
  }

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

  ## To here, catalog is a named vector
  total_count <- sum(catalog)
  K <- length(catalog) # number of mutation types

  if (find_suboptimal & is.null(suboptimal_ref_error)) {
    message("=> Running 'SA' method to getting optimal error as baseline error...")
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
      message("=> The optimal error is 0, resetting it to 1% of total exposure")
      suboptimal_ref_error <- 0.01 * sum(catalog)
    }
  }


  message("=> Running bootstrap...")

  Args_add = list(...)

  if (!is.null(sig_index)) {
    sig = NA
  }
  if (find_suboptimal) {
    threshold.stop = suboptimal_factor * suboptimal_ref_error
  } else {
    threshold.stop = NULL
  }

  Args <- list(
    catalogue_matrix = NA,
    sig = sig,
    sig_index = sig_index,
    sig_db = sig_db,
    db_type = db_type,
    show_index = show_index,
    method = method,
    type = type,
    return_class = "matrix",
    return_error = TRUE,
    rel_threshold = rel_threshold,
    mode = mode,
    true_catalog = catalog,
    threshold.stop = threshold.stop
  )
  Args = c(Args, Args_add)

  res <- replicate(n, {
    sampled <- sample(seq(K), total_count, replace = TRUE, prob = catalog / sum(catalog))
    catalog_mat <- as.integer(table(factor(sampled, levels = seq(K))))
    catalog_mat <- matrix(catalog_mat, ncol = 1)
    rownames(catalog_mat) <- names(catalog)

    Args$catalogue_matrix <- catalog_mat

    suppressMessages(
      base::do.call("sig_fit", args = Args)
    )

    # suppressMessages(
    #   sig_fit(
    #     catalogue_matrix = catalog_mat,
    #     sig = sig,
    #     sig_index = sig_index,
    #     sig_db = sig_db,
    #     db_type = db_type,
    #     show_index = show_index,
    #     method = method,
    #     type = type,
    #     return_class = "matrix",
    #     return_error = TRUE,
    #     rel_threshold = rel_threshold,
    #     mode = mode,
    #     true_catalog = catalog,
    #     unlist(c(list(...), threshold.stop = suboptimal_factor * suboptimal_ref_error))
    #   )
    # )
  })

  message("=> Done")

  expo <- res[1, ]
  expo <- sapply(expo, cbind)
  rownames(expo) <- rownames(res[1, ][[1]])

  errors <- sapply(res[2, ], c)
  names(errors) <- colnames(expo) <- paste0("Rep_", seq(n))

  return(list(expo = expo, errors = errors))
}
