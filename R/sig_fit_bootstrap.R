#' Obtain Bootstrap Distribution of Signature Exposures of a Certain Tumor Sample
#'
#' @inheritParams sig_fit
#' @param catalog a named numeric vector or a numeric matrix with dimension Nx1.
#' N is the number of component, 1 is the sample.
#' @param n the number of bootstrap replicates.
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
                              ...) {
  method <- match.arg(method)
  mode <- match.arg(mode)

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

  res <- replicate(n, {
    sampled <- sample(seq(K), total_count, replace = TRUE, prob = catalog / sum(catalog))
    catalog_mat <- as.integer(table(factor(sampled, levels = seq(K))))
    catalog_mat <- matrix(catalog_mat, ncol = 1)
    rownames(catalog_mat) <- names(catalog)

    suppressMessages(
      sig_fit(
        catalogue_matrix = catalog_mat,
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
        ...
      )
    )
  })

  expo <- res[1, ]
  expo <- sapply(expo, cbind)
  rownames(expo) <- rownames(res[1, ][[1]])

  errors <- sapply(res[2, ], c)
  names(errors) <- colnames(expo) <- paste0("Rep_", seq(n))

  return(list(expo = expo, errors = errors))
}
