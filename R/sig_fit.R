#' Fit Signature Exposures with Linear Combination Decomposition
#'
#' The function performs a signatures decomposition of a given mutational
#' catalogue `V` with known signatures `W` by solving the minimization problem
#' `min(||W*H - V||)` where W and V are known.
#'
#' The method 'NNLS' solves the minimization problem with nonnegative least-squares constraints.
#' The method 'QP' and 'SA' are modified from SignatureEstimation package.
#' See references for details.
#' Of note, when fitting exposures for copy number signatures, only components of
#' feature CN is used.
#'
#' @param catalogue_matrix a numeric matrix `V` with row representing components and
#' columns representing samples, typically you can get `nmf_matrix` from `sig_tally()` and
#' transpose it by `t()`.
#' @inheritParams show_sig_profile
#' @inheritParams show_cosmic_sig_profile
#' @inheritParams get_sig_exposure
#' @inheritParams get_sig_similarity
#' @param sig a `Signature` object obtained either from [sig_extract] or [sig_auto_extract],
#' or just a raw signature matrix with row representing components (motifs) and
#' column representing signatures.
#' @param method method to solve the minimazation problem.
#' 'NNLS' for nonnegative least square; 'QP' for quadratic programming; 'SA' for simulated annealing.
#' @param auto_reduce if `TRUE`, try reducing the input reference signatures to increase
#' the cosine similarity of reconstructed profile to observed profile.
#' @param return_class string, 'matrix' or 'data.table'.
#' @param return_error if `TRUE`, also return sample error (Frobenius norm) and cosine
#' similarity between observed sample profile (asa. spectrum) and reconstructed profile. NOTE:
#' it is better to obtain the error when the type is 'absolute', because the error is
#' affected by relative exposure accuracy.
#' @param rel_threshold numeric vector, a relative exposure lower than this value will be set to 0.
#' Of note, this is a little different from the same parameter in [get_sig_exposure].
#' @param true_catalog used by [sig_fit_bootstrap], user never use it.
#' @param ... control parameters passing to argument `control` in `GenSA` function when use method 'SA'.
#'
#' @return The exposure result either in `matrix` or `data.table` format.
#' If `return_error` set `TRUE`, a `list` is returned.
#' @export
#' @seealso [sig_extract], [sig_auto_extract], [sig_fit_bootstrap], [sig_fit_bootstrap_batch]
#' @references
#' Daniel Huebschmann, Zuguang Gu and Matthias Schlesner (2019). YAPSA: Yet Another Package for Signature Analysis. R package version 1.12.0.
#'
#' Huang X, Wojtowicz D, Przytycka TM. Detecting presence of mutational signatures in cancer with confidence. Bioinformatics. 2018;34(2):330–337. doi:10.1093/bioinformatics/btx604
#'
#' Kim, Jaegil, et al. "Somatic ERCC2 mutations are associated with a distinct genomic signature in urothelial tumors."
#'  Nature genetics 48.6 (2016): 600.
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
#'   H_infer <- sig_fit(V, W, method = "QP")
#'   H_infer
#'   H
#'
#'   H_dt <- sig_fit(V, W, method = "QP", return_class = "data.table")
#'   H_dt
#'
#'   ## Show results
#'   show_sig_fit(H_infer)
#'   show_sig_fit(H_dt)
#'
#'   ## Get clusters/groups
#'   H_dt_rel <- sig_fit(V, W, return_class = "data.table", type = "relative")
#'   z <- get_groups(H_dt_rel, method = "k-means")
#'   show_groups(z)
#' }
#'
#' if (requireNamespace("GenSA", quietly = TRUE)) {
#'   H_infer <- sig_fit(V, W, method = "SA")
#'   H_infer
#'   H
#'
#'   H_dt <- sig_fit(V, W, method = "SA", return_class = "data.table")
#'   H_dt
#'
#'   ## Modify arguments to method
#'   sig_fit(V, W, method = "SA", maxit = 10, temperature = 100)
#'
#'   ## Show results
#'   show_sig_fit(H_infer)
#'   show_sig_fit(H_dt)
#' }
#' @testexamples
#' expect_is(H_infer, "matrix")
#' expect_is(H_dt, "data.table")
sig_fit <- function(catalogue_matrix,
                    sig,
                    sig_index = NULL,
                    sig_db = "legacy",
                    db_type = c("", "human-exome", "human-genome"),
                    show_index = TRUE,
                    method = c("QP", "NNLS", "SA"),
                    auto_reduce = FALSE,
                    type = c("absolute", "relative"),
                    return_class = c("matrix", "data.table"),
                    return_error = FALSE,
                    rel_threshold = 0,
                    mode = c("SBS", "DBS", "ID", "copynumber"),
                    true_catalog = NULL,
                    ...) {
  stopifnot(is.matrix(catalogue_matrix))
  db_type <- match.arg(db_type)
  method <- match.arg(method)

  timer <- Sys.time()
  send_info("Started.")
  on.exit(send_elapsed_time(timer))

  if (is.null(sig_index)) {
    send_info("Signature index not detected.")
    if (inherits(sig, "Signature")) {
      send_success("Signature object detected.")
      sig_matrix <- sig$Signature
    } else if (is.matrix(sig)) {
      send_success("Signature matrix detected.")
      sig_matrix <- sig
    } else {
      send_error("Invalid input for 'sig'.")
      send_stop("Exit.")
    }
  } else {
    send_success("Signature index detected.")
    send_info("Checking signature database in package.")

    sigs_db <- get_sig_db(sig_db)
    sigs <- sigs_db$db

    ## Some extra processing
    if (sig_db == "legacy" & db_type == "human-genome") {
      ## v2 comes from Exome
      sigs <- sig_convert(sig = sigs, from = "human-exome", to = "human-genome")
    } else if (sig_db == "SBS" & db_type == "human-exome") {
      ## v3 comes from WGS (PCAWG)
      ## Should DBS and ID also handle such cases?
      sigs <- sig_convert(sig = sigs, from = "human-genome", to = "human-exome")
    }

    avail_index <- switch(
      sig_db,
      legacy = substring(colnames(sigs), 8),
      SBS = substring(colnames(sigs), 4),
      DBS = substring(colnames(sigs), 4),
      ID = substring(colnames(sigs), 3),
      TSB = substring(colnames(sigs), 4)
    )

    send_info("Checking signature index.")

    msg <- paste(
      paste0("Valid index for db '", sig_db, "':"),
      paste0(avail_index, collapse = " "),
      sep = "\n"
    )
    if (show_index) {
      send_info(msg)
    }

    if (!is.character(sig_index)) {
      sig_index <- as.character(sig_index)
    }

    if ("ALL" %in% sig_index) {
      sig_index <- avail_index
    }

    if (length(sig_index) == 1) {
      if (grepl(",", sig_index)) {
        send_success("Signature indices separated by comma are detected.")
        sig_index <- split_seq(sig_index)
      }
    }

    if (!all(sig_index %in% avail_index)) {
      send_error("Invalid index.")
      send_info(msg)
      send_stop()
    }

    index <- c()
    for (i in sig_index) {
      index <- c(index, colnames(sigs)[avail_index == i])
    }

    sig_matrix <- as.matrix(sigs[, index, drop = FALSE])
  }

  send_success("Database and index checked.")

  mode <- match.arg(mode)
  type <- match.arg(type)
  return_class <- match.arg(return_class)

  if (mode == "copynumber") {
    send_info("Copy number mode detected. Checking and handling component names...")
    ## For copy number signature
    ## Only feature of copy number used for exposure quantification
    ## So the result exposure is estimated copy number segments
    has_cn_w <- grepl("^CN[^C]", rownames(sig_matrix)) | startsWith(rownames(sig_matrix), "copynumber")
    sig_matrix <- sig_matrix[has_cn_w, , drop = FALSE]
    catalogue_matrix <- catalogue_matrix[rownames(sig_matrix), , drop = FALSE]
    send_success("Component names checked and handled.")
  }

  ## Keep sum of signature to 1
  sig_matrix <- apply(sig_matrix, 2, function(x) x / sum(x))
  send_success("Signature normalized.")

  ## Check if V and W have same rows
  send_info("Checking row number for catalog matrix and signature matrix.")
  if (nrow(catalogue_matrix) != nrow(sig_matrix)) {
    send_stop("Catalogue matrix and Signature matrix should have same rows, please check!")
  }
  send_success("Checked.")

  ## If V and W have row names, check the order
  cat_rowname <- rownames(catalogue_matrix)
  sig_rowname <- rownames(sig_matrix)
  if (!is.null(cat_rowname) & !is.null(sig_rowname)) {
    send_info("Checking rownames for catalog matrix and signature matrix.")
    if (!all(sig_rowname == cat_rowname)) {
      send_info("Matrix V and W don't have same orders. Try reordering...")
      if (all(sort(cat_rowname) == sort(sig_rowname))) {
        ## Set catalogue matrix as signature matrix
        catalogue_matrix <- catalogue_matrix[sig_rowname, , drop = FALSE]
      } else {
        stop("The rownames of matrix are not identical, please check your input!")
      }
    }
    send_success("Checked.")
  }

  send_success("Method '", method, "' detected.")
  f_fit <- switch(method,
    NNLS = {
      if (!requireNamespace("nnls", quietly = TRUE)) {
        send_stop("Please install 'nnls' package firstly.")
      }
      decompose_NNLS
    },
    QP = {
      if (!requireNamespace("quadprog", quietly = TRUE)) {
        send_stop("Please install 'quadprog' package firstly.")
      }
      decompose_QP
    },
    SA = {
      if (!requireNamespace("GenSA", quietly = TRUE)) {
        send_stop("Please install 'GenSA' package firstly.")
      }
      decompose_SA
    }
  )
  send_success("Corresponding function generated.")

  send_info("Calling function.")
  # expo <- purrr::map2(as.data.frame(catalogue_matrix), rel_threshold,
  #   f_fit,
  #   sig_matrix,
  #   type = type,
  #   auto_reduce = auto_reduce,
  #   ...
  # )
  expo <- purrr::pmap(list(as.data.frame(catalogue_matrix),
                           rel_threshold,
                           if (is.null(colnames(catalogue_matrix))) {
                             NA_character_
                           } else {
                             colnames(catalogue_matrix)
                           } ),
                      f_fit,
                      sig_matrix,
                      type = type,
                      auto_reduce = auto_reduce,
                      ...
  )
  send_success("Done.")

  send_info("Generating output signature exposures.")
  expo <- dplyr::bind_rows(expo) %>%
    as.matrix()
  rownames(expo) <- colnames(sig_matrix)
  colnames(expo) <- colnames(catalogue_matrix)
  expo[is.na(expo)] <- 0

  expo_mat <- expo

  if (return_class != "matrix") {
    expo <- expo %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Sig") %>%
      tidyr::pivot_longer(cols = -"Sig", names_to = "sample", values_to = "expo") %>%
      tidyr::pivot_wider(id_cols = "sample", names_from = "Sig", values_from = "expo") %>%
      data.table::as.data.table()
  }

  send_success("Done.")
  if (return_error) {
    send_info("Calculating errors (Frobenius Norm).")
    if (!is.null(true_catalog)) {
      ## Make sure component names are same
      if (!is.null(sig_rowname) & !is.null(names(true_catalog))) {
        true_catalog <- true_catalog[sig_rowname]
      }
    }

    ## compute estimation error for each sample/patient (Frobenius norm)
    if (type == "relative") {
      send_warning("When the type is 'relative', the returned error is a little affected by its precision.")
      if (is.null(true_catalog)) {
        errors <- sapply(
          seq(ncol(expo_mat)),
          function(i) {
            FrobeniusNorm(
              catalogue_matrix[, i],
              sig_matrix,
              expo_mat[, i] * sum(catalogue_matrix[, i])
            )
          }
        )
      } else {
        errors <- sapply(
          seq(ncol(expo_mat)),
          function(i) {
            FrobeniusNorm(
              true_catalog,
              sig_matrix,
              expo_mat[, i] * sum(catalogue_matrix[, i])
            )
          }
        )
      }
    } else {
      if (is.null(true_catalog)) {
        errors <- sapply(seq(ncol(expo_mat)), function(i) FrobeniusNorm(catalogue_matrix[, i], sig_matrix, expo_mat[, i]))
      } else {
        errors <- sapply(seq(ncol(expo_mat)), function(i) FrobeniusNorm(true_catalog, sig_matrix, expo_mat[, i]))
      }
    }

    sim <- diag(cosine(sig_matrix %*% expo_mat, catalogue_matrix))
    names(sim) <- names(errors) <- colnames(catalogue_matrix)

    # Set precision
    errors <- round(errors, digits = 3)
    sim <- round(sim, digits = 6)

    send_success("Done.")
    return(list(
      expo = expo,
      errors = errors,
      cosine = sim
    ))
  }

  return(expo)
}


## x: catalogue to decompose
## y: relative exposure threshold
## sig_matrix: reference signature matrix, components X signatures
## type: type of signature contribution to return

decompose_NNLS <- function(x, y, z, sig_matrix, type = "absolute", auto_reduce = FALSE, ...) {
  if (is.na(z)) send_info("Fitting sample: ", z)

  if (sum(x) != 0) {
    ## nnls/lsqnonneg solve nonnegative least-squares constraints problem.
    ## expo <- pracma::lsqnonneg(sig_matrix, x)$x
    expo <- stats::coef(nnls::nnls(sig_matrix, x))
    expo <- expo / sum(expo)

    if (auto_reduce) {
      rec <- (expo %*% t(sig_matrix) * sum(x))[1, ]
      sim <- cosine(rec, x)
      if (sim < 0.99) {
        sim_old <- sim
        # continue to optimize
        send_info("Start optimizing...")
        for (i in seq(0.001, 0.501, 0.01)) {
          expo_low <- expo < i
          send_info("Dropping reference signatures with relative exposure <", i)
          sig_matrix_update <- sig_matrix[, !expo_low, drop = FALSE]
          abs_expo <- decompose_NNLS(x, 0, z, sig_matrix_update, type = "absolute")
          rec_update <- (abs_expo %*% t(sig_matrix_update))[1, ]
          sim_update <- cosine(rec_update, x)
          if (sim_update < sim) {
            break()
          }
          sim <- sim_update
        }
        send_success("Stop optimizing at exposure level: ",
                     i, ", ",
                     sum(expo_low),
                     " signatures dropped.")
        send_info("Cosine similarity with ", round(abs(sim - sim_old), 6), " improved.")

        out_expo <- vector("numeric", length = length(expo))
        ## Correctly assign the exposure
        out_expo[!expo_low] <- abs_expo / sum(abs_expo)
        out_expo <- out_expo / sum(out_expo)

        return(return_expo(out_expo, y, type, total = sum(x)))
      } else {
        send_success("The cosine similarity is very high, just return result.")
      }
    }

  } else {
    expo <- rep(0, ncol(sig_matrix))
  }

  return_expo(expo = expo, y, type, total = sum(x))
}

# m observed turmor profile vector for a single patient/sample, 96 by 1. m is normalized.
# P is same as sig_matrix

decompose_QP <- function(x, y, z, P, type = "absolute", auto_reduce = FALSE, ...) {
  if (is.na(z)) send_info("Fitting sample: ", z)

  if (sum(x) != 0) {
    m <- x / sum(x)
    # N: how many signatures are selected
    N <- ncol(P)
    # G: matrix appearing in the quatric programming objective function
    G <- t(P) %*% P
    # C: matrix constraints under which we want to minimize the quatric programming objective function.
    C <- cbind(rep(1, N), diag(N))
    # b: vector containing the values of b_0.
    b <- c(1, rep(0, N))
    # d: vector appearing in the quatric programming objective function
    d <- t(m) %*% P

    # Solve quadratic programming problem
    out <- quadprog::solve.QP(Dmat = G, dvec = d, Amat = C, bvec = b, meq = 1)

    # Some exposure values are negative, but very close to 0
    # Change these neagtive values to zero and renormalized
    expo <- out$solution
    expo[expo < 0] <- 0
    expo <- expo / sum(expo)

    if (auto_reduce) {
      sig_matrix <- P
      rec <- (expo %*% t(sig_matrix) * sum(x))[1, ]
      sim <- cosine(rec, x)
      if (sim < 0.99) {
        sim_old <- sim
        # continue to optimize
        send_info("Start optimizing...")
        for (i in seq(0.001, 0.501, 0.01)) {
          expo_low <- expo < i
          send_info("Dropping reference signatures with relative exposure <", i)
          sig_matrix_update <- sig_matrix[, !expo_low, drop = FALSE]
          abs_expo <- decompose_QP(x, 0, z, sig_matrix_update, type = "absolute")
          rec_update <- (abs_expo %*% t(sig_matrix_update))[1, ]
          sim_update <- cosine(rec_update, x)
          if (sim_update < sim) {
            break()
          }
          sim <- sim_update
        }
        send_success("Stop optimizing at exposure level: ",
                     i, ", ",
                     sum(expo_low),
                     " signatures dropped.")
        send_info("Cosine similarity with ", round(abs(sim - sim_old), 6), " improved.")

        out_expo <- vector("numeric", length = length(expo))
        ## Correctly assign the exposure
        out_expo[!expo_low] <- abs_expo / sum(abs_expo)
        out_expo <- out_expo / sum(out_expo)

        return(return_expo(out_expo, y, type, total = sum(x)))
      } else {
        send_success("The cosine similarity is very high, just return result.")
      }
    }

  } else {
    expo <- rep(0, ncol(P))
  }

  # return the exposures
  return_expo(expo, y, type, total = sum(x))
}


decompose_SA <- function(x, y, z, P, type = "absolute", auto_reduce = FALSE, ...) {
  if (is.na(z)) send_info("Fitting sample: ", z)

  if (sum(x) != 0) {
    control <- list(...)

    m <- x / sum(x)
    # objective function to be minimized
    # local version of Frobenius norm to simplify and speed-up the objective function
    FrobeniusNorm.local <- function(exposures) {
      estimate <- P %*% exposures
      return(sqrt(sum((m - (estimate / sum(estimate)))^2)))
    }
    # N: how many signatures are selected
    N <- ncol(P)
    # change our suggestion to control GenSA function based on user's requirements
    # https://blog.csdn.net/georgesale/article/details/80631417
    our.control <- list(maxit = 1000, temperature = 10, nb.stop.improvement = 1000, simple.function = TRUE)
    our.control[names(control)] <- control
    # Solve the problem using simulated annealing package GenSA
    sa <- GenSA::GenSA(lower = rep(0.0, N), upper = rep(1.0, N), fn = FrobeniusNorm.local, control = our.control)
    # Normalize the solution
    expo <- sa$par / sum(sa$par)

    if (auto_reduce) {
      sig_matrix <- P
      rec <- (expo %*% t(sig_matrix) * sum(x))[1, ]
      sim <- cosine(rec, x)
      if (sim < 0.99) {
        sim_old <- sim
        # continue to optimize
        send_info("Start optimizing...")
        for (i in seq(0.001, 0.501, 0.01)) {
          expo_low <- expo < i
          send_info("Dropping reference signatures with relative exposure <", i)
          sig_matrix_update <- sig_matrix[, !expo_low, drop = FALSE]
          abs_expo <- decompose_SA(x, 0, z, sig_matrix_update, type = "absolute")
          rec_update <- (abs_expo %*% t(sig_matrix_update))[1, ]
          sim_update <- cosine(rec_update, x)
          if (sim_update < sim) {
            break()
          }
          sim <- sim_update
        }
        send_success("Stop optimizing at exposure level: ",
                     i, ", ",
                     sum(expo_low),
                     " signatures dropped.")
        send_info("Cosine similarity with ", round(abs(sim - sim_old), 6), " improved.")

        out_expo <- vector("numeric", length = length(expo))
        ## Correctly assign the exposure
        out_expo[!expo_low] <- abs_expo / sum(abs_expo)
        out_expo <- out_expo / sum(out_expo)

        return(return_expo(out_expo, y, type, total = sum(x)))
      } else {
        send_success("The cosine similarity is very high, just return result.")
      }
    }

  } else {
    expo <- rep(0, ncol(P))
  }

  # return the exposures
  return_expo(expo, y, type, total = sum(x))
}

## total is used to set the total exposure in a sample
## for method QP and SA
return_expo <- function(expo, y, type = "absolute", total = NULL) {
  rel_expo <- expo
  expo <- expo * total
  expo[rel_expo < y] <- 0
  if (type == "relative") {
    expo <- expo / sum(expo)
  }
  expo <- round(expo, digits = 6)
  expo
}
