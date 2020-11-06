fitComponent <-
  function(dat,
           dist = "norm",
           seed = 123456,
           model_selection = "BIC",
           threshold = 0.1,
           min_prior = 0.001,
           niter = 1000,
           nrep = 1,
           min_comp = 2,
           max_comp = 10,
           cores = 1) {
    control <- new("FLXcontrol")
    control@minprior <- min_prior
    control@iter.max <- niter

    old_opt <- options()
    options(future.rng.onMisue = "ignore")
    on.exit(options(old_opt))
    set.seed(seed, kind = "L'Ecuyer-CMRG")

    stepFlexmix_v2 <- function(..., k = NULL, nrep = 3, verbose = TRUE, drop = TRUE,
                               unique = FALSE, cores = 1, seed = 123456) {
      # doParallel doest not work in Windows
      # doParallel::registerDoParallel(cores = cores)

      if (!requireNamespace("doFuture", quietly = TRUE)) {
        send_info("{.pkg doFuture} is recommended to install for improving computation.")
      } else {
        doFuture::registerDoFuture()
        future::plan("multiprocess", workers = cores)
      }

      MYCALL <- match.call()
      MYCALL1 <- MYCALL
      bestFlexmix <- function(...) {
        z <- new("flexmix", logLik = -Inf)
        logLiks <- rep(NA, length.out = nrep)
        for (m in seq_len(nrep)) {
          if (verbose) {
            cat(" *")
          }
          x <- try(flexmix::flexmix(...))
          if (!is(x, "try-error")) {
            logLiks[m] <- flexmix::logLik(x)
            if (flexmix::logLik(x) > flexmix::logLik(z)) {
              z <- x
            }
          }
        }
        return(list(z = z, logLiks = logLiks))
      }
      z <- list()
      if (is.null(k)) {
        RET <- bestFlexmix(...)
        z[[1]] <- RET$z
        logLiks <- as.matrix(RET$logLiks)
        z[[1]]@call <- MYCALL
        z[[1]]@control@nrep <- nrep
        names(z) <- as.character(z[[1]]@k)
        if (verbose) {
          cat("\n")
        }
      }
      else {
        k <- as.integer(k)
        # logLiks <- matrix(nrow = length(k), ncol = nrep)
        z_list <- foreach(
          n = seq_along(k),
          .packages = "flexmix",
          .export = c("k", "z", "MYCALL1", "nrep", "verbose", "bestFlexmix", "seed", "...")
        ) %dopar% {
          set.seed(seed, kind = "L'Ecuyer-CMRG")
          # loadNamespace("flexmix")
          ns <- as.character(k[n])
          if (verbose) {
            cat(k[n], ":")
          }
          RET <- bestFlexmix(..., k = k[n])
          z[[ns]] <- RET$z
          # logLiks[n, ] <- RET$logLiks
          MYCALL1[["k"]] <- as.numeric(k[n])
          z[[ns]]@call <- MYCALL1
          z[[ns]]@control@nrep <- nrep
          if (verbose) {
            cat("\n")
          }
          list(
            z = z,
            logLiks = RET$logLiks
          )
        }
      }

      z <- lapply(z_list, function(x) x[["z"]])
      z <- purrr::flatten(z)

      # k by nrep matrix
      if (nrep == 1) {
        logLiks <- as.matrix(sapply(z_list, function(x) x[["logLiks"]]))
      } else {
        logLiks <- t(sapply(z_list, function(x) x[["logLiks"]]))
      }

      finite_index <- is.finite(sapply(z, flexmix::logLik))
      logLiks <- logLiks[finite_index, , drop = FALSE]
      z <- z[finite_index]
      rownames(logLiks) <- names(z)
      if (!length(z)) {
        stop("No convergence to a suitable mixture")
      }
      if (drop & (length(z) == 1)) {
        return(z[[1]])
      }
      else {
        z <- return(new("stepFlexmix",
          models = z, k = as.integer(names(z)),
          nrep = as.integer(nrep), logLiks = logLiks, call = MYCALL
        ))
        if (unique) {
          z <- unique(z)
        }
        return(z)
      }
    }


    if (dist == "norm") {
      if (min_comp == max_comp) {
        fit <-
          flexmix::flexmix(
            dat ~ 1,
            model = flexmix::FLXMCnorm1(),
            k = min_comp,
            control = control
          )
      } else {
        fit <-
          stepFlexmix_v2(
            dat ~ 1,
            model = flexmix::FLXMCnorm1(),
            k = min_comp:max_comp,
            nrep = nrep,
            control = control,
            cores = cores,
            seed = 123456
          )
        if (inherits(fit, "stepFlexmix")) {
          fit <- recur_fit_component(
            fit = fit, dist = dist,
            threshold = threshold,
            control = control,
            model_selection = model_selection
          )
        }
      }
    } else if (dist == "pois") {
      if (min_comp == max_comp) {
        fit <-
          flexmix::flexmix(
            dat ~ 1,
            model = flexmix::FLXMCmvpois(),
            k = min_comp,
            control = control
          )
      } else {
        fit <-
          stepFlexmix_v2(
            dat ~ 1,
            model = flexmix::FLXMCmvpois(),
            k = min_comp:max_comp,
            nrep = nrep,
            control = control,
            cores = cores,
            seed = 123456
          )
        if (inherits(fit, "stepFlexmix")) {
          fit <- recur_fit_component(
            fit = fit, dist = dist,
            threshold = threshold,
            control = control,
            model_selection = model_selection
          )
        }
      }
    }
    fit
  }

recur_fit_component <- function(fit, dist, threshold, control, model_selection = "BIC") {
  fits <- fit
  fit <- flexmix::getModel(fits, which = model_selection)
  send_success("Select ", fit@k, " according to ", model_selection, ".")

  mu <- find_mu(fit)
  sub_len <- sum(diff(mu) < threshold)

  if (sub_len > 0) {
    send_info("The model does not pass the threshold for mu difference of adjacent distribution.")
    send_info("Trying to call the optimal model...")
  }

  while (sub_len > 0) {
    K <- fit@k
    fit <- tryCatch(flexmix::getModel(fits, which = as.character(K - sub_len)),
      error = function(e) {
        # If the model is not called by user
        # Call it
        if (dist == "norm") {
          flexmix::flexmix(
            dat ~ 1,
            model = flexmix::FLXMCnorm1(),
            k = K,
            control = control
          )
        } else {
          flexmix::flexmix(
            dat ~ 1,
            model = flexmix::FLXMCmvpois(),
            k = K,
            control = control
          )
        }
      }
    )

    mu <- find_mu(fit)
    sub_len <- sum(diff(mu) < threshold)
  }

  send_success("Finally, select ", fit@k, " after passing threshold ", threshold, ".")
  fit
}

find_mu <- function(fit) {
  if (is.matrix(flexmix::parameters(fit))) {
    mu <- flexmix::parameters(fit)[1, ]
  } else {
    mu <- flexmix::parameters(fit)
  }
  mu <- sort(mu, decreasing = FALSE)
  mu
}

utils::globalVariables(
  c("n")
)
