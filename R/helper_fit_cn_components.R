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
             max_comp = 10) {
    control <- new("FLXcontrol")
    control@minprior <- min_prior
    control@iter.max <- niter

    set.seed(seed)

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
          stepFlexmix(
            dat ~ 1,
            model = flexmix::FLXMCnorm1(),
            k = min_comp:max_comp,
            nrep = nrep,
            control = control
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
          stepFlexmix(
            dat ~ 1,
            model = flexmix::FLXMCmvpois(),
            k = min_comp:max_comp,
            nrep = nrep,
            control = control
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
  message("Select ", fit@k, " according to ", model_selection)

  mu <- find_mu(fit)
  sub_len <- sum(diff(mu) < threshold)

  if (sub_len > 0) {
    message("The model does not pass the threshold for mu difference of adjacent distribution")
    message("Trying to call the optimal model...")
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

  message("Finally, select ", fit@k, " after passing threshold ", threshold)
  message("===================")
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
