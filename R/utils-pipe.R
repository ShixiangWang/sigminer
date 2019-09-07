#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL


# Source code from reporttools (https://github.com/cran/reporttools/blob/master/R/pairwise.fisher.test.r)
pairwise.fisher.test <- function(x, g, p.adjust.method, ...) {
  DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(g)))
  g <- factor(g)

  compare.levels <- function(i, j) {
    ind <- (as.integer(g) %in% c(i, j)) & (is.na(x) == FALSE) & (is.na(g) == FALSE)
    xi <- factor(x[ind], exclude = NULL)
    xj <- factor(g[ind], exclude = NULL)
    tab <- table(xi, xj)
    nonzeromarginal <- (min(apply(tab, 1, sum)) * min(apply(tab, 2, sum)) > 0)
    size <- ((nrow(tab) > 1) * (ncol(tab) > 1) > 0)
    if ((nonzeromarginal == TRUE) & (size == TRUE)) {
      fisher.test(xi, xj, ...)$p.value
    } else {
      NA
    }
  }

  PVAL <- pairwise.table(compare.levels, levels(g), p.adjust.method)
  ans <- list(data.name = DNAME, p.value = PVAL, p.adjust.method = p.adjust.method)
  class(ans) <- "pairwise.htest"
  return(ans)
}


# https://stackoverflow.com/questions/8396577/check-if-character-value-is-a-valid-r-object-name
isValidAndUnreserved <- function(string) {
  make.names(string) == string
}

`%:::%` <- function(pkg, fun) get(fun,
    envir = asNamespace(pkg),
    inherits = FALSE
  )

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
