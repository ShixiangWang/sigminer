
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

`%:::%` <- function(pkg, fun) {
  get(fun,
    envir = asNamespace(pkg),
    inherits = FALSE
  )
}

# From https://gist.github.com/mbannert/e9fcfa86de3b06068c83
rgb2hex <- function(r, g, b) grDevices::rgb(r, g, b, maxColorValue = 255)
col2hex <- function(col, alpha) grDevices::rgb(t(grDevices::col2rgb(col)), alpha = alpha, maxColorValue = 255)


#' Set Color Style for Plotting
#'
#' @param style one of 'default' and 'cosmic'.
#' @param mode only used when the `style` is 'cosmic', can be one of
#' "SBS", "copynumber", "DBS", "ID".
#'
#' @return color values.
#' @export
#'
#' @examples
#' use_color_style("default")
#' use_color_style("cosmic")
use_color_style <- function(style, mode = c("SBS", "copynumber", "DBS", "ID")) {
  # c("red", "cyan", "yellow", "blue", "magenta", "gray50", "orange", "darkgreen", "brown", "black", rainbow(10)[4:10])
  if (style == "default") {
    palette <- c(
      "cyan", "red", "brown", "purple",
      "green", "blue", "black", "magenta",
      "orange", "gray50", "darkgreen", "yellow",
      rainbow(10)[4:10]
    )
  } else {
    mode <- match.arg(mode)

    colors <- switch(mode,
      SBS = list(
        c(3, 189, 239), c(1, 1, 1),
        c(228, 41, 38), c(203, 202, 202),
        c(162, 207, 99), c(236, 199, 197)
      ),
      copynumber = list(
        c(3, 189, 239), c(1, 1, 1),
        c(228, 41, 38), c(203, 202, 202),
        c(162, 207, 99), c(236, 199, 197)
      ),
      DBS = list(
        c(3, 189, 239), c(3, 102, 204),
        c(162, 207, 99), c(1, 102, 1),
        c(255, 153, 153), c(228, 41, 38),
        c(255, 178, 102), c(255, 128, 1),
        c(204, 153, 255), c(76, 1, 153)
      ),
      ID = list(
        c(253, 190, 111), c(255, 128, 2),
        c(176, 221, 139), c(54, 161, 46),
        c(253, 202, 181), c(252, 138, 106),
        c(241, 68, 50), c(188, 25, 26),
        c(208, 225, 242), c(148, 196, 223),
        c(74, 152, 201), c(23, 100, 171),
        c(226, 226, 239), c(182, 182, 216),
        c(134, 131, 189), c(98, 64, 155)
      )
    )

    # https://github.com/AlexandrovLab/SigProfilerPlotting/blob/a60c48c2b74f037b9df5a3f905cc6ada22623341/sigProfilerPlotting/sigProfilerPlotting.py#L96
    # colors = [[3/256,189/256,239/256], [1/256,1/256,1/256],
    #           [228/256,41/256,38/256], [203/256,202/256,202/256], [162/256,207/256,99/256], [236/256,199/256,197/256]]
    #
    # For ID
    # https://github.com/AlexandrovLab/SigProfilerPlotting/blob/a60c48c2b74f037b9df5a3f905cc6ada22623341/sigProfilerPlotting/sigProfilerPlotting.py#L1558
    # colors = [[253/256,190/256,111/256], [255/256,128/256,2/256], [176/256,221/256,139/256], [54/256,161/256,46/256],
    #           [253/256,202/256,181/256], [252/256,138/256,106/256], [241/256,68/256,50/256], [188/256,25/256,26/256],
    #           [208/256,225/256,242/256], [148/256,196/256,223/256], [74/256,152/256,201/256], [23/256,100/256,171/256],
    #           [226/256,226/256,239/256], [182/256,182/256,216/256], [134/256,131/256,189/256], [98/256,64/256,155/256]]
    #
    #
    # For DBS
    # colors = [[3/256,189/256,239/256], [3/256,102/256,204/256],[162/256,207/256,99/256],
    #           [1/256,102/256,1/256], [255/256,153/256,153/256], [228/256,41/256,38/256],
    #           [255/256,178/256,102/256], [255/256,128/256,1/256], [204/256,153/256,255/256],
    #           [76/256,1/256,153/256]]

    palette <- sapply(colors, FUN = function(x) rgb2hex(x[1], x[2], x[3])) %>% as.character()

    if (mode %in% c("copynumber", "SBS")) {
      palette <- c(palette, sapply(c("purple", "brown", "orange"), FUN = col2hex) %>% as.character())
    }
  }
  palette
}


# https://stackoverflow.com/questions/34096162/dplyr-mutate-replace-several-columns-on-a-subset-of-rows
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>%
    dplyr::mutate(...)
  .data
}


FrobeniusNorm <- function(M, P, E) {
  sqrt(sum((M - P %*% E)^2))
}

is.wholenumber <- function(x, tol = .Machine$double.eps) {
  abs(x - round(x)) < tol
}

# Signals ----------------------------------------------------------------

send_info <- function(...) {
  cli::cli_alert_info(
    cli::style_bold(get_timestamp(), cli::col_blue(...))
  )
}

send_success <- function(...) {
  cli::cli_alert_success(
    cli::style_bold(get_timestamp(), cli::col_green(...))
  )
}

send_warning <- function(...) {
  cli::cli_alert_warning(
    cli::style_bold(get_timestamp(), cli::col_yellow(...))
  )
}

send_error <- function(...) {
  cli::cli_alert_danger(
    cli::style_bold(get_timestamp(), cli::col_red(...))
  )
}

send_elapsed_time <- function(timer, ...) {
  timer <- Sys.time() - timer
  send_info(..., round(timer, 3), " ", attr(timer, "units"), " elapsed.")
}

get_timestamp <- function() {
  paste0("[", Sys.time(), "]: ")
}

# https://stackoverflow.com/questions/14469522/stop-an-r-program-without-error
send_stop <- function(...) {
  args <- list(...)
  if (length(args) != 0) {
    send_error(...)
  }
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}
