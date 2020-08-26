
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

file_name <- function(filepath, must_chop = NULL) {
  y <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(filepath))
  if (!is.null(must_chop)) {
    y <- sub(must_chop, "", y)
  }
  y
}

# From https://gist.github.com/mbannert/e9fcfa86de3b06068c83
rgb2hex <- function(r, g, b) grDevices::rgb(r, g, b, maxColorValue = 255)
col2hex <- function(col, alpha) grDevices::rgb(t(grDevices::col2rgb(col)), alpha = alpha, maxColorValue = 255)


#' Set Color Style for Plotting
#'
#' @param style one of 'default' and 'cosmic'.
#' @param mode only used when the `style` is 'cosmic', can be one of
#' "SBS", "copynumber", "DBS", "ID".
#' @param method used to set a more custom palette for different methods.
#'
#' @return color values.
#' @export
#'
#' @examples
#' use_color_style("default")
#' use_color_style("cosmic")
use_color_style <- function(style, mode = c("SBS", "copynumber", "DBS", "ID"), method = "Wang") {
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

    if (mode == "copynumber" & startsWith(method, "T")) {
      colors <- list(
        c(226, 226, 239), c(182, 182, 216), c(134, 131, 189), c(98, 64, 155),
        c(253, 202, 181), c(252, 138, 106), c(241, 68, 50), c(188, 25, 26),
        c(208, 225, 242), c(148, 196, 223), c(74, 152, 201), c(23, 100, 171),
        c(253, 222, 161), c(253, 190, 111), c(255, 128, 2), c(255, 80, 2),
        c(176, 221, 139), c(100, 200, 100), c(54, 161, 46)
      )
    }

    palette <- sapply(colors, FUN = function(x) rgb2hex(x[1], x[2], x[3])) %>% as.character()

    if (mode %in% c("copynumber", "SBS")) {
      palette <- c(palette, sapply(c("purple", "brown", "orange"), FUN = col2hex) %>% as.character())
    }
  }
  palette
}

# https://en.wikipedia.org/wiki/Help:Distinguishable_colors
letter_colors <- c("#FF5005", "#FFFF80", "#990000", "#740AFF", "#E0FF66",
                   "#00998F", "#5EF1F2", "#FF0010", "#426600", "#FFA8BB",
                   "#FFA405", "#003380", "#C20088", "#9DCC00", "#8F7C00",
                   "#94FFB5", "#94FFB5", "#FFCC99", "#2BCE48", "#005C31",
                   "#191919", "#4C005C", "#993F00", "#0075DC", "#F0A3FF",
                   "#000000", "#FFFFFF")


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
