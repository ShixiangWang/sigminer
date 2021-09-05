#' Get Comparison Result between Signature Groups
#'
#' Compare genotypes/phenotypes based on signature groups (samples are assigned to
#' several groups). For categorical
#' type, calculate fisher p value (using [stats::fisher.test]) and count table.
#' In larger than 2 by 2 tables, compute p-values by Monte Carlo simulation.
#' For continuous type, calculate anova p value (using [stats::aov]),
#' summary table and Tukey Honest significant difference (using [stats::TukeyHSD]).
#' The result of this function can be plotted by [show_group_comparison()].
#'
#' @param data a `data.frame` containing signature groups and genotypes/phenotypes
#' (including categorical and continuous type data) want to analyze. User need to
#' construct this `data.frame` by him/herself.
#' @param col_group column name of signature groups.
#' @param cols_to_compare column names of genotypes/phenotypes want to summarize based on groups.
#' @param type a characater vector with length same as `cols_to_compare`,
#' 'ca' for categorical type and 'co' for continuous type.
#' @param NAs default is `NA`, filter `NA`s for categorical columns.
#' Otherwise a value (either length 1 or length same as `cols_to_compare`) fill `NA`s.
#' @param verbose if `TRUE`, print extra information.
#' @author Shixiang Wang <w_shixiang@163.com>
#' @return a `list` contains data, summary, p value etc..
#' @export
#' @examples
#' \donttest{
#' load(system.file("extdata", "toy_copynumber_signature_by_W.RData",
#'   package = "sigminer", mustWork = TRUE
#' ))
#'
#' # Assign samples to clusters
#' groups <- get_groups(sig, method = "k-means")
#'
#' set.seed(1234)
#'
#' groups$prob <- rnorm(10)
#' groups$new_group <- sample(c("1", "2", "3", "4", NA), size = nrow(groups), replace = TRUE)
#'
#' # Compare groups (filter NAs for categorical coloumns)
#' groups.cmp <- get_group_comparison(groups[, -1],
#'   col_group = "group",
#'   cols_to_compare = c("prob", "new_group"),
#'   type = c("co", "ca"), verbose = TRUE
#' )
#'
#' # Compare groups (Set NAs of categorical columns to 'Rest')
#' groups.cmp2 <- get_group_comparison(groups[, -1],
#'   col_group = "group",
#'   cols_to_compare = c("prob", "new_group"),
#'   type = c("co", "ca"), NAs = "Rest", verbose = TRUE
#' )
#' }
get_group_comparison <- function(data, col_group, cols_to_compare,
                                 type = "ca", NAs = NA, verbose = FALSE) {
  if (!all(type %in% c("ca", "co"))) {
    stop("all elements in 'type' must be 'ca' for 'categorical' variable and 'co' for 'continuous' variable.")
  }
  if (!is.data.frame(data)) stop("'data' must be a data.frame object.")
  data <- data.table::as.data.table(data)

  # subset
  data <- data[, c(col_group, cols_to_compare), with = FALSE]
  colnames(data)[1] <- "group"
  data <- data[!is.na(data[["group"]])]
  data[["group"]] <- factor(data[["group"]])

  if (length(levels(data[["group"]])) < 2) {
    stop("At least two factors in ", col_group, " column")
  }

  do_summary <- function(col, type = c("ca", "co"),
                         verbose = FALSE, na = NA) {
    type <- match.arg(type)

    df <- data[, c("group", col), with = FALSE]
    if (!is.na(na) & type == "ca") {
      df[[col]] <- ifelse(is.na(df[[col]]), na, df[[col]])
    } else {
      df <- df[!is.na(df[[col]])]
    }

    if (type == "ca") {
      if (verbose) message("Treat ", col, " as categorical variable.")

      table_df <- table(df[["group"]], df[[col]])

      table_p <- tryCatch(
        {
          if (any(dim(table_df) > 2)) {
            # compute p-values by Monte Carlo simulation, in larger than 2 by 2 tables.
            test <- fisher.test(table_df, simulate.p.value = TRUE)
          } else {
            test <- fisher.test(table_df)
          }
          test[["p.value"]]
        },
        error = function(e) {
          NA
        }
      )

      list(data = df, table = table_df, p_value = table_p, type = "categorical", extra = NA)
    } else {
      if (verbose) message("Treat ", col, " as continuous variable.")

      table_df <- summary(df)

      fit <- tryCatch(
        {
          stats::aov(as.formula(paste0(col, " ~ group")), data = df)
        },
        error = function(e) {
          NA
        }
      )

      if (inherits(fit, "aov")) {
        p_value <- summary(fit)[[1]][["Pr(>F)"]][1] # get anova p value
        extra <- stats::TukeyHSD(fit)[[1]]
      } else {
        p_value <- NA
        extra <- NA
      }
      list(data = df, table = table_df, p_value = p_value, type = "continuous", extra = extra)
    }
  }

  res <- Map(do_summary, cols_to_compare, type, verbose, NAs)
  names(res) <- cols_to_compare
  res
}
