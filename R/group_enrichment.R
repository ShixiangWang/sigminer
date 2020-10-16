#' General Group Enrichment Analysis
#'
#' This function takes a `data.frame` as input, compares proportion of positive
#' cases or mean measure in one subgroup and the remaining samples.
#'
#' @param df a `data.frame`.
#' @param grp_vars character vector specifying group variables to split samples
#' into subgroups (at least 2 subgroups, otherwise this variable will be skipped).
#' @param enrich_vars character vector specifying measure variables to be compared.
#' If variable is not numeric, only binary cases are accepted in the form of
#' `TRUE/FALSE` or `P/N` (P for positive cases and N for negative cases).
#' @param cross logical, default is `TRUE`, combine all situations provided by
#' `grp_vars` and `enrich_vars`. For examples, `c('A', 'B')` and `c('C', 'D')`
#' will construct 4 combinations(i.e. "AC", "AD", "BC" and "BD"). A variable can
#' not be in both `grp_vars` and `enrich_vars`, such cases will be automatically
#' drop. If `FALSE`, use pairwise combinations, see section "examples" for use cases.
#' @param co_method test method for continuous variable, default is 't.test'.
#'
#' @return a `data.table` with following columns:
#' - `grp_var`: group variable name.
#' - `enrich_var`: enrich variable (variable to be compared) name.
#' - `grp1`: the first group name, should be a member in `grp_var` column.
#' - `grp2`: the remaining samples, marked as 'Rest'.
#' - `grp1_size`: sample size for `grp1`.
#' - `grp1_pos_measure`: for binary variable, it stores the proportion of
#' positive cases in `grp1`; for continuous variable, it stores mean value.
#' - `grp2_size`: sample size for `grp2`.
#' - `grp2_pos_measure`: same as `grp1_pos_measure` but for `grp2`.
#' - `measure_observed`: for binary variable, it stores odds ratio;
#' for continuous variable, it stores scaled mean ratio.
#' - `measure_tested`: only for binary variable, it stores
#' estimated odds ratio and its 95% CI from `fisher.test()`.
#' - `p_value`: for binary variable, it stores p value from `fisher.test()`;
#' for continuous variable, it stores value from `wilcox.test()` or `t.test()`.
#' - `type`: one of "binary" and "continuous".
#' - `method`: one of "fish.test", "wilcox.test" and "t.test".
#' @export
#' @seealso [show_group_enrichment]
#'
#' @examples
#' set.seed(1234)
#' df <- dplyr::tibble(
#'   g1 = factor(abs(round(rnorm(99, 0, 1)))),
#'   g2 = rep(LETTERS[1:4], c(50, 40, 8, 1)),
#'   e1 = sample(c("P", "N"), 99, replace = TRUE),
#'   e2 = rnorm(99)
#' )
#'
#' print(str(df))
#' print(head(df))
#'
#' # Compare g1:e1, g1:e2, g2:e1 and g2:e2
#' x1 <- group_enrichment(df, grp_vars = c("g1", "g2"), enrich_vars = c("e1", "e2"))
#' x1
#'
#' # Only compare g1:e1, g2:e2
#' x2 <- group_enrichment(df, grp_vars = c("g1", "g2"),
#'                        enrich_vars = c("e1", "e2"),
#'                        co_method = "wilcox.test",
#'                        cross = FALSE)
#' x2
#'
#' # Visualization
#' p1 <- show_group_enrichment(x1, fill_by_p_value = TRUE)
#' p1
#' p2 <- show_group_enrichment(x1, fill_by_p_value = FALSE)
#' p2
#' p3 <- show_group_enrichment(x1, return_list = TRUE)
#' p3
#' @testexamples
#' expect_is(x1, "data.frame")
#' expect_is(x2, "data.frame")
#' expect_is(p1, "ggplot")
#' expect_is(p2, "ggplot")
#' expect_is(p3, "list")
group_enrichment <- function(df, grp_vars = NULL, enrich_vars = NULL,
                             cross = TRUE,
                             co_method = c("t.test", "wilcox.test")) {
  stopifnot(
    is.character(grp_vars), is.character(enrich_vars),
    is.data.frame(df),
    all(c(grp_vars, enrich_vars) %in% colnames(df))
  )
  co_method <- match.arg(co_method)

  if (cross) {
    comb_df <- expand.grid(grp_vars, enrich_vars) %>%
      t() %>%
      as.data.frame()
  } else {
    comb_df <- rbind(
      matrix(grp_vars, nrow = 1),
      matrix(enrich_vars, nrow = 1)
    ) %>%
      as.data.frame()
  }

  # Remove variables in both grp_vars and enrich_vars
  comb_df <- comb_df[, sapply(comb_df, function(x) length(unique(x)) != 1),
    drop = FALSE
  ]
  purrr::map_df(comb_df, ~ enrich_one(.x[1], .x[2], df, co_method))
}

enrich_one <- function(x, y, df, method = "t.test") {
  # x is group var
  # y is enrich var

  # Check the input
  df[[x]] <- as.character(df[[x]])
  grps <- na.omit(unique(df[[x]]))

  message("Handing pair ", x, ":", y)
  message("============================")

  if (length(grps) < 2) {
    message("Group variable contains <2 groups, skipping...")
    return(data.table::data.table())
  }

  if (is.numeric(df[[y]])) {
    # t.test or wilcox.test
    cmp <- purrr::map_df(grps, function(z, x, y, df, method = "t.test") {
      message("Handing pair ", x, ":", y, " - group: ", z)
      df$Group <- data.table::fifelse(df[[x]] %in% z, z, "Other", na = "Other")
      df$Group <- factor(x = df$Group, levels = c(z, "Other"))
      df$Cval <- df[[y]]

      # Only keep rows with valid values
      df <- df[, c("Group", "Cval")][is.finite(df$Cval), ]

      .test <- if (method == "t.test") stats::t.test else stats::wilcox.test
      cmp <- tryCatch(
        .test(Cval ~ Group, data = df),
        error = function(e) {
          message("  An error occur when testing, it will be skipped. Info:")
          message("    ", e$message)
          list(p.value = NA)
        }
      )

      data_range <- range(df$Cval)
      grp_sum <- df %>%
        dplyr::group_by(.data$Group) %>%
        dplyr::summarise(
          n = dplyr::n(),
          measure = mean(.data$Cval),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          # Scale all measures to range 0-1
          measure_scaled = (.data$measure - data_range[1]) / diff(data_range)
        )

      cmp.tbl <- data.table::data.table(
        grp_var = x,
        enrich_var = y,
        grp1 = z,
        grp2 = "Rest",
        grp1_size = grp_sum$n[1],
        grp1_pos_measure = grp_sum$measure[1],
        grp2_size = grp_sum$n[2],
        grp2_pos_measure = grp_sum$measure[2],
        measure_observed = grp_sum$measure_scaled[1] / grp_sum$measure_scaled[2],
        measure_tested = NA,
        p_value = cmp$p.value,
        type = "continuous",
        method = method
      )
      cmp.tbl
    },
    x = x,
    y = y,
    df = df,
    method = method
    )
  } else {
    # fisher test

    # P and N must in df[[y]]
    if (is.logical(df[[y]])) {
      df$Ctype <- data.table::fifelse(df[[y]], "P", "N", na = "N")
    } else if (is.character(df[[y]])) {
      if (!all(na.omit(df[[y]]) %in% c("P", "N"))) {
        stop("When input is string vector, just 'P' and 'N' are supported!")
      }
    } else {
      stop("Unsupported input for non-numeric variable, convert it to P/N or TRUE/FALSE!")
    }

    cmp <- purrr::map_df(grps, function(z, x, y, df) {
      message("Handing pair ", x, ":", y, " - group: ", z)

      df$Group <- data.table::fifelse(df[[x]] %in% z, z, "Other", na = "Other")
      df$Group <- factor(x = df$Group, levels = c(z, "Other"))
      df$Ctype <- factor(x = df[[y]], levels = c("P", "N"))
      df.tbl <- with(df, table(Group, Ctype))
      df.tbl <- df.tbl[c(z, "Other"), c("P", "N")]

      cmp <- tryCatch(
        fisher.test(df.tbl),
        error = function(e) {
          message("  An error occur when testing, it will be skipped. Info:")
          message("    ", e$message)
          list(
            p.value = NA,
            conf.int = c(NA, NA)
          )
        }
      )

      grp_size <- as.numeric(rowSums(df.tbl))
      grp_frac <- t(apply(df.tbl, 1, function(x) x / sum(x)))

      cmp.tbl <- data.table::data.table(
        grp_var = x,
        enrich_var = y,
        grp1 = z,
        grp2 = "Rest",
        grp1_size = grp_size[1],
        grp1_pos_measure = grp_frac[1, 1],
        grp2_size = grp_size[2],
        grp2_pos_measure = grp_frac[2, 1],
        measure_observed = grp_frac[1, 1] / grp_frac[2, 1],
        measure_tested = paste0(
          round(cmp$estimate, 3),
          " (", paste(round(cmp$conf.int, 3),
            collapse = ","
          ), ")"
        ),
        p_value = cmp$p.value,
        type = "binary",
        method = "fisher.test"
      )
      cmp.tbl
    },
    x = x,
    y = y,
    df = df
    )
  }

  return(cmp)
}
