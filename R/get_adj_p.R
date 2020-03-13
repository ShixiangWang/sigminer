# Get adjust p values -----------------------------------------------------

#' Get Adjust P Values from Group Comparison
#'
#' Setting `aes(label=..p.adj..)` in [ggpubr::compare_means()] does not
#' show adjust p values. The returned result of this function can be combined with [ggpubr::stat_pvalue_manual()] to fix
#' this problem.
#'
#' More info see [ggpubr::compare_means()], [ggpubr::stat_compare_means()] and [stats::p.adjust()].
#' @param data a `data.frame` containing column for groups and column for comparison.
#' @param .col column name for comparison.
#' @param .grp column name for groups.
#' @param p.digits how many significant digits are to be used.
#' @param comparisons Default is `NULL`, use all combination in group column.
#' It can be a list of length-2 vectors. The entries in the vector are either
#' the names of 2 values on the x-axis or the 2 integers that correspond to the
#' index of the groups of interest, to be compared.
#' @param method a character string indicating which method to be used for comparing means.
#' It can be 't.test', 'wilcox.test' etc..
#' @param p.adjust.method correction method, default is 'fdr'. Run `p.adjust.methods` to
#' see all available options.
#' @param ... other arguments passed to [ggpubr::compare_means()]
#' @source https://github.com/kassambara/ggpubr/issues/143
#' @return a `data.frame` containing comparison result
#' @export
#'
#' @examples
#' library(ggpubr)
#' # T-test
#' stat.test <- compare_means(
#'   len ~ dose,
#'   data = ToothGrowth,
#'   method = "t.test",
#'   p.adjust.method = "fdr"
#' )
#' stat.test
#' # Create a simple box plot
#' p <- ggboxplot(ToothGrowth, x = "dose", y = "len")
#' p
#'
#' # Add p values
#' my_comparisons <- list(c("0.5", "1"), c("1", "2"), c("0.5", "2"))
#' p + stat_compare_means(method = "t.test", comparisons = my_comparisons)
#'
#' # Try adding adjust p values
#' # proposed by author of ggpubr
#' # however it does not work
#' p + stat_compare_means(aes(label = ..p.adj..), method = "t.test", comparisons = my_comparisons)
#'
#' # Solution:
#' # calculate adjust p values and their location
#' # then use stat_pvalue_manual() function
#' p_adj <- get_adj_p(ToothGrowth, .col = "len", .grp = "dose")
#' p_adj
#' p + stat_pvalue_manual(p_adj, label = "p.adj")
#'
#' # Show selected comparisons
#' # Of note, p value is ajusted
#' # for three comparisons, but only
#' # two are showed in figure
#' p_adj <- get_adj_p(ToothGrowth,
#'   .col = "len", .grp = "dose",
#'   comparisons = list(c("0.5", "1"), c("1", "2"))
#' )
#' p + stat_pvalue_manual(p_adj, label = "p.adj")
#' @testexamples
#' expect_s3_class(p, "ggplot")
#' expect_s3_class(p_adj, "tbl")
get_adj_p <- function(data, .col, .grp = "Sample", comparisons = NULL,
                      method = "wilcox.test", p.adjust.method = "fdr", p.digits = 3L, ...) {
  # Compute p-values
  comparison.formula <- paste0(.col, "~", .grp) %>%
    as.formula()
  pvalues <- ggpubr::compare_means(
    formula = comparison.formula, data = data,
    method = method,
    p.adjust.method = p.adjust.method,
    ...
  )

  # If a comparison list is provided, extract the comparisons of interest for plotting
  if (!is.null(comparisons)) {
    pvalues <- purrr::map_df(comparisons, ~ pvalues %>% dplyr::filter(group1 == .x[1] & group2 == .x[2]))
  }

  # P-value y coordinates
  y.max <- data %>%
    dplyr::pull(.col) %>%
    max(na.rm = TRUE)
  p.value.y.coord <- rep(y.max, nrow(pvalues))

  step.increase <- (1:nrow(pvalues)) * (y.max / 10)
  p.value.y.coord <- p.value.y.coord + step.increase

  pvalues <- pvalues %>%
    dplyr::mutate(
      y.position = p.value.y.coord,
      p.adj = format.pval(.data$p.adj, digits = p.digits)
    )

  pvalues
}
