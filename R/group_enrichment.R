#
# df <- data.frame(
#   g1 = factor(abs(round(rnorm(99, 0, 1)))),
#   g2 = factor(abs(round(rnorm(99, 0, 1)))),
#   e1 = sample(c("P", "N"), 99, replace = TRUE),
#   e2 = rnorm(99)
# )
group_enrichment <- function(df, grp_vars = NULL, enrich_vars = NULL,
                             co_method = c("t.test", "wilcox.test")) {
  # grp_vars 指代离散的列，只要要两组
  # enrich_vars 指代要对比的列，可以是离散值（只能二分类，用 P/TRUE 表示正分类，N/FALSE表示负分类）或连续值
  # 如果是离散值，使用 fisher 检验，如果是连续值，使用 t/wilcox 检验
  stopifnot(is.character(grp_vars), is.character(enrich_vars),
            is.data.frame(df),
            all(c(grp_vars, enrich_vars) %in% colnames(df)))
  co_method <- match.arg(co_method)

  comb_df <- expand.grid(grp_vars, enrich_vars) %>% t() %>% as.data.frame()
  # Remove variables in both grp_vars and enrich_vars
  comb_df <- comb_df[, sapply(comb_df, function(x) length(unique(x)) != 1),
                     drop = FALSE]
  # purrr::map(comb_df, ~.x)

}

enrich_one <- function(x, y, df, method = "t.test") {
  # x is group var
  # y is enrich var

  # Check the input

  df[[x]] <- as.character(df[[x]])
  grps <- na.omit(unique(df[[x]]))

  if (is.numeric(df[[y]])) {
    # t.test or wilcox.test

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

    cmp = purrr::map_df(grps, function(z, x, y, df){
      df$Group = data.table::fifelse(df[[x]] %in% z, z, "Other", na = "Other")
      df$Group = factor(x = df$Group, levels = c(z, "Other"))
      df$Ctype = factor(x = df[[y]], levels = c("P", "N"))
      df.tbl = with(df, table(Group, Ctype))
      df.tbl = df.tbl[c(z, "Other"), c("P", "N")]
      cmp = fisher.test(df.tbl)

      grp_size <- as.numeric(rowSums(df.tbl))
      grp_frac <- t(apply(df.tbl, 1, function(x) x / sum(x)))

      cmp.tbl = data.table::data.table(
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
                      collapse = ","), ")"),
        p_value = cmp$p.value)
      cmp.tbl
    },
    x = x,
    y = y,
    df = df)

  }
}
