# Generated by roxytest: do not edit by hand!

# File R/get_adj_p.R: @testexamples

test_that("Function get_adj_p() @ L69", {
  
  library(ggpubr)
  # T-test
  stat.test <- compare_means(
    len ~ dose,
    data = ToothGrowth,
    method = "t.test",
    p.adjust.method = "fdr"
  )
  stat.test
  # Create a simple box plot
  p <- ggboxplot(ToothGrowth, x = "dose", y = "len")
  p
  
  # Add p values
  my_comparisons <- list(c("0.5", "1"), c("1", "2"), c("0.5", "2"))
  p + stat_compare_means(method = "t.test", comparisons = my_comparisons)
  
  # Try adding adjust p values
  # proposed by author of ggpubr
  # however it does not work
  p + stat_compare_means(aes(label = ..p.adj..), method = "t.test", comparisons = my_comparisons)
  
  # Solution:
  # calculate adjust p values and their location
  # then use stat_pvalue_manual() function
  p_adj <- get_adj_p(ToothGrowth, .col = "len", .grp = "dose")
  p_adj
  p + stat_pvalue_manual(p_adj, label = "p.adj")
  
  # Show selected comparisons
  # Of note, p value is ajusted
  # for three comparisons, but only
  # two are showed in figure
  p_adj <- get_adj_p(ToothGrowth,
    .col = "len", .grp = "dose",
    comparisons = list(c("0.5", "1"), c("1", "2"))
  )
  p + stat_pvalue_manual(p_adj, label = "p.adj")
  expect_s3_class(p, "ggplot")
  expect_s3_class(p_adj, "tbl")
})

