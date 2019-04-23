context("test-sig_analysis_series")

load(system.file("extdata", "example_cn_list.RData",
  package = "sigminer", mustWork = TRUE
))

segTabs <- data.table::rbindlist(tcga_segTabs, idcol = "sample")


test_that("signature analysis series and result visualization functions work", {
  skip_on_cran() # can not pass windows when using more than 1 cpu

  # focus on copy number signatures for now
  cn <- read_copynumber(segTabs,
    seg_cols = c("chromosome", "start", "end", "segVal"),
    genome_build = "hg19"
  )
  expect_s4_class(cn, "CopyNumber")
  p1 <- draw_cn_distribution(cn)
  expect_true(inherits(p1, "ggplot"))

  cn_prepare <- sig_prepare(cn, cores = 2)
  expect_true(inherits(cn_prepare, "list"))
  expect_true(length(cn_prepare) == 3)

  p2 <- draw_cn_features(cn_prepare$features)
  expect_true(inherits(p2, "ggplot"))

  p3 <- draw_cn_components(cn_prepare$features, cn_prepare$components)
  expect_true(inherits(p3, "ggplot"))

  # estimate rank
  library(NMF)
  cn_estimate <- sig_estimate(cn_prepare$nmf_matrix, cores = 2, nrun = 5, verbose = TRUE)
  expect_true(inherits(cn_estimate, "list"))
  expect_true(length(cn_estimate) == 5)

  # do NMF decomposition
  res <- sig_extract(cn_prepare$nmf_matrix, 3, mode = "copynumber", nrun = 5, cores = 2)
  expect_s4_class(res$nmfObj, "NMFfitX1")

  p1 <- draw_sig_profile(res$nmfObj)
  expect_true(inherits(p1, "ggplot"))

  p2 <- draw_sig_activity(res$nmfObj)
  expect_true(inherits(p2, "ggplot"))

  sig_activity <- sig_get_activity(res$nmfObj)
  expect_true(is.list(sig_activity))

  sig_cor <- sig_get_correlation(sig_activity)
  expect_true(is.list(sig_cor))

  draw_sig_corrplot(sig_cor)

  subtypes <- sig_assign_samples(res$nmfObj)
  expect_true(inherits(subtypes, "data.frame"))
  set.seed(1234)
  subtypes$new_group <- sample(c("1", "2", "3", "4"), size = nrow(subtypes), replace = TRUE)

  subtypes.sum <- sig_summarize_subtypes(subtypes[, -1],
    col_subtype = "nmf_subtypes",
    cols_to_summary = colnames(subtypes[, -1])[-1],
    type = c("co", "ca"), verbose = TRUE
  )
  expect_true(is.list(subtypes.sum))

  p_list <- draw_subtypes_comparison(subtypes.sum)
  expect_true(is.list(p_list))
})

# clean
segTabs <- tcga_segTabs <- NULL
