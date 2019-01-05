context("test-sig_analysis_series")

load(system.file("extdata", "example_cn_list.RData",
                 package = "sigminer", mustWork = TRUE))

segTabs = data.table::rbindlist(tcga_segTabs, idcol = "sample")


test_that("signature analysis series functions work", {
  skip_on_cran()

  # focus on copy number signatures for now
  cn = read_copynumber(segTabs,
                       seg_cols = c("chromosome", "start", "end", "segVal"),
                       genome_build = "hg19")
  expect_s4_class(cn, "CopyNumber")

  cn_prepare =  sig_prepare(cn, cores = 2)
  expect_identical(inherits(cn_prepare, "list"), TRUE)
  expect_identical(length(cn_prepare) == 3, TRUE)

  # estimate rank
  library(NMF)
  cn_estimate = sig_estimate(cn_prepare$nmf_matrix, cores = 2, nrun = 5, verbose = TRUE)
  expect_identical(inherits(cn_estimate, "list"), TRUE)
  expect_identical(length(cn_estimate) == 5, TRUE)

  # do NMF decomposition
  res = sig_extract(cn_prepare$nmf_matrix, 3, mode = "copynumber", nrun = 5, cores = 2)
  expect_s4_class(res$nmfObj, "NMFfitX1")

})

# clean
segTabs <- tcga_segTabs <- NULL
