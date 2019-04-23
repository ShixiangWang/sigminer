context("test-prepare_series")

extdata_dir <- system.file("extdata", package = "sigminer", mustWork = TRUE)
laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools", mustWork = TRUE)

test_that("prepare_* series functions work", {
  skip_on_cran()

  cp <- read_copynumber(extdata_dir, pattern = "txt", genome_build = "hg19")
  expect_s4_class(cp, "CopyNumber")

  cn_prepare <- prepare_copynumber(cp)
  expect_true(inherits(cn_prepare, "list"))
  expect_true(length(cn_prepare) == 3)

  # laml = read_maf(maf = laml.maf)
  # expect_s4_class(laml, "MAF")

  # this function need BSgenome etc. package installed
  # laml.tnm <- prepare_maf(maf = laml, ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
  #                         prefix = 'chr', add = TRUE, useSyn = TRUE)
})

extdata_dir <- laml.maf <- NULL
