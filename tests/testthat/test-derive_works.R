skip_on_cran()

test_that("sig_feed for CopyNumber works", {
  # Load copy number object
  load(system.file("extdata", "toy_copynumber.RData",
    package = "sigminer", mustWork = TRUE
  ))
  # Prepare copy number signature analysis
  cn_prepare <- sig_feed(cn)
  expect_type(cn_prepare, "list")
  cn_prepare <- sig_feed(cn, method = "W")
  expect_type(cn_prepare, "list")
  cn_prepare2 <- sig_feed(cn, cores = 2, method = "W")
  expect_type(cn_prepare2, "list")
  # The result is reproducible
  cn_prepare <- sig_feed(cn, cores = 2, method = "W")
  expect_true(all.equal(cn_prepare, cn_prepare2))
  cn_prepare <- sig_feed(cn, keep_only_matrix = TRUE, method = "W")
  expect_type(cn_prepare, "integer")
})


test_that("sig_feed for MAF works", {
  # Prepare mutational signature analysis
  laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
  laml <- read_maf(maf = laml.maf)
  library(BSgenome.Hsapiens.UCSC.hg19)
  expect_warning(
    mt_prepare <- sig_feed(
      laml,
      ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
      prefix = "chr", add = TRUE, useSyn = TRUE
    )
  )
  expect_type(mt_prepare, "list")
  expect_warning(
    mt_prepare <- sig_feed(
      laml,
      ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
      prefix = "chr", add = TRUE, useSyn = TRUE, keep_only_matrix = TRUE
    )
  )
  expect_type(mt_prepare, "double")
})
