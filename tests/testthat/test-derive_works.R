skip_on_cran()

test_that("derive for CopyNumber works", {
  # Load copy number object
  load(system.file("extdata", "toy_copynumber.RData",
    package = "sigminer", mustWork = TRUE
  ))
  # Prepare copy number signature analysis
  cn_prepare <- derive(cn)
  expect_type(cn_prepare, "list")
  cn_prepare <- derive(cn, keep_only_matrix = TRUE)
  expect_type(cn_prepare, "double")
})


test_that("derive for MAF works", {
  # Prepare mutational signature analysis
  laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
  laml <- read_maf(maf = laml.maf)
  library(BSgenome.Hsapiens.UCSC.hg19)
  expect_warning(
    mt_prepare <- derive(
      laml,
      ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
      prefix = "chr", add = TRUE, useSyn = TRUE
    )
  )
  expect_type(mt_prepare, "list")
  expect_warning(
    mt_prepare <- derive(
      laml,
      ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
      prefix = "chr", add = TRUE, useSyn = TRUE, keep_only_matrix = TRUE
    )
  )
  expect_type(mt_prepare, "double")
})
