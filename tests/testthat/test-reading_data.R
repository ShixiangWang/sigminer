context("test-reading_data")

extdata_dir <- system.file("extdata", package = "sigminer", mustWork = TRUE)
laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools", mustWork = TRUE)

test_that("read copy number or maf work", {
  skip_on_cran()
  expect_message(hello())

  cp <- read_copynumber(extdata_dir, pattern = "txt", genome_build = "hg19")
  expect_s4_class(cp, "CopyNumber")

  laml <- read_maf(maf = laml.maf)
  expect_s4_class(laml, "MAF")

  vt <- read_variation(cp, laml)
  expect_s4_class(vt, "GenomicVariation")
})

extdata_dir <- laml.maf <- cp <- laml <- vt <- NULL
