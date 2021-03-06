# Generated by roxytest: Do not edit by hand!

context("File R/read_maf.R: @testexamples")

test_that("Function read_maf() @ L25", {
  
  laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools", mustWork = TRUE)
  if (!require("R.utils")) {
    message("Please install 'R.utils' package firstly")
  } else {
    laml <- read_maf(maf = laml.maf)
    laml
  }
  expect_s4_class(laml, "MAF")
})

