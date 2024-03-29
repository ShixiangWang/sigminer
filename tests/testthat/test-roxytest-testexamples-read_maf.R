# Generated by roxytest: do not edit by hand!

# File R/read_maf.R: @testexamples

test_that("Function read_maf() @ L35", {
  
  
  laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools", mustWork = TRUE)
  if (!require("R.utils")) {
    message("Please install 'R.utils' package firstly")
  } else {
    laml <- read_maf(maf = laml.maf)
    laml
  
    laml_mini <- laml@data[, list(
      Tumor_Sample_Barcode, Chromosome,
      Start_Position, End_Position,
      Reference_Allele, Tumor_Seq_Allele2
    )]
    laml2 <- read_maf_minimal(laml_mini)
    laml2
  }
  
  expect_s4_class(laml, "MAF")
})

