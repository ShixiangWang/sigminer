skip_on_cran()

test_that("sig_extract works", {
  # Load copy number prepare object
  load(system.file("extdata", "toy_copynumber_prepare.RData",
    package = "sigminer", mustWork = TRUE
  ))
  # Extract copy number signatures
  library(NMF)
  res <- sig_extract(cn_prepare$nmf_matrix, 2, nrun = 1)
  expect_s3_class(res, "Signature")
})
