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

  # Auto-extract copy number signatures
  res <- sig_auto_extract(cn_prepare$nmf_matrix, result_prefix = "Test_copynumber", nrun = 1)
  expect_s3_class(res, "Signature")
  res <- sig_auto_extract(cn_prepare$nmf_matrix, result_prefix = "Test_copynumber", nrun = 1, method = "L1KL")
  expect_s3_class(res, "Signature")
  res <- sig_auto_extract(cn_prepare$nmf_matrix, result_prefix = "Test_copynumber", nrun = 1, method = "L2KL")
  expect_s3_class(res, "Signature")
})
