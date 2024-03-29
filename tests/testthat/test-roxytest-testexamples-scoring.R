# Generated by roxytest: do not edit by hand!

# File R/scoring.R: @testexamples

test_that("Function scoring() @ L69", {
  
  # Load copy number object
  load(system.file("extdata", "toy_copynumber.RData",
    package = "sigminer", mustWork = TRUE
  ))
  
  d <- scoring(cn)
  d
  
  d2 <- scoring(cn, TD_cn_cutoff = 4L)
  d2
  expect_s3_class(d, "data.frame")
  expect_s3_class(d2, "data.frame")
})

