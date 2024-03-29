# Generated by roxytest: do not edit by hand!

# File R/signature_obj_operation.R: @testexamples

test_that("Function sig_names() @ L21", {
  
  ## Operate signature names
  load(system.file("extdata", "toy_mutational_signature.RData",
    package = "sigminer", mustWork = TRUE
  ))
  sig_names(sig2)
  cc <- sig_modify_names(sig2, new_names = c("Sig2", "Sig1", "Sig3"))
  sig_names(cc)
  
  # The older names are stored in tags.
  print(attr(cc, "tag"))
  expect_is(cc, "Signature")
})

