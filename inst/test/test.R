# Load copy number prepare object
load(system.file("extdata", "toy_copynumber_prepare.RData",
                 package = "sigminer", mustWork = TRUE
))

sig = sig_auto_extract(cn_prepare$nmf_matrix)
