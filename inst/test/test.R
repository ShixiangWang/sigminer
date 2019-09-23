load(system.file("extdata", "toy_copynumber_prepare.RData",
                 package = "sigminer", mustWork = TRUE
))
library(NMF)
cn_estimate <- sig_estimate(cn_prepare$nmf_matrix,
                            cores = 1, nrun = 5,
                            verbose = TRUE
)
