load(system.file("extdata", "toy_copynumber_prepare.RData",
                 package = "sigminer", mustWork = TRUE
))
library(NMF)
cn_estimate <- sig_estimate(cn_prepare$nmf_matrix,
                            cores = 1, nrun = 5,
                            verbose = TRUE
)


# Rewriting show cn components --------------------------------------------

# Load copy number prepare object
load(system.file("extdata", "toy_copynumber_prepare.RData",
                 package = "sigminer", mustWork = TRUE
))
show_cn_components(cn_prepare$features, cn_prepare$components)

data = cn_prepare$parameters

