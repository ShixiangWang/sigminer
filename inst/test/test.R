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
zz1 = show_cn_components(cn_prepare$features, cn_prepare$components)

zz2 = show_cn_components2(cn_prepare$parameters, cn_prepare$features)
zz3 = show_cn_components2(cn_prepare$parameters, cn_prepare$features, show_weights = F)
data = cn_prepare$parameters

# 首先按照分布进行取样
# n_obs 用来计算权重
sample_size = 1000
set.seed(1234)
