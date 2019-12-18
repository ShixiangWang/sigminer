# Load copy number object
load(system.file("extdata", "toy_copynumber.RData",
                 package = "sigminer", mustWork = TRUE
))

# Prepare copy number signature analysis
options(sigminer.sex = "male", sigminer.copynumber.max = 10)

cn_prepare <- sig_derive(cn, method = "W")
cn_prepare <- sig_derive(cn, method = "W", cores = 2)
