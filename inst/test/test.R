# Load copy number object
load(system.file("extdata", "toy_copynumber.RData",
                 package = "sigminer", mustWork = TRUE
))
# Prepare copy number signature analysis
system.time(
  cn_prepare <- derive(cn, cores = 4)
)
# 用户   系统   流逝
# 15.206  1.969 42.834

# Load copy number prepare object
load(system.file("extdata", "toy_copynumber_prepare.RData",
                 package = "sigminer", mustWork = TRUE
))
show_cn_components(cn_prepare$features, cn_prepare$components)

parameters(cn_prepare$components$segsize)
parameters(cn_prepare$components$bp10MB)

paras = parameters(cn_prepare$components$segsize)



get_tidy_parameter(cn_prepare$components)
get_tidy_parameter(cn_prepare$components$segsize)
