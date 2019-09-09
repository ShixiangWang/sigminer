load(system.file("extdata", "toy_copynumber.RData",
                 package = "sigminer", mustWork = TRUE
))
# Prepare copy number signature analysis
cn_prepare <- derive(cn, cores = 1)
cn_prepare <- derive.CopyNumber(cn, cores = 2)
cn_prepare <- derive(cn, cores = 4)

registerDoFuture()
plan(multiprocess, workers = 2)
ext <- foreach(file = c("324-fds", "yes")) %dopar% isValidAndUnreserved(file)
