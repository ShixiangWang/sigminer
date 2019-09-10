load(system.file("extdata", "toy_copynumber.RData",
                 package = "sigminer", mustWork = TRUE
))
# Prepare copy number signature analysis
system.time(
  cn_prepare <- derive(cn, cores = 1)
)
system.time(
  cn_prepare <- derive.CopyNumber(cn, cores = 4)
)

cn_prepare <- derive.CopyNumber(cn, cores = 2)


registerDoFuture()
plan(multiprocess, workers = 2)
ext <- foreach(file = c("324-fds", "yes")) %dopar% isValidAndUnreserved(file)


# Load copy number prepare object
load(system.file("extdata", "toy_copynumber_prepare.RData",
                 package = "sigminer", mustWork = TRUE
))
# Extract copy number signatures
res <- sig_extract(cn_prepare$nmf_matrix, 2, nrun = 1)
res
