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

# Load copy number prepare object
load(system.file("extdata", "toy_copynumber_prepare.RData",
                 package = "sigminer", mustWork = TRUE
))
# Extract copy number signatures
res <- sig_extract(cn_prepare$nmf_matrix, 2, nrun = 1)
res

load("/Users/wsx/biosoft/SignatureAnalzyer.052418/INPUT_SignatureAnalyzer/lego96.PAN.SNV.091217.RData")
ttype <- sapply(colnames(lego96.SNV),function(x) strsplit(x,"__")[[1]][1])
lego96.demo <- lego96.SNV[,ttype=="Biliary_AdenoCA"]

tt = sig_auto_extract(t(lego96.demo), result_prefix = "Test", K0 = 25, tol = 1e-7, n_iter = 2e5, cores = 4)
# test copy number
tt = sig_auto_extract(cn_prepare$nmf_matrix, result_prefix = "Test_copynumber", K0 = 25, tol = 1e-7, n_iter = 2e5, cores = 8)
