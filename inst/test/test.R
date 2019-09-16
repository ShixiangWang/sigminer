# Load copy number object
load(system.file("extdata", "toy_copynumber.RData",
                 package = "sigminer", mustWork = TRUE
))
# Prepare copy number signature analysis
cn_prepare3 <- derive(cn)
save(cn_prepare, file = "inst/extdata/toy_copynumber_prepare.RData")


cn_prepare <- derive(cn, cores = 2)

# Load copy number prepare object
load(system.file("extdata", "toy_copynumber_prepare.RData",
                 package = "sigminer", mustWork = TRUE
))

sig = sig_auto_extract(cn_prepare$nmf_matrix)
save(sig, file = "inst/extdata/toy_copynumber_signature.RData")


laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml <- read_maf(maf = laml.maf)
library(BSgenome.Hsapiens.UCSC.hg19)
mt_prepare <- derive(
  laml,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
  prefix = "chr", add = TRUE, useSyn = TRUE
)

sig2 = sig_auto_extract(mt_prepare$nmf_matrix)
save(sig2, file ="inst/extdata/toy_mutational_signature.RData")


#-----------

show_sig_profile(sig)
show_sig_profile(sig, normalize = "column")
show_sig_profile(sig, normalize = "raw")
show_sig_profile(sig$Raw$W)
show_sig_profile(sig$Raw$W, normalize = "raw")


show_sig_profile(sig2, mode = "mutation")


zz = show_cn_components(cn_prepare$features, cn_prepare$components)
zz$parameters
show_sig_profile(sig, params = cn_prepare$parameters, y_expand = 2)
