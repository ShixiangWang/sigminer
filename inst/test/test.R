library(maftools)
?read.maf
laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml <- read.maf(maf = laml.maf)
laml.tnm <- trinucleotideMatrix(maf = laml, ref_genome = 'BSgenome.Hsapiens.UCSC.hg19',
                                prefix = 'chr', add = TRUE, useSyn = TRUE)
laml.sign <- extractSignatures(mat = laml.tnm, n = 3, pConstant = 0.001)
laml.sign$signatures

compareSignatures(laml.sign)


load(system.file("extdata", "toy_copynumber_signature.RData",
                 package = "sigminer", mustWork = TRUE
))
# Load copy number prepare object
load(system.file("extdata", "toy_copynumber_prepare.RData",
                 package = "sigminer", mustWork = TRUE
))
params = get_tidy_parameter(cn_prepare$components)
show_sig_profile(sig, params = params, y_expand = 2)
show_sig_profile(sig, params = params, y_expand = 3, show_cv = T, params_label_angle = 70)



# Optimize signature exposure
load(system.file("extdata", "toy_copynumber_prepare.RData",
                 package = "sigminer", mustWork = TRUE
))
res <- sig_auto_extract(cn_prepare$nmf_matrix, result_prefix = "Test_copynumber", nrun = 1)
