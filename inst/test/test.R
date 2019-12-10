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

str(cn_prepare$features$segsize)


load(system.file("extdata", "toy_segTab.RData",
                 package = "sigminer", mustWork = TRUE
))

segTabs$start = as.character(segTabs$start)
segTabs$end = as.numeric(segTabs$end)

system.time(
  cn <- read_copynumber(segTabs,
                        seg_cols = c("chromosome", "start", "end", "segVal"),
                        genome_build = "hg19", complement = FALSE, verbose = TRUE
  )
)


# Load copy number object
load(system.file("extdata", "toy_copynumber.RData",
                 package = "sigminer", mustWork = TRUE
))
# Prepare copy number signature analysis
system.time(
  cn_prepare3 <- sig_derive(cn, method = "M")
)
system.time(
  cn_prepare2 <- sig_derive(cn, method = "W")
)
# 31s

str(cn_prepare$features$copynumber)

cc = getBPnum(cn_list, chrlen = sigminer::chromsize.hg19)
cc = getOscilation(cn_list)
cc2 = getOscilation(cn_list)

cc = getCentromereDistCounts(cn_list, centromeres.hg19, chromsize.hg19)
cc2 = getCentromereDistCounts2(cn_list, centromeres.hg19, chromsize.hg19)
debug(getCentromereDistCounts2)
