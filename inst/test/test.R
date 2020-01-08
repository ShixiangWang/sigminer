# Load copy number object
load(system.file("extdata", "toy_copynumber.RData",
                 package = "sigminer", mustWork = TRUE
))

# Prepare copy number signature analysis
options(sigminer.sex = "male", sigminer.copynumber.max = 10)

cn_prepare <- sig_tally(cn, method = "W")
cn_prepare <- sig_tally(cn, method = "W", cores = 2)

df <- dplyr::tibble(
  start = c(1, 10, 20, 30, 40, 50, 60, 70, 80),
  end = c(9, 19, 29, 39, 49, 59, 69, 79, 89),
  segVal = c(2, 2, 1, 2, 3, 3, 2, 2, 2)
)

df <- dplyr::tibble(
  start = c(1, 10, 20, 30, 40, 50, 60, 70, 80),
  end = c(9, 19, 29, 39, 49, 59, 69, 79, 89),
  segVal = c(2, 2, 2, 2, 3, 3, 2, 2, 2)
)

df <- dplyr::tibble(
  start = c(1, 10, 20, 30, 40, 50, 60, 70, 80),
  end = c(9, 19, 29, 39, 49, 59, 69, 79, 89),
  segVal = c(1, 2, 1, 2, 3, 3, 2, 2, 2)
)

df <- dplyr::tibble(
  start = c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90),
  end = c(9, 19, 29, 39, 49, 59, 69, 79, 89, 99),
  segVal = c(2, 2, 1, 2, 3, 3, 2, 2, 2, 1)
)

df <- dplyr::tibble(
  start = c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90),
  end = c(9, 19, 29, 39, 49, 59, 69, 79, 89, 99),
  segVal = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2)
)

df <- dplyr::tibble(
  start = c(1, 10, 20, 30, 40, 50, 60, 70, 80, 90),
  end = c(9, 19, 29, 39, 49, 59, 69, 79, 89, 99),
  segVal = c(1, 2, 2, 1, 2, 1, 2, 3, 2, 2)
)


df <- dplyr::tibble(
  start = c(1, 10, 20, 30, 40, 50, 60),
  end = c(9, 19, 29, 39, 49, 59, 69),
  segVal = c(2, 1, 3, 3, 3, 2, 2)
)

df

join_segments(df)
join_segments(df[c(1:2), ])
join_segments(df[c(1), ])
join_segments(df[c(2:4), ])
join_segments(df[c(5:8), ])


library(maftools)
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
laml = read.maf(maf = laml.maf)

library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
laml.tnm = trinucleotideMatrix(maf = laml, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
dim(laml.tnm$nmf_matrix)
library('NMF')
laml.sig = extractSignatures(mat = laml.tnm, n = 6, pConstant = 1e-9)
laml.se = signatureEnrichment(maf = laml, sig_res = laml.sig)
sum(laml.se$cf_sizes$N)

load(system.file("extdata", "toy_copynumber_prepare.RData",
                 package = "sigminer", mustWork = TRUE
))
# Extract copy number signatures
library(NMF)
sig <- sig_extract(cn_prepare$nmf_matrix, 5, nrun = 10)

# Methods 'consensus' and 'samples' are from NMF::predict()
get_groups(sig)
get_groups(sig, method = "exposure")
get_groups(sig, method = "consensus")
get_groups(sig, method = "consensus", match_consensus = TRUE)
get_groups(sig, method = "samples")

debug(get_groups)
undebug(get_groups)

## Test New feature

# Load copy number object
load(system.file("extdata", "toy_copynumber.RData",
                 package = "sigminer", mustWork = TRUE
))
# Use method designed by Wang, Shixiang et al.
cn_prepare <- sig_tally(cn, method = "W", feature_setting = CN.features[1:50])
cn_prepare <- sig_tally(cn, method = "W", feature_setting = CN.features)

cn_prepare$features$BoChr
cn_prepare$features$NChrV
cn_prepare$features$NC50
cn_prepare$components
cn_prepare$parameters
debug(sig_tally)

library(NMF)
sigs <- sig_extract(cn_prepare$nmf_matrix, n_sig = 3, pConstant = 1e-9)
show_sig_profile(sigs, method = "W", normalize = "feature", style = "cosmic")
show_sig_profile(sigs, method = "W", normalize = "feature", style = "cosmic", filters = c("NChrV", "BoChr"))
show_sig_profile(sigs, method = "W", normalize = "feature", style = "cosmic", filters = c("BoChr"))

# Load copy number prepare object
load(system.file("extdata", "toy_copynumber_prepare.RData",
                 package = "sigminer", mustWork = TRUE
))
# Extract copy number signatures
library(NMF)
res <- sig_extract(cn_prepare$nmf_matrix, 5, nrun = 2)
show_sig_profile(res, normalize = "feature", style = "cosmic")
show_sig_profile(res, normalize = "feature", style = "cosmic", filters = c("bp10MB", "bpchrarm"))
show_sig_profile(laml.sig$signatures, mode = "mutation", style = "default",rm_panel_border = F)
show_sig_profile(laml.sig$signatures, mode = "mutation", style = "cosmic")


# Change name of function
#sed -i "" "s/sig_feed/sig_tally/g" `grep "sig_feed" -rl R/*`
#sed -i "" "s/sig_feed/sig_tally/g" `grep "sig_feed" -rl docs/*`
#sed -i "" "s/sig_feed/sig_tally/g" `grep "sig_feed" -rl man/*`
#sed -i "" "s/sig_feed/sig_tally/g" `grep "sig_feed" -rl tests/*`
#sed -i "" "s/sig_feed/sig_tally/g" `grep "sig_feed" -rl _pkgdown.yml`
