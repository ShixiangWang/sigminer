# Load copy number object
load(system.file("extdata", "toy_copynumber.RData",
                 package = "sigminer", mustWork = TRUE
))

# Prepare copy number signature analysis
options(sigminer.sex = "male", sigminer.copynumber.max = 10)

cn_prepare <- sig_derive(cn, method = "W")
cn_prepare <- sig_derive(cn, method = "W", cores = 2)

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
library('NMF')
laml.sig = extractSignatures(mat = laml.tnm, n = 3, pConstant = 1e-9)
laml.se = signatureEnrichment(maf = laml, sig_res = laml.sig)


load(system.file("extdata", "toy_copynumber_prepare.RData",
                 package = "sigminer", mustWork = TRUE
))
# Extract copy number signatures
library(NMF)
sig <- sig_extract(cn_prepare$nmf_matrix, 2, nrun = 10)

# Methods 'consensus' and 'samples' are from NMF::predict()
get_groups(sig)
