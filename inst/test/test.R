## Signature identification for absolute copy number profile

load(system.file("extdata", "toy_segTab.RData",
                 package = "sigminer", mustWork = TRUE
))
cn <- read_copynumber(segTabs,
                      seg_cols = c("chromosome", "start", "end", "segVal"),
                      genome_build = "hg19", complement = FALSE, verbose = TRUE
)

save(cn, file = "inst/extdata/toy_copynumber.RData")

# Load copy number object
load(system.file("extdata", "toy_copynumber.RData",
                 package = "sigminer", mustWork = TRUE
))

# Prepare copy number signature analysis
options(sigminer.sex = "male", sigminer.copynumber.max = 10)

cn_tally_M = sig_tally(cn, method = "M")
cn_tally_W = sig_tally(cn, method = "W")

save(cn_tally_M, file = "inst/extdata/toy_copynumber_tally_M.RData")
save(cn_tally_W, file = "inst/extdata/toy_copynumber_tally_W.RData")


load(file = "inst/extdata/toy_copynumber_tally_W.RData")
library(NMF)
sig = sig_extract(cn_tally_W$nmf_matrix, n_sig = 2, pConstant = 1e-13)
save(sig, file = "inst/extdata/toy_copynumber_signature_by_W.RData")


sig2 = sig_extract(cn_tally_W$nmf_matrix, n_sig = 2, pConstant = 1e-13, optimize = TRUE)

identical(sig$Signature, sig2$Signature)
identical(sig$Signature.norm, sig2$Signature.norm)

sum((apply(sig$Signature[22:28,], 2, function(x) x/sum(x)) %*% sig$Exposure - t(cn_tally_W$nmf_matrix)[22:28, ])^2)
sum((apply(sig2$Signature[22:28,], 2, function(x) x/sum(x)) %*% sig2$Exposure - t(cn_tally_W$nmf_matrix)[22:28, ])^2)

expo1 = get_sig_exposure(sig)
expo2 = get_sig_exposure(sig2)

plot(expo1$Sig1, expo2$Sig1)
plot(expo1$Sig2, expo2$Sig2)

cor(expo1$Sig1, expo2$Sig1)
cor(expo1$Sig2, expo2$Sig2)

## Signature identification for SBS
library(sigminer)
library(NMF)

laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools", mustWork = TRUE)
laml <- read_maf(maf = laml.maf)
if (require("BSgenome.Hsapiens.UCSC.hg19")) {
  mt_tally <- sig_tally(
    laml,
    ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
    prefix = "chr", add = TRUE, useSyn = TRUE
  )
} else {
  message("Please install package 'BSgenome.Hsapiens.UCSC.hg19' firstly!")
}

sig = sig_extract(mt_tally$nmf_matrix, n_sig = 3, pConstant = 1e-13)
get_sig_exposure(sig)


sig2 = sig_extract(mt_tally$nmf_matrix, n_sig = 3, pConstant = 1e-13, optimize = TRUE)

identical(sig$Signature, sig2$Signature)
identical(sig$Signature.norm, sig2$Signature.norm)

sum((apply(sig$Signature, 2, function(x) x/sum(x)) %*% sig$Exposure - t(mt_tally$nmf_matrix))^2)
sum((apply(sig2$Signature, 2, function(x) x/sum(x)) %*% sig2$Exposure - t(mt_tally$nmf_matrix))^2)

expo1 = get_sig_exposure(sig)
expo2 = get_sig_exposure(sig2)

plot(expo1$Sig1, expo2$Sig1)
plot(expo1$Sig2, expo2$Sig2)

cor(expo1$Sig1, expo2$Sig1)
cor(expo1$Sig2, expo2$Sig2)

show_sig_profile(mt_tally$nmf_matrix[1, , drop = FALSE] %>% t(), mode = "SBS",
                 style = "cosmic", normalize = "raw",
                 check_sig_names = FALSE)

show_sig_profile(mt_tally$nmf_matrix[1:2, ] %>% t(), mode = "SBS",
                 style = "cosmic", normalize = "raw",
                 check_sig_names = FALSE)

show_cosmic_sig_profile(sig_index = 1, style = "cosmic")


# Change name of function
# sed -i "" "s/sig_feed/sig_tally/g" `grep "sig_feed" -rl R/*`
# sed -i "" "s/sig_feed/sig_tally/g" `grep "sig_feed" -rl docs/*`
# sed -i "" "s/sig_feed/sig_tally/g" `grep "sig_feed" -rl man/*`
# sed -i "" "s/sig_feed/sig_tally/g" `grep "sig_feed" -rl tests/*`
# sed -i "" "s/sig_feed/sig_tally/g" `grep "sig_feed" -rl _pkgdown.yml`


load(system.file("extdata", "asso_data.RData",
                 package = "sigminer", mustWork = TRUE
))

#?show_sig_feature_corrplot
show_sig_feature_corrplot(tidy_data.seqz.feature, p_val = 1)
show_sig_feature_corrplot(tidy_data.seqz.feature, breaks_count = NULL)
show_sig_feature_corrplot(tidy_data.seqz.feature, breaks_count = NA)


## Quantify exposure

## define raw data
W_prim <- matrix(c(1,2,3,4,5,6),ncol=2)
W_prim_df <- as.data.frame(W_prim)
W_df <- YAPSA:::normalize_df_per_dim(W_prim_df,2) # corresponds to the sigs

W <- as.matrix(W_df)
## 1. Simple case: non-negativity already in raw data

H_df <- as.data.frame(H) # corresponds to the exposures
V <- W %*% H # matrix multiplication
colnames(V) <- paste0("samp", 1:4)
V_df <- as.data.frame(V) # corresponds to the mutational catalogue
YAPSA:::LCD(V_df,W_df)

optimize_sig_exposure(V, W)
optimize_sig_exposure(V[, 1, drop = FALSE], W)
optimize_sig_exposure(V[, 1, drop = FALSE], W, type = "relative")

optimize_sig_exposure(V, W[, 1, drop = FALSE])
optimize_sig_exposure(V, sig_index = 1)


load(system.file("extdata", "toy_mutational_signature.RData",
                 package = "sigminer", mustWork = TRUE
))

get_sig_similarity(sig2)

expo_df <- get_sig_exposure(sig2)
expo_df

laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
laml <- read_maf(maf = laml.maf)
if (require("BSgenome.Hsapiens.UCSC.hg19")) {
  mt_tally <- sig_tally(
    laml,
    ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
    prefix = "chr", add = TRUE, useSyn = TRUE
  )
} else {
  message("Please install package 'BSgenome.Hsapiens.UCSC.hg19' firstly!")
}

expo_mat <- optimize_sig_exposure(t(mt_tally$nmf_matrix), sig_index = c(1, 10, 28), return_class = "data.table")

zz = dplyr::left_join(expo_df, expo_mat, by = "sample")
cor(zz$Sig2, zz$COSMIC_1)

zz[, -1] %>% as.matrix() %>% cor %>%
  pheatmap::pheatmap()


## Converting between genome- and exome-relative representations of signatures
# (Following from the signature extraction example presented above)
# Retrieve genome-derived signatures
library(sigfit)
data("counts_21breast")

mcmc_samples_extr <- extract_signatures(counts = counts_21breast,
                                        nsignatures = 2:7,
                                        iter = 1000,
                                        seed = 1756)

genome_signatures <- retrieve_pars(mcmc_samples_extr[[4]],
                                   par = "signatures")

# Apply exome mutational opportunities
exome_signatures <- convert_signatures(genome_signatures,
                                       opportunities_from = "human-genome",
                                       opportunities_to = "human-exome")

par(mfrow = c(2, 1))
plot_spectrum(genome_signatures$mean[4,],
              name = "Signature D, Genome-relative probabilities")
plot_spectrum(exome_signatures[4,],
              name = "Signature D, Exome-relative probabilities")


#> convert_signatures
function (signatures, opportunities_from = NULL, opportunities_to = NULL)
{
  if (is.null(opportunities_from) & is.null(opportunities_to)) {
    stop("Either 'opportunities_from' or 'opportunities_to' must be provided.")
  }
  ## to_matrix get the mean matrix
  signatures <- to_matrix(signatures)
  opportunities_from <- build_opps_matrix(nrow(signatures),
                                          ncol(signatures), opportunities_from)
  opportunities_to <- build_opps_matrix(nrow(signatures), ncol(signatures),
                                        opportunities_to)
  conv_sigs <- signatures
  for (i in 1:nrow(conv_sigs)) {
    conv_sigs[i, ] <- conv_sigs[i, ]/opportunities_from[i,
                                                        ] * opportunities_to[i, ]
    conv_sigs[i, ] <- conv_sigs[i, ]/sum(conv_sigs[i, ])
  }
  conv_sigs
}

#> sigfit:::build_opps_matrix
function (nsamples, ncat, opps)
{
  if (is.null(opps)) {
    matrix(1, nrow = nsamples, ncol = ncat)
  }
  else {
    if (is.character(opps)) {
      strand <- ncat == 192
      opps <- matrix(rep(human_trinuc_freqs(opps, strand),
                         nsamples), nrow = nsamples, byrow = TRUE)
    }
    else {
      opps <- as.matrix(opps)
    }
    opps/sum(opps) * nsamples
  }
}

#> sigfit:::human_trinuc_freqs
function (type = "human-genome", strand = FALSE)
{
  if (type == "human-genome") {
    freq <- c(1.14e+08, 6.6e+07, 14300000, 91200000, 1.05e+08,
              74600000, 15700000, 1.01e+08, 81700000, 67600000,
              13500000, 79300000, 1.11e+08, 87500000, 12500000,
              1.25e+08, 1.14e+08, 6.6e+07, 14300000, 91200000,
              1.05e+08, 74600000, 15700000, 1.01e+08, 81700000,
              67600000, 13500000, 79300000, 1.11e+08, 87500000,
              12500000, 1.25e+08, 1.14e+08, 6.6e+07, 14300000,
              91200000, 1.05e+08, 74600000, 15700000, 1.01e+08,
              81700000, 67600000, 13500000, 79300000, 1.11e+08,
              87500000, 12500000, 1.25e+08, 1.17e+08, 75700000,
              1.04e+08, 1.41e+08, 73100000, 95500000, 1.15e+08,
              1.13e+08, 64300000, 53600000, 85200000, 82700000,
              1.18e+08, 1.12e+08, 1.07e+08, 2.18e+08, 1.17e+08,
              75700000, 1.04e+08, 1.41e+08, 73100000, 95500000,
              1.15e+08, 1.13e+08, 64300000, 53600000, 85200000,
              82700000, 1.18e+08, 1.12e+08, 1.07e+08, 2.18e+08,
              1.17e+08, 75700000, 1.04e+08, 1.41e+08, 73100000,
              95500000, 1.15e+08, 1.13e+08, 64300000, 53600000,
              85200000, 82700000, 1.18e+08, 1.12e+08, 1.07e+08,
              2.18e+08)
  }
  else if (type == "human-exome") {
    freq <- c(1940794, 1442408, 514826, 1403756, 2277398,
              2318284, 774498, 2269674, 1740752, 1968596, 631872,
              1734468, 1799540, 1910984, 398440, 2024770, 1940794,
              1442408, 514826, 1403756, 2277398, 2318284, 774498,
              2269674, 1740752, 1968596, 631872, 1734468, 1799540,
              1910984, 398440, 2024770, 1940794, 1442408, 514826,
              1403756, 2277398, 2318284, 774498, 2269674, 1740752,
              1968596, 631872, 1734468, 1799540, 1910984, 398440,
              2024770, 1299256, 1166912, 1555012, 1689928, 978400,
              2119248, 2650754, 1684488, 884052, 1173252, 1993110,
              1251508, 1391660, 1674368, 1559846, 2850934, 1299256,
              1166912, 1555012, 1689928, 978400, 2119248, 2650754,
              1684488, 884052, 1173252, 1993110, 1251508, 1391660,
              1674368, 1559846, 2850934, 1299256, 1166912, 1555012,
              1689928, 978400, 2119248, 2650754, 1684488, 884052,
              1173252, 1993110, 1251508, 1391660, 1674368, 1559846,
              2850934)
  }
  else {
    stop("'type' must be either \"human-genome\" or \"human-exome\"")
  }
  if (strand) {
    rep(freq/2, 2)
  }
  else {
    freq
  }
}


# [1] "ACA>AAA" "ACC>AAC" "ACG>AAG" "ACT>AAT" "CCA>CAA" "CCC>CAC" "CCG>CAG" "CCT>CAT" "GCA>GAA" "GCC>GAC" "GCG>GAG"
# [12] "GCT>GAT" "TCA>TAA" "TCC>TAC" "TCG>TAG" "TCT>TAT" "ACA>AGA" "ACC>AGC" "ACG>AGG" "ACT>AGT" "CCA>CGA" "CCC>CGC"
# [23] "CCG>CGG" "CCT>CGT" "GCA>GGA" "GCC>GGC" "GCG>GGG" "GCT>GGT" "TCA>TGA" "TCC>TGC" "TCG>TGG" "TCT>TGT" "ACA>ATA"
# [34] "ACC>ATC" "ACG>ATG" "ACT>ATT" "CCA>CTA" "CCC>CTC" "CCG>CTG" "CCT>CTT" "GCA>GTA" "GCC>GTC" "GCG>GTG" "GCT>GTT"
# [45] "TCA>TTA" "TCC>TTC" "TCG>TTG" "TCT>TTT" "ATA>AAA" "ATC>AAC" "ATG>AAG" "ATT>AAT" "CTA>CAA" "CTC>CAC" "CTG>CAG"
# [56] "CTT>CAT" "GTA>GAA" "GTC>GAC" "GTG>GAG" "GTT>GAT" "TTA>TAA" "TTC>TAC" "TTG>TAG" "TTT>TAT" "ATA>ACA" "ATC>ACC"
# [67] "ATG>ACG" "ATT>ACT" "CTA>CCA" "CTC>CCC" "CTG>CCG" "CTT>CCT" "GTA>GCA" "GTC>GCC" "GTG>GCG" "GTT>GCT" "TTA>TCA"
# [78] "TTC>TCC" "TTG>TCG" "TTT>TCT" "ATA>AGA" "ATC>AGC" "ATG>AGG" "ATT>AGT" "CTA>CGA" "CTC>CGC" "CTG>CGG" "CTT>CGT"
# [89] "GTA>GGA" "GTC>GGC" "GTG>GGG" "GTT>GGT" "TTA>TGA" "TTC>TGC" "TTG>TGG" "TTT>TGT"
