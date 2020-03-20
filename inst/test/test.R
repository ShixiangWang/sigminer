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

library(NMF)
sig = sig_extract(cn_tally_W$nmf_matrix, n_sig = 2, pConstant = 1e-13)
save(sig, file = "inst/extdata/toy_copynumber_signature_by_W.RData")


## Signature identification for SBS
library(sigminer)
library(NMF)

laml.maf <- system.file("extdata", "tcga_laml.maf.gz", package = "maftools", mustWork = TRUE)
laml <- read_maf(maf = laml.maf)
if (require("BSgenome.Hsapiens.UCSC.hg19")) {
  mt_prepare <- sig_tally(
    laml,
    ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
    prefix = "chr", add = TRUE, useSyn = TRUE
  )
} else {
  message("Please install package 'BSgenome.Hsapiens.UCSC.hg19' firstly!")
}

sig = sig_extract(mt_prepare$nmf_matrix, n_sig = 3, pConstant = 1e-13)
get_sig_exposure(sig)


# Change name of function
# sed -i "" "s/sig_feed/sig_tally/g" `grep "sig_feed" -rl R/*`
# sed -i "" "s/sig_feed/sig_tally/g" `grep "sig_feed" -rl docs/*`
# sed -i "" "s/sig_feed/sig_tally/g" `grep "sig_feed" -rl man/*`
# sed -i "" "s/sig_feed/sig_tally/g" `grep "sig_feed" -rl tests/*`
# sed -i "" "s/sig_feed/sig_tally/g" `grep "sig_feed" -rl _pkgdown.yml`


load(system.file("extdata", "asso_data.RData",
                 package = "sigminer", mustWork = TRUE
))


show_sig_feature_corrplot(tidy_data.seqz.feature, p_val = 1)


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
