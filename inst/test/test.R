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


# Test DBS ----------------------------------------------------------------

maf = read_maf("../MatrixGenerator-Test/data/LUSC/tcga_lusc_from_dcc.maf")
mt_tally <- sig_tally(
  maf,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
  use_syn = TRUE, add_trans_bias = TRUE, mode = "DBS"
)



z2 = query[, search_DBS(.SD),
      by = Tumor_Sample_Barcode]

# query %>%
#   dplyr::as_tibble() %>%
#   dplyr::group_by(Tumor_Sample_Barcode) %>%
#   tidyr::nest() %>%
#   dplyr::mutate(data = purrr::map(.data$data, search_DBS)) %>%
#   tidyr::unnest("data") %>%
#   data.table::as.data.table()


# Test INDEL --------------------------------------------------------------
debug(generate_matrix_INDEL)
mt_tally <- sig_tally(
  maf,
  ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
  use_syn = TRUE, mode = "ID"
)
