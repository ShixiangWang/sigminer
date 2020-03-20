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


show_sig_feature_corrplot(tidy_data.seqz.feature, p_val = 0.05)

# > YAPSA::LCD
function (in_mutation_catalogue_df, in_signatures_df, in_per_sample_cutoff = 0)
{
  signatures_matrix <- as.matrix(in_signatures_df)
  out_exposures_df <- data.frame()
  G <- diag(dim(signatures_matrix)[2])
  H <- rep(0, dim(signatures_matrix)[2])
  for (i in seq_len(ncol(in_mutation_catalogue_df))) {
    temp_fractions <- lsei::lsei(a = signatures_matrix, b = in_mutation_catalogue_df[,
                                                                                     i], e = G, f = H)
    temp_exposures_vector <- round(temp_fractions, digits = 6)
    names(temp_exposures_vector) <- names(in_signatures_df)
    rel_exposures_vector <- temp_exposures_vector/sum(temp_exposures_vector)
    deselect_ind <- which(rel_exposures_vector < in_per_sample_cutoff)
    temp_exposures_vector[deselect_ind] <- 0
    out_exposures_df[seq(1, dim(signatures_matrix)[2], 1),
                     i] <- temp_exposures_vector
    rm(temp_fractions)
  }
  colnames(out_exposures_df) <- colnames(in_mutation_catalogue_df)
  rownames(out_exposures_df) <- colnames(in_signatures_df)
  return(out_exposures_df)
}
