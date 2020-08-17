pcawg_cn <- readRDS(file = "../CNSigs/data/pcawg_copynumber.rds")
cn_obj = read_copynumber(pcawg_cn[!Chromosome %in% c("X", "Y")], max_copynumber = 1000L,
                          genome_build = "hg19", complement = FALSE,
                          genome_measure = "called", skip_annotation = TRUE)


# # Load copy number object
# load(system.file("extdata", "toy_copynumber.RData",
#                  package = "sigminer", mustWork = TRUE
# ))

# Use method designed by Wang, Shixiang et al.
cn_tally <- sig_tally(cn_obj, method = "T")

cn_tally$nmf_matrix %>% colnames()
show_catalogue(cn_tally, mode = "copynumber", method = "T", style = "cosmic")

debug(show_sig_profile)
