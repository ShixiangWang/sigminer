# Load copy number object
load(system.file("extdata", "toy_copynumber.RData",
                 package = "sigminer", mustWork = TRUE
))

cn_tally <- sig_tally(cn, method = "T")

cn_tally$nmf_matrix %>% colnames()
show_catalogue(cn_tally, mode = "copynumber", method = "T", style = "cosmic")
