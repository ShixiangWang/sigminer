# Load copy number object
cn <- readRDS("D:/data/pcawg_cn_obj.rds")

cn_tally <- sig_tally(cn, method = "T")

cn_tally$nmf_matrix %>% colnames()
show_catalogue(cn_tally, mode = "copynumber", method = "T", style = "cosmic")
